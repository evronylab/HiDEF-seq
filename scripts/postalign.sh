#!/bin/bash

#file input is unsorted and unfiltered bam file post pbmm2 alignment. 
#Output file is aligned bam file. Script will also output pbmm2filterstats

#DEFINITIONS
#create temporary file for ref positions
TMPrefpos="$(mktemp TMPrefposXXX)"

#create temporary file for filtered sam
TMPfilteredsam="$(mktemp TMPfilteredsamXXX)"

#create temporary file for header
TMPheader="$(mktemp TMPheaderXXX)"

#create temporary file for filtered sam without the header
TMPfilteredsamnohead="$(mktemp TMPfilteredsamnoheadXXX)"

#create temporary file for ref sequences
TMPrefseq="$(mktemp TMPrefseqXXX)"

#create temporary file for finalrefoutput
TMPfinalrefoutput="$(mktemp TMPfinalrefoutputXXX)"

#create temporary file for line QC (make sure # lines in TMPfilteredsamnohead = # lines in TMPrefseq)
TMPlineqc="$(mktemp TMPlineqcXXX)"

#convert $1 (input unsorted post pbmm2 bam file) to sam file (keep header) -> pbmm2filter.awk keeps only ZMWs with 1 forward alignment and 1 reverse alignment -> convert to bam -> sort bam file -> samtools view (keep header) -> baseposition.awk pulls and saves mismatch positions and bases, save as file TMPfilteredsam. 
#baseposition.awk also records full (reference) position of alignment in TMPrefpos. Format: chr:POS-(refmmcounter)
#last cut command saves the job id #
cmd=$(mktemp $TMPoutputFullpath/.XXXXXX)
echo "$samtoolsbin view -h -@4 $1 | awk -v minoverlap=$minoverlap -v pbmm2filterstats=$pbmm2filterstats -f $scriptdir/pbmm2filter.awk | $samtoolsbin view -@4 -b | $samtoolsbin sort -@4 -m 3G | $samtoolsbin view -h -@4 | awk -v TMPrefpos=$TMPrefpos -v TMPfilteredsam=$TMPfilteredsam -f $scriptdir/baseposition.awk > $TMPfilteredsam" > $cmd
chmod +x $cmd
jid1=$(sbatch $slurm_add_options -t 20:00:00 --mem=64G -c 4 --parsable --wrap "$HIDEF_SINGULARITY_WRAPPER $cmd; rm $cmd")

#save header of TMPfilteredsam in temporary file TMPheader
cmd=$(mktemp $TMPoutputFullpath/.XXXXXX)
echo "$samtoolsbin view -H $TMPfilteredsam > $TMPheader" > $cmd
chmod +x $cmd
jid2=$(sbatch $slurm_add_options --dependency=afterok:$jid1 -t 00:30:00 --mem=2G -c 1 --parsable --wrap "$HIDEF_SINGULARITY_WRAPPER $cmd; rm $cmd")

#save SAM data without header from TMPfilteredsam in temporary file TMPfilteredsamnohead
cmd=$(mktemp $TMPoutputFullpath/.XXXXXX)
echo "$samtoolsbin view $TMPfilteredsam > $TMPfilteredsamnohead" > $cmd
chmod +x $cmd
jid3=$(sbatch $slurm_add_options --dependency=afterok:$jid2 -t 05:00:00 --mem=2G -c 1 --parsable --wrap "$HIDEF_SINGULARITY_WRAPPER $cmd; rm $cmd")

#refseq.awk pulls out reference sequences corresponding with the positions listed in TMPrefpos, reformats the file so there is 1 line per reference sequence
cmd=$(mktemp $TMPoutputFullpath/.XXXXXX)
echo "$samtoolsbin faidx $genomereffasta -r $TMPrefpos | awk -v TMPrefseq=$TMPrefseq -f $scriptdir/refseq.awk > $TMPrefseq" > $cmd
chmod +x $cmd
jid4=$(sbatch $slurm_add_options --dependency=afterok:$jid3 -t 20:00:00 --mem=64G -c 4 --parsable --wrap "$HIDEF_SINGULARITY_WRAPPER $cmd; rm $cmd")

#Paste TMPrefseq as new column into BAM file without BAM header (TMPfilteredsamnohead) -> EXTRACTREFSEQscript.awk pulls out reference bases corresponding to reference mismatch positions -> save output as TMPfinalrefoutput
cmd=$(mktemp $TMPoutputFullpath/.XXXXXX)
echo "paste $TMPfilteredsamnohead $TMPrefseq | awk -f $scriptdir/EXTRACTREFSEQscript.awk > $TMPfinalrefoutput" > $cmd
chmod +x $cmd
jid5=$(sbatch $slurm_add_options --dependency=afterok:$jid4 -t 20:00:00 --mem=4G -c 1 --parsable --wrap "$HIDEF_SINGULARITY_WRAPPER $cmd; rm $cmd")

#QC step: save line number of TMPfilteredsamnohead as linesamnoheader -> save line number of TMPrefseq as linerefseq -> save line number of TMPfinalrefoutput as linefinalref -> if line numbers aren't equivalent, output 'FAIL' to TMPlineqc. If line numbers are equivalent, echo 'OK' to TMPlineqc
cmd=$(mktemp $TMPoutputFullpath/.XXXXXX)
echo "linesamnoheader=\$(wc -l $TMPfilteredsamnohead | cut -f1 -d' '); linerefseq=\$(wc -l $TMPrefseq | cut -f1 -d' '); linefinalref=\$(wc -l $TMPfinalrefoutput | cut -f1 -d' '); if [ \"\$linesamnoheader\" -ne \"\$linerefseq\" ] || [ \"\$linesamnoheader\" -ne \"\$linefinalref\" ]; then echo 'FAIL' > $TMPlineqc; else echo 'OK' > $TMPlineqc; fi" > $cmd
chmod +x $cmd
jid6=$(sbatch $slurm_add_options --dependency=afterok:$jid5 -t 20:00:00 --mem=64G -c 4 --parsable --wrap "$HIDEF_SINGULARITY_WRAPPER $cmd; rm $cmd")

#loop to check if lineqc passes (i.e. OK); keeps looping until file TMPlineqc exists
#if line numbers aren't equivalent, exit the script
while [ ! -s $TMPlineqc ]
do
  sleep 10 #sleeps 10 sec
done
sleep 10 # sleeps 10 sec to give time for prior job to write TMPlineqc file
if [[ $(< $TMPlineqc) = "FAIL" ]]; then
	exit 1
fi

#Combine TMPheader and TMPfinalrefoutput into new BAM file.
cat $TMPheader $TMPfinalrefoutput | $samtoolsbin view -b > $finaloutputfile
