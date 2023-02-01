#!/bin/bash

#Usage: process_subreads.sh [config.yaml]
# Prerequisites: ccs, pbmm2, lima, pbbam, samtools >=v1.12, samtools v1.10, SLURM environment
# Uses [config.yaml] file per documentation

##LOAD CONFIGURATION
export configyaml=$(readlink -f $1)
export subreads=`grep "subreads_filename: " $configyaml | cut -d " " -f 2`
export outputfileprefix=`grep "ccs_BAM_prefix: " $configyaml | cut -d " " -f 2`
export outputdrct=`grep "process_subreads_output_path: " $configyaml | cut -d " " -f 2`
export genomereffasta=`grep "fastaref: " $configyaml | cut -d " " -f 2`
export genomerefmmi=`grep "genomemmiindex: " $configyaml | cut -d " " -f 2`
export samtoolsbin_new=`grep "samtools_bin: " $configyaml | cut -d " " -f 2`
export samtoolsbin=`grep "samtools_1.10_bin: " $configyaml | cut -d " " -f 2`
export scriptdir=`grep "HiDEFpipeline: " $configyaml | cut -d " " -f 2`
export condabase_script=`grep "condabase_script: " $configyaml | cut -d " " -f 2`
export conda_pbbioconda_env=`grep "conda_pbbioconda_env: " $configyaml | cut -d " " -f 2`
export slurm_add_options=`grep "slurm_add_options: " $configyaml | cut -d " " -f 2-`
export inputdrct=`dirname $(readlink -f $configyaml)`
#ccschunks is number of chunks used in pbccs.
export ccschunks=`grep "ccschunks: " $configyaml | cut -d " " -f 2`
#minccsrq is minimum predicted quality of ccs consensus. parameter is used in pbccs. We chose 0.99 because this is the default, recommended setting.
export minccsrq=`grep "minccsrq: " $configyaml | cut -d " " -f 2`
#minoverlap is used in pbmm2filter.awk. It's the minimum reciprocal overlap between the forward primary alignment and reverse primary alignment. We decided 90% overlap was sufficient to permit comparison of the forward and reverse alignments.
export minoverlap=`grep "minoverlap: " $configyaml | cut -d " " -f 2`
#if removetempdir is set to true, pipeline will delete the TMPoutput directory. If it's set to anything else, pipeline will not delete the TMPoutput directory
export removetempdir=`grep "removetempdir: " $configyaml | cut -d " " -f 2`

#make output directory and output logs sub-directory
mkdir -p $outputdrct
mkdir -p $outputdrct/logs

#change to input directory
cd $inputdrct

#make temporary output directory
export TMPoutput=$(mktemp -d TMPoutputXXX)
export TMPoutputFullpath=$(realpath $TMPoutput)

#Make temporary barcodes.fasta file from config file for lima demultiplexing, and also load barcode names into bash array
awk '{if($0~/barcodes:/){startprint=1}else if(startprint==1 && $1=="-"){split($2,arr,":");print ">" arr[1];print arr[2]}else if(startprint==1 && $1!="-"){exit}}' $configyaml > $TMPoutputFullpath/barcodes.fasta
export barcodes=$TMPoutputFullpath/barcodes.fasta

mapfile -t barcodenames < <(awk '{if($0~/barcodes:/){startprint=1}else if(startprint==1 && $1=="-"){split($2,arr,":");print arr[1]}else if(startprint==1 && $1!="-"){exit}}' $configyaml)

numbarcodesminusone=`expr ${#barcodenames[@]} - 1`

#Load arrays with sample names
mapfile -t samplenames < <(awk '{if($0~/samplenames:/){startprint=1}else if(startprint==1 && $1=="-"){print $2}else if(startprint==1 && $1!="-"){exit}}' $configyaml)

#switch to temporary output directory
cd $TMPoutput

##STEP ONE: Count ZMWs (for QC purposes)
cmd=$(mktemp $TMPoutputFullpath/.XXXXXX)
echo "source $condabase_script; conda activate $conda_pbbioconda_env; pbindexdump --format cpp $subreads.pbi | grep holeNumber | tr ',' '\n' | sed 's/.*{//;s/}.*//' | uniq | wc -l > $inputdrct/$TMPoutput/subreads.zmwcount.txt" > $cmd
chmod +x $cmd
zmwrawcountjob=$(sbatch $slurm_add_options -t 12:00:00 --mem=64G -c 4 --parsable --wrap "$HIDEF_SINGULARITY_WRAPPER $cmd; rm $cmd")

##STEP TWO: Filter to remove any subread with cx tag != 3 (bad adaptor subreads), often caused by single hairpin artifacts. Then count ZMWs again after filtering.
cmd=$(mktemp $TMPoutputFullpath/.XXXXXX)
echo "$samtoolsbin_new view -b -@8 -d cx:3 $subreads > $inputdrct/$TMPoutput/$outputfileprefix.subreads.cxfiltered.bam; source $condabase_script; conda activate $conda_pbbioconda_env; pbindex $inputdrct/$TMPoutput/$outputfileprefix.subreads.cxfiltered.bam; pbindexdump --format cpp $inputdrct/$TMPoutput/$outputfileprefix.subreads.cxfiltered.bam.pbi | grep holeNumber | tr ',' '\n' | sed 's/.*{//;s/}.*//' | uniq | wc -l > $inputdrct/$TMPoutput/subreads.cxfiltered.zmwcount.txt" > $cmd
chmod +x $cmd
filtercxjob=$(sbatch $slurm_add_options -t 30:00:00 --mem=64G -c 8 --parsable --wrap "$HIDEF_SINGULARITY_WRAPPER $cmd; rm $cmd")

##STEP THREE: PBCCS --by-strand to create forward and reverse consensus sequences
#Creates ccs chunked bam files ($outputfileprefix.ccs.$i.bam)

#declare the array ccsrunjobids, which will contain the slurm job ids for each ccs chunk submitted
declare -a ccsrunjobids=()

#for loop will submit each of the ccs chunk jobs
#for each chunk: make temporary ccs directory (TMPccsoutput) --> activate conda (so you can run ccs) --> run ccs --by-strand, outputting all files to TMPccsoutput--> move ccs_report.txt (which contains info about ccs run), chunked ccs bam file and .bam.pbi file from TMPccsoutput to TMPoutput --> delete TMPccsoutput
for i in `seq 1 $ccschunks | xargs`
do
cmd=$(mktemp $TMPoutputFullpath/.XXXXXX)
echo "mkdir TMPccsoutput$i; cd TMPccsoutput$i; source $condabase_script; conda activate $conda_pbbioconda_env; ccs --by-strand --min-rq $minccsrq --top-passes 0 --report-file $inputdrct/$TMPoutput/$outputfileprefix.ccsreport.chunk${i}.txt --report-json $inputdrct/$TMPoutput/$outputfileprefix.ccsreport.chunk${i}.json --metrics-json $inputdrct/$TMPoutput/$outputfileprefix.ccsmetrics.chunk${i}.json $inputdrct/$TMPoutput/$outputfileprefix.subreads.cxfiltered.bam $outputfileprefix.ccs.$i.bam --chunk $i/$ccschunks -j 8; mv $outputfileprefix.ccs.* $inputdrct/$TMPoutput; cd $inputdrct/$TMPoutput; rm -rf TMPccsoutput$i" > $cmd
chmod +x $cmd
ccsrunjobids[i]=$(sbatch $slurm_add_options -t 24:00:00 --dependency=afterok:$filtercxjob --mem=64G -c 8 --parsable --wrap "$HIDEF_SINGULARITY_WRAPPER $cmd; rm $cmd")
done

##STEP FOUR: merge post ccs chunks using pbmerge
#Also save total number of ZMWs.
#create variable postccsjobids, which contains contents of array ccsrunjobids with : delimiter
#var contains contents of ccsrunjobids with : delimiters and surrounding {}
export var=$(IFS=:; echo "{${ccsrunjobids[*]}}")
#get rid of surrounding { }. below gets rid of last character } and saves the result as postccsjobid
export postccsjobid=${var%?}
#below gets rid of first character { and saves the result as postccsjobids
export postccsjobids="${postccsjobid:1}"

cmd=$(mktemp $TMPoutputFullpath/.XXXXXX)
echo "source $condabase_script; conda activate $conda_pbbioconda_env; pbmerge -o $inputdrct/$TMPoutput/$outputfileprefix.ccs.bam `seq -f $inputdrct/$TMPoutput/$outputfileprefix.ccs.%g.bam 1 $ccschunks | xargs`; pbindexdump --format cpp $inputdrct/$TMPoutput/$outputfileprefix.ccs.bam.pbi | grep holeNumber | tr ',' '\n' | sed 's/.*{//;s/}.*//' | uniq | wc -l > ccs.zmwcount.txt" > $cmd
chmod +x $cmd
mergepostccsjobid=$(sbatch $slurm_add_options --dependency=afterok:$postccsjobids -t 24:00:00 --mem=64G -c 4 --parsable --wrap "$HIDEF_SINGULARITY_WRAPPER $cmd; rm $cmd")

##STEP FIVE: Demultiplexing with lima
#Note: lima automatically detects the number of available threads.
#Outputs: $outputfileprefix.ccs.demux.barcode--barcode.bam for each sample.
cmd=$(mktemp $TMPoutputFullpath/.XXXXXX)
echo "source $condabase_script; conda activate $conda_pbbioconda_env; lima --ccs --same --split-named --min-score 80 --min-end-score 50 --min-ref-span 0.75 --min-scoring-regions 2 $inputdrct/$TMPoutput/$outputfileprefix.ccs.bam $barcodes $inputdrct/$TMPoutput/$outputfileprefix.ccs.demux.bam" > $cmd
chmod +x $cmd
limajobid=$(sbatch $slurm_add_options --dependency=afterok:$mergepostccsjobid -t 12:00:00 --mem=64G -c 8 --parsable --wrap "$HIDEF_SINGULARITY_WRAPPER $cmd; rm $cmd")

##STEP SIX: post pbccs filter (keeps only ZMWs with 1 forward consensus sequence and 1 reverse consensus sequence).
#Run for each demultiplexed sample

#declare array ccsfilterjobids. array entries will be job ids.
declare -A ccsfilterjobids=()

#for loop submits slurm jobs for each demultiplexed sample. Resulting job ids are stored in array ccsfilterjobids
#Before starting, rename each ccs BAM file from $outputfileprefix.ccs.demux.[barcodename]--[barcodename].bam to $outputfileprefix.[samplename].ccs.demux.[barcodename].bam
#Count number of ZMWs per sample. Then convert bam file to sam -> postccsfilter.awk filters out ZMWs that don't meet criteria -> convert to bam file as $outputfileprefix.[samplename].ccs.demux.[barcodename].postccsfilter.bam
#filter also outputs $outputfileprefix.[samplename].ccs.demux.[barcodename].bam.ccsfilterstats for each sample; ccsfilterstats file tells you how many ZMWs are filtered out.
for i in `seq 0 $numbarcodesminusone`
do
	export ccsfilterstats=$outputfileprefix.${samplenames[$i]}.ccs.demux.${barcodenames[$i]}.bam.ccsfilterstats
	cmd=$(mktemp $TMPoutputFullpath/.XXXXXX)
	echo "source $condabase_script; conda activate $conda_pbbioconda_env; mv $outputfileprefix.ccs.demux.${barcodenames[$i]}--${barcodenames[$i]}.bam $outputfileprefix.${samplenames[$i]}.ccs.demux.${barcodenames[$i]}.bam; mv $outputfileprefix.ccs.demux.${barcodenames[$i]}--${barcodenames[$i]}.bam.pbi $outputfileprefix.${samplenames[$i]}.ccs.demux.${barcodenames[$i]}.bam.pbi; pbindexdump --format cpp $outputfileprefix.${samplenames[$i]}.ccs.demux.${barcodenames[$i]}.bam.pbi | grep holeNumber | tr ',' '\n' | sed 's/.*{//;s/}.*//' | uniq | wc -l > ${barcodenames[$i]}.lima.zmwcount.txt; conda deactivate; $samtoolsbin view -h $outputfileprefix.${samplenames[$i]}.ccs.demux.${barcodenames[$i]}.bam | awk -v ccsfilterstats=$ccsfilterstats -f $scriptdir/postccsfilter.awk | $samtoolsbin view -b - > $outputfileprefix.${samplenames[$i]}.ccs.demux.${barcodenames[$i]}.postccsfilter.bam" > $cmd
	chmod +x $cmd
	ccsfilterjobids[${barcodenames[$i]}]=$(sbatch $slurm_add_options --dependency=afterok:$limajobid -t 12:00:00 --mem=32G -c 2 --parsable --wrap "$HIDEF_SINGULARITY_WRAPPER $cmd; rm $cmd")
done

##STEP SEVEN: pbmm2 (unsorted; preset ccs is used for ccs reads)
#declare array pbmm2jobids. This will contain the slurm pbmm2 job ids 
declare -A pbmm2jobids=()

#for loops submits pbmm2 jobs (dependencies correlate with previous job id array ccsfilterjobids)
#activate conda (so you can use pbmm2) -> pbmm2 -> outputs $outputfileprefix.[samplename].ccs.demux.[barcodename].postccsfilter.aligned.bam files.
for i in `seq 0 $numbarcodesminusone`
do
  cmd=$(mktemp $TMPoutputFullpath/.XXXXXX)
  echo "source $condabase_script; conda activate $conda_pbbioconda_env; pbmm2 align --log-level INFO -j 8 $genomerefmmi $inputdrct/$TMPoutput/$outputfileprefix.${samplenames[$i]}.ccs.demux.${barcodenames[$i]}.postccsfilter.bam $inputdrct/$TMPoutput/$outputfileprefix.${samplenames[$i]}.ccs.demux.${barcodenames[$i]}.postccsfilter.aligned.bam --preset CCS" > $cmd
  chmod +x $cmd
	pbmm2jobids[${barcodenames[$i]}]=$(jobid=${ccsfilterjobids[${barcodenames[$i]}]}; sbatch $slurm_add_options --dependency=afterok:$jobid -t 30:00:00 --mem=128G -c 8 --parsable --wrap "$HIDEF_SINGULARITY_WRAPPER $cmd; rm $cmd")
done

##STEP EIGHT: post pbmm2 filter + processing
#Outputs: $outputfileprefix.[samplename].ccs.demux.[barcodename].postccsfilter.aligned.final.bam and $outputfileprefix.[samplename].ccs.demux.[barcodename].aligned.bam.pbmm2filterstats
#declare array pbmm2filter. This will contain the pbmm2filter job ids
declare -A pbmm2filterjobids=()

#for loop submits post pbmm2 script (postalign.sh) for each of the chunked files
for i in `seq 0 $numbarcodesminusone`
do
  export pbmm2filterstats=$outputfileprefix.${samplenames[$i]}.ccs.demux.${barcodenames[$i]}.postccsfilter.aligned.bam.pbmm2filterstats;
  export finaloutputfile=$outputfileprefix.${samplenames[$i]}.ccs.demux.${barcodenames[$i]}.postccsfilter.aligned.final.bam
  
  cmd=$(mktemp $TMPoutputFullpath/.XXXXXX)
  echo "source $scriptdir/postalign.sh $inputdrct/$TMPoutput/$outputfileprefix.${samplenames[$i]}.ccs.demux.${barcodenames[$i]}.postccsfilter.aligned.bam; $samtoolsbin index $finaloutputfile; source $condabase_script; conda activate $conda_pbbioconda_env; pbindex $finaloutputfile" > $cmd
  chmod +x $cmd
  pbmm2filterjobids[${barcodenames[$i]}]=$(jobid=${pbmm2jobids[${barcodenames[$i]}]}; sbatch $slurm_add_options --dependency=afterok:$jobid -t 60:00:00 --mem=4G -c 1 --parsable --wrap "$HIDEF_SINGULARITY_WRAPPER $cmd; rm $cmd")
done

#create variable postpbmm2filterjobids, which contains contents of array pbmm2filterjobids with : delimiter
#var contains contents of pbmm2filterjobids with : delimiters and surrounding {}
export var=$(IFS=:; echo "{${pbmm2filterjobids[*]}}")
#get rid of surrounding { }. below gets rid of last character } and saves the result as postpbmm2filterjobid
export postpbmm2filterjobid=${var%?}
#below gets rid of first character { and saves the result as postpbmm2filterjobids
export postpbmm2filterjobids="${postpbmm2filterjobid:1}"

##STEP NINE: move final output files to output directory, and create file process_subreads.DONE in outputdrct with jobid for dependency in next part of the full pipeline.
cmd=$(mktemp $TMPoutputFullpath/.XXXXXX)
echo "mv *.ccsreport* $outputdrct/logs; mv *.ccsmetrics* $outputdrct/logs; mv *.lima.* $outputdrct/logs; mv *.ccsfilterstats $outputdrct/logs; mv *.final.bam* $outputdrct; mv *.pbmm2filterstats $outputdrct/logs; mv *.zmwcount.txt $outputdrct/logs" > $cmd
chmod +x $cmd
movefilesjobid=$(sbatch $slurm_add_options --dependency=afterok:$postpbmm2filterjobids -t 00:30:00 --parsable --wrap "$HIDEF_SINGULARITY_WRAPPER $cmd; rm $cmd")

echo $movefilesjobid > $outputdrct/process_subreads.DONE

#if $removetempdir is set to true, then delete tmp output directory
if [[ $removetempdir = true ]];
	then
	  cmd=$(mktemp $TMPoutputFullpath/.XXXXXX)
	  echo "cd $inputdrct; rm -rf $TMPoutput" > $cmd
	  chmod +x $cmd
		sbatch $slurm_add_options --dependency=afterok:$movefilesjobid -t 01:00:00 --wrap "$HIDEF_SINGULARITY_WRAPPER $cmd; rm $cmd"
fi
