#post pbmm2 filter keeps only ZMWs with 1 forward primary alignment and 1 reverse primary alignment. Must also align to the same chromosome and overlap must be at least 90% between fwd and rev alignments
BEGIN{
	#need to set first element of line to a blank. otherwise, length(line) outputs an error
	line[1]=""
}
{
#preserve header
if($0 ~ "^@"){
	print $0
} else{
	#start readlength, which records the length of the alignment in the reference space.
	readlength=0;
	#save CIGAR (column 6) numbers in array BASEPOS
	patsplit($6, BASEPOS, /[0-9]+/);
	#save CIGAR non-numerical characters in array SYM
	#I is insertion; D is deletion; X is mismatch; N is skipped region from the reference; S is soft clipping; = is sequence match
	patsplit($6, SYM, /[IDXNS=]/);

	#save ref positions in readlength. The last value of readlength will be the last alignment position. The first alignment position is POS (column 4)
	for(i in BASEPOS){
		if(SYM[i]=="D" || SYM[i]=="X" || SYM[i]=="N" || SYM[i]=="="){ 
			readlength=readlength+BASEPOS[i]
		}
	}

	#set variable zmwname so we can compare different lines' column 1. column 1 is formatted "movie#/ZMW#/ccs/" then "fwd" or "rev"
	zmwname=$1; sub("/fwd","",zmwname); sub("/rev","",zmwname);

	#if the ZMW changes (i.e. is different from the previous ZMW), then do the below. This is also run for the first line
	if(zmwname!=priorzmwname){
		#array line's elements will be all the lines (i.e. alignments) corresponding to a particular ZMW. k is the number of alignments belonging to that particular ZMW
		k=length(line);
		#totalzmw counts the number of unique ZMWs in the file
		totalzmw++;

		#below creates counters to count the number of ZMWs with supplementary alignments. this is for the filter stats file pbmm2filterstats
		#count number of ZMWs with 1, 2, >2 supplementary alignments
		if(numsupp==1){
			onesupp++
		};
		if(numsupp==2){
			twosupp++
		};
		if(numsupp>2){
			greattwosupp++
		};

		#below is the actual filtration step
		#if there is only 1 fwd alignment and 1 rev alignment (impossible to have a supplementary alignment without a primary alignment, so we know these will be primary) and if the fwd chromosome alignment (column 3) equals the rev chromosome alignment
		if(numfwd==1 && numrev==1 && fwdchr==revchr){
			#below covers all the potential cases. overlap formula will differ depending on the first (start) and last positions of the fwd and rev alignments
			if(startrev >= startfwd){
				if(lastrev < lastfwd){
					overlap = lastrev - startrev + 1;
				}
				if(lastrev >= lastfwd){
					overlap = lastfwd - startrev + 1;
				}
			};
			if(startrev < startfwd){
				if(lastrev < lastfwd){
					overlap = lastrev - startfwd + 1;
				}
				if(lastrev >= lastfwd){
					overlap = lastfwd - startfwd + 1;
				}
			};

			#Print zmw lines if overlaps are greater than or equal to the minoverlap (set in process_subreads.sh)
			if(overlap/lengthfwd >= minoverlap && overlap/lengthrev >= minoverlap){
if(fwdstrand!=revstrand){
					#if the fwd and rev alignments meet all these cases, then print their corresponding lines
					for(j=1;j<=k;j++){print line[j]};
					#increase the keptzmw counter by 1. this records the number of ZMWs that are kept
					keptzmw++;
}else{
#Count number of ZMWs that pass all filters but align to the same strand.
samestrand++;
}
			}
		}
			
		##reset variables
		#delete the lines from the array line
		for(j=1;j<=k;j++){
			delete line[j]
		};
		numprim=0;
		numsupp=0;
		numfwd=0; 
		startfwd=0;
		lastfwd=0;
		m=1;
		numrev=0;
		startrev=0;
		lastrev=0;
		overlap=0;
		fwdchr=0;
		revchr=0;
		fwdstrand=0;
		revstrand=0;
	}

	##The following code is performed for every line.
	#if it's a fwd alignment, increase numfwd and record the last alignment position as lastfwd. record the first alignment position as startfwd. record alignment length as lengthfwd. fwdchr is the chromosome alignment
	#Calculation of lastfwd accounts for 1-based coordinates (need to subtract 1 from start position)
	#Also determine alignment strand; '1' = fwd, '2' = rev.

	if($2==0 || $2==2048 ){
		strand=1
	}else{
		strand=2
	}

	if(match($1,"/fwd")!=0){
		numfwd++;
		lastfwd=$4 - 1 + readlength;
		startfwd=$4;
		lengthfwd=readlength;
		fwdchr=$3;
		fwdstrand=strand;
	}
	#do the equivalent for a reverse alignment
	else if(match($1,"/rev")!=0){
		numrev++;
		lastrev=$4 - 1 + readlength;
		startrev=$4;
		lengthrev=readlength;
		revchr=$3;
		revstrand=strand;
	}

	##Count primary and supplementary alignments
	#0 signals primary alignment; #16 signals reverse primary alignment
	if($2 == 0 || $2 == 16){
		numprim++;
	}
	#2048 signals supplementary alignment; #2064 signals reverse supplementary alignment; 
	if($2 == 2048 || $2 == 2064){
		numsupp++;
	}
	#add the current line to array "line" as element m
	line[m]=$0;
	#increase m, so the next line will be added as a the next element in array line
	m++;
	#rename current zmwname as priorzmwname
	priorzmwname=zmwname
}
}
END{
	#this deals w/ case where the final lines have not yet been processed.
	if(length(line)>0){
		k=length(line);
		totalzmw++;
		
		if(numsupp==1){
			onesupp++
		};
		if(numsupp==2){
			twosupp++
		};
		if(numsupp>2){
			greattwosupp++
		};

		if(numfwd==1 && numrev==1 && fwdchr==revchr){
			if(startrev >= startfwd){
				if(lastrev < lastfwd){
					overlap = lastrev - startrev + 1;
				}
				if(lastrev >= lastfwd){
					overlap = lastfwd - startrev + 1;
				}
			};
			if(startrev < startfwd){
				if(lastrev < lastfwd){
					overlap = lastrev - startfwd + 1;
				}
				if(lastrev >= lastfwd){
					overlap = lastfwd - startfwd + 1;
				}
			};

			#Print zmw lines if overlaps are greater than or equal to the minoverlap (set in process_subreads.sh)
			if(overlap/lengthfwd >= minoverlap && overlap/lengthrev >= minoverlap){
if(fwdstrand!=revstrand){
					#if the fwd and rev alignments meet all these cases, then print their corresponding lines
					for(j=1;j<=k;j++){print line[j]};
					#increase the keptzmw counter by 1. this records the number of ZMWs that are kept
					keptzmw++;
}else{
#Count number of ZMWs that pass all filters but align to the same strand.
samestrand++;
}
			}
		}
	}
	
	#Subtract one from totalzmw, since it is too high by 1 due to incrementing in the first line of the file.
	totalzmw=totalzmw-1;
	
	#print out filter stats into file pbmm2filterstats
	printf(keptzmw"/"totalzmw"="keptzmw/totalzmw" ZMWs kept/total initial ZMWs\n") > pbmm2filterstats;
	printf(totalzmw-onesupp-twosupp-greattwosupp"=number of ZMWs with zero supplementary alignments\n") >> pbmm2filterstats;
	printf((totalzmw-onesupp-twosupp-greattwosupp)/totalzmw"=fraction of total ZMWs with zero supplementary alignments\n") >> pbmm2filterstats;
	printf(onesupp"=number of ZMWs with one supplementary alignment\n") >> pbmm2filterstats;
	printf(onesupp/totalzmw"=fraction of total ZMWs with one supplementary alignment\n") >> pbmm2filterstats;
	printf(twosupp"=number of ZMWs with two supplementary alignments\n") >> pbmm2filterstats;
	printf(twosupp/totalzmw"=fraction of total ZMWs with two supplementary alignments\n") >> pbmm2filterstats;
	printf(greattwosupp"=number of ZMWs with more than two supplementary alignments\n") >> pbmm2filterstats;
	printf(greattwosupp/totalzmw"=fraction of total ZMWs with more than two supplementary alignments\n") >> pbmm2filterstats
	printf(samestrand"=number of ZMWs passing all filters but with both alignments aligning to the same strand\n") >> pbmm2filterstats;
	printf(samestrand/totalzmw"=fraction of total ZMWs passing all filters but with both alignments aligning to the same strand\n") >> pbmm2filterstats
}
