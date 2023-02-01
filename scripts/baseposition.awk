#baseposition.awk pulls out the mismatch positions (w.r.t. reference and query) and mismatch sequences (w.r.t. query). 
#baseposition.awk also records full (reference) position of alignment in TMPrefpos. Format: RNAME:POS-(refmmcounter)
BEGIN{
	#make sure array order matches order entries were added into array
	PROCINFO["sorted_in"] = "@ind_num_asc";
}
{
	#print header
if($0 ~ "^@"){
	print $0
} else{
	#save CIGAR numbers in array BASEPOS
	patsplit($6, BASEPOS, /[0-9]+/);
	#save CIGAR non-numerical characters in array SYM
	#I is insertion; D is deletion; X is mismatch; N is skipped region from the reference; S is soft clipping; = is sequence match
	patsplit($6, SYM, /[IDXNS=]/);

	#create position counters
	querymmcounter=0;
	#save ref positions in refmmcounter. The last value of refmmcounter will be the last alignment position. The first alignment position is POS
	refmmcounter=($4-1);

	#reset arrays querymmposarray, refmmposarray, querymmseqarray, querymmqualarray
	delete querymmposarray;
	delete refmmposarray;
	querymmposarrayindex=1;
	refmmposarrayindex=1;
	delete querymmseqarray;
	delete querymmqualarray;

	#for character i in array BASEPOS, add as necessary to querymmcounter and refmmcounter. refer to sam alignment specifications for further details re: whether ISDN=X consumes reference and/or query position
	for(i in BASEPOS){
		if(SYM[i]=="I"){ 
			querymmcounter=querymmcounter+BASEPOS[i]  
		} else if(SYM[i]=="S"){
			querymmcounter=querymmcounter+BASEPOS[i]
		} else if(SYM[i]=="D"){ 
			refmmcounter=refmmcounter+BASEPOS[i]
		} else if(SYM[i]=="N"){ 
			refmmcounter=refmmcounter+BASEPOS[i]
		} else if(SYM[i]=="="){ 
			querymmcounter=querymmcounter+BASEPOS[i];
			refmmcounter=refmmcounter+BASEPOS[i]
		} else if(SYM[i]=="X"){ 
			#if it's a mismatch, we want to save the current querymmcounter and refmmcounter values in arrays querymmposarray and refmmposarray, respectively. If the corresponding ith element of BASEPOS is not 1, we need to save each mismatch position, so we need a for loop.
			for(k=1; k<=BASEPOS[i];k++){
				querymmposarray[querymmposarrayindex]=querymmcounter+k;
				refmmposarray[refmmposarrayindex]=refmmcounter+k;
				querymmposarrayindex++;
				refmmposarrayindex++
			}
			#set querymmcounter and refmmcounter equal to its previous value + the ith element of BASEPOS
			querymmcounter=querymmcounter+BASEPOS[i];
			refmmcounter=refmmcounter+BASEPOS[i]
		};
	}

	#extract corresponding query base (10th column) and quality (11th column) for each mismatch query position in querymmsequarray and querymmqualarray, respectively
	for(i in querymmposarray){
		querymmseqarray[i]=substr($10, querymmposarray[i], 1);
		querymmqualarray[i]=substr($11, querymmposarray[i], 1)
	};
	#Write to temporary file the full position of alignment in reference. TMPrefpos is a variable from process_subreads.sh. Each line will have:
	#chr:firstrefpos-lastrefpos
	print $3 ":" $4 "-" refmmcounter > TMPrefpos;

	#print original line with a trailing tab (tabs separate columns)
	printf "%s\t",$0;

	#add querymmposarray as new column (qp)
	#if there are no mismatch bases, then do not print the tag, since empty 'B' tag arrays are not allowed by SAM.
	if(length(querymmposarray)>0){
		printf "qp:B:I,"; ##See SAM spec section 1.5. 'B' tag arrays need to start with a letter indicating the type of variable. In this case we should use an unsigned 32-bit integer. That means large (up to the number 2^32) positive numbers.
		for(i in querymmposarray){
			if(i!=length(querymmposarray)){
				printf "%s,",querymmposarray[i]
			} else{
				#omit the trailing comma for the last entry in querymmposarray
				printf "%s",querymmposarray[i]
			}
		}
#print the tab separating columns
printf "\t";
	}

	#add querymmseqarray as new column (qn)
	printf "qn:Z:"; ##This needs to be 'Z' because 'B' is only for numeric arrays.
	for(i in querymmseqarray){
		if(i!=length(querymmseqarray)){
			printf "%s,",querymmseqarray[i]
		} else{
			printf "%s",querymmseqarray[i]
		}
	};
	printf "\t";
	#add querymmqualarray as new column (qq)
	printf "qq:Z:"; ##This needs to be 'Z' because 'B' is only for numeric arrays.
	for(i in querymmqualarray){
		if(i!=length(querymmqualarray)){
			printf "%s,",querymmqualarray[i]
		} else{
			printf "%s",querymmqualarray[i]
		}
	};
	printf "\t";
	#add refmmposarray as new column (rp), but not if it is empty, since empty 'B' type tag arrays are not allowed by SAM.
	if(length(refmmposarray)>0){
		printf "rp:B:I,";
		for(i in refmmposarray){
			if(i!=length(refmmposarray)){
				printf "%s,",refmmposarray[i]
			} else{
				printf "%s",refmmposarray[i]
			}
		}
	}
	
#print a line break
printf "\n";

}
}
