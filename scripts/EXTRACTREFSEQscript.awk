#EXTRACTREFSEQscript.awk pulls out reference bases corresponding to reference mismatch positions using the final tag column that contains the entire alignment reference sequence (created from refseq.awk) and replaces it with just the mismatch reference bases as the last tag column.

BEGIN{
	#make sure array order matches order entries were added into array
	PROCINFO["sorted_in"] = "@ind_num_asc"
}
{

if($(NF-1) ~ "rp:"){
#read rp tag array (second to last column) if it exists, that includes the reference positions corresponding to mismatches into array refmmposarray. 
patsplit($(NF-1), refmmposarray, /[0-9]+/);

#extract corresponding ref base for each position in refmmposarray
for(i in refmmposarray){
	#need to account for starting ref POS ($4) by subtracting $4 and adding 1 to account for the 1-based coordinate system.
	refmmseqarray[i]=substr($NF, (refmmposarray[i]-$4+1), 1)
}
}

#print entire line EXCEPT the original final column. We just want to keep the bases corresponding to mismatches; this will be our new final column.
for(i=1;i<NF;i++){
printf "%s\t",$i
}

#add refmmseqarray as new final column
printf "rn:Z:";
for(i in refmmseqarray){
	if(i!=length(refmmseqarray)){
		printf "%s,",refmmseqarray[i]
	} else{
		printf "%s",refmmseqarray[i]
	}
};

printf "\n";

#reset arrays
delete refmmposarray;
delete refmmseqarray;
}
