{
	#postccsfilter.awk runs on ccs.bam files (post ccs). 
	#the filter keeps only ZMWs that produced one forward consensus sequence and one reverse consensus sequence
#print header
if($0 ~ "^@"){
	print $0
} else{
	#file is organized by ZMW, so the ZMWs we wish to keep will have 2 lines (one forward; one reverse)
	#So, keep a line only if it is followed OR preceded by another line from the same ZMW
	#first, set zmwname=first column
	zmwname=$1; 
	#first column contains ZMW, but it also contains /rev and /fwd at the end. This means the first column of lines belonging to the same ZMW will not be identical, so we need to remove the /fwd and /rev. 
	sub("/rev","",zmwname); 
	sub("/fwd","",zmwname); 
	#now, we can compare zmwname to see if lines belong to the same ZMW
	#if the current line ZMW matches the previous line's ZMW, print the current line and the prior line
	#increase keptzmw counter (records number of ZMWs kept)
	if(zmwname==priorzmwname){
		keptzmw++;
		print priorline; 
		print $0
	#if the current line ZMW does not match the previous line's ZMW, then increase the totalzmw counter (counts how many ZMWs are in the file), name the current first column "priorzmwname" and the current line as "priorline"
	} else{
		totalzmw++;
		priorzmwname=zmwname; 
		priorline=$0 
	} 
}
}
END{
	#print number of ZMWs kept and total original number of ZMWs in file "ccsfilterstats"
printf(keptzmw"/"totalzmw "="keptzmw/totalzmw" ZMWs kept/total initial ZMWs\n") > ccsfilterstats
}
