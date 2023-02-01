{
#samtools faidx outputs reference sequences in several lines. We just want one line per ZMW (otherwise the later paste command won't work)
#refseq.awk makes sure each line from the tmpREFpos corresponds with each line in tmpREFseq.
if($0 ~ "^>"){
#if the line starts with ">" it's a new ZMW.
#below if statement ensures that I'm not adding a new line before the very first ZMW. 
if(NR>1){
printf "\n"
}
} else{
printf $0;
}
}
END{
	#need to print an ending line break or wc -l will not count the number of lines properly
	printf "\n"
}
