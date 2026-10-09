# Baseline function preserved verbatim from bin/filterCalls.R at
# 7e37ada2664f07e9534e9dbd695b68a0fbcc22f3 (before observed allele keys).
# Original full-file SHA256: e4e6e67e371c25df9c13b840300df8ab09a9f8421c09cd35903a58b352bd69ab
# Scientific reference only; do not update alongside the implementation.
overlapsAny_bymcols <- function(query, subject, join_mcols = character(), ignore.strand = TRUE, overlap_adjacent_query_insertion = FALSE) {
	
	stopifnot(inherits(query, "GenomicRanges"), inherits(subject, "GenomicRanges"))
	
	nq <- length(query)
	ns <- length(subject)
	
	if(nq == 0){return(logical())} #Empty query
	if(ns == 0){return(rep(FALSE, nq))} #Empty subject
	
	if(length(join_mcols) > 0){
		q_mcols <- mcols(query)
		s_mcols <- mcols(subject)
		
		if(!all(join_mcols %in% colnames(q_mcols)) || !all(join_mcols %in% colnames(s_mcols))){
			stop("All `join_mcols` must exist as metadata columns in both query and subject.")
		}
		
		key_q <- query %>%
			as_tibble %>%
			select(all_of(join_mcols)) %>%
			as.list %>%
			interaction(drop = TRUE)
		
		key_s <- subject %>%
			as_tibble %>%
			select(all_of(join_mcols)) %>%
			as.list %>%
			interaction(drop = TRUE)
		
		keys_all <- factor(c(key_q,key_s))
		id_q <- as.integer(keys_all)[seq_len(nq)]
		id_s <- as.integer(keys_all)[nq + seq_len(ns)]
		
		# Make NA ids unique and disjoint between query and subject so they never overlap.
		if(anyNA(id_q)){
			idx <- id_q %>% is.na %>% which
			id_q[idx] <- -idx
		}
		if(anyNA(id_s)){
			idx <- id_s %>% is.na %>% which
			id_s[idx] <- -(nq + idx)
		}
		
	}else{
		#No join keys: everyone shares the same id => equivalent to no join constraint.
		id_q <- rep.int(1L, nq)
		id_s <- rep.int(1L, ns)
	}
	
	#Change seqnames to contain key
	q2 <- GRanges(
		seqnames = str_c(as.character(seqnames(query)), "|", id_q),
		ranges = ranges(query),
		strand = strand(query)
	)
	
	s2 <- GRanges(
		seqnames = str_c(as.character(seqnames(subject)), "|", id_s),
		ranges = ranges(subject),
		strand = strand(subject)
	) %>%
		GNCList
	
	hits <- findOverlaps(q2, s2, ignore.strand = ignore.strand) %>%
		suppressWarnings #Hide warning about unshared seqlevels between q2 and s2
	
	if(overlap_adjacent_query_insertion == TRUE){
		hits_adj <- findOverlaps(q2, s2, ignore.strand = ignore.strand, maxgap = 0) %>%
			suppressWarnings
		
		qh_adj <- queryHits(hits_adj)
		sh_adj <- subjectHits(hits_adj)
		
		hits_adj <- hits_adj[
			(query[qh_adj] %>% width == 0) & (subject[sh_adj] %>% width != 0)
			]
		
		hits <- GenomicRanges::union(hits, hits_adj)
		rm(hits_adj)
	}
	
	if(length(hits) == 0){return(rep(FALSE, nq))} #No overlap

	#Output result
	tabulate(queryHits(hits), nbins = nq) > 0L
}
