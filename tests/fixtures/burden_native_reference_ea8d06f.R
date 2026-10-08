# Independent original numerical reference from ea8d06f3c0f29a7145e9e56dee1b0bb2424f4b77.
# Exact eight function bodies from bin/calculateBurdens.R, SHA256
# 47bac8c08ba60533f11ac2b79e488376c57bc2faf87e14574b12899615d75211.
# Test fixture only; never sourced by the production pipeline.

sum_RleList <- function(a, b) {
	seqs_union <- union(names(a), names(b))
	
	seqs_union %>%
		map(function(nm){
			seq_in_a <- nm %in% names(a)
			seq_in_b <- nm %in% names(b)
			if(seq_in_a && seq_in_b){
				a[[nm]] + b[[nm]]
			}else if(seq_in_a){
				a[[nm]]
			}else{
				b[[nm]]
			}
		}) %>%
			set_names(seqs_union) %>%
			RleList(compress=FALSE)
}

bc_orientation_is_asymmetric <- function(bc_orientation){
	bc_orientation %>%
		as.character %>%
		str_split("-") %>%
		map_lgl(function(x){length(x) == 2 && x[1] != x[2]})
}

validate_bam.gr.filtertrack <- function(gr){
	x_plus <- gr %>% filter(strand == "+") %>% select(-bc_orientation)
	x_minus <- gr %>% filter(strand == "-") %>% select(-bc_orientation)
	if(!(identical(ranges(x_plus), ranges(x_minus)) & identical(mcols(x_plus), mcols(x_minus)))){
		stop("Mismatched plus and minus strand ranges in bam.gr.filtertrack!")
	}
	invisible(TRUE)
}

calc_duplex_coverage <- function(gr){
	cov <- gr %>% coverage
	if((cov %% 2L != 0L) %>% any %>% any){
		stop("Non-even strand coverage in bam.gr.filtertrack!")
	}
	cov %/% 2L
}

accumulate_bam.gr.filtertracks <- function(state, incoming, orientation_call_types=NULL){
	#Helper function to calculate strand-level coverage for each barcode configuration and aligned read strand.
	calc_by_bc_orientation_strand_coverage <- function(gr, needed){
		#Keep the same column types for empty and populated orientation tables,
		#without converting genomic coordinates and unrelated metadata to a tibble.
		if(!needed){gr <- gr[0]}
		orientations <- tibble(
			bc_orientation = gr$bc_orientation,
			strand = strand(gr) %>% as.factor
		)
		if(!needed){
			return(orientations[0,] %>% mutate(bam.gr.filtertrack.coverage = list()))
		}
		orientations %>%
			distinct(bc_orientation, strand) %>%
			filter(bc_orientation_is_asymmetric(bc_orientation)) %>%
			mutate(
				bam.gr.filtertrack.coverage = map2(
					bc_orientation,
					strand,
					function(x,y){
						gr %>%
							filter(bc_orientation == x, strand == y) %>%
							coverage
					}
				)
			)
	}

	#Helper function to sum strand-level coverage across analysis chunks.
	sum_by_bc_orientation_strand_coverage <- function(a, b){
		full_join(
			a,
			b,
			by = join_by(bc_orientation, strand),
			suffix = c("", ".2")
		) %>%
			mutate(
				bam.gr.filtertrack.coverage = map2(
					bam.gr.filtertrack.coverage,
					bam.gr.filtertrack.coverage.2,
					function(x,y){
						if(is.null(x)){
							y
						}else if(is.null(y)){
							x
						}else{
							sum_RleList(x,y)
						}
					}
				)
			) %>%
			select(-bam.gr.filtertrack.coverage.2)
	}
	
	#Only this environment owns the large coverage lists while chunks are loaded.
	#Do not bind either list to another local variable: that would keep every old
	#category alive while replacements are constructed. Joins use metadata only.
	walk(incoming$bam.gr.filtertrack, validate_bam.gr.filtertrack)
	incoming_metadata <- incoming %>% select(-bam.gr.filtertrack)
	orientation_rows <- if(is.null(orientation_call_types)){
		integer()
	}else{
		incoming_metadata %>%
			mutate(.coverage_row = row_number()) %>%
			filter(SBSindel_call_type != "mutation") %>%
			semi_join(orientation_call_types, by = join_by(call_type, call_class, SBSindel_call_type, filtergroup)) %>%
			pull(.coverage_row)
	}
	first_chunk <- is.null(state$metadata)
	if(first_chunk){
		state$metadata <- incoming_metadata
		state$coverage <- vector("list", nrow(incoming_metadata))
		state$orientation <- vector("list", nrow(incoming_metadata))
		incoming_rows <- seq_len(nrow(incoming_metadata))
	}else{
		#Reproduce the former left join's order and metadata/factor promotion,
		#including repeated keys, without putting coverage into a data mask.
		matched <- left_join(
			state$metadata %>% mutate(.previous_row = row_number()),
			incoming_metadata %>% mutate(.incoming_row = row_number()),
			by = names(state$metadata)
		)
		if(!identical(matched$.previous_row, seq_along(state$coverage))){
			state$coverage <- state$coverage[matched$.previous_row]
			state$orientation <- state$orientation[matched$.previous_row]
		}
		state$metadata <- matched %>% select(-.previous_row, -.incoming_row)
		incoming_rows <- matched$.incoming_row
		#The legacy code checked every incoming category, even an unmatched one.
		for(j in setdiff(seq_len(nrow(incoming_metadata)), incoming_rows)){
			invisible(calc_duplex_coverage(incoming$bam.gr.filtertrack[[j]]))
		}
	}
	for(i in seq_along(incoming_rows)){
		j <- incoming_rows[i]
		cov <- if(is.na(j)) NULL else calc_duplex_coverage(incoming$bam.gr.filtertrack[[j]])
		if(first_chunk){
			state$coverage[[i]] <- cov
		}else{
			state$coverage[[i]] <- sum_RleList(state$coverage[[i]], cov)
		}
		rm(cov)
		orientation <- if(is.na(j)) NULL else calc_by_bc_orientation_strand_coverage(
			incoming$bam.gr.filtertrack[[j]], j %in% orientation_rows
		)
		if(first_chunk){
			state$orientation[[i]] <- orientation
		}else{
			state$orientation[[i]] <- sum_by_bc_orientation_strand_coverage(state$orientation[[i]], orientation)
		}
		rm(orientation)
	}
	invisible(NULL)
}

filtertrack_coverage_result <- function(state){
	state$metadata %>%
		mutate(
			bam.gr.filtertrack.coverage = state$coverage,
			bam.gr.filtertrack.by_bc_orientation_strand.coverage = state$orientation
		)
}

gr_1bp_cov <- function(gr, cov){
	
	stopifnot(all(width(gr) == 1))
	
	n <- length(gr)
	result <- rep.int(0, n) #pre-fill zeros
	chr <- gr %>% seqnames %>% as.character
	pos <- gr %>% start
	
	# indices per chromosome (keeps order in 'gr')
	idx_by_chr <- split(seq_len(n), chr)
	
	for (ch in names(idx_by_chr)) {
		idx <- idx_by_chr[[ch]]
		if (ch %in% names(cov)){
			result[idx] <- cov[[ch]][pos[idx]]
		}
	}
	
	return(result)
}

sum_filtertrack_sensitivity_coverage <- function(bytype, queries, previous = NULL){
	walk(bytype$bam.gr.filtertrack, validate_bam.gr.filtertrack)
	current <- bytype %>% select(-bam.gr.filtertrack)
	current$coverage_start <- vector("list", nrow(current))
	current$coverage_end <- vector("list", nrow(current))
	for(i in seq_len(nrow(current))){
		cov <- calc_duplex_coverage(bytype$bam.gr.filtertrack[[i]])
		query_index <- match(current$call_type[i], queries$call_type)
		if(is.na(query_index)){
			current$coverage_start[[i]] <- numeric()
			current$coverage_end[[i]] <- numeric()
		}else{
			current$coverage_start[[i]] <- gr_1bp_cov(queries$query_start[[query_index]], cov)
			current$coverage_end[[i]] <- gr_1bp_cov(queries$query_end[[query_index]], cov)
		}
		rm(cov)
	}
	if(is.null(previous)){return(current)}
	previous %>%
		left_join(current, by = setdiff(names(current), c("coverage_start", "coverage_end")), suffix = c("", ".2")) %>%
		mutate(
			coverage_start = map2(coverage_start, coverage_start.2, `+`),
			coverage_end = map2(coverage_end, coverage_end.2, `+`)
		) %>%
		select(-coverage_start.2, -coverage_end.2)
}
