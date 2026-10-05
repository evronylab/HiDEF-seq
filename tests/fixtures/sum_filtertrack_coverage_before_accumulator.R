#Benchmark/reference fixture: orientation-gated whole-table accumulation before
#the per-category accumulator. Shared helper dependencies come from calculateBurdens.R.
sum_bam.gr.filtertracks <- function(bam.gr.filtertrack1, bam.gr.filtertrack2=NULL, orientation_call_types=NULL){

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
	
	#Calculate each incoming chunk once. Aggregate coverage and its invariant
	#checks remain required for every row, including rows used only internally.
	incoming <- if(is.null(bam.gr.filtertrack2)) bam.gr.filtertrack1 else bam.gr.filtertrack2
	walk(incoming$bam.gr.filtertrack, validate_bam.gr.filtertrack)
	orientation_rows <- if(is.null(orientation_call_types)){
		integer()
	}else{
		incoming %>%
			mutate(.coverage_row = row_number()) %>%
			filter(SBSindel_call_type != "mutation") %>%
			semi_join(
				orientation_call_types,
				by = join_by(call_type, call_class, SBSindel_call_type, filtergroup)
			) %>%
			pull(.coverage_row)
	}
	incoming <- incoming %>%
		mutate(
			bam.gr.filtertrack.coverage = bam.gr.filtertrack %>% map(calc_duplex_coverage),
			bam.gr.filtertrack.by_bc_orientation_strand.coverage = map2(
				bam.gr.filtertrack,
				seq_along(bam.gr.filtertrack) %in% orientation_rows,
				calc_by_bc_orientation_strand_coverage
			)
		) %>%
		select(-bam.gr.filtertrack)

	if(is.null(bam.gr.filtertrack2)){
		incoming
	}else{
		bam.gr.filtertrack1 %>%
			left_join(
				incoming,
				by = names(.) %>% setdiff(c("bam.gr.filtertrack.coverage", "bam.gr.filtertrack.by_bc_orientation_strand.coverage")),
				suffix = c("",".2")
			) %>%
			mutate(
				bam.gr.filtertrack.coverage = map2(
					bam.gr.filtertrack.coverage,
					bam.gr.filtertrack.coverage.2,
					function(x,y){sum_RleList(x,y)}
				),
				bam.gr.filtertrack.by_bc_orientation_strand.coverage = map2(
					bam.gr.filtertrack.by_bc_orientation_strand.coverage,
					bam.gr.filtertrack.by_bc_orientation_strand.coverage.2,
					function(x,y){sum_by_bc_orientation_strand_coverage(x,y)}
				)
			) %>%
			select(-bam.gr.filtertrack.coverage.2, -bam.gr.filtertrack.by_bc_orientation_strand.coverage.2)
	}
}

