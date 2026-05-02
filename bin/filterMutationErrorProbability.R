#!/usr/bin/env -S Rscript --vanilla

cat("#### Running filterMutationErrorProbability ####\n")

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(qs2))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(tidyverse))

options(warn = 2)

option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = NULL, help = "path to YAML configuration file"),
  make_option(c("-s", "--sample_id_toanalyze"), type = "character", default = NULL, help = "sample_id to analyze"),
  make_option(c("-g", "--chromgroup_toanalyze"), type = "character", default = NULL, help = "chromgroup to analyze"),
  make_option(c("-v", "--filtergroup_toanalyze"), type = "character", default = NULL, help = "filtergroup to analyze"),
  make_option(c("-f", "--files"), type = "character", default = NULL, help = "comma-separated filterCalls qs2 files"),
  make_option(c("-o", "--output_prefix"), type = "character", default = NULL, help = "output prefix")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$config) || is.null(opt$sample_id_toanalyze) || is.null(opt$chromgroup_toanalyze) || is.null(opt$filtergroup_toanalyze) || is.null(opt$files) || is.null(opt$output_prefix)) {
  stop("Missing input parameter(s)!")
}

yaml.config <- suppressWarnings(read.config(opt$config))
filterCallsFiles <- opt$files %>% str_split_1(",") %>% str_trim()
sample_id_toanalyze <- opt$sample_id_toanalyze
chromgroup_toanalyze <- opt$chromgroup_toanalyze
filtergroup_toanalyze <- opt$filtergroup_toanalyze
output_prefix <- opt$output_prefix

filtergroup_config <- yaml.config$filtergroups %>%
  enframe(name = NULL) %>%
  unnest_wider(value) %>%
  filter(filtergroup == filtergroup_toanalyze)

if (nrow(filtergroup_config) != 1) stop("Could not find exactly one filtergroup configuration.")

max_mutation_errors_per_bp <- suppressWarnings(as.numeric(filtergroup_config$max_mutation_errors_per_bp))
if (is.na(max_mutation_errors_per_bp)) max_mutation_errors_per_bp <- NULL

default_qual_bin_edges <- c(10, 20, 30, 40, 50, 60, 70, 80, 89, Inf)
parse_qual_bin_edges <- function(x, strict = FALSE) {
  if (is.null(x) || (length(x) == 1 && is.na(x))) return(default_qual_bin_edges)
  if (is.character(x) && length(x) == 1) x <- str_split_1(x, ",") %>% str_trim()
  parsed <- map(x, ~ {
    v <- suppressWarnings(as.numeric(.x))
    if (!is.na(v)) return(list(value = v, valid = TRUE))
    if (tolower(as.character(.x)) %in% c("inf", "+inf", "infinity", "+infinity")) return(list(value = Inf, valid = TRUE))
    list(value = NA_real_, valid = FALSE)
  })
  vals <- parsed %>% map_dbl("value")
  invalid_entries <- parsed %>% map_lgl(~ !.x$valid)
  if (strict && any(invalid_entries)) {
    stop("Invalid value(s) in filterMutationErrorProbability_qual_bins. Use numeric upper edges and/or Inf.")
  }
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) {
    if (strict) stop("filterMutationErrorProbability_qual_bins could not be parsed.")
    return(default_qual_bin_edges)
  }
  vals <- unique(sort(vals))
  if (strict && any(vals <= 0)) {
    stop("filterMutationErrorProbability_qual_bins must contain positive upper edges.")
  }
  if (!is.infinite(tail(vals, 1))) vals <- c(vals, Inf)
  vals
}

make_qual_bin_labels <- function(edges) {
  lower <- 1
  map_chr(edges, ~ {
    if (is.infinite(.x)) {
      str_c(">=", lower)
    } else {
      label <- str_c(lower, "-", as.integer(.x))
      lower <<- as.integer(.x) + 1L
      label
    }
  })
}

qual_bins_raw <- filtergroup_config$filterMutationErrorProbability_qual_bins[[1]]
if (!is.null(max_mutation_errors_per_bp) && (is.null(qual_bins_raw) || (length(qual_bins_raw) == 1 && (is.na(qual_bins_raw) || qual_bins_raw == "")))) {
  stop("filterMutationErrorProbability_qual_bins must be defined when max_mutation_errors_per_bp is set.")
}
qual_bin_edges <- parse_qual_bin_edges(qual_bins_raw, strict = !is.null(max_mutation_errors_per_bp))
qual_bin_labels <- make_qual_bin_labels(qual_bin_edges)

hp_cache <- new.env(parent = emptyenv())

safe_count_bind <- function(parts, group_cols, value_col) {
  if (length(parts) == 0) {
    out <- as_tibble(setNames(replicate(length(group_cols), character(), simplify = FALSE), group_cols))
    out[[value_col]] <- numeric()
    return(out)
  }
  bind_rows(parts) %>%
    group_by(across(all_of(group_cols))) %>%
    summarize("{value_col}" := sum(n, na.rm = TRUE), .groups = "drop")
}

truncate_indel_len <- function(ref, alt) {
  len <- nchar(alt) - nchar(ref)
  pmax(-5L, pmin(5L, as.integer(len)))
}

safe_mean <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  vals <- suppressWarnings(as.numeric(unlist(x, use.names = FALSE)))
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) return(NA_real_)
  mean(vals)
}

qual_bin_label <- function(x) {
  if_else(
    is.na(x),
    NA_character_,
    as.character(cut(x, breaks = c(-Inf, qual_bin_edges), labels = qual_bin_labels, include.lowest = TRUE, right = TRUE))
  )
}

calc_hp_len <- function(seqname, pos1, bsgenome_obj) {
  key <- str_c(seqname, ":", pos1)
  if (exists(key, envir = hp_cache, inherits = FALSE)) {
    return(get(key, envir = hp_cache, inherits = FALSE))
  }

  chr_len <- seqlengths(bsgenome_obj)[[seqname]]
  if (is.na(chr_len) || pos1 < 1 || pos1 > chr_len) {
    assign(key, 1L, envir = hp_cache)
    return(1L)
  }

  window_start <- max(1L, pos1 - 50L)
  window_end <- min(chr_len, pos1 + 50L)
  seq_window <- as.character(getSeq(bsgenome_obj, seqname, start = window_start, end = window_end))
  center_idx <- pos1 - window_start + 1L
  base <- substr(seq_window, center_idx, center_idx)

  if (base %in% c("N", "")) {
    assign(key, 1L, envir = hp_cache)
    return(1L)
  }

  left <- center_idx
  while (left > 1L && substr(seq_window, left - 1L, left - 1L) == base) left <- left - 1L
  right <- center_idx
  while (right < nchar(seq_window) && substr(seq_window, right + 1L, right + 1L) == base) right <- right + 1L

  hp <- as.integer(min(40L, right - left + 1L))
  assign(key, hp, envir = hp_cache)
  hp
}

make_sbs_channel <- function(reftnc, altbase) {
  str_c(reftnc, ">", str_sub(reftnc, 1, 1), altbase, str_sub(reftnc, 3, 3))
}

sbs_channel_rc <- function(channel) {
  reftnc <- str_sub(channel, 1, 3)
  alttnc <- str_sub(channel, 5, 7)
  reftnc_rc <- reftnc %>% DNAStringSet() %>% reverseComplement() %>% as.character()
  alttnc_rc <- alttnc %>% DNAStringSet() %>% reverseComplement() %>% as.character()
  str_c(reftnc_rc, ">", alttnc_rc)
}

get_passfilter_columns <- function(df) {
  names(df) %>%
    str_subset("passfilter$") %>%
    setdiff(c("mutation_error_probability.passfilter", "mutation_error_source.passfilter"))
}

suppressPackageStartupMessages(library(yaml.config$BSgenome$BSgenome_name, character.only = TRUE, lib.loc = yaml.config$cache_dir))
bsgenome_obj <- get(yaml.config$BSgenome$BSgenome_name)

cat("## Pass 1/2: accumulating mutation error probability inputs...\n")

sbs_channel_parts <- list()
sbs_channel_bq_parts <- list()
sbs_reftnc_cov_parts <- list()
indel_ctx_parts <- list()
indel_ctx_bq_parts <- list()
indel_interrogated_parts <- list()

for (i in seq_along(filterCallsFiles)) {
  cat(" ", i, sep = "")
  x <- qs_read(filterCallsFiles[i])
  calls <- x$calls
  pass_cols <- get_passfilter_columns(calls)

  train <- calls %>%
    filter(
      call_class %in% c("SBS", "indel"),
      SBSindel_call_type == "mismatch-ss",
      if_all(all_of(pass_cols), ~ .x == TRUE)
    ) %>%
    mutate(
      mean_qual = map_dbl(qual, safe_mean),
      mean_qual_opp = map_dbl(qual.opposite_strand, safe_mean),
      mean_qual_both = pmin(mean_qual, mean_qual_opp, na.rm = TRUE),
      mean_qual_both = if_else(is.infinite(mean_qual_both), mean_qual, mean_qual_both),
      qual_bin = qual_bin_label(mean_qual_both)
    )

  sbs_train <- train %>%
    filter(call_class == "SBS") %>%
    mutate(channel = make_sbs_channel(reftnc_template_strand, alt_template_strand))

  if (nrow(sbs_train) > 0) {
    sbs_channel_parts[[length(sbs_channel_parts) + 1]] <- sbs_train %>% count(channel, name = "n")
    sbs_channel_bq_parts[[length(sbs_channel_bq_parts) + 1]] <- sbs_train %>% count(channel, qual_bin, name = "n")
  }

  sbs_cov_tbl <- x$bam.gr.filtertrack.bytype %>%
    filter(call_class == "SBS", SBSindel_call_type == "mismatch-ss") %>%
    pluck("bam.gr.filtertrack.reftnc_both_strands", 1)

  if (!is.null(sbs_cov_tbl) && nrow(sbs_cov_tbl) > 0) {
    sbs_reftnc_cov_parts[[length(sbs_reftnc_cov_parts) + 1]] <- sbs_cov_tbl %>%
      transmute(reftnc = as.character(reftnc), n = as.numeric(count))
  }

  indel_train <- train %>%
    filter(call_class == "indel") %>%
    mutate(
      indel_len = truncate_indel_len(ref_plus_strand, alt_plus_strand),
      hp = map2_int(as.character(seqnames), as.integer(start), ~ calc_hp_len(.x, .y, bsgenome_obj)),
      context = str_c(call_type, "|len", indel_len, "|hp", hp)
    )

  if (nrow(indel_train) > 0) {
    indel_ctx_parts[[length(indel_ctx_parts) + 1]] <- indel_train %>% count(context, name = "n")
    indel_ctx_bq_parts[[length(indel_ctx_bq_parts) + 1]] <- indel_train %>% count(context, qual_bin, name = "n")
  }

  indel_cov_tbl <- x$bam.gr.filtertrack.bytype %>%
    filter(call_class == "indel", SBSindel_call_type == "mismatch-ss")

  if (nrow(indel_cov_tbl) > 0) {
    indel_width <- indel_cov_tbl %>%
      pluck("bam.gr.filtertrack", 1) %>%
      width() %>%
      sum(na.rm = TRUE)
    indel_interrogated_parts[[length(indel_interrogated_parts) + 1]] <- tibble(n = as.numeric(indel_width))
  }

  rm(x, calls, train, sbs_train, indel_train)
  invisible(gc())
}
cat(" DONE\n")

sbs_channel_counts <- safe_count_bind(sbs_channel_parts, c("channel"), "n_channel")
sbs_channel_bq_counts <- safe_count_bind(sbs_channel_bq_parts, c("channel", "qual_bin"), "n_channel_bq")
sbs_reftnc_interrogated <- safe_count_bind(sbs_reftnc_cov_parts, c("reftnc"), "interrogated")
indel_ctx_counts <- safe_count_bind(indel_ctx_parts, c("context"), "n_context")
indel_ctx_bq_counts <- safe_count_bind(indel_ctx_bq_parts, c("context", "qual_bin"), "n_context_bq")
indel_interrogated_total <- if (length(indel_interrogated_parts) == 0) 0 else bind_rows(indel_interrogated_parts) %>% pull(n) %>% sum(na.rm = TRUE)

qc_prob_parts <- list()
qc_source_parts <- list()
qc_context_parts <- list()
qc_context_bin <- tibble()
sbs_spectra_by_qual_bin <- tibble()

make_output_chunk_file <- function(input_path) {
  base <- basename(input_path)
  out <- base %>% str_replace("\\.filterCalls\\.", ".filterMutationErrorProbability.")
  if (identical(out, base)) out <- base %>% str_replace("\\.qs2$", ".filterMutationErrorProbability.qs2")
  out
}

cat("## Pass 2/2: assigning mutation error probabilities and writing updated chunks...\n")

for (i in seq_along(filterCallsFiles)) {
  cat(" ", i, sep = "")
  in_file <- filterCallsFiles[i]
  out_file <- make_output_chunk_file(in_file)

  x <- qs_read(in_file)
  calls <- x$calls %>%
    mutate(
      row_id = row_number(),
      mean_qual = map_dbl(qual, safe_mean),
      mean_qual_opp = map_dbl(qual.opposite_strand, safe_mean),
      mean_qual_both = pmin(mean_qual, mean_qual_opp, na.rm = TRUE),
      mean_qual_both = if_else(is.infinite(mean_qual_both), mean_qual, mean_qual_both),
      qual_bin = qual_bin_label(mean_qual_both)
    )

  eligible <- calls %>%
    filter(call_class %in% c("SBS", "indel"), SBSindel_call_type == "mutation")

  if (is.null(max_mutation_errors_per_bp)) {
    updates <- eligible %>%
      transmute(
        row_id,
        mutation_error_probability_source = "disabled",
        mutation_error_source.passfilter = TRUE,
        call_class,
        context = NA_character_
      )
  } else {
    sbs_rate_table <- sbs_channel_bq_counts %>%
      left_join(sbs_channel_counts, by = "channel") %>%
      group_by(channel) %>%
      mutate(
        n_channel = replace_na(n_channel, 0),
        errors_per_bp = if_else(n_channel > 0, n_channel_bq / n_channel, NA_real_),
        bin_rank = length(qual_bin_labels) - match(qual_bin, qual_bin_labels) + 1,
        context_filtered = case_when(
          is.na(errors_per_bp) ~ FALSE,
          any(errors_per_bp > max_mutation_errors_per_bp, na.rm = TRUE) ~ bin_rank >= min(bin_rank[errors_per_bp > max_mutation_errors_per_bp], na.rm = TRUE),
          TRUE ~ FALSE
        )
      ) %>%
      ungroup() %>%
      select(channel, qual_bin, errors_per_bp, context_filtered)

    indel_rate_table <- indel_ctx_bq_counts %>%
      left_join(indel_ctx_counts, by = "context") %>%
      group_by(context) %>%
      mutate(
        n_context = replace_na(n_context, 0),
        errors_per_bp = if_else(n_context > 0, n_context_bq / n_context, NA_real_),
        bin_rank = length(qual_bin_labels) - match(qual_bin, qual_bin_labels) + 1,
        context_filtered = case_when(
          is.na(errors_per_bp) ~ FALSE,
          any(errors_per_bp > max_mutation_errors_per_bp, na.rm = TRUE) ~ bin_rank >= min(bin_rank[errors_per_bp > max_mutation_errors_per_bp], na.rm = TRUE),
          TRUE ~ FALSE
        )
      ) %>%
      ungroup() %>%
      select(context, qual_bin, errors_per_bp, context_filtered)

    qc_context_bin <- bind_rows(
      sbs_rate_table %>% transmute(call_class = "SBS", context = channel, qual_bin, errors_per_bp, pass = !context_filtered),
      indel_rate_table %>% transmute(call_class = "indel", context, qual_bin, errors_per_bp, pass = !context_filtered)
    )
    sbs_spectra_by_qual_bin <- sbs_channel_bq_counts %>%
      group_by(qual_bin) %>%
      mutate(total_in_bin = sum(n_channel_bq, na.rm = TRUE), frac_in_bin = if_else(total_in_bin > 0, n_channel_bq / total_in_bin, NA_real_)) %>%
      ungroup() %>%
      transmute(qual_bin, channel, n = n_channel_bq, fraction = frac_in_bin)

    sbs_updates <- eligible %>%
      filter(call_class == "SBS") %>%
      mutate(
        channel = make_sbs_channel(reftnc_template_strand, alt_template_strand),
        context = channel
      ) %>%
      left_join(sbs_rate_table %>% select(channel, qual_bin, errors_per_bp, context_filtered), by = c("channel", "qual_bin")) %>%
      mutate(
        context_filtered = replace_na(context_filtered, FALSE),
        mutation_error_probability_source = "sbs_context_bq_cutoff",
        mutation_error_source.passfilter = !context_filtered
      ) %>%
      select(row_id, mutation_error_probability_source, mutation_error_source.passfilter, call_class, context)

    indel_updates <- eligible %>%
      filter(call_class == "indel") %>%
      mutate(
        indel_len = truncate_indel_len(ref_plus_strand, alt_plus_strand),
        hp = map2_int(as.character(seqnames), as.integer(start), ~ calc_hp_len(.x, .y, bsgenome_obj)),
        context = str_c(call_type, "|len", indel_len, "|hp", hp)
      ) %>%
      left_join(indel_rate_table %>% select(context, qual_bin, errors_per_bp, context_filtered), by = c("context", "qual_bin")) %>%
      mutate(
        context_filtered = replace_na(context_filtered, FALSE),
        mutation_error_probability_source = "indel_context_bq_cutoff",
        mutation_error_source.passfilter = !context_filtered
      ) %>%
      select(row_id, mutation_error_probability_source, mutation_error_source.passfilter, call_class, context)

    updates <- bind_rows(sbs_updates, indel_updates)
  }

  calls <- calls %>%
    mutate(
      mutation_error_probability_source = "not_applicable",
      mutation_error_source.passfilter = TRUE
    ) %>%
    left_join(
      updates %>% select(row_id, mutation_error_probability_source, mutation_error_source.passfilter),
      by = "row_id",
      suffix = c("", ".new")
    ) %>%
    mutate(
      mutation_error_probability_source = coalesce(mutation_error_probability_source.new, mutation_error_probability_source),
      mutation_error_source.passfilter = coalesce(mutation_error_source.passfilter.new, mutation_error_source.passfilter)
    ) %>%
    select(-ends_with(".new"))

  if (nrow(updates) > 0) {
    qc_source_parts[[length(qc_source_parts) + 1]] <- updates %>% count(mutation_error_probability_source, name = "n")

    qc_context_parts[[length(qc_context_parts) + 1]] <- updates %>%
      filter(!is.na(context)) %>%
      group_by(call_class, context) %>%
      summarize(
        n_total = n(),
        n_pass = sum(mutation_error_source.passfilter, na.rm = TRUE),
        mean_mutation_error_probability = NA_real_,
        .groups = "drop"
      )
  }

  x$calls <- calls %>% select(-row_id, -mean_qual, -mean_qual_opp, -mean_qual_both, -qual_bin)
  qs_save(x, out_file)

  rm(x, calls, eligible, updates)
  invisible(gc())
}
cat(" DONE\n")

cat("## Writing filterMutationErrorProbability QC outputs...\n")

qc_prefix <- str_c(output_prefix, ".filterMutationErrorProbability.qc")

qc_source <- bind_rows(qc_source_parts) %>%
  {
    if (nrow(.) == 0) tibble(source = character(), n = numeric())
    else . %>% group_by(mutation_error_probability_source) %>% summarize(n = sum(n), .groups = "drop") %>% rename(source = mutation_error_probability_source)
  }

qc_context <- bind_rows(qc_context_parts) %>%
  {
    if (nrow(.) == 0) {
      tibble(call_class = character(), context = character(), n_total = numeric(), n_pass = numeric(), mean_mutation_error_probability = numeric())
    } else {
      . %>%
        group_by(call_class, context) %>%
        summarize(
          n_total = sum(n_total),
          n_pass = sum(n_pass),
          mean_mutation_error_probability = weighted.mean(mean_mutation_error_probability, pmax(n_total, 1), na.rm = TRUE),
          .groups = "drop"
        )
    }
  }

summary_tbl <- tibble(
  sample_id = sample_id_toanalyze,
  chromgroup = chromgroup_toanalyze,
  filtergroup = filtergroup_toanalyze,
  max_mutation_errors_per_bp = if_else(is.null(max_mutation_errors_per_bp), NA_real_, max_mutation_errors_per_bp),
  filterMutationErrorProbability_qual_bins = str_c(qual_bin_labels, collapse = ","),
  sbs_training_channels = nrow(sbs_channel_counts),
  indel_training_contexts = nrow(indel_ctx_counts),
  indel_interrogated_total = indel_interrogated_total
)

write_tsv(summary_tbl, str_c(qc_prefix, ".summary.tsv"))
write_tsv(qc_source %>% arrange(desc(n)), str_c(qc_prefix, ".source.tsv"))
write_tsv(qc_context %>% arrange(call_class, context), str_c(qc_prefix, ".context.tsv"))
write_tsv(qc_context_bin %>% arrange(call_class, context, qual_bin), str_c(qc_prefix, ".context_bin_errors_per_bp.tsv"))
write_tsv(sbs_spectra_by_qual_bin %>% arrange(qual_bin, channel), str_c(qc_prefix, ".sbs_spectra_by_qual_bin.tsv"))

if (nrow(qc_context) > 0) {
  p_context <- ggplot(qc_context, aes(x = reorder(context, mean_mutation_error_probability), y = mean_mutation_error_probability, fill = call_class)) +
    geom_col() +
    coord_flip() +
    labs(x = "context", y = "mean mutation error probability", title = "Mean mutation error probability by context") +
    theme_bw()
  ggsave(str_c(qc_prefix, ".context_mean_probability.png"), p_context, width = 12, height = 10, dpi = 200)
}

if (nrow(qc_context_bin) > 0) {
  p_context_bin <- ggplot(qc_context_bin, aes(x = qual_bin, y = errors_per_bp, color = pass, group = context)) +
    geom_line(alpha = 0.5) +
    geom_point(size = 1) +
    facet_wrap(~ call_class, scales = "free_y") +
    scale_color_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#d95f02")) +
    labs(x = "quality bin", y = "errors per bp", color = "pass", title = "Context-level errors per bp by quality bin") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
  ggsave(str_c(qc_prefix, ".context_bin_errors_per_bp.pdf"), p_context_bin, width = 12, height = 8)
}

if (nrow(sbs_spectra_by_qual_bin) > 0) {
  p_sbs_spectra <- ggplot(sbs_spectra_by_qual_bin, aes(x = qual_bin, y = channel, fill = fraction)) +
    geom_tile() +
    scale_fill_viridis_c(option = "C", na.value = "grey90") +
    labs(x = "quality bin", y = "SBS channel", fill = "fraction", title = "ssDNA SBS spectra by quality bin") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 5))
  ggsave(str_c(qc_prefix, ".sbs_spectra_by_qual_bin.pdf"), p_sbs_spectra, width = 12, height = 10)
}

cat("DONE\n")
