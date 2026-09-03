prep_samples <- function(sample_sheet_path) {
  if (file.exists(sample_sheet_path)) {
    message("Reading sample sheet from ", sample_sheet_path, "...")
    sample_sheet <- read.csv(sample_sheet_path) |>
      dplyr::mutate(time = factor(time, levels = sort(unique(time))),
                    time_lab = factor(paste0("Day ", time)),
                    subject = factor(subject, levels = sort(unique(subject))),
                    subject_lab = factor(paste0("Subject ", subject)),
                    sample_name = paste0("Subject ", subject, ", ", "Day ", time)) # readable sample name for visualizations
    rownames(sample_sheet) <- sample_sheet$sample_id
    print(paste0("Found ", nrow(sample_sheet), " samples in sample sheet."))
    return(sample_sheet)
  }
  else {
    stop("Sample sheet ", sample_sheet_path, " does not exist.")
  }
}

prep_tx2gene <- function(gtf_path) {
  if (file.exists(gtf_path)) {
    message("Reading annotation file from ", gtf_path, "...")
    tx2gene <- rtracklayer::readGFF(gtf_path) |>
      as.data.frame() |>
      dplyr::select(transcript_id,
                    gene_id,
                    gene_name,
                    gene_biotype) |>
      dplyr::filter(!is.na(transcript_id)) |>
      dplyr::distinct() |>
      dplyr::mutate(gene_name = coalesce(gene_name,
                                         gene_id)) # if no gene name present, use gene id
    message("Found ", length(unique(tx2gene$transcript_id)), " unique transcripts and ", length(unique(tx2gene$gene_id)), " unique genes.")
    return(tx2gene)
  }
  else {
    stop("GTF file ", gtf_path, " does not exist.")
  }
}

prep_kallisto <- function(kallisto_path,
                          tx2gene,
                          sample_sheet,
                          transcripts = FALSE) {
  # NOTE: this will not include samples in the file path if they are not in the sample sheet
  print(paste0("Reading abundance.h5 files from ", kallisto_path, "..."))
  files <- file.path(kallisto_path, sample_sheet$sample_id, "abundance.h5") 
  if (length(files) > 0){
    print(paste0("Found ", length(files), " files in ", kallisto_path))
    names(files) <- sample_sheet$sample_id
    txi <- tximport(files,
                    type = "kallisto",
                    tx2gene = tx2gene,
                    txOut = transcripts, # change this to true if you want transcript level data
                    ignoreTxVersion = TRUE)
    return(txi)
  }
  else {
    stop("No kallisto abundance files found at ", kallisto_path, ".")
  }
}

prep_homologs <- function(homolog_path) {
  if (file.exists(homolog_path)) {
    print(paste0("Reading csv from ", homolog_path, "..."))
    homologs <- read.csv(homolog_path, sep = "\t") |>
      dplyr::filter(ref_species == "homo_sapiens") |>
      dplyr::select(ref_gene_stable_id,
                    ref_gene_name,
                    query_gene_stable_id,
                    query_gene_name) |>
      dplyr::rename(human_gene_id = ref_gene_stable_id,
                    human_gene_name = ref_gene_name,
                    dog_gene_id = query_gene_stable_id,
                    dog_gene_name = query_gene_name) |>
      group_by(dog_gene_id) |>
      summarise(across(everything(), ~paste0(unique(.), collapse = ", ")))
    print(paste0("Found ", length(unique(homologs$human_gene_id)), " homologous genes in humans."))
    return(homologs)
  }
  else {
    stop("Homolog file ", homolog_path, " does not exist.")
  }
}

run_model <- function(txi = txi,
                      sample_sheet,
                      design,
                      reduced,
                      test = "Wald",
                      min_count,
                      min_samples) {
  dds <- DESeqDataSetFromTximport(txi,
                                  colData = sample_sheet,
                                  design = design)
  keep <- rowSums(counts(dds) >= min_count) >= min_samples # filtering for low counts
  dds <- dds[keep, ]
  
  if (test == "LRT") {
    dds <- DESeq(dds,
                 test = test,
                 reduced = reduced)
  } else {
    dds <- DESeq(dds,
                 test = test)
  }
  return(dds)
}

process_counts <- function(counts,
                           sample_sheet,
                           tx2gene) {
  genes <- tx2gene |>
    dplyr::select(-transcript_id) |>
    dplyr::distinct()
  counts <- counts |>
    as.data.frame() |>
    rownames_to_column("gene_id") |>
    pivot_longer(cols = starts_with("OR"),
                 names_to = "sample_id",
                 values_to = "count") |>
    left_join(sample_sheet, by = "sample_id") |>
    left_join(genes, by = "gene_id") |>
    dplyr::distinct(gene_id, sample_id, .keep_all = TRUE)
  return(counts)
}

extract_coefs <- function(dds,
                          tx2gene,
                          homologs,
                          min_alpha,
                          min_logfold) {
  # getting a data frame that has log fold change, p val and gene name for every coefficient
  coefficients <- resultsNames(dds)
  df_list <- list()
  for (coef in coefficients[-1]) {
    res <- lfcShrink(dds, coef = coef, type = "apeglm") |>
      as.data.frame() |>
      rownames_to_column("gene_id") |>
      left_join(tx2gene,
                by = "gene_id") |>
      dplyr::select(-transcript_id) |>
      dplyr::distinct(gene_id, .keep_all = TRUE) |>
      dplyr::filter(!is.na(padj) & padj < min_alpha & abs(log2FoldChange) > min_logfold) |> # genes with extreme outliers are set to NA in DESeq2
      dplyr::arrange(padj) |>
      dplyr::mutate(coefficient = coef)
    df_list <- append(df_list, list(res))
  }
  df <- do.call(rbind, df_list) |>
    left_join(homologs, by = c("gene_id" = "dog_gene_id")) |>
    mutate(regulation = case_when(
      (padj < min_alpha & abs(log2FoldChange) > min_logfold & log2FoldChange > 0) ~ "upregulated",
      (padj < min_alpha & abs(log2FoldChange) > min_logfold & log2FoldChange < 0) ~ "downregulated",
      T ~ "nonDE"))
  return(df)
}

run_gsea <- function(dds, main_df, path_data, path_terms){
  set.seed(42)
  df_list <- list()
  coefs <- resultsNames(dds)
  
  for (coef in coefs[-1]) {
    # create ranked list of genes
    ranked_genes <- main_df |>
      filter(coefficient == coef & !is.na(human_gene_name)) |>
      arrange(desc(log2FoldChange)) |>
      pull(log2FoldChange, name = human_gene_id)
    
    df_gsea <- GSEA(ranked_genes,
                    TERM2GENE = path_data,
                    TERM2NAME = path_terms,
                    seed = TRUE) |>
      as.data.frame() |>
      mutate(coefficient = coef)
    
    df_list <- append(df_list, list(df_gsea))
  }
  df_gsea <- do.call(rbind, df_list)
  return(df_gsea)
}