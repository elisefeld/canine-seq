prep_samples <- function(sample_sheet_path) {
  if (file.exists(sample_sheet_path)) {
    message("Reading sample sheet from ", sample_sheet_path, "...")
    sample_sheet <- read.csv(sample_sheet_path) |>
      dplyr::mutate(time_num = as.numeric(time),
                    time = factor(time, levels = sort(unique(time_num))),
                    subject = factor(subject, levels = sort(unique(subject))))
    rownames(sample_sheet) <- sample_sheet$sample_id
    message("Found ", nrow(sample_sheet), " samples in sample sheet.")
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
    if (file.exists(homolog_path)) {
      message("Reading csv from ", homolog_path, "...")
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
        summarise(across(everything(), ~paste0(unique(.), collapse = ", "))) |>
        dplyr::mutate(alt_human_gene_ids = str_extract(human_gene_id, "(.*)"),
               human_gene_id = str_extract(human_gene_id, "^[^,]+"),
               alt_human_gene_names = str_extract(human_gene_name, "(.*)"),
               human_gene_name = str_extract(human_gene_name, "^[^,]+"))
      message("Found ", length(unique(homologs$human_gene_id)), " homologous genes in humans.")
      
      sym_to_entrez <- bitr(unique(homologs$human_gene_name), fromType='SYMBOL', toType='ENTREZID', OrgDb="org.Hs.eg.db") |>
        group_by(SYMBOL) |>
        summarise(across(everything(), ~paste0(unique(.), collapse = ", "))) |>
        dplyr::mutate(alt_entrez_id = str_extract(ENTREZID, "(.*)"),
               entrez_id = str_extract(ENTREZID, "^[^,]+")) |>
        dplyr::select(SYMBOL, entrez_id, alt_entrez_id)
      
      homologs <- homologs |>
        left_join(sym_to_entrez, by = c("human_gene_name" = "SYMBOL"))
      
      tx2gene <- tx2gene |>
        left_join(homologs, by = c("gene_id" = "dog_gene_id")) |>
        dplyr::mutate(human_gene_name = coalesce(human_gene_name, human_gene_id))
      message(sum(is.na(tx2gene$human_gene_id)), " dog transcript IDs do not have a human equivalent.")
      return(tx2gene)
    }
    else {
      message("Homolog file ", homolog_path, " does not exist.")
      return(tx2gene)
    }
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
  message("Reading abundance.h5 files from ", kallisto_path, "...")
  files <- file.path(kallisto_path, sample_sheet$sample_id, "abundance.h5") 
  if (length(files[file.exists(files)]) > 0){
    message("Found ", length(files[file.exists(files)]), " files in ", kallisto_path)
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
    message("Reading csv from ", homolog_path, "...")
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
      summarise(across(everything(), ~paste0(unique(.), collapse = ", "))) |>
      dplyr::mutate(has_alt_human_genes = str_detect(human_gene_id, ","),
             alt_human_gene_ids = str_extract(human_gene_id, "(.*)"),
             human_gene_id = str_extract(human_gene_id, "^[^,]+"),
             alt_human_gene_names = str_extract(human_gene_name, "(.*)"),
             human_gene_name = str_extract(human_gene_name, "^[^,]+"))
    message("Found ", length(unique(homologs$human_gene_id)), " homologous genes in humans.")
    
    sym_to_entrez <- clusterProfiler::bitr(unique(homologs$human_gene_name), fromType='SYMBOL', toType='ENTREZID', OrgDb="org.Hs.eg.db") |>
      group_by(SYMBOL) |>
      summarise(across(everything(), ~paste0(unique(.), collapse = ", "))) |>
      dplyr::mutate(alt_entrez_id = str_extract(ENTREZID, "(.*)"),
             entrez_id = str_extract(ENTREZID, "^[^,]+")) |>
      dplyr::select(SYMBOL, entrez_id, alt_entrez_id)
    
    homologs <- homologs |>
      left_join(sym_to_entrez, by = c("human_gene_name" = "SYMBOL"))
    
    return(homologs)
  }
  else {
    stop("Homolog file ", homolog_path, " does not exist.")
  }
}

run_model <- function(txi,
                      sample_sheet,
                      design,
                      reduced,
                      test = "Wald",
                      min_count,
                      min_samples) {
  dds <- DESeqDataSetFromTximport(txi,
                                  colData = sample_sheet,
                                  design = design)
  message("Filtering out genes that don't have a count above ", min_count, " in ", min_samples, " or more samples")
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
                          #homologs,
                          min_alpha,
                          min_logfold,
                          shrink=TRUE) {
  # getting a data frame that has log fold change, p val and gene name for every coefficient
  coefficients <- resultsNames(dds)
  df_list <- list()
  for (coef in coefficients[-1]) {
      if (shrink == TRUE){
      message("Shrinking log2FoldChange values for coefficient: ", coef)
      res <- lfcShrink(dds, coef = coef, type = "apeglm")
      }
      else{
        res <- results(dds, name = coef)
      }
  
      res <- res |>
      as.data.frame() |>
      rownames_to_column("gene_id") |>
      left_join(tx2gene,
                by = "gene_id") |>
      dplyr::select(-transcript_id) |>
      dplyr::distinct(gene_id, .keep_all = TRUE) |>
      dplyr::filter(!is.na(padj)) |> # genes with extreme outliers are set to NA in DESeq2
      dplyr::arrange(padj) |>
      dplyr::mutate(coefficient = coef)
    df_list <- append(df_list, list(res))
  }
  df <- do.call(rbind, df_list) |>
    #left_join(homologs, by = c("gene_id" = "dog_gene_id")) |>
    dplyr::mutate(regulation = case_when(
      (padj < min_alpha & abs(log2FoldChange) > min_logfold & log2FoldChange > 0) ~ "upregulated",
      (padj < min_alpha & abs(log2FoldChange) > min_logfold & log2FoldChange < 0) ~ "downregulated",
      T ~ "nonDE"),
      coefficient = factor(coefficient, levels = resultsNames(dds))) 
  return(df)
}

run_gsea <- function(dds, main_df, path_data, path_terms){
  set.seed(42)
  df_list <- list()
  coefs <- resultsNames(dds)
  
  for (coef in coefs[-1]) {
    # create ranked list of genes
    ranked_genes <- main_df |>
      dplyr::filter(coefficient == coef & !is.na(human_gene_name)) |>
      arrange(desc(log2FoldChange)) |>
      pull(log2FoldChange, name = human_gene_id)
    
    df_gsea <- GSEA(ranked_genes,
                    TERM2GENE = path_data,
                    TERM2NAME = path_terms,
                    seed = TRUE) |>
      as.data.frame() |>
      dplyr::mutate(coefficient = coef)
    
    df_list <- append(df_list, list(df_gsea))
  }
  df_gsea <- do.call(rbind, df_list)
  return(df_gsea)
}