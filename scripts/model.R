create_dds <- function(data, sample_sheet, design, source = c("tximport", "matrix")) {
    source <- match.arg(source)
    
    if (source == "tximport") {
      dds <- DESeqDataSetFromTximport(data,
                                      colData = sample_sheet,
                                      design = design)
    }
    else {
      dds <- DESeqDataSetFromMatrix(data,
                                    colData = sample_sheet,
                                    design = design)
    }
    return(dds)
}

run_model <- function(dds, min_count, min_samples, test = c("Wald", "LRT"), reduced = NULL){
  test <- match.arg(test)
  message("Filtering out genes that don't have a count above ", min_count, " in ", min_samples, " or more samples")
  keep <- rowSums(counts(dds) >= min_count) >= min_samples # filtering for low counts
  dds <- dds[keep, ]
  
  if (test == "LRT") {
    dds <- DESeq(dds, test = test, reduced = reduced)
  }
  
  else {
    dds <- DESeq(dds, test = test)
  }
  return(dds)
}
                      
                      
                      
extract_results <- function(dds, min_alpha, min_logfold) {
  coefficients <- resultsNames(dds)
  df_list <- list()
  for (coef in coefficients[-1]) {
    message("Processing coefficient: ", coef)
    
    # Unshrunk results
    res_unshrunk <- results(dds, name = coef) |>
      as.data.frame() |>
      rownames_to_column("gene_id") |>
      select(gene_id, log2FoldChange, lfcSE)
    
    # Shrunk results
    res_shrunk <- lfcShrink(dds, coef = coef, type = "apeglm") |>
      as.data.frame() |>
      rownames_to_column("gene_id") |>
      dplyr::mutate(log2FoldChange_shrunk = log2FoldChange,
                    lfcSE_shrunk = lfcSE) |>
      dplyr::select(-log2FoldChange, -lfcSE)
    
    # Combine
    res <- res_unshrunk |>
      dplyr::left_join(res_shrunk, by = "gene_id") |>
      dplyr::filter(!is.na(padj)) |>
      dplyr::arrange(padj) |>
      dplyr::mutate(coefficient = coef)
    
    df_list[[coef]] <- res
  }
  
  df <- do.call(rbind, df_list) |>
    dplyr::mutate(coefficient = factor(coefficient, levels = resultsNames(dds))) |>
    dplyr::rename(no = log2FoldChange,
           yes = log2FoldChange_shrunk) |>
    pivot_longer(cols = c("no", "yes"),
                          names_to = "shrunk",
                          values_to = "log2FoldChange") |>
    dplyr::mutate(regulation = case_when(padj < min_alpha &
                                           abs(log2FoldChange) > min_logfold &
                                           log2FoldChange > 0 ~ "upregulated",
                                         
                                         padj < min_alpha &
                                           abs(log2FoldChange) > min_logfold &
                                           log2FoldChange < 0 ~ "downregulated",
                                         
                                         TRUE ~ "nonDE"),
                  sig = padj < min_alpha & abs(log2FoldChange) > min_logfold)
  return(df)
}