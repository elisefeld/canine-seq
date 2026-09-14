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
                      
                      
                      
extract_results <- function(dds) {
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
    pivot_longer(cols = c("log2FoldChange", "log2FoldChange_shrunk",
                          names_to = shrunk,
                          values_to = log2FC))
  return(df)
}

extract_significance <- function(df,
                                 min_alpha,
                                 min_logfold,
                                 filter_by = c("shrunk", "unshrunk")) {
  filter_by <- match.arg(filter_by)
  filter_col <- if (filter_by == "shrunk") {
    "log2FoldChange_shrunk"
  } else {
    "log2FoldChange"
  }
  
df <- df |>
  dplyr::mutate(regulation = case_when(padj < min_alpha &
                                         abs(.data[[filter_col]]) > min_logfold &
                                         .data[[filter_col]] > 0 ~ "upregulated",
                                       
                                       padj < min_alpha &
                                         abs(.data[[filter_col]]) > min_logfold &
                                         .data[[filter_col]] < 0 ~ "downregulated",
                                       
                                       TRUE ~ "nonDE"),
                sig = padj < min_alpha & abs(.data[[filter_col]]) > min_logfold)
  
}