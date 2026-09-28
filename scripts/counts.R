count_matrix_to_df <- function(count_matrix, sample_sheet) {
  count_df <- count_matrix |>
    as.data.frame() |>
    rownames_to_column("gene_id") |>
    pivot_longer(cols = -1,
                 names_to = "sample_id",
                 values_to = "count") |>
    left_join(sample_sheet, by = "sample_id")
  return(count_df)
}

merge_counts <- function(count_matrix, dds, sample_sheet, annotation = NULL) {
  raw_df <- count_matrix_to_df(count_matrix, sample_sheet)
  
  norm_df <- count_matrix_to_df(counts(dds, normalized = TRUE), sample_sheet) |>
    dplyr::mutate(count_norm = count) |>
    dplyr::select(-count)
    
  count_df <- raw_df |>
    left_join(norm_df) |>
    dplyr::rename(no = count,
                  yes = count_norm) |>
    pivot_longer(cols = c("yes", "no"),
                 names_to = "normalized",
                 values_to = "count")
  
    if (!is.null(annotation)) {
    count_df <- count_df |>
      left_join(annotation, by = "gene_id")
    }
  return(count_df)
}

plot_gene_expression <- function(counts, gene) {
  counts |>
    dplyr::filter(gene_id == gene,
                  normalized == "yes") |>
    ggplot(aes(x = time, y = count, color = subject)) +
    geom_point() +
    geom_smooth(aes(group = subject), method = "loess", se = FALSE) +
    theme_classic() +
    labs(x = "Time (days)", 
         y = "Normalized count",
         title = gene)
}
