create_ranked_list <- function(df, coef){
  ranked_genes <- df |>
    dplyr::filter(coefficient == coef & !is.na(gene_id)) |>
    arrange(desc(log2FoldChange)) |>
    pull(log2FoldChange, name = gene_id)
}

run_gsea <- function(dds, df, collection, subcollection = NULL, match_by){
  set.seed(42)
  df_list <- list()
  coefs <- resultsNames(dds)
  
  path_terms <- msigdbr(species = "Homo sapiens", collection = collection, subcollection = subcollection) |>
    dplyr::select(gs_name, gs_description)
  
  path_data <- msigdbr(species = "Homo sapiens", collection = collection, subcollection = subcollection) |>
    dplyr::select(gs_name, all_of(match_by))
  
  for (coef in coefs[-1]) {
    # create ranked list of genes
    ranked_genes <- create_ranked_list(df, coef)
    
    df_gsea <- GSEA(ranked_genes,
                    TERM2GENE = path_data,
                    TERM2NAME = path_terms,
                    seed = TRUE) |>
      as.data.frame() |>
      dplyr::mutate(coefficient = coef,
                    regulation = case_when((NES > 0) ~ "upregulated",
                                           (NES < 0) ~ "downregulated",
                                           T ~ "nonDE"),
                    pathway = str_trim(str_to_title(str_replace_all(ID, "_", " "))))
    
    df_list <- append(df_list, list(df_gsea))
  }
  df_gsea <- do.call(rbind, df_list)
  df_gsea <- df_gsea |>
    mutate(coefficient = factor(coefficient, levels = resultsNames(dds)))
  return(df_gsea)
}

plot_regulation <- function(df, min_alpha, regulation_type, n, fill, title) {
  plot <- df |>
    dplyr::filter(p.adjust < min_alpha &
                  regulation == regulation_type &
                  str_detect(coefficient, "^subject") == FALSE) |>
    dplyr::mutate(pathway = str_replace_all(pathway, c("Kegg Medicus "="", "Hallmark "=""))) |>
    group_by(coefficient) |>
    slice_max(abs(NES), n=n) |>
    ungroup() |>
      ggplot(aes(x=reorder_within(pathway, -NES, coefficient), y=abs(NES))) +
      geom_col(fill = fill) +
      scale_x_reordered() +
      facet_wrap(~coefficient,
                 scales = 'free_x',
                 nrow = 1) +
      theme_classic() +
      theme(axis.text.x = element_text(angle=45, hjust=1, size = 8)) +
      labs(x = "pathway",
           y = "|NES|",
           title = title) +
      coord_flip()
}

plot_gsea <- function(df, min_alpha, n) {
  up_plot <- plot_regulation(df, min_alpha, regulation_type = "upregulated", n, fill = "red", title = "Upregulated Pathways")
  down_plot <- plot_regulation(df, min_alpha, regulation_type = "downregulated", n, fill = "blue", title = "Downregulated Pathways")
  return(list(up_plot, down_plot))
}