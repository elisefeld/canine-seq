get_msig <- function(species = 'Homo sapiens',
                     collection,
                     subcollection = NULL){
  msig <- msigdbr(species = species,
                  collection = collection,
                  subcollection = subcollection)
  term2gene <- msig |>
    dplyr::select(gs_name, all_of('gene_symbol'))
  term2name <- msig |>
    dplyr::select(gs_name, gs_description)
  msig_list <- list(term2gene, term2name)
  names(msig_list) <- c('term2gene', 'term2name')
  return(msig_list)
}

create_ranked_list <- function(df, coef){
  ranked_genes <- df |>
    dplyr::filter(coefficient == coef & !is.na(gene_id)) |>
    arrange(desc(log2FoldChange)) |>
    pull(log2FoldChange, name = gene_id)
}

run_ora <- function(genes, universe){
  genes <- genes[!is.na(genes)]
  
  enrichGO(gene = genes,
           OrgDb = org.Hs.eg.db,
           keyType = 'ENTREZID',
           ont = 'BP',
           universe = universe,
           readable = TRUE) |>
    clusterProfiler::simplify()
}

run_gsea <- function(coefs,
                     df,
                     species,
                     collection,
                     subcollection = NULL) {
  
  msig <- get_msig(species = species,
                   collection = collection,
                   subcollection = subcollection)
  
  df <- purrr::map_dfr(coefs,
                 run_gsea_coef,
                 df = df,
                 term2gene = msig$term2gene,
                 term2name = msig$term2name)
  
  df <- df |>
    dplyr::mutate(coefficient = factor(coefficient, levels = coefs))
}


run_gsea_coef <- function(coef,
                          df,
                          term2gene,
                          term2name) {
  
  ranked_genes <- create_ranked_list(df, coef)
  
  GSEA(ranked_genes,
       TERM2GENE = term2gene,
       TERM2NAME = term2name,
       seed = TRUE) |>
    as.data.frame() |>
    dplyr::mutate(coefficient = coef,
                  regulation = case_when(NES > 0 ~ 'upregulated',
                                         NES < 0 ~ 'downregulated',
                                         TRUE ~ 'nonDE'),
                  pathway = str_trim(str_to_title(str_replace_all(ID, '_', ' '))))
}

plot_gsea_heatmap <- function(df, min_alpha = 0.05, n = 20){
  
  top_pathways <- df |>
    dplyr::filter(p.adjust < min_alpha) |>
    group_by(pathway) |>
    summarise(min_alpha = min(p.adjust), na.rm = TRUE) |>
    slice_max(min_alpha, n = n) |>
    pull(pathway)
  
  df |>
    dplyr::filter(pathway %in% top_pathways) |>
    mutate(pathway = str_remove_all(pathway,'Kegg Medicus |Hallmark '),
           coefficient = str_remove_all(coefficient, '^time_|_vs_0')) |>
    ggplot(aes(x = coefficient,y = forcats::fct_reorder(pathway, abs(NES), .fun = max), fill = NES)) +
    geom_tile(color = 'white') +
    scale_fill_gradient2(
      low = '#2166AC',
      mid = 'white',
      high = '#B2182B',
      midpoint = 0) +
    theme_classic(base_size = 12) +
    labs(x = 'Time (days)',
         y = NULL,
         fill = 'NES')
}


plot_regulation <- function(df, min_alpha, regulation_type, n, fill, title) {
  plot <- df |>
    dplyr::filter(p.adjust < min_alpha &
                  regulation == regulation_type &
                  str_detect(coefficient, '^subject') == FALSE) |>
    dplyr::mutate(pathway = str_replace_all(pathway, c('Kegg Medicus '='', 'Hallmark '=''))) |>
    group_by(coefficient) |>
    slice_max(abs(NES), n=n) |>
    ungroup() |>
      ggplot(aes(x=reorder_within(pathway, -NES, coefficient), y=abs(NES))) +
      geom_bar(fill = fill, stat = 'identity') +
      scale_x_reordered() +
      facet_wrap(~coefficient,
                 scales = 'free_x',
                 nrow = 1) +
      theme_classic() +
      theme(axis.text.x = element_text(angle=45, hjust=1, size = 8)) +
      labs(x = 'pathway',
           y = '|NES|',
           title = title) +
      coord_flip()
}

plot_gsea <- function(df, min_alpha, n) {
  up_plot <- plot_regulation(df, min_alpha, regulation_type = 'upregulated', n, fill = 'red', title = 'Upregulated Pathways')
  down_plot <- plot_regulation(df, min_alpha, regulation_type = 'downregulated', n, fill = 'blue', title = 'Downregulated Pathways')
  return(list(up_plot, down_plot))
}
