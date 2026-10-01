read_tx2gene <- function(gtf_path) {
  if (file.exists(gtf_path)) {
    message('Reading annotation file from ', gtf_path, '...')
    tx2gene <- rtracklayer::readGFF(gtf_path) |>
      as.data.frame() |>
      dplyr::select(transcript_id, gene_id, gene_name, gene_biotype) |>
      dplyr::filter(!is.na(transcript_id)) |>
      dplyr::distinct() |>
      dplyr::mutate(gene_name = coalesce(gene_name, gene_id)) # if no gene name present, use gene id
    message('Found ', length(unique(tx2gene$transcript_id)), ' unique transcripts and ', length(unique(tx2gene$gene_id)), ' unique genes.')
    return(tx2gene)
  }
  else {
    stop('GTF file ', gtf_path, ' does not exist.')
  }
}

read_homologs <- function(homolog_path){
  if (file.exists(homolog_path)) {
    message('Reading csv from ', homolog_path, '...')
    homologs <- read.csv(homolog_path, sep = '\t') |>
      dplyr::filter(ref_species == 'homo_sapiens') |>
      dplyr::select(ref_gene_stable_id,
                    ref_gene_name,
                    query_gene_stable_id,
                    query_gene_name) |>
      dplyr::rename(human_gene_id = ref_gene_stable_id,
                    human_gene_name = ref_gene_name,
                    dog_gene_id = query_gene_stable_id,
                    dog_gene_name = query_gene_name) |>
      group_by(dog_gene_id) |>
      summarise(across(everything(), ~paste0(unique(.), collapse = ', '))) |>
      dplyr::mutate(has_alt_human_genes = str_detect(human_gene_id, ','),
                    alt_human_gene_ids = str_extract(human_gene_id, '(.*)'),
                    human_gene_id = str_extract(human_gene_id, '^[^,]+'),
                    alt_human_gene_names = str_extract(human_gene_name, '(.*)'),
                    human_gene_name = str_extract(human_gene_name, '^[^,]+'),
                    human_gene_name = coalesce(human_gene_name, human_gene_id))
    message('Found ', length(unique(homologs$human_gene_id)), ' homologous genes in humans.')
    return(homologs)
  }
  else {
    stop('Homolog file ', homolog_path, ' does not exist.')
  }
}

get_entrez_ids <- function(gene_names, from_type = 'SYMBOL', org = 'org.Hs.eg.db') {
  sym_to_entrez <- clusterProfiler::bitr(unique(gene_names), fromType = from_type, toType='ENTREZID', OrgDb = org) |>
    group_by(SYMBOL) |>
    summarise(across(everything(), ~paste0(unique(.), collapse = ', '))) |>
    dplyr::mutate(entrez_id = str_extract(ENTREZID, '^[^,]+')) |> # selects the first entrez id if there are more than one
    dplyr::select(SYMBOL, entrez_id) |>
    dplyr::rename(gene_id = SYMBOL) 
  return(sym_to_entrez)
}
