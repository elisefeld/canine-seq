import_kallisto <- function(kallisto_path,
                            sample_ids,
                            tx2gene,
                            tx_out = FALSE) {
  # NOTE: this will not include samples in the file path if they are not in the sample sheet
  message("Reading abundance.h5 files from ", kallisto_path, "...")
  files <- file.path(kallisto_path,
                     sample_ids,
                     "abundance.h5") 
  if (length(files[file.exists(files)]) > 0){
    message("Found ", length(files[file.exists(files)]), " files in ", kallisto_path)
    names(files) <- sample_ids
    txi <- tximport(files,
                    type = "kallisto",
                    tx2gene = tx2gene,
                    txOut = tx_out, # change this to true if you want transcript level data
                    ignoreTxVersion = TRUE)
    return(txi)
  }
  else {
    stop("No kallisto abundance files found at ", kallisto_path, ".")
  }
}
  
import_feature_counts <- function(fc_path, sample_ids) {
  if (file.exists(fc_path)) {
  message("Reading featurecounts file from ", fc_path, "...")
  counts <- read.table(fc_path, sep = "\t", header = TRUE, check.names = FALSE) |>
    tibble::column_to_rownames("Geneid") |>
    dplyr::select(-(1:5)) |>
    dplyr::rename_with(~ str_remove_all(.x, "runs\\/.*sort\\/|_sorted.bam")) |>
    select(all_of(sample_ids)) |>
    as.matrix()
  message("Counts for ", length(counts), " genes across ", ncol(counts), " samples found.")
  return(counts)
  }
  else {
    stop("Feature counts file ", fc_path, " does not exist.")
  }
}

prep_samples <- function(sample_sheet_path) {
  if (file.exists(sample_sheet_path)) {
    message("Reading sample sheet from ", sample_sheet_path, "...")
    sample_sheet <- read.csv(sample_sheet_path)
    if ("time" %in% names(sample_sheet)) {
      if("subject" %in% names(sample_sheet)){
        sample_sheet <- sample_sheet |>
          dplyr::mutate(across(everything(), ~ factor(.x, levels = sort(unique(.x)))),
                        time_num = as.numeric(as.character(time)))
        rownames(sample_sheet) <- sample_sheet$sample_id
        message("Found ", nrow(sample_sheet), " samples in sample sheet.")
        return(sample_sheet)
      }
      else {
        stop("Sample sheet does not have a 'subject' column")
      }
     
    }
    else {
      stop("Sample sheet does not have a 'time' column.")
    }
  }
  else {
    stop("Sample sheet ", sample_sheet_path, " does not exist.")
  }
}


