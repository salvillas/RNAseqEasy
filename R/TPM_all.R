#' Join TPM values from all samples in a data frame
#'
#' This function reads TPM quantified data from Salmon and joins them into a final data frame.
#' Results are optionally saved in a .txt file
#'
#' @param sampleDir Directory containing sample folders with quant.sf files.
#' @param sample_table Data frame with sample metadata.
#' @param Include A vector of condition(s) to include samples from `sample_table`.
#' It will include samples matching this/these condition(s) and exclude all others.
#' Default to NULL.
#' @param Exclude A vector of condition(s) to exclude samples from `sample_table`.
#' Default to NULL.
#' @param tx2gene A data frame with columns TXNAME and GENEID (transcript-to-gene mapping).
#' @param Save_results Logical. Whether to save the final data frame in a .txt file. Default is TRUE
#' @param output_path Directory where final data frame will be saved.
#'
#' @return No return value. Files are written to disk.
#' @export
TPM_all <- function(sampleDir, sample_table, output_path, Include = NULL, Exclude = NULL, Save_results = TRUE, tx2gene) {
  sample_table_some <- get_sample_subset(sample_table, Include = Include,
                                         Exclude = Exclude)
  sample_table_some <- add_sample_path(sampleDir = sampleDir, sample_table = sample_table_some)
  txi <- load_tximport_data(sampleDir, sample_table_some, tx2gene)
  TPM_df <- as.data.frame(txi$abundance)
  if (Save_results) {
    write.table(TPM_df, file = file.path(output_path, "TPM_all_samples.txt"), quote = FALSE, sep = "\t", dec = ",", row.names = TRUE)
    message("TPM data saved to ", file.path(output_path, "TPM_all_samples.txt"))
  }
  return(TPM_df)
}
