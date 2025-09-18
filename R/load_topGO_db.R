#' Creation of a correct database for topGO
#'
#' This function adapts your annotation file from your organism of interest to the proper formats that topGO requieres for GO enrichment analysis
#'
#' @param Annotation_file Data frame including the annotation data of your organism. This must contain at least one column with the GeneID and another with GO term(s)
#' @param GeneID Name of the column with GeneIDs
#' @param GOterms Name of the column with GOterms representing GeneIDs
#' @return List in the appropiate topGO format
#' @examples
#' # Read data with DEGs from DESeq2
#' DEGs = read.table()
#' # Perform GO enrichment of BP with Tak-1 OPDA vs Tak-1 Mock differentially expressed genes
#' read_table()
#'
#' @export


load_topGO_db <- function(Annotation_file, GeneID, GOterms){
  # First, we rename columns to make it easier
  Annotation_file %>% dplyr::rename(ID = GeneID,
                                    GO = GOterms)
  ## We remove duplicated rows
  Annotation_file <- Annotation_file[!duplicated(Annotation_file),]
  ## Transform to data frame
  Annotation_file <- as.data.frame(Annotation_file)
  geneID2GO <- split(Annotation_file[,GOterms], Annotation_file[,GeneID])
  return(geneID2GO)
}

