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
  Annotation_file %>% rename(ID = GeneID,
                             GO = GOterms)
  ## We remove duplicated rows
  Annotation_file <- Annotation_file[!duplicated(Annotation_file),]
  ## Transform to data frame
  Annotation_file <- as.data.frame(Annotation_file)
  ## We create an empty list
  GO_list <- list()
  ## For loop row by row, assigning GeneIDs and GO terms to objects
  for (row in seq(1, nrow(Annotation_file))){
    ID <- Annotation_file[row,1]
    GO <- Annotation_file[row,2]
    ## GeneIDs will be names of the lists and GOterms the content. If a GeneID was already included, the new GO terms will be appended
    if (ID %in% names(GO_list)){
      GO_list[[ID]] <- append(GO_list[[ID]], GO)
    } else {
      GO_list[ID] <- GO
    }

  }
  geneID2GO <- GO_list
  return(geneID2GO)
}


