#' Prepare genes for topGO analysis
#'
#' This function takes a gene vector or a data.frame with differential expression results
#' and returns a binary factor that shows which genes are in the list of interest.
#'
#' @param DEG Vector of genes or data.frame with differential expression results. Genes must be rownames in the data.frame.
#' @param geneUniverse Vector with all gene names of the organism universe (e.g., `geneID2GO` names).
#' @param padj_threshold Significance threshold for `padj` if `DEG` is a data.frame.
#'
#' @return Binary factor with gene names, to be used in `topGO`.
#' @export
prepare_topGO_genes <- function(DEG, geneUniverse, padj_threshold = 0.05) {
  if (is.vector(DEG)) {
    selected_genes <- DEG
  } else if (is.data.frame(DEG)) {
    if (!"padj" %in% colnames(DEG)) {
      warning("No 'padj' column found. All genes will be considered.")
      selected_genes <- rownames(DEG)
    } else {
      selected_genes <- rownames(DEG[DEG$padj <= padj_threshold, ])
    }
  } else {
    stop("None of the selected genes match the provided gene universe. Check gene IDs and annotation.")
  }

  if (!any(selected_genes %in% geneUniverse)) {
    stop("Gene IDs must correspond with the annotation file provided")
  }
  gene_factor <- factor(as.integer(geneUniverse %in% selected_genes))
  names(gene_factor) <- geneUniverse
  return(gene_factor)
}



#' Run GO enrichment analysis using topGO
#'
#' This function creates a topGOdata object and runs the enrichment test.
#'
#' @param gene_factor Binary factor indicating genes of interest (output of `prepare_topGO_genes()`).
#' @param geneID2GO Named list mapping gene IDs to GO terms.
#' @param ontology GO ontology to use: "BP", "MF", or "CC".
#' @param algorithm Algorithm for topGO test (default: "weight01").
#' @param statistic Statistical test to use (default: "fisher").
#'
#' @return A list with the topGOdata object and the result of the test.
#' @export
run_topGO_analysis <- function(gene_factor, geneID2GO,
                               ontology = "BP",
                               algorithm = "weight01",
                               statistic = "fisher") {
  if (!requireNamespace("topGO", quietly = TRUE)) {
    stop("The 'topGO' package is required but not installed.")
  }

  if (!is.factor(gene_factor)) {
    stop("gene_factor must be a factor.")
  }

  if (!is.list(geneID2GO)) {
    stop("geneID2GO must be a named list.")
  }

  GOdata <- topGO::new("topGOdata",
                       ontology = ontology,
                       allGenes = gene_factor,
                       annot = topGO::annFUN.gene2GO,
                       gene2GO = geneID2GO)

  result <- topGO::runTest(GOdata, algorithm = algorithm, statistic = statistic)

  return(list(GOdata = GOdata, result = result))
}








#' Análisis de enriquecimiento GO con topGO
#'
#' Esta función realiza un análisis de términos GO a partir de una lista de genes.
#'
#' @param gene_list Vector containing genes of interest
#' @param all_genes Vector con todos los genes del universo
#' @param ontology Ontología GO a usar: "BP", "MF" o "CC"
#' @param Name Name for the resulting files (typically the comparison you are testing)
#' @return Data frame with enriched GO terms
#' @examples
#' # Read data with DEGs from DESeq2
#' DEGs = read.table()
#' # Perform GO enrichment of BP with Tak-1 OPDA vs Tak-1 Mock differentially expressed genes
#' read_table()
#'
#' @export




topGO_all <-
function(DEG, Name, Plot_Ancestors, geneID2GO){
  geneNames <- names(geneID2GO)
  if (is.vector(DEG) & DEG[1] %in% geneNames) {
    DEG2 <- factor(as.integer(geneNames_Mpo %in%  DEG))
    names(DEG2) <- geneNames_Mpo
    DEG3 <- new("topGOdata", ontology= "BP", allGenes= DEG2,
                annot= annFUN.gene2GO, gene2GO= geneID2GO_Mpo)
    DEG5 <- runTest(DEG3, algorithm= "weight01", statistic="fisher")
    DEG6 <- usedGO(DEG3)
    DEG8 <- score(DEG5)[score(DEG5)<=0.05]
    DEG9 <- termStat(object= DEG3, whichGO= names(DEG8))
    DEG9$Term <- Term(rownames(DEG9))
    DEG9$pval <- DEG8
    dir.create(file.path(CurrentDir, "topGO"))
    #write.table(DEG9, file= file.path(CurrentDir, paste0("topGO/topGO_weightFS_", Name, "_BP.txt")), sep= "\t", quote= FALSE, row.names= TRUE)
    write.table(DEG9 %>% dplyr::select(pval), file= file.path(CurrentDir, paste0("topGO/topGO_weightFS_", Name, "_BP_simple.txt")), sep= "\t", quote= FALSE, row.names= TRUE, col.names=FALSE)
    DEG9$Term <- factor(DEG9$Term, levels = rev(DEG9$Term))
    DEG9$Enrichment <- (DEG9$Significant / DEG9$Annotated) / (length(DEG)/length(geneID2GO_Mpo))
    DEG9 <- DEG9 %>%
      arrange(Enrichment, Significant) %>%
      mutate(Term = factor(Term, levels= Term)) %>%
      filter(Significant>1, Enrichment>1)
  } else if (is.data.frame(DEG) & rownames(DEG)[1] %in% geneNames) {
    DEG2 <- factor(as.integer(geneNames_Mpo %in%  rownames(subset(DEG, padj <= 0.05))))
    names(DEG2) <- geneNames_Mpo
    DEG3 <- new("topGOdata", ontology= "BP", allGenes= DEG2,
                annot= annFUN.gene2GO, gene2GO= geneID2GO_Mpo)
    DEG5 <- runTest(DEG3, algorithm= "weight01", statistic="fisher")
    DEG6 <- usedGO(DEG3)
    DEG8 <- score(DEG5)[score(DEG5)<=0.05]
    DEG9 <- termStat(object= DEG3, whichGO= names(DEG8))
    DEG9$Term <- Term(rownames(DEG9))
    DEG9$pval <- DEG8
    dir.create(file.path(CurrentDir, "topGO"))
    #write.table(DEG9, file= file.path(CurrentDir, paste0("topGO/topGO_weightFS_", Name, "_BP.txt")), sep= "\t", quote= FALSE, row.names= TRUE)
    write.table(DEG9 %>% dplyr::select(pval), file= file.path(CurrentDir, paste0("topGO/topGO_weightFS_", Name, "_BP_simple.txt")), sep= "\t", quote= FALSE, row.names= TRUE, col.names=FALSE)
    DEG9$Term <- factor(DEG9$Term, levels = rev(DEG9$Term))
    DEG9$Enrichment <- (DEG9$Significant / DEG9$Annotated) / (nrow(DEG)/length(geneID2GO_Mpo))
    DEG9 <- DEG9 %>%
      arrange(Enrichment, Significant) %>%
      mutate(Term = factor(Term, levels= Term)) %>%
      filter(Significant>1, Enrichment>1)
  }

  DEG10 <- ggplot(DEG9, aes(x= Enrichment, y= Term, col= pval))+
    geom_point(aes(size=Significant))+
    scale_color_gradient2(low = "#253494", mid = "#41b6c4", midpoint = 0.025, high = "#edf8b1")+
    scale_size_binned()+
    labs(title = Name)+
    theme_gray()+
    theme(panel.grid.major.y = element_line(color = "#d4d4d4", linetype = "dashed"),
          axis.title.y = element_blank(),
          axis.text.y = element_text(face = "bold"),
          plot.title = element_text(hjust = 0, face = "bold"))
  ggsave(file.path(CurrentDir, paste0("topGO/topGO_weightFS_", Name, "_plot.pdf")), plot = DEG10,
         device = "pdf", width = 10, height = round(8 * nrow(DEG9)/50), units = "in", scale= 1.5, limitsize = FALSE)
  write.table(DEG9, file.path(CurrentDir, paste0("topGO/topGO_weightFS_", Name, "_plot.txt")), quote = FALSE, sep = "\t")
  if (Plot_Ancestors == TRUE) {
    test <- get_ancestors_pval(DEG9)
    #cloud <- wordcloud2(test %>% dplyr::select(Description, pval_add) %>% rename(words = Description, freqs = pval_add),
    #                    color = rep(c("#0C5FA2", "#8FD3BD"),  as.integer(nrow(test)/2)+1), size = 0.3)
    #saveWidget(cloud,file.path(CurrentDir, paste0("topGO/topGO_weightFS_", Name, "_Ancestors_wordcloud.html")),selfcontained = F)
    #webshot::webshot(file.path(CurrentDir, paste0("topGO/topGO_weightFS_", Name, "_Ancestors_wordcloud.html")),file.path(CurrentDir, paste0("topGO/topGO_weightFS_", Name, "_Ancestors_wordcloud.png")),vwidth = 1992, vheight = 1744, delay =10)
    test2 <- test %>% arrange(pval_add) %>%
      mutate(Description=factor(Description, levels=Description))
    DEG11 <- ggplot(test2, aes(x=Description, y=pval_add)) +
      geom_segment( aes(x=Description, xend=Description, y=0, yend=pval_add))+
      #geom_point(size = 4, color= colorRampPalette(c("#edf8b1", "#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0","#225ea8", "#0c2c84"))(length(unique(grouped_terms$General_Description)))) +
      geom_point(aes(size = items), color= "#1d91c0")+
      scale_y_continuous(expand = expand_scale(mult = c(0, .1)))+
      coord_flip()+
      ylab("Abundance")+
      scale_size_area(breaks = c(1, as.integer(max(test2$items)/3), as.integer((max(test2$items)/3)*2),max(test2$items)))+
      theme_minimal()+
      theme(axis.ticks = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_line(linetype = "dotted", colour = "grey60"),
            panel.grid.major.y = element_blank())
    ggsave(file.path(CurrentDir, paste0("topGO/topGO_weightFS_", Name, "_Ancestors_plot.pdf")), plot = DEG11,
           device = "pdf", width = 10, height = round(8 * nrow(test2)/50), units = "in", scale= 1.5, limitsize = FALSE)
  }
  if (length(rownames(DEG9))>1) {
    simMatrix <- calculateSimMatrix(rownames(DEG9),
                                    orgdb="org.At.tair.db",
                                    ont="BP",
                                    semdata = save,
                                    method="Rel")
    scores <- setNames(-log10(DEG9$pval), rownames(DEG9))
    reducedTerms <- reduceSimMatrix(simMatrix,
                                    scores,
                                    threshold=0.5,
                                    orgdb="org.At.tair.db")
    if (nrow(simMatrix)>2) {
      x <- cmdscale(as.matrix(as.dist(1-simMatrix)), eig=TRUE, k=2)$points
      df <- cbind(as.data.frame(x),
                  reducedTerms[match(rownames(x), reducedTerms$go), c("term", "parent", "parentTerm", "score")]) %>% arrange(-score)
      top_terms <- df$term[1:20]
      if (nrow(df)<=10) {
        viz <- fviz_nbclust(df[,1:2], kmeans, method='silhouette', k.max = nrow(df)-1)
      } else {
        viz <- fviz_nbclust(df[,1:2], kmeans, method='silhouette')
      }
      nclust <- as.integer(viz$data[viz$data$y == max(viz$data$y),"clusters"])
      df$clust <- as.factor(kmeans(df[,1:2], centers=nclust)$cluster)
      df <- df %>%
        mutate(label = ifelse(term %in% top_terms, parentTerm, ""))
      scatterplot <- ggplot(df, ggplot2::aes(x=V1, y=V2, color=clust)) +
        geom_point(aes(size= score), alpha=.5) +
        scale_color_manual(guide="none", values = as.vector(pal_npg(palette = "nrc")(10))) +
        scale_size_continuous(guide="none", range=c(0, 25)) +
        scale_x_continuous(name="") +
        scale_y_continuous(name="") +
        theme_minimal() +
        theme(axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              panel.grid = element_blank())+
        ggrepel::geom_label_repel(aes(label=label),
                                  data=subset(df, parent == rownames(df)),
                                  box.padding=grid::unit(1, "lines"), size=10, max.overlaps = 100)
      ggsave(file.path(CurrentDir, paste0("topGO/topGO_weightFS_", Name, "_Scatterplot.pdf")), plot = scatterplot,
             device = "pdf", width = 18, height = 10, units = "in", scale= 1.5, limitsize = FALSE)

    }
    pdf(file.path(CurrentDir, paste0("topGO/topGO_weightFS_",Name, "_Treemap.pdf")),
        height= 8, width= 10)
    treemap::treemap(reducedTerms, index=c("parentTerm", "term"), vSize="score",
                     type="index", title="",
                     palette=rep(as.vector(pal_npg(palette = "nrc")(10)), length(unique(reducedTerms$parent))),
                     #palette=gg_color_hue(length(unique(reducedTerms$parent))),
                     fontcolor.labels=c("#FFFFFFDD", "#00000080"), bg.labels=0,
                     border.col="#00000080")
    dev.off()
    return(df)
  }
}
