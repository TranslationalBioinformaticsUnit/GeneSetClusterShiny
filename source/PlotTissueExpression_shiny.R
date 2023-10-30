#' PlotTissueExpression
#'
#' Plots a the tissue expression total median expression per tissue and the percentage of each cluster in each tissue
#' @import patchwork
#' @import ggplot2
#' @import reshape2
#' @import RColorBrewer
#'
#' @param Object a pathway object
#' @param all boolean to add all tissues
#'
#' @return plot
#'
#' @export
#'
setGeneric(name="PlotTissueExpression",
           def=function(Object, all = FALSE)
           {
             standardGeneric("PlotTissueExpression")
           }
)

#' PlotTissueExpression
#'
#' @param Object a pathway object
#' @param all boolean to add all tissues
#'
#' @return plot
#'
#' @examples
#' Great.files <- c(system.file("extdata", "MM10.GREAT.KO.uGvsMac.bed.tsv",
#'                              package = "GeneSetCluster"),
#' system.file("extdata", "MM10.GREAT.KO.uGvsMac.bed_BCKGRND.tsv", package = "GeneSetCluster"),
#' system.file("extdata", "MM10.GREAT.WT.uGvsMac.bed.tsv", package = "GeneSetCluster"),
#' system.file("extdata", "MM10.GREAT.WT.uGvsMac.bed_BCKGRND.tsv", package = "GeneSetCluster"))
#' Great.files.bckgrnd <- Great.files[grepl("BCKGRND", Great.files)]
#'
#'
#' Great.bckgnrd.Object1 <- LoadGeneSets(file_location = Great.files.bckgrnd,
#'                                       groupnames= c("KO", "WT"),
#'                                       P.cutoff = 0.05,
#'                                       Mol.cutoff = 5,
#'                                       Source = "Great",
#'                                       Great.Background = TRUE,
#'                                       type = "Canonical_Pathways",
#'                                     topranks = 20,
#'                                    structure = "SYMBOL",
#'                                    Organism = "org.Mm.eg.db",
#'                                    seperator = ",")
#' man.Great.Object1 <- ManageGeneSets(Object = Great.bckgnrd.Object1,
#'                                    keep.type =c("Disease Ontology",
#'                                    "GO Biological Process" ),
#'                                    exclude.type="")
#' man.Great.Object2 <- CombineGeneSets(Object = man.Great.Object1)
#' man.Great.Object3 <- ClusterGeneSets(Object = man.Great.Object2,
#'                                      clusters = 5,
#'                                      method = "kmeans")
#'
#' man.Great.Object4 <- TissueExpressionPerGeneSet(man.Great.Object3)
#' PlotTissueExpression(man.Great.Object4, all = F)


setMethod(f="PlotTissueExpression",
          signature = "PathwayObject",
          definition = function(Object, all=FALSE)
  {

    if (is.null(Object@dfTissue) | anyNA(Object@dfTissue)) {
      stop("First you have to run TissueExpressionPerGeneSet to obtain tissue expression information.")
    }

    df <- Object@dfTissue
    df$tissue <- rownames(df)
    plot_melt_mtx = melt(df, id="tissue")
    # plot_melt_mtx = melt(Object2@dfTissue)
    colnames(plot_melt_mtx)[1:3] = c("cluster", "Dataset", "value")
    plot_cluster_size <- aggregate(value ~ cluster, data = plot_melt_mtx, FUN = sum)
    colnames(plot_cluster_size) = c("cluster", "value")

    df_tissue_z <- apply(Object@dfTissue, 2, function(x) (x-mean(x))/sd(x))
    df_tissue_z[is.na(df_tissue_z)] <- 0
    plot_melt_mtx_z <- melt(df_tissue_z)
    colnames(plot_melt_mtx_z)[1:3] <- c("cluster", "Dataset", "value")
    plot_cluster_size_z <- aggregate(value ~ cluster, data = plot_melt_mtx_z, FUN = sum)
    colnames(plot_cluster_size_z) <- c("cluster", "value")

    if (all == FALSE) {
      plot_cluster_size <- plot_cluster_size[order(plot_cluster_size[,"value"], decreasing = T),]
      keep <- plot_cluster_size$cluster[1:20] # top 20
      plot_cluster_size <- plot_cluster_size[plot_cluster_size$cluster %in% keep,]

      plot_cluster_size_z <- plot_cluster_size_z[plot_cluster_size_z$cluster %in% keep,]
      plot_cluster_size_z <- plot_cluster_size_z[rownames(plot_cluster_size),]
      plot_melt_mtx_z <- plot_melt_mtx_z[plot_melt_mtx_z$cluster %in% keep,]
    }

    sorted_labels <- paste(sort(levels(plot_cluster_size$cluster), decreasing = F))
    sorted_labels <- paste(sort(plot_cluster_size$cluster, decreasing = F))

    plot_cluster_size$cluster <- factor(plot_cluster_size$cluster,levels = sorted_labels)
    plot_cluster_size_z$cluster <- factor(plot_cluster_size_z$cluster,levels = sorted_labels)
    plot_melt_mtx_z$cluster <- factor(plot_melt_mtx_z$cluster,levels = sorted_labels)
    plot_melt_mtx_z$tot_Ncells <- ave(plot_melt_mtx_z$value, plot_melt_mtx_z$cluster, FUN=sum)
    plot_melt_mtx_z$percent <- plot_melt_mtx_z$value*100/plot_melt_mtx_z$tot_Ncells

    cluster_size <- plot_cluster_size
    cluster_size_z <- plot_cluster_size_z
    melt_mtx_z <- plot_melt_mtx_z
    cluster_size[order(cluster_size$value, decreasing = T),]

    my_pal <- brewer.pal(ncol(Object@dfTissue), name = "Paired")

    #Plot total proportion
    p1 <- ggplot(plot_cluster_size, aes(x = value, y = reorder(cluster, value)))

    #Plot total z-score proportion
    p2 <- ggplot(plot_cluster_size_z, aes(x = value, y = reorder(cluster, value)))

    #Plot Barplot by cluster
    p3 <- ggplot(plot_melt_mtx_z, aes(x=reorder(cluster, value), y=percent, fill=Dataset))

    p1 <- p1 + geom_bar(position="dodge", stat="identity",fill = "grey60", width = 0.85) +
      theme_bw() + xlab("Sum of NES values") +
      geom_text(aes(y=cluster,x=1.1,label=round(value, digits = 2),hjust="bottom"), size=4) + ylab("")+
      theme(axis.title = element_text(size = 16), axis.text.x = element_text(size = 12),
            axis.text.y = element_blank(), axis.ticks.y=element_blank())

    p2 <- p2 + geom_bar(position="dodge", stat="identity",fill = "grey60", width = 0.85) +
      theme_bw() + xlab("Sum of Z-score") +
      geom_text(aes(y=cluster,x=1.1,label=round(value, digits = 2),hjust="bottom"), size=4) + ylab("")+
      theme(axis.title = element_text(size = 16), axis.text.x = element_text(size = 12),
            axis.text.y = element_blank(), axis.ticks.y=element_blank())

    p3 <- p3 + geom_bar(position="dodge", stat="identity", width = 0.7) + theme_bw() + coord_flip() +
      scale_fill_manual(values = my_pal) +
      ylab("Percentage of cluster Z-score in tissue expression") + xlab("Tissue") +
      theme(legend.position="top", axis.title = element_text(size = 16),
            axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 14))

    #Combine plots
    finalPlot <- p3 + p2 + p1 + plot_layout(widths = c(4, 2, 2))

    message("[You may want to plot the cluster tissue expression using PlotClusterTissueExpression next.]")

    return(finalPlot)
}
)
