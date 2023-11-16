#' PlotGeneSets
#'
#' Plots a heatmap with the distances per cluster
#' @import ComplexHeatmap
#' @import RColorBrewer
#' @import colorRamp2
#' @importFrom grDevices colorRampPalette
#'
#' @param Object a pathway object
#' @param fontsize a numeric with the fontsize for plotting
#' @param legend  add a legend to the plot
#' @param annotation.mol should the genes from the genes set be added to the plot.
#' @param main is the plot title
#' @param RR.max is the maximum distance size to be added, it cutsoff the max to this. Usefull if the distance score gets very high.
#' @param cluster.order is a user defined order of the clusters in the heatmap.
#'
#' @return plot
#'
#' @export
#'
setGeneric(name="PlotGeneSets",
           def=function(Object, fontsize = 5,
                        legend = T,
                        annotation.mol=F,
                        main="",
                        RR.max = "",
                        cluster.order = "",
                        doORA = T,
                        wordcloud = T)
           {
             standardGeneric("PlotGeneSets")
           }
)

#' PlotGeneSets
#'
#' @param Object a PathwayObject
#' @param PathwayObject a PathwayObject
#' @param fontsize a numeric with the fontsize for plotting
#' @param legend  add a legend to the plot
#' @param annotation.mol should the genes from the genes set be added to the plot.
#' @param main is the plot title
#' @param RR.max is the maximum distance size to be added, it cutsoff the max to this. Usefull if the distance score gets very high.
#' @param cluster.order is a user defined order of the clusters in the heatmap.
#'
#' @return plot
#' @export
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
#'  PlotGeneSets(Object = man.Great.Object3, fontsize =5,
#'               legend = TRUE,
#'               annotation.mol=FALSE,
#'               main="man.Great.Object3",
#'               RR.max = 50)
setMethod(f="PlotGeneSets",
          signature="PathwayObject",
          definition=function(Object, fontsize = 5,
                              legend = T,
                              annotation.mol=F,
                              main="",
                              RR.max = "",
                              cluster.order = "",
                              doORA = T,
                              wordcloud = T)
          {
            Data.RR <- Object@Data.RR
            if(!RR.max == "")
            {
              RR.max <- as.numeric(as.character(RR.max))
              for(rows.i in 1:nrow(Data.RR))
              {
                idx <- Data.RR[rows.i,] > RR.max
                names(idx) <- NULL
                Data.RR[rows.i,idx] <- RR.max
              }
            }
            if(!cluster.order[1]== "")
            {
              if(!length(cluster.order) == max(as.numeric(as.character(Object@plot$aka2$Cluster))))
              {
                message("Make sure youre the order of supplied clusters is the same length as the number of clusters")

              }else{
                order.x <- vector()
                for(order.i in cluster.order)
                {
                  order.x <- c(order.x, which(order.i == as.numeric(as.character(Object@plot$aka2$Cluster))))
                }
                Object@plot$aka2 <- Object@plot$aka2[order.x,]

                Data.RR <- Data.RR[order.x, order.x]
              }
            }

            #performing ora
            clus <- Object@Data[[1]]$cluster
            names(clus) <- Object@Data[[1]]$Pathways
            clus<-gsub("Cluster_", "", clus)
            print(length(unique(Object@Data[[1]]$cluster)))
            keywords_ora <- ORAperCluster(Object, doORA, length(unique(Object@Data[[1]]$cluster)), clus)

            go_id_list <- list()
            for (z in 1:length(unique(clus)))
            {
              go_id_list[[z]] <- names(clus[clus==z])
            }


            if (wordcloud == T)
            {
              #performing keyword erichment to generate wordcloud
              if (checkGO(Object) == FALSE) {
                message("No GO terms have been detected in the pathways. The semantic enrichment word cloud will not be generated.")
                wordcloud = FALSE
              } else {
                #adapted from anno_word_cloud_from_GO function
                env_tdm_GO <- readRDS(system.file("extdata", "tdm_GO.rds", package = "simplifyEnrichment"))
                names(go_id_list) <- as.character(1:length(go_id_list))

                # keyword enrichment
                message(paste0("Performing keyword enrichment for "), length(go_id_list), " group(s) of pathways.")
                term <- lapply(go_id_list, function(x, min_stat=5) {
                  df <- keywordEnrichment(x, env_tdm_GO)
                  df <- df[df$p <= min_stat, , drop = FALSE]
                  data.frame(df[, 1], -log10(df$p))
                })
              }

              #plot labels
              clus2 <- Object@Data[[1]]$cluster
              clus2<-gsub("Cluster_", "", clus2)
              names(clus2) <- Object@Data[[1]]$RR_name
              go_id_list2 <- list()
              for (z in 1:length(unique(clus2)))
              {
                go_id_list2[[z]] <- names(clus2[clus2==z])
              }

              names(go_id_list2) <- as.character(1:length(go_id_list2))
              align_to <- go_id_list2
              for (i in 1:length(align_to)){
                for (j in 1:length(align_to[[i]]))
                  align_to[[i]][j] <- as.numeric(which(colnames(Data.RR) == align_to[[i]][j]))
              }

              align_to <- lapply(align_to, as.numeric)

              annot_label <- rowAnnotation(keywords = anno_word_cloud(align_to, term),
                                           annotation_name_align=T)

            } else {
              annot_label <- NULL
            }

            # rowNumbers <- rowAnnotation(foo = anno_mark(at = unlist(lapply(align_to, function(x) round(mean(x)))),
            #                                             labels = rep(paste0("Group ", 1:length(align_to))),
            #                                             link_width = unit(0,"mm")))

            if (length(unique(Object@Data[[1]]$Groups)) < 3)
            {
              my_palgroups <- brewer.pal(3, name = "Set3")
              my_palgroups <- my_palgroups[1:length(unique(Object@Data[[1]]$Groups))]
              my_palcluster <- brewer.pal(3, name = "Paired")
              my_palcluster <- my_palcluster[1:length(unique(Object@Data[[1]]$Groups))]
            } else {
              my_palgroups <- brewer.pal(length(unique(Object@Data[[1]]$Groups)), name = "Set3")
              my_palcluster <- brewer.pal(length(unique(Object@Data[[1]]$cluster)), name = "Paired")
            }

            groups <- my_palgroups[1:length(unique(Object@Data[[1]]$Groups))]
            names(groups) <- unique(Object@Data[[1]]$Groups)
            cluster <- my_palcluster[1:length(unique(Object@Data[[1]]$cluster))]
            names(cluster) <- gsub("Cluster_","",unique(Object@Data[[1]]$cluster))
            #names(cluster) <- unique(Object@Data[[1]]$cluster)
            colors <- list(Group = groups, Cluster = cluster)

            color_fun <- colorRamp2(seq(min(Data.RR), max(Data.RR), length = 9), brewer.pal(9, "Reds"))
            annot_top <- HeatmapAnnotation(df=Object@plot$aka2, show_legend = T, col = colors)
            row_annot <- rowAnnotation(df=Object@plot$aka2, show_legend = F, col = colors)

            if(annotation.mol==T)
            {
              if (wordcloud == T)
              {
                plot <- Heatmap(matrix = Data.RR, cluster_rows = F, cluster_columns = F, show_row_names = T, show_column_names = F ,
                                name="Similarity", col=color_fun, top_annotation = annot_top, left_annotation = row_annot, right_annotation = annot_label)
              } else {
                plot <- Heatmap(matrix = Data.RR, cluster_rows = F, cluster_columns = F, show_row_names = T, show_column_names = F ,
                                name="Similarity", col=color_fun, top_annotation = annot_top, left_annotation = row_annot)
              }
            } else {

              if(wordcloud == T)
              {
                plot <- Heatmap(matrix = Data.RR, cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F ,
                                name="Similarity", col=color_fun, top_annotation = annot_top, left_annotation = row_annot, right_annotation = annot_label)
              } else {


                plot <- Heatmap(matrix = Data.RR, cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F ,
                                name="Similarity", col=color_fun, top_annotation = annot_top, left_annotation = row_annot)
              }
            }

            if (doORA == TRUE) {
              legend <- Legend(labels = rep(paste0(1:length(keywords_ora), "- ", keywords_ora)),
                               title = "\nORA", legend_gp = gpar(fill = 1:length(keywords_ora)),
                               nr=8,  title_position = "leftcenter-rot")

              mergedplot <- draw(plot, annotation_legend_list=legend, annotation_legend_side = "bottom")
            } else {
              mergedplot <- plot
            }

            return(mergedplot)

}
)

