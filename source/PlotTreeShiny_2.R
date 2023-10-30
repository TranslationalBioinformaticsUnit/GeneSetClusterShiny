#' PlotTreeShiny_2
#'
#' Plots a tree plot with the distances per cluster
#' @import ggplot2
#' @import ggnewscale
#' @import ggtree
#' @import GO.db
#'
#' @param Object a pathway object
#' @param fontsize a numeric with the fontsize for plotting nodenames
#' @param main is the plot title
#' @param nodenames boolean to add names to the nodes
#'
#' @return plot
#'
#' @export
#'
setGeneric(name="PlotTreeShiny_2",
           def=function(Object, fontsize = 3,
                        main="",
                        nodenames= TRUE)
           {
             standardGeneric("PlotTreeShiny_2")
           }
)

#' PlotTreeShiny_2
#'
#' @param Object a pathway object
#' @param fontsize a numeric with the fontsize for plotting rownames/colnames
#' @param main is the plot title
#' @param nodenames boolean to add names to the nodes
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
#'                                      method = "Hierarchical")
#'  PlotTreeShiny_2(Object = man.Great.Object3, fontsize = 3,
#'               main= "man.Great.Object3",
#'               nodenames = TRUE)


setMethod(f="PlotTreeShiny_2",
          signature="PathwayObject",
          definition=function(Object, fontsize = 3,
                              main="",
                              nodenames = TRUE)
  {

            
    
    maximun<-200
    npathways<-dim(Object@Data[[1]])[1]
    if (npathways>maximun) {
      print("Size greater than 200. No treeplot perform")
      NULL
    }else{
      metadata<-as.data.frame(Object@Data[[1]][,c("RR_name","Ratio")])

      clus.x<-hclust(dist(t(Object@Data.RR)), method = "ward.D2")
      
      p<-ggtree(clus.x) %<+% metadata + geom_tippoint(aes(size=Ratio)) + theme_tree2() + ggtitle(main) + theme(plot.title = element_text(hjust=0.5, face="bold"))
      offset<-0.05

      for (i in 1:dim(Object@plot$aka2)[2]){
        p<-gheatmap(p,Object@plot$aka2[i], offset=offset, colnames_position="top", colnames_offset_y=0.05, width=0.05, font.size=fontsize) + scale_fill_manual(values=Object@plot$aka3[[i]], name=names(Object@plot$aka3)[i])
        #calculate offset
        if ((npathways>20) & (npathways<60)){
          offset<-npathways+30
        }else if(npathways<=20){
          offset<-npathways+5
        }else{
          offset<-npathways+10
        }
        p<-p + new_scale_fill()
      }
      p
    }

  }
)