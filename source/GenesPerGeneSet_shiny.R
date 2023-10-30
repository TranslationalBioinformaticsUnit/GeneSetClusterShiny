#' GenesPerGeneSet_shiny
#'
#' Extracts a data.frame from the object with every gene in per cluster
#' @param Object a PathwayObject
#'
#' @return dataframe of pathwayobjects
#' @export
#'

setGeneric(name="GenesPerGeneSet_shiny",
           def=function(Object)
           {
             standardGeneric("GenesPerGeneSet_shiny")
           }
)

#' GenesPerGeneSet_shiny
#'
#' @param Object a PathwayObject
#' @param PathwayObject  a PathwayObject
#'
#' @return dataframe of pathwayobjects
#' @export
#'
#' @examples
#' IPA.files <- c(system.file("extdata",
#'                            "MM10.IPA.KO.uGvsMac.Canonical_pathways.xls",
#'                             package = "GeneSetCluster"),
#'              system.file("extdata",
#'                             "MM10.IPA.WT.uGvsMac.Canonical_pathways.xls",
#'                              package = "GeneSetCluster"),
#'              system.file("extdata",
#'                              "MM10.IPA.KO.uGvsMac.Functional_annotations.xls",
#'                              package = "GeneSetCluster"),
#'              system.file("extdata",
#'                              "MM10.IPA.WT.uGvsMac.Functional_annotations.xls",
#'                              package = "GeneSetCluster"))
#' canonical.files <- IPA.files[grep("Canonical", IPA.files)]
#'
#' IPA.object1 <- LoadGeneSets(file_location = canonical.files,
#'                          groupnames= c("KO", "WT"),
#'                          P.cutoff = 1.3,
#'                          Mol.cutoff = 5,
#'                          Source = "IPA",
#'                          type = "Canonical_Pathways",
#'                          structure = "SYMBOL",
#'                          seperator = ",")
#'IPA.object2 <- CombineGeneSets(Object = IPA.object1)
#'IPA.object3 <- ClusterGeneSets(Object = IPA.object2,
#'                               clusters = 7,
#'                               method = "kmeans")
#' GenesPerGeneSet_shiny(Object =IPA.object3 )
setMethod(f="GenesPerGeneSet_shiny",
          signature="PathwayObject",
          definition=function(Object)
          {
            message("[=============================================================]")
            message("[<<<<               GenesPerGeneSet_shiny START                >>>>>]")
            message("---------------------------------------------------------------")
            message(paste("Extracting genes for all ",max(Object@Data[[1]]$cluster), " clusters", sep=""))
            clus.mol.ls <- list()
            #for(clus.i in 1:max(Object@Data[[1]]$cluster))
            for(clus.i in 1:max(Object@plot$aka2[["Cluster"]]))
            {
              clus.x <- Object@Data[[1]][Object@Data[[1]]$cluster == paste("Cluster_",clus.i, sep=""),]
              clus.x$Molecules <- as.character(clus.x$Molecules)
              clus.mol.ls[[clus.i]] <- unique(as.vector(strsplit2(clus.x$Molecules, split = Object@metadata$seperator[1])))
            }
            unique.mol <- unique(do.call(what = c, args = clus.mol.ls))

            mol.unique.df <- as.data.frame(matrix(0, nrow = length(unique.mol), ncol = length(clus.mol.ls)))
            rownames(mol.unique.df) <- unique.mol
            colnames(mol.unique.df) <-  paste("Cluster_", 1:length(clus.mol.ls), sep="")
            for(clus.i in 1:max(Object@plot$aka2[["Cluster"]]))
            {
              mol.unique.df[clus.mol.ls[[clus.i]],clus.i] <- 1
            }

            message("--------------------------------------------------------------")
            message("[<<<<               GenesPerGeneSet_shiny ends                >>>>>]")
            message("[============================================================]")

            return(mol.unique.df)          }
)
