#' TissueExpressionPerGeneSet
#'
#' Extracts every gen of every cluster from the object. Using the median gene expression from the GTEx database build a data.frame with every tissue per cluster.
#'
#' @import jsonlite
#' @import httr
#' @import reshape2
#' @import dplyr
#' @import utils
#' @import pbapply
#'
#' @param Object a PathwayObject
#' @param localDatabase data frame with the GTEx database information
#'
#' @return dataframe of pathwayobjects with tissue expression per cluster
#'

setGeneric(name="TissueExpressionPerGeneSet",
           def=function(Object, localDatabase=NULL, dic=NULL)
           {
             standardGeneric("TissueExpressionPerGeneSet")
           }
)

#' TissueExpressionPerGeneSet
#'
#' @param Object a PathwayObject
#' @param localDatabase data frame with the GTEx database information
#'
#' @return dataframe of pathwayobjects with tissue expression per cluster
#' @export
#'
#' @examples
#' #' IPA.files <- c(system.file("extdata",
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
#' IPA.object2 <- CombineGeneSets(Object = IPA.object1)
#' IPA.object3 <- ClusterGeneSets(Object = IPA.object2,
#'                               clusters = 7,
#'                               method = "kmeans")
#' IPA.object4 <- TissueExpressionPerGeneSet(Object = IPA.object3, threads = 8)
setMethod(f="TissueExpressionPerGeneSet",
          signature="PathwayObject",
          definition=function(Object, localDatabase=NULL, dic=NULL)
          {
          message("[=========================================================]")
          message("[<<<<         ObtainTissueExpression START           >>>>>]")
          message("-----------------------------------------------------------")

          options <- c("SYMBOL", "ENTREZID", "ENSEMBLID", "ENTREZ", "ENSEMBL")
          entrezOptions <- c("ENTREZID", "ENTREZ")
          ensemblOptions <- c("ENSEMBL", "ENSEMBLID")

          if (!toupper(Object@metadata$structure[1]) %in% options) {
            message("Genes are not included as genes Symbol, EntrezID or ENSEMBLID.\nSo TissuePerGeneSet will NOT be performed.")
            message("Please make sure that you are using one of the followings: GENE SYMBOL, ENTREZID or ENSEMBLID.")
            stop()
          }

          #data("dic", package = "GeneSetCluster")

          mol.unique.df <- GenesPerGeneSet_shiny(Object)
          genes <- rownames(mol.unique.df)

          if (substr(genes[1], 1, 3) == "ENS") {
            dic.selected <- dic[which(dic$GTEx.median.TPM.Name %in% toupper(genes)),]
            genes.selected <- dic.selected$GTEx.median.TPM.Name
            touse <- mol.unique.df[which(toupper(genes) %in% dic$GTEx.median.TPM.Name),]

            #label number of genes per cluster
            ngenes <- apply(touse, 2, function(x) length(which(x==1)))
            colnames(touse) <- paste0(names(ngenes), " (", ngenes, ")")

            touse$ENSID <- rownames(touse)
            meltedtouse <- melt(touse, id="ENSID")
            meltedtouse <- meltedtouse[which(meltedtouse$value==1), c("variable", "ENSID")]
            colnames(meltedtouse) <- c("Term", "Gene")

          } else if (toupper(Object@metadata$structure[1]) %in% entrezOptions) {
            usingOrg <- obtainOrg(Object)
            resSymbol <- AnnotationDbi::select(usingOrg, keys=genes, columns='SYMBOL', keytype='ENTREZID')

            dic.selected <- dic[which(dic$GTEx.median.TPM.Description %in% toupper(resSymbol$SYMBOL)),]
            genes.selected <- dic.selected$GTEx.median.TPM.Name
            colnames(dic.selected) <- c("ENSEMBLID", "SYMBOL")
            dicreference <- merge(resSymbol, dic.selected, by="SYMBOL")
            rownames(dicreference) <- dicreference$ENSEMBLID
            touse <- dicreference[genes.selected,]

            mol.unique.df$ENTREZID = rownames(mol.unique.df)
            touse = merge(mol.unique.df, touse, by="ENTREZID")
            rownames(touse) = touse$ENSEMBLID
            #touse <- subset(touse, select=-c(ENTREZID, SYMBOL, GTEx.median.TPM.Description))
            touse <- subset(touse, select=-c(ENTREZID, SYMBOL))

            #label number of genes per cluster
            ngenes <- apply(touse, 2, function(x) length(which(x==1)))
            colnames(touse)[1:ncol(touse)-1] <- paste0(names(ngenes), " (", ngenes, ")")[1:ncol(touse)-1]

            meltedtouse <- melt(touse, id="ENSEMBLID")
            meltedtouse <- meltedtouse[which(meltedtouse$value==1), c("variable", "ENSEMBLID")]
            colnames(meltedtouse) <- c("Term", "Gene")

          } else {
            mol.unique.df$GTEx.median.TPM.Description <- toupper(rownames(mol.unique.df))
            dic.selected <- dic[which(dic$GTEx.median.TPM.Description %in% toupper(genes)),]
            genes.selected <- dic.selected$GTEx.median.TPM.Name

            touse <- merge(mol.unique.df, dic.selected, by="GTEx.median.TPM.Description")
            touse <- subset(touse, select = -GTEx.median.TPM.Description)

            #label number of genes per cluster
            ngenes <- apply(touse, 2, function(x) length(which(x==1)))[1:ncol(touse)-1]
            newnames <- paste0(names(ngenes), " (", ngenes, ")")
            colnames(touse) <- c(newnames, colnames(touse[ncol(touse)]))

            meltedtouse <- melt(touse, id="GTEx.median.TPM.Name")
            meltedtouse <- meltedtouse[which(meltedtouse$value==1), c("variable", "GTEx.median.TPM.Name")]
            colnames(meltedtouse) <- c("Term", "Gene")
          }

          if (length(genes.selected) == 0) {
            stop("The genes used to create the object are not correct. Please make sure they are introduced either as GENE SYMBOL, ENTREZID or ENSEMBLID.")
          }

          if (is.null(localDatabase)) {
            #GTEx API REST query of median Gene Expression: database -> gtex_v8
            message("Performing GTEx API REST query to obtain the GTEx database.")

            url <- "http://127.0.0.1:8001/GTExdatabase"
            urlNames <- "http://127.0.0.1:8001/GTExdatabaseNames"
            tryCatch(
              {
                response <- httr::GET(url)
                GTEx.info <- fromJSON(rawToChar(response$content))

                responseNames <- httr::GET(urlNames)
                GTEx.infoNames <- fromJSON(rawToChar(responseNames$content))

                for (i in 1:length(GTEx.info))
                {
                  names(GTEx.info[[i]]) <- GTEx.infoNames[[i]]
                }

                message("Performing GSEA...\n")
                results <- pblapply(GTEx.info, function(x) clusterProfiler::GSEA(geneList=x, TERM2GENE=meltedtouse, minGSSize=5, maxGSSize=10000, verbose=F))

              }, error=function(e) {
                stop("The API service is not currently working... Use the local database instead, you can download it on: INCLUDE LINk")
              }
            )


          } else {
            tryCatch({
              message("Performing GSEA...\n")
              results <- pblapply(localDatabase, function(x) clusterProfiler::GSEA(geneList=x, TERM2GENE=meltedtouse, minGSSize=5, maxGSSize=10000, verbose=F))
            # 12 mins GREAT ds
            }, error=function(e){
              stop("Seems that the database provided is not correct. Please make sure you are using the correct one. You can download it on INDLUDE LINK")
            })
          }

          tissue.df <- data.frame(matrix(0, nrow = length(results), ncol = length(unique(Object@Data[[1]]$cluster))))
          rownames(tissue.df) <- names(results)
          colnames(tissue.df) <- colnames(touse)[1:ncol(tissue.df)]

          for (i in 1:length(results))
          {
            tissue.df[i, results[[i]]@result$ID] <- results[[i]]@result$NES
          }

          Object@dfTissue <- tissue.df

          message("\n")
          message("[=========================================================]")
          message("[<<<<         ObtainTissueExpression END             >>>>>]")
          message("-----------------------------------------------------------")
          message("[You may want to plot the results using PlotTissueExpression next.]")

          return(Object)
        }
)
