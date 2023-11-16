### Functions --------------------------------------------------------------------
# These functions are not meant to be invoked directly by the user.
# See the PlotPathwayCluster method instead.
obtainDFmetadata <- function(df, mat_cor) {
  ngroups <- unique(df$Groups)
  info <- vector("list", length(ngroups))

  for (i in 1:length(ngroups)) {
    info[[i]] <- df[which(df$Groups==ngroups[i]), "Pathways"]
  }

  df_metadata = as.data.frame(matrix(0, nrow = nrow(mat_cor), ncol = length(ngroups)))
  rownames(df_metadata) = colnames(mat_cor)
  colnames(df_metadata) = ngroups

  for (i in 1:nrow(mat_cor)) {
    for (j in 1:length(ngroups)) {
      if (rownames(mat_cor)[i] %in% info[[j]]) {
        df_metadata[i, j] <- 1
      }
    }
  }

  return(df_metadata)
}


obtainDefCluster <- function(mat_using){

  candidates_clustering <- c()
  definitive_clusters <- list()
  j <- 1
  num_sim = 0.6
  for (i in 1:(ncol(mat_using)-1)) {
    if (i < j-1 | names(mat_using[i,])[j] %in% candidates_clustering) {  ## Avoid repeating same cluster
      next
    }

    if (i == ncol(mat_using)-1) {
      if (mat_using[i,i+1] >= num_sim) {
        candidates_clustering <- c(names(mat_using[i,])[i], names(mat_using[i,])[i+1])
        definitive_clusters <- append(definitive_clusters, list(candidates_clustering))
      }
    } else {

      if (mat_using[i,i+1] >= num_sim) {
        candidates_clustering <- c(names(mat_using[i,])[i])

        for (j in (i+1):(ncol(mat_using))) {

          if (j==ncol(mat_using)) {
            candidates_clustering <- c(candidates_clustering, names(mat_using[i,])[j])
            definitive_clusters <- append(definitive_clusters, list(candidates_clustering))
            break
          }

          if(mat_using[i,j] >= num_sim) {
            candidates_clustering <- c(unique(candidates_clustering), names(mat_using[i,])[j])
            if (j == (ncol(mat_using)-2)) {
              append(definitive_clusters, list(candidates_clustering))
            }

          } else if (mat_using[i+1,j] >= num_sim){ # The lower this value is, the larger group of clusters is obtained
            candidates_clustering <- c(unique(candidates_clustering), names(mat_using[i,])[j])
            if (j == (ncol(mat_using)-2)) {
              append(definitive_clusters, list(candidates_clustering))
            }

          } else {
            definitive_clusters <- append(definitive_clusters, list(candidates_clustering))
            break
          }
        }
      }
    }
  }

  return(definitive_clusters)
}



# This function has been adapted from simplifyEnrichment package (Gu, Z. & Hubschmann, D. 2022).
#' @import simplifyEnrichment
#' @import GetoptLong
keywordEnrichment <- function(term_id, tdm, min_bg = 5, min_term = 2) {
  GO_EXCLUDE_WORDS <- c("via", "protein", "factor", "side", "type", "specific", "biological", "process")

  tdm2 <- tdm[slam::row_sums(tdm) >= 5, ]
  l <- colnames(tdm2) %in% term_id

  n <- nrow(tdm2)
  n_term <- numeric(n)
  n_bg <- numeric(n)
  p <- numeric(n)
  for(i in seq_len(n)) {
    if(interactive() && se_opt$verbose) {
      if(i %% 100 == 0 || i == n) {
        message(strrep("\r", 100), appendLF = FALSE)
        message(GetoptLong::qq("performing keyword enrichment, @{i}/@{n}..."), appendLF = FALSE)
      }
    }
    v <- as.vector(tdm2[i, ])
    s11 <- sum(v & l)
    if(s11 < 2) {
      next
    }
    s12 <- sum(!v & l)
    s21 <- sum(v & !l)
    s22 <- sum(!v & !l)

    n_term[i] <- s11
    n_bg[i] <- s11 + s21

    p[i] <- fisher.test(cbind(c(s11, s21), c(s12, s22)), alternative = "greater")$p.value

  }
  if(interactive() && se_opt$verbose) {
    message("")
  }

  df <- data.frame(keyword = rownames(tdm2), n_term = n_term, n_bg = n_bg, p = p)
  df <- df[df$n_term >= 2, , drop = FALSE]
  df$padj <- p.adjust(df$p)
  df[order(df$padj, df$p), , drop = FALSE]

  df <- df[which(!df$keyword %in% GO_EXCLUDE_WORDS),]

  return(df)
}


#' @import AnnotationDbi
#' @import GO.db
obtainGOdescription <- function(go_id) {
  go_descriptions <- AnnotationDbi::select(GO.db, keys = go_id, columns = "TERM")
  return(go_descriptions)
}


#' @import simplifyEnrichment
#' @import ggwordcloud
wordcloud_generation <- function(go_id_list, min_stat=5) {
  env_tdm_GO <- readRDS(system.file("extdata", "tdm_GO.rds", package = "simplifyEnrichment"))
  df <- keywordEnrichment(go_id_list, env_tdm_GO)
  df <- df[df$p <= min_stat, , drop = FALSE]
  df <- data.frame(df[, 1], -log10(df$p))
  if (nrow(df) > 10) {
    df <- df[1:10,]
  }

  wordplot <- ggplot(df, aes(label = df[,1], size = df[,2])) +
                             # color = factor(sample.int(10, nrow(df), replace = TRUE)))) + #to add colour
                    geom_text_wordcloud(shape="square") +
                    theme_minimal() +
                    labs(title = "")+
                    theme(
                      plot.title = element_blank(),
                      legend.title = element_blank(),
                      legend.background = element_blank(),
                      legend.box.spacing = unit(0, "mm")
                    )

  return(wordplot)
}

obtainOrg <- function(Object) {
  mm <- c("MUS MUSCULUS" ,"MM", "ORG.MM.EG.DB")
  hs <- c("HOMO SAPIENS", "HS", "ORG.HS.EG.DB")
  if (toupper(Object@metadata$organism[1]) %in% hs) {
    usingOrg <- org.Hs.eg.db
  } else if (toupper(Object@metadata$organism[1]) %in% mm) {
    usingOrg <- org.Mm.eg.db
  } else {
    message("Not recognised organism. Will use Homo sapiens data base as default.")
    usingOrg <- org.Hs.eg.db
  }
  return(usingOrg)
}

#' @import AnnotationDbi
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import clusterProfiler
#' @import stringr
#' @export
ORAperCluster <- function(Object, doORA, clusters, clus, completeRes = FALSE) {
  options <- c("SYMBOL", "ENTREZID", "ENSEMBLID", "ENTREZ", "ENSEMBL")
  entrezOptions <- c("ENTREZID", "ENTREZ")
  ensemblOptions <- c("ENSEMBL", "ENSEMBLID")

  # internal check of parameters
  if (!toupper(Object@metadata$structure[1]) %in% options | doORA == FALSE) {
    if (!toupper(Object@metadata$structure[1]) %in% options) {
      message("Genes are not included as genes Symbol, EntrezID or ENSEMBLID.\nORA of the genes involved in the cluster pathway will NOT be performed.")
    }
    keywords <- rep(".", times=clusters)
    return(keywords)
  }

  #organism
  usingOrg <- obtainOrg(Object)

  #obtain genes EntrezID for ORA
  info_df <- data.frame(Object@Data[[1]]$Pathways, Object@Data[[1]]$Molecules)
  keywords <- vector()
  message(paste0("Performing ORA for ", clusters, " clusters and choosing the most overrepresented pathway..."))
  for (i in 1:clusters) {
    paths_of_clust <- names(clus[clus == i])
    genes2check <- unique(unlist(str_split(paste0(info_df[which(info_df[,1] %in% paths_of_clust),2], collapse = ","), Object@metadata$seperator[1])))

    if (toupper(Object@metadata$structure[1]) == "SYMBOL") {
      if (table(genes2check==toupper(genes2check))[1] < length(genes2check)/2 | table(genes2check==tolower(genes2check))[1] < length(genes2check)/2 | all(genes2check==toupper(genes2check)) | all(genes2check==tolower(genes2check))){ # To correct capital letters of gene symbol
        genes2check <- tolower(genes2check)
        genes2check <- paste(toupper(substr(tolower(genes2check), 1, 1)), substr(tolower(genes2check), 2, nchar(tolower(genes2check))), sep = "")
      }
      entrezID <- AnnotationDbi::select(usingOrg, keys=genes2check, columns='ENTREZID', keytype='SYMBOL')
      # entrezID <- AnnotationDbi::select(usingOrg, keys=toupper(genes2check), columns='ENTREZID', keytype='SYMBOL')
      entrezID <- entrezID$ENTREZID
    } else if (toupper(Object@metadata$structure[1]) %in% entrezOptions) {
      entrezID <- genes2check
    } else if (toupper(Object@metadata$structure[1]) %in% ensemblOptions) {
      entrezID <- AnnotationDbi::select(usingOrg, keys=genes2check, columns='ENTREZID', keytype='ENSEMBL')
      entrezID <- entrezID$ENTREZID
    }

    ora <- clusterProfiler::enrichGO(na.omit(entrezID), OrgDb=usingOrg, keyType="ENTREZID",
                                     ont="BP", pvalueCutoff=1, pAdjustMethod="BH", qvalueCutoff=1)

    if (completeRes == TRUE) {
      keywords[i] <- ora
    } else {

      if (nchar(ora@result$Description[1]) > 50) {
        keywords[i] <- paste0(substr(ora@result$Description[1], 1, 40), "... ", ora@result$ID[1])
      } else {
        keywords[i] <- ora@result$Description[1]
      }
    }
  }
  return(keywords)
}

#' @import GO.db
getterm <- function(goid){
  termid <- GO.db::GOTERM[[goid]]
  if(is.null(termid)) {
    return("NA")
  } else {
    return(AnnotationDbi::Term(termid))
  }
}

convertNA <- function(x) {
  na_count <- sum(x == "NA")
  ifelse(x == "NA", paste0("NA.", seq_len(na_count)), x)
}

checkGO <- function(Object) {
  if (substr(Object@Data[[1]]$Pathways[1], 1, 3) == "GO:") {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

parallelization.message <- function(threads) {
  message(paste0("Parallelizing processes, ", threads, " cores have been selected. "), sep="")

  if (threads > detectCores()) {
    message("More cores have been selected than detected.")
    threads <- detectCores()
    message(paste0("The maximum number of cores will be used:", threads, " cores."))
  }

  return(threads)
}

scale_values <- function(x, ...) {(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}

#' @import seriation
optimalDist <- function(mat_sym, all_methods=FALSE, iterations=3, use_method=NULL, min_length=NULL) {
  if (is.null(use_method))
  {
    if (all_methods) {
      message("Selected all methods... Running 32 mehtods this could take some time... Check INCLUDE SOMETHING for more information.")
    } else {
      message("Running 26 distance methods to choose the optimal one. Check INCLUDE SOMETHING for more information.")
    }

    if (ncol(mat_sym) <= 300 | all_methods == TRUE) { # For big datasets the most time consuming methods gonna be avoid (as default)
      methods_list <- list_seriation_methods('dist')[c(-2, -3, -4, -5, -16, -17, -33, -34, -36, -42)] #for small datasets
      methods_repeat <- list_seriation_methods('dist')[c(1, 15, 16, 17, 18, 19, 20, 21, 22, 28, 29, 30, 31, 39, 40)]

    } else {
      methods_list <- list_seriation_methods('dist')[c(-1, -2, -3, -4, -5, -16, -17, -21, -30, -31, -33, -34, -36, -39, -40, -42)] #for big datasets
      methods_repeat <- list_seriation_methods('dist')[c(15, 16, 17, 18, 19, 39, 40)]
    }

    if (is.null(min_length))
    {
      min_length <- selectNpathways(mat_sym)
    } else {
      min_length <- min_length
    }
    res_list <- list()
    res_order <- list()
    for (meth_sel in methods_list) {

      if (meth_sel %in% methods_repeat)
      {
        best_res <- 0
        for (z in 1:iterations)
        {
          i <- dist2order(mat_sym, meth_sel)
          i_def <- i
          res <- obtainDefCluster(mat_sym[i,i])
          go_id_list <- res[which(lapply(res, function(x) length(x)) >= min_length)]
          score <- mean(unlist(lapply(go_id_list, length)))

          if (mean(score) > mean(best_res) & length(best_res) > 1)
          {
            best_res <- score
            i_def <- i
          }
        }
        cluster_num <- best_res
      } else {
        i <- dist2order(mat_sym, meth_sel)
        res <- obtainDefCluster(mat_sym[i,i])
        go_id_list <- res[which(lapply(res, function(x) length(x)) >= min_length)]
        cluster_num <- unlist(lapply(go_id_list, length))
        i_def <- i
      }

      res_list <- append(res_list, list(c(length(cluster_num), max(cluster_num), mean(cluster_num), median(cluster_num))))
      res_order <- append(res_order, list(i_def))

    }
    names(res_list) <- methods_list
    names(res_order) <- methods_list

    valid_methods <- names(which(unlist(lapply(res_list, function(x) x[1])>1))) # discard those ones with 1 or 0 clusters
    for (z in 1:length(res_list))
    {
      best_method <- sort(unlist(lapply(res_list, function(x) x[3])), decreasing = T)[z] # selected by mean of length of all the clusters
      if (names(best_method) %in% valid_methods)
      {
        break()
      }
    }
    order_method <- res_order[names(best_method)]
    message(paste0("The best estimated method with ", iterations, " iterations is: ", names(best_method), "\n Using ", names(best_method), " to order the correlation matrix..."))

    return(list(order_method[[1]], names(best_method), iterations))

  } else {
    if (use_method %in% list_seriation_methods('dist'))
    {
      message(paste0("Selected ", use_method, " method and ", iterations, " iterations", sep=""))
      methods_repeat <- list_seriation_methods('dist')[c(1, 15, 16, 17, 18, 19, 20, 21, 22, 28, 29, 30, 31, 39, 40)]
      if (!use_method %in% methods_repeat) {
        iterations <- 1
      }
      best_score <- 0

      if (is.null(min_length))
      {
        min_length <- selectNpathways(mat_sym)
      } else {
        min_length <- min_length
      }
      res_list <- list()
      res_order <- list()
      for (y in 1:iterations)
      {
        i <- dist2order(mat_sym, method = use_method)
        i_def <- i
        res <- obtainDefCluster(mat_sym[i,i])
        go_id_list <- res[which(lapply(res, function(x) length(x)) >= min_length)]
        score <- mean(unlist(lapply(go_id_list, length)))

        if (mean(score) > mean(best_score) & length(best_score) > 1)
        {
          best_score <- score
          i_def <- i
        }
      }
      return(list(i_def, use_method, iterations))

    } else {
      message("The selected method is not correct or it is not among the possible methods. Please check INCLUDE SOMETHING for more information.")
      stop()
    }
  }
}


#' @import seriation
dist2order <- function(mat, method) {
  d_corr <- as.dist(1 - mat)
  s <- seriation::seriate(d_corr, method = method)
  order <- seriation::get_order(s)
  return(order)
}

selectNpathways <- function(mat_sym){
  if (ncol(mat_sym) <= 50)
  {
    min_length <- 4
  } else if (ncol(mat_sym) >= 300)
  {
    min_length <- 10
  } else {
    min_length <- 6
  }

  return(min_length)
}

