library(BiocManager)
library(shiny)
library(shinydashboard)
library(shinyBS)
library(shinythemes)
library(DT)
library(ggplot2)
library(htmltools)
library(shinyjs)
library(ggnewscale)
library(ggtree)
library(GO.db)
library(reshape2)
library(clusterProfiler)
library(powerjoin)
library(imputeTS)
library(clustree)
library(cluster)
library(factoextra)
library(GGally)
library(shinyWidgets)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(dplyr)
library(limma)
library(stringr)
library(shinyalert)
library(jsonlite)
library(doParallel)
library(parallel)
library(httr)
library(utils)
library(readxl)
library(pbapply)
library(RColorBrewer)
library(patchwork)
library(gridExtra)
library(grid)
library(pheatmap)

library(ggplotify)
library(cowplot)
library(simplifyEnrichment)
library(GetoptLong)
library(ggwordcloud)
library(ComplexHeatmap)
library(colorRamp2)


#Load R scripts
source("source/PlotTree.R")
source("source/PlotGeneSets.R")


source("source/PathwayObject.R")
source("source/LoadGeneSets_shiny.R")
source("source/HighlightGeneSets_shiny.R")
source("source/TissueExpressionPerGeneSet_shiny.R")
source("source/PlotTissueExpression_shiny.R")
source("source/GenesPerGeneSet_shiny.R")
source("source/OptimalGeneSets_shiny.R")
source("source/internalFunctions.R")
source("source/CombineGeneSets.R")
source("source/ManageGeneSets.R")
source("source/ClusterGeneSets.R")



###########
#VERSION 15
###########

load("databases/hpoDatabase2.RData")
load("databases/mpDatabase2.RData")


underlineCellsCols <- function(rows, cols){
  stopifnot(length(rows) == length(cols))
  c(
    "function(row, data, num, index){",
    sprintf("  var rows = [%s];", paste0(rows-1, collapse = ",")),
    sprintf("  var cols = [%s];", paste0(cols, collapse = ",")),
    "  for(var i = 0; i < rows.length; ++i){",
    "      $('td:eq(' + cols[i] + ')', row)",
    "        .css({'text-decoration': 'underline'});",
    "  }",
    "}"
  )
}


getcharacter<-function(list, seperator){
  nseperator<-str_count(list,seperator)
  return(nseperator)
}


calculateORA<-function(object, cluster){
  genes<-as.vector(unlist(sapply(strsplit(object@Data[[1]][,"Molecules"][object@plot$aka2$Cluster==cluster],object@metadata$seperator[1]),unique)))
  ora<-enrichGO(gene = genes, OrgDb = object@metadata[1,"organism"], keyType = object@metadata[1,"structure"], ont="BP", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable=TRUE)
  #ora<-enrichGO(gene = genes[1:20], OrgDb = object@metadata[1,"organism"], keyType = object@metadata[1,"structure"], ont="BP", pvalueCutoff = 0.05, qvalueCutoff = 0.02, readable=TRUE, minGSSize=10, maxGSSize = 50)
  oratop5<-data.frame(ora[1:20,-c(8)], rep(paste("Cluster_",cluster, sep=""), nrow(ora[1:20,])))
  rm(ora, genes)
  oratop5[,c(5,6,7)]<-round(oratop5[,c(5,6,7)],3)
  colnames(oratop5)[9]<-"Cluster"
  return(oratop5)
}


topORA<-function(ora, cluster, top){
  return(head(ora[ora$Cluster==paste("Cluster_",cluster,sep=""),],top))
}

getterm<- function(goid){
  termid<-GOTERM[[goid]]
  if(is.null(termid)) {
    return(paste(goid,"NA","NA", sep="__"))
  }else{
    return(paste(goid,Term(termid), Definition(termid), sep="__"))
  }
}



# Define server function  
server <- function(input, output, session) {
  

  output$tissueIntro<-renderUI({tags$p("Press the button below if you want to calculate the tissue enrichment analysis (the tissues were selected from", HTML("<a href='https://gtexportal.org/home/' target='_blank'>here</a>"), "):")})

  ###########################
  ###########################
  #upload input files
  tablegroup<-reactiveValues(data=NULL)
  observeEvent(input$files, {
    tablegroup$data<-cbind(input$files[,1], as.character(c(1:dim(input$files)[1])))
    colnames(tablegroup$data)<-c("name","group")
    output$filesInfo<-renderDataTable(datatable(tablegroup$data, rownames=FALSE, editable=TRUE, options=list(pageLength=25, dom='t')))
  })
  
  ###########################
  ###########################
  #modify table and add group names
  observeEvent(input$filesInfo_cell_edit, {
    #get values
    info = input$filesInfo_cell_edit
    i = as.numeric(info$row)
    k = as.character(info$value)
    #write values to reactive
    tablegroup$data[i,2]<-k
  })
  
  ###########################
  ###########################
  #ANNOTATIONS GO
  annotations<-reactive({})
  
  dataInfo_previousvalue<-reactiveVal("none")
  oraInfo_previousvalue<-reactiveVal("none")
  genesInfo_previousvalue<-reactiveVal("none")
  #links to quickGO
  observe({
    #Data table
    if (!is.null(input$dataInfo_cell_clicked$col)){
      if (input$dataInfo_cell_clicked$col==0){
        #if(input$dataInfo_cell_clicked$value!=""){
        if((input$dataInfo_cell_clicked$value!="")&&(input$dataInfo_cell_clicked$value!=dataInfo_previousvalue())){
          dataInfo_previousvalue(input$dataInfo_cell_clicked$value)
          js$browseURL(paste("https://www.ebi.ac.uk/QuickGO/term/",input$dataInfo_cell_clicked$value, sep=""))
        }
      }
    }
    #ORA table
    if (!is.null(input$oraInfo_cell_clicked$col)){
      if (input$oraInfo_cell_clicked$col==0){
        #if(input$oraInfo_cell_clicked$value!=""){
        if((input$oraInfo_cell_clicked$value!="")&&(input$oraInfo_cell_clicked$value!=oraInfo_previousvalue())){
          oraInfo_previousvalue(input$oraInfo_cell_clicked$value)
          js$browseURL(paste("https://www.ebi.ac.uk/QuickGO/term/",input$oraInfo_cell_clicked$value, sep=""))
        }
      }
    }
    #Gene table
    if (!is.null(input$genesInfo_cell_clicked$col)){
      if (input$genesInfo_cell_clicked$col==0){
        #if(input$genesInfo_cell_clicked$value!=""){
        if((input$genesInfo_cell_clicked$value!="")&&(input$genesInfo_cell_clicked$value!=genesInfo_previousvalue())){
          genesInfo_previousvalue(input$genesInfo_cell_clicked$value)
          js$browseURL(paste("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",input$genesInfo_cell_clicked$value, sep=""))
        }
      }
    }
  })
  
  ###########################
  ###########################
  #CLEAR
  observeEvent(input$reset, {
    shinyjs::runjs("location.reload(true);")
  })
  
  ###########################
  ###########################
  #RUN EXAMPLE
  runexample<- reactiveVal(FALSE)
  
  observeEvent(input$setExample, {
    runexample(TRUE)
    updateRadioButtons(session,"source", selected = "Great")
    updateRadioButtons(session,"gene.id", selected = "SYMBOL")
    updateRadioButtons(session,"organism", selected = "org.Mm.eg.db")
    
    tablegroup$data<-cbind(c("MM10.GREAT.KO.uGvsMac.bed.tsv", "MM10.GREAT.WT.uGvsMac.bed.tsv"), c("KO", "WT"))
    colnames(tablegroup$data)<-c("name","group")
    output$filesInfo<-renderDataTable(datatable(tablegroup$data, rownames=FALSE, editable=TRUE, options=list(pageLength=25, dom='t')))
    
    click("printResults")
  })

  ###########################
  ###########################
  #CALCULATE RESULTS (including duplicates, clustering=hierarchival and number of cluster=5)
  treeplotoutput<-reactive({})
  heatmapoutput<-reactive({})
  tissueoutput<-reactive({})
  data.object<-reactive({})
  ora.object<-reactive({})
  ora.objectALL<-reactive({})
  
  observeEvent(input$printResults, {
    
    tocontinue=TRUE
    if(((is.null(input$source))||(is.null(input$files))||(is.null(input$gene.id))||(is.null(input$organism)))&&(!runexample())){
      shinyalert(title = "Missing inputs", text = "Please, all the inputs are mandatory.", type = "error")
    }else{
      showModal(modalDialog("Load gene sets (1/5)", footer=NULL))
      #1-if input=great load gene set with specifc parameters and perform manage.object
      if ((input$source=="Great") || (input$source=="GSEA") || (is.null(input$source))){
        
        if(is.null(input$source)){
          #run example
          load.object <- LoadGeneSets_shiny(file_location = files<-c("example/MM10.GREAT.KO.uGvsMac.bed.tsv",
                                                                     "example/MM10.GREAT.WT.uGvsMac.bed.tsv"),
                                            groupnames= c("KO","WT"),
                                            P.cutoff = 0.05,
                                            Mol.cutoff = 5,
                                            Source = "Great",
                                            Great.Background = FALSE, #TRUE
                                            type = "Canonical_Pathways",
                                            topranks = "",
                                            structure = "SYMBOL",
                                            Organism = "org.Mm.eg.db",
                                            seperator = ",") #,
        }else{
          load.object <- LoadGeneSets_shiny(file_location = input$files[,4],
                                            groupnames= tablegroup$data[,2],
                                            P.cutoff = 0.05,
                                            Mol.cutoff = 5,
                                            Source = input$source,
                                            Great.Background = FALSE, #TRUE
                                            type = "Canonical_Pathways",
                                            topranks = "",
                                            structure = input$gene.id,
                                            Organism = input$organism,
                                            seperator = ",") #,
        }
        
        
        
        if(is.character(load.object)){
          removeModal()
          #one of the group has 0 pathways
          shinyalert("Error", paste("Group ", load.object, " has 0 pathways passing cutoff (p<0.05)"), type = "error")
          tocontinue=FALSE
        }else{
          showModal(modalDialog("Manage gene sets (2/5).", footer=NULL))
          manage.object <- ManageGeneSets(Object = load.object,
                                          keep.type =c("Disease Ontology",
                                                       "GO Biological Process" ),
                                          exclude.type="")
        }
        
      
      }
      #1-if input=ipa load gene set with specifc parameters and don't perform manage.object
      else if ((input$source=="IPA")){
        
        load.object <- LoadGeneSets_shiny(file_location = input$files[,4],
                                              groupnames= tablegroup$data[,2],
                                              P.cutoff = 1.3,
                                              Mol.cutoff = 5,
                                              Source = input$source,
                                              type = "Canonical_Pathways",
                                              structure = input$gene.id,
                                              Organism = input$organism,
                                              seperator = ",")
        
        if(is.character(load.object)){
          removeModal()
          #one of the group has 0 pathways
          shinyalert("Error", paste("Group ", load.object, " has 0 pathways passing cutoff (p<0.05)"), type = "error")
          tocontinue=FALSE
        }else{
          manage.object<-load.object
        }
        
      }
      if(tocontinue){
  
        rm(load.object)
        
        #check seperator
        possibleSeparators<-c(",","/",";")
        nelements<-length(manage.object@Data)
        for (i in (1:nelements)){
          counts<-unlist(lapply(possibleSeparators,getcharacter,list=manage.object@Data[[i]]["Molecules"]))
          manage.object@metadata$seperator[i]<-possibleSeparators[which(counts==max(counts))]
        }
        if (length(unique(manage.object@metadata$seperator))!=1){
          removeModal()
          shinyalert("Error", "Please, use the same gene-seperator (, ; /) in all data ", type = "error")
          
        }else{
          
          showModal(modalDialog("Combine gene sets (3/5)", footer=NULL))
          #2-perform combine gene set
          combine.object <- CombineGeneSets(Object = manage.object)
          rm(manage.object)
          
          #2.1-calculate optimal number of cluster
          optimal.object<-OptimalGeneSets_2(object = combine.object, method = "silhouette", max_cluster= 24, cluster_method = "Hierarchical", main= "")
          nclust<-optimal.object$data
          optimalCluster<-as.numeric(nclust$clusters[which.max(nclust$y)])
          updateNumericInput(session, inputId="cluster", value=optimalCluster)
          
          showModal(modalDialog(paste("Cluster gene sets, optimal k=",optimalCluster," (4/5)", sep=""), footer=NULL))
          #3-perform cluster gene set
          cluster.object <- ClusterGeneSets(Object = combine.object,
                                            clusters = optimalCluster,
                                            method = "Hierarchical")
          
          rm(combine.object)
          data.object<<-reactive({cluster.object})
          #plot treeplot
          #treeplotoutput<<-reactive({PlotTreeShiny_2(Object = cluster.object, fontsize = 3, main= "", nodenames = FALSE)})
          treeplotoutput<<-reactive({PlotTree(Object = cluster.object, clusters = optimalCluster, doORA = F, nodenames = F, wordcloud = T)})
          
          if(is.null(treeplotoutput())){
            output$treeIntro<-renderText("Treeplot not performed because size greater than 200.")
          }
          output$treeplot<-renderPlot(treeplotoutput())
          #plot heatmap
          heatmapoutput<<-reactive({PlotGeneSets(cluster.object, doORA = F, wordcloud = T)})
          output$heatmap<-renderPlot(heatmapoutput())
          
          #download data and plots
          output$downloadData <- renderUI({
            downloadButton("download.data","Download data")
          })
          #download plots using different format
          output$downloadDataFormats <- renderUI({
            dropdownButton(inputId="download.plot.format",label="Download image", circle=FALSE, radioButtons("format", label="", choices = c("jpg","png","pdf"), selected = NULL), downloadButton("download.plots","Download plots"))
          })
          
          #put check box dependinf the number of cluster
          output$checkboxCluster <- renderUI({
            clusterchoices<-paste("Cluster_", names(cluster.object@plot$aka3$Cluster), sep="")
            checkboxGroupInput("optionCluster", label = "", choices=clusterchoices, selected = clusterchoices, inline=TRUE)
          })
          
          #change names of the cluster with Cluster_1
          cluster.object@Data[[1]]$cluster<-paste(rep("Cluster_",dim(cluster.object@Data[[1]]["cluster"])[1]), cluster.object@Data[[1]]$cluster, sep = "")
          if(!is.character(cluster.object@Data[[1]]$Ratio[1])){
            cluster.object@Data[[1]]$Ratio<-round(cluster.object@Data[[1]]$Ratio,3)
          }
          
          #calculate annotations
          if(is.null(annotations())){
            resultAnnotations<-as.data.frame(str_split_fixed(as.data.frame(unlist(lapply(cluster.object@Data[[1]][,"Pathways"],getterm)))[,1],"__",3))
            colnames(resultAnnotations)<-c("ID", "Term", "Definition")
            resultAnnotationsFinal<-data.frame(resultAnnotations[,c(1,2)], cluster.object@Data[[1]][c("Type","Groups","cluster","Pval","Ratio")], resultAnnotations[,3])
            colnames(resultAnnotationsFinal)[c(4,5,8)]<-c("Group", "Cluster", "Definition")
            annotations<<-reactive({resultAnnotationsFinal})
          }
          
          #4-ORA per cluster using enrichGO of clusterprofiler
          showModal(modalDialog("ORA per cluster (5/5)", footer=NULL))
          for (i in (1:max(as.numeric(cluster.object@plot$aka2$Cluster)))){
            #print(paste("cluster",i))
            genes<-as.vector(unlist(sapply(strsplit(cluster.object@Data[[1]][,"Molecules"][cluster.object@plot$aka2$Cluster==i],cluster.object@metadata$seperator[1]),unique)))
            #count duplication per gene
            genesTimes<-as.data.frame(genes) %>% group_by_all() %>% count
            totalPathways<-sum(data.object()@plot$aka2$Cluster==i)
            genesTimes$n<-round(genesTimes$n/totalPathways,3)
            colnames(genesTimes)[1]<-"gene"
            if (i==1){
              genesFinal<-genesTimes
              colnames(genesFinal)[i+1]<- paste("Cluster_",i, " (n=", totalPathways,")", sep="")
            }else{
              genesFinal<-power_full_join(genesFinal, genesTimes, by="gene")
              colnames(genesFinal)[i+1]<- paste("Cluster_",i, " (n=", totalPathways,")", sep="")
            }
          }
          #calulate ORA in parallel
          oraResults<-do.call("rbind", lapply(1:max(as.numeric(cluster.object@plot$aka2$Cluster)), calculateORA, object=cluster.object))
          #oraResults<-calculateORA(object=cluster.object, cluster=1)
          ora.objectALL<<-reactive({oraResults})
          oraResultsTop<-do.call("rbind", lapply(1:max(as.numeric(cluster.object@plot$aka2$Cluster)), topORA, ora=oraResults, top=as.integer(input$top)))
          ora.object<<-reactive({oraResultsTop})
          #output$oraInfo<-renderDataTable(datatable(as.data.frame(oraResultsTop), rownames = FALSE, options = list(searching=TRUE, pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0))), columnDefs = list(list(targets = c(7), visible = FALSE))))
          #                                %>% formatStyle(colnames(oraResultsTop)[1], color = "blue"))
          output$oraInfo<-renderDataTable(datatable(as.data.frame(oraResultsTop), rownames = FALSE, options = list(searching=TRUE, pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                          %>% formatStyle(colnames(oraResultsTop)[1], color = "blue"))
          
          genesFinal<-na_replace(genesFinal,0)
          output$genesInfo <-renderDataTable(datatable(genesFinal, rownames = FALSE, filter="top", options = list(dom='lrtip', pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                             %>% formatStyle(colnames(genesFinal)[1], color = "blue"))
          
          rm(genesTimes,genes)
          removeModal()
          
          #hide tissue results
          hideTab(inputId = "tabs", target = "Tissue")
          #if(input$organism=="org.Mm.eg.db"){
          #  #for mouse data hide tissue enrichment analysis
          #  hideTab(inputId = "tabs2", target = "Tissue enrichment")
          #}
          
          shinyjs::toggle("main1")
          shinyjs::toggle("main2")
        }
      }

    }

    
  }) #end of printResults
  
  
  ###########################
  ###########################
  #DATABASE OPTIONS
  output$databases <- renderUI({
    if(!is.null(data.object())){
      #print(data.object()@metadata[1,"organism"])
      if (data.object()@metadata[1,"organism"]=="org.Hs.eg.db") {
        updateSelectizeInput(session = session, inputId = "database", choices=c("","Human Phenotype Ontology (HPO)"), server=TRUE)
        selectizeInput("database", "Human databases:", choices=NULL, selected = NULL, multiple=FALSE)  
      }else if (data.object()@metadata[1,"organism"]=="org.Mm.eg.db") {
        updateSelectizeInput(session = session, inputId = "database", choices=c("","Mammalian Phenotype (MP)"), server=TRUE)
        selectizeInput("database", "Mouse databases:", choices=NULL, selected = NULL, multiple=FALSE)  
      }
    }
  })
  
  output$databaseOptions <- renderUI({
    
    if(!is.null(input$database)){ 
      if (input$database=="Human Phenotype Ontology (HPO)") {
        updateSelectizeInput(session = session, inputId = "hpo", choices=c("",unique(hpoDatabase$hpo_id)), server=TRUE)
        selectizeInput("hpo", HTML("Choose HPO <a href='https://hpo.jax.org/app/browse/term/HP:0000118' target='_blank'>(info)</a>:"), choices=NULL, selected = NULL, multiple=FALSE)  
      }else if (input$database=="Mammalian Phenotype (MP)") {
        updateSelectizeInput(session = session, inputId = "mp", choices=c("",unique(mpDatabase$mp_id)), server=TRUE)
        selectizeInput("mp", HTML("Choose MP <a href='https://www.informatics.jax.org/vocab/mp_ontology' target='_blank'>(info)</a>:"), choices=NULL, selected = NULL, multiple=FALSE)  
      }
    }
  })
  
  
  observe({
    if((input$hpo!="")&&(!is.null(input$hpo))){
      output$databaseOptionTitle<-renderText({hpoDatabase$hpo_name[hpoDatabase$hpo_id==input$hpo][1]})
      output$databaseGenes<-renderText({hpoDatabase$gene_symbol[hpoDatabase$hpo_id==input$hpo]})
    }
    if((input$mp!="")&&(!is.null(input$mp))){
      output$databaseOptionTitle<-renderText({mpDatabase$mp_term[mpDatabase$mp_id==input$mp][1]})
      output$databaseGenes<-renderText({mpDatabase$gene_symbol[mpDatabase$mp_id==input$mp]})
    }
  })
  
  observeEvent(input$calculateHighlight, {
    if((input$database=="")|((input$database=="Mammalian Phenotype (MP)")&&(input$mp==""))|((input$database=="Human Phenotype Ontology (HPO)")&&(input$hpo==""))){
      shinyalert(title = "Wrong database", text = "Please, select a database and a gene set.", type = "error")
    }else{
      showModal(modalDialog("Calculate highlight genes", footer=NULL))
      if((input$hpo!="")&&(!is.null(input$hpo))){
        genes<-hpoDatabase$gene_symbol[hpoDatabase$hpo_id==input$hpo]
      }
      if((input$mp!="")&&(!is.null(input$mp))){
        genes<-mpDatabase$gene_symbol[mpDatabase$mp_id==input$mp]
      }
      highlight.object<-HighlightGeneSets_shiny(data.object(), highligt.genes = genes, name="")
      data.object<<-reactive({highlight.object})
      resultHighlight<-(t(as.data.frame(str_split(unique(paste(data.object()@Data[[1]]$cluster, round(data.object()@Data[[1]]$Highlight.mean, 3)))," "))))
      colnames(resultHighlight)<-c("Cluster", "Mean")
      resultHighlight<-resultHighlight[order(resultHighlight[,1]),]
      output$highlightInfo<-renderDataTable(datatable(resultHighlight, rownames = FALSE, options = list(scrollY=400, dom='t')))
      removeModal()
    }
  })
  
  ###########################
  ###########################
  #FILTER TOP ORA
  observeEvent(input$top, {
    if(!is.null(data.object())){
      oraResultsTop<-do.call("rbind", lapply(1:max(as.numeric(data.object()@plot$aka2$Cluster)), topORA, ora=ora.objectALL(), top=as.integer(input$top)))
      ora.object<<-reactive({oraResultsTop})
      oraFiltered<-as.data.frame(ora.object())
      oraFiltered<-oraFiltered[oraFiltered$Cluster %in% input$optionCluster,]
      #output$oraInfo<-renderDataTable(datatable(oraFiltered, rownames = FALSE, options = list(searching=TRUE, pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0))), columnDefs = list(list(targets = c(7), visible = FALSE))))
      #                                %>% formatStyle(colnames(oraFiltered)[1], color = "blue"))
      output$oraInfo<-renderDataTable(datatable(oraFiltered, rownames = FALSE, options = list(searching=TRUE, pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                      %>% formatStyle(colnames(oraFiltered)[1], color = "blue"))
      
    }
    
  })
  
  ###########################
  ###########################
  #FILTER BY CLUSTER -> ORA and table of GO
  observe({

    if((is.null(input$optionCluster)) && (!is.null(data.object()))){
      #if any option is selected print NULL
      tableFiltered<-annotations()
      tableFiltered<-tableFiltered[tableFiltered$Cluster %in% c("Cluster_X"),]
      output$dataInfo <-renderDataTable(datatable(tableFiltered, rownames = FALSE, options = list(autoWidth=TRUE, scrollX=TRUE, scrollY=400, searching=TRUE, pageLength=25, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                        %>% formatStyle(colnames(tableFiltered)[1], color = "blue"))


      oraFiltered<-as.data.frame(ora.object())
      oraFiltered<-oraFiltered[oraFiltered$Cluster %in% c("Cluster_X"),]
      #output$oraInfo<-renderDataTable(datatable(oraFiltered, rownames = FALSE, options = list(searching=TRUE, pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0))), columnDefs = list(list(targets = c(7), visible = FALSE))))
      #                                  %>% formatStyle(colnames(oraFiltered)[1], color = "blue"))
      output$oraInfo<-renderDataTable(datatable(oraFiltered, rownames = FALSE, options = list(searching=TRUE, pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                      %>% formatStyle(colnames(oraFiltered)[1], color = "blue"))

    }
  })
  
  observeEvent(input$optionCluster,{
    
    tableFiltered<-annotations()
    tableFiltered<-tableFiltered[tableFiltered$Cluster %in% input$optionCluster,]
    output$dataInfo <-renderDataTable(datatable(tableFiltered, rownames = FALSE, options = list(autoWidth=TRUE, scrollX=TRUE, scrollY=400, searching=TRUE, pageLength=25, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                        %>% formatStyle(colnames(tableFiltered)[1], color = "blue"))
    
    oraFiltered<-as.data.frame(ora.object())
    oraFiltered<-oraFiltered[oraFiltered$Cluster %in% input$optionCluster,]
    #output$oraInfo<-renderDataTable(datatable(oraFiltered, rownames = FALSE, options = list(searching=TRUE, pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0))), columnDefs = list(list(targets = c(7), visible = FALSE))))
    #                                    %>% formatStyle(colnames(oraFiltered)[1], color = "blue"))
    output$oraInfo<-renderDataTable(datatable(oraFiltered, rownames = FALSE, options = list(searching=TRUE, pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                    %>% formatStyle(colnames(oraFiltered)[1], color = "blue"))
    
    
  })
  

  ###########################
  ###########################
  #RECALCULATING CLUSTERING
  observeEvent(input$recalculateClustering, {
    showModal(modalDialog("Recalculating...", footer=NULL))
    new.object<-ClusterGeneSets(Object = data.object(),
                    clusters = input$cluster,
                    method = "Hierarchical")
    
    #change names of the cluster with Cluster_1
    new.object@Data[[1]]$cluster<-paste(rep("Cluster_",dim(new.object@Data[[1]]["cluster"])[1]), new.object@Data[[1]]$cluster, sep = "")
    
    data.object<<-reactive({new.object})
    #plot new treeplot
    #treeplotoutput<<-reactive({PlotTreeShiny_2(Object = data.object(), fontsize = 3, main= "", nodenames = FALSE)})
    treeplotoutput<<-reactive({PlotTree(Object = data.object(), clusters = input$cluster, doORA = F, nodenames = F, wordcloud = T)})
    if(is.null(treeplotoutput)){
      output$treeIntro<-renderText("Treeplot not performed because size greater than 200.")
    }
    output$treeplot<-renderPlot(treeplotoutput())
    #plot new heatmap
    #heatmapoutput<<-reactive({PlotGeneSets(Object = data.object(), fontsize = 5, main= "", legend=TRUE, annotation.mol=FALSE, RR.max=50)})
    heatmapoutput<<-reactive({PlotGeneSets(Object = data.object(), doORA = F, wordcloud = T)})
    output$heatmap<-renderPlot(heatmapoutput())
    #show new data results
    annotationsUpdate<-annotations()
    annotationsUpdate$Cluster<-data.object()@Data[[1]]$cluster
    annotations<<-reactive({annotationsUpdate})
    
    
    #recalculate ORA per cluster
    for (i in (1:max(as.numeric(data.object()@plot$aka2$Cluster)))){
      genes<-as.vector(unlist(sapply(strsplit(data.object()@Data[[1]][,"Molecules"][data.object()@plot$aka2$Cluster==i],data.object()@metadata$seperator[1]),unique)))
      #count duplication per gene
      genesTimes<-as.data.frame(genes) %>% group_by_all() %>% count
      totalPathways<-sum(data.object()@plot$aka2$Cluster==i)
      genesTimes$n<-round(genesTimes$n/totalPathways,3)
      colnames(genesTimes)[1]<-"gene"
      if (i==1){
        #auxfinal<-aux
        genesFinal<-genesTimes
        colnames(genesFinal)[i+1]<- paste("Cluster_",i, " (n=", totalPathways,")", sep="")
      }else{
        #auxfinal<-rbind(auxfinal, aux)
        genesFinal<-power_full_join(genesFinal, genesTimes, by="gene")
        colnames(genesFinal)[i+1]<- paste("Cluster_",i, " (n=", totalPathways,")", sep="")
      }
    }
    
    oraResults<-do.call("rbind", lapply(1:max(as.numeric(data.object()@plot$aka2$Cluster)), calculateORA, object=data.object()))
    #oraResults[,c(5,6,7)]<-round(oraResults[,c(5,6,7)],3)
    ora.objectALL<<-reactive({oraResults})
    
    oraResultsTop<-do.call("rbind", lapply(1:max(as.numeric(data.object()@plot$aka2$Cluster)), topORA, ora=oraResults, top=as.integer(input$top)))
    ora.object<<-reactive({oraResultsTop})

    genesFinal<-na_replace(genesFinal,0)
    output$genesInfo <-renderDataTable(datatable(genesFinal, rownames = FALSE, filter="top", options = list(pageLength=10, rowCallback = JS(underlineCellsCols(c(1),c(0)))))
                                        %>% formatStyle(colnames(genesFinal)[1], color = "blue"))
    #put check box of cluster
    output$checkboxCluster <- renderUI({
      clusterchoices<-paste("Cluster_", names(data.object()@plot$aka3$Cluster), sep="")
      checkboxGroupInput("optionCluster", label = "", choices=clusterchoices, selected = clusterchoices, inline=TRUE)
    })
    
    #reset tissue results
    shinyjs::show("calculateTissueEnrichment")
    hideTab(inputId = "tabs", target = "Tissue")
    output$tissueOutput<-renderDataTable(NULL)
    output$tissueIntro<-renderUI({tags$p("Press the button below if you want to calculate the tissue enrichment analysis (the tissues were selected from", HTML("<a href='https://gtexportal.org/home/' target='_blank'>here</a>"), "):")})
    removeModal()
    
  })
  
  
  ###########################
  #CALCULATE TISSUE ENRICHMENT CLUSTERING
  observeEvent(input$calculateTissueEnrichment, {
    showModal(modalDialog("Calculate tissue enrichment", footer=NULL))
    load("databases/TissueLocalDatabase15.RData")
    load("databases/dic.rda")
    shinyjs::hide("calculateTissueEnrichment")
    tissue.object <- TissueExpressionPerGeneSet(data.object(), localDatabase = localDatabase15, dic=dic)
    tissue.object@dfTissue<-round(tissue.object@dfTissue,3)
    tissue.object@dfTissue<-na_replace(tissue.object@dfTissue,0)
    data.object<<-reactive({tissue.object})

    tissueData<-data.frame(rownames(tissue.object@dfTissue), tissue.object@dfTissue)
    colnames(tissueData)[1]<-"Tissue"
    output$tissueOutput<-renderDataTable(datatable(tissueData, rownames = FALSE, options = list(searching=TRUE, pageLength=10)))
    shinyjs::hide("calculateTissueEnrichment")

    output$tissueIntro<-renderUI({tags$p("Tissue enrichment results:")})
    showTab(inputId = "tabs", target = "Tissue")
    #plot tissue plots
    tissueoutput<<-reactive({PlotTissueExpression(tissue.object)})
    output$tissuePlot<-renderPlot(tissueoutput())

    rm(localDatabase15)
    rm(dic)
    removeModal()
    
  })
  
  
  
  ###########################
  ###########################
  #DOWNLOAD
  output$downloadTemplate <- downloadHandler(
    filename = "template.xls",
    content = function(file) {
      file.copy("source/template.xls",file)
    }
  )
  
  output$downloadUserGuide <- downloadHandler(
    filename = "UserGuide.pdf",
    content = function(file) {
      file.copy("source/UserGuide.pdf",file)
    }
  )
  
  output$download.plots<-downloadHandler(
    file = function() {
      paste('plots-', Sys.Date(), '.zip', sep = '')
    },
    content = function(file) {

      if((is.null(treeplotoutput())) & (is.null(tissueoutput()))){
        filenames<-c(paste("heatmap_plot.",input$format, sep=""))
        # Save the heatmap plot
        ggsave(plot = as.ggplot(heatmapoutput()), filename = filenames[1], width = 30, unit = "cm", device = input$format)

      }else if ((is.null(treeplotoutput())) & (!is.null(tissueoutput()))){
        filenames<-c(paste("heatmap_plot.",input$format, sep=""), paste("tissue_plot.",input$format, sep=""))
        # Save the heatmap plot
        ggsave(plot = as.ggplot(heatmapoutput()), filename = filenames[1], width = 30, unit = "cm", device = input$format)
        # Save the tissue plot
        ggsave(plot = tissueoutput(), filename = filenames[2], width = 30, unit = "cm", device = input$format)

      }else if ((!is.null(treeplotoutput())) & (is.null(tissueoutput()))){
        filenames<-c(paste("heatmap_plot.",input$format, sep=""), paste("tree_plot.",input$format, sep=""))
        # Save the heatmap plot
        ggsave(plot = as.ggplot(heatmapoutput()), filename = filenames[1], width = 30, unit = "cm", device = input$format)
        # Save the tree plot
        ggsave(plot = treeplotoutput(), filename = filenames[2], width = 30, unit = "cm", device = input$format)

      }else{
        filenames<-c(paste("heatmap_plot.",input$format, sep=""), paste("tree_plot.",input$format, sep=""),paste("tissue_plot.",input$format, sep=""))
        # Save the heatmap plot
        ggsave(plot = as.ggplot(heatmapoutput()), filename = filenames[1], width = 30, unit = "cm", device = input$format)
        # Save the tree plot
        ggsave(plot = treeplotoutput(), filename = filenames[2], width = 30, unit = "cm", device = input$format)
        # Save the tissue plot
        ggsave(plot = tissueoutput(), filename = filenames[3], width = 30, unit = "cm", device = input$format)
      }

      # Create the zip file
      zip::zipr(zipfile = file, files = filenames)
    }
  )
  
  
  output$download.data<-downloadHandler(
    filename = "data.csv",
    content = function(file) {
      if (input$source=="IPA"){
        write.csv(data.object()@Data[[1]][,c(1:9)], file, row.names=FALSE) # add parentheses to data arg if reactive
      }else{
        write.csv(data.object()@Data[[1]][,c(1:7,10)], file, row.names=FALSE) # add parentheses to data arg if reactive
      }
    }
  )

} # server



