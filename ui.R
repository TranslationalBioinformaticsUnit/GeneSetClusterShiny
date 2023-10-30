library(shiny)
library(shinythemes)
library(DT)
library(shinydashboard)
library(shinyWidgets)
library(shinyBS)
library(shinyjs)
library(ggplot2)
library(gridExtra)
library(htmltools)



###########
#VERSION 13
###########

#data
options(shiny.maxRequestSize=50*1024^3)

js <- "$(document).ready(function(){
    var $treeplot = $('#hello li > a[data-value=treeplot_val]').parent(); 
    $treeplot.removeClass('active').addClass('hide');
    $('#printPlot').on('click', function(){
      $treeplot.removeClass('hide');
    });
  });
"
js_code <- "
shinyjs.browseURL = function(url) {
  window.open(url,'_blank');
}
"

# Define UI
ui <- fluidPage(theme = shinytheme("flatly"),
                tags$head(
                  tags$script(HTML(js))
                ),
                useShinyjs(),
                extendShinyjs(text = js_code, functions = 'browseURL'),
                navbarPage(
                  #"GeneSetCluster",
                  title=div("", img(src = "GeneSetClusterLogo_small.png", id = "logo", height = "40px", width="80px", style = "position: relative; margin:-5px 0px; display:right-align;")),
                  id="navbar",
                  position = c("fixed-top"),
                  #header = tagList(useShinydashboard()),
                  tabPanel("Main",
                           fluidRow(
                             br(),
                             br(),
                             br(),
                             br(),
                             column(width=3,
                                   downloadButton("downloadTemplate", "Download template", style = "font-size: 12px; padding: 6px; background-color: #7EBC72; border-color: #7EBC72; color: white;"),
                                   actionButton(inputId="setExample", label="Run example", style = "font-size: 12px; padding: 6px; background-color: #70B7D7; border-color: #70B7D7; color: white;"),
                                   br(),
                                   radioButtons("source", label = h4("Source", bsButton("infosource", label="", icon=icon("info"), style="info", size="extra-small")), choices=c("GREAT"="Great","IPA","GSEA", "Template"="GSEA"), selected = character(0), inline=TRUE),
                                   bsPopover(id="infosource", title="Source", content="GREAT, IPA or GSEA inputs", placement = "right", trigger="hover"),
                                   radioButtons("gene.id", label = h4("Gene ID", bsButton("infogeneid", label="", icon=icon("info"), style="info", size="extra-small")), choices=c("Ensembl ID"="ENSEMBL", "Symbol"="SYMBOL", "Entrez ID"="ENTREZID"), selected = character(0), inline=TRUE),
                                   bsPopover(id="infogeneid", title="Gene ID", content="Gene identification used", placement = "right", trigger="hover"),
                                   radioButtons("organism", label = h4("Organism", bsButton("infoorganism", label="", icon=icon("info"), style="info", size="extra-small")), choices=c("Homo sapiens"="org.Hs.eg.db", "Mus musculus"="org.Mm.eg.db"), selected = character(0), inline=TRUE),
                                   bsPopover(id="infoorganism", title="Organism", content="Organism of your data", placement = "right", trigger="hover"),
                                   fileInput("files", h4("Upload files",bsButton("infofiles", label="", icon=icon("info"), style="info", size="extra-small")), multiple=TRUE, accept = c("text/csv","text/comma-separated-values,text/plain",".csv", ".txt", ".xls", ".xlsx")),
                                   bsPopover(id="infofiles", title="Files input", content=".txt or .csv files", placement = "right", trigger="hover"),
                                   DTOutput("filesInfo"),
                                   br(),
                                   actionButton(inputId="printResults", label="Run analysis", icon=icon("search", lib = "glyphicon")),actionButton(inputId="reset", label="Reset", icon=icon("refresh", lib = "glyphicon")),
                             ),
                             div(id="main1",
                               column(width=9,
                                     fluidRow(
                                       column(10,
                                         tabsetPanel(id="tabs", type="tabs",
                                                     tabPanel("Heatmap",
                                                          plotOutput("heatmap", height = "500px"),
                                                     ),
                                                     tabPanel("Treeplot",
                                                          textOutput("treeIntro"),
                                                          plotOutput("treeplot", height = "500px"),
                                                     ),
                                                     tabPanel("Tissue",
                                                          plotOutput("tissuePlot", height = "500px"),
                                                     )
                                         ),
                                         
                                       ),
                                       column(2,
                                         #div(style="display:inline-block",uiOutput("download"),width=6),
                                         div(style="display:inline-block",uiOutput("downloadData"), br(), width=6),
                                         div(style="display:inline-block",uiOutput("downloadDataFormats"), br(), width=6),
                                         h4("Clustering"),
                                         numericInput("cluster","Set number of clusters:",value=NULL,min=2, width="200px"),
                                         actionButton(inputId="recalculateClustering", label="Recalculate", icon=icon("search", lib = "glyphicon"))
                                         #),
                                       ),
                                     ),
                                  ),#column main panel
                              ) %>% shinyjs::hidden(),
                            ),
                            div(id="main2",
                              fluidRow(
                                 column(width=12,
                                     hr(style = "border-top: 1px solid #000000;"),
                                     uiOutput("checkboxCluster"),
                                     tabsetPanel(id="tabs2", type="pills",
                                                 tabPanel("Data",
                                                      br(),
                                                      column(9,
                                                        DTOutput("dataInfo"),
                                                      ),
                                                      column(3,
                                                        #selectInput("database", "Databases:", choices=c("","Human Phenotype Ontology (HPO)"), selected = NULL),
                                                        uiOutput("databases"),
                                                        uiOutput("databaseOptionInformation"),
                                                        uiOutput("databaseOptions"),
                                                        textOutput("databaseOptionTitle"),
                                                        verbatimTextOutput("databaseGenes"),
                                                        actionButton(inputId="calculateHighlight", label="Calculate", icon=icon("search", lib = "glyphicon")),
                                                        br(),
                                                        DTOutput("highlightInfo"),
                                                      )
                                                 ),
                                                 tabPanel("ORA",
                                                      div(style="display:inline-block;vertical-align:center;horizontal-align:center;", class="row-fluid", strong("Top"),),
                                                      div(style="display:inline-block; padding-left:10px; width:auto", class="row-fluid", selectInput("top", "", choices=c(5,10,15,20), selected = c(5), width="60px")),
                                                      div(style="display:inline-block;vertical-align:center;horizontal-align:center;padding-left:10px", class="row-fluid", strong("per cluster:")),
                                                      br(),
                                                      br(),
                                                      DTOutput("oraInfo"),
                                                 ),
                                                 tabPanel("Genes",
                                                          br(),
                                                          DTOutput("genesInfo"),
                                                 ),
                                                 tabPanel("Tissue enrichment",
                                                      br(),
                                                      #p("Press the button below if you want to calculate the tissue enrichment analysis:"),
                                                      uiOutput("tissueIntro"),
                                                      actionButton(inputId="calculateTissueEnrichment", label="Calculate", icon=icon("search", lib = "glyphicon")),
                                                      DTOutput("tissueOutput"),
                                                 ),
                                     ),
                                ),
                              ),
                            ) %>% shinyjs::hidden(),
                  ), # Navbar 1, tabPanel
                  #tabPanel("Video", br(), br(), br(), br(),"Videos page"),
                  tabPanel("Help",  br(), br(), br(), br(), "Please, click here to download the user guide:",
                           downloadButton("downloadUserGuide", "User guide", style = "font-size: 12px; padding: 6px; background-color: #E4B5BB; border-color: #E4B5BB; color: black;")),
                ) # navbarPage
) # fluidPage



