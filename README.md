# GeneSetClusterShiny
Shiny app for gene set cluster tool

### Install packages
```
list.of.packages = c("BiocManager","shiny","shinydashboard","shinyBS","shinythemes","DT","ggplot2","htmltools","shinyjs","ggnewscale","ggtree","GO.db","reshape2"
,"clusterProfiler","powerjoin","imputeTS","clustree","cluster","factoextra","GGally","shinyWidgets","org.Hs.eg.db","org.Mm.eg.db","dplyr","limma","stringr","shinyalert","jsonlite","doParallel","parallel","httr","utils","readxl","pbapply","RColorBrewer","patchwork","gridExtra","pheatmap")
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
```

### Run shiny app
```
shiny::runGitHub(repo = "GeneSetClusterShiny", username = "TranslationalBioinformaticsUnit", ref = "main")
```
