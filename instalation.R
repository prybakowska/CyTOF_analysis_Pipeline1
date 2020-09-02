
if (!requireNamespace("BiocManager", quietly = TRUE)){
  BiocManager::install(version="3.11")
}
  
if(!require("CATALYST")){
  BiocManager::install("CATALYST")
  library(CATALYST)
}

if(!require("flowDensity")){
  BiocManager::install("flowDensity")
  library(flowDensity)
}

library(flowCore)
library(FlowSOM)

# BiocManager::install("SingleCellExperiment")
library(SingleCellExperiment)

library(ggplot2)


if(!require("remotes")){
  install.packages("remotes")
  library(remotes)
}

if(!require("CytoNorm")){
  remotes::install_github("saeyslab/CytoNorm")
  library(CytoNorm)
}

if(!require("flowAI")){
  BiocManager::install("flowAI")
  library(flowAI)
}

# if(!require("devel")){
#   BiocManager::install("devel")
# }

if(!require("flowCut")){
  BiocManager::install("flowCut")
  library(flowCut)
}

