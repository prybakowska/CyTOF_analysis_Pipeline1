if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!require("CATALYST")){
  BiocManager::install("CATALYST")
  library(CATALYST)
}
#
if(!require("flowDensity")){
  BiocManager::install("flowDensity")
  library(flowDensity)
}

library(flowCore)
library(FlowSOM)
library(SingleCellExperiment)
library(ggplot2)


if(!require("remotes")){
  install.packages("remotes")
  # library(remotes)
}

if(!require("CytoNorm")){
  remotes::install_github("saeyslab/CytoNorm")
  library(CytoNorm)
}

if(!require("flowAI")){
  BiocManager::install("flowAI")
  library(flowAI)
}

BiocManager::install(version='devel')

if(!require("flowCut")){
  BiocManager::install("flowCut")
  library(flowCut)
}
