
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!require("devtools")){
  install.packages("devtools")
  library(devtools)
}

if(!require("flowCore")){
  BiocManager::install("flowCore")
  library(flowCore)
}

if(!require("SummarizedExperiment")){
  BiocManager::install("SummarizedExperiment")
  library(SummarizedExperiment)
}

if(!require("SingleCellExperiment")){
  BiocManager::install("SingleCellExperiment")
  library(SingleCellExperiment)
}

if(!require("openCyto")){
  BiocManager::install("openCyto")
  library(openCyto)
}

if (!require("CytoML")) {
  BiocManager::install("CytoML")
  library(CytoML)
}

if(!require("flowDensity")){
  BiocManager::install("flowDensity")
  library(flowDensity)
}

if(!require("FlowSOM")){
  BiocManager::install("FlowSOM")
  library(FlowSOM)
}

if(!require("CATALYST")){
  BiocManager::install("CATALYST")
  library(CATALYST)
}

if(!require("ggplot2")){
  install.packages("ggplot2")
  library(ggplot2)
}

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

# Cairo

if(!require("flowCut")){
  BiocManager::install("flowCut")
  library(flowCut)
}

if(!require("stringr")){
  install.packages("stringr")
  library(stringr)
}

if(!require("cytutils")){
  install_github("ismmshimc/cytutils")
  library(cytutils)
}

if(!require("pheatmap")){
  install.packages(pheatmap)
  library(pheatmap)
}

if(!require(cytofclean)){
  devtools::install_github("JimboMahoney/cytofclean", 
                           dependencies = TRUE)
  library("cytofclean")
}

if(!require("uwot")){
  install.packages("uwot")
  library(uwot)
}

if(!require("tidyverse")){
  install.packages("tidyverse")
  library(tidyverse)
}

if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}

if (!require("ggpubr")) {
  install.packages("ggpubr")
  library(ggpubr)
}

if (!require("flowWorkspace")) {
  BiocManager::install("flowWorkspace")
  library(flowWorkspace)
}


if (!require("pals")) {
  install.packages("pals")
  library(pals)
}

if (!require("scales")) {
  install.packages("scales")
  library(scales)
}


