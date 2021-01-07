
# todo installation double l 


if (!requireNamespace("BiocManager", quietly = TRUE)){
  BiocManager::install(version="3.11")
}
  
if(!require("FlowRepositoryR")){
  BiocManager::install("FlowRepositoryR")
  library(FlowRepositoryR)
}

if(!require("CATALYST")){
  BiocManager::install("CATALYST")
  library(CATALYST)
}

if(!require("flowDensity")){
  BiocManager::install("flowDensity")
  library(flowDensity)
}

if(!require("flowCore")){
  BiocManager::install("flowCore")
  library(flowCore)
}

if(!require("FlowSOM")){
  BiocManager::install("FlowSOM")
  library(FlowSOM)
}

if(!require("SingleCellExperiment")){
  BiocManager::install("SingleCellExperiment")
  library(SingleCellExperiment)
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

# if(!require("devel")){
#   BiocManager::install("devel")
# }

# TODOD check if yiou can install in the version 4 in r oit worked for Katrin by github

if(!require("flowCut")){
  BiocManager::install("flowCut")
  library(flowCut)
}

if(!require("stringr")){
  install.packages("stringr")
  library(stringr)
}

if(!require("devtools")){
  install.packages("devtools")
  library(devtools)
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

#TODO: install cytoml and flowrep and correct flowcut from github
# 
# if (!require("CytoML")) {
#   BiocManager::install("CytoML")
#   library(CytoML)
# }

# library("devtools")
# devtools::install_github("RGLab/CytoML")

# library(CytoML)
if (!require("CytoML")) {
devtools::install_github("RGLab/CytoML", 
                         ref = "8ef06308392c83395b26f6143b149ebbda5164d5" )
library(CytoML)
}

if (!require("pals")) {
  install.packages("pals")
  library(pals)
}


