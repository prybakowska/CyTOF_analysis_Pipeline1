#install CATALYST package 

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

# set data directory and subdirectories and locate and unzip the data 

# !!! TODO !!!
# this can be done wit flowrepositry package 
zipFile <- "/home/paulina/Downloads/FlowRepository_FR-FCM-ZYND_files.zip"
projectDir <- "/home/paulina/Documents/CyTOF_workflow/data"
dataDir <- file.path(projectDir, "RawFiles")
if(!dir.exists(dataDir))(dir.create(dataDir))

unzip(zipFile, exdir = dataDir)

dir <- "/home/paulina/Documents/CyTOF_workflow/data"
rawDir <- file.path(dir, "RawFiles")
dir.exists(rawDir)

# ------------------------------------------------------------------------------
# Bead normalization -----------------------------------------------------------
#--------------------------------------------------------------------------------

# set bead normalization directory where normalized fcs files will be saved
beadNormDir <- file.path(dir, "BeadNorm")
if(!dir.exists(beadNormDir)) dir.create(beadNormDir)

# define sample to which all the files should be normalized and read in flowframe
ff_ref <- read.FCS(file.path(rawDir, "181129_RUN1_01.FCS"))

# Define which files will be normalized
files <- list.files(file.path(rawDir), pattern = ".FCS$")

for (file in files){
  print(paste0("   ", Sys.time()))
  print(paste0("   Normalizing ", file))
  
  ff <- read.FCS(file.path(rawDir, file))
  
  # at this point for further analysis you can select only the markers the 
  # markers neccessary for the analysis, this will reduce the size of your data
  channels_to_keep <- c(grep("Time|Event_length|Pd|Ir|Ce140|Center|Offset|Width|Residual",
                             colnames(ff)),
                        grep("_", get_markers(ff, colnames(ff))))
  channels_to_keep <- colnames(ff)[sort(unique(channels_to_keep))]
  
  # prepare the data for bead normalization 
  dat <- prepData(ff[, channels_to_keep])
  
  # normalize the data and remove beads
  dat_norm <- normCytof(x = dat,
                        beads = "dvs",
                        remove_beads = TRUE,
                        norm_to = ff_ref,
                        k = 80,
                        plot = TRUE)
  
  # convert back to .fcs files and save 
  ff <- sce2fcs(dat_norm$data)
  write.FCS(ff, file.path(beadNormDir, sub_dir, gsub(".FCS","_beadNormDir.fcs", file)))
  
  # plot and save diagnostic plots 
  dat_norm$scatter
  ggsave(filename = file.path(beadNormDir, sub_dir, gsub(".FCS","_beadGate.png", file)))
  dat_norm$lines
  ggsave(filename = file.path(beadNormDir, sub_dir, gsub(".FCS","_beadLines.png", file)))
}








