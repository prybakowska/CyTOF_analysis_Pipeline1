#install CATALYST package 

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

BiocManager::install(version='devel')

if(!require("flowCut")){
  BiocManager::install("flowCut")
  library(flowCut)
}

# set data directory and subdirectories and locate and unzip the data 

# !!! TODO !!!
# this can be done wit flowrepositry package 
zipFile <- "/home/paulina/Downloads/FlowRepository_FR-FCM-ZYND_files.zip"
projectDir <- "/home/paulina/Documents/CyTOF_workflow/data"
dataDir <- file.path(projectDir, "RawFiles")
if(!dir.exists(dataDir))(dir.create(dataDir))

unzip(zipFile, exdir = dataDir)

dir <- "/home/paulina/Documents/CyTOF_workflow/data"
raw_data_dir <- file.path(dir, "RawFiles")
dir.exists(raw_data_dir)

# ------------------------------------------------------------------------------
# Bead normalization -----------------------------------------------------------
#-------------------------------------------------------------------------------

# set directory where normalized fcs files will be saved
bead_norm_dir <- file.path(dir, "BeadNorm")
if(!dir.exists(bead_norm_dir)) dir.create(bead_norm_dir)

# define sample to which all the files should be normalized and read in flowframe
ff_ref <- read.FCS(file.path(raw_data_dir, "181129_RUN1_01.FCS"))

# Define which files will be normalized
files <- list.files(file.path(raw_data_dir), pattern = ".FCS$")

for (file in files){
  print(paste0("   ", Sys.time()))
  print(paste0("   Normalizing ", file))
  
  ff <- read.FCS(file.path(raw_data_dir, file))
  
  # at this point for further analysis you can select only the markers
  # necessary for the analysis, this will reduce the size of your data
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
  write.FCS(ff, file.path(bead_norm_dir, sub_dir, gsub(".FCS","_beadNorm.fcs", file)))
  
  # plot and save diagnostic plots 
  dat_norm$scatter
  ggsave(filename = file.path(bead_norm_dir, sub_dir, gsub(".FCS","_beadGate.png", file)))
  dat_norm$lines
  ggsave(filename = file.path(bead_norm_dir, sub_dir, gsub(".FCS","_beadLines.png", file)))
}


# ------------------------------------------------------------------------------
# Signal Cleaning --------------------------------------------------------------
#-------------------------------------------------------------------------------

# set directory where cleaned fcs files will be saved
clean_dir <- file.path(dir, "Cleaned")
if(!dir.exists(clean_dir)) dir.create(clean_dir)

# Define which files will be normalized
files <- list.files(file.path(bead_norm_dir), pattern = "_beadNorm.fcs$")

for (file in files) {
  
  # read fcs files and arcshinh transform 
  ff <- read.FCS(file.path(bead_norm_dir,file), transformation = FALSE)
  
  # find mass channels
  channels_to_transform <- find_mass_ch(ff)
  
  #transfrom mass channels
  ff_t <- transform(ff, transformList(channels_to_transform, cytofTransform))
  
  # norm_not_na <- which(apply(ff_t@exprs, 1, function(x){all(!is.na(x))}))
  # ff_t <- ff_t[norm_not_na, ]
  
  # clean Flow Rate 
  cleaned_data <- clean_flow_rate(flow_frame = ff_t, out_dir = clean_dir, plot_rate = TRUE) 
  
  channels_to_clean <- c(grep("Event_length", colnames(ff)),
                         grep("Ce140Di", colnames(ff)))
  
  flowCut_res <- flowCut(cleaned_data$FRnewFCS[,-channels_to_clean], 
                         Segment = 1000, 
                         MaxPercCut = 0.5,
                         FileID = gsub("_normalised", "", file),
                         Plot = "All",
                         Directory = file.path(clean_dir, sub_dir),
                         UseOnlyWorstChannels = TRUE,
                         AllowFlaggedRerun = TRUE)
  
  #subdir <- gsub("_normalised", "", file)
  # selection <- norm_not_na
  selection <- cleaned_data$goodCellIDs
  if(length(flowCut_res$ind) > 0) {
    selection <- selection[-flowCut_res$ind]
  }
  
  # Write out selected subset of untransformed data
  write.FCS(ff[selection, ],
            file = file.path(clean_dir,sub_dir, gsub("_normalised", "", file))) 
  
}





