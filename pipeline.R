
source('~/Documents/CyTOF_workflow/CytofPipeline1/instalation.R')
source('~/Documents/CyTOF_workflow/CytofPipeline1/functions.R')


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

# set directory 
dir <- "/home/paulina/Documents/CyTOF_workflow/data"

# # set directory where bead-normalized fcs files will be saved
bead_norm_dir <- file.path(dir, "BeadNorm")
if(!dir.exists(bead_norm_dir)) dir.create(bead_norm_dir)

# set files input directory
raw_data_dir <- file.path(dir, "RawFiles")

# define sample to which all the files should be normalized and read in flowframe
ff_ref <- read.FCS(file.path(raw_data_dir, "181129_RUN1_01.FCS"))

# define which files will be normalized
files <- list.files(file.path(raw_data_dir), pattern = ".FCS$")

for (file in files){
  print(paste0("   ", Sys.time()))
  print(paste0("   Normalizing ", file))
  
  # read flow frame
  ff <- read.FCS(file.path(raw_data_dir, file))
  
  ff <- bead_normalize()
  
  
  
  
  
  
  
  
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

# set directory 
dir <- "/home/paulina/Documents/CyTOF_workflow/data"

# set directory where cleaned fcs files will be saved
clean_dir <- file.path(dir, "Cleaned")
if(!dir.exists(clean_dir)) dir.create(clean_dir)

# set files input directory
bead_norm_dir <- file.path(dir, "BeadNorm")

# Define which files will be normalized
files <- list.files(file.path(bead_norm_dir), pattern = "_beadNorm.fcs$")

for (file in files) {
  
  # read fcs file
  ff <- read.FCS(file.path(bead_norm_dir, file), transformation = FALSE)
  
  # norm_not_na <- which(apply(ff_t@exprs, 1, function(x){all(!is.na(x))}))
  # ff_t <- ff_t[norm_not_na, ]
  
  # clean Flow Rate and signal instability
  ff <- clean_flow_rate(flow_frame = ff, out_dir = clean_dir, 
                                to_plot = TRUE)
  
  # clean Signal 
  ff <- clean_signal(flow_frame = ff, to_plot = "All", out_dir = clean_dir)
 
  # Write FCS files
  write.FCS(ff,
            file = file.path(clean_dir, gsub("_beadNorm","_cleaned", file))) 
  
}





