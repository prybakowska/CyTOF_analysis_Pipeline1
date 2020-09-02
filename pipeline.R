
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

# define full pathway to files thatwill be normalized
files <- list.files(file.path(raw_data_dir), pattern = ".FCS$", full.names = T)

# create reference sample to which all the files will be normalized 
ref_sample <- create_ref(fcs_file = files[1], beads = "dvs")

# Normalize file by file in the loop, saving new file with each loop
for (file in files){
  
  # read fcs file
  ff <- read.FCS(file, transformation = FALSE, truncate_max_range = FALSE)
  
  # bead normalize the files
  ff_norm <- bead_normalize(flow_frame = ff, keep_all_markers = FALSE, 
                            out_dir = bead_norm_dir, norm_to_ref = ref_sample, 
                            to_plot = TRUE)
  
  # write FCS files
  write.FCS(ff_norm, filename = file.path(bead_norm_dir, 
                                 gsub(".FCS","_beadNorm.fcs", basename(file)))) 

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
files <- list.files(file.path(bead_norm_dir), pattern = "_beadNorm.fcs$", 
                    full.names = TRUE)

# Clean file by file in the loop, saving new file with each loop
for (file in files) {
  
  # read fcs file
  ff <- read.FCS(filename = file, transformation = FALSE)
  
  # norm_not_na <- which(apply(ff_t@exprs, 1, function(x){all(!is.na(x))}))
  # ff_t <- ff_t[norm_not_na, ]
  
  # clean Flow Rate and signal instability
  ff <- clean_flow_rate(flow_frame = ff, out_dir = clean_dir, to_plot = TRUE)
  
  # clean Signal 
  ff <- clean_signal(flow_frame = ff, to_plot = "All", out_dir = clean_dir, 
                     arcsine_transform = TRUE)
 
  # Write FCS files
  write.FCS(ff,
            file = file.path(clean_dir, gsub("_beadNorm","_cleaned", basename(file)))) 
}





