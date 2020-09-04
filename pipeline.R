
source('~/Documents/CyTOF_workflow/CytofPipeline1/instalation.R')
source('~/Documents/CyTOF_workflow/CytofPipeline1/functions.R')


# set working directory

dir <- "/home/paulina/Documents/CyTOF_workflow/data"
setwd(dir)


# # !!! TODO !!!
# # this can be done wit flowrepositry package 
# zipFile <- "/home/paulina/Downloads/FlowRepository_FR-FCM-ZYND_files.zip"
# projectDir <- "/home/paulina/Documents/CyTOF_workflow/data"
# dataDir <- file.path(projectDir, "RawFiles")
# if(!dir.exists(dataDir))(dir.create(dataDir))
# 
# unzip(zipFile, exdir = dataDir)
# 
# dir <- "/home/paulina/Documents/CyTOF_workflow/data"
# raw_data_dir <- file.path(dir, "RawFiles")
# dir.exists(raw_data_dir)

# ------------------------------------------------------------------------------
# Bead normalization -----------------------------------------------------------
#-------------------------------------------------------------------------------

# set input directory (pathway to the files that are going to be normalized)
raw_data_dir <- file.path(dir, "RawFiles")

# set and create the directory where bead-normalized fcs files will be saved
bead_norm_dir <- file.path(dir, "BeadNorm")
if(!dir.exists(bead_norm_dir)) dir.create(bead_norm_dir)

# define full pathway to the files that you want to normalize
files <- list.files(file.path(raw_data_dir), pattern = ".FCS$", full.names = T)

# create reference sample to which all the files will be normalized 
ref_sample <- create_ref(fcs_file = files[1], beads = "dvs")

# Normalize file by file in the loop, saving new file with each loop execution
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

# set input directory (pathway to the files that are going to be cleaned)
bead_norm_dir <- file.path(dir, "BeadNorm")

# set and create the directory where cleaned fcs files will be saved
clean_dir <- file.path(dir, "Cleaned")
if(!dir.exists(clean_dir)) dir.create(clean_dir)

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

# ------------------------------------------------------------------------------
# Files inspection -------------------------------------------------------------
#-------------------------------------------------------------------------------

# Set input directory 
clean_dir <- file.path(dir, "Cleaned")

# Define files for visualization
files <- list.files(file.path(clean_dir), pattern = "_cleaned.fcs$", 
                    full.names = TRUE)

# Define batch id and sample id for each file
file_batch_id <- str_match(basename(files), "(RUN[0-9]*)_[0-9]*_.*.fcs")[,2]
file_id <- str_match(basename(files), "RUN[0-9]*_([0-9]*)_.*.fcs")[,2]

plot_marker_quantiles(fcs_files = files, file_batch_id = file_batch_id, 
                      file_id = file_id, arcsine_transform = TRUE, 
                      out_dir = clean_dir, markers_to_plot = NULL)


# ------------------------------------------------------------------------------
# Files outliers detection -----------------------------------------------------
#-------------------------------------------------------------------------------

# Set input directory 
clean_dir <- file.path(dir, "Cleaned")

# Define files for visualization
files <- list.files(file.path(clean_dir), pattern = "_cleaned.fcs$", 
                    full.names = TRUE)

# Define batch_id for each file 
file_batch_id <- stringr::str_match(basename(files), 
                                    "(RUN[0-9]*)_[0-9]*_.*.fcs")[,2]

# Define out_dir for diagnostic plots
quality_dir <- file.path(dir, "Quality_control")

file_quality_check(fcs_files = files, file_batch_id = file_batch_id, 
                   out_dir = quality_dir,
                   phenotyping_markers = c("Ir","CD", "HLA", "IgD", "Pt"), 
                   arcsine_transform = TRUE)


# ------------------------------------------------------------------------------
# Files outliers detection -----------------------------------------------------
#-------------------------------------------------------------------------------

# Set input directory 
clean_dir <- file.path(dir, "Cleaned")

# Define files for debarcoding
files <- list.files(file.path(clean_dir), pattern = "_cleaned.fcs$", 
                    full.names = TRUE)

# Define out_dir for diagnostic plots
debarcode_dir <- file.path(dir, "Debarcoded")

bad_quality_file <- list.files(dir, recursive = TRUE, 
                               pattern = "AOF_sample_scores.RDS")

debarcode_files <- function(fcs_files = files, bad_quality_files = NULL, 
                            out_dir = debarcode_dir){
   
  if(!is.null(bad_quality_files)){
    file_scores <- readRDS(bad_quality_file)
  }
  
  for (file in fcs_files){
    print(paste0("   ", Sys.time()))
    print(paste0("   Debarcoding ", file))
    ff <- read.FCS(file)
    
    file_label <- basename(file)
    if(!dir.exists(out_dir)) dir.create(out_dir)
    
    
    
    debarcode(ff = ff, 
              sample_key = CATALYST::sample_key,
              output_dir = debarcode_tmp_dir,
              ref = NULL,
              min_threshold = 0.18)
  }
}

























