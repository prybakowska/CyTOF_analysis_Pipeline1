
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
# Files debarcoding ------------------------------------------------------------
#-------------------------------------------------------------------------------

# Set input directory 
clean_dir <- file.path(dir, "Cleaned")

# Define files for debarcoding
files <- list.files(file.path(clean_dir), pattern = "_cleaned.fcs$", 
                    full.names = TRUE)

# Define out_dir for diagnostic plots
debarcode_dir <- file.path(dir, "Debarcoded")

file_scores <- readRDS(list.files(dir, recursive = TRUE, 
                               pattern = "AOF_sample_scores.RDS"))

good_files <- file_scores$file_names[file_scores$quality == "good"]
fcs_files_clean <- files[basename(files) %in% good_files]

file_batch_id <- stringr::str_match(basename(fcs_files_clean), 
                                    "(RUN[0-9]*)_[0-9]*_.*.fcs")[,2]

barcodes_list <- list("RUN1" = rownames(sample_key)[c(1:16)], 
                      "RUN2" = rownames(sample_key)[c(5:20)])

debarcode_files(fcs_files = fcs_files_clean, 
                out_dir = debarcode_dir, min_threshold = TRUE, 
                barcodes_used = barcodes_list, file_batch_id = file_batch_id)


#TODO check bad quality files
# ------------------------------------------------------------------------------
# Files aggregation ------------------------------------------------------------
#-------------------------------------------------------------------------------

# Set input directory 
debarcode_dir <- file.path(dir, "Debarcoded")

# Define files for debarcoding
files <- list.files(file.path(debarcode_dir), pattern = "_debarcoded.fcs$", 
                    full.names = TRUE, recursive = T)

# Define out_dir for aggtegated files
aggregate_dir <- file.path(dir, "Aggregated")
if(!dir.exists(aggregate_dir))(dir.create(aggregate_dir))

# Bring meta data 
md <- read.csv(file.path(dir, "meta_data.csv"))

# assign barcodes names the to barcodes 
md$barcode_name <- paste0(rownames(sample_key)[md$BARCODE])

# assign new sample names specifying patient id and its batch name 
md$fcs_new_name <- paste0(md$PATIENT_ID, "_", md$RUN, ".fcs")  

# aggregate files batch by batch 
for (i in seq_len(nrow(md))){
  
  patterns <- as.character(md[i, c("barcode_name", "RUN")])
  
  files_to_agg <- grep(pattern = patterns[2],
                       grep(pattern = patterns[1], files, value = TRUE), 
                       value = TRUE)
  
  print(paste0("Creating ", md[[i, "fcs_new_name"]]))
  
  aggregate_files(fcs_files = files_to_agg,
                  outputFile = file.path(aggregate_dir, md[[i, "fcs_new_name"]]))
}


# ------------------------------------------------------------------------------
# Files gating -----------------------------------------------------------------
#-------------------------------------------------------------------------------

# open cytofcelan GUI and select the file sthat you want to gate, 
# de-select bead removal 

cytofclean::cytofclean_GUI()

# Set input directory 
cytof_clean_dir <- file.path(dir, "Aggregated", "CyTOFClean")

# Set output directory 
gate_dir <- file.path(dir, "Gated")

if (!dir.exists(gate_dir)) { 
  dir.create(gate_dir)
}

files <- list.files(cytof_clean_dir, pattern = ".fcs$")
gate(fcs_files = files, cd45_ch = "Pr141Di", viability_ch = "Pt195Di",
     out_dir = gate_dir)


# ------------------------------------------------------------------------------
# Normalization using reference sample -----------------------------------------
#-------------------------------------------------------------------------------

# Set input directory 
gate_dir <- file.path(dir, "Gated")

# Define files for debarcoding
files_ref <- list.files(file.path(gate_dir), pattern = "REF.*_gated.fcs$", 
                    full.names = TRUE, recursive = T)

#
labels_ref <- stringr::str_match(basename(files_ref), 
                                   ".*_(RUN[0-9]*).*.fcs")[,2]

ff <- read.FCS(files[1])
channels <- grep("Pd|Rh", grep("Di", colnames(ff), value = T), 
                 value = T, invert = T) #TODO channels should not have pd 

# Define out_dir for aggtegated files
norm_dir <- file.path(dir, "CytoNormed")
if(!dir.exists(norm_dir))(dir.create(norm_dir))

png(file.path(norm_dir, "001_099_normalization.png"),
    width = length(channels) * 300,
    height = (length(files_ref) * 2 + 1) * 300)
model <- QuantileNorm.train(files = files_ref,  labels = labels_ref, 
                            channels = channels, 
                            transformList = transformList(channels, 
                                                          cytofTransform), 
                            nQ = 2, quantileValues = c(0.01, 0.99), goal = "mean", 
                            plot = TRUE)
dev.off()

saveRDS(object = model, file = file.path(norm_dir, "001_099_model.RDS"))


files <- list.files(file.path(gate_dir), pattern = "_gated.fcs$", 
                        full.names = TRUE, recursive = T)

labels <- stringr::str_match(basename(files), 
                             ".*_(RUN[0-9]*).*.fcs")[,2]

QuantileNorm.normalize(model = model, files = files, labels = labels, 
                       transformList = transformList(channels, 
                                                     cytofTransform),
                       transformList.reverse = transformList(channels, 
                                                             cytofTransform.reverse), 
                       outputDir = norm_dir)

# show how to see the batch effect 

files_before_norm <- list.files(gate_dir, pattern = ".fcs", full.names = T)
files_after_norm <- list.files(norm_dir, pattern = ".fcs", full.names = T)

plot_batch(files_before_norm = files_before_norm, files_after_norm = files_after_norm,
           out_dir = norm_dir, batch_pattern = "RUN[0-9]*" )

# TODO add option for plotting cytokines!!!!








