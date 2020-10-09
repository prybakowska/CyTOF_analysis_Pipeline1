
source('~/Documents/CyTOF_workflow/CytofPipeline1/instalation.R')
source('~/Documents/CyTOF_workflow/CytofPipeline1/functions.R')


# set working directory
dir <- "/home/paulina/Documents/CyTOF_workflow/data_cyt"
setwd(dir)

# ------------------------------------------------------------------------------
# Download data from flowrepository --------------------------------------------
#-------------------------------------------------------------------------------

# connect to flowrepository data
ds <- flowRep.get("FR-FCM-ZZJ7") # TODOD change the number to your flowrepository 
summary(ds)

# Download the flowrepository data
ds <- download(object = ds, dirpath = file.path(dir, "RawFiles")) 

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
ref_sample <- create_ref(fcs_files = files, beads = "dvs", out_dir = bead_norm_dir)

# Normalize file by file in the loop, saving new file with each loop execution
for (file in files){
  
  # read fcs file
  ff <- read.FCS(file, transformation = FALSE, truncate_max_range = FALSE)
  
  # bead normalize the files
  ff_norm <- bead_normalize(flow_frame = ff, keep_all_markers = FALSE, 
                            out_dir = bead_norm_dir, norm_to_ref = ref_sample, 
                            to_plot = TRUE, k = 80,
                            markers_to_keep = c("CD", "HLA", "IgD", "TCR", "Ir", 
                                                "Viability","IL",
                                                "TNF", "TGF", "MIP", "MCP", "Granz"))
  
  # write FCS files
  write.FCS(ff_norm, filename = file.path(bead_norm_dir, 
                                 gsub(".FCS","_beadNorm.fcs", basename(file)))) 
}

# ------------------------------------------------------------------------------
# Visualized files after bead normalization  -----------------------------------
#-------------------------------------------------------------------------------

# Set input directory 
raw_data_dir <- file.path(dir, "RawFiles")
bead_norm_dir <- file.path(dir, "BeadNorm")

# Define files for visualization

# before normalization 
files_b <- list.files(file.path(raw_data_dir), pattern = ".FCS$", 
                    full.names = TRUE)

# after normalization 
files_a <- list.files(file.path(bead_norm_dir), pattern = "_beadNorm.fcs$", 
                    full.names = TRUE)

# Define batch id and sample id for each file
batch_pattern <- str_match(basename(files_b), ".*(day[0-9]*).*.FCS")[,2]
# file_id <- str_match(basename(files_b), ".*day[0-9]*_[0-9]*.*.FCS$")[,2]

plot_marker_quantiles(files_after_norm = files_a, files_before_norm = files_b, 
                      batch_pattern = batch_pattern, arcsine_transform = TRUE,
                      file_id = NULL,
                      markers_to_plot = c("CD", "HLA", "IgD", "IL", "TNF",
                                          "TGF", "GR"), 
                      out_dir = bead_norm_dir)

# ------------------------------------------------------------------------------
# Signal Cleaning --------------------------------------------------------------
#-------------------------------------------------------------------------------

# set input directory (pathway to the files that are going to be cleaned)
bead_norm_dir <- file.path(dir, "BeadNorm")

# set and create the directory where cleaned fcs files will be saved
clean_dir <- file.path(dir, "Cleaned")
if(!dir.exists(clean_dir)) dir.create(clean_dir)

# Define which files will be cleaned
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
  
  #TODO think which parametesr from flocut to put as user defined
 
  # Write FCS files
  write.FCS(ff,
            file = file.path(clean_dir, gsub("_beadNorm","_cleaned", basename(file)))) 
}

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
                                    "(day[0-9]*)_[0-9]*_.*.fcs")[,2]

# Define out_dir for diagnostic plots
quality_dir <- file.path(dir, "Quality_control")

file_quality_check(fcs_files = files, file_batch_id = file_batch_id, 
                   out_dir = quality_dir,
                   phenotyping_markers = c("Ir","CD", "HLA", "IgD", "Pt"), 
                   arcsine_transform = TRUE, sd = 3)
#TODO put the flowsom parameters in to the user selection
#TODO make the heatmap plot with the scores


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

# Read in files scores
file_scores <- readRDS(list.files(dir, recursive = TRUE, 
                               pattern = "AOF_sample_scores.RDS"))

# Select good quality files
good_files <- file_scores$file_names[file_scores$quality == "good"]
fcs_files_clean <- files[basename(files) %in% good_files]

# Define file batch ID for each file
file_batch_id <- stringr::str_match(basename(fcs_files_clean), 
                                    "(day[0-9]*)_[0-9]*_.*.fcs")[,2]

#TODO check if you can use meta data to read barcodes instead of doing this 
# Define which barcodes were used in each batch 
barcodes_list <- list("day1" = rownames(sample_key)[c(2:6, 8:13, 15)], 
                      "day2" = rownames(sample_key)[c(6:17)],
                      "day3" = rownames(sample_key)[c(8:19)])

# Deabarcode files 
debarcode_files(fcs_files = fcs_files_clean, 
                out_dir = debarcode_dir, min_threshold = TRUE, 
                barcodes_used = barcodes_list, file_batch_id = file_batch_id)
#TODO plot less cells in the barcoding file  
#TODO user defined threshold and make infromatrion about flagging 

# TODO check bad quality files and maybe put heatmap scores 
# ------------------------------------------------------------------------------
# Files aggregation ------------------------------------------------------------
# ------------------------------------------------------------------------------

# Set input directory 
debarcode_dir <- file.path(dir, "Debarcoded")

# Define files for debarcoding
files <- list.files(file.path(debarcode_dir), pattern = "_debarcoded.fcs$", 
                    full.names = TRUE, recursive = T)

# Define out_dir for aggregated files
aggregate_dir <- file.path(dir, "Aggregated")
if(!dir.exists(aggregate_dir))(dir.create(aggregate_dir))

# Bring metadata 
md <- read.csv(file.path(dir, "meta_data.csv"))

# assign barcodes names the to barcodes 
md$barcode_name <- paste0(rownames(sample_key)[md$BARCODE])

# assign new sample names specifying patient id and its batch name 
md$fcs_new_name <- paste0(md$ID, "_", md$STIM, "_", md$BATCH, ".fcs")  

# TODO aggregate files put this part in teh function so it is shorter
# aggregate files batch by batch 
for (i in seq_len(nrow(md))){
  
  patterns <- as.character(md[i, c("barcode_name", "BATCH")])
  
  files_to_agg <- grep(pattern = patterns[2],
                       grep(pattern = patterns[1], files, value = TRUE), 
                       value = TRUE)
  
  print(paste0("Creating ", md[[i, "fcs_new_name"]]))
  
  aggregate_files(fcs_files = files_to_agg,
                  outputFile = file.path(aggregate_dir, md[[i, "fcs_new_name"]]),
                  write_agg_file = TRUE)
}


# ------------------------------------------------------------------------------
# Files gating -----------------------------------------------------------------
#-------------------------------------------------------------------------------

# open cytofcelan GUI and select the file that you want to gate, 
# de-select bead removal 
cytofclean::cytofclean_GUI()

# Set input directory 
cytof_clean_dir <- file.path(dir, "Aggregated", "CyTOFClean")

# Set output directory 
gate_dir <- file.path(dir, "Gated")
if (!dir.exists(gate_dir)) { 
  dir.create(gate_dir)
}

# List files for gating 
files <- list.files(cytof_clean_dir, pattern = ".fcs$")

#TODO remove all the CD45 gating
# Gate the files
gate(fcs_files = files, viability_ch = "Pt195Di",
     out_dir = gate_dir)


# ------------------------------------------------------------------------------
# Normalization using reference sample -----------------------------------------
#-------------------------------------------------------------------------------

# Set input directory 
gate_dir <- file.path(dir, "Gated")

# Define files for debarcoding
files_ref <- list.files(file.path(gate_dir), pattern = "p1_REF.*_gated.fcs$", 
                    full.names = TRUE, recursive = T)

#
labels_ref <- stringr::str_match(basename(files_ref), 
                                   ".*_(day[0-9]*).*.fcs")[,2]

ff <- read.FCS(files_ref[1])
channels <- grep("Pd|Rh|140", grep("Di", colnames(ff), value = T), 
                 value = T, invert = T) #TODO channels should not have pd 

# Define out_dir for aggregated files
norm_dir <- file.path(dir, "CytoNormed")
if(!dir.exists(norm_dir))(dir.create(norm_dir))

png(file.path(norm_dir, "005_095_normalization.png"),
    width = length(channels) * 300,
    height = (length(files_ref) * 2 + 1) * 300)
model <- QuantileNorm.train(files = files_ref,  labels = labels_ref, 
                            channels = channels, 
                            transformList = transformList(channels, 
                                                          cytofTransform), 
                            nQ = 2, quantileValues = c(0.05, 0.95), goal = "mean", 
                            plot = TRUE)
dev.off()

saveRDS(object = model, file = file.path(norm_dir, "005_095_model.RDS"))


files <- list.files(file.path(gate_dir), pattern = "_gated.fcs$", 
                        full.names = TRUE, recursive = T)

labels <- stringr::str_match(basename(files), 
                             ".*_(day[0-9]*).*.fcs")[,2]

QuantileNorm.normalize(model = model, files = files, labels = labels, 
                       transformList = transformList(channels, 
                                                     cytofTransform),
                       transformList.reverse = transformList(channels, 
                                                             cytofTransform.reverse), 
                       outputDir = norm_dir)

# show how to see the batch effect 
batch_pattern <- "day[0-9]*"

files_before_norm <- list.files(gate_dir, pattern = ".fcs", full.names = T)
batch <- str_match(files_before_norm, "day[0-9]*")[,1]
files_before_norm <- files_before_norm[order(factor(batch))]

files_after_norm <- list.files(norm_dir, pattern = ".fcs", full.names = T)
batch <- str_match(files_after_norm, "day[0-9]*")[,1]
files_after_norm <- files_after_norm[order(factor(batch))]

plot_batch(files_before_norm = files_before_norm, 
           files_after_norm = files_after_norm,
           out_dir = norm_dir, batch_pattern = batch_pattern, 
           clustering_markers = c("CD", "IgD", "HLA"),
           functional_markers = c("IL", "Gran", "TNF", "TGF", "MIP", "MCP"))


# ------------------------------------------------------------------------------
# Runa UMAP y map manual labels ------------------------------------------------
#-------------------------------------------------------------------------------

# Set input directory 
gate_dir <- file.path(dir, "Gated")

# Set output directory 
analysis_dir <- file.path(dir, "Analysis")

if (!dir.exists(analysis_dir)) { 
  dir.create(analysis_dir)
}

files <- list.files(gate_dir, pattern = ".fcs$", full.names = TRUE)
batch_pattern <- str_match(basename(files), ".*(day[0-9]*).*.fcs")[,2]


plot_UMAP <- function(fcs_files = files, clustering_markers = c("CD", "HLA", "IgD"),
                      out_dir = analysis_dir, batch_pattern = "day[0-9]*"){
  
  set.seed(1)
  ff_agg <- AggregateFlowFrames(fileNames = fcs_files,
                                cTotal = length(fcs_files) * 1000,
                                verbose = TRUE,
                                writeMeta = FALSE,
                                writeOutput = FALSE,
                                outputFile = file.path(out_dir, paste0("norm_REF_aggregated_", ".fcs")))
  
  ff_agg <- transform(ff_agg,
                      transformList(grep("Di", colnames(ff_agg), value = TRUE), 
                                    cytofTransform))
  
  markers <- get_markers(ff_agg, colnames(ff_agg))
  
  cl_markers <- paste(clustering_markers, collapse="|")
  cl_markers <- grep(cl_markers, markers, value = T)
  
  ff_agg@exprs[, names(cl_markers)] <- apply(ff_agg@exprs[, names(cl_markers)], 2, function(x){
    q <- quantile(x, 0.9999)
    x[x > q] <- q
    x
  })
  
  samp <- length(fcs_files)
  
  set.seed(1)
  ff_samp <- ff_agg@exprs[sample(nrow(ff_agg@exprs), samp*500), ]
  
  
  #TODO use the ce of nocicka do show mds or umap but by sample ypu have a code in TODO doc in google
  set.seed(123)
  dimred_res <- uwot::umap(X = ff_samp[, names(cl_markers)], 
                           n_neighbors = 15, scale = TRUE)
  
  dimred_df <- data.frame(dim1 = dimred_res[,1], dim2= dimred_res[,2],
                          ff_samp[, names(cl_markers)])
  
  dimred_df$file_id <- ff_samp[,"File2"]
  dimred_df$batch <- NA
  
  for (i in 1:length(fcs_files)){
    
    file <- fcs_files[i]
    batch <- stringr::str_match(file, batch_pattern)[,1]
    dimred_df[dimred_df[, "file_id"] == i, "batch"] <- batch
  }
  
  p <- ggplot(dimred_df,  aes_string(x = "dim1", y = "dim2")) +
    geom_point(aes(color = batch), size = 0.8, position="jitter") +
    theme_bw() +
    # ggtitle()+
    # scale_color_manual(values=col )+
    scale_color_gradientn("Eu151Di", 
                           colours = grDevices::colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))+
    theme(legend.position = "bottom")
  
  p 
  
}






