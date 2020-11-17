# TODO in the files change the name for 165 as you do not haave this marker


source('~/Documents/CyTOF_workflow/CytofPipeline1/instalation.R')
source('~/Documents/CyTOF_workflow/CytofPipeline1/functions.R')

# set working directory
dir <- "/home/paulina/Documents/CyTOF_workflow/data_cyt"
setwd(dir)

# ------------------------------------------------------------------------------
# Download data from flowrepository --------------------------------------------
#-------------------------------------------------------------------------------

# connect to flowrepository data
ds <- flowRep.get("FR-FCM-Z3YR") 
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
  ff_norm <- bead_normalize(flow_frame = ff, 
                            out_dir = bead_norm_dir, norm_to_ref = ref_sample, 
                            to_plot = TRUE, k = 80,
                            markers_to_keep = c("CD", "HLA", "IgD", "TCR", "Ir", 
                                                "Viability","IL", "IFNa",
                                                "TNF", "TGF", "MIP", "MCP", "Granz"))
  
  # save normalized FCS files
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

plot_marker_quantiles(files_after_norm = files_a, files_before_norm = files_b, 
                      batch_pattern = batch_pattern, arcsine_transform = TRUE,
                      markers_to_plot = c("CD", "HLA", "IgD", "IL", "TNF",
                                          "TGF", "GR", "IF"), 
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
                     arcsine_transform = TRUE, 
                     non_used_bead_ch = "140")
  # TOSOFIE TOKATRIEN shoudl we take more parameters from flowcut as user defined? 

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
                                    "(day[0-9]*).*.fcs")[,2]

# Define out_dir for diagnostic plots
quality_dir <- file.path(dir, "Quality_control")


#ToDO adapt code to the new function

file_quality_check(fcs_files = files, file_batch_id = file_batch_id, 
                   out_dir = quality_dir,
                   phenotyping_markers = c("Ir","CD", "HLA", "IgD", "Pt"), 
                   arcsine_transform = TRUE, 
                   sd = 3)

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
                                    "(day[0-9]*).*.fcs")[,2]

# Read in meta data 
md <- read.csv(file.path(dir, "meta_data.csv"))

# read in barcode key 
sample_key <- CATALYST::sample_key

# Extract information about barcodes used in each batch
barcodes_list <- list()
for (batch in unique(file_batch_id)){
  idx <- md[md[,"BATCH"] == batch, "BARCODE"]
  barcodes_list[[batch]] <- rownames(sample_key)[idx] 
}

# Alternatively define barcodes list manually 
# barcodes_list <- list("day1" = rownames(sample_key)[c(2:6, 8:13, 15)],
#                       "day2" = rownames(sample_key)[c(6:17)],
#                       "day3" = rownames(sample_key)[c(8:19)])

# Deabarcode files 
debarcode_files(fcs_files = fcs_files_clean, 
                out_dir = debarcode_dir, min_threshold = TRUE, 
                barcodes_used = barcodes_list, file_batch_id = file_batch_id, 
                less_than_th = TRUE, barcode_key = sample_key)

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
out_dir <- file.path(dir, "Gated")
if (!dir.exists(gate_dir)) { 
  dir.create(gate_dir)
}

# List files for gating 
files <- list.files(cytof_clean_dir, pattern = ".fcs$", full.names = TRUE)

# Gate the files and plot the gating strategy for each file 
n_plots <- 2  
png(file.path(gate_dir,
              paste0("gating.png")),
    width = n_plots * 300, height = length(files) * 300)
layout(matrix(1:(length(files) * n_plots), ncol = n_plots, byrow = TRUE))

for (file in files){
  
  ff <- read.FCS(file, transformation = FALSE)
  
  ff <- gate_intact_cells(flow_frame = ff, file_name = NULL)
  ff <- gate_live_cells(flow_frame = ff, viability_channel = "Pt195Di",
                        out_dir = gate_dir)
  
  write.FCS(ff, file.path(out_dir,
                      gsub(".fcs", "_gated.fcs", basename(file))))
}

dev.off()


# ------------------------------------------------------------------------------
# Normalization using reference sample -----------------------------------------
#-------------------------------------------------------------------------------

# Set input directory 
gate_dir <- file.path(dir, "Gated")

# Define reference samples
files_ref <- list.files(file.path(gate_dir), pattern = "REF.*_gated.fcs$", 
                    full.names = TRUE, recursive = T)

# Define batch labels for each files
labels_ref <- stringr::str_match(basename(files_ref), 
                                   ".*_(day[0-9]*).*.fcs")[,2]
# Define channels to be normalized
ff <- read.FCS(files_ref[1])
channels <- grep("Pd|Rh|140", grep("Di", colnames(ff), value = T), 
                 value = T, invert = T) #TODO channels should not have pd 

# Define out_dir for normalized files
norm_dir <- file.path(dir, "CytoNormed")
if(!dir.exists(norm_dir))(dir.create(norm_dir))

# Build the normalization model using reference samples and plot quantiles 
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

# save the model
saveRDS(object = model, file = file.path(norm_dir, "005_095_model.RDS"))

# Define path to the files for normalization
files <- list.files(file.path(gate_dir), pattern = "_gated.fcs$", 
                        full.names = TRUE, recursive = T)

# Define batch labels for each files, note that they need to corresponds to 
# reference labels 
labels <- stringr::str_match(basename(files), 
                             ".*_(day[0-9]*).*.fcs")[,2]

# Normalize files 
QuantileNorm.normalize(model = model, files = files, labels = labels, 
                       transformList = transformList(channels, 
                                                     cytofTransform),
                       transformList.reverse = transformList(channels, 
                                                             cytofTransform.reverse), 
                       outputDir = norm_dir)


# ------------------------------------------------------------------------------
# Plot batch effect ------------------------------------------------------------
#-------------------------------------------------------------------------------

# Define batch pattern
batch_pattern <- "day[0-9]*"

# Define files before normalization and order them according to the batch 
# Set input directory for files before CytoNorm normalization 
gate_dir <- file.path(dir, "Gated")
files_before_norm <- list.files(gate_dir, pattern = ".fcs", full.names = T)
batch <- str_match(files_before_norm, "day[0-9]*")[,1]
files_before_norm <- files_before_norm[order(factor(batch))]

# Define files after normalization and order them according to the batch 
# Set input directory for files after CytoNorm normalization
norm_dir <- file.path(dir, "CytoNormed")
files_after_norm <- list.files(norm_dir, pattern = ".fcs", full.names = T)
batch <- str_match(files_after_norm, "day[0-9]*")[,1]
files_after_norm <- files_after_norm[order(factor(batch))]

# Plot batch effect 
plot_batch(files_before_norm = files_before_norm, 
           files_after_norm = files_after_norm,
           out_dir = norm_dir, batch_pattern = batch_pattern, 
           clustering_markers = c("CD", "IgD", "HLA"),
           functional_markers = c("IL", "Gran", "TNF", "TGF", "MIP", "MCP"))

# ------------------------------------------------------------------------------
# Run UMAP ---------------------------------------------------------------------
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

# Build UMAP on aggregated files
UMAP <- UMAP(fcs_files = files, clustering_markers = c("CD", "HLA", "IgD"),
             functional_markers = c("IL", "TNF", "TGF", "Gr", "IF", "MIP", "MCP1"),
             out_dir = analysis_dir, batch_pattern = "day[0-9]*", 
             arcsine_transform = TRUE)

saveRDS(UMAP, file.path(analysis_dir, "UMAP.RDS"))

# Plot UMAP between two donors using RSQ stimulation and day 1 batch

# Add additional column to be able to do facet in ggplot
UMAP$indyvidual <- str_match(UMAP$sample_name, "p1|p2|REF")[,1]

# filter UMAP data frame for the needed data
df <- UMAP %>% 
  filter(batch == "day1" & 
           indyvidual %in% c("p1", "p2"))

df <- df %>% filter(sample_name == grep("RSQ", df[,"sample_name"], 
                                        value = TRUE))

# select markers to plot 
markers_to_plot <- grep("CD|HLA|IgD",
                           colnames(df), value = TRUE)


# create a list to store the singel plots
plots <- list()

# plot the map for each markers
for(m in markers_to_plot){
 
  p <- ggplot(df,  aes_string(x = "dim1", y = "dim2", color = m)) +
    geom_point(size = 4) +
    facet_wrap("indyvidual") +
    scale_color_gradientn(markers_to_plot[markers_to_plot == m], 
                          colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))  +
    theme(panel.background = element_rect(fill = "white", colour = "black",
                                          size = 2, linetype = "solid"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.subtitle = element_text(color="black", size=26, hjust = 0.95, face = "bold"),
          axis.text = element_text(size = 24, colour = "black"),
          axis.title = element_text(size = 20, colour = "black"), 
          strip.text.x = element_text(size = 23, color = "black"),
          strip.background = element_rect(fill = "white"), 
          legend.text = element_text(size = 14), 
          legend.title = element_text(size = 20),
          legend.position = "right") 
  
  plots[[m]] <- p
}

# plot the plots
# ggarrange(plotlist = plots, nrow=length(markers_to_plot), ncol= 1)
ggarrange(plotlist = plots, ncol = 3, nrow = 8)

