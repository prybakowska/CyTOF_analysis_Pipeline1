
# set your working directory to the folder where the files were downloaded 
# using setwd()

#  execute
source('installation.R')
source('functions.R')

# create data folder where all analysis will be stored
if(!dir.exists("data")) dir.create("data")

# set it as a main directory 
dir <- file.path(getwd(), "data")

# ------------------------------------------------------------------------------
# Bead normalization -----------------------------------------------------------
#-------------------------------------------------------------------------------

# set input directory (pathway to the files that are going to be normalized)
raw_data_dir <- file.path(dir, "RawFiles")

# set and create the directory where bead-normalized fcs files will be saved
bead_norm_dir <- file.path(dir, "BeadNorm")
if(!dir.exists(bead_norm_dir)) dir.create(bead_norm_dir)

# define full pathway to the files that you want to normalize
files <- list.files(file.path(raw_data_dir), 
                    pattern = ".FCS$", 
                    full.names = T)

# create baseline file to which all the files will be normalized 
ref_sample <- baseline_file(fcs_files = files, 
                            beads = "dvs", 
                            out_dir = bead_norm_dir)

# Normalize file by file in the loop, saving new file with each loop execution
for (file in files){
  
  # read fcs file
  ff <- flowCore::read.FCS(file, transformation = FALSE, 
                 truncate_max_range = FALSE)
  
  # bead normalize the files
  ff_norm <- bead_normalize(flow_frame = ff, 
                            out_dir = bead_norm_dir, 
                            norm_to_ref = ref_sample, 
                            to_plot = TRUE, 
                            k = 80,
                            markers_to_keep = c("CD", "HLA", "IgD", "TCR", "Ir", 
                                                "Viability","IL", "IFNa",
                                                "TNF", "TGF", "MIP", "MCP", "Granz"))
  
  # save normalized FCS files
  flowCore::write.FCS(ff_norm, filename = file.path(bead_norm_dir, 
                                 gsub(".FCS","_beadNorm.fcs", basename(file), 
                                      ignore.case = TRUE))) 
}

# ------------------------------------------------------------------------------
# Visualized files after bead normalization  -----------------------------------
#-------------------------------------------------------------------------------

# Set input directory 
raw_data_dir <- file.path(dir, "RawFiles")
bead_norm_dir <- file.path(dir, "BeadNorm")

# Define files for visualization

# before normalization 
files_b <- list.files(raw_data_dir, 
                      pattern = ".FCS$", 
                      ignore.case = T,
                      full.names = TRUE)

# after normalization 
files_a <- list.files(bead_norm_dir, 
                      pattern = "_beadNorm.fcs$", 
                      ignore.case = T,
                      full.names = TRUE)

# Define batch id and sample id for each file
batch_pattern <- str_match(basename(files_b), "(?i).*(day[0-9]*).*.FCS")[,2]

plot_marker_quantiles(files_after_norm = files_a, 
                      files_before_norm = files_b, 
                      batch_pattern = batch_pattern, 
                      arcsine_transform = TRUE,
                      remove_beads = TRUE,
                      bead_channel = "140", 
                      uncommon_prefix = "_beadNorm.fcs|.FCS", 
                      markers_to_plot = c("CD", "HLA", "IgD", "IL", "TNF",
                                          "TGF", "GR", "IFNa"),
                      manual_colors = c("darkorchid4", "darkorange", "darkgreen"),
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
files <- list.files(bead_norm_dir,
                    ignore.case = TRUE, 
                    pattern = "_beadNorm.fcs$", 
                    full.names = TRUE)

# Clean file by file in the loop, saving new file with each loop
for (file in files) {
  
  # read fcs file
  ff <- flowCore::read.FCS(filename = file, 
                           transformation = FALSE)
  
  # clean flow rate 
  ff <- clean_flow_rate(flow_frame = ff, 
                        out_dir = clean_dir, 
                        to_plot = TRUE,
                        data_type = "MC")
  
  # clean signal 
  ff <- clean_signal(flow_frame = ff,
                     to_plot = "All",
                     out_dir = clean_dir,
                     Segment = 1000,
                     arcsine_transform = TRUE,
                     data_type = "MC",
                     non_used_bead_ch = "140")

  # Write FCS files
  flowCore::write.FCS(ff,
            file = file.path(clean_dir, gsub("_beadNorm","_cleaned", basename(file)))) 
}

# ------------------------------------------------------------------------------
# Files outliers detection -----------------------------------------------------
#-------------------------------------------------------------------------------

# Set input directory 
clean_dir <- file.path(dir, "Cleaned")

# Define files for visualization
files <- list.files(clean_dir, 
                    pattern = "_cleaned.fcs$", 
                    full.names = TRUE)

# Define batch_id for each file 
file_batch_id <- stringr::str_match(basename(files), 
                                    "(day[0-9]*).*.fcs")[,2]

# Define out_dir for diagnostic plots
quality_dir <- file.path(dir, "Quality_control")

file_quality_check(fcs_files = files, 
                   file_batch_id = file_batch_id, 
                   out_dir = quality_dir,
                   phenotyping_markers = c("Ir","CD", "HLA", "IgD", "Pt"), 
                   arcsine_transform = TRUE, 
                   nClus = 10,
                   sd = 3)

# ------------------------------------------------------------------------------
# Files debarcoding ------------------------------------------------------------
#-------------------------------------------------------------------------------

# Set input directory 
clean_dir <- file.path(dir, "Cleaned")

# Define files for debarcoding
files <- list.files(clean_dir, 
                    pattern = "_cleaned.fcs$", 
                    full.names = TRUE)

# Define out_dir for diagnostic plots
debarcode_dir <- file.path(dir, "Debarcoded")

# Read in files scores
file_scores <- readRDS(list.files(dir, 
                                  recursive = TRUE, 
                                  full.names = TRUE,
                                  pattern = "Quality_AOF_score.RDS"))

# Select good quality files
good_files <- file_scores$file_names[file_scores$quality == "good"]
fcs_files_clean <- files[basename(files) %in% good_files]

# Define file batch ID for each file
file_batch_id <- stringr::str_match(basename(fcs_files_clean), 
                                    "(day[0-9]*).*.fcs")[,2]

# Read in metadata 
md <- utils::read.csv(file.path(dir, "RawFiles", "meta_data.csv"))

# read in barcode key 
sample_key <- CATALYST::sample_key

# Extract information about barcodes used in each batch
barcodes_list <- list()
for (batch in unique(file_batch_id)){
  idx <- md[md[,"BATCH"] == batch, "BARCODE"]
  barcodes_list[[batch]] <- rownames(sample_key)[idx] 
}

# Debarcode files 
debarcode_files(fcs_files = fcs_files_clean, 
                out_dir = debarcode_dir, 
                min_threshold = TRUE, 
                barcodes_used = barcodes_list, 
                file_batch_id = file_batch_id, 
                less_than_th = TRUE, 
                barcode_key = sample_key)

# ------------------------------------------------------------------------------
# Files aggregation and file name deconvolution --------------------------------
# ------------------------------------------------------------------------------

# Set input directory 
debarcode_dir <- file.path(dir, "Debarcoded")

# Define files for debarcoding
files <- list.files(debarcode_dir, 
                    pattern = "_debarcoded.fcs$", 
                    full.names = TRUE, recursive = T)

# Define out_dir for aggregated files
aggregate_dir <- file.path(dir, "Aggregated")
if(!dir.exists(aggregate_dir))(dir.create(aggregate_dir))

# Bring metadata 
md <- utils::read.csv(file.path(dir, "RawFiles", "meta_data.csv"))

# assign barcodes names the to barcodes 
md$barcode_name <- paste0(rownames(sample_key)[md$BARCODE])

# assign new sample names specifying patient id and its batch name 
md$fcs_new_name <- paste0(md$ID, "_", md$STIM, "_", md$BATCH, ".fcs")  

# aggregate files batch by batch 
for (i in seq_len(nrow(md))){
  
  patterns <- as.character(md[i, c("barcode_name", "BATCH")])
  
  files_to_agg <- grep(pattern = patterns[2],
                       grep(pattern = patterns[1], 
                            files, value = TRUE), 
                       value = TRUE)
  
  print(paste0("Creating ", md[[i, "fcs_new_name"]]))
  
  aggregate_files(fcs_files = files_to_agg,
                  outputFile = md[[i, "fcs_new_name"]],
                  out_dir = aggregate_dir,
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
files <- list.files(aggregate_dir, 
                    pattern = ".fcs$", 
                    full.names = TRUE)

# Gate the files and plot the gating strategy for each file 
n_plots <- 3  
png(file.path(gate_dir, paste0("gating.png")),
    width = n_plots * 300, 
    height = length(files) * 300)
layout(matrix(1:(length(files) * n_plots), 
              ncol = n_plots, 
              byrow = TRUE))

for (file in files){
  
  ff <- flowCore::read.FCS(filename = file, 
                           transformation = FALSE)
  
  ff <- gate_intact_cells(flow_frame = ff, 
                          file_name = basename(file))
  
  ff <- gate_singlet_cells(flow_frame = ff,
                           channels = "Event_length",
                           file_name = basename(file))
  
  ff <- gate_live_cells(flow_frame = ff, 
                        viability_channel = "Pt195Di",
                        out_dir = gate_dir)
  
  flowCore::write.FCS(ff, file.path(gate_dir,
                                    gsub(".fcs", "_gated.fcs", basename(file))))
}

dev.off()

# ------------------------------------------------------------------------------
# Normalization using reference sample -----------------------------------------
#-------------------------------------------------------------------------------

# Set input directory 
gate_dir <- file.path(dir, "Gated")

# Define reference samples
files_ref <- list.files(gate_dir, 
                        pattern = "REF.*_gated.fcs$", 
                        full.names = TRUE, 
                        recursive = T)

# Define batch labels for each files
labels_ref <- stringr::str_match(basename(files_ref), 
                                   ".*_(day[0-9]*).*.fcs")[,2]
# Define channels to be normalized
ff <- read.FCS(files_ref[1])
channels <- grep("Pd|Rh|140", 
                 grep("Di", colnames(ff), 
                      value = T), 
                 value = T, invert = T) 

# Define out_dir for normalized files
norm_dir <- file.path(dir, "CytoNormed")
if(!dir.exists(norm_dir))(dir.create(norm_dir))

# Build the normalization model using reference samples and plot quantiles 
png(file.path(norm_dir, "005_095_normalization.png"),
    width = length(channels) * 300,
    height = (length(files_ref) * 2 + 1) * 300)
model <- CytoNorm::QuantileNorm.train(files = files_ref,
                                      labels = labels_ref, 
                                      channels = channels, 
                                      transformList = transformList(channels, 
                                                                    CytoNorm::cytofTransform), 
                                      nQ = 2, 
                                      limit = c(0,8),
                                      quantileValues = c(0.05, 0.95), 
                                      goal = "mean", 
                                      plot = TRUE)
dev.off()

# save the model
saveRDS(object = model, 
        file = file.path(norm_dir, "005_095_model.RDS"))

# Define path to the files for normalization
files <- list.files(file.path(gate_dir), 
                    pattern = "_gated.fcs$", 
                    full.names = TRUE, recursive = T)

# Define batch labels for each files, note that they need to corresponds to 
# reference labels 
labels <- stringr::str_match(basename(files), 
                             ".*_(day[0-9]*).*.fcs")[,2]

# Normalize files 
CytoNorm::QuantileNorm.normalize(model = model, 
                                 files = files, 
                                 labels = labels, 
                                 transformList = transformList(channels, 
                                                               CytoNorm::cytofTransform),
                                 transformList.reverse = transformList(channels, 
                                                                       CytoNorm::cytofTransform.reverse), 
                                 outputDir = norm_dir)

# ------------------------------------------------------------------------------
# Plot batch effect ------------------------------------------------------------
#-------------------------------------------------------------------------------

# Define batch pattern
batch_pattern <- "day[0-9]*"

# Define files before normalization and order them according to the batch 
# Set input directory for files before CytoNorm normalization 
gate_dir <- file.path(dir, "Gated")
files_before_norm <- list.files(gate_dir, 
                                pattern = ".fcs", 
                                full.names = T)
batch <- str_match(files_before_norm, "day[0-9]*")[,1]
files_before_norm <- files_before_norm[order(factor(batch))]

# Define files after normalization and order them according to the batch 
# Set input directory for files after CytoNorm normalization
norm_dir <- file.path(dir, "CytoNormed")
files_after_norm <- list.files(norm_dir, 
                               pattern = ".fcs", 
                               full.names = T)
batch <- stringr::str_match(files_after_norm, "day[0-9]*")[,1]
files_after_norm <- files_after_norm[order(factor(batch))]

# Plot batch effect 
plot_batch(files_before_norm = files_before_norm, 
           files_after_norm = files_after_norm,
           out_dir = norm_dir, 
           batch_pattern = batch_pattern, 
           clustering_markers = c("CD", "IgD", "HLA"),
           manual_colors = c("darkorchid4", "darkorange", "chartreuse4"))

batch_pattern <- "day[0-9]*"
plot_marker_quantiles(files_after_norm = files_after_norm, 
                      files_before_norm = files_before_norm, 
                      batch_pattern = batch_pattern, 
                      arcsine_transform = TRUE,
                      markers_to_plot = c("CD", "HLA", "IgD", "IL", "TNF",
                                          "TGF", "GR", "IFNa", "MCP", "MIP"),
                      manual_colors = c("darkorchid4", "darkorange", "darkgreen"),
                      out_dir = norm_dir)

# Use FlowSOM to extract cell frequency and MSI

# Create a list with files before and after normalization
all_files <- list("before" = files_before_norm,
                 "after" = files_after_norm)

# perform FlowSOM clustering and extract cell frequency and msi per cluster and metacluster
mx <- extract_pctgs_msi_per_flowsom(file_list = all_files, 
                                    nCells = 50000, 
                                    phenotyping_markers =  c("CD", "HLA", "IgD"),
                                    functional_markers = c("MIP", "MCP", "IL", "IFNa", "TNF", "TGF", "Gr"),
                                    xdim = 10, 
                                    ydim = 10, 
                                    n_metaclusters = 35, 
                                    out_dir = norm_dir, 
                                    arcsine_transform = TRUE)

# perform dimensional reduction and plot the data 
before <- mx[["before"]]
after <- mx[["after"]]

# set title for each plot
title_gg <- c("Clusters frequency" = "cl_pctgs",
              "Metaclusters frequency" = "mcl_pctgs",
              "Clusters MSI" = "mfi_cl",
              "Metaclusters MSI" = "mfi_mcl")

# create the list to store the plots
plots <- list()
for (name in c("cl_pctgs","mcl_pctgs","mfi_cl","mfi_mcl")){
  print(name)
  
  # set the number of closest neighborhoods 
  n <- 14
  
  # process files before normalization
  df_b <- before[[name]]
  
  # select the columns for which MSI is higher than 0.2 SD
  if(grepl("mfi", name)){
    id_cols <-  which(apply(df_b, 2, sd) > 0.2)
    df_b <- df_b[,id_cols]
  }
  
  # build UMAP for files before the normalization
  set.seed(654)
  df_b_umap <- data.frame(umap(df_b, n_neighbors = n, scale = T))
  
  # process files after normalization
  df_a <- after[[name]]
  
  # select the columns for which MSI is higher than 0.2 SD
  if(grepl("mfi", name)){
    id_cols <-  which(apply(df_a, 2, sd) > 0.2)
    df_a <- df_a[,id_cols]
  }
  
  # build UMAP for files after the normalization
  set.seed(654)
  df_a_umap <- data.frame(umap(df_a, n_neighbors = n, scale = T))
  
  # extract rownames to use the for ggplot annotation
  rnmes <- c(rownames(df_b), rownames(df_a))
  
  #join two UMAP data frames
  dr <- data.frame(rbind(df_b_umap, df_a_umap), check.names = F)
  
  # extract annotation for plotting
  colnames(dr) <- c("dim1", "dim2")
  dr$day <- str_match(rnmes, "day[0-9]*")[,1]
  dr$reference <- ifelse(grepl("REF", rnmes),"ref", "")
  dr$sample <- ifelse(grepl("p1", rnmes),"p1", ifelse(grepl("p2", rnmes), "p2", "ref"))
  dr$stimulation <- str_match(rnmes, "RSQ|LPS|IMQ|CPG|UNS")
  dr$normalization <- ifelse(grepl("Norm", rnmes),"Normalized", "Raw")
  dr$normalization <- factor(dr$normalization, levels = c("Raw", "Normalized"))
  
  # plot
  gg <- ggplot(dr, aes(x = dim1, y = dim2))+
    geom_point(data=dr, aes_string(x="dim1", y="dim2", fill = "day", shape = "sample", color = "day"), 
               size = 3)+
    facet_wrap(~normalization)+
    ggtitle(names(title_gg)[which(title_gg%in% name)])+
    scale_shape_manual(values = c(22, 21, 24))+
    scale_fill_manual(values = c("darkorchid4", "darkorange", "chartreuse4"))+
    scale_color_manual(values = c("darkorchid4", "darkorange", "chartreuse4"))+
    theme(panel.background = element_rect(fill = "white", colour = "black",
                                          size = 1, linetype = "solid"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "right",
          legend.key=element_blank(),
          title = element_text(size = 10),
          strip.text = element_blank(), 
          strip.background = element_rect(fill = "white", colour = "white"))
  
  gg <- ggplotGrob(gg)
  plots[[name]] <- gg
}
gg_a <- ggarrange(plotlist = plots,
                  ncol = 2,
                  nrow = 2)

ggplot2::ggsave(filename ="batch_effect_frequency_MSI.png",
                device = "png",
                path = norm_dir,
                plot = gg_a,
                units = "cm",
                width = 19,
                height = 10, dpi = 300)

# ------------------------------------------------------------------------------
# Run UMAP ---------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Set input directory 
norm_dir <- file.path(dir, "CytoNormed")

# Set output directory 
analysis_dir <- file.path(dir, "Analysis")

if (!dir.exists(analysis_dir)) { 
  dir.create(analysis_dir)
}

files <- list.files(norm_dir, 
                    pattern = ".fcs$", 
                    full.names = TRUE)
batch_pattern <- stringr::str_match(basename(files), ".*(day[0-9]*).*.fcs")[,2]

# Build UMAP on aggregated files
UMAP_res <- UMAP(fcs_files = files, 
                 clustering_markers = c("CD", "HLA", "IgD"),
                 functional_markers = c("IL", "TNF", "TGF", "Gr", "IFNa", "MIP", "MCP1"),
                 out_dir = analysis_dir,
                 batch_pattern = "day[0-9]*", 
                 arcsine_transform = TRUE, 
                 cells_total = 5000)

saveRDS(UMAP_res, 
        file.path(analysis_dir, "UMAP.RDS"))

# get manual labels for UMAP

# copy "gating_strategy.wsp" file from RawFiles to Analysis folder so the workspace 
# is in the same place as generated aggregated file called: "aggregated_for_UMAP_analysis.fcs"

file.copy(file.path(dir, "RawFiles","gating_strategy.wsp"), analysis_dir)

# open flowJo workspace
wsp <- CytoML::open_flowjo_xml(
  file.path(analysis_dir, paste0("gating_strategy.wsp")))

# parse the flowJo workspace
gates <- CytoML::flowjo_to_gatingset(wsp,
                                     "All Samples",
                                     sampNloc = "sampleNode")
# get cell population gate names
gate_names <- flowWorkspace::gs_get_pop_paths(gates, path = "auto")

# select which gates to plot
celltypes_of_interest <- gate_names[-c(1, 3, 10, 14, 15, 18, 19, 20, 26, 28)]

# get gating matrix for analyzed cells
gatingMatrix <- flowWorkspace::gh_pop_get_indices_mat(gates, gate_names)

# get the vector of the cell names for the sellected gates
gating_labels <- manual_labels(gatingMatrix, celltypes_of_interest)

# Add gating names column to be able to plot them on UMAP
UMAP_res$cell_labels <- gating_labels

# Plot UMAP between two donors using RSQ stimulation and day 1 batch
# Add additional column to be able to do facet in ggplot
UMAP_res$indyvidual <- stringr::str_match(UMAP_res$sample_name, "p1|p2|REF")[,1]

# filter UMAP data frame for the needed data
df <- UMAP_res %>% 
  dplyr::filter(batch == "day1" & 
           indyvidual %in% c("p1", "p2") &
           sample_name %in% grep("RSQ", UMAP_res[,"sample_name"], 
                                 value = TRUE))
# select markers to plot 
# markers_to_plot <- grep("CD|HLA|IgD",
#                            colnames(df), value = TRUE)

markers_to_plot <- grep("TNF|MIP1",
                        colnames(df), value = TRUE)

# scale markers between 0-1 for easier visualization 

df[,markers_to_plot] <- apply(df[,markers_to_plot], 2, rescale)

# create a list to store the single plots
plots <- list()

# plot the map for each markers
for(m in markers_to_plot){
  
  set.seed(20)
  p <- df %>% dplyr::sample_n(size = 5000) %>%
    ggplot(aes_string(x = "dim1", y = "dim2", color = m)) +
    geom_point(size = 1.5) +
    ggtitle(m)+
    facet_wrap("indyvidual") +
    scale_color_gradientn(markers_to_plot[markers_to_plot == m], 
                          colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))  +
    theme(panel.background = element_rect(fill = "white", colour = "black",
                                          size = 1, linetype = "solid"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(color="black", size=10, hjust = 0.5, face = "bold", vjust = -1),
          # axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          strip.text.x = element_blank(),
          strip.background = element_rect(fill = "white"), 
          legend.text = element_text(size = 16), 
          legend.title = element_text(size = 20),
          legend.position = "none", 
          plot.margin = unit(c(0, 0.5, 0, 0.5), "lines")) 
  
  plots[[m]] <- p
}

# print the plots and save them 
gg_a <- ggpubr::ggarrange(plotlist = plots, 
                          ncol = 2, #1, 3
                          nrow = 1) #2, 8
ggplot2::ggsave(filename = "marker_expressions_in_UMAP_cyt.png", device = "png", #marker_expressions_in_UMAP_cyt.png
                path = analysis_dir, 
                plot = gg_a,
                width = 18, #18
                height = 5, # 21
                units = "cm", 
                dpi = 300)

# select number of colors equal to the number of cell populations
n <- length(unique(df$cell_labels))
myColors <- pals::glasbey(n = n)
#  assign population names to the color
names(myColors) <- unique(df$cell_labels)
#  for better visualization change color for Unknown population to white
myColors["Unknown"] <- "white"
#  set seed to get reproducible results
set.seed(21)
# plot manual labels on UMAP, subset the number of cells for easier interpretation 
df %>% dplyr::sample_n(size = 8000) %>%
  ggplot(aes_string(x = "dim1", y = "dim2", fill = "cell_labels")) +
  geom_point(size = 2, pch = 21) +
  scale_fill_manual(values = myColors)+
  guides(fill = guide_legend(ncol = 1)) +
  theme(panel.background = element_rect(fill = "white", colour = "black",
                                        size = 1, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.subtitle = element_text(color="black", size=16, hjust = 0.95, face = "bold"),
        axis.text = element_text(size = 10, colour = "black"),
        axis.title = element_blank(), 
        legend.text = element_text(size = 8), 
        legend.title = element_blank(),
        legend.position = "right", 
        legend.key = element_blank(), 
        legend.spacing.y = unit(0.3, 'cm'),
        legend.key.size = unit(0.25, "cm")) 

# save the plot 
ggsave(filename = "manual_labels_UMAP.png", device = "png", 
       path = analysis_dir, width = 12, height = 7.5, units = "cm", dpi = 300)


# plot density 
df <- UMAP_res %>% 
  dplyr::filter(batch == "day1" & 
                  indyvidual %in% c("p1", "p2") &
                  sample_name %in% grep("RSQ", UMAP_res[,"sample_name"], 
                                        value = TRUE))

set.seed(20)
gg_dens <-  df %>%
  ggplot(aes_string(x = "dim1", y = "dim2")) +
  geom_point(shape=16, size=0.25, show.legend = FALSE) +
  stat_density_2d(geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)), colour = "black",
                  bins = 10, show.legend = FALSE)+
  scale_fill_distiller(palette = "Spectral", direction = -1)+
  facet_wrap("indyvidual", nrow = 1) +
  xlim(-10, NA)+
  theme(panel.background = element_rect(fill = "white", colour = "black",
                                        size = 1, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        strip.text.x = element_blank(),
        strip.background = element_rect(fill = "white"), 
        legend.position = "none", 
        plot.margin = unit(c(0, 0.5, 0, 0.5), "lines"))

ggplot2::ggsave(filename = "density_plot_UMAP.png", device = "png", #marker_expressions_in_UMAP_cyt.png
                path = analysis_dir, 
                plot = gg_dens,
                width = 9, #18
                height = 4.6, # 21
                units = "cm", 
                dpi = 300)








