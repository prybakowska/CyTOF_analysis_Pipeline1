#' @description Finds all the mass channels
#' @param flow_frame Untranfosrmed flow frame
#' @param channels Pattern for non-mass channels
#' @param ... Additional arguments to pass to grep
#' @return Logical vector with TRUE values for mass channels


check <- function(x) tryCatch(if(class(x) == 'logical') 1 else 1, 
                              error=function(e) 0)


find_mass_ch <- function(flow_frame, 
                         channels = "Time|Event_length|Center|Offset|Width|Residual",
                         ...){
  non_mass_ch <- grep(c(channels), 
       colnames(flow_frame), 
       invert = TRUE, ...)
  return(non_mass_ch)
}


#' @description Cleans the flow rate using functions from flowAI package.
#' @param flow_frame Flow frame.
#' @param to_plot Logical if to plot cleaning results.
#' @param out_dir Character, pathway to where the plots should be saved, 
#' only if to_plot = TRUE.
#' @return Cleand, untransformed flow frame

clean_flow_rate <- function(flow_frame, to_plot = TRUE, 
                            out_dir = getwd()) {
  
  flow_frame@exprs[, "Time"] <- flow_frame@exprs[, "Time"]/100
  
  FlowRateData <- flowAI:::flow_rate_bin(flow_frame, 
                                         timeCh = "Time", 
                                         timestep = 0.01)
  
  FlowRateQC <- flowAI:::flow_rate_check(flow_frame, 
                                         FlowRateData, 
                                         alpha = 0.01, 
                                         use_decomp = TRUE) 
  
  if (to_plot == TRUE){
    
    out_dir <- file.path(out_dir, "FlowRateCleaning")
    if(!dir.exists(out_dir)){
      dir.create(out_dir)
    }
    
    png(file.path(out_dir,
                  gsub(".fcs", "_flowAI.png", 
                       basename(flow_frame@description$FILENAME))),
        width = 800,
        height = 600)
    p <- flowAI:::flow_rate_plot(FlowRateQC)
    print(p)
    dev.off()
  }
  
  flow_frame_cl <- flow_frame[FlowRateQC$goodCellIDs,]
  flow_frame_cl@exprs[,"Time"] <- flow_frame_cl@exprs[,"Time"]*100

  return(flow_frame_cl)
}

#' @description Cleans the flow rate using functions from flowAI package
#' @param flow_frame Flow frame
#' @param channels_to_clean Character vector of the channels that needs to be cleaned
#' @param to_plot Characer variable that indicates if plots should be genarated.
#' The default is "All", which generates plot for all channels. Other options are
#' "Flagged Only", plots the channels that were spotted with flowcut as incorrect
#' and "None", does not plots anything.
#' @param out_dir Character, pathway to where the plots should be saved, 
#' only if argument to_plot = TRUE
#' @param arcsine_transform Logical, if the data should be transformed with 
#' arcsine transformation and cofactor 5.
#' @param ... Additional arguments to pass to flowcut.
#' @return Cleand, untransformed flow frame.

clean_signal <- function(flow_frame, channels_to_clean = NULL, to_plot = "All", 
                         out_dir = getwd(), arcsine_transform = TRUE,  ...){
  
  channels_to_transform <- find_mass_ch(flow_frame, value = FALSE)
  if (arcsine_transform == TRUE){
    ff_t <- transform(flow_frame, 
                      transformList(colnames(flow_frame)[channels_to_transform], 
                                    CytoNorm::cytofTransform))
  } else {
    ff_t <- flow_frame
  }
  
  if (!is.null(channels_to_clean)){
    
    ch_to_clean <- which(colnames(flow_frame) %in% channels_to_clean)
    
    if(!("TIME" %in% toupper(colnames(flow_frame)[ch_to_clean]))){
      
      ind_Time <- grep("TIME", toupper(colnames(flow_frame)))
      channels <- unique(sort(c(ch_to_clean, ind_Time)))
    }
    
  } else {
    
    ind_Time <- grep("TIME", toupper(colnames(flow_frame)))
    ch_to_clean <- c(ind_Time, channels_to_transform)
    ind_140 <- grep("140", toupper(colnames(flow_frame)))
    channels <- ch_to_clean[ch_to_clean != 14]
  
    }
  
  out_dir <- file.path(out_dir, "SignalCleaning")
  if(!dir.exists(out_dir)){
    dir.create(out_dir)
  }
  
  cleaned_data <- flowCut(ff_t, 
                          Segment = 1000, 
                          MaxPercCut = 0.5, 
                          Channels = channels,
                          FileID = gsub("_beadNorm", "_flowCutCleaned", 
                                        basename(flow_frame@description$FILENAME)),
                          Plot = to_plot,
                          Directory = out_dir,
                          UseOnlyWorstChannels = TRUE,
                          AllowFlaggedRerun = TRUE,
                          AlwaysClean = TRUE,
                          ...)
  
  ff_t_clean <- cleaned_data$frame
  
  if (arcsine_transform == TRUE){
    ff_clean <- transform(ff_t_clean, 
                          transformList(colnames(flow_frame)[channels_to_transform], 
                                        cytofTransform.reverse))
  } else {
    ff_clean <- ff_t_clean
  }
  
  return(ff_clean)
}

#' @description Creates the reference flowframe which bead´s values will be used as reference 
#'for other files during the normalization.
#' @param fcs_files Path to selected fcs file which will be used as reference 
#'file upon normalization.
#' @param beads The same as in CATALYAST package, character variable:
#'"dvs" (for bead masses 140, 151, 153 ,165, 175) or 
#'"beta" (for bead masses 139, 141, 159, 169, 175) or a numeric vector of masses.

create_ref <- function(fcs_file, beads = NULL){
  
  if (!file.exists(fcs_file[1])){
    stop("incorrect file path, the fcs file does not exist")
  }
  
  print(paste("preparing reference sample from ", 
              basename(fcs_file)))
  dat <- prepData(fcs_file) 
  dat_norm <- normCytof(x = dat,
                        beads = beads,
                        remove_beads = TRUE,
                        norm_to = NULL,
                        k = 80,
                        plot = FALSE, verbose = FALSE)
  ff_ref <- sce2fcs(dat_norm$beads)
  return(ff_ref)
}


#' @description Performs bead normalization usinf function from CATALYST package.
#' @param flow_frame Untranfosrmed flow frame.
#' @param markers_to_keep Character vector, marker names to be plotted, can be 
#' full marker name e.g. "CD45" or "CD" if all CD-markers needs to be plotted
#' @param to_plot Logical if to plot bead gate and bead normalization for each 
#' file.
#' @param out_dir Pathway to where the plots should be saved, 
#' only if to_plot = TRUE, the default is working directory.
#' @return Bead normalized flow frame without the beads.

bead_normalize <- function(flow_frame, keep_all_markers = TRUE, 
                           markers_to_keep = c("PD","CD", "HLA", "IgD", "TCR", 
                                               "BAFF", "Ir", "Pt195"), 
                           beads = "dvs", norm_to_ref = NULL, 
                           remove_beads = TRUE, to_plot = TRUE, 
                           out_dir = getwd(), ...){
  
  if (keep_all_markers == FALSE){
    
    if (is.null(markers_to_keep)){
      stop("pattern for markers needs to be specify")
    }
    
    matches <- paste(markers_to_keep, collapse="|")
    
    m_to_keep <- grep(matches, get_markers(flow_frame, colnames(flow_frame)), 
                             ignore.case = TRUE, value = FALSE)
    
    non_mass_ch <- grep("Time|length|Ce140|Center|
                        Offset|Width|Residual", colnames(ff),
                        ignore.case = TRUE, value = FALSE)
    
    channels_to_keep <- c(m_to_keep, non_mass_ch)
    channels_to_keep <- colnames(flow_frame)[sort(unique(channels_to_keep))]
    
    flow_frame <- flow_frame[, channels_to_keep]
  }
  
  if (is.null(norm_to_ref)){
    warning("the reference file is not defined. Each file will be normalized 
       to its own bead mean but not across all files")
  }
  
  dat <- prepData(flow_frame) 
  
  # normalize the data and remove beads
  dat_norm <- normCytof(x = dat,
                        beads = beads,
                        remove_beads = remove_beads,
                        norm_to = norm_to_ref,
                        k = 80,
                        plot = TRUE)
 
  # convert back to .fcs files and save 
  ff <- sce2fcs(dat_norm$data)
  
  filename <- basename(flow_frame@description$FILENAME)
  
  if (to_plot == TRUE){
    
    # plot and save diagnostic plots 
    plot_dir <- file.path(out_dir, "Plots_BeadNormalization")
    if(!dir.exists(plot_dir))(dir.create(plot_dir))
    
    p <- dat_norm$scatter
    ggsave(filename = file.path(plot_dir, gsub(".FCS","_beadGate.png", filename)), 
           plot = p)
   
    p <- dat_norm$lines
    ggsave(filename = file.path(plot_dir, gsub(".FCS","_beadLines.png", filename)), 
           plot = p)
   
  }
  
  return(ff)
}

#' @description Calculates quantiles for selected markers and plots them as 
#' diagnostic plot 
#' @param fcs_files Character, full path to fcs_files.
#' @param file_batch_id Character vector, batch label for each fcs_file, 
#'the order needs to be the same as in fcs_files.
#' @param file_id Character vector, file label (shorter file name) for each 
#'fcs_file, the order needs to be the same as in fcs_files.
#' @param arcsine_transform Logical, if the data should be transformed with 
#' arcsine transformation and cofactor 5.
#' @param markers_to_plot character vector, marker names to be plotted, can be 
#' full marker name e.g. "CD45" or "CD" if all CD-markers needs to be plotted  

# Plot diagnostic plot for markers across each run 
plot_marker_quantiles <- function(fcs_files, file_batch_id, file_id,
                                  arcsine_transform = TRUE, 
                                  markers_to_plot = NULL, out_dir = getwd()){
  
  if (!file.exists(fcs_files[1])){
    stop("incorrect file path, the fcs file does not exist")
  }
  
  if(check(norm_markers) == 0){
    
    o <- capture.output(ff_tmp <- read.FCS(file.path(fcs_files[1])))
    
    if (!is.null(markers_to_plot)){
      
      if(!is.character(markers_to_plot)){
        stop ("markers are not a character vector")
      }
      
      matches <- paste(markers_to_plot, collapse="|")
      
      norm_markers <- grep(matches,
                           get_markers(ff_tmp, find_mass_ch(ff_tmp, 
                                                            value = TRUE)), 
                           value = TRUE)
    }
    norm_markers <- find_mass_ch(ff_tmp, value = TRUE)
    norm_markers <- get_markers(ff_tmp, norm_markers)
  }
  
  quantile_values <-  c(0.05, 0.25, 0.5, 0.75, 0.95)
  
  files <- fcs_files
  
  quantiles <- expand.grid(File = files,
                           Marker = norm_markers,
                           Quantile = quantile_values,
                           Value = NA)
  for (file in files) {
    print(file)
    
    o <- capture.output(ff <- read.FCS(file.path(file)))   
    
    if(arcsine_transform == TRUE){
      ff <- transform(ff, transformList(grep("Di", colnames(ff), value = TRUE),
                                        arcsinhTransform(a = 0, b = 1/5, c = 0)))
    }
    
    for (marker in names(norm_markers)) {
      quantiles_res <- quantile(exprs(ff)[, marker],
                                quantile_values)
      for (i in seq_along(quantiles_res)) {
        quantiles <- quantiles %>%
          dplyr::mutate(Value = replace(Value, 
                                        File == file & 
                                          Marker == norm_markers[marker] &
                                          Quantile == quantile_values[i],
                                        quantiles_res[i]))
        
      }
    }
  }
  
  quantiles$Batch <- file_batch_id
  quantiles$Sample <- file_id
  
  p <- quantiles %>%
    ggplot(aes(x = Sample,
               y = Value,
               color = Batch)) +
    geom_point(aes(alpha = ifelse(Quantile == "0.5", 1, 0))) +
    geom_line(aes(alpha = ifelse(Quantile != "0.5", 1, 0),
                  group = interaction(Batch, Quantile),
                  size = ifelse(Quantile %in% c("0.05", "0.95"), 0.5,
                                ifelse(Quantile == "0.5", 0, 1))),
              alpha = 0.5) +
    scale_size_identity() +
    scale_alpha_identity() +
    facet_wrap( ~ Marker + Batch, ncol = 2) +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "bottom")
  # legend.justification=c(1,0), legend.position=c(1,-0.15)) +
  # ggtitle(paste("Marker distribution across aliquots for each batch"))
  
  ggsave(filename = paste("Marker distribution across aliquots and batches.pdf"), 
         plot = p, 
         path = file.path(out_dir), 
         width = 5, height = 50, limitsize = F)
  
}


fsom_aof <- function(fcs_files, 
                     phenotyping_markers,
                     nCells = length(fcs_files)*10000,
                     xdim = 10,
                     ydim = 10,
                     nClus = 10,
                     out_dir, 
                     pattern = NULL, ...){
  
  
  if(check(phenotyping_channels) == 0){
    o <- capture.output(ff_tmp <- read.FCS(file.path(files[1])))
    markers <- get_markers(ff_tmp, colnames(ff_tmp))
    phenotyping_channels <- grep(paste(phenotyping_markers, 
                                       collapse = ("|")), markers, value = TRUE)
  }
  
  fsom <- prepareFlowSOM(file = fcs_files,
                         colsToUse = names(phenotyping_channels),
                         seed = 1, 
                         plot = FALSE,
                         nCells = nCells,
                         transformList = transformList(names(phenotyping_channels), 
                                                       cytofTransform),
                         FlowSOM.params = list(xdim = xdim, 
                                               ydim = ydim, 
                                               nClus = nClus, 
                                               scale = FALSE))
  
  myCol <- c("tomato", "violet", "grey50", "slateblue1", "yellow","turquoise2",
             "yellowgreen", "skyblue", "wheat2","steelblue", "blue2", "navyblue",
             "orange", "violetred", "red4", "springgreen2",  "palegreen4",
             "tan", "tan2", "tan3", "brown", "grey70", "grey30")
  
  pdf(file.path(out_dir, paste0("FlowSOM_clustering.pdf")))
  PlotStars(fsom$FlowSOM, main = "FlowSOM clustering", 
            backgroundValues = fsom$metaclustering, 
            backgroundColor = myCol)
  dev.off()
  
  return(fsom)
  
}

aof_scoring <- function(fcs_files = files,
                        phenotyping_channels,
                        fsom,
                        out_dir,
                        list_for_scores = NULL){
  
  aof_scores <- matrix(NA,
                       nrow = length(fcs_files), 
                       ncol = length(phenotyping_channels),
                       dimnames = list(fcs_files,  
                                       names(phenotyping_channels)))
  
  for(file in fcs_files){ 
    print(paste("calculating AOF", file))
    File_ID <- which(fcs_files == file)
    idx <- which(fsom$FlowSOM$data[,"File"] == File_ID)
    fcs_data <- fsom$FlowSOM$data[idx,]
    MC <- fsom$metaclustering[fsom$FlowSOM$map$mapping[idx, 1]]
    
    aof_tmp <- greedyCytometryAof(fcs_data = fcs_data, 
                                  y = MC, 
                                  channel_names = names(phenotyping_channels),
                                  width = 0.05, 
                                  cofactor = NULL, 
                                  verbose = TRUE) 
    aof_scores[file, ] <- aof_tmp$Aof
  }
  
  aof_scores_scaled <- scale(aof_scores)
  
  # 
  # aof_scores_scaled <- apply(aof_scores_scaled,2, function(x){
  #   missing <- which(is.na(x))
  #   x[missing] <- 0
  #   x
  # })
  
  aof_scores_scaled <- pmax(aof_scores_scaled, 0)^2
  sample_scores <- apply(aof_scores_scaled, 1, sum, na.rm = TRUE)
  
  df <- as.data.frame(sample_scores)
  
  # list_scores <- list("aof_scores" = aof_scores,
  #                     "aof_scores_scaled" = aof_scores_scaled)
  # 
  # # aof <- c("aof_scores", "aof_scores_scaled")
  # print(paste("plotting AOF scores for "))
  # 
  # for (name in names(list_scores)) {
  #   
  #   filename <- file.path(out_dir, paste0(name, ".pdf"))
  #   pheatmap(list_scores[[name]],
  #            #scale = "column", 
  #            cluster_rows = FALSE, 
  #            cluster_cols = FALSE,
  #            display_numbers = TRUE,
  #            #display_numbers = round(aof_scores, 3),
  #            labels_col = phenotyping_channels,
  #            labels_row = basename(rownames(list_scores[[name]])),
  #            filename = filename, 
  #            number_format = "%.1f", 
  #            width = 10)
  
  
  
  #saveRDS(list_for_scores, file.path(output_dir, "aof_scores_all.RDS"))
  return(df)
  
}










  