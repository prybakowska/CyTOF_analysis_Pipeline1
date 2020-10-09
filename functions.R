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
#' only if argument to_plot = TRUE, default is set to working directory.
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
#' only if argument to_plot = TRUE, default is set to working directory.
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

create_ref <- function(fcs_files, beads = NULL, to_plot = TRUE, 
                       out_dir = getw(), k = 80){
  
  # if (!file.exists(fcs_file[1])){
  #   stop("incorrect file path, the fcs file does not exist")
  # }
  
  # print(paste("preparing reference sample from ", 
  #             basename(fcs_file)))
  # dat <- prepData(fcs_file) 
  # dat_norm <- normCytof(x = dat,
  #                       beads = beads,
  #                       remove_beads = TRUE,
  #                       norm_to = NULL,
  #                       k = 80,
  #                       plot = FALSE, 
  #                       verbose = FALSE, 
  #                       transform = FALSE)
  # ff_ref <- sce2fcs(dat_norm$beads)
  # return(ff_ref)
  
  ff <- AggregateFlowFrames(fileNames = fcs_files, cTotal = length(fcs_files)*25000)
  
  dat <- prepData(ff) 
  
  dat_norm <- normCytof(x = dat,
                        beads = beads,
                        remove_beads = TRUE,
                        norm_to = NULL,
                        k = k,
                        plot = to_plot, 
                        verbose = FALSE, 
                        transform = FALSE)
  ff_ref <- sce2fcs(dat_norm$beads)
  rm(ff)
  
  if (to_plot == TRUE){
    
    # plot and save diagnostic plots 
    plot_dir <- file.path(out_dir, "Plots_BeadNormalization")
    if(!dir.exists(plot_dir))(dir.create(plot_dir))
    
    p <- dat_norm$scatter
    ggsave(filename = file.path(plot_dir,"RefBeadGate.png"), 
           plot = p, limitsize = FALSE)
    
    p <- dat_norm$lines
    ggsave(filename = file.path(plot_dir,"RefBeadLines.png"), 
           plot = p, limitsize = FALSE)
  }
  return(ff_ref)
}


#' @description Performs bead normalization usinf function from CATALYST package.
#' @param flow_frame Untranfosrmed flow frame.
#' @param markers_to_keep Character vector, marker names to be plotted, can be 
#' full marker name e.g. "CD45" or "CD" if all CD-markers needs to be plotted
#' @param to_plot Logical if to plot bead gate and bead normalization for each 
#' file.
#' @param out_dir Character, pathway to where the plots should be saved, 
#' only if argument to_plot = TRUE, default is set to working directory.
#' @return Bead normalized flow frame without the beads.

bead_normalize <- function(flow_frame, keep_all_markers = TRUE, 
                           markers_to_keep = NULL, 
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
    
    non_mass_ch <- grep("Time|length|Ce140|151|153|165|175|Center|Offset|Width|
                        |Residual", 
                        colnames(flow_frame),
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
                        plot = TRUE, 
                        transform = FALSE)
 
  # convert back to .fcs files and save 
  f <- sce2fcs(dat_norm$data)
  
  filename <- basename(flow_frame@description$FILENAME)
  
  if (to_plot == TRUE){
    
    # plot and save diagnostic plots 
    plot_dir <- file.path(out_dir, "Plots_BeadNormalization")
    if(!dir.exists(plot_dir))(dir.create(plot_dir))
    
    p <- dat_norm$scatter
    ggsave(filename = file.path(plot_dir, gsub(".FCS","_beadGate.png", filename)), 
           plot = p, limitsize = FALSE)
   
    p <- dat_norm$lines
    ggsave(filename = file.path(plot_dir, gsub(".FCS","_beadLines.png", filename)), 
           plot = p, limitsize = FALSE)
   
  }
  
  return(f)
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
plot_marker_quantiles <- function(files_after_norm, files_before_norm, 
                                  batch_pattern,
                                  file_id = NULL, 
                                  arcsine_transform = TRUE, 
                                  markers_to_plot = NULL, out_dir = getwd()){
  
  
  fcs_files <- c(files_after_norm, files_before_norm)
  tmp <- c(paste0(files_after_norm, "_YES"), paste0(files_before_norm, "_NO"))
  
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
    } else {
      # TODO CHEK if else should be here
      norm_markers <- find_mass_ch(ff_tmp, value = TRUE)
      norm_markers <- get_markers(ff_tmp, norm_markers)
    }
    
  }
  
  quantile_values <-  c(0.05, 0.25, 0.5, 0.75, 0.95)
  quantiles <- expand.grid(File = tmp,
                           Marker = norm_markers,
                           Quantile = quantile_values,
                           Value = NA)
  quantiles <- cbind(quantiles, "Batch" = stringr::str_match(
    basename(as.character(quantiles$File)), batch_pattern)[,1])
  
  quantiles$Normalization <- gsub(".*.fcs_|.*.FCS_", "", quantiles$File) 
  quantiles$File <- gsub("_YES|_NO", "", quantiles$File)
  
  for (file in fcs_files) {
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
  
  # quantiles$Batch <- file_batch_id
  quantiles$Sample <- gsub("Norm_", "", 
                           gsub("_CC_gated.fcs|_gated.fcs|_gated.fcs|_beadNorm.fcs|.FCS",
                                "", basename(as.character(quantiles$File))))
  
  ncols <- length(unique(quantiles$Batch))
  p <- quantiles %>% filter(Normalization == "YES") %>%
    ggplot(aes(x = Sample,
               y = Value,
               color = Batch)) +
    geom_point(data = quantiles %>% filter(Normalization == "NO"), 
               aes(alpha = ifelse(Quantile == "0.5", 1, 0)), color = "grey48") +
    geom_line(data = quantiles %>% filter(Normalization == "NO"), 
              aes(alpha = ifelse(Quantile != "0.5", 1, 0),
                  group = interaction(Batch, Quantile),
                  size = ifelse(Quantile %in% c("0.05", "0.95"), 0.5,
                                ifelse(Quantile == "0.5", 0, 1))),
              alpha = 0.5, color = "grey48") +
    
    geom_point(aes(alpha = ifelse(Quantile == "0.5", 1, 0))) +
    geom_line(aes(alpha = ifelse(Quantile != "0.5", 1, 0),
                  group = interaction(Batch, Quantile),
                  size = ifelse(Quantile %in% c("0.05", "0.95"), 0.5,
                                ifelse(Quantile == "0.5", 0, 1))),
              alpha = 0.5) +
    scale_size_identity() +
    scale_alpha_identity() +
    facet_wrap(~ Marker + Batch, ncol = ncols, scales = "free_x") +
    # facet_grid(~ Marker + Batch,  space = "free_x") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "bottom")
  # legend.justification=c(1,0), legend.position=c(1,-0.15)) +
  # ggtitle(paste("Marker distribution across aliquots for each batch"))
  
  ggsave(filename = paste("Marker distribution across aliquots and batches.pdf"), 
         plot = p, 
         path = file.path(out_dir), 
         width = length(fcs_files)*0.25, height = length(norm_markers)*4, limitsize = F)
  
}

#' @description Builds FlowSOM tree for the files scoring 
#' @param fcs_files Character, full path to fcs_files.
#' @param phenotyping_markers Character vector, marker names to be used for phenotyping,
#' can be full marker name e.g. "CD45" or "CD" if all CD-markers needs to be plotted
#' @param nCells Numeric, the total number of cells, to use for FlowSOM clustering. 
#' This number is determined by total number of fcs files, as a defult 1000 cells 
#' is used per file
#' @param xdim Numeric, parameter to pass to FlowSOM, width of the SOM grid
#' @param ydim Numeric, parameter to pass to FlowSOM, geight of the SOM grid
#' @param nClus Numeric, exact number of clusters for metaclustering 
#' @param out_dir Character, pathway to where the FlowSOM clustering plot should
#'  be saved, default is set to working directory.
#' @param pattern ?? check
#' @param arcsine_transform Logical, if the data should be transformed with 
#' arcsine transformation and cofactor 5.

fsom_aof <- function(fcs_files, 
                     phenotyping_markers,
                     nCells = length(fcs_files)*10000,
                     xdim = 10,
                     ydim = 10,
                     nClus = 10,
                     out_dir, 
                     # pattern = NULL,
                     arcsine_transform, ...){
  
  
  if(check(phenotyping_channels) == 0){
    o <- capture.output(ff_tmp <- read.FCS(file.path(files[1])))
    markers <- get_markers(ff_tmp, colnames(ff_tmp))
    phenotyping_channels <- grep(paste(phenotyping_markers, 
                                       collapse = ("|")), markers, value = TRUE)
  }
  
  if(arcsine_transform  == TRUE){
    trans <- transformList(names(phenotyping_channels), cytofTransform)
  } else {
    trans <- NULL
  }
  
  fsom <- prepareFlowSOM(file = fcs_files,
                         colsToUse = names(phenotyping_channels),
                         seed = 1, 
                         plot = FALSE,
                         nCells = nCells,
                         transformList = trans,
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

aof_scoring <- function(fcs_files,
                        phenotyping_markers,
                        fsom,
                        out_dir,
                        list_for_scores = NULL){
  
  if(check(phenotyping_channels) == 0){
    
    o <- capture.output(ff_tmp <- read.FCS(file.path(files[1])))
    markers <- get_markers(ff_tmp, colnames(ff_tmp))
    phenotyping_channels <- grep(paste(phenotyping_markers, 
                                       collapse = ("|")), markers, value = TRUE)
  }
  
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
                                  cofactor = 5, 
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

#' @param out_dir Character, pathway to where the plots should be saved, 
#' only if argument to_plot = TRUE, default is set to working directory

plot_aof_all_files <- function(scores, out_dir, sd) {
  
  df_scores <- do.call(rbind, scores)
  df_scores$file_names <- basename(rownames(df_scores))
  
  scores_median <- median(df_scores$sample_scores)
  scores_MAD <- mad(df_scores$sample_scores)
  
  df_scores$quality <- ifelse(df_scores$sample_scores > 
                                (scores_median + sd * scores_MAD),"bad","good")
  
  bad_scores <- sum(df_scores$quality == "bad")
  
  colors <- c("bad" = "red", "good" = "darkblue", "threshold= " = "orange")  
  
  max_score <- max(df_scores$sample_scores)
  max_pctgs <- max_score + (max_score * 0.1)
  
  p <- ggplot(df_scores, aes(x = file_names, y = sample_scores, color = quality)) +
    geom_point(size = 4) +
    scale_colour_manual(values = colors) +
    ylim(-0.5, max_pctgs) + 
    geom_hline(yintercept = scores_median + sd * scores_MAD, 
               linetype = "dashed", color = "darkgreen", size = 1)+
    # geom_text(aes(x= 0, y= scores_median + sd * scores_MAD, label = 
    #                 paste("threshold= ", round(scores_median + sd * scores_MAD, digits = 2)), 
    #               vjust = -1, hjust = -0.2), color = "darkgreen", size = 6) + 
    annotate(geom="text", x = mean(as.numeric(as.factor(df_scores$file_names))), 
             y= max_score - 0.05*max_score, label=paste("N bad = ", bad_scores),
             color="red", size = 5) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.text = element_text(size = 11),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          # panel.border = fill = "black",
          axis.line = element_line(colour = "black"),
          legend.position = "none") +
    scale_x_discrete(breaks = df_scores$file_names[df_scores$quality == "bad"])
  p
  
  ggsave(filename = "AOF_scores.pdf", plot = p, path = file.path(out_dir))
  
  # bad_files <- sample_scores[which(sample_scores$vector == "red"),]
  
  saveRDS(df_scores, file.path(out_dir, "AOF_sample_scores.RDS"))
  write.csv(df_scores, file = file.path(out_dir, "AOF_excluded_files.csv"))
  
  return(df_scores)
}


#' @description performs sample quality score 
#' @param fcs_files Character, full path to fcs_files.
#' @param file_batch_id Character vector, batch label for each fcs_file, 
#'the order needs to be the same as in fcs_files.
#' @param out_dir Character, pathway to where the plots should be saved, 
#' only if argument to_plot = TRUE, default is set to working directory
#' @param phenotyping_markers Character vector, marker names to be used for 
#' flowsom clustering inlcuding DNA marker Iridium and viability staining if avaiable,
#' can be full marker name e.g. "CD45" or pattern "CD" if 
#' all CD-markers needs to be plotted, default is set to NULL, so all the mass 
#' channels will be used 
#' @param arcsine_transform Logical, if the data should be transformed with 
#' arcsine transformation and cofactor 5.
#' @param sd numeric, number of standard deviation allowed for outlier detection. 
#' Default = 3.

file_quality_check <- function(fcs_files, file_batch_id, 
                               out_dir, phenotyping_markers = NULL, 
                               arcsine_transform = TRUE, sd = 3){
  
  if(!dir.exists(out_dir)) dir.create(out_dir)
  
  scores <- list()
  
  for (batch in unique(file_batch_id)){
    
    files <- fcs_files[file_batch_id == batch]
    
    fsom <- fsom_aof(fcs_files = files, phenotyping_markers = phenotyping_markers, 
                     out_dir = out_dir, arcsine_transform = arcsine_transform)
    
    scores[[batch]] <- aof_scoring(fcs_files = files, phenotyping_markers = phenotyping_markers, 
                                   fsom = fsom, out_dir = out_dir)
    
    final_score <- plot_aof_all_files(scores = scores, out_dir = out_dir, 
                                      sd = sd)
  }
}


#' @description performs sample debarcoding 
#' @param fcs_files Character, full path to fcs_files.
#' @param bad_quality_files Character, full path to the .RDS file that shows 
#' which files should be removed due to bad quality events. Default NULL.
#' @param out_dir Character, pathway to where the plots should be saved, 
#' only if argument to_plot = TRUE, default is set to working directory
#' @param min_threshold if the minimal treshold for barcoding should be applied. 
#' Default is set to 0.18.
#' @param barcodes_used character, the names of the barcodes that were used, eg.
#' barcode 1 is the same as A1. 
#' @param file_batch_id Character vector, batch label for each fcs_file, 
#'the order needs to be the same as in fcs_files.

debarcode_files <- function(fcs_files, 
                            out_dir = getwd(), min_threshold = TRUE, 
                            barcodes_used = NULL, file_batch_id){
  
  if(anyDuplicated(fcs_files) != 0){
    stop("names of fcs files are duplicated")
  }
  
  # if(!is.null(bad_quality_files)){
  #     file_scores <- readRDS(bad_quality_file)
  #     selection <- file_scores$file_names[file_scores$quality == "good"]
  #     
  #     print(paste("following files were discarded from further analysis",
  #          basename(fcs_files[!(basename(fcs_files) %in% selection)])))
  #     
  #     fcs_files_clean <- fcs_files[basename(fcs_files) %in% selection]
  # } else {
  #   fcs_files_clean <- fcs_files
  # }
  
  for (file in fcs_files){
    print(paste0("   ", Sys.time()))
    print(paste0("   Debarcoding ", file))
    ff <- read.FCS(file)
     
    file_id <- which(file == fcs_files)
    batch_id <- file_batch_id[file_id]
    
    if(!dir.exists(out_dir)) dir.create(out_dir)
    
    if(!is.null(barcodes_used)){
      # s_key <- sample_key[rownames(sample_key) %in% barcodes_used,]
      if(is.list(barcodes_used)){
        s_key <- sample_key[rownames(sample_key) %in% barcodes_list[[batch_id]],]
      } else {
        s_key <- sample_key[rownames(sample_key) %in% barcodes_used,]
      }
      
    } else {
      s_key <- sample_key
    }

    dat <- prepData(ff)
    dat <- assignPrelim(dat, s_key)
    rownames(dat)[rowData(dat)$is_bc]
    # table(colData(dat)$bc_id)
    dat <- estCutoffs(dat)
    
    if (min_threshold == TRUE){
      id <- metadata(dat)$sep_cutoffs < 0.18
      metadata(dat)$sep_cutoffs[id] <- 0.18 
    }
    
    id <- is.na(metadata(dat)$sep_cutoffs)
    metadata(dat)$sep_cutoffs[id] <- 1 
  
    p <- plotYields(dat, which = rownames(s_key))
   
    pdf(file.path(out_dir, paste(gsub(".fcs", "_yields.pdf", basename(file)))))
    for (name in names(p)){
      print(p[[name]])
    }
    dev.off()
    
    dat <- applyCutoffs(dat)
    p <- plotEvents(dat)
    
    pdf(file.path(out_dir, paste(gsub(".fcs", "_debarcode_quality.pdf", 
                                      basename(file)))))
    for (name in names(p)){
      print(p[[name]])
    }
    dev.off()
    
    dat <- dat[, dat$bc_id !=0]
    fs <- sce2fcs(dat, split_by = "bc_id")
    
    # for (i in 1:length(fs)){
    # 
    #   ff <- fs[i]
    #   rname <- paste0("_",rownames(ff@phenoData), ".fcs")
    
    tmp_dir <- file.path(out_dir, batch_id)
    if(!dir.exists(tmp_dir)) dir.create(tmp_dir)
    
    file_name <- gsub("_cleaned.fcs|.fcs", "", basename(file))
    
    write.flowSet(fs, outdir = tmp_dir, 
                  filename = paste0(rownames(fs@phenoData), "_", file_name, 
                                    "_debarcoded.fcs")) 
    
    # }
  }
}

#' @description performs sample aggregation
#' @param fcs_files Character, full path to fcs_files.
#' @param channels_to_keep Character vector with channel names to be kept. 
#' Default NULL
#' @param outputFile character, the names of the file that will be given
#' @param maxcells numeric, maximum cells to randomly aggregate from each file

aggregate_files <- function(fcs_files, 
                            channels_to_keep = NULL,
                            outputFile = "aggregate.fcs", maxcells = NULL, 
                            write_agg_file = FALSE){
  
  nFiles <- length(fcs_files)
  flowFrame <- NULL
  
  for (i in seq_len(nFiles)) {
    f <- flowCore::read.FCS(fcs_files[i])
    
    if(!is.null(maxcells)){
      c <- sample(seq_len(nrow(f)), min(nrow(f), maxcells))
      f <- f[c,]
    }
    
    m <- matrix(rep(i, nrow(f)))
    m2 <- m + stats::rnorm(length(m), 0, 0.1)
    m <- cbind(m, m2)
    colnames(m) <- c("File", "File_scattered")
    prev_agg <- length(grep("File[0-9]*$", colnames(f)))
    if (prev_agg > 0) {
      colnames(m) <- paste0(colnames(m), prev_agg + 1)
    }
    if(is.null(channels_to_keep)){
      f <- flowCore::fr_append_cols(f, m)
    } else {
      f <- flowCore::fr_append_cols(f[ , channels_to_keep], m)
    }
    if (is.null(flowFrame)) {
      flowFrame <- f
      flowFrame@description$`$FIL` <- gsub(".*/", "",
                                           outputFile)
      flowFrame@description$FILENAME <- gsub(".*/", "",
                                             outputFile)
    }
    else {
      f@exprs[, "Time"] <- f@exprs[, "Time"] + max(flowFrame@exprs[,"Time"]) + 1000
      flowCore::exprs(flowFrame) <- rbind(flowCore::exprs(flowFrame), 
                                          flowCore::exprs(f))
    }
  }
  flowFrame@description[[paste("flowCore_$P", ncol(flowFrame) - 
                                 1, "Rmin", sep = "")]] <- 0
  flowFrame@description[[paste("flowCore_$P", ncol(flowFrame) - 
                                 1, "Rmax", sep = "")]] <- nFiles + 1
  flowFrame@description[[paste("flowCore_$P", ncol(flowFrame), 
                               "Rmin", sep = "")]] <- 0
  flowFrame@description[[paste("flowCore_$P", ncol(flowFrame), 
                               "Rmax", sep = "")]] <- nFiles + 1
  flowFrame@description[[paste("$P", ncol(flowFrame) - 1, 
                               "B", sep = "")]] <- 32
  flowFrame@description[[paste("$P", ncol(flowFrame), "B", 
                               sep = "")]] <- 32
  flowFrame@description$FIL <- gsub(".*/", "", outputFile)
  
  #write.FCS.corrected(flowFrame, filename = outputFile, endian = "big")
  
  if(write_agg_file == TRUE){
     flowCore::write.FCS(flowFrame, filename = outputFile, endian = "big")
  }

  return(flowFrame)
}



gate <- function(fcs_files, cd45_ch = "Pr141Di", viability_ch = "Pt195Di",
                 out_dir = getwd()) {
  
  n_plots <- 2  
  png(file.path(gate_dir,
                paste0("gating.png")),
      width = n_plots * 300, height = length(files) * 300)
  layout(matrix(1:(length(files) * n_plots), ncol = n_plots, byrow = TRUE))
  
  for(file in files){
    
    print(paste0("Preprocessing ", file))
    print(Sys.time())
    
    ff <- read.FCS(file.path(cytof_clean_dir, file), transformation = FALSE)
    ff_t <- transform(ff, 
                      transformList(colnames(ff)[grep("Di", colnames(ff))], 
                                    CytoNorm::cytofTransform))
    
    selection <- matrix(TRUE,
                        nrow = nrow(ff),
                        ncol = 3,
                        dimnames = list(NULL,
                                        c("intact", 
                                          "life", 
                                          "cd45")))
    
    tr <- list()
    for(m in c("Ir193Di", "Ir191Di")){
      
      tr[[m]] <- c(deGate(ff_t, m,
                          tinypeak.removal = 0.8, upper = FALSE, use.upper = TRUE,
                          alpha = 0.05, percentile = .95, verbose = F, count.lim = 3), #alpha
                   deGate(ff_t, m,
                          tinypeak.removal = 0.8, upper = TRUE, use.upper = TRUE,
                          alpha = 0.1, percentile = .95, verbose = F, count.lim = 3)) #alpha
    }
    
    for(m in c("Ir193Di", "Ir191Di")){
      selection[ff_t@exprs[,m] < tr[[m]][1], "intact"] <- FALSE
      selection[ff_t@exprs[,m] > tr[[m]][2], "intact"] <- FALSE
    }
    
    percentage <- (sum(selection)/length(selection))*100
    plotDens(ff_t, c("Ir193Di", "Ir191Di"), 
             main = paste0(file," ( ", format(round(percentage, 2), nsmall = 2), "% )"))
    
    abline(h = c(tr[["Ir191Di"]]))
    abline(v = c(tr[["Ir193Di"]]))
    points(ff_t@exprs[!selection[,"intact"], c("Ir193Di", "Ir191Di")], pch = ".")
    
    subset <- selection[, "intact"]
    selection[, "life"] <- subset
    
    v_ch <- grep(viability_ch, colnames(ff), value = T)
    
    tr <- list()
    for(m in c("Ir191Di", v_ch)){
      if (m == v_ch) {
        upper = TRUE
        alpha = 0.1
        tr[[m]] <- deGate(ff_t[subset,], m,
                          tinypeak.removal = 0.8, upper = upper, use.upper = TRUE,
                          alpha = alpha, percentile = .80, verbose = F, count.lim = 3)
        
      } else {
        alpha = 0.05
        tr[[m]] <- c(deGate(ff_t[subset,], m,
                            tinypeak.removal = 0.2, upper = FALSE, use.upper = TRUE,
                            alpha = 0.03, percentile = .95, verbose = F, count.lim = 3), #0.017
                     deGate(ff_t[subset,], m,
                            tinypeak.removal = 0.2, upper = TRUE, use.upper = TRUE,
                            alpha = alpha, percentile = .95, verbose = F, count.lim = 3)) 
        
      }
    }
    
    for(m in c(v_ch, "Ir191Di")){
      if (m == v_ch) {
        selection[ff_t@exprs[,m] > tr[[m]][1], "life"] <- FALSE 
      } else {
        selection[ff_t@exprs[,m] < tr[[m]][1], "life"] <- FALSE
        selection[ff_t@exprs[,m] > tr[[m]][2], "life"] <- FALSE  
      }
    }
    percentage <- (sum(selection[,"life"]/sum(subset)))*100  
    plotDens(ff_t[subset,], c(v_ch, "Ir191Di"), 
             main = paste0(file," ( ", format(round(percentage, 2), nsmall = 2), "% )"),
             xlim = c(0, 8), ylim = c(0, 8))
    
    abline(h = tr[["Ir191Di"]])
    abline(v = tr[[v_ch]])
    
    points(ff_t@exprs[subset & !selection[,"life"], c(v_ch, "Ir191Di")], pch = ".") 
    
    # subset <- selection[, "life"]
    # selection[, "cd45"] <- subset
    # 
    # m <- grep(cd45_ch, colnames(ff), value = T)
    # 
    # tr <- deGate(ff_t[subset,], m,
    #              tinypeak.removal = 0.1, upper = FALSE, use.upper = TRUE,
    #              alpha = 0.05, percentile = .95, verbose = F, count.lim = 3)
    # 
    # if (tr < 1.7) {
    #   tr <- 1.7
    # }
    
    # selection[ff_t@exprs[,m] < tr[1], "cd45"] <- FALSE 
    # percentage <- (sum(selection[,"cd45"]/sum(subset)))*100
    # plotDens(ff_t[subset,], c("Ir191Di", m), 
    #          main = paste0(file," ( ", format(round(percentage, 2), nsmall = 2), "% )"), 
    #          xlim = c(0, 8), 
    #          ylim = c(0, 8))
    # abline(h = tr)
    # points(ff_t@exprs[subset & !selection[,"cd45"], c("Ir191Di", m)], pch = ".")
    
    write.FCS(ff[selection[,"life"], ],
              file.path(out_dir,
                        gsub(".fcs", "_gated.fcs", file)))
    
    saveRDS(selection,
            file.path(out_dir,
                      gsub(".fcs", "_gatingMatrix.RDS", file)))
    
  }
  dev.off()
}

#TODO put option for trsanforming or not 
plot_batch <- function(files_before_norm , files_after_norm, out_dir = getwd(), 
                       clustering_markers = "CD|HLA|IgD|PD|BAFF|TCR", 
                       functional_markers = NULL, 
                       batch_pattern = "RUN[0-9]*"){
  
  # ff_agg  <- file.path(fSOM_dir, paste0("norm_REF_aggregated_", p, ".fcs"))
  # if (file.exists(ff_agg)){
  #   ff_agg <- read.FCS(ff_agg)
  # } else {
  
  files_list <- list("files_before_norm" = files_before_norm, 
                     "files_after_norm" = files_after_norm)
  
  plots <- list()
  for (name in names(files_list)) {
    
    set.seed(1)
    ff_agg <- AggregateFlowFrames(fileNames = files_list[[name]],
                                  cTotal = length(files_list[[name]]) * 1000,
                                  verbose = TRUE,
                                  writeMeta = FALSE,
                                  writeOutput = FALSE,
                                  outputFile = file.path(out_dir, paste0("norm_REF_aggregated_", ".fcs")))
    # }
    
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
    
    set.seed(1)
    samp <- length(files_list[[name]])
    ff_samp <- ff_agg@exprs[sample(nrow(ff_agg@exprs), samp*500), ]
    
    
    #TODO use the ce of nocicka do show mds or umap but by sample ypu have a code in TODO doc in google
    set.seed(123)
    dimred_res <- uwot::umap(X = ff_samp[, names(cl_markers)], 
                             n_neighbors = 15, scale = TRUE)
    
    dimred_df <- data.frame(dim1 = dimred_res[,1], dim2= dimred_res[,2],
                            ff_samp[, names(cl_markers)])
    
    dimred_df$file_id <- ff_samp[,"File2"]
    dimred_df$batch <- NA
    
    for (i in 1:length(files_list[[name]])){
      
      file <- files_list[[name]][i]
      batch <- stringr::str_match(file, batch_pattern)[,1]
      dimred_df[dimred_df[, "file_id"] == i, "batch"] <- batch
    }
    
    p <- ggplot(dimred_df,  aes_string(x = "dim1", y = "dim2", color = "batch")) +
      geom_point(aes(color = batch), size = 0.8, position="jitter") +
      theme_bw() +
      ggtitle(name)+
      # scale_color_manual(values=col )+
      # scale_color_gradientn(norm_markers[m], 
      #                       colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))+
      theme(legend.position = "bottom")
    
    p 
    plots[[name]] <- p
    
    # if (!(is.null(functional_markers))){
    #   
    #   names <- gsub("_CC_gated.fcs|gated.fcs|.fcs", "", basename(files_list[[name]]))
    #   labels <- str_match(basename(files_list[[name]]), batch_pattern)[,1]
    #   
    #   plot_aggregate(names = names, files = files_list[[name]], 
    #                  markers = functional_markers, arcsine_transform = TRUE, 
    #                  labels = labels, flow_frame_agg = ff_agg, 
    #                  output_image = file.path(out_dir, paste0(name, "_aggregate.png")))
    # }
  }
  
  pdf(file.path(norm_dir, "batch_effect"), width = length(plots)*6, height = 6)
  gridExtra::grid.arrange(grobs = plots, ncol = 2)
  dev.off()
  
  
  if (!(is.null(functional_markers))){
    
    # batch_id <- stringr::str_match(basename(files_after_norm), 
    #                                batch_pattern)[,1]
    # file_id <- gsub("_CC_gated.fcs|.fcs|_gated.fcs|Norm_", "", basename(files_after_norm))
    
    markers <- c(clustering_markers, functional_markers)
    plot_marker_quantiles(files_after_norm = files_after_norm, 
                          files_before_norm = files_before_norm,
                          batch_pattern = batch_pattern, 
                          # file_id = file_id, 
                          arcsine_transform = TRUE, 
                          markers_to_plot = markers, out_dir = out_dir)
    
  }
}

plot_aggregate <- function(names = NULL, 
                           files = files, 
                           markers = NULL, 
                           labels = NULL, 
                           output_image = "aggregate.png", 
                           flow_frame_agg = NULL,
                           arcsine_transform = TRUE,
                           max_nm = 6 * 5000, #6 * 100000
                           ymargin = c(-2,10),
                           agg_name = "File2",
                           colors_man = NULL){ #c(-5,10) for manual transformations, #c(-100,300) for FlowJO transformations
  
  if (!is.null(flow_frame_agg)){
    ff <- AggregateFlowFrames(fileNames = files, maxcells = 1000, 
                              cTotal = length(files)*1000, writeOutput = FALSE, 
                              outputFile = FALSE, writeMeta = FALSE)
  }
  
  matches <- paste(markers, collapse="|")
  m <- grep(matches, get_markers(ff, colnames(ff)), 
            ignore.case = TRUE, value = TRUE)
  ch <- get_channels(ff, m)
  
  if(arcsine_transform == TRUE){
    
    ff <- transform(ff, transformList(ch, cytofTransform))
  }
  
  data <- ff@exprs
  file_values <- data[,agg_name]
  #file_values_scattered <- data[,gsub("File","File_scattered",agg_name)]
  
  subset <- sample(seq_len(nrow(data)), min(max_nm, nrow(data)))
  if(is.null(markers)) {
    data <- data[subset,]
  } else {
    data <- data[subset, ch]
  }
  
  file_values <- file_values[subset]
  file_values_scattered <- file_values + stats::rnorm(length(file_values), 
                                                      0, 0.1)
  channels <- colnames(data)
  nrows_in_plot <- floor(sqrt(length(channels)))
  ncols_in_plot <- ceiling(length(channels)/nrows_in_plot)
  png(output_image, width = 800*ncols_in_plot, height = 1000*nrows_in_plot)
  par(mar = c(30.1, 12.1, 2.1, 2.1))
  layout(matrix(seq_len(nrows_in_plot * ncols_in_plot), nrow = nrows_in_plot, 
                byrow = TRUE))
  if (is.null(labels)) {
    colors <- "#00000055"
  } else {
    labels <- factor(labels)
    color_palette <- colorRampPalette(RColorBrewer::brewer.pal(9,
                                                               "Set1"))
    colors <- paste0(color_palette(length(levels(labels))),
                     "55")
    if (is.null(colors_man)) {
      colors <- colors[labels][file_values]
    } else {
      colors <- colors_man[labels][file_values]
    }
  }
  xlabels <- names#labels 
  for (channel in channels) {
    print(FlowSOM::get_markers(ff, channel))
    plot(0,
         type = "n", xaxt = "n",
         xlab = "", ylab = FlowSOM::get_markers(ff, channel),
         cex.lab = 5,
         ylim = ymargin, 
         xlim = c(0, length(files) + 1))
    abline(v = seq_len(length(files)), col = "lightgrey")
    axis(side = 1, at = seq_len(length(files)), labels = xlabels, 
         las = 2, cex.axis = 1.5)
    points(file_values_scattered, data[, channel], pch = ".", 
           col = colors, cex = 3)
  }
  dev.off()
}



  