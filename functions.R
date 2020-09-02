#' @description Finds all the mass channels
#' @param flow_frame Untranfosrmed flow frame
#' @param channels Pattern for non-mass channels
#' @param ... Additional arguments to pass to grep
#' @return Logical vector with TRUE values for mass channels

find_mass_ch <- function(flow_frame, 
                         channels = "Time|Event_length|Center|Offset|Width|Residual",
                         ...){
  non_mass_ch <- grep(c(channels), 
       colnames(ff), 
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

#'@description Creates the reference flowframe which bead´s values will be used as reference 
#'for other files during the normalization.
#'@param fcs_files Path to selected fcs file which will be used as reference 
#'file upon normalization.
#'@param beads The same as in CATALYAST package, character variable:
#'"dvs" (for bead masses 140, 151, 153 ,165, 175) or 
#'"beta" (for bead masses 139, 141, 159, 169, 175) or a numeric vector of masses.

create_ref <- function(fcs_file, beads = NULL){
  
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
#' @param channels_to_keep Vector of the channels that needs to be kept for 
#' further analysis.
#' @param to_plot Logical if to plot bead gate and bead normalization for each 
#' file.
#' @param out_dir Pathway to where the plots should be saved, 
#' only if to_plot = TRUE, the default is working directory.
#' @return Bead normalized flow frame without the beads.

bead_normalize <- function(flow_frame, keep_all_markers = TRUE, 
                           markers_to_keep = "PD|CD|HLA|IgD|TCR|BAFF|Ir|Pt195", 
                           beads = "dvs", norm_to_ref = NULL, 
                           remove_beads = TRUE, to_plot = TRUE, 
                           out_dir = getwd(), ...){
  
  if (keep_all_markers == FALSE){
    
    if (is.null(markers_to_keep)){
      stop("pattern for markers needs to be specify")
    }
    
    m_to_keep <- grep(markers_to_keep, 
                            get_markers(flow_frame, colnames(flow_frame)), 
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














  