#' Finds all the mass channels
#' @param flow_frame untranfosrmed flow frame
#' @param channels pattern for non-mass channels
#' @return logical vector with TRUE values for mass channels

find_mass_ch <- function(flow_frame, channels = "Time|Event_length|Center|Offset|Width|Residual"){
  non_mass_ch <- grep(c(channels), 
       colnames(ff), 
       value = TRUE,
       invert = TRUE)
  return(non_mass_ch)
}


#' Cleans the flow rate using functions from flowAI package
#' @param flow_frame untranfosrmed flow frame
#' @param to_plot logical if to plot cleaning results
#' @param out_dir pathway to where the plots should be saved, 
#' only if to_plot = TRUE
#' @return cleand, untransformed flow frame

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
                  gsub(".fcs", "_flowAI.png", basename(ff@description$FILENAME))),
        width = 800,
        height = 600)
    p <- flowAI:::flow_rate_plot(FlowRateQC)
    print(p)
    dev.off()
  }
  
  flow_frame_cl <- FlowRateQC$FRnewFCS
  flow_frame_cl@exprs[,"Time"] <- flow_frame_cl@exprs[,"Time"]*100

  return(flow_frame_cl)
}

#' Cleans the flow rate using functions from flowAI package
#' @param flow_frame untranfosrmed flow frame
#' @param channels_to_clean vector of the channels that needs to be cleaned
#' @param to_plot a characer variable that indicates if plots should be genarated.
#' The default is "All", which generates plot for all channels. Other options are
#' "Flagged Only", plots the channels that were spotted with flowcut as incorrect
#' and "None", does not plots anything.
#' @param out_dir pathway to where the plots should be saved, 
#' only if to_plot = TRUE
#' @return cleand, untransformed flow frame

clean_signal <- function(flow_frame, channels_to_clean = NULL, to_plot = "All", 
                         out_dir = getwd()){
  
  channels_to_transform <- find_mass_ch(flow_frame)
  ff_t <- transform(flow_frame, transformList(channels_to_transform, 
                                              CytoNorm::cytofTransform))
  
  if (!is.null(channels_to_clean)){
    
    ch_to_clean <- colnames(flow_frame)[colnames(flow_frame) %in% 
                                          channels_to_clean]
    if(any(ch_to_clean != "Time" | ch_to_clean != "time")){
      
      ind <- grep("TIME", toupper(colnames(flow_frame)))
      ch_to_clean <- c(ch_to_clean, colnames(flow_frame)[ind]) 
      ch_ids <- colnames(flow_frame) %in% ch_to_clean
    }
    
  } else {
    ind <- grep("TIME", toupper(colnames(flow_frame)))
    ch_to_clean <- c(colnames(flow_frame)[ind], channels_to_transform)
    ch_ids <- colnames(flow_frame) %in% ch_to_clean
  }
  
  out_dir <- file.path(out_dir, "SignalCleaning")
  if(!dir.exists(out_dir)){
    dir.create(out_dir)
  }
  
  cleaned_data <- flowCut(ff_t[,ch_ids], 
                          Segment = 1000, 
                          MaxPercCut = 0.5,
                          FileID = gsub("_beadNorm", "_flowCutCleaned", file),
                          Plot = to_plot,
                          Directory = out_dir,
                          UseOnlyWorstChannels = TRUE,
                          AllowFlaggedRerun = TRUE)
  
  ff_t_clean <- cleaned_data$frame
  ff_clean <- transform(ff_t_clean, transformList(channels_to_transform, 
                                                  cytofTransform.reverse))
  
  return(ff_clean)
}


#' performs bead normalization usinf function from CATALYST package
#' @param flow_frame untranfosrmed flow frame
#' @param channels_to_keep vector of the channels that needs to be kept for 
#' further analysis
#' @param to_plot logical if to plot bead gate and bead normalization for each 
#' file
#' @param out_dir pathway to where the plots should be saved, 
#' only if to_plot = TRUE, the default is working directory
#' @return bead normalized flow frame without the beads

bead_normalize <- function(flow_frame, keep_all_markers = TRUE, 
                           marker_pattern = "PD|CD|HLA|IgD|TCR|BAFF|Ir|Pt195", 
                           beads = "dvs",
                           to_plot = TRUE, out_dir = getwd()){
  
  if (keep_all_markers == FALSE){
    
    if (is.null(marker_pattern)){
      stop("pattern for markers needs to be specify")
    }
 
    markers_to_keep <- grep(marker_pattern, 
                            get_markers(flow_frame, colnames(flow_frame)), 
                             ignore.case = TRUE, value = FALSE)
    non_mass_ch <- grep("Time|length|Ce140|Center|
                        Offset|Width|Residual", colnames(ff),
                        ignore.case = TRUE, value = FALSE)
    
    channels_to_keep <- c(markers_to_keep, non_mass_ch)
    channels_to_keep <- colnames(flow_frame)[sort(unique(channels_to_keep))]
    
    flow_frame <- flow_frame[, channels_to_keep]
  }
  
  dat <- prepData(flow_frame) 
  
  # normalize the data and remove beads
  dat_norm <- normCytof(x = dat,
                        beads = beads,
                        remove_beads = TRUE,
                        norm_to = ff_ref,
                        k = 80,
                        plot = TRUE)
  
  # convert back to .fcs files and save 
  ff <- sce2fcs(dat_norm$data)
  # write.FCS(ff, file.path(bead_norm_dir, sub_dir, gsub(".FCS","_beadNorm.fcs", file)))
  
  if (to_plot == TRUE){
    # plot and save diagnostic plots 
    dat_norm$scatter
    ggsave(filename = file.path(out_dir, gsub(".FCS","_beadGate.png", file)))
    dat_norm$lines
    ggsave(filename = file.path(out_dir, gsub(".FCS","_beadLines.png", file)))
  }
  
  return(ff)
}














  