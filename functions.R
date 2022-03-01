
check <- function(x) tryCatch(if(class(x) == 'logical') 1 else 1, 
                              error=function(e) 0)


#' @description Finds all the mass channels
#' 
#' @param flow_frame Untransformed flow frame
#' @param channels Pattern for non-mass channels
#' @param ... Additional arguments to pass to grep
#' 
#' @return Logical vector with TRUE values for mass channels

find_mass_ch <- function(flow_frame, 
                         channels = "Time|Event_length|Center|Offset|Width|Residual|SSC|FSC|File_scattered",
                         ...){
  non_mass_ch <- grep(c(channels), 
       colnames(flow_frame), 
       invert = TRUE, ...)
  return(non_mass_ch)
}


#' flow_rate_bin_addapted
#' 
#' @param x flow frame
#' @param second_fraction the fraction of the seconds used in the data 
#' @param timeCh Time channel 
#' @param timestep 
#' 
#' @return timeFlowData, the cell assignment to the bins

flow_rate_bin_addapted <- function (x, second_fraction = 0.1, timeCh = timeCh, timestep = timestep) 
{
  xx <- exprs(x)[, timeCh]
  idx <- c(1:nrow(x))
  endsec <- ceiling(timestep * max(xx))
  lenx <- length(xx)
  secbegin <- as.numeric(gsub("(.*)(\\.)(.{0}).*", "\\1\\2\\3", xx[1]))
  tbins <- seq(secbegin, endsec/timestep, by = as.numeric(second_fraction)/timestep)
  if (tail(tbins, n=1) < endsec/timestep){
    tbins <- c(tbins, tail(tbins, n=1) + 10)
  }
  if (secbegin == 0){
    secbegin2 <- 0
  } else {
    secbegin2 <- as.numeric(gsub("(.*)(\\.)(.{1}).*", "\\1\\2\\3", xx[1]/100))
  }
  
  secbin <- seq(secbegin2, endsec, by = as.numeric(second_fraction))
  minbin <- round(secbin/60, 3)
  nrBins <- length(tbins) - 1
  tbCounts <- c(0, hist(xx, tbins, plot = FALSE)$counts)
  expEv <- lenx/(nrBins)
  binID <- do.call(c, mapply(rep, x = 1:length(tbCounts), 
                             times = tbCounts, SIMPLIFY = FALSE))
  if (length(idx) != length(binID)) 
    stop("length of cell ID not equal length of bin ID")
  timeFlowData <- list(frequencies = cbind(tbins, minbin, 
                                           secbin, tbCounts), 
                       cellBinID = data.frame(cellID = idx, 
                                              binID = binID), info = data.frame(second_fraction = second_fraction, 
                                                                                expFrequency = expEv, bins = nrBins))
  return(timeFlowData)
}


#' clean_flow_rate
#' 
#' @description Cleans the flow rate using functions from flowAI package.
#' 
#' @param flow_frame Untransformed flow frame
#' @param to_plot Logical if to plot cleaning results, default is set to TRUE.
#' @param out_dir Character, pathway to where the plots should be saved, 
#' only if argument to_plot = TRUE, default is set to working directory.
#' @param alpha numeric, as in flow_auto_qc {flowAI}. The statistical 
#' significance level used to accept anomalies. The default value is 0.01.
#' @param data_type Character, if MC (mass cytometry) of FC (flow cytometry) data 
#' are analyzed
#' 
#' @return Clean, untransformed flow frame and save the plot _beadNorm_flowAI.png
#'  in out_dir 
clean_flow_rate <- function(flow_frame, to_plot = TRUE, 
                            out_dir = getwd(), alpha = 0.01, data_type = "MC") {
  
  if (data_type == "MC"){
    time_division <- 100
    timestep <- 0.01
  } else if (data_type == "FC") {
    time_division <- 1
    word <- which(grepl("TIMESTEP", names(keyword(flowSet(ff)[[1]])), 
                        ignore.case = TRUE))
    timestep <- as.numeric(keyword(flowSet(ff)[[1]])[[word[1]]])
  } else{
    stop("type of data MC or FC needs to be specify")
  }
  
  
  flow_frame@exprs[, "Time"] <- flow_frame@exprs[, "Time"]/time_division
  
  
  FlowRateData <- flow_rate_bin_addapted(flow_frame, 
                                         timeCh = "Time", 
                                         timestep = timestep)
  
  FlowRateQC <- flowAI:::flow_rate_check(x = flow_frame, 
                                         FlowRateData = FlowRateData, 
                                         alpha = alpha, 
                                         use_decomp = TRUE) 
  
  if (to_plot == TRUE){
    
    out_dir <- file.path(out_dir, "FlowRateCleaning")
    if(!dir.exists(out_dir)){
      dir.create(out_dir)
    }
    
    png(file.path(out_dir,
                  gsub(".fcs", "_flowAI.png", 
                       basename(flow_frame@description$FILENAME), ignore.case = TRUE)),
        width = 800,
        height = 600)
    if(data_type == "MC"){
      FlowRateQC$res_fr_QC[,1] <- timestep
    }
    p <- plot_flowrate(FlowRateQC, data_type = data_type)
    print(p)
    dev.off()
  }
  
  flow_frame_cl <- flow_frame[FlowRateQC$goodCellIDs,]
  flow_frame_cl@exprs[,"Time"] <- flow_frame_cl@exprs[,"Time"]*time_division
  
  return(flow_frame_cl)
}


#' plot_flowrate
#' 
#' @description plots flow rate for .fcs files
#'
#' @param data_type 
#' @param FlowRateQC list obtained using flowAI:::flow_rate_check function
#'
#' @return xgraph plot 
plot_flowrate <- function (FlowRateQC, data_type = "MC") 
{
  if (data_type == "MC"){
    lab <- "Time (10 * Seconds)"
  } else {
    lab <- "Time (Seconds)"
  }
  second_fraction <- FlowRateQC$res_fr_QC$second_fraction
  num_obs = FlowRateQC$res_fr_QC$num_obs
  frequencies = as.data.frame(FlowRateQC$frequencies)
  anoms = as.data.frame(FlowRateQC$anoms)
  anoms_points = as.data.frame(cbind(sec_anom = frequencies$secbin[anoms$index], 
                                     count_anom = anoms$anoms))
  xgraph <- ggplot2::ggplot(frequencies, aes_string(x = "secbin", y = "tbCounts")) + 
    theme_bw() + theme(panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), text = element_text(size = 30)) + 
    geom_line(colour = "darkblue")
  xgraph <- xgraph + ggplot2::labs(x = lab, y = paste0("Number of events per 1/", 
                                                       1/second_fraction, " of a second"))
  if (!is.null(anoms_points)) {
    xgraph <- xgraph + ggplot2::geom_point(data = anoms_points, aes_string(x = "sec_anom", 
                                                                           y = "count_anom"), 
                                           color = "green4", size = 5, 
                                           shape = 1, stroke = 3)
  }
  return(xgraph)
}


#' clean_signal
#' 
#' @description Cleans the signal using flowCut package
#' 
#' @param flow_frame Flow frame, if unstransformed arcsine_transform should be 
#' kept as default, TRUE.
#' @param channels_to_clean Character vector of the channels that needs to be cleaned
#' @param to_plot Characer variable that indicates if plots should be genarated.
#' The default is "All", which generates plot for all channels. Other options are
#' "Flagged Only", plots the channels that were spotted with flowcut as incorrect
#' and "None", does not plots anything.
#' @param Segment As in flowCut, an integer value that specifies the 
#' number of events in each segment to be analyzed. 
#' Each segment is defaulted to 1000 events.
#' @param out_dir Character, pathway to where the plots should be saved, 
#' only if argument to_plot = TRUE, default is set to working directory.
#' @param arcsine_transform Logical, if the data should be transformed with 
#' arcsine transformation and cofactor 5, default is set to TRUE
#' @param non_used_bead_ch Character vector, bead channels that does not contain 
#' any marker information, thus do not need to be cleaned and used for further analysis
#' @param ... Additional arguments to pass to flowcut.
#' @param MaxPercCut As in flowCut, numeric between 0-1 the maximum percentage of 
#' event tha will be removed form the data. 
#' @param UseOnlyWorstChannels as in flowCut, logical, automated detection of the 
#' worst channel that will be used for cleaninig 
#' @param AllowFlaggedRerun as in flowCut, logical, specify if flowCut will run
# second time in case the file was flagged
#' @param AlwaysClean as in flowCut, logicle. The file will be cleaned even if it has a
#' relatively stable signal. The segments that are 7 SD away from the mean of all
#' segments are removed 
#' 
#' @return Cleaned, untransformed flow frame if arcsine_transform argument
#' set to TRUE, otherwise transformed flow frame is returned. Save plots with prefix
#' "flowCutCleaned.png" to out_dir if parameter to_plot set to "All" or "Flagged Only".
clean_signal <- function(flow_frame, 
                         channels_to_clean = NULL, 
                         to_plot = "All", 
                         Segment = 1000,
                         out_dir = getwd(), 
                         arcsine_transform = TRUE, 
                         non_used_bead_ch = NULL, 
                         MaxPercCut = 0.5,
                         UseOnlyWorstChannels = TRUE,
                         AllowFlaggedRerun = TRUE,
                         AlwaysClean = TRUE,
                         data_type = "MC",
                         ...){
  
  channels_to_transform <- find_mass_ch(flow_frame, value = FALSE)
  
  if (arcsine_transform == TRUE){
    
    if(data_type == "MC"){
      ff_t <- flowCore::transform(flow_frame, 
                                  transformList(colnames(flow_frame)[channels_to_transform], 
                                                CytoNorm::cytofTransform))
    } else if (data_type == "FC"){
      ff_t <- flowCore::transform(flow_frame, 
                                  transformList(colnames(flow_frame)[channels_to_transform], 
                                                arcsinhTransform(a = 0, b = 1/150, c = 0))) 
      
    } else {
      stop("specify data type MC or FC")
    }
    
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
    
    if (!is.null(non_used_bead_ch)) {
      non_bead_ch <- "140"
    } else {
      non_bead_ch <- paste(non_used_bead_ch, collapse="|")
    }

    ind_Time <- grep("TIME", colnames(flow_frame), value = T, ignore.case = T)
    ch_to_clean <- c(ind_Time, find_mass_ch(flow_frame, value = TRUE))
    ind_nonbeads <- grep(non_bead_ch, colnames(flow_frame), value = TRUE)
    channels <- ch_to_clean[!(ch_to_clean %in% ind_nonbeads)] 
    channels <- grep(paste(channels, collapse = "|"), colnames(flow_frame))
  }
  
  out_dir <- file.path(out_dir, "SignalCleaning")
  if(!dir.exists(out_dir)){
    dir.create(out_dir)
  }
  
  cleaned_data <- flowCut::flowCut(f = ff_t, 
                                   Segment = Segment, 
                                   MaxPercCut = MaxPercCut, 
                                   Channels = channels,
                                   FileID = gsub("_beadNorm", "_flowCutCleaned", 
                                                 basename(flow_frame@description$FILENAME)),
                                   Plot = to_plot,
                                   Directory = out_dir,
                                   UseOnlyWorstChannels = UseOnlyWorstChannels,
                                   AllowFlaggedRerun = AllowFlaggedRerun,
                                   AlwaysClean = AlwaysClean)
  
  ff_t_clean <- cleaned_data$frame
  
  if (arcsine_transform == TRUE){
    ff_clean <- flowCore::transform(ff_t_clean, 
                                    transformList(colnames(flow_frame)[channels_to_transform], 
                                                  cytofTransform.reverse))
  } else {
    ff_clean <- ff_t_clean
  }
  
  return(ff_clean)
}

#' baseline_file
#' 
#' @description Creates the reference flowframe for which bead´s mean values will 
#' be computed and use during the normalization.
#'
#' @param fcs_files Path to fcs files to be normalized
#' @param beads The same as in normCytof from CATALYAST package, character variable:
#'"dvs" (for bead masses 140, 151, 153 ,165, 175) or 
#'"beta" (for bead masses 139, 141, 159, 169, 175)
#' @param to_plot Logicle variable that indicates if plots should be genarated, 
#' default set to FALSE
#' @param out_dir Character, pathway to where the plots should be saved, 
#' only if argument to_plot = TRUE, default is set to working directory.
#' @param k the same as in normCytof from CATALYST package, integer width of the 
#' median window used for bead smoothing (affects visualizations only!).
#' @param ... Additional arguments to pass to normCytof
#' @param ncells number of cells to be aggregated per each file, 
#' @param seed numeric, set to obtain reproducible results, default 654
#' default is set to 25000, so around 250 beads can be aggregated per each file 
#' 
#' @return returns reference, aggregated flow frame 

baseline_file <- function(fcs_files, beads = "dvs", to_plot = FALSE, 
                       out_dir = getw(), k = 80, ncells = 25000, seed = 2, ...){

  set.seed(seed)
  ff <- FlowSOM::AggregateFlowFrames(fileNames = fcs_files, 
                            cTotal = length(fcs_files)*ncells)
  
  dat <- CATALYST::prepData(ff) 

  dat_norm <- CATALYST::normCytof(x = dat,
                        beads = beads,
                        remove_beads = TRUE,
                        norm_to = NULL,
                        k = k,
                        plot = to_plot, 
                        verbose = FALSE, 
                        transform = FALSE, 
                        ...)
  ff_ref <- CATALYST::sce2fcs(dat_norm$beads)
  rm(ff)
  
  if (to_plot == TRUE){
    
    # plot and save diagnostic plots 
    plot_dir <- file.path(out_dir, "Plots_BeadNormalization")
    if(!dir.exists(plot_dir))(dir.create(plot_dir))
    
    p <- dat_norm$scatter
    ggplot2::ggsave(filename = file.path(plot_dir,"RefBeadGate.png"), 
           plot = p, limitsize = FALSE)
    
    p <- dat_norm$lines
    ggplot2::ggsave(filename = file.path(plot_dir,"RefBeadLines.png"), 
           plot = p, limitsize = FALSE)
  }
  return(ff_ref)
}

#' bead_normalize
#' 
#' @description Performs bead normalization using function from CATALYST package.
#' normalized fcs files and plots are stored in out_dir directory
#'
#' @param flow_frame Untranfosrmed flow frame.
#' @param markers_to_keep Character vector, marker names to be kept after 
#' the normalization, can be full marker name e.g. "CD45" or "CD". Additionally, 
#' non_mass channels like Time, Event_lenght, Gaussian parameter and palladium 
#' barcoding channels are kept
#' If NULL (default) all markers will be normalized and kept in flowframe, 
#' selection of the markers will reduce file volume and speedup the analysis. 
#' @param beads character, as in normCytof, "dvs" (for bead masses 140, 151, 153 ,165, 175) 
#' or "beta" (for bead masses 139, 141, 159, 169, 175) or a numeric vector of masses. 
#' Default is set to "dvs"
#' @param norm_to_ref flow frame, created by baseline_file function to which input data 
#' will be normalized, default is set to NULL
#' @param to_plot Logical if to plot bead gate and bead normalization lines for each 
#' file.
#' @param out_dir Character, pathway to where the plots should be saved, 
#' only if argument to_plot = TRUE, default is set to working directory.
#' @param k the same as in normCytof from CATALYST package, integer width of the 
#' median window used for bead smoothing (affects visualizations only!).
#' @param remove_beads 
#' @param ... Additional arguments to pass to normCytof
#'
#' @return Bead normalized flow frame without the beads. Save plots to out_dir
#'  if argument to_plot set to TRUE

bead_normalize <- function(flow_frame,  
                           markers_to_keep = NULL, 
                           beads = "dvs", 
                           norm_to_ref = NULL, 
                           remove_beads = TRUE, 
                           to_plot = TRUE, 
                           out_dir = getwd(), 
                           k = 80, 
                           ...){
  
  if (!is.null(markers_to_keep)){

    matches <- paste(markers_to_keep, collapse="|")
    
    m_to_keep <- grep(matches, FlowSOM::GetMarkers(flow_frame, colnames(flow_frame)), 
                             ignore.case = TRUE, value = FALSE)
    
    non_mass_ch <- grep("Time|length|Ce140|151|153|165|175|Center|Offset|Width|
                        |Residual|Pd", 
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
  
  dat <- CATALYST::prepData(flow_frame) 
  
  # normalize the data and remove beads
  dat_norm <- CATALYST::normCytof(x = dat,
                                  beads = beads,
                                  remove_beads = remove_beads,
                                  norm_to = norm_to_ref,
                                  k = k,
                                  plot = TRUE, 
                                  transform = FALSE, 
                                  ...)
 
  # convert back to .fcs files and save 
  f <- CATALYST::sce2fcs(dat_norm$data)
  
  filename <- basename(flow_frame@description$FILENAME)
  
  if (to_plot == TRUE){
    
    # plot and save diagnostic plots 
    plot_dir <- file.path(out_dir, "Plots_BeadNormalization")
    if(!dir.exists(plot_dir))(dir.create(plot_dir))
    
    p <- dat_norm$scatter
    ggplot2::ggsave(filename = file.path(plot_dir, gsub(".FCS|.fcs","_beadGate.png", filename)), 
           plot = p, limitsize = FALSE)
   
    p <- dat_norm$lines
    ggplot2::ggsave(filename = file.path(plot_dir, gsub(".FCS|.fcs","_beadLines.png", filename)), 
           plot = p, limitsize = FALSE)
   
  }
  
  f@description$FILENAME <- basename(flow_frame@description$FILENAME)
  f@description$FIL <- basename(flow_frame@description$FILENAME)
  
  return(f)
}


#' gate_out_beads 
#'
#' @description removes beads from the files that contain them e.g files before
#' bead normalization
#'
#' @param bead_channel character, the mass for bead channel that is exclusively used for 
#' beads identification, no marker is present at this channel, default 140.
#' @param flow_frame Flow frame, if unstransformed arcsine_transform should be 
#' kept as default, TRUE.
#'
#' @return flow frame without beads removed 
gate_out_beads <- function(bead_channel, 
                           flow_frame){
  ch <- grep(pattern = bead_channel, x = colnames(flow_frame), value = TRUE)
  ids <- flow_frame[,ch] > 0
  # calculate threshold 
  th <- deGate(obj = flow_frame[ids,], channel = ch)
  #remove beads 
  cells_to_remove <- flow_frame@exprs[, ch] < th
  flow_frame <- flow_frame[cells_to_remove,]
  return(flow_frame)
}

#' plot_marker_quantiles
#' 
#' @description Calculates quantiles for selected markers and plots them as 
#' diagnostic plot 
#'
#' @param files_before_norm Character, full path to the unnormalized fcs_files.
#' @param files_after_norm Character, full path to the normalized fcs_files.
#' @param batch_pattern Character, batch pattern to be match in the fcs file name
#' @param uncommon_prefix Character vetor or string, uncommon prefix in the basenames
#' of the fcs files. The file names needas to be exactly matched so, uncommon prefix 
#' needs to be removed. If NULL prefix like 
#' "Norm|_CC_gated.fcs|_gated.fcs|_beadNorm.fcs|.FCS|.fcs" will be removed. 
#' Default is set to NULL.
#' @param bead_channel character, the mass for bead channel that is exclusively used for 
#' beads identification, no marker is present at this channel, default 140.
#' @param remove_beads logical, if beads needs to be removed. This needs to be 
#' set to TRUE if files contains the beads e.g before beads normalization, 
#' default is set to TRUE
#' @param arcsine_transform Logical, if the data should be transformed with 
#' arcsine transformation and cofactor 5.
#' @param markers_to_plot character vector, marker names to be plotted, can be 
#' full marker name e.g. "CD45" or "CD" if all CD-markers needs to be plotted 
#' @param manual_colors character, vector of the colors to be used, 
#' the number of colors needs to be equal to the length of batch_pattern
#' @param out_dir Character, pathway to where the plots should be saved, 
#' default is set to working directory.
#' 
#' @return save the plots to out_dir with the name 
#' "Marker distribution across aliquots and batches.pdf"

plot_marker_quantiles <- function(files_before_norm,
                                  files_after_norm, 
                                  batch_pattern,
                                  remove_beads = FALSE,
                                  bead_channel = "140", 
                                  uncommon_prefix = NULL, 
                                  arcsine_transform = TRUE, 
                                  markers_to_plot = NULL, 
                                  manual_colors = NULL, 
                                  out_dir = getwd()){
  
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
                           FlowSOM::GetMarkers(ff_tmp, find_mass_ch(ff_tmp, 
                                                            value = TRUE)), 
                           value = TRUE, ignore.case = F)
    } else {
      norm_markers <- find_mass_ch(ff_tmp, value = TRUE)
      norm_markers <- FlowSOM::GetMarkers(ff_tmp, norm_markers)
    }
  }
  
  quantile_values <-  c(0.01, 0.25, 0.5, 0.75, 0.99)
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
      ff <- flowCore::transform(ff, transformList(grep("Di", colnames(ff), value = TRUE),
                                        arcsinhTransform(a = 0, b = 1/5, c = 0)))
    }
    
    norm <- quantiles$Normalization[(which(quantiles$File == file)[1])]
    
    if(norm == "NO" & remove_beads){
      ff <- gate_out_beads(bead_channel = bead_channel, flow_frame = ff)
    }
    
    for (marker in names(norm_markers)) {
      quantiles_res <-stats::quantile(exprs(ff)[, marker],
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
  
  if(is.null(uncommon_prefix)){
    quantiles$Sample <- gsub(pattern = "Norm_", replacement = "", ignore.case = TRUE,
                             x = gsub(pattern = "_CC_gated.fcs|_gated.fcs|_beadNorm.fcs|.FCS|.fcs",
                                  replacement = "", ignore.case = TRUE,
                                  x = basename(as.character(quantiles$File))))
  } else {
    uncommon_prefix <- paste(uncommon_prefix, collapse = ("|"))
    quantiles$Sample <- gsub(pattern = "Norm_", replacement = "", ignore.case = TRUE,
                             x = gsub(pattern = uncommon_prefix, replacement = "", 
                                      ignore.case = TRUE,
                                      x =  basename(as.character(quantiles$File))))
  }
  
  ncols <- length(unique(quantiles$Batch))
  p <- quantiles %>% dplyr::filter(Normalization == "YES") %>%
    ggplot2::ggplot(aes(x = Sample,
               y = Value,
               color = Batch)) +
    geom_point(data = quantiles %>% dplyr::filter(Normalization == "NO"), 
               aes(alpha = ifelse(Quantile == "0.5", 2, 0)), color = "grey31") +
    geom_line(data = quantiles %>% dplyr::filter(Normalization == "NO"), 
              aes(alpha = ifelse(Quantile != "0.5", 2, 0),
                  group = interaction(Batch, Quantile),
                  size = ifelse(Quantile %in% c("0.05", "0.95"), 1,
                                ifelse(Quantile == "0.5", 0, 2))),
              alpha = 0.5, color = "grey31") +
    ylab(label = levels(quantiles$Marker))+
    geom_point(aes(alpha = ifelse(Quantile == "0.5", 2, 0))) +
    geom_line(aes(alpha = ifelse(Quantile != "0.5", 2, 0),
                  group = interaction(Batch, Quantile),
                  size = ifelse(Quantile %in% c("0.05", "0.95"), 1,
                                ifelse(Quantile == "0.5", 0, 2))),
              alpha = 0.5) +
    scale_size_identity() +
    scale_alpha_identity() +
    facet_wrap(~ Marker + Batch, ncol = ncols, scales = "free_x") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "bottom")
  
  if (!is.null(manual_colors)){
    p <- p + ggplot2::scale_colour_manual(values = c(manual_colors))
  }
  
  ggplot2::ggsave(filename = "Marker_distribution_across_aliquots_and_batches.pdf", 
                 plot = p, 
                 path = file.path(out_dir), 
                 width = length(fcs_files)*0.25, height = length(norm_markers)*4, limitsize = F)
}

#' fsom_aof
#' 
#' @description Builds FlowSOM tree for the files scoring 
#'
#' @param fcs_files Character, full path to fcs_files.
#' @param phenotyping_markers Character vector, marker names to be used for clustering,
#' can be full marker name e.g. "CD45" or "CD" if all CD-markers needs to be plotted
#' @param nCells Numeric, the total number of cells, to use for FlowSOM clustering. 
#' This number is determined by total number of fcs files, as a defult 1000 cells 
#' is used per file
#' @param xdim Numeric, parameter to pass to FlowSOM, width of the SOM grid
#' @param ydim Numeric, parameter to pass to FlowSOM, geight of the SOM grid
#' @param nClus Numeric, exact number of clusters for metaclustering 
#' @param out_dir Character, pathway to where the FlowSOM clustering plot should
#'  be saved, default is set to working directory.
#' @param batch Character, Character, aqusition batch for each fcs file..
#' Pass to FlowSOM plot name, defult is set to NULL
#' @param arcsine_transform Logical, if the data should be transformed with 
#' arcsine transformation and cofactor 5. Default set to TRUE.
#' @param seed numeric, set to obtain reproducible results, default 654
#' 
#' @return fsom object

fsom_aof <- function(fcs_files, 
                     phenotyping_markers,
                     nCells = length(fcs_files)*10000,
                     xdim = 10,
                     ydim = 10,
                     nClus = 10,
                     out_dir, 
                     batch = NULL,
                     arcsine_transform = TRUE, 
                     seed = 1){
  
  
  if(check(phenotyping_channels) == 0){
    o <- capture.output(ff_tmp <- flowCore::read.FCS(file.path(files[1])))
    markers <- FlowSOM::GetMarkers(ff_tmp, colnames(ff_tmp))
    phenotyping_channels <- grep(paste(phenotyping_markers, 
                                       collapse = ("|")), markers, value = TRUE)
  }
  
  if(arcsine_transform  == TRUE){
    trans <- transformList(names(phenotyping_channels), CytoNorm::cytofTransform)
  } else {
    trans <- NULL
  }
  
  fsom <- CytoNorm::prepareFlowSOM(file = fcs_files,
                                   colsToUse = names(phenotyping_channels),
                                   seed = seed, 
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
  
  if(max(as.numeric(fsom$metaclustering)) > length(myCol)){
    backgroundColors <- NULL
  } else {
    backgroundColors <- myCol
  }
  
  if(!is.null(batch)){
    filename <- paste0(batch, "_FlowSOM_clustering.pdf")
  } else {
    filename <- "FlowSOM_clustering.pdf"
  }
  
  # pdf(file.path(out_dir, filename), width = 14, height = 10)
  fsomPlot <- FlowSOM::PlotStars(fsom = fsom, 
                                 title = "FlowSOM clustering",
                                 backgroundValues = fsom$metaclustering, 
                                 maxNodeSize = 3,
                                 backgroundColors = backgroundColors)
  fsomTsne <- FlowSOM::PlotDimRed(fsom = fsom, plotFile = NULL, seed = seed, cTotal = 20000,  
                                  title = "tSNE visualization of FlowSOM metaclusters")
  
  figure <- ggarrange(fsomPlot, fsomTsne,
                      # labels = c("FlowSOM clustering", "tsne"),
                      ncol = 2, nrow = 1)
  
  ggplot2::ggsave(filename = filename, plot = figure, device = "pdf", path = out_dir,
         width =24, height = 10)
  # dev.off()
  
  return(fsom)
}

#' scaled_aof_score
#' 
#' @description Calculates scaled AOF and sample quality AOF scores 
#'
#' @param aof_scores Matrix, array, Aof scores obtained using function 
#' greedyCytometryAof from cytutils package 
#' @param out_dir Character, pathway to where the FlowSOM clustering plot should
#' @param aof_channels Character, channels or markers for which aof was calculated
#' @param batch Character, aqusition batch for each fcs file.
#' Pass to FlowSOM plot name, defult is set to NULL
#' be saved, default is set to working directory.
#' 
#' @return returns data frame with the scaled AOF scores

scaled_aof_score <- function(aof_scores, out_dir, aof_channels, batch = NULL){
  phenotyping_channels <- aof_channels
  aof_scores_scaled <- scale(aof_scores)
  aof_scores_scaled <- pmax(aof_scores_scaled, 0)^2
  sample_scores <- apply(aof_scores_scaled, 1, sum, na.rm = TRUE)
  
  df <- as.data.frame(sample_scores)
  
  list_scores <- list("aof_scores_per_marker" = aof_scores,
                      "scaled_AOF" = aof_scores_scaled)
  
  for (name in names(list_scores)) {
    
    if(!is.null(batch)){
      filename <- file.path(out_dir, paste0(batch, "_", name, ".pdf"))
      main <- paste0(batch, "_", name) 
    } else {
      filename <- file.path(out_dir, paste0(name, ".pdf"))
      main <- name
    }
    
    pheatmap::pheatmap(list_scores[[name]],
                       cluster_rows = FALSE,
                       cluster_cols = FALSE,
                       color = colorRampPalette(brewer.pal(n = 9, name =
                                                             "YlGnBu"))(100),
                       display_numbers = TRUE,
                       labels_col = phenotyping_channels,
                       labels_row = basename(rownames(list_scores[[name]])),
                       filename = filename,
                       main = main,
                       number_format = "%.1f",
                       fontsize_number = 8,
                       number_color = "black",
                       width = 10)
  }
  
  if(is.null(batch)){
    saveRDS(list_scores, file.path(out_dir, "AOF_scores_and_Scaled_AOF_scores.RDS"))
  } else {
    saveRDS(list_scores, file.path(out_dir, paste0(batch, "_AOF_scores_and_Scaled_AOF_scores.RDS")))
  }
  
  return(df)
}

#' aof_scoring
#' 
#' @description Calculates AOF (Average Overlap Frequency) scores 
#'
#' @param fcs_files Character, full path to fcs_files.
#' @param phenotyping_markers Character vector, marker names to be used for 
#' clustering, can be full marker name e.g. "CD45" or "CD" if all CD-markers 
#' needs to be plotted
#' @param fsom FlowSOM object as generated by fsom_aof
#' @param out_dir Character, pathway to where the AOF scores and plots should
#' be saved, default is set to working directory.
#' @param batch Character, aqusition batch for each fcs file.
#' Pass to AOF plot names, defult is set to NULL
#' 
#' @return returns data frame with the scaled AOF scores

aof_scoring <- function(fcs_files,
                        phenotyping_markers,
                        fsom,
                        out_dir,
                        batch = NULL){
  
  if(check(phenotyping_channels) == 0){
    
    o <- capture.output(ff_tmp <- read.FCS(file.path(files[1])))
    markers <- FlowSOM::GetMarkers(ff_tmp, colnames(ff_tmp))
    phenotyping_channels <- grep(paste(phenotyping_markers, 
                                       collapse = ("|")), markers, value = TRUE)
    
    if(length(grep("Ir", phenotyping_channels)) > 1){
      phenotyping_channels <- phenotyping_channels[-(grep("Ir", 
                                                          phenotyping_channels)[2])]
    }
  }
  
  aof_scores <- matrix(NA,
                       nrow = length(fcs_files), 
                       ncol = length(phenotyping_channels),
                       dimnames = list(fcs_files,  
                                       names(phenotyping_channels)))
  
  for(file in fcs_files){ 
    print(paste("calculating AOF", file))
    File_ID <- which(fcs_files == file)
    idx <- which(fsom$data[,"File"] == File_ID)
    fcs_data <- fsom$data[idx,]
    MC <- fsom$metaclustering[fsom$map$mapping[idx, 1]]
    
    aof_tmp <- cytutils::greedyCytometryAof(fcs_data = fcs_data, 
                                            y = MC, 
                                            channel_names = names(phenotyping_channels),
                                            width = 0.05, 
                                            cofactor = 5, 
                                            verbose = TRUE) 
    aof_scores[file, ] <- aof_tmp$Aof
  }
  
  scaled_aof_score(aof_scores = aof_scores, 
                   out_dir = out_dir,
                   aof_channels = phenotyping_channels, 
                   batch = batch)
}

#' file_outlier_detecion
#' 
#' @description Detects outlier files based on sample AOF score, generates plot 
#' for outlier and .csv file which indicates which fcs files could be discarded 
#' from further analysis
#'
#' @param scores list of scaled scores per aqusition batch or data frame of scaled scores,
#' both generated by scaled_aof_score or aof_scoring function 
#' @param out_dir Character, pathway to where the plot and .csv files with files
#' quality scores should be saved, default is set to getwd(). 
#' @param sd how many standard deviation should be use to detect outliers  
#' only if argument to_plot = TRUE, default is set to working directory
#' 
#' @return plots Quality AOF scores for all files and save .RDS and .csv Quality 
#' scores for further analysis, files are saved in out_dir 

file_outlier_detecion <- function(scores, out_dir = getwd(), sd) {
  
  if(!inherits(scores, "data.frame") & !inherits(scores, "list")){
    stop("df scores are neither data frame nor list of the data frames")
  }
  
  if(inherits(scores, "list")){
    df_scores <- do.call(rbind, scores)
  } else {
    df_scores <- scores
  }
  
  df_scores$file_names <- basename(rownames(df_scores))
  
  scores_median <- median(df_scores$sample_scores)
  scores_MAD <- mad(df_scores$sample_scores)
  
  df_scores$quality <- ifelse(df_scores$sample_scores > 
                                (scores_median + sd * scores_MAD),"bad","good")
  
  bad_scores <- sum(df_scores$quality == "bad")
  
  colors <- c("bad" = "red", "good" = "darkblue", "threshold= " = "orange")  
  
  max_score <- max(df_scores$sample_scores)
  max_pctgs <- max_score + (max_score * 0.1)
  
  p <- ggplot2::ggplot(df_scores, aes(x = file_names, y = sample_scores, color = quality)) +
    geom_point(size = 4) +
    scale_colour_manual(values = colors) +
    ylim(-0.5, max_pctgs) + 
    geom_hline(yintercept = scores_median + sd * scores_MAD, 
               linetype = "dashed", color = "darkgreen", size = 1)+
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
  
  ggplot2::ggsave(filename = "Quality_AOF_score.png", plot = p, path = file.path(out_dir))
  
  saveRDS(df_scores, file.path(out_dir, "Quality_AOF_score.RDS"))
  write.csv(df_scores, file = file.path(out_dir, "Quality_AOF_score.csv"))
}

#' file_quality_check
#' 
#' @description wraper function to perform sample quality scoring, it will create 
#' Quality control folder where all the quality plots and files will be stored.
#'
#' @param fcs_files Character, full path to fcs_files.
#' @param file_batch_id Character vector, batch label for each fcs_file, 
#'the order needs to be the same as in fcs_files.
#' @param out_dir Character, pathway to where the plots should be saved, 
#' only if argument to_plot = TRUE, default is set to working directory
#' @param phenotyping_markers Character vector, marker names to be used for 
#' flowsom clustering including DNA marker Iridium and viability staining if available,
#' can be full marker name e.g. "CD45" or pattern "CD" if 
#' all CD-markers needs to be plotted, default is set to NULL, so all the mass 
#' channels will be used 
#' @param arcsine_transform Logical, if the data should be transformed with 
#' arcsine transformation and cofactor 5.
#' @param sd numeric, number of standard deviation allowed for file outlier
#' detection, default = 3.
#' @param nClus numeric, as in FlowSOM, numer of metaclusters to be obtained
#' @param ... arguments to be passed to fsom_aof function for FlowSOM parameter 
#' adjustment
#' 
#' @return plots Quality AOF scores for all files and save .RDS and .csv Quality 
#' scores for further analysis, files are saved in out_dir 

file_quality_check <- function(fcs_files, 
                               file_batch_id = NULL, 
                               out_dir = getwd(), 
                               phenotyping_markers = NULL, 
                               arcsine_transform = TRUE, 
                               sd = 3, 
                               nClus = 10, 
                               ...){
  
  if(!dir.exists(out_dir)) dir.create(out_dir)
  
  if (!is.null(file_batch_id)) {
    scores <- list()
    for (batch in unique(file_batch_id)){
      print(batch)
      
      files <- fcs_files[file_batch_id == batch]
      fsom <- fsom_aof(fcs_files = files, 
                       phenotyping_markers = phenotyping_markers, 
                       out_dir = out_dir, 
                       arcsine_transform = arcsine_transform,
                       nClus = nClus,
                       batch = batch, ...)
    
      scores[[batch]] <- aof_scoring(fcs_files = files, 
                                     phenotyping_markers = phenotyping_markers,
                                     fsom = fsom, out_dir = out_dir, batch = batch)
    }
    
  } else {
    files <- fcs_files
    fsom <- fsom_aof(fcs_files = files, phenotyping_markers = phenotyping_markers, 
                     out_dir = out_dir, arcsine_transform = arcsine_transform, 
                     nClus = nClus,
                     batch = NULL)
    
    scores <- aof_scoring(fcs_files = files, 
                          phenotyping_markers = phenotyping_markers,
                          fsom = fsom, out_dir = out_dir, batch = NULL)
  }
  
  final_score <- file_outlier_detecion(scores = scores, out_dir = out_dir, 
                                       sd = sd)
}

#' debarcode_files
#' 
#' @description performs sample debarcoding 
#'
#' @param fcs_files Character, full path to fcs_files.
#' @param file_batch_id Character vector, batch label for each fcs_file, 
#' the order needs to be the same as in fcs_files.
#' @param out_dir Character, pathway to where the plots should be saved, 
#' only if argument to_plot = TRUE, default is set to working directory
#' @param min_threshold logicle, if the minimal threshold for barcoding should be applied. 
#' @param threshold numeric, value for the minimum threshold for debarcoding, 
#' default is set to 0.18 
#' @param to_plot Logical, if plots for yields and debarcoding quality shoudl we drawn
#' @param barcodes_used character, the names of the barcodes that were used, eg.
#' barcode 1 is the same as A1. If NULL ae the barcodes contained in sample_key 
#' from CATALYST package will be used 
#' @param less_than_th logicle, if file names for which debarcoding threshold lower than
#' set in treshold parameter should be saved in out_dir, default FALSE
#' @param barcode_key matrix as in CATALYST::assignPrelim, the debarcoding scheme. 
#' A binary matrix with sample names as row names and numeric masses as column names 
#' OR a vector of numeric masses corresponding to barcode channels. 
#' When the latter is supplied, 'assignPrelim' will create a scheme of the
#' appropriate format internally.
#' 
#' @return save debarcoded fcs files in out_dir. If parameter to_plot set 
#' to TRUE save plots for yields and debarcode_quality in out_dir. If less_than_th 
#' set to TRUE save file names with lower than the value set in threshold parameter 
#' are saved into file called "files_with_lower_debarcoding_threshold.RDS" in out_dir 

debarcode_files <- function(fcs_files, 
                            file_batch_id,
                            out_dir = getwd(), 
                            min_threshold = TRUE, 
                            threshold = 0.18, 
                            to_plot = TRUE, 
                            barcodes_used = NULL, 
                            less_than_th = FALSE, 
                            barcode_key = NULL){
  
  if(anyDuplicated(fcs_files) != 0){
    stop("names of fcs files are duplicated")
  }
  
  for (file in fcs_files){
    print(paste0("   ", Sys.time()))
    print(paste0("   Debarcoding ", file))
    ff <- flowCore::read.FCS(file, transformation = FALSE)
     
    file_id <- which(file == fcs_files)
    batch_id <- file_batch_id[file_id]
    
    if(!dir.exists(out_dir)) dir.create(out_dir)
    
    if(!is.null(barcodes_used)){
      if(is.list(barcodes_used)){
        s_key <- barcode_key[rownames(barcode_key) %in% barcodes_list[[batch_id]],]
      } else {
        s_key <- barcode_key[rownames(barcode_key) %in% barcodes_used,]
      }
      
    } else {
      s_key <- barcode_key
    }

    dat <- CATALYST::prepData(ff)
    dat <- CATALYST::assignPrelim(dat, bc_key = s_key)
    rownames(dat)[rowData(dat)$is_bc]
    # table(colData(dat)$bc_id)
    dat <- CATALYST::estCutoffs(dat)
    
    less_than_th <- c()
    if (min_threshold == TRUE){
      if(any(metadata(dat)$sep_cutoffs < threshold)){
        warning(paste0("cutoff lower than 0.18 has been detected for ", basename(file), 
                      ", cutoff will be set to 0.18"))
        less_than_th <- c(less_than_th, basename(file))
      }
      
      id <- metadata(dat)$sep_cutoffs < threshold
      metadata(dat)$sep_cutoffs[id] <- threshold 

    } else {
      if(any(metadata(dat)$sep_cutoffs < threshold)){
        warning(paste0("cutoff lower than ", threshold, " detected for ", basename(file))) 
        less_than_th <- c(less_than_th, basename(file))
      }
    }
    
    id <- is.na(metadata(dat)$sep_cutoffs)
    metadata(dat)$sep_cutoffs[id] <- 1 
  
    if (to_plot == TRUE){
      p <- CATALYST::plotYields(dat, which = rownames(s_key))
      
      pdf(file.path(out_dir, paste(gsub(".fcs", "_yields.pdf", basename(file)))))
      for (name in names(p)){
        print(p[[name]])
      }
      dev.off()
    }
    
    dat <- CATALYST::applyCutoffs(dat)
    
    if (to_plot == TRUE){
      p <- CATALYST::plotEvents(dat, n = 500)
      
      pdf(file.path(out_dir, paste(gsub(".fcs", "_debarcode_quality.pdf", 
                                        basename(file)))))
      for (name in names(p)){
        print(p[[name]])
      }
      dev.off()
    }
    
    dat <- dat[, dat$bc_id !=0]
    fs <- CATALYST::sce2fcs(dat, split_by = "bc_id")
  
    tmp_dir <- file.path(out_dir, batch_id)
    if(!dir.exists(tmp_dir)) dir.create(tmp_dir)
    
    file_name <- gsub("_cleaned.fcs|.fcs", "", basename(file))
    
    flowCore::write.flowSet(fs, outdir = tmp_dir, 
                  filename = paste0(rownames(fs@phenoData), "_", file_name, 
                                    "_debarcoded.fcs")) 
  }
  
  if(less_than_th == TRUE){
    saveRDS(less_than_th, file.path(out_dir, "files_with_lower_debarcoding_threshold.RDS"))
  }
}

#' aggregate_files
#' 
#' @description performs aggregation of debarcoded files
#' 
#' @param fcs_files Character, full path to fcs_files.
#' @param channels_to_keep Character vector with channel names to be kept. 
#' Default NULL
#' @param outputFile Character, the names of the file that will be given
#' @param maxcells Numeric, maximum cells to randomly aggregate from each file, 
#' default is set to NULL
#' @param write_agg_file Logicle, if the fcs files should be saved, if TRUE 
#' files will be saved in getwd(). Default set to FALSE
#' @param out_dir Character, pathway to where the files should be saved, 
#' only if argument to_plot = TRUE, default is set to working directory
#' 
#' @return aggregated flow frame 

aggregate_files <- function(fcs_files, 
                            channels_to_keep = NULL,
                            outputFile = "aggregate.fcs", 
                            maxcells = NULL, 
                            write_agg_file = FALSE,
                            out_dir = getwd()){
  
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
  
  if(write_agg_file == TRUE){
    
     flowCore::write.FCS(x = flowFrame, filename = file.path(out_dir, outputFile), endian = "big")
  }

  return(flowFrame)
}

#' gate_intact_cells
#' 
#' @description Performs gating of intact cells using flowDensity package
#'
#' @param flow_frame Character, full path to fcs_file.
#' @param file_name Character, the file name used only for plotting, if NULL
#' the file name stored in keyword GUID.original will be used, default is set to NULL
#' @param tinypeak_removal1 numeric from 0-1, as in deGate to exclude/include 
#' tiny peaks in the head of the density distribution curve for both Iridium 
#' channels
#' @param tinypeak_removal2 the same as tinypeak_removal1 but for the tail 
#' in the density distribution curve
#' @param alpha1 numeric, 0-1, as in deGate specify the significance of change 
#' in the slope being detected at the head of the density distribution curve
#' @param alpha2 the same as in alpha1 but for the tail of the density distribution curve
#' @param arcsine_transform Logical, if the data should be transformed 
#' with arcsine transformation and cofactor 5.
#' @param ... 
#' 
#' @return flow frame with intact cells

gate_intact_cells <- function(flow_frame, 
                              file_name = NULL,
                              tinypeak_removal1 = 0.8,
                              tinypeak_removal2 = 0.8,
                              alpha1 = 0.05,
                              alpha2 = 0.1, 
                              arcsine_transform = TRUE, ...){
  
  ff <- flow_frame
  
  if (is.null(file_name)){
    file_name <- ff@description$GUID.original
  } else {
    file_name 
  }
  
  if(arcsine_transform == TRUE){
    
    ff_t <- flowCore::transform(ff, 
                                flowCore::transformList(colnames(ff)[grep("Di", colnames(ff))], 
                                                        CytoNorm::cytofTransform))
  } else {
    ff_t <- ff
  }
  
  selection <- matrix(TRUE,
                      nrow = nrow(ff),
                      ncol = 1,
                      dimnames = list(NULL,
                                      c("intact")))
  
  tr <- list()
  for(m in c("Ir193Di", "Ir191Di")){
    
    tr[[m]] <- c(flowDensity::deGate(ff_t, m,
                                     tinypeak.removal = tinypeak_removal1, 
                                     upper = FALSE, use.upper = TRUE,
                                     alpha = alpha1, verbose = F, count.lim = 3, ...), 
                 flowDensity::deGate(ff_t, m,
                                     tinypeak.removal = tinypeak_removal2, 
                                     upper = TRUE, use.upper = TRUE,
                                     alpha = alpha2, verbose = F, count.lim = 3, ...)) 
  }
  
  for(m in c("Ir193Di", "Ir191Di")){
    selection[ff_t@exprs[,m] < tr[[m]][1], "intact"] <- FALSE
    selection[ff_t@exprs[,m] > tr[[m]][2], "intact"] <- FALSE
  }
  
  percentage <- (sum(selection)/length(selection))*100
  flowDensity::plotDens(ff_t, c("Ir193Di", "Ir191Di"), 
                        main = paste0(basename(file_name)," ( ", format(round(percentage, 2), 
                                                                        nsmall = 2), "% )"))
  
  abline(h = c(tr[["Ir191Di"]]))
  abline(v = c(tr[["Ir193Di"]]))
  points(ff_t@exprs[!selection[,"intact"], c("Ir193Di", "Ir191Di")], pch = ".")
  
  ff <- ff[selection[,"intact"], ]
  
  return(ff)
}


#' remove_mad_outliers
#'
#' @description detects outliers in the selected channel(s) using MAD 
#' (mean absolute deviation). 
#'
#' @param flow_frame 
#' @param channels character, channel names used for gating, default is set to 
#' "Event_length"
#' @param n_mad numeric, how many MAD should be use to detect outliers 
#' @param mad_f function used to compute deviation, default set to "mad" 
#' @param plot logicle, if to plot the data, default TRUE
#' @param center 
#' @param main character, title of the plot, default set to ""
#' @param ... other arguments to pass plotDens
#'
#' @return matrix with the selected cells 
remove_mad_outliers <- function(flow_frame, 
                                channels = "Event_length", 
                                n_mad = 2,
                                mad_f = mad,
                                plot = TRUE,
                                center = "center",
                                main = "",
                                ...){
  boundaries <- matrix(NA,
                       nrow = 5,
                       ncol = length(channels),
                       dimnames = list(c("median", "center", "mad", "l_lim", "u_lim"),
                                       channels))
  for (channel in channels) {
    x <- flow_frame@exprs[, channel]
    boundaries["median", channel] <- median(x)
    boundaries["center", channel] <- density(x)$x[which.max(density(x)$y)]
    boundaries["mad", channel] <- mad_f(x,
                                        center = boundaries[center, channel] )
    boundaries["l_lim", channel] <- boundaries[center, channel] - n_mad * boundaries["mad", channel]
    boundaries["u_lim", channel] <- boundaries[center, channel] + n_mad * boundaries["mad", channel]
  }
  
  selection <- rep(TRUE, nrow(flow_frame))
  for (channel in channels) {
    selection <- selection & (flow_frame@exprs[, channel] > boundaries["l_lim", channel])
    selection <- selection & (flow_frame@exprs[, channel] < boundaries["u_lim", channel])
  }
  percentage <- (sum(selection)/length(selection))*100
  if (plot) {
    flowDensity::plotDens(flow_frame, 
                          c(channels, "Ir191Di"), 
                          main = paste0(main, " ( ", format(round(percentage, 2), 
                                                            nsmall = 2), "% )"),
                          ...)
    if(length(channels) == 2) {
      points(flow_frame@exprs[!selection, channels], col = "red", pch = ".")
      abline(v = boundaries[c("l_lim", "u_lim"), channels[1]], col = "grey")
      abline(h = boundaries[c("l_lim", "u_lim"), channels[2]], col = "grey")
    } else if(length(channels) == 1) {
      points(flow_frame@exprs[!selection, c(channels, "Ir191Di")], pch = ".")
      abline(v = boundaries[c("l_lim", "u_lim"), channels[1]], col = "grey")
    }
  }
  
  return(selection)
}


#' gate_singlet_cells
#'
#' @param flow_frame flow frame
#' @param channels character, channels name to be used for gating, default is 
#' to Event_length
#' @param arcsine_transform Logical, if the data should be transformed with 
#' arcsine transformation and cofactor 5.
#' @param file_name Character, the file name used only for plotting, if NULL
#' the file name stored in keyword GUID.original will be used, default is set to NULL
#' @param n_mad numeric, number of MADs to detect outliers
#' @param ... arguments to pass to plotDens
#'
#' @return flow frame with singlets 
gate_singlet_cells <- function(flow_frame, 
                               channels = "Event_length", 
                               arcsine_transform = TRUE,
                               file_name = NULL,
                               n_mad = 2,
                               ...){
  
  if (is.null(file_name)){
    file_name <- flow_frame@description$FIL
  } else {
    file_name 
  }
  
  if(arcsine_transform == TRUE){
    
    flow_frame_t <- flowCore::transform(flow_frame, 
                                        flowCore::transformList(colnames(flow_frame)[grep("Di", colnames(flow_frame))], 
                                                                CytoNorm::cytofTransform))
  } else {
    flow_frame_t <- flow_frame
  }
  
  selection <- matrix(TRUE,
                      nrow = nrow(flow_frame),
                      ncol = 1,
                      dimnames = list(NULL,
                                      c("singlets")))  
  
  selection[, "singlets"] <- remove_mad_outliers(flow_frame = flow_frame_t, 
                                                 channels = channels,
                                                 main = paste("Singlets", file_name),
                                                 n_mad = n_mad,
                                                 xlim = c(0, 100), ylim = c(0, 8), ...)
  
  flow_frame <- flow_frame[selection[,"singlets"], ]
  
  return(flow_frame)
  
}

#' gate_live_cells
#' 
#' @description Performs gating of live cells using flowDensity package
#' 
#' @param flow_frame Character, full path to fcs_file.
#' @param file_name Character, the file name used only for plotting, if NULL
#' the file name stored in keyword GUID.original will be used, default is set to NULL
#' @param viability_channel Character, the channel name used for viability staining
#' @param tinypeak_removal_viability, numeric from 0-1, as in deGate to exclude/include 
#' tiny peaks in the tail of the density ditribution curve for both viability channel 
#' @param tinypeak_removal_Iridium the same as tinypeak_removal_viablity but for
#' the head and tail of the density ditribution curve in Iridium channel
#' @param alpha_viability numeric, 0-1, as in deGate specify the significance of change 
#' in the slope of viability channel
#' @param alpha_Iridium the same as in alpha_viability but for the Iridium
#' @param arcsine_transform Logical, if the data should be transformed with 
#' arcsine transformation and cofactor 5.
#' @param ... arguments to pass to plotDens function
#' 
#' @return flow frame with live cells

gate_live_cells <- function(flow_frame, 
                            file_name = NULL,
                            viability_channel,
                            tinypeak_removal_viability = 0.8,
                            alpha_viability = 0.1,
                            tinypeak_removal_Iridium = 0.8,
                            alpha_Iridium = 0.05,
                            arcsine_transform = TRUE, ... ){
  
  ff <- flow_frame
  
  if (is.null(file_name)){
    file_name <- ff@description$FIL
  } else {
    file_name 
  }
  
  if(arcsine_transform == TRUE){
    
    ff_t <- flowCore::transform(ff, 
                                transformList(colnames(ff)[grep("Di", colnames(ff))], 
                                              CytoNorm::cytofTransform))
  } else {
    ff_t <- ff
  }
  
  selection <- matrix(TRUE,
                      nrow = nrow(ff),
                      ncol = 1,
                      dimnames = list(NULL,
                                      c("live")))
  
  
  v_ch <- grep(viability_channel, colnames(ff), value = T)
  
  tr <- list()
  for(m in c("Ir191Di", v_ch)){
    if (m == v_ch) {
      upper = TRUE
      alpha = alpha_viability
      tr[[m]] <- flowDensity::deGate(ff_t, m,
                                     tinypeak.removal = tinypeak_removal_viability, 
                                     upper = upper, use.upper = TRUE,
                                     alpha = alpha, verbose = F, count.lim = 3)
      
    } else {
      alpha = alpha_Iridium
      tr[[m]] <- c(flowDensity::deGate(ff_t, m,
                                       tinypeak.removal = tinypeak_removal_Iridium, 
                                       upper = FALSE, use.upper = TRUE,
                                       alpha = alpha,  verbose = F, count.lim = 3), 
                   flowDensity::deGate(ff_t, m,
                                       tinypeak.removal = tinypeak_removal_Iridium, 
                                       upper = TRUE, use.upper = TRUE,
                                       alpha = alpha, verbose = F, count.lim = 3)) 
      
    }
  }
  
  for(m in c(v_ch, "Ir191Di")){
    if (m == v_ch) {
      selection[ff_t@exprs[,m] > tr[[m]][1], "live"] <- FALSE 
    } else {
      selection[ff_t@exprs[,m] < tr[[m]][1], "live"] <- FALSE
      selection[ff_t@exprs[,m] > tr[[m]][2], "live"] <- FALSE  
    }
  }
  percentage <- (sum(selection)/length(selection))*100
  flowDensity::plotDens(ff_t, c(v_ch, "Ir191Di"), 
                        main = paste0(file_name," ( ", format(round(percentage, 2), nsmall = 2), "% )"),
                        xlim = c(0, 8), ylim = c(0, 8), ...)
  
  abline(h = tr[["Ir191Di"]])
  abline(v = tr[[v_ch]])
  
  points(ff_t@exprs[!selection[,"live"], c(v_ch, "Ir191Di")], pch = ".") 
  
  ff <- ff[selection[,"live"], ]
  
  return(ff)
  
}

#' plot_batch
#' 
#' @description Plots batch effect using UMAP and clustering markers
#'
#' @param files_before_norm Character, full path to the unnormalized fcs_files.
#' @param files_after_norm Character, full path to the normalized fcs_files.
#' @param out_dir Character, pathway to where the plots should be saved, 
#' default is set to working directory.
#' @param clustering_markers Character vector, marker names to be used for clustering,
#' can be full marker name e.g. "CD45" or "CD" if all CD-markers needs to be plotted. 
#' These markers are used for UMAP builduing and plotting
#' @param arcsine_transform Logical, if the data should be transformed with 
#' arcsine transformation and cofactor 5, default is set to TRUE
#' @param batch_pattern Character, bach pattern to be match in the fcs file name
#' @param manual_colors character, vector of the colors to be used, 
#' the number of colors needs to be equal to the length of batch_pattern
#' @param seed seed to be set to obtain reproducible results, 
#' default is set to 789
#' @param cells_total number of cells to plot per each file 
#'
#' @return save plots for batch effect in the out_dir 

plot_batch <- function(files_before_norm , 
                       files_after_norm, 
                       out_dir = getwd(), 
                       clustering_markers = "CD|HLA|IgD|PD|BAFF|TCR", 
                       arcsine_transform = TRUE,
                       batch_pattern = "RUN[0-9]*",
                       manual_colors = NULL, 
                       cells_total = 1000,
                       seed = 789){
  
  files_list <- list("files_before_norm" = files_before_norm, 
                     "files_after_norm" = files_after_norm)
  
  plots <- list()
  for (name in names(files_list)) {
    
    set.seed(seed)
    ff_agg <- FlowSOM::AggregateFlowFrames(fileNames = files_list[[name]],
                                           cTotal = length(files_list[[name]]) * cells_total,
                                           verbose = TRUE,
                                           writeMeta = FALSE,
                                           writeOutput = FALSE,
                                           outputFile = file.path(out_dir, paste0("aggregated_for_batch_plotting.fcs")))
    
    if (arcsine_transform == TRUE){
      ff_agg <- flowCore::transform(ff_agg,
                                    transformList(grep("Di", colnames(ff_agg), value = TRUE), 
                                                  cytofTransform))
    }
    
    markers <- FlowSOM::GetMarkers(ff_agg, colnames(ff_agg))
    
    cl_markers <- paste(clustering_markers, collapse="|")
    cl_markers <- grep(cl_markers, markers, value = T)
    
    ff_agg@exprs[, names(cl_markers)] <- apply(ff_agg@exprs[, names(cl_markers)], 
                                               2, function(x){
      q <- quantile(x, 0.9999)
      x[x > q] <- q
      x
    })
    
    set.seed(seed)
    samp <- length(files_list[[name]])
    ff_samp <- ff_agg@exprs[sample(nrow(ff_agg@exprs), samp*cells_total), ]
    
    set.seed(seed)
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
    
    p <- ggplot2::ggplot(dimred_df,  aes_string(x = "dim1", y = "dim2", color = "batch")) +
      geom_point(aes(color = batch), size = 3, position="jitter") +
      ggtitle(name)+
      guides(color = guide_legend(override.aes = list(size = 5)))+
      theme(panel.background = element_rect(fill = "white", colour = "black",
                                            size = 2, linetype = "solid"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.subtitle = element_text(color="black", size=26, 
                                         hjust = 0.95, face = "bold"),
            axis.text = element_blank(),
            axis.title = element_blank(), 
            axis.ticks = element_blank(),
            strip.text.x = element_text(size = 23, color = "black"),
            strip.background = element_rect(fill = "white"), 
            legend.text = element_text(size = 18), 
            legend.title = element_text(size = 22),
            legend.position = "bottom", 
            # legend.key.size = unit(3,"point"),
            legend.key = element_blank()) 
  
    if (!is.null(manual_colors)){
      p <- p+scale_color_manual(values = manual_colors)
    }  
    
    plots[[name]] <- p
    
  }
  
  png(file.path(norm_dir, "batch_effect.png"), 
      width = length(plots)*1500, 
      height = 1500, res = 300)
  gridExtra::grid.arrange(grobs = plots, ncol = 2)
  dev.off()
}


#' prepare_data_for_plotting
#' 
#' @description construct data frame for plotting cell frequencies and MSI 
#' per clusters and metaclusters obtaind from FlowSOM clustering 
#' 
#' @param frequency_msi_list list containing matrices with cell frequency and msi
#' obstained in step extract_pctgs_msi_per_flowsom
#' @param matrix_type the name of the matrix to be plotted
#' @param seed numeric set to obtain reproducible results, default 654
#' @param n_neighbours The size of local neighborhood in UMAP analysis 
#'
#' @return data frame for plotting 
prepare_data_for_plotting <- function(frequency_msi_list,
                                      matrix_type,
                                      seed = 654,
                                      n_neighbours = 14){
  print(matrix_type)
  
  # set the number of closest neighborhoods 
  n <- n_neighbours
  
  # process files before normalization
  df_b <- frequency_msi_list[["before"]][[matrix_type]]
  
  # select the columns for which MSI is higher than 0.2 SD
  if(grepl("mfi", matrix_type)){
    id_cols <-  which(apply(df_b, 2, sd) > 0.2)
    df_b <- df_b[,id_cols]
  }
  
  # build UMAP for files before the normalization
  set.seed(seed)
  df_b_umap <- data.frame(umap(df_b, n_neighbors = n, scale = T))
  
  # process files after normalization
  df_a <- frequency_msi_list[["after"]][[matrix_type]]
  
  # select the columns for which MSI is higher than 0.2 SD
  if(grepl("mfi", matrix_type)){
    id_cols <-  which(apply(df_a, 2, sd) > 0.2)
    df_a <- df_a[,id_cols]
  }
  
  # build UMAP for files after the normalization
  set.seed(seed)
  df_a_umap <- data.frame(umap(df_a, n_neighbors = n, scale = T))
  
  # extract rownames to use the for ggplot annotation
  rnmes <- c(rownames(df_b), rownames(df_a))
  
  #join two UMAP data frames
  dr <- data.frame(rbind(df_b_umap, df_a_umap), check.names = F)
  colnames(dr) <- c("dim1", "dim2")
  return(dr)
}

#' UMAP
#' 
#' @description Builds UMAP on aggregated flow frame 
#'
#' @param fcs_files Character, full path to fcs_files.
#' @param clustering_markers Character vector, marker names to be used for clustering,
#' can be full marker name e.g. "CD45" or "CD" if all CD-markers needs to be plotted
#' @param functional_markers Character vector, marker names for functional markers 
#' e.g cytokines, phosphorylated proteins etc. The functional markers will be plotted 
#' in the file "Marker distribution across aliquots and batches"
#' @param out_dir Character, pathway to where the plots should be saved, 
#' default is set to working directory.
#' @param batch_pattern Character, bach pattern to be match in the fcs file name
#' @param arcsine_transform arcsine_transform Logical, if the data should be transformed with 
#' arcsine transformation and cofactor 5, default is set to TRUE
#' @param cells_total numeric, number of cells taken from each file to buil UMAP 
#' @param seed numeric set to obtain reproducible results, default is set to 1
#' 
#' @return data frame with UMAP coordinates
 
UMAP <- function(fcs_files, 
                 clustering_markers = c("CD", "HLA", "IgD"),
                 functional_markers = c("IL", "TNF", "TGF", "Gr", "IF"),
                 out_dir = getwd(), 
                 batch_pattern = "day[0-9]*", 
                 arcsine_transform = TRUE, 
                 cells_total = 1000, 
                 seed = 1){
  
  set.seed(seed)
  ff_agg <- FlowSOM::AggregateFlowFrames(fileNames = fcs_files,
                                         cTotal = length(fcs_files) * cells_total,
                                         verbose = TRUE,
                                         writeMeta = TRUE,
                                         writeOutput = TRUE,
                                         outputFile = file.path(out_dir, paste0("aggregated_for_UMAP_analysis.fcs")))
  
  if (arcsine_transform == TRUE){
    ff_agg <- flowCore::transform(ff_agg,
                                  transformList(grep("Di", colnames(ff_agg), 
                                                     value = TRUE), 
                                                CytoNorm::cytofTransform))
  }
  
  markers <- FlowSOM::GetMarkers(ff_agg, colnames(ff_agg))
  
  cl_markers <- paste(clustering_markers, collapse="|")
  cl_markers <- grep(cl_markers, markers, value = T)
  
  ff_agg@exprs[, names(cl_markers)] <- apply(ff_agg@exprs[, names(cl_markers)], 2, function(x){
    q <- quantile(x, 0.9999)
    x[x > q] <- q
    x
  })
  
  set.seed(seed)
  dimred_res <- uwot::umap(X = ff_agg@exprs[, names(cl_markers)], 
                           n_neighbors = 15, scale = TRUE)
  
  if (!is.null(functional_markers)){
    f_markers <- paste(functional_markers, collapse="|")
    f_markers <- grep(f_markers, markers, value = T)
    all_markers <- c(cl_markers, f_markers)
  } else {
    all_markers <- cl_markers
  }
  
  df <- data.frame(dim1 = dimred_res[,1], dim2= dimred_res[,2],
                   ff_agg@exprs[, names(all_markers)])
  
  colnames(df)[3:ncol(df)] <- gsub("_|-", "", 
                                   sub("^[0-9]*[A-Za-z]*|^[0-9]*[A-Za-z]*_","", 
                                       all_markers))
  
  df$file_id <- ff_agg@exprs[,"File2"]
  df$batch <- NA
  df$sample_name <- NA
  
  for (i in 1:length(fcs_files)){
    
    file <- fcs_files[i]
    batch <- stringr::str_match(file, batch_pattern)[,1]
    df[df[, "file_id"] == i, "batch"] <- batch
    sample_id <- gsub(".fcs|_gated.fcs|_CC_gated.fcs","", basename(file))
    df[df[, "file_id"] == i, "sample_name"] <- sample_id
  }
  
  return(df)
  
}

#' manual_labels
#' 
#' @description Get the vector of manual labels for each cell
#'
#' @param manual_matrix matrix with TRUE and FALSE values for each cell and 
#' each population name
#' @param cell_types character, the cells of interest
#'
#' @return vector of manual labels for each cells 

manual_labels <- function (manual_matrix, cell_types) 
{
  if (is.list(manual_matrix)) {
    manual_matrix <- do.call(rbind, manual_matrix)
  }
  manual <- rep("Unknown", nrow(manual_matrix))
  for (cellType in cell_types) {
    manual[manual_matrix[, cellType]] <- cellType
  }
  manual <- factor(manual, levels = c("Unknown", cell_types))
  return(manual)
}

#'@description Calculates breaks for flow frame splitting
make_breaks <- function(event_per_flow_frame, total_events){
  breaks <- split_flowFrames(seq_len(total_events),
                             event_per_flow_frame)
  
  names(breaks) <- seq_along(breaks)
  
  return(list("breaks"=breaks, "events_per_flowframe"=event_per_flow_frame))
}

#'@description Calculates beginning and end of each flow frame
split_flowFrames <- function(vec, seg.length) {
  starts=seq(1, length(vec), by=seg.length)
  ends  = starts + seg.length - 1
  ends[ends > length(vec)]=length(vec)
  
  lapply(seq_along(starts), function(i) vec[starts[i]:ends[i]])
}

#' split_big_flowFrames
#' 
#' @description Split the big flow frames into smaller ones
#'
#' @param flow_frame flow frame
#' @param event_per_flow_frame numeric, the number of events to be split to 
#' small flow frames, default is set to 500000  
#' @param min_cell_per_fcs numeric, minimal number of cells in flow frame 
#' to save fcs file, default 20000
#' @param out_dir Character, pathway to where the files should be saved, 
#' default is set to working directory.
#' 
#' @return save splitted fcs files in out_dir 

split_big_flowFrames <- function(flow_frame, 
                                 event_per_flow_frame = 500000, 
                                 out_dir = getwd(), 
                                 min_cell_per_fcs = 20000){
  
  total_events <- nrow(flow_frame)
  res_breaks <- make_breaks(event_per_flow_frame, total_events)
  
  for (i in as.numeric((names(res_breaks$breaks)))) {
    id <- res_breaks$breaks[[i]]
    
    if (i < 10) {
      num <- paste0("0", i)       
    }
    
    if(nrow(flow_frame[id, ]) > min_cell_per_fcs){
      write.FCS(flow_frame[id, ], 
              file.path(out_dir, gsub(".fcs|.FCS", 
                                      paste0("_", num, ".fcs"), flow_frame@description$ORIGINALGUID)))    
    }
      
  }    
}


#' extract_pctgs_msi_per_flowsom
#' 
#' @description performs FlowSOM clustering and extracts cluster and metacluster 
#' frequency and MSI 
#'
#' @param file_list list, pathway to the files before and after normalization 
#' @param nCells Numeric, number of cells to be cluster per each file, 
#' default is set to 50 000
#' @param phenotyping_markers Character vector, marker names to be used for clustering,
#' can be full marker name e.g. "CD45" or "CD" if all CD-markers needs to be plotted
#' @param functional_markers Character vector, marker names to be used for 
#' functional markers, can be full marker name e.g. "IL-6" or "IL" if all IL-markers needs to be plotted
#' @param xdim Numeric, parameter to pass to FlowSOM, width of the SOM grid, 
#' default is set to 10
#' @param ydim Numeric, parameter to pass to FlowSOM, geight of the SOM grid,
#' default is set to 10
#' @param n_metaclusters Numeric, exact number of clusters for metaclustering 
#' in FlowSOM, default is set to 35
#' @param out_dir Character, pathway to where the FlowSOM clustering plot should
#' be saved, default is set to working directory.
#' @param seed numeric, set to obtain reproducible results, default is set to 789
#' @param arcsine_transform arcsine_transform Logical, if the data should 
#' be transformed with arcsine transformation and cofactor 5, default is set to TRUE
#'  
#' @return list of four matrices that contain calculation for 
#' cl_pctgs (cluster percentages), mcl_pctgs (metaclusters percentages), 
#' cl_msi (cluster MSIs for selected markers), mcl_msi (metaclusters MSI
#' for selected markers)

extract_pctgs_msi_per_flowsom <- function(file_list, 
                                          nCells = 50000, 
                                          phenotyping_markers = c("CD", "HLA", "IgD"), 
                                          functional_markers = NULL,
                                          xdim = 10, 
                                          ydim = 10, 
                                          n_metaclusters = 35,
                                          out_dir = getwd(), 
                                          seed = 789, 
                                          arcsine_transform = TRUE) {
  
  res <- list()
  for (f in names(file_list)){

    nCells <- length(file_list[[f]]) * 50000
    print(paste("aggregating files for", f, "normalization"))
    set.seed(seed)
    ff_agg <- FlowSOM::AggregateFlowFrames(fileNames = file_list[[f]],
                                  cTotal = nCells,
                                  writeOutput = F,
                                  outputFile = file.path(out_dir, paste0(f, "_flowsom_agg.fcs")))
  
    if(arcsine_transform == TRUE){
      ff_aggt <- flowCore::transform(ff_agg, 
                                     transformList(colnames(ff_agg)[grep("Di", colnames(ff_agg))], 
                                                   CytoNorm::cytofTransform))
    } else {
      ff_aggt <- ff_agg
    }
   
    markers <- FlowSOM::GetMarkers(ff_agg, colnames(ff_agg))
    phenotyping_channels <- grep(paste(phenotyping_markers, 
                                       collapse = ("|")), markers, value = TRUE)
    functional_channels <- grep(paste(functional_markers, 
                                      collapse = ("|")), markers, value = TRUE)
    
    # Define parameters for FlowSOM analysis
    xdim <- xdim
    ydim <- ydim
    nClus <- n_metaclusters
    s <- seed
    
    print(paste("building FlowSOM for", f, "normalization"))
    fsom <- FlowSOM::FlowSOM(ff_aggt,
                             colsToUse = names(phenotyping_channels),
                             scale = FALSE,
                             nClus = nClus,
                             seed = s,
                             xdim = xdim,
                             ydim = ydim)
    
    fsomPlot <- FlowSOM::PlotStars(fsom, backgroundValues = fsom$metaclustering)
    fsomTsne <- FlowSOM::PlotDimRed(fsom = fsom, cTotal = 5000, seed = s)
    
    figure <- ggarrange(fsomPlot, fsomTsne,
                        # labels = c("FlowSOM clustering", "tsne"),
                        ncol = 2, nrow = 1)
    
    ggplot2::ggsave(filename = paste0(f, "_FlowSOM.pdf"), plot = figure, device = "pdf", path = out_dir,
           width =24, height = 10)
    
    # Define matrices for frequency (pctgs) calculation and MSI (msi). These calculation is performed 
    # for clusters (cl) and metaclusters (mcl)
    cl_pctgs <- matrix(data = NA, nrow = length(file_list[[f]]),
                       ncol = xdim * ydim,
                       dimnames = list(basename(file_list[[f]]), 1:(xdim*ydim)))
    
    mcl_pctgs <- matrix(data = NA, nrow = length(file_list[[f]]),
                        ncol = nClus,
                        dimnames = list(basename(file_list[[f]]), 1:nClus))
    mfi_cl_names <- apply(expand.grid(paste0("Cl", seq_len(fsom$map$nNodes)),
                                      FlowSOM::GetMarkers(ff_agg, c(phenotyping_channels,functional_channels))),
                          1, paste, collapse = "_")
    mfi_mc_names <- apply(expand.grid(paste0("MC", 1:nClus),
                                      FlowSOM::GetMarkers(ff_agg, c(phenotyping_channels,functional_channels))),
                          1, paste, collapse = "_")
    cl_msi <- matrix(NA,
                     nrow = length(file_list[[f]]),
                     ncol = fsom$map$nNodes * length(names(c(phenotyping_channels,
                                                             functional_channels))),
                     dimnames = list(basename(file_list[[f]]), mfi_cl_names))
    mcl_msi <- matrix(NA,
                      nrow = length(file_list[[f]]),
                      ncol =  length(mfi_mc_names),
                      dimnames = list(basename(file_list[[f]]), mfi_mc_names))
    
    print(paste("calculating frequency and msi for:", f, "normalization"))
    
    for (i in unique(fsom$data[,"File2"])){
      
      file <- basename(file_list[[f]][i])
      
      id <- which(fsom$data[,"File2"] == i)
      fsom_subset <- FlowSOM::FlowSOMSubset(fsom = fsom, ids = id)
      
      cl_counts <- rep(0, xdim * ydim)
      counts_tmp <- table(FlowSOM::GetClusters(fsom_subset))
      cl_counts[as.numeric(names(counts_tmp))] <- counts_tmp
      
      cl_pctgs[file,] <- (cl_counts/sum(cl_counts, na.rm = T))*100
      
      mcl_counts <- tapply(cl_counts, fsom$metaclustering, sum)
      mcl_pctgs[file,] <- tapply(cl_pctgs[file,], fsom$metaclustering, sum)
      
      cluster_mfis <- FlowSOM::GetClusterMFIs(fsom_subset)
      cl_msi[file,] <- as.numeric(cluster_mfis[,names(c(phenotyping_channels,functional_channels))])
      mcluster_mfis <- as.matrix(FlowSOM::GetMetaclusterMFIs(list(FlowSOM = fsom_subset,
                                                     metaclustering = fsom$metaclustering)))
      mcl_msi[file,] <- as.numeric(mcluster_mfis[,names(c(phenotyping_channels,functional_channels))])
      
    }
    
    # impute 0 values for NAs
    mfi_cl_imp <- apply(cl_msi, 2,
                        function(x){
                          missing <- which(is.na(x))
                          x[missing] <- 0
                          x
                        })
    
    mfi_mcl_imp <- apply(mcl_msi, 2,
                         function(x){
                           missing <- which(is.na(x))
                           x[missing] <- 0
                           x
                         })
    
    # store the matrices in the list for convenient plotting 
    all_mx <- list("cl_pctgs" = cl_pctgs,
                   "mcl_pctgs" = mcl_pctgs,
                   "cl_msi" = mfi_cl_imp,
                   "mcl_msi" = mfi_mcl_imp)
    
    saveRDS(all_mx, file.path(out_dir, paste0(f, "_calculated_features.RDS")))
    res[[f]] <- all_mx
  }
  return(res)
}









  