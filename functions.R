find_mass_ch <- function(flow_frame, channels = "Time|Event_length|Center|Offset|Width|Residual"){
  non_mass_ch <- grep(c(channels), 
       colnames(ff), 
       value = TRUE,
       invert = TRUE)
  return(non_mass_ch)
}

clean_flow_rate <- function(flow_frame, out_dir = NULL, plot_rate = TRUE) {
  
  flow_frame@exprs[, "Time"] <- flow_frame@exprs[, "Time"]/100
  
  FlowRateData <- flowAI:::flow_rate_bin(flow_frame, 
                                         timeCh = "Time", 
                                         timestep = 0.01)
  
  FlowRateQC <- flowAI:::flow_rate_check(flow_frame, 
                                         FlowRateData, 
                                         alpha = 0.01, 
                                         use_decomp = TRUE) 
  
  
  if (plot_rate == TRUE){
    png(file.path(out_dir,
                  gsub(".fcs", "_flowAI.png", basename(ff@description$FILENAME))),
        width = 800,
        height = 600)
    p <- flowAI:::flow_rate_plot(FlowRateQC)
    print(p)
    dev.off()
  }
return(FlowRateQC)
}
  