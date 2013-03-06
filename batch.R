rppa.batch.sdc <- function(slideList, normalizeTo=NA, surface.normalization=T, plotHeatmaps=T, saveDir=NA, positive.control="IgG 400",
                           swap=F, horizontal.line=T, error.bars=T, scales="free", sample.subset=NA, reference=NA,...)
{
  #requirements
  require(foreach)
  require(plyr)
  require(ggplot2)
  
  #processing serial dilution curve and surface normalization method
  sdc <- function(slide)
  {
    cat(paste("Processing", attr(slide, "title"), "..."))
    if(!is.na(saveDir)){
      #heatmap
      if(plotHeatmaps){
        png(paste(attr(slide, "title"), "- Heatmap.png"), width=1024, height=768)
        rppa.plot.heatmap(slide, title=attr(slide, "title"))
        dev.off()
      }
      
      #write raw data
      write.csv2(slide, paste(attr(slide, "title"), "Raw Data.csv"))
    }
  
    if(surface.normalization){
      cat("Applying surface normalization...")
      if(!is.na(saveDir)) png(paste(attr(slide, "title"), "- Surface.png"), width=1024, height=768)  
      
      slide <- rppa.surface.normalization(slide, positive.control)
      
      if(!is.na(saveDir)) dev.off()
      
      #heatmap after surface normalization
      if(plotHeatmaps && !is.na(saveDir)){
        png(paste(attr(slide, "title"), "- Heatmap After Surface Normalization.png"), width=1024, height=768)
        rppa.plot.heatmap(slide, title=attr(slide, "title"))
        dev.off()
        
        png(paste(attr(slide, "title"), "- Heatmap Showing Surface Effects.png"), width=1024, height=768)
        rppa.plot.heatmap(slide, fill="surface", title=attr(slide, "title"))
        dev.off()
      }
    }
    
    #serial dilution curve plot
    if(!is.na(saveDir)) png(paste(attr(slide, "title"), "- Serial Dilution Curve Fit.png"), width=1024, height=768)
    cat("Serial dilution curve fit...")
    
    #backup title
    slideTitle <- attr(slide, "title")
    slideAntibody <- attr(slide, "antibody")
    
    #call serial dilution curve
    slide <- rppa.serialDilution(slide, ...)
    
    #reset title
    attr(slide, "title") <- slideTitle
    attr(slide, "antibody") <- slideAntibody
  
    if(!is.na(saveDir)){
      #close device for sdc plot
      dev.off()

      #write out estimates found through serial dilution curve fit
      write.csv2(slide, paste(attr(slide, "title"), "Serial Dilution Curve Estimates.csv"))
    } 
    return(slide)
  }
  
  #check if data should be saved
  if(!is.na(saveDir))
  {
    keepWd <- getwd()
    setwd(saveDir)
  }
  
  #serial dilution curve on all slides except normalization slide
  result.sdc <- foreach(slide=slideList) %do% sdc(slide)
  
  #method for protein concentration normalization
  proteinConcNorm <- function(slide.sdc, normalizeTo.sdc){  
    slide.normalized <- rppa.proteinConc.normalize(slide.sdc, normalizeTo.sdc)
    slide.normalized$Slide <- attr(slide.sdc, "antibody")
    
    if(!is.na(saveDir)){
      write.csv2(slide.normalized, paste(attr(slide.sdc, "title"), " - PConcEst normalized by", attr(normalizeTo.sdc, "antibody"), ".csv"))
    }
    return(slide.normalized)
  }
  
  #process normalization slide
  if(!is.na(normalizeTo))
  {
    normalizeTo.sdc <- sdc(normalizeTo)
    cat(paste("Normalizing slides to", attr(normalizeTo, "title"), "..."))
    result.normalized <- foreach(slide.sdc=result.sdc) %do% proteinConcNorm(slide.sdc, normalizeTo.sdc) 
  }
  
  if(!is.na(saveDir)) setwd(keepWd)
  cat("Done!")
  
  if(exists("result.normalized")) result <- result.normalized
  else result <- result.sdc
  
  #method for protein concentration plot
  proteinConc <- function(pConc){
    if(!is.na(saveDir)) png(paste(attr(pConc, "title"), "- Protein Concentration Estimates.png"), width=1024, height=768)
    rppa.proteinConc.plot(pConc, attr(pConc, "title"), swap, horizontal.line, error.bars, scales, sample.subset, reference)
    if(!is.na(saveDir)) dev.off()
  }
  
  #plot protein concentration 
  foreach(pConc=result) %do% proteinConc(pConc)
  
  rppa.proteinConc.overview(ldply(result))
  
  return(result)
}