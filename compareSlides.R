rppa.compare.heatmap <- function(spotsA, spotsB)
{
  #check if dimensions are compatible
  if(nrow(spotsA) == nrow(spotsB) && max(spotsA$Block) == max(spotsB$Block) && max(spotsA$Column) == max(spotsB$Column) && max(spotsA$Row) == max(spotsB$Row))
  {
    spotsA$Signal <- log2(spotsA$Signal) - log2(spotsB$Signal) 
    
    if(!is.null(attr(spotsA, "title")) && !is.null(attr(spotsB, "title")))
      title <- paste("Comparison of ", attr(spotsA, "title"), " and ", attr(spotsB, "title"))
    else title <- NA
    
    rppa.plot.heatmap(spotsA, title=title, discreteColorA="blue", discreteColorB="red", discreteColorC="white")
  }
  
  else{
    require(tcltk)
    tkmessageBox(message = "The slide dimensions are not compatible!", icon = "error", type = "ok")
    return()
  }
}

rppa.compare.serialDilution <- function(slideList, normalizeBy=NA, specific.dilution=NA){  
  require(manipulate)
  
  sampleNames <- unique(as.vector(laply(slideList, function(x){return (levels(x$SampleName))})))
  
  manipulate(
    rppa.compare.serialDilution.format(slideList, normalizeBy, specific.dilution, Reference, normalize.depositions, unify.depositions, title="", normalize.each.cellLine, compare),
    Reference=picker(as.list(sampleNames)),
    compare=picker("CellLine", "SampleName"),
    normalize.each.cellLine=checkbox(TRUE, "Normalize celllines to reference sample individually"),
    normalize.depositions=checkbox(TRUE, "Normalize depositions"),
    unify.depositions=checkbox(FALSE, "Unify depositions")
  )
}

rppa.compare.serialDilution.format <- function(slideList, normalizeBy=NA, specific.dilution=NA, sampleReference=NA, normalize.depositions=F, 
                                        unify.depositions=F, title="",
                                        normalize.each.cellLine=F, compare="SampleName"){
  require(plyr)
  
  for(slide in slideList)
  {
    if(is.null(attr(slide, "title"))){
      cat("One or more slides are without title! Please use rppa.set.title to assign a title before comparing multiple slides.")
      return()
    }
  }
  
  data.protein.conc <- ldply(slideList, rppa.serialDilution.format, sampleReference=sampleReference, normalize.depositions=normalize.depositions,
                             unify.depositions=unify.depositions, plot.serial.dilution=F, 
                             normalize.each.cellLine=normalize.each.cellLine, compare=compare, produce.plot=F)
  
  if(length(normalizeBy) > 1)
  {
    if(is.na(specific.dilution))
    {
      normalization.data <- rppa.serialDilution.format(normalizeBy, sampleReference=sampleReference, normalize.depositions=normalize.depositions,
      unify.depositions=unify.depositions, plot.serial.dilution=F, 
      normalize.each.cellLine=normalize.each.cellLine, compare=compare, produce.plot=F)
    
      relative.error <- data.protein.conc$x.err / data.protein.conc$x.weighted.mean + normalization.data$x.err / normalization.data$x.weighted.mean 
      data.protein.conc$x.weighted.mean <- data.protein.conc$x.weighted.mean / normalization.data$x.weighted.mean
      data.protein.conc$x.err <- relative.error * data.protein.conc$x.weighted.mean
    }
    else{
      normalization.data <- subset(normalizeBy, DilutionFactor == specific.dilution)
      if(!unify.depositions && !normalize.depositions)
      {
        normalization.data <- ddply(normalization.data, .(SampleName,CellLine,Deposition), summarise, correction.factor=mean(Signal, na.rm=T))
      }
      else{
        normalization.data <- ddply(normalization.data, .(SampleName,CellLine), summarise, correction.factor=mean(Signal/Deposition, na.rm=T))
      }
      
      if(is.na(sampleReference))
      {
        normalization.data$correction.factor <- normalization.data$correction.factor / mean(normalization.data$correction.factor, na.rm=T)
      }
      
      else
      {
        if(!normalize.each.cellLine){
        normalization.data$correction.factor <- normalization.data$correction.factor / mean(subset(normalization.data, SampleName==sampleReference)$correction.factor, na.rm=T)
        }
        else{
        normalization.data <- ddply(normalization.data, .(CellLine), function(x, sampleReference){ 
            x$correction.factor <- x$correction.factor / mean(subset(x, SampleName==sampleReference)$correction.factor, na.rm=T)
            return(x)
            }, sampleReference=sampleReference)
        }
      }
      
      if(!unify.depositions && !normalize.depositions)
      {
        data.protein.conc <- merge(data.protein.conc, normalization.data, by=c("SampleName", "CellLine", "Deposition"), all.x=T)
      }
      else{
        data.protein.conc <- merge(data.protein.conc, normalization.data, by=c("SampleName", "CellLine"), all.x=T)
      }
      
      data.protein.conc <- within(data.protein.conc, {
        x.weighted.mean <- x.weighted.mean / correction.factor
        x.err <- x.err / correction.factor
      })
      
    }
  }
  
  rppa.serialDilution.plot(data.protein.conc, plot.serial.dilution=F, title=title, normalize.depositions=normalize.depositions, unify.depositions=unify.depositions, 
                          sampleReference, normalize.each.cellLine=normalize.each.cellLine, compare=compare)
  
  return(data.protein.conc)
}
