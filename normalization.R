rppa.median.normalization <- function(spots){
  if(!is.null(attr(spots, "median.normalized"))) cat("Warning: This slide has been median normalized before!")
  
  spots$Signal = spots$Signal / median(spots$Signal, na.rm=T)
  spots$FG = spots$FG / median(spots$FG, na.rm=T)
  spots$BG = spots$BG / median(spots$BG, na.rm=T)
  attr(spots, "median.normalized") <- TRUE
  
  return(spots)
}

rppa.slides.normalize <- function(slideList, normalization.slide, specific.dilution=NA, median.normalization=T, combineTitle=T)
{
  require(plyr)
  cat("Warning: It usually makes more sense to normalize protein amount estimates (e.g. serial dilution curve)")
  
  return(llply(slideList, rppa.slide.normalize, normalization.slide, specific.dilution, median.normalization, combineTitle))
}

rppa.slide.normalize <- function(spotsA, spotsB, specific.dilution=NA, median.normalization=F, combineTitle=T)
{
  cat("Warning: It usually makes more sense to normalize protein amount estimates (e.g. serial dilution curve)")
  
  if(!is.null(attr(spotsA, "normalized"))) cat("Warning: This slide has been normalized before!")
  
  if(!is.null(attr(spotsA, "title")) && !is.null(attr(spotsB, "title")) && combineTitle)
  {
    attr(spotsA, "title") <- paste(attr(spotsA, "title"), " normalized with ", attr(spotsB, "title"))
  }
  
  if(median.normalization)
  {
    spotsA <- rppa.median.normalization(spotsA)
    spotsB <- rppa.median.normalization(spotsB)
  }
  
  if(!is.na(specific.dilution))
  {
    spotsB <- subset(spotsB, DilutionFactor == specific.dilution)
    spotsB <- ddply(spotsB, .(SampleName,CellLine,Deposition), summarise, Signal=mean(Signal, na.rm=T), FG=mean(FG,na.rm=T), BG=mean(BG, na.rm=T))
    spots.temp<- merge(spotsA, spotsB, by=c("SampleName", "CellLine", "Deposition"), all.x=T)
    
    require(doBy)
    
    spots.temp <- orderBy(~Block+Row+Column, data=spots.temp)
    
    #just make sure both are ordered in the same way
    spotsA <- orderBy(~Block+Row+Column, data=spotsA)
    
    spotsA$Signal <- spots.temp$Signal.x / spots.temp$Signal.y
    spotsA$FG <- spots.temp$FG.x / spots.temp$FG.y
    spotsA$BG <- spots.temp$BG.x / spots.temp$BG.y
  }
    
  else{
    spotsA$Signal = spotsA$Signal / spotsB$Signal
    spotsA$FG = spotsA$FG / spotsB$FG
    spotsA$BG = spotsA$BG / spotsB$BG
  }
  attr(spotsA, "normalized") <- TRUE
  
  return(spotsA)
}