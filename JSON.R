rppa.load <- function(slideIndex, baseUrl="http://localhost:8080/RPPAscanner/spotExport/")
{
  require(RJSONIO)
  require(plyr)
  
  spots <- ldply(fromJSON(paste(baseUrl, "exportAsJSON/", slideIndex, sep=""), simplify=T, nullValue=NA))
  spots$FG <- as.double(spots$FG)
  spots$BG <- as.double(spots$BG)
  spots$Signal <- as.double(spots$Signal)
  spots$Block <- as.integer(spots$Block)
  spots$Row <- as.integer(spots$Row)
  spots$Column <- as.integer(spots$Column)
  spots$CellLine <- as.factor(spots$CellLine)
  spots$Treatment <- as.factor(spots$Treatment)
  spots$Inducer <- as.factor(spots$Inducer)
  spots$LysisBuffer <- as.factor(spots$LysisBuffer)
  
  if(length(fromJSON(paste(baseUrl, "exportShiftsAsJSON/", slideIndex, sep=""))) > 0)
  {
    shifts <- ldply(fromJSON(paste(baseUrl, "exportShiftsAsJSON/", slideIndex, sep=""), simplify=T, nullValue="NA"))
    spots <- merge(spots, shifts, by="Block", all.x=T)
  }
  
  depositionPattern <- scan(paste(baseUrl, "getDepositionPattern/", slideIndex, sep=""), what="integer")
  depositionPattern <- gsub("\\[", "", depositionPattern)
  depositionPattern <- gsub("\\]", "", depositionPattern)
  depositionPattern <- as.integer(strsplit(depositionPattern, ",")[[1]])
  
  spots$Deposition <- spots$Column %% length(depositionPattern)
  spots$Deposition[spots$Deposition==0] <- length(depositionPattern)
  spots$Deposition <- depositionPattern[spots$Deposition]
  
  return(spots)
}

rppa.filter.diameter <- function(spots)
{
  spots$Signal[spots$Diameter >= 250] <- NA
  return(spots)
}

rppa.filter.neg.values <- function(spots)
{
  spots$Signal[(spots$FG-spots$BG <= 0)] <- NA
  return(spots)
}

rppa.filter.flagged <- function(spots)
{
  spots$Signal[spots$Flag != 0] <- NA
  return(spots)
}