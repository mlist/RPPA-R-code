rppa.load <- function (slideIndex, baseUrl = "http://localhost:8080/MIRACLE/spotExport/", filter.bad.signals=T, apply.shifts=T) 
{
  require(RJSONIO)
  require(plyr)
  
  #read the data from database
  spots <- ldply(fromJSON(paste(baseUrl, "exportAsJSON/", slideIndex, 
                                sep = ""), simplify = T, nullValue = "NA"))
  
  #replace "NA" with proper NA
  replace.na <- colwise(function(col) { col[col=="NA"] <- NA; return(col) })
  spots <- replace.na(spots)
  
  #reformat column types
  spots$FG <- as.double(spots$FG)
  spots$BG <- as.double(spots$BG)
  spots$Signal <- as.double(spots$Signal)
  spots$Diameter <- as.double(spots$Diameter)
  spots$Flag <- as.double(spots$Flag)
  spots$Block <- as.integer(spots$Block)
  spots$Row <- as.integer(spots$Row)
  spots$Column <- as.integer(spots$Column)
  spots$CellLine <- as.factor(spots$CellLine)
  spots$Treatment <- as.factor(spots$Treatment)
  spots$Inducer <- as.factor(spots$Inducer)
  spots$LysisBuffer <- as.factor(spots$LysisBuffer)
  spots$SampleName <- as.factor(spots$SampleName)
  spots$SampleType <- as.factor(spots$SampleType)
  spots$TargetGene <- as.factor(spots$TargetGene)
  spots$DilutionFactor <- as.double(spots$DilutionFactor)
  
  #add shifts
  if (length(fromJSON(paste(baseUrl, "exportShiftsAsJSON/", 
                            slideIndex, sep = ""))) > 0) {
    shifts <- ldply(fromJSON(paste(baseUrl, "exportShiftsAsJSON/", 
                                   slideIndex, sep = ""), simplify = T, nullValue = NA))
    spots <- merge(spots, shifts, by = "Block", all.x = T)
  }
  else{
    spots$vshift <- 0
    spots$hshift <- 0
  }
  #apply shifts
  spots <- rppa.vshift(spots)
  spots <- rppa.hshift(spots)
  
  #add depositions
  depositionPattern <- scan(paste(baseUrl, "getDepositionPattern/", 
                                  slideIndex, sep = ""), what = "integer")
  depositionPattern <- gsub("\\[", "", depositionPattern)
  depositionPattern <- gsub("\\]", "", depositionPattern)
  depositionPattern <- as.integer(strsplit(depositionPattern, 
                                           ",")[[1]])
  spots$Deposition <- spots$Column%%length(depositionPattern)
  spots$Deposition[spots$Deposition == 0] <- length(depositionPattern)
  spots$Deposition <- depositionPattern[spots$Deposition]
  
  #filter bad signals
  if(filter.bad.signals)
  {
    spots <- rppa.filter.diameter(spots)
    spots <- rppa.filter.flagged(spots)
    spots <- rppa.filter.neg.values(spots)
  }
  
  return(spots)
}

rppa.filter.diameter <- function(spots)
{
  spots$Diameter <- as.double(spots$Diameter)
  spots$Signal[spots$Diameter >= 250] <- NA
  return(spots)
}

rppa.filter.neg.values <- function(spots)
{
  spots$FG <- as.double(spots$FG)
  spots$BG <- as.double(spots$BG)
  spots$Signal[(spots$FG-spots$BG <= 0)] <- NA
  return(spots)
}

rppa.filter.flagged <- function(spots)
{
  spots$Flag <- as.double(spots$Flag)
  spots$Signal[spots$Flag != 0] <- NA
  return(spots)
}