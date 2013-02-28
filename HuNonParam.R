rppa.nonparam <- function(spots, nrep=1, ...){
  
  spots <- subset(spots, SpotClass=="Sample")
  
  #convert input table so that each dilution is in one column
  spots.c <- rppa.serialDilution.format(spots)
  
  #extract number of different dilutions that are not NA
  numOfDilutions <- length(unique(spots$DilutionFactor[!is.na(spots$DilutionFactor)]))
  
  #calculate matrix of dilutions
  spots.m <- rppa.serialDilution.dilutionMatrix(spots.c, numOfDilutions, highestDilutionFirst=F)
  
  nonpa <- getnonpest(spots.m,3);
  
  #combine estimates with signal information
  spots.result <- cbind(spots.c[,1:(ncol(spots.c)-numOfDilutions)], x.weighted.mean=tabus$pass2, x.err=NA)
  
  spots.summarize <- rppa.serialDilution.summarize(spots.result, ...)
  spots.summarize$concentrations <- 2^spots.summarize$x.weighted.mean
  spots.summarize$upper <- 0
  spots.summarize$lower <- 0
  
  spots.summarize <- spots.summarize[,!(colnames(spots.summarize) %in% c("x.weighted.mean", "x.err"))]
  attr(spots.summarize, "title") <- attr(spots, "title")
  attr(spots.summarize, "antibody") <- attr(spots, "antibody")
  
  return(spots.summarize)
}