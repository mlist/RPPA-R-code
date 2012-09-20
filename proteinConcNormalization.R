rppa.proteinConc.normalize <- function(slideA, slideB, normalize.with.median.first = T, 
                                       target.column="Slide", normalize.per.deposition=F, output.all=F)
{   
  #check if target column is free on both slides
  if(output.all && !is.null(slideA[[target.column]]))
  {
    cat("The target column of slideA is not available.")
    return(NA)
  }  
  if(output.all && !is.null(slideB[[target.column]]))
  {
    cat("The target column of slideB is not available.")
    return(NA)
  }
  slideA[[target.column]] <- attr(slideA, "title")
  slideB[[target.column]] <- attr(slideB, "title")
  
  sub.normalize <- function(slideA){
    slideA$upper <- slideA$upper / median(slideA$concentrations, na.rm=T)
    slideA$lower <- slideA$lower / median(slideA$concentrations, na.rm=T)
    slideA$concentrations <- slideA$concentrations / median(slideA$concentrations, na.rm=T)
    return(slideA)
  }
  
  if(normalize.with.median.first)
  {
    if(!normalize.per.deposition){
      slideA <- sub.normalize(slideA)
      slideB <- sub.normalize(slideB)
    }
    else{
      slideA <- ddply(slideA, .(Deposition), sub.normalize)
      slideB <- ddply(slideB, .(Deposition), sub.normalize)
    }
  }
  
  result <- slideA
  result$concentrations <- slideA$concentrations / slideB$concentrations
  #calculate percentage error, build sum and calculate real error on the new value.
  result$upper <- ((( slideA$upper - slideA$concentrations) /slideA$concentrations ) + (( slideB$upper - slideB$concentrations) / slideB$concentrations )) * result$concentrations 
  result$lower <- (( (slideA$concentrations - slideA$lower) / slideA$concentrations ) + ( (slideB$concentrations - slideB$lower) /slideB$concentrations )) * result$concentrations
  
  result$upper <- result$upper + result$concentrations
  result$lower <- result$concentrations - result$lower
  result[[target.column]] <- paste(slideA[[target.column]], "normalized by", slideB[[target.column]])
  if(output.all) result <- rbind(slideA, result, slideB)
    
  return(result)
}

rppa.specific.dilution <- function(spots, dilution=0.25, deposition=4, ...)
{
  #if(!is.null(spots$Inducer)) spots$Inducer <- gsub(" [0-9]+[.][0-9] mM", "", spots$Inducer )
  
  spots.subset <- subset(spots, DilutionFactor == dilution & Deposition == deposition & !is.na(Signal))
  spots.subset$x.weighted.mean <- spots.subset$Signal
  spots.subset$x.err <- 0
  spots.summarize <- rppa.serialDilution.summarize(spots.subset, ...)
  spots.summarize$x.err <- spots.summarize$sem
  
  if(length(spots.summarize$x.err[is.na(spots.summarize$x.err)])>0) cat("WARNING: some samples have only one valid value and no standard error could be computed!")
  
  
  spots.summarize$concentrations <- spots.summarize$x.weighted.mean
  
  spots.summarize$upper <- spots.summarize$x.weighted.mean + spots.summarize$x.err
  spots.summarize$lower <- spots.summarize$x.weighted.mean - spots.summarize$x.err
  
  spots.summarize <- spots.summarize[,!(colnames(spots.summarize) %in% c("sem", "x.weighted.mean", "x.err", "Deposition"))]
  
  attr(spots.summarize, "title") <- attr(spots, "title")
  attr(spots.summarize, "antibody") <- attr(spots, "antibody")
  
  return(spots.summarize)
}

rppa.normalize.to.ref.sample <- function(data.protein.conc, sampleReference, each.A=F, each.B=F, specific.A, specific.B, method="mean")
{
  require(plyr)
  
  toRefSample <- function(data.protein.conc){ 
    if(is.null(specific.A) && is.null(specific.B)) my.subset <- subset(data.protein.conc, Sample == sampleReference)
    else if(is.null(specific.B)){
      if(!is.na(specific.A)) my.subset <- subset(data.protein.conc, Sample == sampleReference & A == specific.A)
      else my.subset <- subset(data.protein.conc, Sample == sampleReference & is.na(A))
    } 
    else if(is.null(specific.A)){
      if(!is.na(specific.B)) my.subset <- subset(data.protein.conc, Sample == sampleReference & B == specific.B)
      else  my.subset <- subset(data.protein.conc, Sample == sampleReference & is.na(B))
    } 
    else{
      if(!is.na(specific.A) && !is.na(specific.B)) my.subset <- subset(data.protein.conc, Sample == sampleReference & A == specific.A & B == specific.B)
      else if(is.na(specific.A) & !is.na(specific.B)) my.subset <- subset(data.protein.conc, Sample == sampleReference & B == specific.B & is.na(A))
      else if(is.na(specific.B) & !is.na(specific.A)) my.subset <- subset(data.protein.conc, Sample == sampleReference & A == specific.A & is.na(B))
      else my.subset <- subset(data.protein.conc, Sample == sampleReference & is.na(A) & is.na(B))  
    } 
       
    if(method == "mean")  meanOfRefSample <- mean(my.subset$concentrations, na.rm=T)
    else if(method == "mean")  meanOfRefSample <- median(my.subset$concentrations, na.rm=T)
       
    data.protein.conc <- within(data.protein.conc, {
      concentrations <- concentrations / meanOfRefSample  
      upper <- upper / meanOfRefSample
      lower <- lower / meanOfRefSample
    }, meanOfRefSample=meanOfRefSample)
  }
  
  if(each.A && each.B){  
    data.protein.conc <- ddply(data.protein.conc, .(A, B), toRefSample)
  }
  else if(each.A){
    data.protein.conc <- ddply(data.protein.conc, .(A), toRefSample)
  }
  else if(each.B){
    data.protein.conc <- ddply(data.protein.conc, .(B), toRefSample)
  }
  else 
  {
    data.protein.conc <- toRefSample(data.protein.conc)
  }
  
  return(data.protein.conc)
}