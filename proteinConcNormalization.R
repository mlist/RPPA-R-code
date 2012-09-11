rppa.proteinConc.normalize <- function(slideA, slideB, normalize.with.median.first = T, normalize.per.deposition=T, output.all=F)
{   
  sub.normalize <- function(slideA){
    slideA$x.err <- slideA$x.err / median(slideA$x.weighted.mean, na.rm=T)
    slideA$x.weighted.mean <- slideA$x.weighted.mean / median(slideA$x.weighted.mean, na.rm=T)
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
  result$x.weighted.mean <- slideA$x.weighted.mean / slideB$x.weighted.mean
  #calculate percentage error, build sum and calculate real error on the new value.
  result$x.err <- (( slideA$x.err /slideA$x.weighted.mean ) + ( slideB$x.err / slideB$x.weighted.mean )) * result$x.weighted.mean
  
  result$Slide <- paste(slideA$Slide, "normalized by", slideB$Slide)
  if(output.all) result <- rbind(slideA, result, slideB)
    
  return(result)
}

rppa.normalize.depositions <- function(data.protein.conc)
{
  data.protein.conc <- within(data.protein.conc, {
    x.weighted.mean <- x.weighted.mean / as.numeric(Deposition)
    x.err <- x.err / as.numeric(Deposition)
  })
  
  return(data.protein.conc)
}

rppa.normalize.depositions.mean <- function(data.protein.conc)
{
  require(plyr)
  
  depos.means <- ddply(subset(data.protein.conc, select=c("x.weighted.mean", "Deposition")), .(Deposition), summarise, deposition.mean = mean(x.weighted.mean))  
  
  data.protein.conc <- merge(data.protein.conc, depos.means)
  
  data.protein.conc <- within(data.protein.conc, {
    x.weighted.mean <- x.weighted.mean / deposition.mean
    x.err <- x.err / as.numeric(deposition.mean)
  })
  
  return(data.protein.conc)
}

rppa.unify.depositions.mean <- function(data.protein.conc){
  require(plyr)
  
  if(!is.null(data.protein.conc$A) && !is.null(data.protein.conc$B))
    data.protein.conc <- ddply(data.protein.conc, .(Sample, A, B, Fill), summarize, x.weighted.mean=mean(x.weighted.mean), x.err=mean(x.err))
  
  else if(!is.null(data.protein.conc$A))
    data.protein.conc <- ddply(data.protein.conc, .(Sample, A, Fill), summarize, x.weighted.mean=mean(x.weighted.mean), x.err=mean(x.err))
  
  else if(!is.null(data.protein.conc$B))
    data.protein.conc <- ddply(data.protein.conc, .(Sample, B, Fill), summarize, x.weighted.mean=mean(x.weighted.mean), x.err=mean(x.err))
  
  else data.protein.conc <- ddply(data.protein.conc, .(Sample, Fill), summarize, x.weighted.mean=mean(x.weighted.mean), x.err=mean(x.err))
  
  return(data.protein.conc)
}

rppa.unify.depositions.median <- function(data.protein.conc){
  require(plyr)
  
  if(!is.null(data.protein.conc$A) && !is.null(data.protein.conc$B))
    data.protein.conc <- ddply(data.protein.conc, .(Sample, A, B, Fill), summarize, x.weighted.mean=median(x.weighted.mean), x.err=mean(x.err))
  
  else if(!is.null(data.protein.conc$A))
    data.protein.conc <- ddply(data.protein.conc, .(Sample, A, Fill), summarize, x.weighted.mean=median(x.weighted.mean), x.err=mean(x.err))
  
  else if(!is.null(data.protein.conc$B))
    data.protein.conc <- ddply(data.protein.conc, .(Sample, B, Fill), summarize, x.weighted.mean=median(x.weighted.mean), x.err=mean(x.err))
  
  else data.protein.conc <- ddply(data.protein.conc, .(Sample, Fill), summarize, x.weighted.mean=median(x.weighted.mean), x.err=mean(x.err))
  
  return(data.protein.conc)
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
       
    if(method == "mean")  meanOfRefSample <- mean(my.subset$x.weighted.mean, na.rm=T)
    else if(method == "mean")  meanOfRefSample <- median(my.subset$x.weighted.mean, na.rm=T)
       
    data.protein.conc <- within(data.protein.conc, {
      x.weighted.mean <- x.weighted.mean / meanOfRefSample  
      x.err <- x.err / meanOfRefSample
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