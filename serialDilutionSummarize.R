rppa.serialDilution.summarize <- function(data.protein.conc, method="mean", 
                                          select.columns.sample=c("CellLine"), 
                                          select.columns.A=c("LysisBuffer"),
                                          select.columns.B=c("Inducer"), select.columns.fill=c("Treatment"))
{
  data.protein.conc$x.err.percent <- data.protein.conc$x.err / data.protein.conc$x.weighted.mean
  
  
  if(length(select.columns.sample) > 1) Sample <- data.frame(Sample=apply(data.protein.conc[,select.columns.sample], 1 , paste, collapse=" | "))
  else Sample <- data.frame(Sample=data.protein.conc[,select.columns.sample])
  
  newColumns <- Sample
  
  if(!is.na(select.columns.A))
  {
    if(length(select.columns.A) > 1) A <- apply(data.protein.conc[,select.columns.A], 1 , paste, collapse=" | ")
    else A <- data.protein.conc[,select.columns.A]
    
    newColumns = cbind(newColumns, A)
  }
  
  if(!is.na(select.columns.B))
  {
    if(length(select.columns.B) > 1) B <- apply(data.protein.conc[,select.columns.B], 1 , paste, collapse=" | ")
    else B <- data.protein.conc[,select.columns.B]
    
    newColumns = cbind(newColumns, B)
  }
  
  if(!is.na(select.columns.fill))
  {
    if(length(select.columns.fill) > 1) Fill <- apply(data.protein.conc[,select.columns.fill], 1, paste, collapse=" | ")
    else Fill <- data.protein.conc[,select.columns.fill]
    
    newColumns = cbind(newColumns, Fill)
    fillAttribute <- paste(select.columns.fill, collapse=" | ")
  }
  else{
    cat("Fill attribute cannot be empty in the current version, using Depositions per default.")
    Fill <- data.protein.conc[,"Deposition"]
    newColumns <- cbind(newColumns, Fill)
    fillAttribute <- "Deposition"
  }
  
  result <- cbind(newColumns, data.protein.conc[,c("Deposition", "x.weighted.mean", "x.err.percent")])
  
  if(!is.na(select.columns.A) && !is.na(select.columns.B))
    result <- ddply(result, .(Sample, A, B, Fill, Deposition), summarise, x.weighted.mean=mean(x.weighted.mean, na.rm=T), x.err=sum(x.err.percent), sem=sqrt(var(x.weighted.mean,na.rm=TRUE)/length(na.omit(x.weighted.mean))))
  
  else if(!is.na(select.columns.A) && is.na(select.columns.B))
    result <- ddply(result, .(Sample, A, Fill, Deposition), summarise, x.weighted.mean=mean(x.weighted.mean, na.rm=T), x.err=sum(x.err.percent), sem=sqrt(var(x.weighted.mean,na.rm=TRUE)/length(na.omit(x.weighted.mean))))
  
  else if(!is.na(select.columns.B) && is.na(select.columns.A))
    result <- ddply(result, .(Sample, B, Fill, Deposition), summarise, x.weighted.mean=mean(x.weighted.mean, na.rm=T), x.err=sum(x.err.percent), sem=sqrt(var(x.weighted.mean,na.rm=TRUE)/length(na.omit(x.weighted.mean))))
  
  else result <- ddply(result, .(Sample, Fill, Deposition), summarise, x.weighted.mean=mean(x.weighted.mean, na.rm=T), x.err=sum(x.err.percent), sem=sqrt(var(x.weighted.mean,na.rm=TRUE)/length(na.omit(x.weighted.mean))))
  
  result$x.err <- result$x.weighted.mean * result$x.err
  
  #make sure A, B and sample are factors
  result$Sample <- as.factor(result$Sample)
  if(!is.null(result$A)) result$A <- as.factor(result$A)
  if(!is.null(result$B)) result$B <- as.factor(result$B)
  result$Fill <- as.factor(result$Fill)
  
  if(fillAttribute != "Deposition")
  {
    if(length(unique(result$Deposition)) > 1)
    {
      spots.normalized.depositions <- rppa.normalize.depositions(result)
      spots.unified <- rppa.unify.depositions.mean(spots.normalized.depositions)  
      return(spots.unified)
    }
  }
  
  return(result)
}

rppa.normalize.depositions <- function(data.protein.conc)
{
  data.protein.conc <- within(data.protein.conc, {
    x.weighted.mean <- x.weighted.mean / as.numeric(Deposition)
    x.err <- x.err / as.numeric(Deposition)
    sem <- sem / as.numeric(Deposition)
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
    sem <- sem / as.numeric(deposition.mean)
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
