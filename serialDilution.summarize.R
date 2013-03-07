rppa.serialDilution.summarize <-
function(data.protein.conc, method="mean", 
                                          select.columns.sample=c("SampleName"), 
                                          select.columns.A=c("LysisBuffer"),
                                          select.columns.B=c("CellLine"), select.columns.fill=c("Deposition"))
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
      spots.normalized.depositions <- within(result, {
        x.weighted.mean <- x.weighted.mean / as.numeric(Deposition)
        x.err  <- x.err / as.numeric(Deposition) 
      })
      spots.unified <- ddply(spots.normalized.depositions, .(Sample, Fill, A, B), summarise, 
            x.weighted.mean = mean(x.weighted.mean, na.rm=T),
            x.err = mean(x.err, na.rm=T))
      return(spots.unified)
    }
  }
  
  return(result)
}
