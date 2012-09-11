rppa.proteinConc.summarize <- function(data.protein.conc, method="mean", 
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
    Fill <- data.protein.conc[,"Deposition"]
    newColumns <- cbind(newColumns, Deposition)
    fillAttribute <- "Deposition"
  }
  
  result <- cbind(newColumns, data.protein.conc[,c("Deposition", "x.weighted.mean", "x.err.percent")])
  
  if(!is.na(select.columns.A) && !is.na(select.columns.B))
    result <- ddply(result, .(Sample, A, B, Fill, Deposition), summarise, x.weighted.mean=mean(x.weighted.mean), x.err=sum(x.err.percent))
  
  else if(!is.na(select.columns.A) && is.na(select.columns.B))
    result <- ddply(result, .(Sample, A, Fill, Deposition), summarise, x.weighted.mean=mean(x.weighted.mean), x.err=sum(x.err.percent))
  
  else if(!is.na(select.columns.B) && is.na(select.columns.A))
    result <- ddply(result, .(Sample, B, Fill, Deposition), summarise, x.weighted.mean=mean(x.weighted.mean), x.err=sum(x.err.percent))
  
  else result <- ddply(result, .(Sample, Fill, Deposition), summarise, x.weighted.mean=mean(x.weighted.mean), x.err=sum(x.err.percent))
  
  result$x.err <- result$x.weighted.mean * result$x.err
  
  #make sure A, B and sample are factors
  result$Sample <- as.factor(result$Sample)
  if(!is.null(result$A)) result$A <- as.factor(result$A)
  if(!is.null(result$B)) result$B <- as.factor(result$B)
  result$Fill <- as.factor(result$Fill)
  
  if(fillAttribute != "Deposition")
  {
      spots.normalized.depositions <- rppa.normalize.depositions.mean(result)
      spots.unified <- rppa.unify.depositions.median(spots.normalized.depositions)  
      return(spots.unified)
  }
  
  return(result)
}

rppa.proteinConc.manipulate <- function(data.protein.conc, title="")
{
  require(manipulate)
  
  manipulate({
      data.protein.conc.copy <- data.protein.conc
      
      if(normalize.to.ref.sample)
      {
        each.A <- (normalize.each == "each A" || normalize.each == "A and B")
        each.B <- (normalize.each == "each B" || normalize.each == "A and B")
        
        if(normalize.each != "specify"){
          specific.A.copy <- NULL
          specific.B.copy <- NULL
        }
        else{
          specific.A.copy <- specific.A
          specific.B.copy <- specific.B
        }
        
        data.protein.conc.copy <- rppa.normalize.to.ref.sample(data.protein.conc.copy, reference, each.A, each.B, specific.A.copy, specific.B.copy)    
      }
      
      rppa.proteinConc.plot(data.protein.conc.copy, title, swap, horizontal.line, error.bars)
      
    }, swap = checkbox(FALSE, "Swap category orientation"),
       horizontal.line = checkbox(FALSE, "Draw horizontal line through 1"),
       error.bars = checkbox(TRUE, "Plot error bars"),       
       normalize.to.ref.sample=checkbox(FALSE, "Normalize to reference"),
       reference=picker(as.list(c(levels(data.protein.conc$Sample), NA))),
       normalize.each=picker("overall", "each A", "each B", "A and B", "specify"),
       specific.A=picker(as.list(c(levels(data.protein.conc$A), NA))),
       specific.B=picker(as.list(c(levels(data.protein.conc$B), NA)))     
  )  
}

rppa.proteinConc.plot <- function(data.protein.conc, title="", swap=F, horizontal.line=T, error.bars=T){
  
  require(ggplot2)
  require(gridExtra)
  
  #check if data has been summarized
  if(is.null(data.protein.conc$Sample))
  {
    return("You have to summarize the data before you can plot it! rppa.proteinConc.summarize()")
  }
  
  #plot protein concentrations  
  limits <- aes(ymax = x.weighted.mean + x.err, ymin= x.weighted.mean - x.err)
  dodge <- position_dodge(width=0.9)
  
  if(!is.null(data.protein.conc$Deposition)){  
      p <- qplot(Sample, x.weighted.mean, data=data.protein.conc, 
                 main=title, 
                 ylab="Estimated Protein Concentration (Relative Scale)",xlab="Sample", geom="bar", fill=Deposition, position="dodge")
  }
  else if(!is.null(data.protein.conc$Fill))
  {
    p <- qplot(Sample, x.weighted.mean, data=data.protein.conc, 
               main=title, 
               ylab="Estimated Protein Concentration (Relative Scale)",xlab="Sample", geom="bar", fill=Fill, position="dodge")   
    p <- p + guides(fill=guide_legend(title=NULL))
  }
  else { 
      p <- qplot(Sample, x.weighted.mean, data=data.protein.conc, 
                 main=title, 
                 ylab="Estimated Protein Concentration (Relative Scale)",xlab="Sample", geom="bar")
  }
  
  if(!is.null(data.protein.conc$B) && !is.null(data.protein.conc$A)){
    if(swap) p <- p + facet_grid(B~A)
    else p <- p + facet_grid(A~B)
  }
  
  else if(!is.null(data.protein.conc$A))
    p <- p + facet_wrap(~A)
  
  else if(!is.null(data.protein.conc$B))
    p <- p + facet_wrap(~B)
  
  if(error.bars) p <- p + geom_errorbar(limits, width=0.25, position=dodge)
  
  if(horizontal.line)  p <- p + geom_hline(aes(yintercept=1))
  
  p <- p + opts(axis.text.x = theme_text(angle=-45, hjust=0, vjust=1))
  p <- p + opts(plot.margin = unit(c(1,2,1,1), "cm"))
  
  print(p)
}
