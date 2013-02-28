rppa.proteinConc.plot <-
function(data.protein.conc, title="", swap=F, horizontal.line=T, error.bars=T, scales, sample.subset=NA, reference=NA, ...){
  
  require(ggplot2)
  require(gridExtra) 
  
  #subset samples and reorder 
  if(!is.na(sample.subset))
  {
    data.protein.conc <- subset(data.protein.conc, Sample %in% sample.subset)
    data.protein.conc$Sample <- factor(data.protein.conc$Sample, sample.subset)
  }
  #normalize data
  if(!is.na(reference)){
    data.protein.conc <- rppa.normalize.to.ref.sample(data.protein.conc, reference, ...)    
  }
      
  #plot protein concentrations  
  limits <- aes(ymax = upper, ymin= lower)
  dodge <- position_dodge(width=0.9)
  
  if(!is.null(data.protein.conc$Deposition)){  
      p <- qplot(Sample, concentrations, data=data.protein.conc, 
                 main=title, stat="identity", 
                 ylab="Estimated Protein Concentration (Relative Scale)",xlab="Sample", geom="bar", fill=Deposition, position="dodge")
  }
  else if(!is.null(data.protein.conc$Fill))
  {
    p <- qplot(Sample, concentrations, data=data.protein.conc, 
               main=title, stat="identity",
               ylab="Estimated Protein Concentration (Relative Scale)",xlab="Sample", geom="bar", fill=Fill, position="dodge")   
    p <- p + guides(fill=guide_legend(title=NULL))
  }
  else { 
      p <- qplot(Sample, concentrations, data=data.protein.conc, 
                 main=title, stat="identity",
                 ylab="Estimated Protein Concentration (Relative Scale)",xlab="Sample", geom="bar")
  }
  
  if(!is.null(data.protein.conc$B) && !is.null(data.protein.conc$A)){
    if(swap) p <- p + facet_grid(B~A, scales=scales)
    else p <- p + facet_grid(A~B, scales=scales)
  }
  
  else if(!is.null(data.protein.conc$A))
    p <- p + facet_wrap(~A, scales=scales)
  
  else if(!is.null(data.protein.conc$B))
    p <- p + facet_wrap(~B, scales=scales)
  
  if(error.bars) p <- p + geom_errorbar(limits, width=0.25, position=dodge)
  
  if(horizontal.line)  p <- p + geom_hline(aes(yintercept=1))
  
  p <- p + theme(axis.text.x = element_text(angle=-45, hjust=0, vjust=1))
  p <- p + theme(plot.margin = unit(c(1,2,1,1), "cm"))
  
  print(p)
}
