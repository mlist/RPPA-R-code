rppa.proteinConc <- function(data.protein.conc, title="")
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
  
  #plot protein concentrations  
  limits <- aes(ymax = upper, ymin= lower)
  dodge <- position_dodge(width=0.9)
  
  if(!is.null(data.protein.conc$Deposition)){  
      p <- qplot(Sample, concentrations, data=data.protein.conc, 
                 main=title, 
                 ylab="Estimated Protein Concentration (Relative Scale)",xlab="Sample", geom="bar", fill=Deposition, position="dodge")
  }
  else if(!is.null(data.protein.conc$Fill))
  {
    p <- qplot(Sample, concentrations, data=data.protein.conc, 
               main=title, 
               ylab="Estimated Protein Concentration (Relative Scale)",xlab="Sample", geom="bar", fill=Fill, position="dodge")   
    p <- p + guides(fill=guide_legend(title=NULL))
  }
  else { 
      p <- qplot(Sample, concentrations, data=data.protein.conc, 
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
  
  p <- p + theme(axis.text.x = element_text(angle=-45, hjust=0, vjust=1))
  p <- p + theme(plot.margin = unit(c(1,2,1,1), "cm"))
  
  print(p)
}
