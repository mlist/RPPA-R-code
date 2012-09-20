rppa.plot.heatmap <- function(spots, log=NA, fill="Signal", plotNA=T, palette=NA, 
                              discreteColorA=NA, discreteColorB=NA, discreteColorC=NA, title=NA){
  
  require(ggplot2)
  require(gridExtra)
  
  if(nrow(spots[!is.na(spots[,fill]),]) == 0){
    cat("There is no information on that property!")
    fill <- "Signal"
  } 
  
  if(is.na(title)){
   if(is.null(attr(spots, "title"))) title <- "" 
   else title <- attr(spots, "title") 
  }
  #fix row and column nums
  spots$Row <- as.integer(spots$Row)
  spots$Column <- as.integer(spots$Column)
  
  #transform continuous into descrete
  spots$Deposition <- as.factor(spots$Deposition)
  spots$Dilution <- as.factor(spots$Dilution)
  
  #reverse order of levels so that first factor will be purest concentration 
  spots$Dilution <- factor(spots$Dilution, levels=rev(levels(spots$Dilution)))
  
  if(!is.na(log))
  {
    if(log=="log2") spots$Signal <- log2(spots[["Signal"]])  
    else if(log=="log10") spots$Signal <- log10(spots[["Signal"]])
  }
  
  p <- qplot(Column, Row, data=spots, main=title)
  
  
  if(plotNA) 
  {
    p <- p + geom_tile(data=subset(spots, !is.na(Signal)), line=0, aes_string(fill = fill));
    p <- p + geom_tile(data=subset(spots, is.na(Signal)), line=0, fill="black");
  }
  else {
    p <- p + geom_tile(data=spots, line=0, aes_string(fill = fill));
  }

  if(fill!="Signal")
  {
    empty <- spots[is.na(spots[[fill]]),]
    spots <- spots[!is.na(spots[[fill]]),]
    
    if(dim(empty)[1] != 0) p <- p + geom_tile(data=empty, line=0, fill="grey");
  }
  
  #how many blocks per row?
  blocksPerRow <- 12
  if(!is.null(attr(spots, "blocksPerRow"))) {
    blocksPerRow <- attr(spots, "blocksPerRow") 
  }
  
  p <- p + coord_cartesian(ylim=c(max(spots$Row)+0.5,0.5));
  p <- p + facet_wrap(~Block, ncol=blocksPerRow);
  p <- p + scale_x_continuous(expand=c(0,0), breaks=seq(1, max(spots$Column), 3)) + scale_y_continuous(expand=c(0,0));
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                panel.margin=unit(0.1, "lines"), panel.margin=unit(0, "lines"), 
                plot.margin=unit(c(1, 1, 0.5, 0.5), "lines"),   
                plot.title=element_text(size=18), strip.background=element_rect(fill="grey90", colour="grey50"))
  
  #colorbrewer color palette
  if(!is.na(palette) && is.factor(spots[[fill]])) p <- p + scale_fill_brewer(palette = palette) 
  
  else if(!is.na(discreteColorA) && !is.na(discreteColorB) && !is.factor(spots[[fill]]))
  {
    if(!is.na(discreteColorC)){
      p <- p + scale_fill_gradient2(low = discreteColorA, high = discreteColorB, mid=discreteColorC);
    }
    else
    {
      p <- p + scale_fill_gradient(low = discreteColorA, high = discreteColorB);
    }
  }
  
  print(p + scale_y_reverse());
}


rppa.heatmap <- function(spots){  
  
  require(manipulate)
  
  manipulate(
  rppa.plot.heatmap(spots, log, fill, plotNA, palette, discreteColorA, discreteColorB, NA, NA)
, log = picker("none", "log2", "log10")
, fill = picker("Signal", "FG", "BG","Deposition", "CellLine", "LysisBuffer", "DilutionFactor", "Inducer", "SpotType", "SpotClass", "SampleName", "SampleType", "TargetGene")
, palette = picker("Set1", "Set2", "Set3", "Accent", "Dark2", "Paired", "Pastel1", "Pastel2")
, plotNA = checkbox(TRUE, "Plot NA values")
, discreteColorA = picker("darkblue", "red", "blue", "steelblue", "magenta", "yellow", "white", "green")
, discreteColorB = picker("red", "darkblue", "blue", "steelblue", "magenta", "yellow", "white", "green")
  )
}

