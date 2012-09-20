rppa.slide.plot <- function(spots, select.columns.sample="CellLine", 
                            select.columns.A="LysisBuffer",
                            select.columns.B="Inducer", select.columns.fill="Deposition")
{
  require(ggplot2)
  spots <- subset(spots, !is.na(DilutionFactor))
  
  limits <- aes(ymax = mean + sd, ymin= mean - sd)
  
  spots$Deposition <- as.factor(spots$Deposition)
  q <- qplot(aes_string(select.columns.sample), Signal, data=spots, geom="bar", fill=aes_string(select.columns.fill), position="dodge", stat="summary", fun.y="mean") + facet_grid(aes_string(select.columns.A)~aes_string(select.columns.B))
  print(q)
}

#normalize depositions
rppa.normalize.depositions <- function(spots)
{
  spots <- within(spots, {
    Signal <- Signal / as.numeric(Deposition)
    FG <- FG / as.numeric(Deposition)
    BG <- BG / as.numeric(Deposition)
  })
  
  return(spots)
}

#normalize for a specific sample
rppa.normalize.sample <- function(spots, sampleReference, signalColumn="Signal", normalize.each.cellLine=T){
  
  toRefSample <- function(spots){
    meanOfRefSample <- mean(subset(spots, SampleName == sampleReference)[[signalColumn]], na.rm=T)
    spots[[signalColumn]] <- spots[[signalColumn]] / meanOfRefSample
    return(spots)
  }
  
  if(normalize.each.cellLine){  
    spots <- ddply(spots, .(CellLine), toRefSample)
  }
  else 
  {
    spots <- toRefSample(spots)
  }
  
  return(spots)
}

rppa.slide.cv.plot <- function(spots, plotSampleNumber=T){
  require(ggplot2)
  require(scales)
  require(gridExtra)
  
  #plot CVs
  q <- qplot(as.factor(Deposition), Signal, data=subset(spots, !is.na(DilutionFactor) & !is.na(Deposition)),  
        xlab="Depositions",
        ylab="CV of Signal",
        geom="bar", 
        fill=SampleName,
        stat="summary", 
        fun.y="rppa.calculate.cv") 
  q <- q + facet_grid(SampleName~DilutionFactor) + opts(legend.title=theme_blank())
  q <- q + scale_y_continuous(labels=percent)
  
  #plot sample size
  smplsize <- qplot(as.factor(Deposition), Signal, data=subset(spots, !is.na(DilutionFactor) & !is.na(Deposition)), 
             main="Sample Size", 
             xlab="Depositions",
             ylab="Samples",
             geom="bar",
             position="dodge",
             fill=SampleName,
             stat="summary", 
             fun.y="length") 
  smplsize <- smplsize + facet_grid(~DilutionFactor) + opts(legend.position="none")
  
  sidebysideplot <- grid.arrange(q, smplsize, heights=c(3/4, 1/4))
  
  print(sidebysideplot)
}

rppa.slide.concentration.plot <- function(spots, title="", select.columns.A="CellLine", select.columns.B="Inducer", x.log2=T, y.log2=T, only.samples=T){
  require(ggplot2)
  spots <- subset(spots, !is.na(Signal))
  if(only.samples) spots <- subset(spots, SpotClass=="Sample")
  spots$Concentration <- spots$DilutionFactor * spots$Deposition
  spots$Deposition <- as.factor(spots$Deposition)
  spots$DilutionFactor <- as.factor(spots$DilutionFactor)
  spots$A <- spots[,select.columns.A]
  spots$B <- spots[,select.columns.B]
  
  q <- qplot(Concentration, Signal, data=spots, main=title, color=Deposition, shape=DilutionFactor) 
  q <- q + geom_point(position = position_jitter(width=0.1))
  q <- q + stat_smooth(aes(group=1), method="loess")
  
  if(!is.null(spots$A) && !is.null(spots$B))
  {
    q <- q + facet_grid(A~B)
  }
  else if(!is.null(spots$A))
  {
    q <- q + facet_grid(~A)
  }
  else if(!is.null(spots$B))
  {
    q <- q + facet_grid(~B)
  }
  
  if(x.log2)
  {
    q <- q + scale_x_continuous(trans="log2", name="log2(Concentration)")
  }
  if(y.log2)
  {
    q <- q + scale_y_continuous(trans="log2", name="log2(Signal)")
  }
  print(q)
}

rppa.calculate.cv <- function(x) { y <- sd(x) / mean(x); return(y)}