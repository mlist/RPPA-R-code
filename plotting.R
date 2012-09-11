rppa.slide.plot <- function(spots, specific.dilution=NA)
{
  require(ggplot2)
  spots <- subset(spots, !is.na(SampleName) & !is.na(DilutionFactor))
  
  if(!is.na(specific.dilution))
  {
    spots <- subset(spots, DilutionFactor==specific.dilution)
  
  
  spots.stats <- ddply(subset(spots, DilutionFactor==specific.dilution), .(CellLine,SampleName,Deposition), summarise, 
        mean=mean(Signal, na.rm=T), 
        median=median(Signal, na.rm=T), 
        sd=sd(Signal, na.rm=T), 
        mad=mad(Signal, na.rm=T))
  }
  
  limits <- aes(ymax = mean + sd, ymin= mean - sd)
  
  spots$Deposition <- as.factor(spots$Deposition)
  q <- qplot(DilutionFactor, Signal, data=spots, geom="bar", fill=Deposition, position="dodge", stat="summary", fun.y="mean") + facet_grid(CellLine~SampleName)
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

rppa.slide.concentration.plot <- function(spots, title=""){
  require(ggplot2)
  spots <- subset(spots, !is.na(SampleName) & !is.na(CellLine))
  spots$realConcentration <- spots$DilutionFactor * spots$Deposition
  spots$Deposition <- as.factor(spots$Deposition)
  spots$DilutionFactor <- as.factor(spots$DilutionFactor)
  
  q <- qplot(realConcentration, Signal, data=spots, main=title, color=Deposition, shape=DilutionFactor) 
  q <- q + geom_point(position = position_jitter(width=0.1))
  q <- q + facet_grid(SampleName~CellLine, scales="free_y")
  q <- q + stat_smooth(aes(group=1), method="loess")
  print(q)
}

rppa.calculate.cv <- function(x) { y <- sd(x) / mean(x); return(y)}