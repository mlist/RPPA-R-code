rppa.slide.plot <- function(spots, specific.dilution=NA)
{
  require(ggplot2)
  spots <- subset(spots, !is.na(SampleName) & !is.na(DilutionFactor))
  
  if(!is.na(specific.dilution))
  {
    spots <- subset(spots, DilutionFactor==specific.dilution)
  }
  
  spots$Deposition <- as.factor(spots$Deposition)
  q <- qplot(DilutionFactor, Signal, data=spots, geom="bar", fill=Deposition, position="dodge", stat="summary", fun.y="mean") + facet_grid(CellLine~SampleName)
  print(q)
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