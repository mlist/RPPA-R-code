rppa.surface.normalization <- function(spots, positive.control="IgG 400", deposition=4)
{
  require(gam)
  require(mgcv)
  require(ggplot2)
  spots$BlockColumn <- spots$Block %% 12
  spots$BlockColumn[spots$BlockColumn==0] <- 12
  spots$BlockRow <- ((spots$Block-1) %/% 12)+1
  
  spots.pos.ctrl <- spots
  spots.pos.ctrl[(spots.pos.ctrl$SpotType!=positive.control | spots.pos.ctrl$Deposition!=deposition),"Signal"] <- NA
  spots.gam <- gam(Signal ~ s(BlockRow, BlockColumn), data=subset(spots, SpotType==positive.control & Deposition == 4))
  vis.gam(spots.gam)
  spots$surface<- predict.gam(spots.gam, spots)
  spots$Signal <- spots$Signal / (spots$surface / mean(spots$surface, na.rm=T))
  return(spots)
}
