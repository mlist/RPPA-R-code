rppa.proteinConc.overview <-
  function(data.protein.conc, title="", subset.sample=NA)
  { 
    require(manipulate)
    
    manipulate({
      data.protein.conc.copy <- data.protein.conc
      
      normalizeFill <- function(data.protein.conc){
        ddply(data.protein.conc, .(Fill, A, B), transform, 
              concentrations = concentrations / mean(concentrations, na.rm=T),
              upper = upper / mean(concentrations, na.rm=T),
              lower = lower / mean(concentrations, na.rm=T)) 
      }
      
      data.protein.conc.copy <- ddply(normalizeFill(data.protein.conc), .(Sample, Slide, A, B), summarise,
                                        concentrations = mean(concentrations, na.rm=T),
                                        upper = max(upper, na.rm=T),
                                        lower = min(lower, na.rm=T)) 
        each.A <- F
        each.B <- F
        specific.A.copy <- NULL
        specific.B.copy <- NULL
        each.fill <- T
      
      if(duplicate.na)
      {
        foreach(property=c("A","B")) %do%{
          data.protein.conc.copy <- foreach(i=1:nrow(data.protein.conc.copy), .combine=rbind) %do%  {
            if(is.na(data.protein.conc.copy[i,property])){
              foreach(A=levels(data.protein.conc.copy[[property]]), .combine=rbind) %do% {
                currentRow <- data.protein.conc.copy[i,]
                currentRow[[property]] <- A
                return(currentRow)
              }
            }
            else return(data.protein.conc.copy[i,])
          }
        }
      }
      
      if(!normalize.to.ref.sample) reference <- NA
      
      rppa.proteinConc.plot(data.protein.conc.copy, title, swap, horizontal.line, error.bars, scales.free, subset.sample, reference, slideAsFill=T, each.A, each.B, specific.A.copy, specific.B.copy, each.fill)
      
    }, swap = checkbox(FALSE, "Swap category orientation"),
               horizontal.line = checkbox(FALSE, "Draw horizontal line through 1"),
               error.bars = checkbox(TRUE, "Plot error bars"),   
               duplicate.na = checkbox(TRUE, "Duplicate NA values"),
               normalize.to.ref.sample=checkbox(FALSE, "Normalize to reference"),
               reference=picker(as.list(if(!is.na(subset.sample)) subset.sample else c(levels(data.protein.conc$Sample), NA))),
               scales.free=picker("fixed", "free_y", "free_x", "free")      
    )  
  }

rppa.batch.correlation <- function(batch.result)
{
  #create all combinations
  sampleCompare <- expand.grid(unique(batch.result[[1]]$Sample), unique(batch.result[[1]]$Sample))
  sampleCompare$ttest.p <- NA
  sampleCompare$Slide <- NA
  sampleCompare$A <- NA
  sampleCompare$B <- NA
  
  #for each of the slides calculate the t-test all against all samples
  pvalues <- foreach(slide=batch.result) %dopar%
  {
      foreach(A = levels(slide$A)) %do%{
          foreach(B = levels(slide$B)) %do%{
            slide <- subset(slide, A==A & B==B)
            foreach(currRow=1:nrow(sampleCompare), .combine=rbind) %do% {
              sampleCompare[currRow,"ttest.p"] <- t.test(slide[slide[,"Sample"]==sampleCompare[currRow,1], "concentrations"], slide[slide[,"Sample"]==sampleCompare[currRow,2], "concentrations"])$p.value  
              sampleCompare$Slide <- slide$Slide[1]
              sampleCompare$A <- slide$A[1]
              sampleCompare$B <- slide$B[1]
              return(sampleCompare[currRow,])
          }
        }
      }
  }
  
  q <- qplot(x=Var1, y=Var2, data=ldply(pvalues), fill=log(ttest.p), geom="bar", stat="identity", fill=Slide)
  q <- q + theme(axis.text.x = element_text(angle=-45, hjust=0, vjust=1))
  q <- q + scale_fill_gradient2(midpoint = -8)
  print(q)
  
  return(pvalues)
}