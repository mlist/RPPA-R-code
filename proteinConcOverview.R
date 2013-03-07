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
        data.protein.conc.copy <- rppa.duplicate.nas(data.protein.conc.copy)
      }
      
      if(!normalize.to.ref.sample) reference <- NA
      
      rppa.proteinConc.plot(data.protein.conc.copy, title, swap, horizontal.line, error.bars, scales.free, subset.sample, reference, slideAsFill=T, each.A, each.B, specific.A.copy, specific.B.copy, each.fill)
      
    }, swap = checkbox(FALSE, "Swap category orientation"),
               horizontal.line = checkbox(FALSE, "Draw horizontal line through 1"),
               error.bars = checkbox(TRUE, "Plot error bars"),   
               duplicate.na = checkbox(TRUE, "Duplicate NA values"),
               normalize.to.ref.sample=checkbox(FALSE, "Normalize to reference"),
               reference=picker(as.list(if(!is.na(subset.sample)[1]) subset.sample else c(levels(data.protein.conc$Sample), NA))),
               scales.free=picker("fixed", "free_y", "free_x", "free")      
    )  
  }