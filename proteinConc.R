rppa.proteinConc <-
  function(data.protein.conc, title="", subset.sample=NA)
  { 
    require(manipulate)
    
    manipulate({
      data.protein.conc.copy <- data.protein.conc
      
      if(normalize.fill)
      {
        data.protein.conc.copy <- ddply(data.protein.conc, .(Fill, A, B), transform, 
                                        concentrations = concentrations / mean(concentrations, na.rm=T),
                                        upper = upper / mean(concentrations, na.rm=T),
                                        lower = lower / mean(concentrations, na.rm=T)) 
      }
      
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
        
      }
      else{
        each.A <- F
        each.B <- F
        specific.A.copy <- NULL
        specific.B.copy <- NULL
        reference <- NA
      }
      
      rppa.proteinConc.plot(data.protein.conc.copy, title, swap, horizontal.line, error.bars, scales.free, subset.sample, reference, each.A, each.B, specific.A.copy, specific.B.copy)
      
    }, swap = checkbox(FALSE, "Swap category orientation"),
               horizontal.line = checkbox(FALSE, "Draw horizontal line through 1"),
               error.bars = checkbox(TRUE, "Plot error bars"),   
               normalize.fill=checkbox(FALSE, "Normalize fill"),
               normalize.to.ref.sample=checkbox(FALSE, "Normalize to reference"),
               reference=picker(as.list(if(!is.na(subset.sample)) subset.sample else c(levels(data.protein.conc$Sample), NA))),
               normalize.each=picker( "A and B", "overall", "each A", "each B", "specify"),
               specific.A=picker(as.list(c(levels(data.protein.conc$A), NA))),
               specific.B=picker(as.list(c(levels(data.protein.conc$B), NA))),  
               scales.free=picker("fixed", "free_y", "free_x", "free")      
    )  
  }