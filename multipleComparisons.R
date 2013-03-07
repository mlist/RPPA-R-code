rppa.tukeyHSD <- function(slide)
{
  foreach(A = unique(slide$A), .combine=rbind) %do%{
    foreach(B = unique(slide$B), .combine=rbind) %do%{
      slide.subset <- subset(slide, A==A & B==B)
      tukey.df <- as.data.frame(TukeyHSD(aov(concentrations ~ Sample, data=slide.subset), ordered=T)$Sample)
      tukey.df$Samples <- row.names(tukey.df)
      tukey.df$A <- A
      tukey.df$B <- B
      tukey.df$slide <- slide$Slide[1]
      return(tukey.df)
    }
  }
}

rppa.dunnett <- function(slide, referenceSample="OTC#3")
{
  require(multcomp)
  foreach(currentA = unique(slide$A), .combine=rbind) %do%{
    foreach(currentB = unique(slide$B), .combine=rbind) %do%{
      slide.subset <- subset(slide, A==currentA & B==currentB)
      slide.subset$Sample <- relevel(slide.subset$Sample, ref=referenceSample)
      slide.aov <- aov(concentrations ~ Sample, data=slide.subset)
      dunnett.df <- with(summary(glht(slide.aov, linfct=mcp(Sample="Dunnett")))$test, 
                         {  data.frame(estimates=coefficients, stderror=sigma, pvalues=pvalues) } )
      dunnett.df$Samples <- row.names(dunnett.df)
      dunnett.df$A <- currentA
      dunnett.df$B <- currentB
      dunnett.df$slide <- slide$Slide[1]
      return(dunnett.df)
    }
  }
}