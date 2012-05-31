rppa.boxPlotComparison <- function(spots, freeScales="free") {
  require(ggplot2)
  title <- attr(spots, "title")
  spots <- subset(spots, !is.na(DilutionFactor))
  p <- qplot(SampleName, Signal, data=spots, main=title)
  p <- p + geom_boxplot(aes(fill = factor(Deposition)))
  p <- p + facet_grid(DilutionFactor~CellLine, scales=freeScales)
  p <- p + opts(axis.text.x=theme_text(angle=-45))
  p <- p + labs(fill = "Deposition")
  print(p);
}