rppa.compare.heatmap <- function(spotsA, spotsB)
{
  #check if dimensions are compatible
  if(nrow(spotsA) == nrow(spotsB) && max(spotsA$Block) == max(spotsB$Block) && max(spotsA$Column) == max(spotsB$Column) && max(spotsA$Row) == max(spotsB$Row))
  {
    spotsA$Signal <- log2(spotsA$Signal) - log2(spotsB$Signal) 
    
    if(!is.null(attr(spotsA, "title")) && !is.null(attr(spotsB, "title")))
      title <- paste("Comparison of ", attr(spotsA, "title"), " and ", attr(spotsB, "title"))
    else title <- NA
    
    rppa.plot.heatmap(spotsA, title=title, discreteColorA="blue", discreteColorB="red", discreteColorC="white")
  }
  
  else{
    cat("The slide dimensions are not compatible!")
    return()
  }
}