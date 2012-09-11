rppa.reorder <- function(spots, which.column=NA, new.order=NA, add.na=NA){
  
  if(is.na(which.column))
  {
    cat("Which factor do you want to reorder?\n")
    cat(colnames(spots))
    which.column <- as.character(readline())
    
    if(is.null(spots[[which.column]])) return ("invalid factor column selected!")
  }
  
  selected.column <- spots[[which.column]]
  
  if(is.na(add.na))
  {
    cat("Do you want to have <NA> included as factor? yes / no")
    if(readline() == "yes") selected.column <- addNA(selected.column)
  }
  
  if(is.na(new.order))
  {
    cat("Please enter the new order as a vector of indices, example 2,3,1 for order 2, 3, 1. No whitespaces!\n")
    cat("Current order:")
    cat(levels(selected.column))
    
    new.order <- as.numeric(strsplit2(readline(), ","))
    if(!is.vector(new.order)) return ("invalid vector entered!")
    if(length(new.order) != length(levels(selected.column))) return ("vector needs to have the same length as levels exist.") 
  }
  
  spots[[which.column]] <- reorder.factor(selected.column, new.order=new.order)
  
  return(spots)
}