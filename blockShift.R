rppa.vshift.vector <-  function (v, by) 
{
  if (by == 0) return(v)
  
  else if (by >= 0) return(c(NA + 1:by, v[1:(length(v) - by)]))
  
  else return(c(v[(-by + 1):(length(v) - by)]))
}

rppa.hshift <- function(spots)
{
  if(!is.null(attr(spots, "hshifted")))
  {
    cat("This slide has already been vshifted! Do you really want to continue? (yes/no)")
    answer <- readline()
  
    if(answer != "yes") return()
  }
  
  range <- min(spots$Block): max(spots$Block)
  
  for(b in range)
  {
    blockB <- subset(spots, Block == b)
    
    spots[spots$Block==b,] <- unsplit(
      lapply(split(blockB, blockB$Row), function(x)
      {
        by <- unique(x$hshift)
        x <- x[with(x, order(Column)),]
        x$Signal <- rppa.vshift.vector(x$Signal, by)
        x$FG <- rppa.vshift.vector(x$FG, by)
        x$BG <- rppa.vshift.vector(x$BG, by)
        return(x)
      }), blockB$Row)
  }
  
  attr(spots, "hshifted") <- TRUE
  
  return(spots)
}

rppa.vshift <- function(spots, blocks=NA, rows=NA, by=NA)
{
  if(!is.null(attr(spots, "vshifted"))){
    
    cat("This slide has already been vshifted! Do you really want to continue?")
    answer <- readline()
    
    if(answer != "yes") return()
  }
  
  if(is.na(blocks[1])) range <- min(spots$Block):max(spots$Block)
  else range <- blocks
  
  for(b in range)
  {
    blockB <- subset(spots, Block==b);
          
    spots[spots$Block==b,] <- unsplit(
      lapply(split(blockB, blockB$Column), function(x, rows, by)
      {
        if(is.na(by))
          by <- unique(x$vshift)
        
        x <- x[with(x, order(Row)),]  
        
        if(is.na(rows[1])){
          x$Signal <- rppa.vshift.vector(x$Signal, by)
          x$FG <- rppa.vshift.vector(x$FG, by)
          x$BG <- rppa.vshift.vector(x$BG, by)
        }
        
        else{
          for(field in c("Signal", "FG", "BG"))
            x[rows,field] <-  rppa.vshift.vector(x[rows, field], by) 
        }
        
        return(x);
      }, rows=rows, by=by)  
      , blockB$Column)
  }
  attr(spots, "vshifted") <- TRUE
  return(spots)
}