rppa.set.title <- function(spots, title)
{
  attr(spots, "title") <- title
  
  return (spots)
}

rppa.get.title <- function(spots)
{
  attr(spots, "title")
}

rppa.set.blocksPerRow <- function(spots, blocksPerRow)
{
  attr(spots, "blocksPerRow") <- blocksPerRow
  
  return(spots)
}

rppa.get.blocksPerRow <- function(spots)
{
  attr(spots, "blocksPerRow")
}