rppa.set.title <- function(spots, title)
{
  attr(spots, "title") <- title
  
  return (spots)
}

rppa.get.title <- function(spots)
{
  attr(spots, "title")
}