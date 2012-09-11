rppa.set.attributes <- function(spots)
{
  cat("Title of the Slide / Experiment?")
  title <- readline()
  
  cat("Antibody?")
  antibody <- readline()
  
  cat("Blocks per row?")
  blocksPerRow <- readline()
  
  attr(spots, "title") <- title
  attr(spots, "blocksPerRow") <- as.numeric(blocksPerRow)
  attr(spots, "antibody") <- antibody
  
  cat("Do you want to apply vertical and horizontal shifts? yes/no")
  if(readline() == "yes")
  {
    spots <- rppa.vshift(spots)
    spots <- rppa.hshift(spots)
  }
  return(spots)
}

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

rppa.set.antibody <- function(spots, antibody)
{
  attr(spots, "antibody") <- antibody
  
  return(spots)
}

rppa.get.antibody <- function(spots)
{
  attr(spots, "antibody")
}