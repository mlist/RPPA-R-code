rppa.superCurve.create.sample.names <- function(spots, select.columns.sample, select.columns.A, select.columns.B, select.columns.fill){
  if(length(select.columns.sample) > 1 )
    Sample <- as.data.frame(apply(spots[,select.columns.sample], 1, paste, collapse=" | "))
  else Sample <- as.data.frame(spots[,select.columns.sample])
  
  if(!is.null(select.columns.A) && !is.na(select.columns.A))
  {
    if(length(select.columns.A) > 1)
      Sample <- cbind(Sample, apply(spots[,select.columns.A], 1, paste, collapse=" | "))
    else Sample <- cbind(Sample, spots[,select.columns.A])
  }  
  
  if(!is.null(select.columns.B) && !is.na(select.columns.B))
  {
    if(length(select.columns.B) > 1)
      Sample <- cbind(Sample, apply(spots[,select.columns.B], 1, paste, collapse=" | "))
    else Sample <- cbind(Sample, spots[,select.columns.B])
  }
  
  if(!is.null(select.columns.fill) && !is.na(select.columns.fill))
  {
    if(length(select.columns.fill) > 1)
      Sample <- cbind(Sample, apply(spots[,select.columns.fill], 1, paste, collapse=" | "))
    else Sample <- cbind(Sample, spots[,select.columns.fill])
  }
  
  return(Sample)
}

rppa.superCurve.parse.data <- function(Sample, spots)
{}
  blocksPerRow <- attr(spots, "blocksPerRow")
  Sample <- apply(Sample, 1, paste, collapse=" # ")
  Sample[spots$SpotType!="Sample"] <- "Control"
  Mean.Net <- spots$Signal
  Mean.Total <- spots$FG
  Vol.Bkg <- spots$BG
  Main.Row <- (spots$Block %/% blocksPerRow) + 1
  Main.Col <- (spots$Block %% blocksPerRow) 
  
  #columns zero are actually columns blocksPerRow
  Main.Row[Main.Col==0] <- Main.Row[Main.Col==0] - 1
  Main.Col[Main.Col==0] <- blocksPerRow  
  Sub.Row <- spots$Row
  Sub.Col <- spots$Column
  
  parsed.data <- data.frame(Main.Row, Main.Col, Sub.Row, Sub.Col, Sample, Mean.Net, Mean.Total, Vol.Bkg)
  
  return(parsed.data)
}

rppa.superCurve.create.rppa <- function(parsed.data, spots)
{
  new.rppa = new("RPPA")
  new.rppa@data <- parsed.data
  new.rppa@file <- attr(spots, "title")
  new.rppa@antibody <- attr(spots, "antibody")
  
  return(new.rppa)
}

rppa.superCurve.create.series <- function(parsed.data, spots)
{
  #we also need to assign dilution series, for now based on depositions
  series <- parsed.data[,c("Main.Row", "Main.Col", "Sub.Col", "Sub.Row")]
  
  #assuming four different individual dilution series in one block top: left/right, bottom: left/right 
  series$Sub.Col <- series$Sub.Col %/% ((max(series$Sub.Col) / 2)+1)
  series$Sub.Row <- series$Sub.Row %/% ((max(series$Sub.Row) / 2)+1)
  
  series <- apply(series, 1, paste, collapse=" ")
  series <- as.factor(series)
  
  return(series)
}

rppa.superCurve.create.df <- function(new.fit, select.columns.A, select.columns.B, select.columns.fill)
{
  new.fit <- as.data.frame(new.fit@concentrations)
  colnames(new.fit) <- c("x.weighted.mean")
  new.fit$x.weighted.mean <- 2^(new.fit$x.weighted.mean)
  new.fit$x.err <- 
  
  new.cols <- strsplit2(row.names(new.fit), " # ")
  new.cols <- as.data.frame(new.cols)
  
  col.temp <- c("Sample")
  if(!is.na(select.columns.A)) col.temp <- c(col.temp, "A")
  if(!is.na(select.columns.B)) col.temp <- c(col.temp, "B")
  if(!is.na(select.columns.fill)) col.temp <- c(col.temp, "Fill")
  
  colnames(new.cols) <- col.temp
  
  #fix NAs
  new.cols[new.cols=="NA"] <- NA
  new.cols <- apply(new.cols, 2, factor)
  
  new.df <- cbind(new.fit, new.cols)
  
  return(new.df)
}

rppa.superCurve <- function(spots, select.columns.sample=c("CellLine"), 
                            select.columns.A="LysisBuffer", select.columns.B="Inducer", 
                            select.columns.fill="Treatment", return.fit.only=F, model="logistic", method="nlrob", ci=T, interactive=T){
  
  require(limma)
  require(SuperCurve)
  
  #check for necessary attributes title, antibody, 
  if(is.null(attr(spots, "title"))) return("Please set attribute 'title' first!")
  if(is.null(attr(spots, "antibody"))) return("Please set attribute 'antibody' first!")
  if(is.null(attr(spots, "blocksPerRow")))return("Please set attribute 'blocksPerRow' first!")
  
  #correct inducer format
  if(length(unique(spots$Inducer)) > 1)
    spots$Inducer <- gsub(" [0-9]+[.][0-9] mM", "", spots$Inducer )
  
  #correct dilution factors
  spots$DilutionFactor <- as.double(spots$DilutionFactor)
  
  #create data object for SuperCurve package  
  
  #create sample name from all selected columns  
  Sample <- rppa.superCurve.create.sample.names(spots, select.columns.sample, select.columns.A, select.columns.B, select.columns.fill)
  parsedData <- rppa.superCurve.parse.data(Sample, spots)
  
  #put the information in a RPPA data object
  new.rppa <- rppa.superCurve.create.rppa(parsedData, spots)
  
  #we need the dilution factors as log2 to a reference point, which we will choose to be undiluted 1.0
  steps <- round(log2(spots$DilutionFactor)) + log2(spots$Deposition)
  #steps[is.na(steps)] <- 0
  
  series <- rppa.superCurve.create.series(parsedData, spots)
  
  new.design <- RPPADesign(new.rppa, steps=steps, series=new.rppa@data$Sample, controls=c("Control"), center=T)
  
  if(interactive){
    image(new.design)
    cat("Here you can see the dilution steps that are assumed to be correct. Press enter to continue.")
    readline()
  }
  
  new.fit <- RPPAFit(new.rppa, new.design, "Mean.Net", ci=ci, method=method, model=model)
  if(return.fit.only) return(new.fit)
  
  plot(new.fit)

  new.df <- rppa.superCurve.create.df(new.fit, select.columns.A, select.columns.B, select.columns.fill)
  #new.df$Slide <- attr(spots, "title")
  
  return(new.df)
}