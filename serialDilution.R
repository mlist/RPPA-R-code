rppa.serialDilution.pairColumns <- function(spots, startColumn=1){
  
  pairedData <- data.frame(x=numeric(0), y=numeric(0))
  
  for(i in startColumn:(ncol(spots)-1))
  {
    pairedData <- rbind(pairedData, cbind(x=spots[,i+1], y=spots[,i]))  
  }
  
  return(pairedData)
  
}

rppa.serialDilution.filter <- function(data, sensible.min, sensible.max)
{
  #filter NA values
  data <- subset(data, !is.na(x) & !is.na(y))
  
  #filter bad log odds ratios
  ratio <- log(abs(data$y)/abs(data$x)) 
  ratio.median <- median(ratio, na.rm=T)
  ratio.mad <- mad(ratio,na.rm=T)
  filter <- abs(ratio - ratio.median)/2 > ratio.mad
  
  data <- data[!filter,] 
  
  #filter non-sensible values
  data <- subset(data, x > sensible.min & x < sensible.max & y > sensible.min & y < sensible.max)
  
  return(data)
}

rppa.serialDilution.compute <- function(spots.c, spots.m, slide.name="unnamed slide", title=NA, normalize.depositions=F, unify.depositions=F, 
                                        sampleReference=NA, plot.serial.dilution=T
                                        , normalize.each.cellLine=F, compare="SampleName", produce.plot=T){
  
  #pair columns for serial dilution plot
  pairedData <- rppa.serialDilution.pairColumns(spots.c, 4)
  
  #starting values for non-linear model fit
  a <- max(5, min(spots.m, na.rm=T))
  M <- max(spots.m, na.rm =T) 
  D0 <- 2
  D <- D0
  minimal.err <- 5
  sensible.min <- 5 
  sensible.max <- 1.e5
  
  #filter non-sensible values
  pairedData <- rppa.serialDilution.filter(pairedData, sensible.min, sensible.max)
  
  #nls model fit
  fit <- nls(y ~ a +1/((1/(x -a) -c)/D+c), data=pairedData, start=list(a=a,D=D,c=1/M),alg="port", lower=list(minimal.err,1,0),weights=1/(minimal.err+abs(as.numeric(x))))
  
  #calculate fitted data
  fittedData <- data.frame(x=pairedData$x, y=predict(fit, data.frame(x=pairedData$x)))
  
  a <- summary(fit)$parameter[1]
  D <- summary(fit)$parameter[2]
  c <- summary(fit)$parameter[3]
  d.a <- summary(fit)$parameter[4]
  d.D <- summary(fit)$parameter[5]
  d.c <- summary(fit)$parameter[6]
  M <- a+1/summary(fit)$parameter[3]
  
  #estimate protein concentrations
  serialDilutionResult <- rppa.serialDilution.protein.con(D0=D0, D=D,c=c,a=a,d.a=d.a, d.D=d.D, d.c=d.c,data.dilutes=spots.m)
  
  data.protein.conc <- rppa.serialDilution.normalize(spots.c, serialDilutionResult,
                               fittedData, title, pairedData, dilutionFactor=D, normalize.depositions, 
                               unify.depositions, sampleReference, plot.serial.dilution, 
                               normalize.each.cellLine, compare)
  
  data.protein.conc$Slide <- slide.name
  
  if(produce.plot)
  {
    rppa.serialDilution.plot(data.protein.conc, fittedData, title, pairedData, 
                           dilutionFactor=D, normalize.depositions, unify.depositions, 
                           sampleReference, plot.serial.dilution, normalize.each.cellLine, compare)
  }
  
  return(data.protein.conc)
}

rppa.serialDilution.normalize <- function(spots.c, serialDilutionResult, fitted, title, pairedData, dilutionFactor, normalize.depositions=F, unify.depositions=F, sampleReference=NA, 
                                          plot.serial.dilution=T, normalize.each.cellLine=F, compare="SampleName"){
  data.protein.conc <- cbind(spots.c[,1:3], serialDilutionResult)
  
  #normalize depositions
  if(normalize.depositions)
  {
    data.protein.conc <- within(data.protein.conc, {
      x.weighted.mean <- x.weighted.mean / as.numeric(Deposition)
      x.err <- x.err / as.numeric(Deposition)
    })
  }
  
  if(unify.depositions)
  {
    data.protein.conc <- ddply(data.protein.conc, .(SampleName, CellLine), summarize, x.weighted.mean=mean(x.weighted.mean), x.err=mean(x.err))
  }
  
  #normalize to reference sample
  if(!is.na(sampleReference))
  {
    toRefSample <- function(data.protein.conc){
      meanOfRefSample <- mean(subset(data.protein.conc, SampleName == sampleReference)$x.weighted.mean, na.rm=T)
      data.protein.conc <- within(data.protein.conc, {
        x.weighted.mean <- x.weighted.mean / meanOfRefSample  
        x.err <- x.err / meanOfRefSample
      }, meanOfRefSample=meanOfRefSample)
    }
    
    if(normalize.each.cellLine){  
      data.protein.conc <- ddply(data.protein.conc, .(CellLine), toRefSample)
    }
    else 
    {
      data.protein.conc <- toRefSample(data.protein.conc)
    }
  }
  
  return (data.protein.conc)
}

rppa.serialDilution.plot <- function(data.protein.conc, fitted=NA, title=NA, pairedData=NA, dilutionFactor=NA, normalize.depositions=F, unify.depositions=F, sampleReference=NA, plot.serial.dilution=T, normalize.each.cellLine=F, compare="SampleName"){
  
  require(ggplot2)
  require(gridExtra)
          
  if(is.na(fitted) || is.na(pairedData))
  {
    plot.serial.dilution <- FALSE
  }
  
  #plot serial dilution curve
  if(plot.serial.dilution)
  {
    serialDilutionPlot <- ggplot(pairedData, aes(x=x, y=y)) + opts(title=paste("Serial Dilution Curve Fit, estimated dilution factor ", round(dilutionFactor, 2))) + xlab("Signal at next dilution step") + ylab("Signal") + geom_point() + geom_line(data=fitted, color="blue") + geom_abline(intercept=0, slope=1, color="red")
  }
  
  #plot protein concentrations  
  limits <- aes(ymax = x.weighted.mean + x.err, ymin= x.weighted.mean - x.err)
  dodge <- position_dodge(width=0.9)
  
  if(!unify.depositions){  
    if(compare=="SampleName")
    {
      p <- qplot(SampleName, x.weighted.mean, data=data.protein.conc, 
               main=title, 
               ylab="Estimated Protein Concentration (Relative Scale)",xlab="Sample", geom="bar", fill=Deposition, position="dodge")
    }
    else
    {
      p <- qplot(CellLine, x.weighted.mean, data=data.protein.conc, 
                 main=title, 
                 ylab="Estimated Protein Concentration (Relative Scale)",xlab="CellLine", geom="bar", fill=Deposition, position="dodge")
    }
  }
  else { 
    if(compare=="SampleName")
    {
      p <- qplot(SampleName, x.weighted.mean, data=data.protein.conc, 
               main=title, 
               ylab="Estimated Protein Concentration (Relative Scale)",xlab="Sample", geom="bar")
    }
    else
    {
      p <- qplot(CellLine, x.weighted.mean, data=data.protein.conc, 
                 main=title, 
                 ylab="Estimated Protein Concentration (Relative Scale)",xlab="CellLine", geom="bar")
    }
  }
  
  if(compare=="SampleName")
  { 
    if(length(unique(data.protein.conc$Slide)) > 1)
    {
      p <- p + facet_grid(Slide~CellLine)
    }
    else
    {
      p <- p + facet_wrap(~CellLine)
    }
    p <- p + opts(axis.text.x=theme_text(angle=-45))
  }
  else if(compare=="CellLine")
  {
    if(length(unique(data.protein.conc$Slide)) > 1)
    {
      p <- p + facet_grid(Slide~SampleName)
    }
    else
    {
      p <- p + facet_wrap(~SampleName)
    }
  }
  
  p <- p + geom_errorbar(limits, width=0.25, position=dodge)
  p <- p + geom_hline(aes(yintercept=1))
  
  if(plot.serial.dilution) {
    sidebysideplot <- grid.arrange(p, serialDilutionPlot, heights=c(3/4, 1/4))
    
    print(sidebysideplot)
  }
  else { print(p) }
  
  return(data.protein.conc)
}

rppa.selectReference <- function(spots)
{
  require(tcltk)
  tt<-tktoplevel()
  tl<-tklistbox(tt,height=10,selectmode="single",background="white")
  tkgrid(tklabel(tt,text="Please select a reference sample!"))
  tkgrid(tl)
  
  sampleNames <- levels(spots$SampleName)
  sampleChoice <<- "Mock"
  
  for (i in (1:length(sampleNames)))
  {
    tkinsert(tl,"end",sampleNames[i])
  }
  tkselection.set(tl,0)  
  
  OnOK <- function()
  {
    sampleChoice <<- sampleNames[as.numeric(tkcurselection(tl))+1]
    tkdestroy(tt)
  }
  OK.but <-tkbutton(tt,text="   OK   ",command=OnOK)
  tkgrid(OK.but)
  tkfocus(tt)
  tkwait.window(tt)
  return(sampleChoice)
}

rppa.serialDilution <- function(spots)
{  
  require(manipulate)
  
  manipulate(
    rppa.serialDilution.format(spots, Reference, normalize.depositions, unify.depositions, plot.serial.dilution, normalize.each.cellLine, compare),
    Reference=picker(as.list(levels(spots$SampleName))),
    compare=picker("CellLine", "SampleName"),
    normalize.each.cellLine=checkbox(TRUE, "Normalize celllines to reference sample individually"),
    normalize.depositions=checkbox(TRUE, "Normalize depositions"),
    unify.depositions=checkbox(FALSE, "Unify depositions"),
    plot.serial.dilution=checkbox(TRUE, "Plot serial dilution curve")
  )
}

rppa.serialDilution.format <- function(spots, sampleReference=NA, normalize.depositions=F, 
                                       unify.depositions=F, plot.serial.dilution=T,
                                       normalize.each.cellLine=F, compare="SampleName", produce.plot=T) {  
  require(reshape)
  
  #sampleReference <- rppa.selectReference(spots)
  
  title <- attr(spots, "title")
  
  #filter NA values
  spots <- spots[!is.na(spots$Dilution),]
  
  #transform continuous into descrete
  spots$Dilution <- as.factor(spots$Dilution)
  
  #reverse order of levels so that first column will be purest concentration 
  spots$Dilution <- factor(spots$Dilution, levels=rev(levels(spots$Dilution)))
  
  #cast into table
  spots.c <- cast(spots, CellLine + SampleName + Deposition ~ Dilution, value="Signal", add.missing=TRUE, fun.aggregate="median", na.rm=T)
  
  #matrix without identifying columns
  spots.m <- as.matrix(spots.c[,4:ncol(spots.c)] )  
  
  slide.name <- "unnamed slide"
  if(!is.null(attr(spots, "title"))) { slide.name <- attr(spots, "title") }
  
  return(rppa.serialDilution.compute(spots.c, spots.m, slide.name, title, normalize.depositions, 
                     unify.depositions, sampleReference, plot.serial.dilution, 
                                   normalize.each.cellLine, compare, produce.plot))
}

rppa.serialDilution.protein.con <- function (D0,D,c,a,d.D,d.c, d.a, data.dilutes,r=1.2,minimal.err=5) {
  #D0 = dilution.factor # this is a preset value in diluteion experiments, typical D0=10, 3, or 2.
  #D fitted dilution factor. Ideally, D = D0 ^ gamma, where gamma is a parameter in Sips model
  # k 1:ncol(data.dilutes)    # index of dilute dilution steps in each dilution sereies
  # Np 1:nrow(data.dilutes)   # index of samples
  x.weighted.mean= rep(NA,nrow(data.dilutes))
  x.err = x.weighted.mean
  xflag = x.weighted.mean    # takes values of 0,1,2, which means under detection, OK, saturated
  K = ncol(data.dilutes) 		# number of total dilution steps for a sample
  igamma = log(D0)/log(D)	#where gamma is 1/gamma, a parameter in Sips model
  M =min(1e9,1/c+a)			#when M is too large, take 1e9.
  
  x.saturation.level=   D0^(K-1)/((1/( M/r - a)- 1/(M-a)))^igamma 
  x.nodetection.level = D0^(1-1)/((1/( r*a - a)- 1/(M-a)))^igamma
  
  for (Np in 1:nrow(data.dilutes)){ # for each sample
    x=rep(NA,K); w=x; xL=x; xH=x;  #initialization
    y = data.dilutes[Np,]
    if((y[K]> M/r) && length(y[y<M/r])<2) {#condition to call saturation
      xflag[Np] = 2 
      x.weighted.mean[Np] = x.saturation.level # Use M/r value
      x.err[Np] = NA
    } else {
      if((y[1]<r*a) & length(y[y>r*a])<2) {#condition to call undetected
        
        xflag[Np] = 0
        x.weighted.mean[Np] = x.nodetection.level # Use r*a value
        x.err[Np] = NA
      } else {
        
        y[y>M/r] = NA # for removing signals near saturation 
        y[y<a*r] = NA # for removing signals near bg noise
        
        for (k in 1:K){# for each signal in a dilution series
          y[k] =max(min(M/1.01,y[k]), a+minimal.err) # limit y[k] to be within a+minimal.err and M/1.01
          x[k] =   D0^(k-1) /(1/(y[k]-a)- c)^igamma #estimated protein concentration prior dilution
          #estimate the derivitives
          de.x.over.de.a = igamma * D0^(k-1)*(1/(y[k]-a)- c)^(-igamma-1)/(y[k]-a)^2
          
          de.x.over.de.c = igamma * D0^(k-1)*(1/(y[k]-a)- c)^(-igamma-1)
          de.x.over.de.D = x[k] *log(1/(y[k]-a)- c) * igamma/D/log(D)/D0^(k-1)
          w[k] = (de.x.over.de.a * d.a)^2 + ( de.x.over.de.c * d.c)^2 + (de.x.over.de.D * d.D)^2
        }
        w = w[!is.na(x)] # removing signals near saturation or bg noise
        x = x[!is.na(x)] # removing signals near saturation or bg noise
        
        if(length(x) > 0 ) {
          x.range = 3* max(1, median(x)*0.01,mad(x)) 
          x.f = (abs(x-median(x)) < x.range) # removing outliers
          #c(x,w)
          x=x[x.f]
          w=w[x.f] # removing outliers
          w= 1/w
          x.weighted.mean[Np] = sum (x*w) /sum(w)
          x.err[Np]=1/sqrt(sum(w))
        }
      }#end of else saturation
    }# end of else below detection 
    
  }#end of for each sample Np
  #return value:
  cbind(x.weighted.mean, x.err,xflag)
}#end of function
