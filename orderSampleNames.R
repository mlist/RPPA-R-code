rppa.reorderFactors <- function(spots, factor="SampleName")
{
  require(tcltk)
  
  tt<-tktoplevel()
  
  tl<-tklistbox(tt,height=10,selectmode="single",background="white")
  tkgrid(tklabel(tt,text="Please change the order."))
  tkgrid(tl)
  spotNames <- levels(spots[[factor]])
  for (i in (1:length(spotNames)))
  {
    tkinsert(tl,"end",spotNames[i])
  }
  tkselection.set(tl,0) 
  
  OKSelection <- function()
  {
    tkdestroy(tt)
  }
  
  UpSelection <- function()
  {
    currentIndex <- as.integer(tkcurselection(tl))
    
    if(currentIndex != 0)
    {
      temp <- spotNames[currentIndex] 
      tkdelete(tl, currentIndex)
      tkinsert(tl, (currentIndex-1), spotNames[currentIndex+1])
             
      spotNames[currentIndex] <<- spotNames[currentIndex+1]
      spotNames[currentIndex+1] <<- temp
    }
    
    tkselection.set(tl, currentIndex-1)
  }
  
  DownSelection <- function()
  {
    currentIndex <- as.integer(tkcurselection(tl))
    
    if(currentIndex != length(spotNames)-1)
    {
      temp <- spotNames[currentIndex+1] 
      tkdelete(tl, currentIndex)
      tkinsert(tl, (currentIndex+1), spotNames[currentIndex+1])
      
      spotNames[currentIndex+1] <<- spotNames[currentIndex+2]
      spotNames[currentIndex+2] <<- temp
    }
    
    tkselection.set(tl, currentIndex+1)
    
  }
  
  UpSelection.but <- tkbutton(tt,text="Up",command=UpSelection)
  DownSelection.but <- tkbutton(tt, text="Down", command=DownSelection)
  OK.but <-tkbutton(tt,text="   OK   ",command=OKSelection)
                              
  tkgrid(UpSelection.but)
  tkgrid(DownSelection.but)
  tkgrid(tklabel(tt,text="    "))
  tkgrid(OK.but)
  tkfocus(tt)
  tkwait.window(tt)
  
  spots[[factor]] <- factor(spots[[factor]], levels=spotNames)
  
  return(spots)
}