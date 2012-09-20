rppa.selectFactors <- function(spots, factor="SampleName")
{
  require(tcltk)
  
  tt<-tktoplevel()
  
  tl<-tklistbox(tt,height=10,selectmode="multiple",background="white")
  tkgrid(tklabel(tt,text="Please select factors to include."))
  tkgrid(tl)
  spotNames <- levels(spots[[factor]])
  for (i in (1:length(spotNames)))
  {
    tkinsert(tl,"end",spotNames[i])
  }
  tkselection.set(tl,0) rm
  
  OKSelection <- function()
  {
    spotNames <<- as.integer(tkcurselection(tl))
    
    tkdestroy(tt)
  }
  
  OK.but <-tkbutton(tt,text="   OK   ",command=OKSelection)
  
  tkgrid(tklabel(tt,text="    "))
  tkgrid(OK.but)
  tkfocus(tt)
  tkwait.window(tt)
  spotNames <- spotNames + 1
  
  return(spots[spots$SampleName %in% levels(spots$SampleName)[spotNames],])
}