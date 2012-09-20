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