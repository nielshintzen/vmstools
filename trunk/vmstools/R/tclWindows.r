callNumberPeak <- function(){
tt <- tktoplevel()
peaks <- tclVar(5)

f1 <- tkframe(tt)
tkpack(f1, side='top')
tkpack(tklabel(f1, text='peaks: '), side='left')
tkpack(tkentry(f1, textvariable=peaks), side='left')

tkpack(tkbutton(tt, text='Next', command=function() tkdestroy(tt)),
        side='right', anchor='s')

tkwait.window(tt)
return(as.numeric(tclvalue(peaks)))}




callPeakValue <- function(pks){

#-Put default values for peaks
for(iPeaks in 1:pks){
  if(iPeaks == 1) entry1 <- tclVar("-10")
  if(iPeaks == 2) entry2 <- tclVar("-5")
  if(iPeaks == 3) entry3 <- tclVar("0")
  if(iPeaks == 4) entry4 <- tclVar("5")
  if(iPeaks == 5) entry5 <- tclVar("10")
}

#-Create input window
tt <- tktoplevel()
tkwm.title(tt,"Value of peaks")
for(iPeaks in 1:pks){
  if(iPeaks==1) box1 <- tkentry(tt, textvariable=entry1)
  if(iPeaks==2) box2 <- tkentry(tt, textvariable=entry2)
  if(iPeaks==3) box3 <- tkentry(tt, textvariable=entry3)
  if(iPeaks==4) box4 <- tkentry(tt, textvariable=entry4)
  if(iPeaks==5) box5 <- tkentry(tt, textvariable=entry5)
}
#-Create input rows
tkgrid(tklabel(tt,text="value of peaks"),columnspan=pks)
for(iPeaks in 1:pks){
  if(iPeaks==1) tkgrid(tklabel(tt,text=paste("peak",iPeaks)), box1)
  if(iPeaks==2) tkgrid(tklabel(tt,text=paste("peak",iPeaks)), box2)
  if(iPeaks==3) tkgrid(tklabel(tt,text=paste("peak",iPeaks)), box3)
  if(iPeaks==4) tkgrid(tklabel(tt,text=paste("peak",iPeaks)), box4)
  if(iPeaks==5) tkgrid(tklabel(tt,text=paste("peak",iPeaks)), box5)
}

done <- tclVar(0)
eqvar <- tclVar(0)

#-Create submit button
submit.but <- tkbutton(tt, text="submit",command=function()tclvalue(done)<-1)

tkgrid(submit.but)
tkbind(tt, "<Destroy>", function()tclvalue(done)<-2)
tkwait.variable(done)

if(tclvalue(done)=="2") stop("aborted")
tkdestroy(tt)
valPeaks <- numeric()
for(iPks in 1:pks) valPeaks <- paste(valPeaks,tclvalue(get(paste("entry",iPks,sep=""))))
return(valPeaks)}

