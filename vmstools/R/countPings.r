countPings <- function(tacsat=tacsat1f,agg.level='trip')
{
#This function counts the number of pings in the tacsat data for each trip or each ICES statistical rectangle by trip
# and then merges them with that same file.

if(agg.level=="trip"){
npings.by.trip <-  aggregate(tacsat$FT_REF,list(tacsat$FT_REF),length)# The number of pings on each trip
dimnames(npings.by.trip)[[2]] <- c("FT_REF","npings")
tacsat.with.pings <- merge(tacsat,npings.by.trip,all=T)
}

if(agg.level=="statsq"){
#Make sure there is an ices statistical rectangle on tacsat
tacsat$LE_RECT<- ICESrectangle(tacsat) 
npings.by.trip.statsq <-  aggregate(tacsat$FT_REF,list(tacsat$FT_REF,tacsat$LE_RECT),length)# The number of pings in each statistical rectangle on each trip
dimnames(npings.by.trip.statsq)[[2]] <- c("FT_REF","LE_RECT","npings")
tacsat.with.pings <- merge(tacsat,npings.by.trip.statsq,all=T)
}




tacsat.with.pings
}
