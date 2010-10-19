countPings <- function(tacsat=tacsat1f)
{
#This function counts the number of pings in the tacsat data for each trip and then merges them with that same file.
npings.by.trip <-  aggregate(tacsat$FT_REF,list(tacsat$FT_REF),length)
dimnames(npings.by.trip)[[2]] <- c("FT_REF","npings")
tacsat.with.pings <- merge(tacsat,npings.by.trip,all=T)
tacsat.with.pings
}
