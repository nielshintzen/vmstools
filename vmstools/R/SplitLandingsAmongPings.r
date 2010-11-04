
SplitLandingsAmongPings <- function(tacsat=tacsat, eflalo2=eflalo2, which.species="COD")
{
#This function uses the output from count pings to distribute the catches of ONE chosen species recorded in eflalo over the pings in tacsat.
#It then adds on the total weight and total landing values of all species in eflalo2

#First sum logbook landings over the trips



dn <- dimnames(eflalo2)[[2]]
yp <- NULL
for(ss in which.species){yp <- c(yp,grep(ss,dn)) }

data <- aggregate(list(kgs=eflalo2[,yp[1]],cash=eflalo2[,yp[1]]), list(FT_REF=eflalo2$FT_REF), sum,na.rm=T)

kgs <- eflalo2[,yp[1]][match(tacsat$FT_REF,eflalo2$FT_REF)]
tacsat$kgs.by.ping <- kgs/tacsat$npings
cash <- eflalo2[,yp[2]][match(tacsat$FT_REF,eflalo2$FT_REF)]
tacsat$cash.by.ping <- cash/tacsat$npings
#All species
#KGS
yp <- grep("_KG_",dn)
eflalo2$kgs <- apply(eflalo2[,yp],1,sum,na.rm=T)
tacsat$total.kgs <- eflalo2$kgs[ match(tacsat$FT_REF,eflalo2$FT_REF)]
#VALUES
yp <- grep("_EURO_",dn)
eflalo2$cash <- apply(eflalo2[,yp],1,sum,na.rm=T)
tacsat$total.cash <- eflalo2$cash[ match(tacsat$FT_REF,eflalo2$FT_REF)]
tacsat

}

#Example

#tacsat1f <- SplitLandingsAmongPings(tacsat=tacsat1f,eflalo2=eflalo2,which.species="COD")