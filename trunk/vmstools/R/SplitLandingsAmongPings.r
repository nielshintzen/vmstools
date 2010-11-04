
SplitLandingsAmongPings <- function(tacsat=tacsat, eflalo2=eflalo2, which.species="COD")
{
#This function uses the output from count pings to distribute the catches of ONE chosen species recorded in eflalo over the pings in tacsat.
#It then adds on the total weight and total landing values of ALL species in eflalo2

#Get index for the chosen species
dn <- dimnames(eflalo2)[[2]]
yp <- NULL
for(ss in which.species){yp <- c(yp,grep(ss,dn)) }

#Sum logbook landings over trip id

eflalo.by.species <- aggregate(list(kgs=eflalo2[,yp[1]],cash=eflalo2[,yp[2]]), list(FT_REF=eflalo2$FT_REF), sum,na.rm=T)

#Match landings to tacsat using trip id

kgs <- eflalo.by.species$kgs[match(tacsat$FT_REF,eflalo.by.species$FT_REF)]

#Divide the landings among the pings

tacsat$kgs.by.ping <- kgs/tacsat$npings

#Match cash value to tacsat using trip id

cash <- eflalo.by.species$cash[match(tacsat$FT_REF,eflalo.by.species$FT_REF)]

#Divide cash values among the pings

tacsat$cash.by.ping <- cash/tacsat$npings

#All species
#KGS
nyp <- grep("_KG_",dn)
eflalo2$kgs <- apply(eflalo2[,nyp],1,sum,na.rm=T)
eflalo.totals  <- aggregate(list(kgs=eflalo2$kgs),list(FT_REF=eflalo2$FT_REF),sum,na.rm=T)
total.kgs <- eflalo.totals$kgs[ match(tacsat$FT_REF,eflalo.totals$FT_REF)]
tacsat$total.kgs.by.ping <- total.kgs/tacsat$npings
#VALUES
#Identify cash columns
nyp <- grep("_EURO_",dn)
#Sum cash for all species
eflalo2$cash <- apply(eflalo2[,nyp],1,sum,na.rm=T)
#Sum over trip
eflalo.totals  <- aggregate(list(cash=eflalo2$cash),list(FT_REF=eflalo2$FT_REF),sum,na.rm=T)
total.cash <- eflalo.totals$cash[ match(tacsat$FT_REF,eflalo.totals$FT_REF)]
tacsat$total.cash.by.ping <- total.cash/tacsat$npings

tacsat

}

#Example

#tacsat1f <- SplitLandingsAmongPings(tacsat=tacsat1f,eflalo2=eflalo2,which.species="COD")