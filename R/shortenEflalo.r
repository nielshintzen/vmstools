shortenEflalo <- function(data=eflalo2,which.species = c("ANE","BIB","BLL","COD","DAB","HAD","HER","MAC","NEP","PLE","SOL","WHG") )
{
#Eflalo is an unwieldy format so this is handy if you want to select a few species
dn <- dimnames(data)[[2]]
yp <- NULL
for(ss in which.species){yp <- c(yp,grep(ss,dn)) }
short.eflalo2 <- data[,c(1:26,yp)]
short.eflalo2
}