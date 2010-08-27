mergeEflalo2Tacsat <- function(eflalo2,tacsat){

                          #-Add datim format of time
                          D_DATIM   <- as.POSIXct(paste(eflalo2$FT_DDAT, eflalo2$FT_DTIME, sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M:%S") 
                          L_DATIM   <- as.POSIXct(paste(eflalo2$FT_LDAT, eflalo2$FT_LTIME, sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M:%S") 
                          SI_DATIM  <- as.POSIXct(paste(tacsat$SI_DATE,  tacsat$SI_TIME,   sep=" "), tz="GMT", format="%d/%m/%Y  %H:%M:%S")
                          
                          Ta        <- cbind(tacsat,SI_DATIM)
                          Ef        <- cbind(eflalo2,D_DATIM,L_DATIM)
                          
                          #-Order the data by time
                          Ef        <- orderBy(~VE_REF+D_DATIM+L_DATIM,data=Ef)
                          Ta        <- orderBy(~VE_REF+SI_DATIM,data=Ta)
                          
                          #-Split the Eflalo and Tacsat data by vessel
                          splitEf   <- split(Ef,Ef$VE_REF)
                          splitTa   <- split(Ta,Ta$VE_REF)
                          
                          #-link vessels in tacsat to eflalo
                          tacefmatch <- pmatch(sort(unique(Ta$VE_REF)),sort(unique(Ef$VE_REF)))
                          
                          #-loop over all the vessels in tacsat
                          for(i in 1:length(tacefmatch)){
                            pm    <- tacefmatch[i]
                            eftim <- splitEf[[pm]][which(duplicated(splitEf[[pm]]$FT_REF)==F),c("D_DATIM","L_DATIM","FT_REF")]
                            dtime <- eftim[,1]
                            ltime <- eftim[,2]
                            stime <- splitTa[[i]]$SI_DATIM
                            tripn <- eftim[,3]
                            
                            smdtime <- t(outer(stime,dtime,"-"))
                            gtltime <- outer(ltime,stime,"-")
                          
                            #-Find first point where tacsat time is greater or equal to departure time and smaller than arrival time
                            st <- apply(smdtime,1,function(x){which(x>=0)[1]})
                            en <- apply(gtltime,1,function(x){rev(which(x>=0))[1]})
                            
                            #-Make sure that values are within the interval of departure and arrival time
                            subse <- which(is.na(st<=en)==F & (st<=en)==T)
                            st <- st[subse]
                            en <- en[subse]
                            
                            #-Assign Tacsat data with FT_REF from Eflalo2 dataset where they link
                            if(length(st)!=1){
                              idx   <- unlist(mapply(seq,st,en))
                              reps  <- unlist(lapply(mapply(seq,st,en),length))
                              splitTa[[i]]$FT_REF      <- 0
                              splitTa[[i]]$FT_REF[idx] <- rep(tripn[subse],reps)
                            } 
                            if(length(st)==1){
                              splitTa[[i]]$FT_REF <- 0
                              splitTa[[i]]$FT_REF[seq(st,en)] <- rep(tripn[subse],length(seq(st,en)))
                            }
                            if(length(st)==0){
                              splitTa[[i]]$FT_REF <- 0
                            }
                          }
                          Ta$FT_REF <- unlist(lapply(splitTa,function(x){return(x$FT_REF)}))
                              
                      return(Ta)}

