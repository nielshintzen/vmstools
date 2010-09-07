ICESrectangle <- function(dF){
                    rectChar1n2 <- as.integer( 2 * (dF[,"SI_LATI"]-35.5) )

                    rectChar3 <- ifelse(dF[,"SI_LONG"]<=-40,"A"
                                ,ifelse(dF[,"SI_LONG"]<=-30,"B"
                                ,ifelse(dF[,"SI_LONG"]<=-20,"C"
                                ,ifelse(dF[,"SI_LONG"]<=-10,"D"
                                ,ifelse(dF[,"SI_LONG"]<= 0 ,"E"
                                ,ifelse(dF[,"SI_LONG"]<=10,"F"
                                ,ifelse(dF[,"SI_LONG"]<=20,"G"
                                ,ifelse(dF[,"SI_LONG"]<=30,"H"
                                ,"I"))))))))

                    rectChar4 <- as.integer( dF[,"SI_LONG"] %% 10)

                    rectID <- paste(rectChar1n2,rectChar3,rectChar4,sep='')
                return(rectID)}

