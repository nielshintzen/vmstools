makeStackINLA <- function(form,data,projMats,spatFields,idxs,distris,splitVars=NULL,tag="Fit"){
              #- Get types of variables
              expl      <- unlist(strsplit(ac(form)[2],"+",fixed=TRUE))
              vars      <- strsplit(ac(form)[3],"+",fixed=TRUE)[[1]]; vars <- gsub("\n","",vars,fixed=T); vars <- gsub(" ","",vars,fixed=T)
              spats     <- grep("spde",vars)
              rands     <- grep("iid",vars); randsVal <- gsub(")","",gsub("(","",gsub(",model=\"iid\")","",gsub("f(","",gsub(" ","",vars[rands],fixed=T),fixed=T),fixed=T),fixed=T),fixed=T)
              interc    <- grep("Intercept",vars)
              fixs      <- (2:length(vars))[which(!(2:length(vars)) %in% c(spats,rands,interc))]

              #- Add intercepts
              idxsmat   <- matrix(NA,nrow=length(idxs),ncol=length(idxs))
              diag(idxsmat) <- 1
              lens      <- unlist(lapply(idxs,length))
              intercepts<- numeric()
              for(i in 1:ncol(idxsmat))
                intercepts <- cbind(intercepts,unlist(mapply(rep,idxsmat[,i],lens)))
              colnames(intercepts) <- paste0("Intercept",1:length(idxs))

              #- Split variables if necessary
              if(!is.null(splitVars)){
                splitVarsRands <- which(splitVars %in% unique(gsub("2","",gsub("1","",randsVal))))
                splitVarsVal   <- splitVars[splitVarsRands]
                if(length(splitVarsRands)>0)
                  splitVars <- splitVars[-splitVarsRands]
                Xsplit  <- data.frame(rep(0,nrow(data)))
                for(i in 1:length(splitVars))
                  Xsplit <- cbind(Xsplit,as.data.frame(t(mapply(rep,data[,splitVars[i]],times=length(idxs)))))
                Xsplit  <- Xsplit[,-1]
                colnames(Xsplit) <- c(t(outer(splitVars,1:length(idxs),paste0)))
                idxcol  <- anf(substr(colnames(Xsplit),nchar(colnames(Xsplit)),nchar(colnames(Xsplit))))
                for(i in 1:length(idxs))
                  Xsplit[-idxs[[i]],which(idxcol==i)] <- NA
                varsSplit <- vars[fixs][-as.vector(sapply(splitVars,grep,gsub(" ","",vars[fixs])))]
                X         <- cbind(intercepts,Xsplit,as.data.frame(data[,gsub(" ","",varsSplit)]))
              }

              #- Combine the data
              if(is.null(splitVars))
                X       <- cbind(intercepts,as.data.frame(data[,gsub(" ","",vars[fixs])]))
              Y         <- matrix(NA,nrow=nrow(data),ncol=length(idxs))
              colnames(X)[2:(2+length(fixs)-1)] <- gsub(" ","",vars[fixs])
              for(i in 1:ncol(idxsmat)){
                if(distris[i]=="binomial"){
                  Y[idxs[[i]][which(data[,expl]>0)], i] <- 1
                  Y[idxs[[i]][which(data[,expl]==0)],i] <- 0
                }
                if(distris[i]%in%c("gaussian","gamma","beta"))
                  Y[idxs[[i]],i]  <- data[idxs[[i]],expl]
              }
              Y <- as.matrix(Y)
              colnames(Y) <- NULL

              #- Prep the stack fit
              datPrep   <- list(Y)
              names(datPrep)  <- eval(expl)
              obsPrep   <- c(projMats,as.list(rep(1,length(rands)+1)))
              if(length(randsVal)>0){
                if(length(randsVal)>1){
                  if(!is.null(splitVars)){
                    if(length(grep(splitVarsVal,randsVal))>0){
                      Rs  <- as.data.frame(t(mapply(rep,data[,splitVarsVal],times=length(idxs))),stringsAsFactors=F)
                      for(i in 1:length(idxs)){
                        Rs[-idxs[[i]],i] <- NA
                        Rs[,i] <- an(af(Rs[,i]))
                        colnames(Rs) <- randsVal
                        rownames(Rs) <- 1:nrow(Rs)
                      }
                      Rs  <- as.list(Rs)
                    }
                  } else {
                      Rs    <- as.list(data[,randsVal])
                  }
                } else {
                  Rs      <- list(data[,randsVal])
                }
                effPrep   <- c(spatFields,list(X),Rs)
              } else {
                effPrep   <- c(spatFields,list(X))
              }
              names(effPrep) <- c(unlist(lapply(spatFields,function(x){names(x)[1]})),"X",randsVal)


              #- Get the stack
              StackFit  <- inla.stack(
                              tag  = tag,
                              data = datPrep,
                    	        A = obsPrep,
                    	        effects = effPrep)
            return(StackFit)}
