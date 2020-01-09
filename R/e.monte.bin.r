########################################
### Monte Carlo Estimation (echebin) ###
########################################

e.monte.bin <- function(rin,cas,pop,n.sim,K,cluster.type){

	multi <- rmultinom(n.sim,round(sum(cas)),prob=pop)
	lambda <- NULL

	for(i in 1:n.sim){
		ctl <- pop - multi[,i]

		if(cluster.type == "high"){
			temp <- e.main(multi[,i]/ctl,rin,length(cas))
			x <- multi[,i]/ctl
			locs <- temp$locs
			peaks <- temp$peaks
			separates <- temp$separates
			c_separates <- c(0,cumsum(separates))
			parents <- temp$parents

			if(K < 1) temp <- e.scan.pop(x=x,locs=locs,peaks=peaks,c_separates=c_separates,parents=parents,K=K,pop=pop)
			else temp <- e.scan(x=x,locs=locs,peaks=peaks,c_separates=c_separates,parents=parents,K=K)
			reg_data <- temp$reg_data

			if(!is.null(reg_data)){
				popg <- sum(pop)
				casg <- sum(multi[,i])
				ctlg <- sum(ctl)

				casz <- apply(array(multi[,i][reg_data],c(nrow(reg_data),ncol(reg_data))),1,sum,na.rm=TRUE)
				ctlz <- apply(array(ctl[reg_data],c(nrow(reg_data),ncol(reg_data))),1,sum,na.rm=TRUE)
				popz <- apply(array(pop[reg_data],c(nrow(reg_data),ncol(reg_data))),1,sum,na.rm=TRUE)

				temp1 <- casz*log(casz/popz)
				temp2 <- ctlz*log(ctlz/popz)
				temp3 <- (casg-casz)*log((casg-casz)/(popg-popz))
				temp4 <- (ctlg-ctlz)*log((ctlg-ctlz)/(popg-popz))
				temp5 <- -casg*log(casg/popg)
				temp6 <- -ctlg*log(ctlg/popg)
				temp <- cbind(temp1,temp2,temp3,temp4,temp5,temp6)
				log.lambda <- apply(temp,1,sum,na.rm=TRUE)
				log.lambda[which(casz/popz < casg/popg)] <- 0

				lambda <- c(lambda,max(log.lambda))
			}
			else lambda <- c(lambda,0)
		}

		if(cluster.type == "low"){
			temp <- e.main(-multi[,i]/ctl,rin,length(cas))
			x <- -multi[,i]/ctl
			locs <- temp$locs
			peaks <- temp$peaks
			separates <- temp$separates
			c_separates <- c(0,cumsum(separates))
			parents <- temp$parents

			if(K < 1) temp <- e.scan.pop(x=x,locs=locs,peaks=peaks,c_separates=c_separates,parents=parents,K=K,pop=pop)
			else temp <- e.scan(x=x,locs=locs,peaks=peaks,c_separates=c_separates,parents=parents,K=K)
			reg_data <- temp$reg_data

			if(!is.null(reg_data)){
				popg <- sum(pop)
				casg <- sum(multi[,i])
				ctlg <- sum(ctl)

				casz <- apply(array(multi[,i][reg_data],c(nrow(reg_data),ncol(reg_data))),1,sum,na.rm=TRUE)
				ctlz <- apply(array(ctl[reg_data],c(nrow(reg_data),ncol(reg_data))),1,sum,na.rm=TRUE)
				popz <- apply(array(pop[reg_data],c(nrow(reg_data),ncol(reg_data))),1,sum,na.rm=TRUE)

				temp1 <- casz*log(casz/popz)
				temp2 <- ctlz*log(ctlz/popz)
				temp3 <- (casg-casz)*log((casg-casz)/(popg-popz))
				temp4 <- (ctlg-ctlz)*log((ctlg-ctlz)/(popg-popz))
				temp5 <- -casg*log(casg/popg)
				temp6 <- -ctlg*log(ctlg/popg)
				temp <- cbind(temp1,temp2,temp3,temp4,temp5,temp6)
				log.lambda <- apply(temp,1,sum,na.rm=TRUE)
				log.lambda[which(casz/popz > casg/popg)] <- 0

				lambda <- c(lambda,max(log.lambda))
			}
			else lambda <- c(lambda,0)
		}
		#cat(paste(i,", "),sep="")
		#if(i%%10 == 0) cat("\n")
	}
	cat("\n\n")
	return(lambda)
}
