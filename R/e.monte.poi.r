########################################
### Monte Carlo Estimation (echepoi) ###
########################################

e.monte.poi <- function(rin,cas,pop,ex,n.sim,K,cluster.type){

	multi <- rmultinom(n.sim,round(sum(cas)),prob=ex)
	lambda <- NULL

	for(i in 1:n.sim){

		if(cluster.type == "high"){
			temp <- e.main(multi[,i]/ex,rin,length(cas))
			x <- multi[,i]/ex
			locs <- temp$locs
			peaks <- temp$peaks
			separates <- temp$separates
			c_separates <- c(0,cumsum(separates))
			parents <- temp$parents

			if(K < 1) temp <- e.scan.pop(x=x,locs=locs,peaks=peaks,c_separates=c_separates,parents=parents,K=K,pop=pop)
			else temp <- e.scan(x=x,locs=locs,peaks=peaks,c_separates=c_separates,parents=parents,K=K)
			reg_data <- temp$reg_data

			if(!is.null(reg_data)){
				cg <- sum(multi[,i])
				eg <- sum(ex)
				cz <- apply(array(multi[,i][reg_data],c(nrow(reg_data),ncol(reg_data))),1,sum,na.rm=TRUE)
				ez <- apply(array(ex[reg_data],c(nrow(reg_data),ncol(reg_data))),1,sum,na.rm=TRUE)

				temp <- cbind(cz*log(cz/ez),(cg-cz)*log((cg-cz)/(eg-ez)))
				log.lambda <- apply(temp,1,sum,na.rm=TRUE)
				log.lambda[which(cz<ez)] <- 0
			
				lambda <- c(lambda,max(log.lambda))
			}
			else lambda <- c(lambda,0)
		}

		if(cluster.type == "low"){
			temp <- e.main(-multi[,i]/ex,rin,length(cas))
			x <- -multi[,i]/ex
			locs <- temp$locs
			peaks <- temp$peaks
			separates <- temp$separates
			c_separates <- c(0,cumsum(separates))
			parents <- temp$parents

			if(K < 1) temp <- e.scan.pop(x=x,locs=locs,peaks=peaks,c_separates=c_separates,parents=parents,K=K,pop=pop)
			else temp <- e.scan(x=x,locs=locs,peaks=peaks,c_separates=c_separates,parents=parents,K=K)
			reg_data <- temp$reg_data

			if(!is.null(reg_data)){
				cg <- sum(multi[,i])
				eg <- sum(ex)
				cz <- apply(array(multi[,i][reg_data],c(nrow(reg_data),ncol(reg_data))),1,sum,na.rm=TRUE)
				ez <- apply(array(ex[reg_data],c(nrow(reg_data),ncol(reg_data))),1,sum,na.rm=TRUE)

				temp <- cbind(cz*log(cz/ez),(cg-cz)*log((cg-cz)/(eg-ez)))
				log.lambda <- apply(temp,1,sum,na.rm=TRUE)
				log.lambda[which(cz>ez)] <- 0

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