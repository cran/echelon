###########################################
### Echelon scan based on Binomial model ###
###########################################

echebin<-function(echelon.obj,cas,ctl,K=length(cas)/2,n.sim=99,
									cluster.type="high",cluster.legend.pos="bottomleft",
									dendrogram=TRUE,cluster.info=FALSE,coo=NULL,...){


##############################
### Check of echelon class ###
##############################

	if(class(echelon.obj) != "echelon")
		stop(paste("The class 'echelon' is incorrect\n\n"))

	x <- echelon.obj$x
	rin <- echelon.obj$rin
	locs <- echelon.obj$locs
	peaks <- echelon.obj$peaks
	c_separates <- echelon.obj$c_separates
	parents <- echelon.obj$parents
	pop <- cas + ctl

####################
### Echelon scan ###
####################

	if(is.null(cas)||is.null(ctl)) stop("The binomial model requires two arguments 'cas' and 'ctl'\n")
	if(K <= 0) K <- floor(length(x)/2)
	else if(K < 1) temp <- e.scan.pop(x,locs,peaks,c_separates,parents,K,pop)
	else if(K > length(x)) stop("Please check the argument 'K'. It is 'K' > length (x).\n")
	else{
		K <- floor(K)
		temp <- e.scan(x,locs,peaks,c_separates,parents,K)
	}
	reg_data <- temp$reg_data
	if(is.null(reg_data)) stop("There are no clusters! Try again changing the argument 'x' or 'K'\n")


#############################
### Statistic calculation ###
#############################

	popg <- sum(pop)
	casg <- sum(cas)
	ctlg <- sum(ctl)

	casz <- apply(array(cas[reg_data],c(nrow(reg_data),ncol(reg_data))),1,sum,na.rm=TRUE)
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

	if(cluster.type == "high") log.lambda[which(casz/popz < casg/popg)] <- 0
	if(cluster.type == "low") log.lambda[which(casz/popz > casg/popg)] <- 0


########################
### Cluster decision ###
########################

	temp <- e.cluster.decision(reg_data,log.lambda)
	cluster_reg <- temp$cluster_reg
	cluster_log.lambda <- temp$cluster_log.lambda


##############################
### Monte Carlo Estimation ###
##############################

	p_rank <- NULL
	if(n.sim > 0){
		n.sim <- floor(n.sim)
		cat(paste("Starting",n.sim,"Monte Carlo replications...\n"),sep="")
		sim_lambda <- e.monte.bin(rin,cas,pop,n.sim,K,cluster.type)
		for(i in 1:nrow(cluster_reg)) p_rank <- c(p_rank,floor(rank(c(cluster_log.lambda[i],sim_lambda)*-1)[1]))
	}
	else n_sim <- 0


##########################
### Echelon dendrogram ###
##########################

	if(dendrogram){
		temp <- e.cluster.dendrogram(echelon.obj,n.sim,cluster.legend.pos,cluster_reg,p_rank,para=list(...))
		coord <- temp$coord
	}
	else coord <- NULL


###################
### Cluster map ###
###################

	if(!is.null(coo)){
		if(nrow(coo) != length(x)) stop("length(x) and nrow(coo) must have same size\n\n")
		e.cluster.map(x,c_separates,locs,coo,rin,p_rank,cluster_reg,n.sim,cluster.type)
	}


###############
### Out put ###
###############

	if(cluster.info){
		cat("------------- CLUSTERS DETECTED -------------\n")
		cat(paste("Number of locations ......: ",length(x)," region\n",sep=""))

		if(K >= 1) cat(paste("Limit length of cluster ..: ",K," regions\n",sep=""))
		else cat(paste("Limit length of cluster ..: ",K*100," percent of population\n",sep=""))

		cat(paste("Total cases ..............: ",sum(cas),"\n",sep=""))
		if(!is.null(pop)) cat(paste("Total population .........: ",sum(pop),"\n",sep=""))
		if(cluster.type == "high") cat(paste("Scan for Area with .......: High Rates\n"))
		if(cluster.type == "low") cat(paste("Scan for Area with .......: Low Rates\n"))
		cat(paste("Number of Replications ...: ",n.sim,"\n",sep=""))
		cat("Model ....................: Binomial\n\n")

		cat("---------------------------------------------\n")
	}

	cat(paste("MOST LIKELY CLUSTER -- ",length(cluster_reg[1,][!is.na(cluster_reg[1,])])," regions\n Cluster regions included : ",sep=""))
	cat(echelon.obj$reg_name[cluster_reg[1,][!is.na(cluster_reg[1,])]],sep=", ")
	MLC <- list(regionsID=cluster_reg[1,][!is.na(cluster_reg[1,])])

	pop_inZ <- sum(pop[cluster_reg[1,]],na.rm=TRUE)
	cat(paste("\n Population ..............: ",pop_inZ,sep=""))
	MLC <- c(MLC,list(pop_inZ=pop_inZ))

	cas_inZ <- sum(cas[cluster_reg[1,]],na.rm=TRUE)
	cas_outZ <- casg-cas_inZ
	ex_inZ <- sum(pop[cluster_reg[1,]]*casg/popg,na.rm=TRUE)
	ex_outZ <- casg-ex_inZ

	cat(paste("\n Number of cases .........: ",cas_inZ,sep=""))
	cat(paste("\n Expected cases ..........: ",round(ex_inZ,digits=4),sep=""))
	cat(paste("\n Observed / expected .....: ",round(cas_inZ/ex_inZ,digits=4),sep=""))
	cat(paste("\n Relative risk ...........: ",round((cas_inZ/ex_inZ)/(cas_outZ/ex_outZ),digits=4),sep=""))
	cat(paste("\n Log likelihood ratio ....: ",round(cluster_log.lambda[1],digits=4),"",sep=""))
	MLC <- c(MLC,list(cas_inZ=cas_inZ,ex_inZ=ex_inZ,LLR=cluster_log.lambda[1]))

	if(n.sim != 0){
		cat(paste("\n Monte Carlo rank ........: ",p_rank[1],"/",n.sim+1,"",sep=""))
		cat(paste("\n P-value .................: ",p_rank[1]/(n.sim+1),"",sep=""))
		cat("\n\n")
		MLC <- c(MLC,list(p=p_rank[1]/(n.sim+1)))
	}
	clusters <- MLC

	if(cluster.info) cat("----------------------------------------------\n")

	if(nrow(cluster_reg) != 1){
	  if(cluster.info) cat("SECONDARY CLUSTERS\n")

		if(nrow(cluster_reg) > 5) len2C <- 5
		else len2C <- nrow(cluster_reg)

		for(i in 2:len2C){
			secondC <- NULL
			if(cluster.info){
			  cat(paste(i," -- ",length(cluster_reg[i,][!is.na(cluster_reg[i,])])," regions\n Cluster regions included : ",sep=""))
				cat(echelon.obj$reg_name[cluster_reg[i,][!is.na(cluster_reg[i,])]],sep=", ")
			}
			secondC <- list(regionsID=cluster_reg[i,][!is.na(cluster_reg[i,])])

			if(!is.null(pop)){
				pop_inZ <- sum(pop[cluster_reg[i,]],na.rm=TRUE)
				if(cluster.info) cat(paste("\n Population ..............: ",pop_inZ,sep=""))
				secondC <- c(secondC,list(pop_inZ=pop_inZ))
			}

			cas_inZ <- sum(cas[cluster_reg[i,]],na.rm=TRUE)
			cas_outZ <- casg-cas_inZ
			ex_inZ <- sum(pop[cluster_reg[i,]]*casg/popg,na.rm=TRUE)
			ex_outZ <- casg-ex_inZ

			if(cluster.info){
			  cat(paste("\n Number of cases .........: ",cas_inZ,sep=""))
				cat(paste("\n Expected cases ..........: ",round(ex_inZ,digits=4),sep=""))
				cat(paste("\n Observed / expected .....: ",round(cas_inZ/ex_inZ,digits=4),sep=""))
				cat(paste("\n Relative risk ...........: ",round((cas_inZ/ex_inZ)/(cas_outZ/ex_outZ),digits=4),sep=""))
				cat(paste("\n Log likelihood ratio ....: ",round(cluster_log.lambda[i],digits=4),"",sep=""))
			}
			secondC <- c(secondC,list(cas_inZ=cas_inZ,ex_inZ=ex_inZ,LLR=cluster_log.lambda[i]))

			if(n.sim > 0){
			  if(cluster.info){
			    cat(paste("\n Monte Carlo rank ........: ",p_rank[i],"/",n.sim+1,"",sep=""))
					cat(paste("\n P-value .................: ",p_rank[i]/(n.sim+1),"",sep=""))
				}
			  secondC <- c(secondC,list(p=p_rank[i]/(n.sim+1)))
			}
			if(i == 2) clusters <- list(clusters,secondC)
			else clusters[[i]] <- secondC
			if(cluster.info) cat("\n\n")
		}
		if(cluster.info) cat("----------------------------------------------\n")

		if(nrow(cluster_reg) >5){
		  if(cluster.info) cat("Display only the top 5 clusters. See object 'clusters' for more details\n\n")
			for(i in 6:nrow(cluster_reg)){
				secondC <- NULL
				secondC <- list(regionsID=cluster_reg[i,][!is.na(cluster_reg[i,])])
				if(!is.null(pop)){
					pop_inZ <- sum(pop[cluster_reg[i,]],na.rm=TRUE)
					secondC <- c(secondC,list(pop_inZ=pop_inZ))
				}
				cas_inZ <- sum(cas[cluster_reg[i,]],na.rm=TRUE)
				ex_inZ <- sum(pop[cluster_reg[i,]]*casg/popg,na.rm=TRUE)
				secondC <- c(secondC,list(cas_inZ=cas_inZ,ex_inZ=ex_inZ,LLR=cluster_log.lambda[i]))
				if(n.sim != 0){
					secondC <- c(secondC,list(p=p_rank[i]/(n.sim+1)))
				}
				clusters[[i]] <- secondC
			}
		}
	}
	else clusters <- list(clusters, "not detected")

	result <- c(clusters=list(clusters),list(scanned.regions=reg_data))
	if(n.sim > 0) result <- c(result,list(simulated.LLR=sim_lambda))
	if(dendrogram) result <- c(result,list(coord=coord,regions.value=echelon.obj$regions.value,
															regions.name=echelon.obj$reg_name[locs]))
	invisible(result)
}
