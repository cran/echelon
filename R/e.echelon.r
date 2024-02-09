#######################
### Echelon program ###
#######################

echelon<-function(x,nb,dendrogram=TRUE,name=NULL,
							main=NULL,ylab=NULL,yaxes=TRUE,ylim=NULL,xaxes=FALSE,xdper=c(0,1),dmai=NULL,
							col=1,lwd=1,symbols=4,cex.symbols=1,col.symbols=4,
							ens=TRUE,adj.ens=1,cex.ens=0.8,col.ens=1,
							profiles=FALSE){


	T<-length(x)

	nas<-length(x[is.na(x)])
	x[is.na(x)]<-min(x[!is.na(x)])

	if(is.null(name)) reg_name <- 1:T
	else reg_name <- as.vector(name)


####################################
### Neighbor information setting ###
####################################
  if(any(class(nb) == "nb")){
		if(length(nb) != T) stop("the number of regions in nb and x must have same size\n")
		max_po <- max(sapply(nb,length))
		rin <- array(NA, dim=c(T, max_po))
		for(i in 1:T){
			if(any(nb[[i]] == 0)) next
			else rin[i,] <- c(nb[[i]], rep(NA, length = max_po - sapply(nb,length)[i]))
		}
	}

	else if(is.matrix(nb) && nrow(nb) == ncol(nb) && all(diag(nb) == 0)){
		if(nrow(nb) != T) stop("the number of regions in nb and x must have same size\n")
		max_po <- max(apply(nb, 1, function(y) length(which(y!=0))))
 		rin <- array(NA, dim=c(nrow(nb), max_po))
 		for(i in 1:nrow(rin)){
			temp <- which(nb[i,] != 0)
			rin[i,] <- c(temp, rep(NA, length = max_po - length(temp)))
		}
	}
	else{
		if(nrow(nb) != T) stop("the number of regions in nb and x must have same size\n")
		rin <- as.matrix(nb[,-1])
	}


#######################
### Neighbor checks ###
#######################

	for(i in 1:T){
		if(any(rin[i,][!is.na(rin[i,])]>T)) stop(paste("nb error!\nID",i,"includes a larger number than",T,"\n"))
		if(any(duplicated(rin[i,][!is.na(rin[i,])]))) stop(paste("nb error!\nID",i,"includes same numbers\n"))
		for(j in rin[i,][!is.na(rin[i,])])
			if(!any(i==rin[j,][!is.na(rin[j,])])) stop(paste("nb error!\nThere is an incorrect relationship between ID",i,"and",j,"\n"))
	}


########################
### Echelon analysis ###
########################

	e.result <- e.main(x,rin,T)

	len_en<-length(e.result$peaks)
	locs<-e.result$locs
	peaks<-e.result$peaks
	separates<-e.result$separates
	c_separates<-c(0,cumsum(separates))
	parents<-e.result$parents
	tops<-e.result$tops
	bottoms<-e.result$bottoms
	progeny<-e.result$progeny
	family<-e.result$family
	pare_locs<-e.result$pare_locs

	number<-rep(0,times=length(peaks))
	number[which(peaks==1)]<-c(1:sum(peaks))
	number[number==0]<-c((sum(peaks)+1):length(peaks))

	Echelons<-NULL
	for(i in 1:len_en){
		temp<-locs[(c_separates[which(i==number)]+1):c_separates[which(i==number)+1]]
		a<-list(paste(reg_name[temp],"(",x[temp],")",sep=""))
		names(a)<-paste("EN",i,sep="")
		Echelons<-c(Echelons,a)
	}

	parent2<-parents
	for(i in 1:len_en){
		if(i!=number[i]) parent2[which(parents==i)]<-number[i]
	}

	e_order<-peaks
	for(i in 1:len_en){
		if(e_order[i]!=1){
			o_temp<-which(parents==i)
			if(all(e_order[o_temp]==e_order[o_temp][1])) e_order[i]<-e_order[o_temp][1]+1
			else e_order[i]<-max(e_order[o_temp])
		}
	}
	e_order[is.na(e_order)]<-1

	e_max<-x[tops]

	e_min<-x[bottoms]

	e_length<-NULL
	temp<-parents
	temp[temp==0]<-1
	for(i in 1:len_en){
		if(temp[i]==1) e_length<-c(e_length,x[tops][i]-min(x))
		else e_length<-c(e_length,x[tops][i]-x[tops[temp]][i])
	}

	e_level<-rep(0,times=len_en)
	for(i in 1:len_en){
		l_temp<-parents[i]
		while(l_temp!=0){
			e_level[i]<-e_level[i]+1
			l_temp<-parents[l_temp]
		}
	}

	info<-data.frame(EN=number,Order=e_order,Parent=parent2,Maxval=e_max,Minval=e_min,Length=e_length,Cells=separates,Progeny=progeny,Family=family,Level=e_level)
	info<-info[order(info$EN),]
	info[[1]]<-NULL
	row.names(info)<-c(1:len_en)



########################
### Echelon profiles ###
########################

	if(profiles){
		temp<-e.profile(peaks,parents,separates)
		profiles_list <- temp$outmeasure
		limb <- temp$eachlimb
	}
	else{
		limb <- FALSE
		profiles_list <- NULL
	}

##########################
### Echelon dendrogram ###
##########################

	if(dendrogram){
		e.dendrogram.axis(main,ylab,yaxes,ylim,xaxes,xdper,dmai,peaks,x)
		temp <- e.dendrogram(peaks,locs,x,separates,c_separates,number,parents,
													pare_locs,progeny,symbols,col.symbols,cex.symbols,
													lwd,col,ens,adj.ens,col.ens,cex.ens,limb)
		coord <- temp$coordinate
		x_pos <- temp$x_pos
		current.dendrogram <- recordPlot()
	}
	else{
		coord <- NULL
		x_pos <- NULL
		current.dendrogram <- NULL
	}


###############
### Out put ###
###############
	cat(nrow(info), "echelons are created\n")
	cat("See object 'Table' and 'Echelons' for more details\n\n")
	if(profiles) print(list(Profiles=profiles_list))

	k<-list(Echelons=Echelons,Table=info,coord=coord,
	        regions.value=x[locs], regions.name=reg_name[locs], reg_name=reg_name,
	        current.dendrogram=current.dendrogram,x_pos=x_pos,x=x,rin=rin,
	        locs=locs,peaks=peaks,separates=separates,c_separates=c_separates,
	        number=number,parents=parents,pare_locs=pare_locs,progeny=progeny,
	        profiles=profiles_list)

	class(k)<-"echelon"
	invisible(k)
}

