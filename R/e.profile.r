#####################
## Echelon profile ##
#####################

e.profile<-function(peaks,parents,separates){
	if(length(peaks)==1){
		di<-1
		sc<-separates
		bu<-0
		st<-1
		eachlimb<-1
	}
	else{
		di<-NULL
		sc<-NULL
		bu<-NULL
		st<-NULL
		eachlimb<-NULL
		genlimb<-NULL
		k2<-0
		mezasu<-0
		limb<-NULL
		while(!sum(di)==length(peaks)){
			for(j in mezasu){
				for(i in which(peaks!=0)){
					if(parents[i]==j) limb<-c(limb,i)
				}
				for(l in which(parents==j)[which(!is.element(which(parents==j),mezasu))][which(peaks[which(parents==j)[which(!is.element(which(parents==j),mezasu))]]==0)]){
					for(i in which(peaks!=0)[which(!is.element(which(peaks!=0),limb))]){
						kalimb<-l
						ii<-i
						while(parents[ii]!=l){
							if(is.element(parents[ii],mezasu)){
								kalimb<-NULL
								break
							}
							ii<-parents[ii]
							if(ii==0){
								kalimb<-NULL
								break
							}
							kalimb<-c(kalimb,ii)
						}
						if(length(kalimb)>k2){
							k2<-length(kalimb)
							genlimb<-kalimb
						}
						else if(length(kalimb)==k2){
							genlimb<-rbind(genlimb,kalimb)
						}
					}
					if(k2!=0){
						genlimb<-unique(genlimb,MARGIN=1)
						limb<-c(limb,nlimb(genlimb,k2))
					}
					k2<-0
					genlimb<-NULL
			}}
			di<-c(di,length(limb))
			sc<-c(sc,sum(separates[limb]))
			eachlimb<-c(eachlimb,limb)
			k3<-NULL
			for(i in 1:length(limb)){
				k3<-c(k3,which(parents==limb[i]))
			}
			bu<-c(bu,length(k3[which(!is.element(k3,limb))]))
			st<-c(st,length(limb[which(peaks[limb]==1)]))
			mezasu<-limb[which(peaks[limb]==0)]
			limb<-NULL
	}}
	outmeasure<-rbind(cumsum(di)/length(peaks),cumsum(sc)/sum(separates),bu/sum(peaks),cumsum(st)/sum(peaks))
	rownames(outmeasure)<-c("Divergence","Scope","Bunching","Stacking")
	k2<-NULL
	for(i in 1:length(di)){
		k2<-c(k2,paste("cycle ",i,sep=""))
	}
	colnames(outmeasure)<-k2
	
	if(!length(peaks)==1){
		k2<-NULL
		for(i in 1:length(di)){
			k2<-c(k2,rep(i,times=di[i]))
		}
		eachlimb<-cbind(eachlimb,k2)
		eachlimb<-eachlimb[sort.list(eachlimb[,1]),][,2]
	}
	
	dev.new()
	plot(c(1:length(di)),outmeasure[1,],col=2,type="b",main="Echelon profiles",xlab="cycle",ylab="",ylim=c(0,1),xaxt="n")
	axis(side=1,at=c(1:length(di)),labels=c(1:length(di)))
	lines(c(1:length(di)),outmeasure[2,],col=3,type="b",pch=2)
	lines(c(1:length(di)),outmeasure[3,],col=4,type="b",pch=3)
	lines(c(1:length(di)),outmeasure[4,],col=6,type="b",pch=0)
	legend("right",legend=c("Divergence","Scope","Bunching","Stacking"), col=c(2:4,6),pch=c(1:3,0),cex=0.9,pt.cex=0.9)
	abline(h=c(0.2,0.4,0.6,0.8),lty=3)
	dev.set(dev.prev())
	list(outmeasure=outmeasure,di=di,sc=sc,bu=bu,st=st,eachlimb=eachlimb)
}

nlimb<-function(x,k2){
  if(k2==1) outlimb<-as.vector(x)
	else{
		outlimb<-NULL
  	k<-0
  	nuki<-NULL
  	for(i in 1:nrow(x)){
    	for(j in c(1:nrow(x))[-i]){
      	if(any(x[i,]==x[j,])) k<-1
    	}
    	if(k==0){
      	outlimb<-c(outlimb,x[i,])
      	nuki<-c(nuki,i)
    	}
    	k<-0
  	}
  	if(!is.null(nuki)) x<-x[-nuki,]
  	for(i in k2:1){
    	if(any(duplicated(x[,i]))) outlimb<-c(outlimb,unique(x[,i][which(duplicated(x[,i]))]))
  	}
  }
  return(outlimb)
}
