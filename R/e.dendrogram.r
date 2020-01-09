##########################
### Echelon dendrogram ###
##########################

e.dendrogram<-function(peaks,locs,x,separates,c_separates,number,parents,pare_locs,progeny,
		symbols,col.symbols,cex.symbols,
		lwd,col,ens,adj.ens,col.ens,cex.ens,limb){
	
	x_pos<-rep(0,times=length(peaks))
	x_pos[1]<-1
	if(length(peaks)!=1){
		sumi<-1
		k<-1
		ima<-1
		while(length(sumi)!=length(peaks)){
			imabura<-which(parents==parents[ima])
			if(all(is.element(imabura,sumi))){
				ima<-parents[ima]
				sumi<-c(sumi,ima)
			}
			else{
				for(i in 1:length(imabura)){
					if(!any(imabura[i]==sumi)){
						ima<-imabura[i]
						break
					}
				}
				if(peaks[ima]==1){
					k<-k+1
					x_pos[ima]<-k
					sumi<-c(sumi,ima)
				}
				else{
					ima<-which(parents==ima)[1]
				}
			}
		}
		x_pos<-x_pos-1
	}
	else x_pos<--0.15

	for(i in 1:length(x_pos)){
		if(parents[i]==0&&progeny[i]==0){
			lines(c(x_pos[i],x_pos[i]),c(x[locs[(c_separates[i]+1)]],min(x)),lwd=lwd,col=col,lty=ifelse(limb,limb[i],1))
		}

		else if(x_pos[i]!=-1){
			lines(c(x_pos[i],x_pos[i]),c(x[locs[(c_separates[i]+1)]],x[pare_locs[i]]),lwd=lwd,col=col,lty=ifelse(limb,limb[i],1))
		}
		
		else{
			if(parents[i]!=0) a<-c(x[locs[(c_separates[i]+1)]],x[pare_locs[i]])
			else a<-c(x[locs[(c_separates[i]+1)]],min(x))
			child<-which(parents==i)
			k<-c(min(x_pos[child]),max(x_pos[child]))
			lines(c(median(k),median(k)),a,lwd=lwd,col=col,lty=ifelse(limb,limb[i],1))
			lines(k,c(a[1],a[1]),lwd=lwd,col=col,lty=1)
			x_pos[i]<-median(k)[1]
		}
	}
	

	if(!is.null(symbols)){  #symbols
		temp<-NULL
		for(i in 1:length(separates)) 
			temp<-c(temp,rep(x_pos[i],times=separates[i]))
		coordinate<-array(c(temp,x[locs]),dim=c(length(temp),2))
		points(temp,x[locs],pch=symbols,col=col.symbols,cex=cex.symbols)
	}


	if(ens){  #Echelon numbers
		temp<-NULL
		for(i in 1:length(x_pos)) temp<-c(temp,x[locs[(c_separates[i]+1)]])
		text(x_pos,temp,labels=paste("EN",number,sep=""),adj=adj.ens,col=col.ens,cex=cex.ens)
	}


	if(is.numeric(limb)){  ##profile legend
		temp<-c("1st limb","2nd limb","3rd limb","4th limb","5th limb","6th limb-")
		if(max(limb)<7) legend("bottomright",temp[1:max(limb)],lty=c(1:max(limb)))
		else legend("bottomright",temp[1:6],lty=c(1:6))
	}

	list(coordinate=coordinate,x_pos=x_pos)
}
