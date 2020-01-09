########################################
### Echelon dendrogram with clusters ###
########################################

e.cluster.dendrogram <- function(echelon.obj,n.sim,cluster.legend.pos,cluster_reg,p_rank,para){

	if(!is.null(echelon.obj$coord) && length(para) == 0){
		x <- echelon.obj$x
		coord <- echelon.obj$coord
		locs <- echelon.obj$locs
		x_pos <- echelon.obj$x_pos
		separates <- echelon.obj$separates
		replayPlot(echelon.obj$current.dendrogram)
	}

	else{
		main <- ifelse(length(para$main)==0, numeric(0), para$main)
		ylab <- ifelse(length(para$ylab)==0, numeric(0), para$ylab)
		yaxes <- ifelse(length(para$yaxes)==0, TRUE, para$yaxes)

		if(length(para$ylim)==0) ylim <- NULL
		else ylim <- para$ylim

		xaxes <- ifelse(length(para$xaxes)==0, FALSE, para$xaxes)

		if(length(para$xdper)==0) xdper <- c(0,1)
		else xdper <- para$xdper

		if(length(para$dmai)==0) dmai <- NULL
		else dmai <- para$dmai

		lwd <- ifelse(length(para$lwd)==0, 1, para$lwd)
		symbols <- ifelse(length(para$symbols)==0, 4, para$symbols)
		cex.symbols <- ifelse(length(para$cex.symbols)==0, 0.5, para$cex.symbols)
		col2 <- ifelse(length(para$col.symbols)==0, 4, para$col.symbols)
		ens <- ifelse(length(para$ens)==0, TRUE, para$ens)
		adj.ens <- ifelse(length(para$adj.ens)==0, 1, para$adj.ens)
		cex.ens <- ifelse(length(para$cex.ens)==0, 0.8, para$cex.ens)
		col3 <- ifelse(length(para$col.ens)==0, 1, para$col.ens)
		col1 <- ifelse(length(para$col)==0, 1, para$col)

		peaks <- echelon.obj$peaks
		locs <- echelon.obj$locs
		x <- echelon.obj$x
		separates <- echelon.obj$separates
		c_separates <- echelon.obj$c_separates
		number <- echelon.obj$number
		parents <- echelon.obj$parents
		pare_locs <- echelon.obj$pare_locs
		progeny <- echelon.obj$progeny

		e.dendrogram.axis(main,ylab,yaxes,ylim,xaxes,xdper,dmai,peaks,x)
		temp <- e.dendrogram(peaks,locs,x,separates,c_separates,number,parents,
													pare_locs,progeny,symbols,col.symbols=col2,cex.symbols,
													lwd,col=col1,ens,adj.ens,col.ens=col3,cex.ens,limb=FALSE)
		coord <- temp$coordinate
		x_pos <- temp$x_pos
	}


	if(nrow(cluster_reg) <= 5) l <- nrow(cluster_reg)
	else l <- 5
	cluster.col <- c(2:6)
	leg.txt <- NULL
	leg.col <- NULL
	for(j in 1:l){
		k <- NULL
		k2 <- NULL
		for(i in 1:length(separates)) k <- c(k,rep(x_pos[i],times=separates[i]))
			temp <- cluster_reg[j,][!is.na(cluster_reg[j,])]
		for(i in 1:length(temp))  k2 <- c(k2,which(temp[i] == locs))
			points(k[k2],x[locs][k2],pch=10,cex=1.5,col=cluster.col[j])
		if(!is.null(p_rank)) leg.txt <- c(leg.txt,paste(j,"- p-value:",p_rank[j]/(n.sim+1)))
		else leg.txt <- c(leg.txt,paste(j,"- p-value:---"))
		leg.col <- c(leg.col,cluster.col[j])
	}
	legend(cluster.legend.pos,legend=leg.txt,text.col=leg.col)

	list(coord = coord,regions.value=x[locs])
}
