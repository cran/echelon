#-------------------------------------------------------------------------------
#  Clsuter map
#-------------------------------------------------------------------------------
e.cluster.map <- 
  function(x, c_separates, locs, coo, rin, p_rank, cluster_reg, n.sim, cluster.type
  ){

  dev.new()
  par(mai=c(0.1, 0.1, 0.5, 0.1))

  if(cluster.type == "high") plot(coo, ylab = "", xlab = "", axes = F, main = "High rate clusters")
  else if(cluster.type == "low") plot(coo, ylab = "", xlab = "", axes = F, main = "Low rate clusters")

  cluster.col <- c(2:6)
  for(i in seq_along(coo[,1])){
    for(j in seq_along(rin[i,])){
      if(!is.na(rin[i,j])){
        temp <- as.integer(rin[i,j])
        lines(x = c(coo[i, 1], coo[temp, 1]), y = c(coo[i, 2], coo[temp, 2]), lwd = 0.5, lty = 3)
        for(k in 1:nrow(cluster_reg)){
          if(length(which(!is.na(match(c(i, temp), cluster_reg[k,])))) == 2) 
            lines(x = c(coo[i, 1], coo[temp, 1]), y = c(coo[i, 2], coo[temp, 2]), lwd = 2, col = cluster.col[k])
        }
      }
      else break
    }
  }

  if(nrow(cluster_reg) <= 5) l <- nrow(cluster_reg)
  else l <- 5
  leg.txt <- NULL
  leg.col <- NULL
  for(i in 1:l){
    temp <- cluster_reg[i,][!is.na(cluster_reg[i,])]
    points(coo[temp, 1], coo[temp, 2], pch = 21, bg = "white", col = cluster.col[i], cex = 2)
    text(coo[temp, 1], coo[temp, 2], labels = i, col = cluster.col[i], cex = 1.5)
    if(!is.null(p_rank)) leg.txt <- c(leg.txt,paste(i, "- p-value:", p_rank[i]/(n.sim+1)))
    else leg.txt <- c(leg.txt,paste(i, "- p-value:---"))
    leg.col <- c(leg.col, cluster.col[i])
  }
  legend("bottomleft", legend = leg.txt, text.col = leg.col)

  dev.set(dev.prev())
}
