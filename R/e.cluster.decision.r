#-------------------------------------------------------------------------------
# Cluster decision
#-------------------------------------------------------------------------------
e.cluster.decision <- function(reg_data, log.lambda){
	
  log.lambda2 <- log.lambda[order(log.lambda, decreasing = TRUE)]
  reg_data2 <- reg_data[order(log.lambda, decreasing = TRUE), , drop = FALSE]

  k <- NULL
  cluster_reg <- reg_data2[1, , drop = FALSE]
	
  if(nrow(reg_data2) != 1){
    for(i in 2:nrow(reg_data2)){
      if(any(is.element(cluster_reg[which(!is.na(cluster_reg))], reg_data2[i,][which(!is.na(reg_data2[i,]))]))) k <- c(k, i)
      else cluster_reg <- rbind(cluster_reg, reg_data2[i,])
    }
    if(!is.null(k)) cluster_log.lambda <- log.lambda2[-k]
    else cluster_log.lambda <- log.lambda2
  }
  else cluster_log.lambda <- log.lambda2

  if(all((cluster_log.lambda == 0))) stop("No clusters found! Please Try again by modifying the argument 'x', 'K' or 'cluster.type'\n")

  log0 <- which(cluster_log.lambda == 0)
  if(length(log0) != 0){
    cluster_reg <- cluster_reg[-log0, , drop = FALSE]
    cluster_log.lambda <- cluster_log.lambda[-log0]
  }

  list(cluster_reg = cluster_reg, cluster_log.lambda = cluster_log.lambda)
}
