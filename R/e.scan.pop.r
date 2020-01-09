######################################
### Echelon scan (population base) ###
######################################

e.scan.pop <- function(x,locs,peaks,c_separates,parents,K,pop){

	e_num <- 1
	Kout <- 0
	ng <- sum(pop)

	while(e_num <= length(peaks)){
		if(peaks[e_num] == 1){
			areas <- locs[(c_separates[e_num]+1):c_separates[e_num+1]]

			if(sum(pop[areas]) > ng*K){
				i<-1
				temp <- NULL
				while(sum(pop[temp]) <= ng*K){
					t_check <- x[areas[i]]
					t_check2 <- which(x[areas] == t_check)
					if(sum(pop[temp])+sum(pop[areas[t_check2]]) <= ng*K){
						temp <- c(temp,areas[t_check2])
						i <- length(temp)+1
					}
					else{
					  if(is.null(temp)) Kout <-1
						break
					}
				}
				if(!is.null(temp)) areas <- temp
			}

			if(Kout == 0){
				temp <- t(array(areas,c(length(areas),length(areas))))
				temp[upper.tri(temp)] <- NA

				if(any(duplicated(x[areas]))){
					del_tie <- which(duplicated(x[areas]))-1
					temp <- temp[-del_tie,,drop=FALSE]
				}

				if(e_num == 1){
					reg_data <- temp
					reg_data_col <- ncol(reg_data)
				}
				else if(ncol(temp) == reg_data_col) reg_data <- rbind(reg_data,temp)
				else if(ncol(temp) > reg_data_col){
					reg_data <- cbind(reg_data,array(NA,dim=c(nrow(reg_data),ncol(temp)-reg_data_col)))
					reg_data <- rbind(reg_data,temp)
					reg_data_col <- ncol(reg_data)
				}
				else{
					temp <- cbind(temp,array(NA,dim=c(nrow(temp),reg_data_col-ncol(temp))))
					reg_data <- rbind(reg_data,temp)
				}
			}
			else Kout <- 0
		}

		else{
			areas <- locs[(c_separates[e_num]+1):c_separates[e_num+1]]
			now_eches <- e_num
			temp <- NULL

			while(any(peaks[now_eches] == 0)){
				child_eches <- which(is.element(parents,now_eches))
				for(i in 1:length(child_eches)) temp <- c(temp,locs[(c_separates[child_eches[i]]+1):c_separates[child_eches[i]+1]])

				e_top <- x[locs[(c_separates[e_num]+1):c_separates[e_num+1]]]
				if(sum(pop[temp])+sum(pop[areas[which(e_top[1]==e_top)]]) > ng*K){
					Kout <- 1
					break
				}

				now_eches <- child_eches
			}

			if(Kout == 0){
				if(sum(pop[temp])+sum(pop[areas]) > ng*K){
					i <- 1
					temp2 <- NULL
					while(sum(pop[temp])+sum(pop[temp2]) <= ng*K){
						t_check <- x[areas[i]]
						t_check2 <- which(x[areas] == t_check)
						if(sum(pop[temp])+sum(pop[temp2])+sum(pop[areas[t_check2]]) <= ng*K){
							temp2 <- c(temp2,areas[t_check2])
							i <- length(temp2)+1
						}
						else break
					}
					areas <- temp2
				}

				temp2 <- t(array(temp,c(length(temp),length(areas))))
				temp <- t(array(areas,c(length(areas),length(areas))))
				temp[upper.tri(temp)] <- NA
				temp <- cbind(temp2,temp)

				if(any(duplicated(x[areas]))){
					del_tie <- which(duplicated(x[areas]))-1
					temp <- temp[-del_tie,,drop=FALSE]
				}

				if(ncol(temp) == reg_data_col) reg_data <- rbind(reg_data,temp)
				else if(ncol(temp) > reg_data_col){
					reg_data <- cbind(reg_data,array(NA,dim=c(nrow(reg_data),ncol(temp)-reg_data_col)))
					reg_data <- rbind(reg_data,temp)
					reg_data_col <- ncol(reg_data)
				}
				else{
					temp <- cbind(temp,array(NA,dim=c(nrow(temp),reg_data_col-ncol(temp))))
					reg_data <- rbind(reg_data,temp)
				}
			}
			else Kout <- 0
		}
		e_num <- e_num+1
	}

	list(reg_data=reg_data)
}
