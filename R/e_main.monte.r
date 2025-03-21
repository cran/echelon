#-------------------------------------------------------------------------------
# Main and scan program for Monte Carlo
#-------------------------------------------------------------------------------
e.main.monte <- function(x, rin, K, Kmin, par1, par2, type){
	
  if(type == 11){
    cas <- x
    ex <- par1
    x <- x/ex
  }
  else if(type == 12){
    cas <- x
    ex <- par1
    x <- -x/ex
  }
  else if(type == 21){
    cas <- x
    ctl <- par2 - cas
    x <- x/ctl
  }
  else if(type == 22){
    cas <- x
    ctl <- par2 - cas
    x <- -x/ctl
  }
  dat <- x


#-------------------------------------------------------------------------------
# Echelon analysis
#-------------------------------------------------------------------------------
  total_eche_locs <- numeric(length(dat))
  include_eche_locs <- NULL
  eche_separates <- NULL
  eche_separates_include <- NULL
  eche_peaks <- NULL
  pare_temp <- NULL

  while(any(!is.na(dat))){
    now_val <- max(dat,na.rm = TRUE)
    now_loc <- which(dat == now_val)[1]
		
    t_dummy <- now_loc

    now_eche_locs <- now_loc
    p_dummy <- 1

    while(!is.na(t_dummy)){
      now_nei <- as.vector(rin[now_loc,])
      now_nei <- now_nei[!is.na(now_nei)]
      if(any(is.element(now_nei, now_loc))) now_nei <- now_nei[-which(is.element(now_nei, now_loc))]
      now_nei <- unique(now_nei)

      if(any(is.element(total_eche_locs, now_nei))){
        temp <- rep(seq_along(eche_separates_include), times = eche_separates_include)
        del_eche_separates_include <- unique(temp[is.element(include_eche_locs, now_nei)])
        now_loc <- c(include_eche_locs[is.element(temp, del_eche_separates_include)], now_loc)

        include_eche_locs <- include_eche_locs[!is.element(include_eche_locs, now_loc)]
        include_eche_locs <- c(include_eche_locs, now_loc)

        if(p_dummy == 1) eche_separates_include <- eche_separates_include[-del_eche_separates_include]
        else eche_separates_include <- eche_separates_include[-c(del_eche_separates_include, length(eche_separates_include))]
        eche_separates_include <- c(eche_separates_include, length(now_loc))

        now_nei <- as.vector(rin[now_loc,])
        now_nei <- now_nei[!is.na(now_nei)]
        if(any(is.element(now_nei, now_loc))) now_nei <- now_nei[-which(is.element(now_nei, now_loc))]
        now_nei <- unique(now_nei)

        p_dummy <- p_dummy + length(del_eche_separates_include)
      }

      if(any(x[now_nei] == x[t_dummy])){
        now_loc <- c(now_loc, now_nei[which(x[now_nei] == x[t_dummy])])
        now_eche_locs <- c(now_eche_locs, now_nei[which(x[now_nei] == x[t_dummy])])
				
        if(p_dummy != 1){
          include_eche_locs <- c(include_eche_locs, now_nei[which(x[now_nei] == x[t_dummy])])
          eche_separates_include[length(eche_separates_include)] <- eche_separates_include[length(eche_separates_include)] + length(now_nei[which(x[now_nei] == x[t_dummy])])
        }
      }
      else t_dummy <- NA
    }

    i_dummy <- now_eche_locs
    s_dummy <- 0

    while(s_dummy == 0){
      now_nei <- as.vector(rin[now_loc,])
      now_nei <- now_nei[!is.na(now_nei)]
      if(any(is.element(now_nei, now_loc))) now_nei <- now_nei[-which(is.element(now_nei, now_loc))]
      now_nei <- unique(now_nei)

      if(length(now_nei) == 0){
        dat[now_eche_locs] <- NA
        total_eche_locs[which(total_eche_locs==0)[1] + seq_along(now_eche_locs) - 1] <- now_eche_locs
        if(p_dummy == 1) eche_peaks <- c(eche_peaks, 1)
        else eche_peaks <- c(eche_peaks, 0)
        eche_separates <- c(eche_separates, length(now_eche_locs))
        now_eche_locs <- NULL
        pare_temp <- c(pare_temp, 0)
        break
      }

      now_nei_max_val <- max(dat[now_nei])
      now_nei_max_locs <- now_nei[which(now_nei_max_val == dat[now_nei])]
      t_dummy <- now_nei_max_val

      while(!is.na(t_dummy)){
        sub_nei <- c(now_nei,as.vector(rin[now_nei_max_locs,]))
        sub_nei <- sub_nei[!is.na(sub_nei)]
        sub_nei <- sub_nei[-which(is.element(sub_nei, c(now_loc, now_nei_max_locs)))]
        sub_nei <- unique(sub_nei)

        if(any(x[sub_nei] == t_dummy)){
          now_nei <- c(now_nei, sub_nei[x[sub_nei] == t_dummy])
          now_nei_max_locs <- c(now_nei_max_locs, sub_nei[x[sub_nei] == t_dummy])
        }
				
        else t_dummy <- NA
      }

      if(length(sub_nei) == 0){
        now_eche_locs <- c(now_eche_locs, now_nei_max_locs)
        dat[now_eche_locs] <- NA
        total_eche_locs[which(total_eche_locs==0)[1] + seq_along(now_eche_locs) - 1] <- now_eche_locs
        if(p_dummy == 1){
          eche_peaks <- c(eche_peaks, 1)
          include_eche_locs <- c(include_eche_locs, now_eche_locs)
          eche_separates_include <- c(eche_separates_include, length(now_eche_locs))
        }
        else{
          eche_peaks <- c(eche_peaks, 0)
          include_eche_locs <- c(include_eche_locs, now_eche_locs[-c(seq_along(i_dummy))])
          eche_separates_include[length(eche_separates_include)] <- eche_separates_include[length(eche_separates_include)] + length(now_eche_locs[-c(seq_along(i_dummy))])
        }
        eche_separates <- c(eche_separates, length(now_eche_locs))
        now_eche_locs <- NULL
        pare_temp <- c(pare_temp, 0)
        break
      }

      else if(now_nei_max_val >= max(x[sub_nei])){
        now_eche_locs <- c(now_eche_locs, now_nei_max_locs)
        now_loc <- unique(c(now_loc, now_eche_locs))
      }

      else{
        total_eche_locs[which(total_eche_locs==0)[1] + seq_along(now_eche_locs) - 1] <- now_eche_locs
        dat[now_eche_locs] <- NA
        pare_temp <- c(pare_temp, now_nei_max_locs[1])

        if(p_dummy == 1){
          eche_peaks <- c(eche_peaks, 1)
          include_eche_locs <- c(include_eche_locs, now_eche_locs)
          eche_separates_include <- c(eche_separates_include, length(now_eche_locs))
        }
        else{
          eche_peaks <- c(eche_peaks, 0)
          include_eche_locs <- c(include_eche_locs, now_eche_locs[-c(seq_along(i_dummy))])
          eche_separates_include[length(eche_separates_include)] <- eche_separates_include[length(eche_separates_include)] + length(now_eche_locs[-c(seq_along(i_dummy))])
        }
        eche_separates <- c(eche_separates, length(now_eche_locs))
        now_eche_locs <- NULL
        s_dummy <- 1
      }
    }
  }

  eche_parent <- rep(seq_along(eche_separates), times = eche_separates)[match(pare_temp, total_eche_locs)]
  eche_parent[is.na(eche_parent)] <- 0
	

#-------------------------------------------------------------------------------
# Echelon scan
#-------------------------------------------------------------------------------
  e_num <- 1
  Kout <- 0
  nowrow <- 1
  nowcol <- 1
  c_separates <- c(0, cumsum(eche_separates))

# K>=1
  if(K >= 1){
    reg_data <- matrix(NA, length(x), K)
    while(e_num <= length(eche_peaks)){
      if(eche_peaks[e_num] == 1){
        areas <- total_eche_locs[(c_separates[e_num]+1):c_separates[e_num+1]]

        if(length(areas) > K){
          i <- 1
          temp <- NULL

          while(length(temp) < K){
            t_check <- x[areas[i]]
            t_check2 <- which(x[areas] == t_check)
            if(length(temp) + length(t_check2) <= K){
              temp <- c(temp, areas[t_check2])
              i <- length(temp) + 1
            }
            else{
              if(is.null(temp)) Kout <- 1
              break
            }
          }
          if(!is.null(temp)) areas <- temp
        }

        if(Kout == 0){
          temp <- t(array(areas,c(length(areas), length(areas))))
          temp[upper.tri(temp)] <- NA

          if(any(duplicated(x[areas]))){
            del_tie <- which(duplicated(x[areas])) - 1
            temp <- temp[-del_tie, , drop = FALSE]
          }

          reg_data[nowrow:(nowrow + nrow(temp) - 1), 1:ncol(temp)] <- temp
          nowrow <- nowrow + nrow(temp)
          if(ncol(temp) > nowcol) nowcol <- ncol(temp)
        }
        else Kout <- 0
      }

      else{
        areas <- total_eche_locs[(c_separates[e_num]+1):c_separates[e_num+1]]
        now_eches <- e_num
        temp <- NULL
        ii <- 1

        while(any(eche_peaks[now_eches] == 0)){
          child_eches <- which(is.element(eche_parent, now_eches))
          for(i in seq_along(child_eches)){
            temp2 <- total_eche_locs[(c_separates[child_eches[i]]+1):c_separates[child_eches[i]+1]]
            temp[ii:(ii + length(temp2) - 1)] <- temp2
            ii <- ii + length(temp2)
          }

          e_top <- x[total_eche_locs[(c_separates[e_num]+1):c_separates[e_num+1]]]
          if((length(temp) + length(which(e_top[1] == e_top))) > K){
            Kout <- 1
            break
          }
          now_eches <- child_eches
        }

        if(Kout == 0){
          if(length(temp) + length(areas) > K){
            i <- 1
            temp2 <- NULL
            while(length(temp) + length(temp2) <= K){
              t_check <- x[areas[i]]
              t_check2 <- which(x[areas] == t_check)
              if(length(temp) + length(temp2) + length(t_check2) <= K){
                temp2 <- c(temp2, areas[t_check2])
                i <- length(temp2) + 1
              }
              else break
            }
            areas <- temp2
          }

          temp2 <- t(array(temp, c(length(temp), length(areas))))
          temp <- t(array(areas, c(length(areas), length(areas))))
          temp[upper.tri(temp)] <- NA
          temp <- cbind(temp2, temp)

          if(any(duplicated(x[areas]))){
            del_tie <- which(duplicated(x[areas])) - 1
            temp <- temp[-del_tie, , drop = FALSE]
          }

          reg_data[nowrow:(nowrow + nrow(temp) - 1), 1:ncol(temp)] <- temp
          nowrow <- nowrow + nrow(temp)
          if(ncol(temp) > nowcol) nowcol <- ncol(temp)
        }
        else Kout <- 0
      }
      e_num <- e_num + 1
    }
    if(is.na(reg_data[length(x),1])) reg_data <- reg_data[-(nowrow:length(x)),]
    if(K > nowcol) reg_data <- reg_data[, -((nowcol+1):K)]
  }

# K<1
  else if(K < 1){
    reg_data <- matrix(NA, length(x), length(x))
    ng <- sum(par2)

    while(e_num <= length(eche_peaks)){
      if(eche_peaks[e_num] == 1){
        areas <- total_eche_locs[(c_separates[e_num]+1):c_separates[e_num+1]]

        if(sum(par2[areas]) > ng*K){
          i <- 1
          temp <- NULL
          while(sum(par2[temp]) <= ng*K){
            t_check <- x[areas[i]]
            t_check2 <- which(x[areas] == t_check)
            if(sum(par2[temp]) + sum(par2[areas[t_check2]]) <= ng*K){
              temp <- c(temp, areas[t_check2])
              i <- length(temp) + 1
            }
            else{
              if(is.null(temp)) Kout <- 1
              break
            }
          }
          if(!is.null(temp)) areas <- temp
        }

        if(Kout == 0){
          temp <- t(array(areas, c(length(areas), length(areas))))
          temp[upper.tri(temp)] <- NA

          if(any(duplicated(x[areas]))){
            del_tie <- which(duplicated(x[areas])) - 1
            temp <- temp[-del_tie, , drop = FALSE]
          }

          reg_data[nowrow:(nowrow + nrow(temp) - 1), 1:ncol(temp)] <- temp
          nowrow <- nowrow + nrow(temp)
          if(ncol(temp) > nowcol) nowcol <- ncol(temp)
        }
        else Kout <- 0
      }

      else{
        areas <- total_eche_locs[(c_separates[e_num]+1):c_separates[e_num+1]]
        now_eches <- e_num
        temp <- NULL
        ii <- 1

        while(any(eche_peaks[now_eches] == 0)){
          child_eches <- which(is.element(eche_parent, now_eches))
          for(i in seq_along(child_eches)){
            temp2 <- total_eche_locs[(c_separates[child_eches[i]]+1):c_separates[child_eches[i]+1]]
            temp[ii:(ii + length(temp2) - 1)] <- temp2
            ii <- ii + length(temp2)
          }

          e_top <- x[total_eche_locs[(c_separates[e_num]+1):c_separates[e_num+1]]]
          if(sum(par2[temp]) + sum(par2[areas[which(e_top[1] == e_top)]]) > ng*K){
            Kout <- 1
            break
          }
          now_eches <- child_eches
        }

        if(Kout == 0){
          if(sum(par2[temp]) + sum(par2[areas]) > ng*K){
            i <- 1
            temp2 <- NULL
            while(sum(par2[temp]) + sum(par2[temp2]) <= ng*K){
              t_check <- x[areas[i]]
              t_check2 <- which(x[areas] == t_check)
              if(sum(par2[temp]) + sum(par2[temp2]) + sum(par2[areas[t_check2]]) <= ng*K){
                temp2 <- c(temp2, areas[t_check2])
                i <- length(temp2) + 1
              }
              else break
            }
            areas <- temp2
          }

          temp2 <- t(array(temp, c(length(temp), length(areas))))
          temp <- t(array(areas, c(length(areas), length(areas))))
          temp[upper.tri(temp)] <- NA
          temp <- cbind(temp2, temp)

          if(any(duplicated(x[areas]))){
            del_tie <- which(duplicated(x[areas])) - 1
            temp <- temp[-del_tie, , drop = FALSE]
          }

          reg_data[nowrow:(nowrow + nrow(temp) - 1), 1:ncol(temp)] <- temp
          nowrow <- nowrow + nrow(temp)
          if(ncol(temp) > nowcol) nowcol <- ncol(temp)
        }
        else Kout <- 0
      }
      e_num <- e_num + 1
    }
    if(is.na(reg_data[length(x),1])) reg_data <- reg_data[-(nowrow:length(x)),]
    if(length(x) > nowcol) reg_data <- reg_data[, -((nowcol+1):length(x))]
  }

  if(!is.null(reg_data)){
    if(is.vector(reg_data)) reg_data <- t(reg_data)
    if(ncol(reg_data) < Kmin) reg_data <- NULL
    else if(Kmin != 1) reg_data <- reg_data[which(!is.na(reg_data[,Kmin])),]
    if(is.vector(reg_data)) reg_data <- t(reg_data)


#-------------------------------------------------------------------------------
# Statistic calculation
#-------------------------------------------------------------------------------
    if(!is.null(reg_data)){

# Poisson model
      if(type == 11 || type == 12){
        cg <- sum(cas)
        eg <- sum(ex)
        cz <- rowSums(array(cas[reg_data], dim(reg_data)), na.rm = TRUE)
        ez <- rowSums(array(ex[reg_data], dim(reg_data)), na.rm = TRUE)

        if(type == 11) temp <- which(cz > ez)
        else if(type == 12) temp <- which(cz < ez)

        cz <- cz[temp]
        if(length(cz) == 0) maxLLR <- 0
        else{
          ez <- ez[temp]
          term1 <- cz*log(cz/ez)
          term1[is.nan(term1)] <- 0
          term2 <- (cg-cz)*log((cg-cz)/(eg-ez))
          term2[is.nan(term2)] <- 0
          maxLLR <- max(term1 + term2)
        }
      }

# Binomial model
      if(type == 21 || type == 22){
        casg <- sum(cas)
        ctlg <- sum(ctl)
        casz <- rowSums(array(cas[reg_data], dim(reg_data)), na.rm = TRUE)
        ctlz <- rowSums(array(ctl[reg_data], dim(reg_data)), na.rm = TRUE)

        if(type == 21) temp <- which(casz/(casz + ctlz) > casg/(casg + ctlg))
        else if(type == 22) temp <- which(casz/(casz + ctlz) < casg/(casg + ctlg))

        casz <- casz[temp]
        if(length(casz) == 0) maxLLR <- 0
        else{
          ctlz <- ctlz[temp]
          popz <- casz + ctlz
          popg <- casg + ctlg
          term1 <- casz*log(casz/popz)
          term1[is.nan(term1)] <- 0
          term2 <- ctlz*log(ctlz/popz)
          term2[is.nan(term2)] <- 0
          term3 <- (casg-casz)*log((casg-casz)/(popg-popz))
          term3[is.nan(term3)] <- 0
          term4 <- (ctlg-ctlz)*log((ctlg-ctlz)/(popg-popz))
          term4[is.nan(term4)] <- 0
          term5 <- -casg*log(casg/popg)
          term5[is.nan(term5)] <- 0
          term6 <- -ctlg*log(ctlg/popg)
          term6[is.nan(term6)] <- 0
          maxLLR <- max(term1 + term2 + term3 + term4 + term5 + term6)
        }
      }

# Normal model

    }
    else maxLLR <- -1
  }
  else maxLLR <- -1

  return(maxLLR)
}

