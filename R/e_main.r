#-------------------------------------------------------------------------------
# Main program
#-------------------------------------------------------------------------------
e.main <- function(x, rin){
	
  dat <- x

  total_eche_locs <- numeric(length(dat))
  include_eche_locs <- NULL
  eche_separates <- NULL
  eche_separates_include <- NULL
  eche_peaks <- NULL
  eche_tops <- NULL
  eche_bottoms <- NULL
  pare_temp <- NULL
  eche_pro <- NULL
  eche_family <- NULL
  eche_family_include <- NULL

  while(any(!is.na(dat))){
    now_val <- max(dat, na.rm = TRUE)
    now_loc <- which(dat == now_val)[1]
    eche_tops <- c(eche_tops, now_loc)
		
    t_dummy <- now_loc

    now_eche_locs <- now_loc
    p_dummy <- 1
    f_dummy <- 1

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

        f_dummy <- f_dummy + sum(eche_family_include[del_eche_separates_include])
        eche_family_include <- eche_family_include[-del_eche_separates_include]

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
    eche_family_include <- c(eche_family_include, f_dummy)		

    s_dummy <- 0

    while(s_dummy == 0){
      now_nei <- as.vector(rin[now_loc,])
      now_nei <- now_nei[!is.na(now_nei)]
      if(any(is.element(now_nei, now_loc))) now_nei <- now_nei[-which(is.element(now_nei, now_loc))]
      now_nei <- unique(now_nei)

      if(length(now_nei) == 0){
        dat[now_eche_locs] <- NA
        total_eche_locs[which(total_eche_locs == 0)[1] + seq_along(now_eche_locs) - 1] <- now_eche_locs
        if(p_dummy == 1) eche_peaks <- c(eche_peaks, 1)
        else eche_peaks <- c(eche_peaks, 0)
        eche_separates <- c(eche_separates, length(now_eche_locs))
        eche_bottoms <- c(eche_bottoms, now_eche_locs[length(now_eche_locs)])
        now_eche_locs <- NULL
        eche_pro <- c(eche_pro, p_dummy-1)
        eche_family <- c(eche_family, f_dummy)
        pare_temp <- c(pare_temp, 0)
        break
      }

      now_nei_max_val <- max(dat[now_nei])
      now_nei_max_locs <- now_nei[which(now_nei_max_val == dat[now_nei])]
      t_dummy <- now_nei_max_val

      while(!is.na(t_dummy)){
        sub_nei <- c(now_nei, as.vector(rin[now_nei_max_locs,]))
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
        total_eche_locs[which(total_eche_locs == 0)[1] + seq_along(now_eche_locs) - 1] <- now_eche_locs

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
        eche_bottoms <- c(eche_bottoms, now_eche_locs[length(now_eche_locs)])
        now_eche_locs <- NULL
        eche_pro <- c(eche_pro, p_dummy - 1)
        eche_family <- c(eche_family, f_dummy)
        pare_temp <- c(pare_temp, 0)
        break
      }

      else if(now_nei_max_val >= max(x[sub_nei])){
        now_eche_locs <- c(now_eche_locs, now_nei_max_locs)
        now_loc <- unique(c(now_loc, now_eche_locs))
      }

      else{
        total_eche_locs[which(total_eche_locs == 0)[1] + seq_along(now_eche_locs) - 1] <- now_eche_locs
        dat[now_eche_locs] <- NA
        pare_temp <- c(pare_temp, now_nei_max_locs[1])

        if(p_dummy == 1){
          eche_peaks <- c(eche_peaks, 1)
          include_eche_locs <- c(include_eche_locs, now_eche_locs)
          eche_separates_include <- c(eche_separates_include, length(now_eche_locs))
          eche_pro <- c(eche_pro, 0)
        }
        else{
          eche_peaks <- c(eche_peaks, 0)
          include_eche_locs <- c(include_eche_locs, now_eche_locs[-c(seq_along(i_dummy))])
          eche_separates_include[length(eche_separates_include)] <- eche_separates_include[length(eche_separates_include)] + length(now_eche_locs[-c(seq_along(i_dummy))])
          eche_pro <- c(eche_pro, p_dummy - 1)
        }
        eche_separates <- c(eche_separates, length(now_eche_locs))
        eche_bottoms <- c(eche_bottoms, now_eche_locs[length(now_eche_locs)])
        eche_family <- c(eche_family, f_dummy)
        now_eche_locs <- NULL
        s_dummy <- 1
      }
    }
  }

  eche_parent <- rep(seq_along(eche_separates), times = eche_separates)[match(pare_temp, total_eche_locs)]
  eche_parent[is.na(eche_parent)] <- 0
	
  list(locs = total_eche_locs, peaks = eche_peaks, separates = eche_separates, parents = eche_parent,
       tops = eche_tops, bottoms = eche_bottoms, progeny = eche_pro, family = eche_family, pare_locs = pare_temp)
}
