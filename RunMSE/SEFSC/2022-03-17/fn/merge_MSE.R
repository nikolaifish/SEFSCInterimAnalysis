# Author: Quang Huynh
# Modified for compatability with newer version of SAMtool by Nikolai Klibansky

# Merge batches of MSE output
merge_MSE <- function(...) {
  dots <- list(...)
  
  slots_identical <- function(slotname, x = dots, is_logical = FALSE) {
    res <- lapply(x, getElement, slotname)
    is_identical <- all(vapply(res[-1], identical, logical(1), res[[1]]))
    if(is_logical) {
      return(is_identical)
    } else return(unique(do.call(c, res)))
  }
  
  slots_identical("Name")
  nyears <- slots_identical("nyears")
  proyears <- slots_identical("proyears")
  nsim <- slots_identical("nsim")
  nallyears <- nyears+proyears
  
  stopifnot(slots_identical("OM", is_logical = TRUE))
  stopifnot(slots_identical("Obs", is_logical = TRUE))
  stopifnot(slots_identical("SSB_hist", is_logical = TRUE))
  stopifnot(slots_identical("CB_hist", is_logical = TRUE))
  stopifnot(slots_identical("FM_hist", is_logical = TRUE))
  #stopifnot(slots_identical("Hist", is_logical = TRUE))
  
  nMPs <- vapply(dots, getElement, numeric(1), "nMPs")
  
  
  # These objects are all arrays (nsim, nMPs, proyears)
  slotvec <- c("SB_SBMSY", "F_FMSY", "N", "B", "SSB", "VB", "FM", "Catch", "Removals", "Effort",
               "TAC", "TAE")
  res <- list()
  for(i in 1:length(slotvec)) {
    new_mat <- array(NA, dim = c(nsim, sum(nMPs), proyears))
    for(j in 1:length(dots)) {
      array_ij <- getElement(dots[[j]], slotvec[i])
      if(j == 1) new_mat[, 1:nMPs[1], ] <- array_ij
      if(j > 1) new_mat[, (sum(nMPs[1:(j-1)]) + 1):(sum(nMPs[1:(j-1)]) + nMPs[j]), ] <- array_ij
    }
    res[[i]] <- new_mat
  }
  
  
  # These objects are basically lists of arrays (nlevels, nsim, nMPs, proyears)
  slotvec2 <- c("SPR", "BioEco")
  res2 <- list()
  for(i in 1:length(slotvec2)) {
    for(j in 1:length(dots)) {
      list_ij <- getElement(dots[[j]], slotvec2[i])
      if(slotvec2[i]=="RefPoint"){
        list_ij <- list_ij[c("MSY","FMSY","SSBMSY")]
      }
      new_mat <- array(NA, dim = c(nsim, sum(nMPs), proyears))
      for(k in 1:length(list_ij)) {
        array_ijk <- list_ij[[k]]
        if(j == 1) new_mat[, 1:nMPs[1], ] <- array_ijk
        if(j > 1) new_mat[, (sum(nMPs[1:(j-1)]) + 1):(sum(nMPs[1:(j-1)]) + nMPs[j]), ] <- array_ijk
        list_ij[[k]] <- new_mat
      }
      res2[[i]] <- list_ij
    }
  }
  
  # RefPoint is a lists of arrays and lists which vary in shape
  RefPointlist <- list()
  RefPointvec <- c("MSY", "FMSY", "SSBMSY","F_SPR")
  
  for(i in 1:length(RefPointvec)) {
    RefPoint_i <- RefPointvec[i]
    if(RefPoint_i!="F_SPR"){
      new_mat <- array(NA, dim = c(nsim, sum(nMPs), nallyears))
      for(j in 1:length(dots)) {
        array_ij <- getElement(slot(dots[[j]],"RefPoint"), RefPointvec[i])
        if(j == 1) new_mat[, 1:nMPs[1], ] <- array_ij
        if(j > 1) new_mat[, (sum(nMPs[1:(j-1)]) + 1):(sum(nMPs[1:(j-1)]) + nMPs[j]), ] <- array_ij
      }
    }else{
      new_mat <- array(NA, dim = c(nsim, sum(nMPs), 9, nallyears))
      for(j in 1:length(dots)) {
        array_ij <- getElement(slot(dots[[j]],"RefPoint"), RefPointvec[i])
        if(j == 1) new_mat[, 1:nMPs[1], , ] <- array_ij
        if(j > 1) new_mat[, (sum(nMPs[1:(j-1)]) + 1):(sum(nMPs[1:(j-1)]) + nMPs[j]), , ] <- array_ij
      }
    }
    RefPointlist[[RefPoint_i]] <- new_mat
  }
  # These should be these same for all MPs
  RefPointlist$Dynamic_Unfished <- getElement(slot(dots[[1]],"RefPoint"), "Dynamic_Unfished")
  RefPointlist$ByYear <- getElement(slot(dots[[1]],"RefPoint"), "ByYear")
  
  # PPD is a list of Data class objects (nMPs)
  PPDlist <- list()
  ct <- 0
  for(j in 1:length(dots)) {
    list_j <- getElement(dots[[j]], "PPD")
    for(k in 1:length(list_j)){
      ct <- ct+1
      PPDlist[[ct]] <- list_j[[k]]
    }
  }
  
  ## Create MSE Object ####
  MSEout <- new("MSE",
                Name = slots_identical("Name"),
                nyears = nyears,
                proyears = proyears,
                nMPs = length(slots_identical("MPs")),
                MPs = slots_identical("MPs"),
                nsim = nsim,
                OM = dots[[1]]@OM,
                Obs = dots[[1]]@Obs,
                SB_SBMSY = res[[1]],
                F_FMSY = res[[2]],
                N = res[[3]],
                B = res[[4]],
                SSB = res[[5]],
                VB = res[[6]],
                FM = res[[7]],
                SPR = res2[[1]],
                Catch = res[[8]],
                Removals = res[[9]],
                Effort = res[[10]],
                TAC = res[[11]],
                TAE = res[[12]],
                BioEco = res2[[2]],
                RefPoint = RefPointlist,
                SSB_hist = dots[[1]]@SSB_hist,
                CB_hist = dots[[1]]@CB_hist,
                FM_hist = dots[[1]]@FM_hist,
                Hist = dots[[1]]@Hist,
                PPD = PPDlist,
                Misc = list(Data = do.call(c, lapply(dots, function(x) x@Misc$Data))))
  
  # Store MSE info
  attr(MSEout, "version") <- packageVersion("DLMtool")
  attr(MSEout, "date") <- date()
  attr(MSEout, "R.version") <- R.version
  
  MSEout
}