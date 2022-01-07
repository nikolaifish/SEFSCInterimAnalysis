# "perfect" MP suggested by Quang Huynh to investigate deviation of HCR_MSY from MSY reference points
MP_PerfectMSY <- make_MP(Perfect, HCR_MSY)

# Averaged Index and Buffered Index MPs with (numbers indicate assessment interval)
iMP_avg_5 <- make_interim_MP(.Assess = "SCA_Pope", .HCR = "HCR_MSY", MSY_frac = 1,
                             assessment_interval = 5,
                             type = "mean", type_par = 3)
iMP_avg_10 <- make_interim_MP(SCA_Pope, .HCR = "HCR_MSY", MSY_frac = 1,
                              assessment_interval = 10,
                              type = "mean", type_par = 3)

iMP_buffer_5 <- make_interim_MP(.Assess = "SCA_Pope", .HCR = "HCR_MSY", MSY_frac = 1,
                                assessment_interval = 5,
                                type = "buffer", type_par = 1)
iMP_buffer_10 <- make_interim_MP(.Assess = "SCA_Pope", .HCR = "HCR_MSY", MSY_frac = 1,
                                 assessment_interval = 10,
                                 type = "buffer", type_par = 1)

# Fixed TAC MPs (numbers indicate assessment interval) and annual assessment MPs
SCA_5 <- SCA_10 <- SCA_1 <- make_MP(.Assess = "SCA_Pope", .HCR = "HCR_MSY", MSY_frac = 1,
                                    diagnostic = "min")

# # Merge batches of MSE output (modified for newer version of SAMtool by Nikolai)
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
  stopifnot(slots_identical("Hist", is_logical = TRUE))

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


# projection_MP <- eval(bquote(function(x, Data, reps = 1, assessment_interval, HCR = HCR_MSY, SCA_arg = list(I_type = "VB", CAA_multiplier = 20)) {
#   dependencies <- .(SAMtool:::get_dependencies("SCA_Pope"))
#
#   current_yr <- Data@Year[length(Data@Year)]
#   run_SCA <- current_yr == Data@LHYear
#   if (current_yr > Data@LHYear) {
#     run_SCA <- current_yr == Data@Misc[[x]]$next_assess_yr
#   }
#   if (!run_SCA) get_projected_TAC <- TRUE else get_projected_TAC <- FALSE
#
#   if (run_SCA) {
#     # Return TAC = UMSY * VB_current when run_SCA = TRUE
#     SCA_formals <- list(x = x, Data = Data)
#     do_Assessment <- do.call(SCA_Pope, c(SCA_formals, SCA_arg))
#     Assess_output <- Assess_diagnostic(x, Data, do_Assessment, include_assessment = FALSE)
#
#     if (do_Assessment@conv) {
#       Rec <- do.call(match.fun(HCR), list(Assessment = do_Assessment, reps = reps))
#       pro <- projection(do_Assessment, FMort = do_Assessment@UMSY, p_years = 52, p_sim = 1, obs_error = c(0, 0), process_error = 0)
#       Rec@Misc <- c(list(project_catch = pro@Catch[1, ], last_assess_yr = current_yr, next_assess_yr = current_yr + assessment_interval), Assess_output)
#
#       get_projected_TAC <- FALSE
#
#     } else {
#
#       if (current_yr == Data@LHYear || length(Data@Misc[[x]]) == 2) {
#         Rec <- new("Rec")
#         Rec@TAC <- TACfilter(rep(NA, reps))
#         Rec@Misc <- c(list(next_assess_yr = current_yr + 1), Assess_output)
#
#         get_projected_TAC <- FALSE
#       } else {
#         get_projected_TAC <- TRUE
#         next_assess_yr <- current_yr + 1
#       }
#     }
#   }
#
#   if (get_projected_TAC) {
#     Rec <- new("Rec")
#     TAC_used <- Data@Misc[[x]]$project_catch[current_yr - Data@Misc[[x]]$last_assess_yr + 1]
#     if(is.infinite(TAC_used) || is.na(TAC_used)) stop("Error in TAC during interim")
#     Rec@TAC <- TACfilter(TAC_used)
#     if (run_SCA) {
#       Rec@Misc <- c(list(project_catch = Data@Misc[[x]]$project_catch, last_assess_yr = Data@Misc[[x]]$last_assess_yr,
#                          next_assess_yr = ifelse(exists("next_assess_yr"), next_assess_yr, Data@Misc[[x]]$next_assess_yr)),
#                     Assess_output)
#
#     } else {
#       Rec@Misc <- list(project_catch = Data@Misc[[x]]$project_catch, last_assess_yr = Data@Misc[[x]]$last_assess_yr,
#                        next_assess_yr = ifelse(exists("next_assess_yr"), next_assess_yr, Data@Misc[[x]]$next_assess_yr),
#                        diagnostic = Data@Misc[[x]]$diagnostic)
#     }
#   }
#   return(Rec)
#
# }))
# class(projection_MP) <- "MP"
# environment(projection_MP) <- asNamespace("SAMtool")
#
#
# make_projection_MP <- function(...) {
#   fn <- projection_MP
#   dots <- list(...)
#   arg_ind <- pmatch(names(dots), names(formals(fn)))
#   formals(fn)[arg_ind] <- dots
#   class(fn) <- "MP"
#   return(fn)
# }
#
# # Projection MPs
# pMP_5 <- make_projection_MP(assessment_interval = 5)
# pMP_10 <- make_projection_MP(assessment_interval = 10)
