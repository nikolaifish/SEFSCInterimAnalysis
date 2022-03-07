# Author: Quang Huynh

projection_MP <- eval(bquote(function(x, Data, reps = 1, assessment_interval, HCR = HCR_MSY, SCA_arg = list(I_type = "VB", CAA_multiplier = 20)) {
  dependencies <- .(SAMtool:::get_dependencies("SCA_Pope"))

  current_yr <- Data@Year[length(Data@Year)]
  run_SCA <- current_yr == Data@LHYear
  if (current_yr > Data@LHYear) {
    run_SCA <- current_yr == Data@Misc[[x]]$next_assess_yr
  }
  if (!run_SCA) get_projected_TAC <- TRUE else get_projected_TAC <- FALSE

  if (run_SCA) {
    # Return TAC = UMSY * VB_current when run_SCA = TRUE
    SCA_formals <- list(x = x, Data = Data)
    do_Assessment <- do.call(SCA_Pope, c(SCA_formals, SCA_arg))
    Assess_output <- Assess_diagnostic(x, Data, do_Assessment, include_assessment = FALSE)

    if (do_Assessment@conv) {
      Rec <- do.call(match.fun(HCR), list(Assessment = do_Assessment, reps = reps))
      pro <- projection(do_Assessment, FMort = do_Assessment@UMSY, p_years = 52, p_sim = 1, obs_error = c(0, 0), process_error = 0)
      Rec@Misc <- c(list(project_catch = pro@Catch[1, ], last_assess_yr = current_yr, next_assess_yr = current_yr + assessment_interval), Assess_output)

      get_projected_TAC <- FALSE

    } else {
 
      if (current_yr == Data@LHYear || length(Data@Misc[[x]]) == 2) {
        Rec <- new("Rec")
        Rec@TAC <- TACfilter(rep(NA, reps))
        Rec@Misc <- c(list(next_assess_yr = current_yr + 1), Assess_output)

        get_projected_TAC <- FALSE
      } else {
        get_projected_TAC <- TRUE
        next_assess_yr <- current_yr + 1
      }
    }
  }

  if (get_projected_TAC) {
    Rec <- new("Rec")
    TAC_used <- Data@Misc[[x]]$project_catch[current_yr - Data@Misc[[x]]$last_assess_yr + 1]
    if(is.infinite(TAC_used) || is.na(TAC_used)) stop("Error in TAC during interim")
    Rec@TAC <- TACfilter(TAC_used)
    if (run_SCA) {
      Rec@Misc <- c(list(project_catch = Data@Misc[[x]]$project_catch, last_assess_yr = Data@Misc[[x]]$last_assess_yr,
                         next_assess_yr = ifelse(exists("next_assess_yr"), next_assess_yr, Data@Misc[[x]]$next_assess_yr)),
                    Assess_output)

    } else {
      Rec@Misc <- list(project_catch = Data@Misc[[x]]$project_catch, last_assess_yr = Data@Misc[[x]]$last_assess_yr,
                       next_assess_yr = ifelse(exists("next_assess_yr"), next_assess_yr, Data@Misc[[x]]$next_assess_yr),
                       diagnostic = Data@Misc[[x]]$diagnostic)
    }
  }
  return(Rec)

}))
class(projection_MP) <- "MP"
environment(projection_MP) <- asNamespace("SAMtool")
