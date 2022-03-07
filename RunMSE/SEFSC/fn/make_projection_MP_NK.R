# Author: Quang Huynh
# Modified and totally wrecked: Nikolai Klibansky
# I was trying to modify this code to make it usable with the current version of openMSE, but there are some settings
# I don't readily understand and Quang says he can probably code it up pretty quickly. So I'm giving up


make_projection_MP <- function (.Assess = "SCA",
                                .HCR = "HCR_MSY",
                                AddInd = 1,
                                assessment_interval = 5,
                                diagnostic = c("min", "full","none"),
                                ...) {
  diagnostic <- match.arg(diagnostic)
  if (is.character(.Assess)) {
    Assess_char <- .Assess
    .Assess <- as.symbol(.Assess)
  }
  else {
    .Assess <- substitute(.Assess)
    Assess_char <- as.character(.Assess)
  }
  if (is.character(.HCR)) {
    .HCR <- as.symbol(.HCR)
  }
  else {
    .HCR <- substitute(.HCR)
  }
  if (!inherits(eval(.Assess), "Assess")) {
    stop(paste(.Assess, "does not belong to class 'Assess'. Use: avail('Assess') to find eligible objects."))
  }
  if (!inherits(eval(.HCR), "HCR")) {
    stop(paste(.HCR, "does not belong to class 'HCR.' Use: avail('HCR') to find eligible objects."))
  }
  dots <- list(...)
  dots$AddInd <- AddInd

  if(!is.null(dots$lag)){
    lag <- dots$lag
  }else{
    lag <- 0
  }

  dots_in_Assess <- dots[match(names(formals(eval(.Assess))),
                               names(dots), nomatch = 0)]
  dots_in_HCR <- dots[match(names(formals(eval(.HCR))), names(dots),
                            nomatch = 0)]
  Assess_call <- as.call(c(.Assess, x = quote(x), Data = quote(Data),
                           dots_in_Assess))
  HCR_call <- as.call(c(.HCR, Assessment = quote(do_Assessment),
                        reps = quote(reps), dots_in_HCR))
  MP_body <- bquote({
    dependencies <- .(SAMtool:::get_dependencies(Assess_char, dots))
    ny <- length(Data@Year)      # Number of years in available Data object
    ny_Assess <- ny - .(lag)     # Number of year of Data available for assessment
    Current_Yr <- Data@Year[ny]  # Current year
    first_assess_yr <- Current_Yr == Data@LHYear # Is the current year the last year of available data?
    run_assessment <- first_assess_yr || Current_Yr == Data@Misc[[x]]$interim$next_assess_yr # Should the assessment be run?
    #run_interim <- !run_assessment # If the assessment should not be run this year, then run the interim analysis.


    # Current_Yr <- Data@Year[length(Data@Year)]
    # run_assessment <- Current_Yr == Data@LHYear
    # if (Current_Yr > Data@LHYear) {
    #   run_assessment <- Current_Yr == Data@Misc[[x]]$next_assess_yr
    # }
    if (!run_assessment) {
      get_projected_TAC <- TRUE
      }else{
        get_projected_TAC <- FALSE
      }

    if (run_assessment) {
      # Return TAC = UMSY * VB_current when run_assessment = TRUE
      # SCA_formals <- list(x = x, Data = Data)
      # do_Assessment <- do.call(SCA_Pope, c(SCA_formals, SCA_arg))
      do_Assessment <- .(Assess_call) # Run the assessment
      # Assess_output <- Assess_diagnostic(x, Data, do_Assessment, include_assessment = FALSE)

      if (do_Assessment@conv) {
        if (do_Assessment@conv) {
          Rec <- .(HCR_call)
          if (diagnostic != "none") {
            Rec@Misc <- SAMtool:::Assess_diagnostic(x, Data, do_Assessment,
                                                    include_assessment = .(diagnostic == "full"))
          }

        Rec <- do.call(match.fun(HCR), list(Assessment = do_Assessment, reps = reps))
        pro <- projection(do_Assessment, FMort = do_Assessment@UMSY, p_years = 52, p_sim = 1, obs_error = c(0, 0), process_error = 0)
        Rec@Misc <- c(list(project_catch = pro@Catch[1, ], last_assess_yr = Current_Yr, next_assess_yr = Current_Yr + assessment_interval), Assess_output)

        get_projected_TAC <- FALSE

      } else {

        if (Current_Yr == Data@LHYear || length(Data@Misc[[x]]) == 2) {
          Rec <- new("Rec")
          Rec@TAC <- TACfilter(rep(NA, reps))
          Rec@Misc <- c(list(next_assess_yr = Current_Yr + 1), Assess_output)

          get_projected_TAC <- FALSE
        } else {
          get_projected_TAC <- TRUE
          next_assess_yr <- Current_Yr + 1
        }
      }
    }

    if (get_projected_TAC) {
      Rec <- new("Rec")
      TAC_used <- Data@Misc[[x]]$project_catch[Current_Yr - Data@Misc[[x]]$last_assess_yr + 1]
      if(is.infinite(TAC_used) || is.na(TAC_used)) stop("Error in TAC during interim")
      Rec@TAC <- TACfilter(TAC_used)
      if (run_assessment) {
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

  })
  custom_MP <- eval(call("function", as.pairlist(alist(x = 1, Data = , reps = 1)), MP_body))
  formals(custom_MP)$assessment_interval <- assessment_interval
  # formals(custom_MP)$AddInd <- AddInd
  # formals(custom_MP)$type <- type
  # formals(custom_MP)$type_par <- type_par
  # formals(custom_MP)$diagnostic <- diagnostic
  return(structure(custom_MP, class = "MP"))



  # fn <- projection_MP
  # dots <- list(...)
  # arg_ind <- pmatch(names(dots), names(formals(fn)))
  # formals(fn)[arg_ind] <- dots
  # class(fn) <- "MP"
  # return(fn)
}
