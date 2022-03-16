make_projection_MP_lag <- function (.Assess = "SCA",
                                    .HCR = "HCR_MSY",
                                    assessment_interval = 5,
                                    Ftarget = expression(Assessment@FMSY),
                                    proj_args = list(process_error = 1, p_sim = 1),
                                    diagnostic = c("min", "full", "none"),
                                    ...)
{
  diagnostic <- match.arg(diagnostic)
  if (is.character(.Assess)) {
    Assess_char <- .Assess
    .Assess <- as.symbol(.Assess)
  }
  else {
    .Assess <- substitute(.Assess)
    Assess_char <- as.character(.Assess)
  }
  if (!inherits(eval(.Assess), "Assess")) {
    stop(paste(.Assess, "does not belong to class 'Assess'. Use: avail('Assess') to find eligible objects."))
  }
  dots <- list(...)

  if(!is.null(dots$lag)){
    lag <- dots$lag
  }else{
    lag <- 0
  }

  dots_in_Assess <- dots[match(names(formals(eval(.Assess))),
                               names(dots), nomatch = 0)]
  Assess_call <- as.call(c(.Assess, x = quote(x), Data = quote(Data),
                           dots_in_Assess))
  proj_args <- proj_args[match(names(formals(projection)),
                               names(proj_args), nomatch = 0)]
  proj_args$p_years <- assessment_interval + 1
  proj_call <- as.call(c(substitute(projection), Assessment = quote(Assessment),
                         Ftarget = quote(eval(Ftarget)), proj_args))
  MP_body <- bquote({
    dependencies <- .(SAMtool:::get_dependencies(Assess_char, dots))
    ny <- length(Data@Year)
    ny_Assess <- ny- .(lag)           # Number of year of Data available for assessment
    Current_Yr <- Data@Year[ny]
    first_assess_yr <- Current_Yr == Data@LHYear
    run_assessment <- first_assess_yr || Current_Yr == Data@Misc[[x]]$proj$next_assess_yr
    get_projection <- !run_assessment
    Rec <- new("Rec")
    if (run_assessment) {
      Assessment <- .(Assess_call)
      if (diagnostic != "none") {
        Rec@Misc <- SAMtool:::Assess_diagnostic(x, Data, Assessment,
                                      include_assessment = .(diagnostic == "full"))
      }
      if (!is.null(Assessment@info$Misc))
        Rec@Misc <- c(Rec@Misc, Assessment@info$Misc)
      if (Assessment@conv) {
        Ftarget <- eval(Ftarget)
        do_project <- .(proj_call)
        Rec@Misc$proj <- list(TACvec = apply(do_project@Catch,
                                             2, median) %>% structure(names = Current_Yr +
                                                                        1 + 0:assessment_interval), next_assess_yr = Current_Yr +
                                assessment_interval)
        Rec@TAC <- rep(Rec@Misc$proj$TACvec[1], reps)
      }
      else {
        Rec@Misc$proj <- list(TACvec = NA_real_, next_assess_yr = Current_Yr +
                                1)
        Rec@TAC <- rep(NA_real_, reps)
      }
    }
    if (get_projection) {
      Rec@Misc <- Data@Misc[[x]]
      Rec@TAC <- rep(Data@Misc[[x]]$proj$TACvec[as.character(Current_Yr +
                                                               1)], reps)
    }
    return(Rec)
  })
  custom_MP <- eval(call("function", as.pairlist(alist(x = 1,
                                                       Data = , reps = 1)), MP_body))
  formals(custom_MP)$Ftarget <- Ftarget
  formals(custom_MP)$assessment_interval <- assessment_interval
  formals(custom_MP)$diagnostic <- diagnostic
  return(structure(custom_MP, class = "MP"))
}
