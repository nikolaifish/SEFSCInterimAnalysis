Assess_diagnostic_NK <- function (x, Data, Assessment, include_assessment = TRUE)
{
  Year <- Data@Year[length(Data@Year)]
  if (inherits(Assessment, "Assessment")) {
    msg <- ifelse(is.character(Assessment@opt), Assessment@opt,
                  Assessment@opt$message)
    hess <- ifelse(is.character(Assessment@SD), FALSE, Assessment@SD$pdHess)
    maxgrad <- ifelse(is.character(Assessment@SD), 1e+10,
                      max(abs(Assessment@SD$gradient.fixed)))
    iter <- ifelse(is.character(Assessment@opt), NA, Assessment@opt$iterations)
    fn_eval <- ifelse(is.character(Assessment@opt), NA,
                      Assessment@opt$evaluations[1])

    # Stuff added by Nikolai Klibansky (2022-03-18)
    conv <- slot(Assessment,"conv")
    R0 <- slot(Assessment,"R0")
    h <- slot(Assessment,"h")
    MSY <- slot(Assessment,"MSY")
    FSY <- slot(Assessment,"FMSY")
    SSBMSY <- slot(Assessment,"SSBMSY")
    F_FMSY <- tail(slot(Assessment,"F_FMSY"),1)
    SSB_SSBMSY <- tail(slot(Assessment,"SSB_SSBMSY"),1)

    dg <- list(hess = hess, msg = msg, maxgrad = maxgrad,
               iter = iter, fn_eval = fn_eval, Year = Year,
               conv = conv,
               R0 = R0, h = h,
               MSY = MSY, FSY = FSY, SSBMSY = SSBMSY,
               F_FMSY = F_FMSY, SSB_SSBMSY = SSB_SSBMSY
    )
  }
  else {
    dg <- list(hess = FALSE, msg = msg, maxgrad = 1e+10,
               iter = NA, fn_eval = NA, Year = Year,
               conv = NA,
               R0 = NA, h = NA,
               MSY = NA, FSY = NA, SSBMSY = NA,
               F_FMSY = NA, SSB_SSBMSY = NA
    )
  }
  if (length(Data@Misc) == 0)
    Data@Misc <- vector("list", length(Data@Mort))
  diagnostic <- Data@Misc[[x]]$diagnostic
  len_diag <- length(diagnostic)
  diagnostic[[len_diag + 1]] <- dg
  output <- list(diagnostic = diagnostic)
  if (include_assessment) {
    if (inherits(Assessment, "Assessment")) {
      Assessment@info <- Assessment@obj <- list()
      Assessment@SD <- character(0)
      Assessment@N_at_age <- Assessment@C_at_age <- Assessment@Obs_C_at_age <- Assessment@Selectivity <- array(dim = c(0,
                                                                                                                       0))
    }
    Assessment_report <- Data@Misc[[x]]$Assessment_report
    len_Assess <- length(Assessment_report)
    if (len_Assess > 0) {
      Assessment_report[[len_Assess + 1]] <- Assessment
    }
    else Assessment_report <- list(Assessment)
    output$Assessment_report <- Assessment_report
  }
  return(output)
}
