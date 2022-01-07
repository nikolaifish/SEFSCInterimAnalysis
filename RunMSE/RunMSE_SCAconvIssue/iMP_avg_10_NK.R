iMP_avg_10_NK <- function (x = 1, Data, reps = 1, assessment_interval = 10, AddInd = "VB", 
          type = "mean", type_par = 3, diagnostic = "min") 
{
  dependencies <- "Data@Cat, Data@Ind, Data@Mort, Data@L50, Data@L95, Data@CAA, Data@vbK, Data@vbLinf, Data@vbt0, Data@wla, Data@wlb, Data@MaxAge, Data@steep, Data@CV_Ind, Data@sigmaR, Data@CV_Cat, Data@Mort, Data@CV_Mort, Data@steep, Data@CV_steep, Data@vbLinf, Data@vbK, Data@vbt0, Data@wla, Data@wlb, Data@MaxAge, Data@L50, Data@L95"
  interim_get_index <- SAMtool:::interim_get_index
  Assess_diagnostic <- SAMtool:::Assess_diagnostic
  
  
  ny <- length(Data@Year)
  Current_Yr <- Data@Year[ny]
  first_assess_yr <- Current_Yr == Data@LHYear
  run_assessment <- first_assess_yr || Current_Yr == Data@Misc[[x]]$interim$next_assess_yr
  run_interim <- !run_assessment
  print(paste("iMP_avg_10_NK: Current_Yr=",Current_Yr))
  print(paste("iMP_avg_10_NK: run_assessment=",run_assessment))
  if (run_assessment) {
    do_Assessment <- SCA_Pope_NK(x = x, Data = Data, AddInd = "VB")
    if (diagnostic != "none") {
      Assess_diag_output <- Assess_diagnostic(x, Data, 
                                              do_Assessment, include_assessment = FALSE)
    }
    if (do_Assessment@conv) {
      Rec <- HCR_MSY(Assessment = do_Assessment, reps = reps, 
                     MSY_frac = 1)
      if (diagnostic != "none") 
        Rec@Misc$diagnostic <- Assess_diag_output$diagnostic
      if (!is.null(do_Assessment@info$Misc)) 
        Rec@Misc <- c(Rec@Misc, do_Assessment@info$Misc)
      Rec@Misc$interim <- list(Cref = Rec@TAC, Iref = do_Assessment@Index[ny, 
                                                                          1], next_assess_yr = Current_Yr + assessment_interval)
      if (type == "buffer") {
        Rec@Misc$interim$s <- log(do_Assessment@Obs_Index[, 
                                                          1]/do_Assessment@Index[, 1]) %>% sd(na.rm = TRUE)
      }
      else {
        Rec@Misc$interim$s <- 0
      }
      return(Rec)
    }
    else {
      next_assess_yr <- Current_Yr + 1
      if (first_assess_yr || is.null(Data@Misc[[x]]$interim)) {
        Rec <- new("Rec")
        Rec@TAC <- rep(NA_real_, reps)
        if (diagnostic != "none") 
          Rec@Misc$diagnostic <- Assess_diag_output$diagnostic
        if (!is.null(do_Assessment@info$Misc)) 
          Rec@Misc <- c(Rec@Misc, do_Assessment@info$Misc)
        Rec@Misc$interim <- list(next_assess_yr = next_assess_yr)
        return(Rec)
      }
      else {
        run_interim <- TRUE
      }
    }
  }
  if (!run_assessment || run_interim) {
    print("iMP_avg_10_NK: running interim")
    Cref <- Data@Misc[[x]]$interim$Cref
    Iref <- Data@Misc[[x]]$interim$Iref
    I_y <- switch(type, buffer = interim_get_index(AddInd[1],x, Data, ny),
                  mean = local({
                  if (is.null(type_par)) type_par <- 3
                  interim_get_index(AddInd[1], x, Data, seq(ny - type_par + 1, ny)) %>% mean(na.rm = TRUE)
                                                   }),
                  loess = local({
                    I_df <- data.frame(Year = Data@Year, Ind = interim_get_index(AddInd[1], x, Data, 1:ny))
                    if (is.null(type_par)) type_par <- formals(loess)$span
                    fit <- loess(Ind ~ Year, data = I_df, span = type_par)
                    fit$fitted[length(fit$fitted)]
                    }),
                  none = interim_get_index(AddInd[1], x, Data, ny)
                  )
    print("iMP_avg_10_NK: got I_y")
    myI_y<<-I_y
    if (type == "buffer") {
      if (is.null(type_par)) {
        b <- 1
      }
      else {
        b <- type_par
      }
      s <- Data@Misc[[x]]$interim$s
    }
    else {
      b <- s <- 0
    }
    Rec <- new("Rec")
    Rec@TAC <- Cref * (I_y + b * s)/(Iref + b * s)
    Rec@Misc <- Data@Misc[[x]]
    if (exists("next_assess_yr", inherits = FALSE)) 
      Rec@Misc$interim$next_assess_yr <- next_assess_yr
    if (exists("Assess_diag_output", inherits = FALSE)) 
      Rec@Misc$diagnostic <- Assess_diag_output$diagnostic
  }
  return(Rec)
}
iMP_avg_10_NK <- structure(iMP_avg_10_NK, class = "MP")
