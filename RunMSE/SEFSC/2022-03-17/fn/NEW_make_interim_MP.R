NEW_make_interim_MP<-function (.Assess = "SCA", .HCR = "HCR_MSY", AddInd = "VB", 
                              assessment_interval = 5, type = c("buffer", "mean", "loess", "none"), 
                              type_par = NULL, diagnostic = c("min", "full", "none"), 
                              const=1, eq="default", ...) 
{
  
  type <- match.arg(type)
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
  dots_in_Assess <- dots[match(names(formals(eval(.Assess))), 
                               names(dots), nomatch = 0)]
  dots_in_HCR <- dots[match(names(formals(eval(.HCR))), names(dots), 
                            nomatch = 0)]
  Assess_call <- as.call(c(.Assess, x = quote(x), Data = quote(Data), 
                           dots_in_Assess))
  HCR_call <- as.call(c(.HCR, Assessment = quote(do_Assessment), 
                        reps = quote(reps), dots_in_HCR))
  MP_body <- bquote({
    dependencies <- .(get_dependencies(Assess_char, dots))
    ny <- length(Data@Year)
    Current_Yr <- Data@Year[ny]
    first_assess_yr <- Current_Yr == Data@LHYear
    run_assessment <- first_assess_yr || Current_Yr == Data@Misc[[x]]$interim$next_assess_yr
    run_interim <- !run_assessment
    if (run_assessment) {
      do_Assessment <- .(Assess_call)
      if (diagnostic != "none") {
        Assess_diag_output <- Assess_diagnostic(x, Data, 
                                                do_Assessment, include_assessment = .(diagnostic == 
                                                                                        "full"))
      }
      if (do_Assessment@conv) {
        Rec <- .(HCR_call)
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
          Rec@Misc$interim <- list(Cref = NA_real_, Iref = NA_real_, 
                                   next_assess_yr = next_assess_yr)
          return(Rec)
        }
        else {
          run_interim <- TRUE
        }
      }
    }
    if (!run_assessment || run_interim) {
      Cref <- Data@Misc[[x]]$interim$Cref
      Iref <- Data@Misc[[x]]$interim$Iref
      I_y <- switch(type, buffer = Assess_I_hist(AddInd[1], Data, x, ny)$I_hist, 
                    mean = local({
                      if (is.null(type_par)) type_par <- 3
                      Assess_I_hist(AddInd[1], Data, x, seq(ny - type_par + 1, ny))$I_hist %>% mean(na.rm = TRUE)
                    }), 
                    loess = local({
                      I_df <- data.frame(Year = Data@Year, Ind = Assess_I_hist(AddInd[1], Data, x, 1:ny)$I_hist)
                      if (is.null(type_par)) type_par <- formals(loess)$span
                      fit <- loess(Ind ~ Year, data = I_df, span = type_par)
                      fit$fitted[length(fit$fitted)]
                    }), 
                    none = Assess_I_hist(AddInd[1], Data, x, ny)$I_hist)
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
      
      `%!in%` <- Negate(`%in%`)
      if(eq %!in% c("default","change","index")){ eq="default"}
      Rec <- new("Rec") ##########################################################################################
      if(eq=="default"){
        Rec@TAC <- const * Cref * (I_y + b * s)/(Iref + b * s)
      }
      if(eq=="change"){ ##############################################################
        # Rec@TAC <- Cref * (I_y + b * s)/(Iref + b * s)
        Rec@TAC <- ( ( (Cref * (I_y + b * s)/(Iref + b * s)) - Cref) * const) + Cref
      }
      if(eq=="index"){ ##############################################################
        # Rec@TAC <- Cref * (I_y + b * s)/(Iref + b * s)
        Rec@TAC <- Cref * (I_y + b * s)/((Iref*const) + b * s)
      }
      
      Rec@Misc <- Data@Misc[[x]]
      if (exists("next_assess_yr", inherits = FALSE)) 
        Rec@Misc$interim$next_assess_yr <- next_assess_yr
      if (exists("Assess_diag_output", inherits = FALSE)) 
        Rec@Misc$diagnostic <- Assess_diag_output$diagnostic
    }
    return(Rec)
  })
  custom_MP <- eval(call("function", as.pairlist(alist(x = 1, 
                                                       Data = , reps = 1)), MP_body))
  formals(custom_MP)$assessment_interval <- assessment_interval
  formals(custom_MP)$AddInd <- AddInd
  formals(custom_MP)$type <- type
  formals(custom_MP)$type_par <- type_par
  formals(custom_MP)$diagnostic <- diagnostic
  return(structure(custom_MP, class = "MP"))
}
class(NEW_make_interim_MP)<-"MP"
environment(NEW_make_interim_MP) <- asNamespace("SAMtool")