SCA_NK <- function (x = 1, Data, AddInd = "B", SR = c("BH", "Ricker", 
                                            "none"), vulnerability = c("logistic", "dome"), catch_eq = c("Baranov", 
                                                                                                         "Pope"), comp = c("age", "length"), comp_dist = c("multinomial", 
                                                                                                                                                           "lognormal"), comp_multiplier = c(50, 50), rescale = "mean1", 
          max_age = Data@MaxAge, start = NULL, prior = list(), fix_h = TRUE, 
          fix_F_equilibrium = TRUE, fix_omega = TRUE, fix_tau = TRUE, 
          tv_M = c("none", "walk", "DD"), M_bounds = NULL, refyear = 1, 
          LWT = list(), early_dev = c("comp_onegen", "comp", "all"), 
          late_dev = "comp50", integrate = FALSE, silent = TRUE, opt_hess = FALSE, 
          n_restart = ifelse(opt_hess, 0, 1), control = list(iter.max = 2e+05, 
                                                             eval.max = 4e+05), inner.control = list(), ...) 
{
  print("running SCA_NK")
  Assess_I_hist <- SAMtool:::Assess_I_hist
  make_prior <- SAMtool:::make_prior
  tiny_comp <- SAMtool:::tiny_comp
  logit <- SAMtool:::logit
  optimize_TMB_model <- SAMtool:::optimize_TMB_model
  SCA_dynamic_SSB0 <- SAMtool:::SCA_dynamic_SSB0
  ref_pt_SCA <- SAMtool:::ref_pt_SCA

  dependencies <- "Data@Cat, Data@Ind, Data@Mort, Data@L50, Data@L95, Data@CAA, Data@vbK, Data@vbLinf, Data@vbt0, Data@wla, Data@wlb, Data@MaxAge"
  dots <- list(...)
  start <- lapply(start, eval, envir = environment())
  max_age <- as.integer(min(max_age, Data@MaxAge))
  n_age <- max_age + 1
  vulnerability <- match.arg(vulnerability)
  catch_eq <- match.arg(catch_eq)
  comp <- match.arg(comp, several.ok = TRUE)
  comp_dist <- match.arg(comp_dist)
  if (length(comp_multiplier) == 1) 
    comp_multiplier <- rep(comp_multiplier, 2)
  SR <- match.arg(SR)
  tv_M <- match.arg(tv_M)
  if (is.character(early_dev)) 
    early_dev <- match.arg(early_dev)
  if (is.numeric(early_dev)) 
    stopifnot(early_dev < length(Data@Year))
  if (any(names(dots) == "yind")) {
    yind <- eval(dots$yind)
  }
  else {
    yind <- which(!is.na(Data@Cat[x, ]))[1]
    yind <- yind:length(Data@Cat[x, ])
  }
  Year <- Data@Year[yind]
  C_hist <- Data@Cat[x, yind]
  if (any(is.na(C_hist) | C_hist < 0)) 
    warning("Error. Catch time series is not complete.")
  n_y <- length(C_hist)
  M <- rep(Data@Mort[x], n_age)
  a <- Data@wla[x]
  b <- Data@wlb[x]
  Linf <- Data@vbLinf[x]
  K <- Data@vbK[x]
  t0 <- Data@vbt0[x]
  La <- Linf * (1 - exp(-K * (c(0:max_age) - t0)))
  SD_La <- La * Data@LenCV[x]
  Wa <- a * La^b
  A50 <- min(0.5 * max_age, iVB(t0, K, Linf, Data@L50[x]))
  A95 <- max(A50 + 0.5, iVB(t0, K, Linf, Data@L95[x]))
  mat_age <- c(0, 1/(1 + exp(-log(19) * (c(1:max_age) - A50)/(A95 - 
                                                                A50))))
  mat_age <- mat_age/max(mat_age)
  if (any(comp == "age")) {
    Data <- SAMtool:::expand_comp_matrix(Data, "CAA")
    print(paste0("x=",x))
    #print(paste0("yind=",paste(yind,collapse=" ")))
    #print(paste0("n_age=",n_age))
    myCAA <<- Data@CAA

    CAA_hist <- Data@CAA[x, yind, 1:n_age]

    if (max_age < Data@MaxAge) 
      CAA_hist[, n_age] <- rowSums(Data@CAA[x, yind, n_age:(Data@MaxAge + 
                                                              1)], na.rm = TRUE)
  }
  else {
    CAA_hist <- matrix(0, n_y, n_age)
  }
  CAA_n_nominal <- rowSums(CAA_hist, na.rm = TRUE)
  if (comp_multiplier[1] <= 1) {
    CAA_n_rescale <- comp_multiplier[1] * CAA_n_nominal
  }
  else {
    CAA_n_rescale <- pmin(CAA_n_nominal, comp_multiplier[1])
  }
  if (any(comp == "length")) {
    CAL_hist <- Data@CAL[x, yind, ]
    CAL_sum <- colSums(CAL_hist, na.rm = TRUE)
    CAL_cdf <- cumsum(CAL_sum)
    CAL_ind <- which(CAL_sum > 0)[1]:which.max(CAL_cdf)[1]
    CAL_hist <- CAL_hist[, CAL_ind]
    CAL_mids <- Data@CAL_mids[CAL_ind]
    CAL_bins <- Data@CAL_bins[CAL_ind]
    PLA <- generate_PLA(La, SD_La, CAL_bins, CAL_mids)
  }
  else {
    CAL_hist <- matrix(0, n_y, 1)
    PLA <- matrix(1, n_age, 1)
    CAL_bins <- CAL_mids <- 1
  }
  CAL_n_nominal <- rowSums(CAL_hist, na.rm = TRUE)
  if (comp_multiplier[2] <= 1) {
    CAL_n_rescale <- comp_multiplier[2] * CAL_n_nominal
  }
  else {
    CAL_n_rescale <- pmin(CAL_n_nominal, comp_multiplier[2])
  }
  if (early_dev == "all") {
    est_early_rec_dev <- rep(1, n_age - 1)
    est_rec_dev <- rep(1, n_y)
  }
  else if (early_dev == "comp") {
    est_early_rec_dev <- rep(0, n_age - 1)
    if (any(comp == "age")) {
      est_rec_dev <- ifelse(1:n_y < which(CAA_n_nominal > 
                                            0)[1], 0, 1)
    }
    else {
      est_rec_dev <- ifelse(1:n_y < which(CAL_n_nominal > 
                                            0)[1], 0, 1)
    }
  }
  else if (early_dev == "comp_onegen") {
    if (any(comp == "age")) {
      istart <- which(CAA_n_nominal > 0)[1] - n_age
    }
    else {
      istart <- which(CAL_n_nominal > 0)[1] - n_age
    }
    if (istart < 0) {
      early_start <- n_age + istart
      est_early_rec_dev <- ifelse(2:n_age < early_start, 
                                  0, 1) %>% rev()
      est_rec_dev <- rep(1, n_y)
    }
    else {
      est_early_rec_dev <- rep(0, n_age - 1)
      est_rec_dev <- ifelse(1:n_y < istart, 0, 1)
    }
  }
  else if (is.numeric(early_dev)) {
    if (early_dev > 1) {
      est_early_rec_dev <- rep(0, n_age - 1)
      est_rec_dev <- ifelse(1:n_y < early_dev, 0, 11)
    }
    else {
      istart <- early_dev - 1
      est_early_rec_dev <- rep(c(1, 0), c(istart, n_age - 
                                            istart - 1))
      est_rec_dev <- rep(1, n_y)
    }
  }
  if (tv_M == "DD") 
    est_early_rec_dev <- rep(0, n_age - 1)
  if (is.character(late_dev) && late_dev == "comp50") {
    if (any(comp == "age")) {
      comp_ldev <- colSums(CAA_hist, na.rm = TRUE)/max(colSums(CAA_hist, 
                                                               na.rm = TRUE))
    }
    else {
      comp_ldev <- colSums(CAL_hist, na.rm = TRUE)/max(colSums(CAL_hist, 
                                                               na.rm = TRUE))
    }
    comp_mode <- which.max(comp_ldev)[1]
    comp50_ind <- which(comp_ldev[1:comp_mode] <= 0.5)
    comp50_ind <- comp50_ind[length(comp50_ind)]
    if (all(comp != "age")) 
      comp50_ind <- ceiling(LinInterp(La, 0:n_age, Data@CAL_mids[comp50_ind]))
    late_dev <- ifelse(is.na(comp50_ind), 0, comp50_ind)
  }
  if (is.numeric(late_dev) && late_dev > 0) {
    if (late_dev > length(est_rec_dev)) 
      late_dev <- length(est_rec_dev)
    ind_late <- (length(est_rec_dev) - late_dev + 1):length(est_rec_dev)
    est_rec_dev[ind_late] <- ifelse(SR == "none", max(est_rec_dev[-ind_late]), 
                                    0)
  }
  if (rescale == "mean1") 
    rescale <- 1/mean(C_hist)
  Ind <- lapply(AddInd, Assess_I_hist, Data = Data, x = x, 
                yind = yind)
  I_hist <- vapply(Ind, getElement, numeric(n_y), "I_hist")
  I_sd <- vapply(Ind, getElement, numeric(n_y), "I_sd") %>% 
    pmax(0.05)
  I_units <- vapply(Ind, getElement, numeric(1), "I_units")
  I_vul <- vapply(AddInd, function(xx) {
    if (xx == "B") {
      return(rep(1, n_age))
    }
    else if (xx == "SSB") {
      return(mat_age)
    }
    else if (xx == "VB") {
      return(rep(0, n_age))
    }
    else {
      return(Data@AddIndV[x, suppressWarnings(as.numeric(xx)), 
                          1:n_age])
    }
  }, numeric(n_age))
  nsurvey <- ncol(I_hist)
  if (!is.list(LWT)) {
    if (!is.null(LWT) && length(LWT) != nsurvey) 
      stop("LWT needs to be a vector of length ", nsurvey)
    LWT <- list(Index = LWT)
    LWT$CAA <- LWT$CAL <- LWT$Catch <- 1
  }
  else {
    if (is.null(LWT$Index)) 
      LWT$Index <- rep(1, nsurvey)
    if (is.null(LWT$CAA)) 
      LWT$CAA <- 1
    if (is.null(LWT$CAL)) 
      LWT$CAL <- 1
    if (is.null(LWT$Catch)) 
      LWT$Catch <- 1
  }
  prior <- make_prior(prior, nsurvey, ifelse(SR == "BH", 1, 
                                             2), msg = FALSE)
  if (is.null(M_bounds)) {
    if (tv_M == "none") {
      M_bounds <- c(0, 10000)
    }
    else {
      M_bounds <- c(0.75, 1.25) * mean(M)
    }
  }
  data <- list(model = "SCA", C_hist = C_hist, rescale = rescale, 
               I_hist = I_hist, I_sd = I_sd, I_units = I_units, I_vul = I_vul, 
               abs_I = rep(0, nsurvey), nsurvey = nsurvey, CAA_hist = apply(CAA_hist, 
                                                                            1, tiny_comp) %>% t(), CAA_n = CAA_n_rescale, CAL_hist = apply(CAL_hist, 
                                                                                                                                           1, tiny_comp) %>% t(), CAL_n = CAL_n_rescale, LWT = c(LWT$Index, 
                                                                                                                                                                                                 LWT$CAA, LWT$CAL, LWT$Catch), n_y = n_y, n_age = n_age, 
               n_bin = ncol(PLA), weight = Wa, PLA = PLA, mat = mat_age, 
               vul_type = vulnerability, SR_type = SR, comp_dist = comp_dist, 
               catch_eq = catch_eq, est_early_rec_dev = est_early_rec_dev, 
               est_rec_dev = est_rec_dev, yindF = as.integer(0.5 * 
                                                               n_y), tv_M = tv_M, M_bounds = M_bounds, use_prior = prior$use_prior, 
               prior_dist = prior$pr_matrix)
  if (data$n_bin == 1) 
    data$CAL_hist <- t(data$CAL_hist)
  params <- list()
  if (!is.null(start)) {
    if (!is.null(start$R0) && is.numeric(start$R0)) 
      params$R0x <- log(start$R0[1] * rescale)
    if (!is.null(start$h) && is.numeric(start$h)) {
      if (SR == "BH") {
        h_start <- (start$h[1] - 0.2)/0.8
        params$transformed_h <- logit(h_start)
      }
      else if (SR == "Ricker") {
        params$transformed_h <- log(start$h[1] - 0.2)
      }
    }
    if (!is.null(start$M) && is.numeric(start$M)) 
      params$log_M0 <- log(start$M)
    if (catch_eq == "Baranov" && !is.null(start$F_equilibrium) && 
        is.numeric(start$F_equilibrium)) {
      params$F_equilibrium <- start$F_equilibrium
    }
    if (catch_eq == "Pope" && !is.null(start$U_equilibrium) && 
        is.numeric(start$U_equilibrium)) {
      params$F_equilibrium <- start$U_equilibrium
    }
    if (!is.null(start$vul_par) && is.numeric(start$vul_par)) {
      if (start$vul_par[1] > 0.75 * max_age) 
        stop("start$vul_par[1] needs to be less than 0.75 * Data@MaxAge (see help).")
      if (vulnerability == "logistic") {
        if (length(start$vul_par) < 2) 
          stop("Two parameters needed for start$vul_par with logistic vulnerability (see help).")
        if (start$vul_par[1] <= start$vul_par[2]) 
          stop("start$vul_par[1] needs to be greater than start$vul_par[2] (see help).")
        params$vul_par <- c(logit(start$vul_par[1]/max_age/0.75), 
                            log(start$vul_par[1] - start$vul_par[2]))
      }
      if (vulnerability == "dome") {
        if (length(start$vul_par) < 4) 
          stop("Four parameters needed for start$vul_par with dome vulnerability (see help).")
        if (start$vul_par[1] <= start$vul_par[2]) 
          stop("start$vul_par[1] needs to be greater than start$vul_par[2] (see help).")
        if (start$vul_par[3] <= start$vul_par[1] || 
            start$vul_par[3] >= max_age) {
          stop("start$vul_par[3] needs to be between start$vul_par[1] and Data@MaxAge (see help).")
        }
        if (start$vul_par[4] <= 0 || start$vul_par[4] >= 
            1) 
          stop("start$vul_par[4] needs to be between 0-1 (see help).")
        params$vul_par <- c(logit(start$vul_par[1]/max_age/0.75), 
                            log(start$vul_par[1] - start$vul_par[2]), 
                            logit(1/(max_age - start$vul_par[1])), logit(start$vul_par[4]))
      }
    }
    if (!is.null(start$F) && is.numeric(start$F)) {
      Fstart <- numeric(n_y)
      Fstart_ind <- data$yindF + 1
      Fstart[Fstart_ind] <- log(start$F[Fstart_ind])
      Fstart[-Fstart_ind] <- log(start$F[-Fstart_ind]/Fstart[Fstart_ind])
      params$log_F_dev <- Fstart
    }
    if (!is.null(start$omega) && is.numeric(start$omega)) 
      params$log_omega <- log(start$omega)
    if (!is.null(start[["tau"]]) && is.numeric(start[["tau"]])) 
      params$log_tau <- log(start[["tau"]])
    if (!is.null(start[["tau_M"]]) && is.numeric(start[["tau_M"]])) 
      params$log_tau_M <- log(start[["tau_M"]])
  }
  if (is.null(params$R0x)) {
    params$R0x <- ifelse(is.null(Data@OM$R0[x]), log(mean(data$C_hist)) + 
                           4, log(1.5 * rescale * Data@OM$R0[x]))
  }
  if (is.null(params$transformed_h)) {
    h_start <- ifelse(!fix_h && is.na(Data@steep[x]), 0.9, 
                      Data@steep[x])
    if (SR == "BH") {
      h_start <- (h_start - 0.2)/0.8
      params$transformed_h <- logit(h_start)
    }
    else if (SR == "Ricker") {
      params$transformed_h <- log(h_start - 0.2)
    }
    else {
      params$transformed_h <- 0
    }
  }
  if (is.null(params$log_M0)) 
    params$log_M0 <- log(M) %>% mean()
  if (is.null(params$logit_M_walk)) 
    params$logit_M_walk <- rep(0, n_y)
  if (is.null(params$F_equilibrium)) 
    params$F_equilibrium <- 0
  if (is.null(params$vul_par)) {
    if (any(comp == "age")) {
      comp_ldev <- colSums(CAA_hist, na.rm = TRUE)/max(colSums(CAA_hist, 
                                                               na.rm = TRUE))
    }
    else {
      comp_ldev <- colSums(CAL_hist, na.rm = TRUE)/max(colSums(CAL_hist, 
                                                               na.rm = TRUE))
    }
    comp_mode <- which.max(comp_ldev)[1]
    if (is.na(Data@LFC[x]) && is.na(Data@LFS[x]) || Data@LFC[x] > 
        Linf || Data@LFS[x] > Linf) {
      if (all(comp == "length")) {
        comp_mode <- ceiling(LinInterp(La, 0:max_age, 
                                       Data@CAL_mids[comp_mode]))
        if (!length(comp_mode)) 
          comp_mode <- 1
      }
      if (vulnerability == "logistic") 
        params$vul_par <- c(logit(comp_mode/max_age/0.75), 
                            log(1))
      if (vulnerability == "dome") {
        params$vul_par <- c(logit(comp_mode/max_age/0.75), 
                            log(1), logit(1/(max_age - comp_mode)), logit(0.5))
      }
    }
    else {
      A5 <- min(iVB(t0, K, Linf, Data@LFC[x]), comp_mode - 
                  1)
      Afull <- min(iVB(t0, K, Linf, Data@LFS[x]), 0.5 * 
                     max_age)
      A5 <- min(A5, Afull - 0.5)
      A50_vul <- mean(c(A5, Afull))
      if (vulnerability == "logistic") 
        params$vul_par <- c(logit(Afull/max_age/0.75), 
                            log(Afull - A50_vul))
      if (vulnerability == "dome") {
        params$vul_par <- c(logit(Afull/max_age/0.75), 
                            log(Afull - A50_vul), logit(0.1/(max_age - 
                                                               Afull)), logit(0.5))
      }
    }
  }
  if (is.na(params$vul_par[1])) 
    params$vul_par[1] <- 1
  if (is.null(params$log_F_dev)) {
    Fstart <- numeric(n_y)
    Fstart[data$yindF + 1] <- log(0.75 * mean(M))
    params$log_F_dev <- Fstart
  }
  if (is.null(params$log_omega)) {
    sigmaC <- max(0.01, sdconv(1, Data@CV_Cat[x]), na.rm = TRUE)
    params$log_omega <- log(sigmaC)
  }
  if (is.null(params[["log_tau"]])) {
    tau_start <- ifelse(is.na(Data@sigmaR[x]), 0.6, Data@sigmaR[x])
    params$log_tau <- log(tau_start)
  }
  if (is.null(params[["log_tau_M"]])) 
    params$log_tau_M <- log(0.05)
  params$log_early_rec_dev <- rep(0, n_age - 1)
  params$log_rec_dev <- rep(0, n_y)
  LH <- list(LAA = La, SD_LAA = SD_La, CAL_mids = CAL_mids, 
             WAA = Wa, Linf = Linf, K = K, t0 = t0, a = a, b = b, 
             A50 = A50, A95 = A95)
  info <- list(Year = Year, data = data, params = params, 
               LH = LH, control = control, inner.control = inner.control)
  map <- list()
  if (catch_eq == "Baranov" && any(info$data$C_hist <= 0)) {
    ind <- info$data$C_hist <= 0
    info$params$log_F_dev[ind] <- -20
    map_logF <- length(params$log_F_dev)
    map_logF[ind] <- NA
    map_logF[!ind] <- 1:sum(!ind)
    map$log_F_dev <- factor(map_logF)
  }
  else if (catch_eq == "Pope") {
    map$log_F_dev <- factor(rep(NA, n_y))
  }
  if (fix_h && !prior$use_prior[2]) 
    map$transformed_h <- factor(NA)
  if (!prior$use_prior[3]) 
    map$log_M0 <- factor(NA)
  if (tv_M != "walk") 
    map$logit_M_walk <- factor(rep(NA, n_y))
  if (fix_F_equilibrium) 
    map$F_equilibrium <- factor(NA)
  if (fix_omega) 
    map$log_omega <- factor(NA)
  if (fix_tau) 
    map$log_tau <- factor(NA)
  map$log_tau_M <- factor(NA)
  if (any(!est_early_rec_dev)) 
    map$log_early_rec_dev <- factor(ifelse(est_early_rec_dev, 
                                           1:sum(est_early_rec_dev), NA))
  if (any(!est_rec_dev)) 
    map$log_rec_dev <- factor(ifelse(est_rec_dev, 1:sum(est_rec_dev), 
                                     NA))
  if (vulnerability == "dome") 
    map$vul_par <- factor(c(1, 2, NA, 3))
  random <- NULL
  if (integrate) 
    random <- c("log_early_rec_dev", "log_rec_dev", "logit_M_walk")
  obj <- MakeADFun(data = info$data, parameters = info$params, 
                   hessian = TRUE, map = map, random = random, DLL = "SAMtool", 
                   inner.control = inner.control, silent = silent)
  if (catch_eq == "Pope") {
    high_U <- try(obj$report(c(obj$par, obj$env$last.par[obj$env$random]))$penalty > 
                    0, silent = TRUE)
    if (!is.character(high_U) && !is.na(high_U) && high_U) {
      Recruit <- try(Data@Rec[x, ], silent = TRUE)
      if (is.numeric(Recruit) && length(Recruit) == n_y && 
          any(!is.na(Recruit))) {
        log_rec_dev <- log(Recruit/mean(Recruit, na.rm = TRUE))
        log_rec_dev[is.na(est_rec_dev) | is.na(log_rec_dev) | 
                      is.infinite(log_rec_dev)] <- 0
        info$params$log_rec_dev <- log_rec_dev
        obj <- MakeADFun(data = info$data, parameters = info$params, 
                         hessian = TRUE, map = map, random = random, 
                         DLL = "SAMtool", inner.control = inner.control, 
                         silent = silent)
      }
      while (obj$par["R0x"] < 30 && obj$report(c(obj$par, 
                                                 obj$env$last.par[obj$env$random]))$penalty > 
             0) {
        obj$par["R0x"] <- obj$par["R0x"] + 1
      }
    }
  }
  mod <- optimize_TMB_model(obj, control, opt_hess, n_restart)
  opt <- mod[[1]]
  SD <- mod[[2]]
  report <- obj$report(obj$env$last.par.best)
  Yearplusone <- c(Year, max(Year) + 1)
  YearEarly <- (Year[1] - n_age + 1):(Year[1] - 1)
  YearDev <- c(YearEarly, Year)
  YearR <- c(YearDev, max(YearDev) + 1)
  R <- c(rev(report$R_early), report$R)
  Dev <- structure(c(rev(report$log_early_rec_dev), report$log_rec_dev), 
                   names = YearDev)
  report$dynamic_SSB0 <- SCA_dynamic_SSB0(obj, data = info$data, 
                                          params = info$params, map = map) %>% structure(names = Yearplusone)
  nll_report <- ifelse(is.character(opt), ifelse(integrate, 
                                                 NA, report$nll), opt$objective)
  Assessment <- new("Assessment", Model = "SCA", Name = Data@Name, 
                    conv = SD$pdHess, B0 = report$B0, R0 = report$R0, N0 = report$N0, 
                    SSB0 = report$E0, VB0 = report$VB0, B = structure(report$B, 
                                                                      names = Yearplusone), B_B0 = structure(report$B/report$B0, 
                                                                                                             names = Yearplusone), SSB = structure(report$E, 
                                                                                                                                                   names = Yearplusone), SSB_SSB0 = structure(report$E/report$E0, 
                                                                                                                                                                                              names = Yearplusone), VB = structure(report$VB, 
                                                                                                                                                                                                                                   names = Yearplusone), VB_VB0 = structure(report$VB/report$VB0, 
                                                                                                                                                                                                                                                                            names = Yearplusone), R = structure(R, names = YearR), 
                    N = structure(rowSums(report$N), names = Yearplusone), 
                    N_at_age = report$N, Selectivity = matrix(report$vul, 
                                                              nrow = length(Year), ncol = n_age, byrow = TRUE), 
                    Obs_Catch = structure(C_hist, names = Year), Obs_Index = structure(I_hist, 
                                                                                       dimnames = list(Year, paste0("Index_", 1:nsurvey))), 
                    Obs_C_at_age = CAA_hist, Catch = structure(report$Cpred, 
                                                               names = Year), Index = structure(report$Ipred, dimnames = list(Year, 
                                                                                                                              paste0("Index_", 1:nsurvey))), C_at_age = report$CAApred, 
                    Dev = Dev, Dev_type = "log-Recruitment deviations", 
                    NLL = structure(c(nll_report, report$nll_comp, report$prior, 
                                      report$penalty), names = c("Total", paste0("Index_", 
                                                                                 1:nsurvey), "CAA", "CAL", "Catch", "Dev", "M_dev", 
                                                                 "Prior", "Penalty")), info = info, obj = obj, opt = opt, 
                    SD = SD, TMB_report = report, dependencies = dependencies)
  if (tv_M != "walk") 
    Assessment@NLL <- Assessment@NLL[names(Assessment@NLL) != 
                                       "M_dev"]
  if (all(comp != "age")) 
    Assessment@NLL <- Assessment@NLL[names(Assessment@NLL) != 
                                       "CAA"]
  if (all(comp != "length")) 
    Assessment@NLL <- Assessment@NLL[names(Assessment@NLL) != 
                                       "CAL"]
  if (catch_eq == "Pope") 
    Assessment@NLL <- Assessment@NLL[names(Assessment@NLL) != 
                                       "Catch"]
  if (SR != "none") 
    Assessment@h <- report$h
  if (catch_eq == "Baranov") {
    Assessment@FMort <- structure(report$F, names = Year)
  }
  else {
    Assessment@U <- structure(report$U, names = Year)
  }
  if (Assessment@conv) {
    SE_Early <- as.list(SD, "Std. Error")$log_early_rec_dev %>% 
      rev()
    SE_Main <- as.list(SD, "Std. Error")$log_rec_dev
    SE_Dev <- structure(c(SE_Early, SE_Main), names = YearDev)
    if (any(is.na(SE_Dev))) {
      Dev <- Dev[seq(which(!is.na(SE_Dev))[1], length(YearDev))]
      SE_Dev <- SE_Dev[seq(which(!is.na(SE_Dev))[1], length(YearDev))]
      SE_Dev[is.na(SE_Dev)] <- 0
    }
    ref_pt <- lapply(seq_len(ifelse(tv_M == "walk", n_y, 
                                    1)), ref_pt_SCA, obj = obj, report = report)
    if (catch_eq == "Baranov") {
      report$FMSY <- vapply(ref_pt, getElement, numeric(1), 
                            "FMSY")
    }
    else {
      report$UMSY <- vapply(ref_pt, getElement, numeric(1), 
                            "UMSY")
    }
    report$MSY <- vapply(ref_pt, getElement, numeric(1), 
                         "MSY")
    report$VBMSY <- vapply(ref_pt, getElement, numeric(1), 
                           "VBMSY")
    report$RMSY <- vapply(ref_pt, getElement, numeric(1), 
                          "RMSY")
    report$BMSY <- vapply(ref_pt, getElement, numeric(1), 
                          "BMSY")
    report$EMSY <- vapply(ref_pt, getElement, numeric(1), 
                          "EMSY")
    refyear <- eval(refyear)
    if (catch_eq == "Baranov") {
      Assessment@FMSY <- report$FMSY[refyear]
      Assessment@F_FMSY <- structure(report$F/report$FMSY[refyear], 
                                     names = Year)
    }
    else {
      Assessment@UMSY <- report$UMSY[refyear]
      Assessment@U_UMSY <- structure(report$U/report$UMSY[refyear], 
                                     names = Year)
    }
    Assessment@MSY <- report$MSY[refyear]
    Assessment@BMSY <- report$BMSY[refyear]
    Assessment@SSBMSY <- report$EMSY[refyear]
    Assessment@VBMSY <- report$VBMSY[refyear]
    Assessment@B_BMSY <- structure(report$B/report$BMSY[refyear], 
                                   names = Yearplusone)
    Assessment@SSB_SSBMSY <- structure(report$E/report$EMSY[refyear], 
                                       names = Yearplusone)
    Assessment@VB_VBMSY <- structure(report$VB/report$VBMSY[refyear], 
                                     names = Yearplusone)
    Assessment@Dev <- Dev
    Assessment@SE_Dev <- SE_Dev
    Assessment@TMB_report <- report
    catch_eq_fn <- function(Ftarget) {
      projection_SCA(Assessment, Ftarget = Ftarget, p_years = 1, 
                     p_sim = 1, obs_error = list(array(1, c(1, 1, 
                                                            nsurvey)), matrix(1, 1, 1)), process_error = matrix(1, 
                                                                                                                1, 1)) %>% slot("Catch") %>% as.vector()
    }
    Assessment@forecast <- list(per_recruit = ref_pt[[refyear]][[7]], 
                                catch_eq = catch_eq_fn)
  }
  return(Assessment)
}
