Project_NK <- function (Hist = NULL, MPs = NA, parallel = FALSE, silent = FALSE, 
                        extended = FALSE, checkMPs = TRUE) 
{
  popdynOneTScpp <- MSEtool:::popdynOneTScpp
  calcRecruitment <- MSEtool:::calcRecruitment
  tiny <- MSEtool:::tiny
  CalcMPDynamics <- MSEtool:::CalcMPDynamics
  updateData <- MSEtool:::updateData
  CheckMPs <- MSEtool:::CheckMPs
  # readline("press any key to continue")
  if (class(Hist) != "Hist") 
    stop("Must provide an object of class `Hist`")
  OM <- Hist@OM
  set.seed(OM@seed)
  nsim <- OM@nsim
  nyears <- OM@nyears
  proyears <- OM@proyears
  interval <- OM@interval
  maxF <- OM@maxF
  pstar <- OM@pstar
  reps <- OM@reps
  control <- OM@cpars$control
  OM@cpars$control <- NULL
  if (checkMPs) 
    MPs <- CheckMPs(MPs = MPs, silent = silent)
  nMP <- length(MPs)
  if (nMP < 1) 
    stop("No valid MPs found", call. = FALSE)
  isrunning <- snowfall::sfIsRunning()
  if (!parallel & isrunning) 
    snowfall::sfStop()
  if (parallel) {
    if (!isrunning) 
      setup()
    Export_customMPs(MPs)
  }
  if (length(interval) != nMP) 
    interval <- rep(interval, nMP)[1:nMP]
  if (!all(interval == interval[1])) {
    if (!silent) 
      message("Variable management intervals:")
    df <- data.frame(MP = MPs, interval = interval)
    for (i in 1:nrow(df)) {
      message(df$MP[i], "has management interval:", df$interval[i])
    }
  }
  MSY_y <- Hist@Ref$ByYear$MSY
  FMSY_y <- Hist@Ref$ByYear$FMSY
  SSBMSY_y <- Hist@Ref$ByYear$SSBMSY
  BMSY_y <- Hist@Ref$ByYear$BMSY
  VBMSY_y <- Hist@Ref$ByYear$VBMSY
  F01_YPR_y <- Hist@Ref$ByYear$F01_YPR
  Fmax_YPR_y <- Hist@Ref$ByYear$Fmax_YPR
  F_SPR_y <- Hist@Ref$ByYear$F_SPR
  SPR_target <- F_SPR_y %>% dimnames() %>% getElement(2) %>% 
    substr(3, 4) %>% as.numeric()
  SPR_target <- SPR_target/100
  MSY_y <- array(MSY_y, dim = c(nsim, nyears + proyears, nMP)) %>% 
    aperm(c(1, 3, 2))
  FMSY_y <- array(FMSY_y, dim = c(nsim, nyears + proyears, 
                                  nMP)) %>% aperm(c(1, 3, 2))
  SSBMSY_y <- array(SSBMSY_y, dim = c(nsim, nyears + proyears, 
                                      nMP)) %>% aperm(c(1, 3, 2))
  BMSY_y <- array(BMSY_y, dim = c(nsim, nyears + proyears, 
                                  nMP)) %>% aperm(c(1, 3, 2))
  VBMSY_y <- array(VBMSY_y, dim = c(nsim, nyears + proyears, 
                                    nMP)) %>% aperm(c(1, 3, 2))
  F01_YPR_y <- array(F01_YPR_y, dim = c(nsim, nyears + proyears, 
                                        nMP)) %>% aperm(c(1, 3, 2))
  Fmax_YPR_y <- array(Fmax_YPR_y, dim = c(nsim, nyears + proyears, 
                                          nMP)) %>% aperm(c(1, 3, 2))
  F_SPR_y <- array(F_SPR_y, dim = c(nsim, length(SPR_target), 
                                    nyears + proyears, nMP)) %>% aperm(c(1, 4, 2, 3))
  Data <- Hist@Data
  Data_Misc <- Data@Misc
  Data@Misc <- list()
  MSElist <- list(Data)[rep(1, nMP)]
  SB_SBMSY_a <- array(NA, dim = c(nsim, nMP, proyears))
  F_FMSYa <- array(NA, dim = c(nsim, nMP, proyears))
  Ba <- array(NA, dim = c(nsim, nMP, proyears))
  SSBa <- array(NA, dim = c(nsim, nMP, proyears))
  VBa <- array(NA, dim = c(nsim, nMP, proyears))
  FMa <- array(NA, dim = c(nsim, nMP, proyears))
  Ca <- array(NA, dim = c(nsim, nMP, proyears))
  CaRet <- array(NA, dim = c(nsim, nMP, proyears))
  TACa <- array(NA, dim = c(nsim, nMP, proyears))
  Effort <- array(NA, dim = c(nsim, nMP, proyears))
  SPReqa <- array(NA, dim = c(nsim, nMP, proyears))
  SPRdyna <- array(NA, dim = c(nsim, nMP, proyears))
  Cost_out <- array(NA, dim = c(nsim, nMP, proyears))
  Rev_out <- array(NA, dim = c(nsim, nMP, proyears))
  LatEffort_out <- array(NA, dim = c(nsim, nMP, proyears))
  TAE_out <- array(NA, dim = c(nsim, nMP, proyears))
  StockPars <- Hist@SampPars$Stock
  FleetPars <- Hist@SampPars$Fleet
  ObsPars <- Hist@SampPars$Obs
  ImpPars <- Hist@SampPars$Imp
  StockPars$N <- Hist@AtAge$Number
  StockPars$Z <- Hist@AtAge$Z.Mortality
  StockPars$FM <- Hist@AtAge$F.Mortality
  StockPars$FMret <- Hist@AtAge$Fret.Mortality
  StockPars$CB <- Hist@AtAge$Removals
  StockPars$CBret <- Hist@AtAge$Landings
  StockPars$Biomass <- Hist@AtAge$Biomass
  StockPars$SSB <- Hist@AtAge$SBiomass
  StockPars$VBiomass <- Hist@AtAge$VBiomass
  n_age <- StockPars$n_age
  nareas <- StockPars$nareas
  RealData <- Hist@OM@cpars$Data
  ReferencePoints <- Hist@Ref$ReferencePoints
  LatentEff <- Hist@Misc$BioEco$LatentEff
  RevCurr <- Hist@Misc$BioEco$RevCurr
  CostCurr <- Hist@Misc$BioEco$CostCurr
  Response <- Hist@Misc$BioEco$Response
  CostInc <- Hist@Misc$BioEco$CostInc
  RevInc <- Hist@Misc$BioEco$RevInc
  N_P_mp <- array(NA, dim = c(nsim, n_age, nMP, proyears, 
                              nareas))
  B_P_mp <- array(NA, dim = c(nsim, n_age, nMP, proyears, 
                              nareas))
  SB_P_mp <- array(NA, dim = c(nsim, n_age, nMP, proyears, 
                               nareas))
  VB_P_mp <- array(NA, dim = c(nsim, n_age, nMP, proyears, 
                               nareas))
  Catch_P_mp <- array(NA, dim = c(nsim, n_age, nMP, proyears, 
                                  nareas))
  Removals_P_mp <- array(NA, dim = c(nsim, n_age, nMP, proyears, 
                                     nareas))
  FM_P_mp <- array(NA, dim = c(nsim, n_age, nMP, proyears, 
                               nareas))
  FMret_P_mp <- array(NA, dim = c(nsim, n_age, nMP, proyears, 
                                  nareas))
  mm <- 1
  for (mm in 1:nMP) {
    Data_MP <- MSElist[[mm]]
    Data_MP@Misc <- Data_Misc
    if (!silent) 
      message(mm, "/", nMP, " Running MSE for ", MPs[mm])
    checkNA <- rep(0, OM@proyears)
    upyrs <- seq(from = 1, to = proyears, by = interval[mm])
    L5_P <- FleetPars$L5_y
    LFS_P <- FleetPars$LFS_y
    Vmaxlen_P <- FleetPars$Vmaxlen_y
    SLarray_P <- FleetPars$SLarray_real
    V_P <- FleetPars$V_real
    LR5_P <- FleetPars$LR5_y
    LFR_P <- FleetPars$LFR_y
    Rmaxlen_P <- FleetPars$Rmaxlen_y
    retA_P <- FleetPars$retA_real
    retL_P <- FleetPars$retL_real
    Fdisc_P <- FleetPars$Fdisc_array1
    DR_P <- FleetPars$DR_y
    LatentEff_MP <- LatentEff
    N_P <- array(NA, dim = c(nsim, n_age, proyears, nareas))
    Biomass_P <- array(NA, dim = c(nsim, n_age, proyears, 
                                   nareas))
    VBiomass_P <- array(NA, dim = c(nsim, n_age, proyears, 
                                    nareas))
    SSN_P <- array(NA, dim = c(nsim, n_age, proyears, nareas))
    SSB_P <- array(NA, dim = c(nsim, n_age, proyears, nareas))
    FM_P <- array(NA, dim = c(nsim, n_age, proyears, nareas))
    FM_Pret <- array(NA, dim = c(nsim, n_age, proyears, 
                                 nareas))
    Z_P <- array(NA, dim = c(nsim, n_age, proyears, nareas))
    CB_P <- array(NA, dim = c(nsim, n_age, proyears, nareas))
    CB_Pret <- array(NA, dim = c(nsim, n_age, proyears, 
                                 nareas))
    SAYRL <- as.matrix(expand.grid(1:nsim, 1:n_age, nyears, 
                                   1:nareas))
    SAYRt <- as.matrix(expand.grid(1:nsim, 1:n_age, 1 + 
                                     nyears, 1:nareas))
    SAYR <- as.matrix(expand.grid(1:nsim, 1:n_age, 1, 1:nareas))
    SYt <- SAYRt[, c(1, 3)]
    SAYt <- SAYRt[, 1:3]
    SR <- SAYR[, c(1, 4)]
    SA1 <- SAYR[, 1:2]
    S1 <- SAYR[, 1]
    SY1 <- SAYR[, c(1, 3)]
    SAY1 <- SAYRt[, 1:3]
    SYA <- as.matrix(expand.grid(1:nsim, 1, 1:n_age))
    SY <- SYA[, 1:2]
    SA <- SYA[, c(1, 3)]
    SAY <- SYA[, c(1, 3, 2)]
    S <- SYA[, 1]
    y <- 1
    if (!silent) {
      pb <- txtProgressBar(min = 1, max = proyears, style = 3, 
                           width = min(getOption("width"), 50))
    }
    NextYrN <- lapply(1:nsim, function(x) popdynOneTScpp(nareas, 
                                                         StockPars$maxage, Ncurr = StockPars$N[x, , nyears, 
                                                         ], Zcurr = StockPars$Z[x, , nyears, ], plusgroup = StockPars$plusgroup))
    N_P[, , 1, ] <- aperm(array(unlist(NextYrN), dim = c(n_age, 
                                                         nareas, nsim, 1)), c(3, 1, 4, 2))
    Biomass_P[SAYR] <- N_P[SAYR] * StockPars$Wt_age[SAY1]
    VBiomass_P[SAYR] <- Biomass_P[SAYR] * V_P[SAYt]
    SSN_P[SAYR] <- N_P[SAYR] * StockPars$Mat_age[SAY1]
    SSB_P[SAYR] <- N_P[SAYR] * StockPars$Fec_Age[SAY1]
    SSBcurr <- apply(SSB_P[, , 1, ], c(1, 3), sum)
    recdev <- StockPars$Perr_y[, nyears + n_age]
    rec_area <- sapply(1:nsim, calcRecruitment, SRrel = StockPars$SRrel, 
                       SSBcurr = SSBcurr, recdev = recdev, hs = StockPars$hs, 
                       aR = StockPars$aR, bR = StockPars$bR, R0a = StockPars$R0a, 
                       SSBpR = StockPars$SSBpR, SSB0 = StockPars$SSB0)
    N_P[, 1, y, ] <- t(rec_area)
    Ntemp <- lapply(1:nsim, function(x) movestockCPP(nareas, 
                                                     StockPars$maxage, StockPars$mov[x, , , , y + nyears], 
                                                     N_P[x, , y, ]))
    N_P[, , y, ] <- array(unlist(Ntemp), dim = c(n_age, 
                                                 nareas, nsim)) %>% aperm(c(3, 1, 2))
    Biomass_P[SAYR] <- N_P[SAYR] * StockPars$Wt_age[SAY1]
    VBiomass_P[SAYR] <- Biomass_P[SAYR] * V_P[SAYt]
    SSN_P[SAYR] <- N_P[SAYR] * StockPars$Mat_age[SAY1]
    SSB_P[SAYR] <- N_P[SAYR] * StockPars$Fec_Age[SAY1]
    StockPars$N_P <- N_P
    runMP <- applyMP_NK(Data = Data_MP, MPs = MPs[mm], reps = reps, 
                     silent = TRUE)
    MPRecs <- runMP[[1]][[1]]
    Data_p <- runMP[[2]]
    Data_p@TAC <- MPRecs$TAC
    LastSpatial <- array(FleetPars$MPA[nyears, ], dim = c(nareas, 
                                                          nsim))
    LastAllocat <- rep(1, nsim)
    LastTAC <- LastCatch <- apply(StockPars$CBret[, , nyears, 
    ], 1, sum)
    TACused <- apply(Data_p@TAC, 2, quantile, p = pstar, 
                     na.rm = T)
    if (length(MPRecs$TAC) > 0) {
      checkNA[y] <- sum(is.na(TACused))
      TACused[is.na(TACused)] <- LastTAC[is.na(TACused)]
      TACused[TACused < tiny] <- tiny
      TACa[, mm, y] <- TACused
    }
    RevPC <- RevCurr/LastCatch
    PMargin <- 1 - CostCurr/(RevPC * LastCatch)
    Profit <- (RevPC * LastCatch) - CostCurr
    HistEffort <- rep(1, nsim)
    Effort_pot <- HistEffort + Response * Profit
    Effort_pot[Effort_pot < 0] <- tiny
    if (!all(is.na(LatentEff_MP))) {
      LastTAE <- histTAE <- HistEffort/(1 - LatentEff_MP)
    }
    else {
      LastTAE <- histTAE <- rep(NA, nsim)
    }
    MPCalcs <- CalcMPDynamics(MPRecs, y, nyears, proyears, 
                              nsim, Biomass_P, VBiomass_P, LastTAE, histTAE, LastSpatial, 
                              LastAllocat, LastTAC, TACused, maxF, LR5_P, LFR_P, 
                              Rmaxlen_P, retL_P, retA_P, L5_P, LFS_P, Vmaxlen_P, 
                              SLarray_P, V_P, Fdisc_P, DR_P, FM_P, FM_Pret, Z_P, 
                              CB_P, CB_Pret, Effort_pot, StockPars, FleetPars, 
                              ImpPars, control = control)
    TACa[, mm, y] <- MPCalcs$TACrec
    LastSpatial <- MPCalcs$Si
    LastAllocat <- MPCalcs$Ai
    LastTAE <- MPCalcs$TAE
    LastTAC <- MPCalcs$TACrec
    Effort[, mm, y] <- MPCalcs$Effort
    CB_P <- MPCalcs$CB_P
    CB_Pret <- MPCalcs$CB_Pret
    FM_P <- MPCalcs$FM_P
    FM_Pret <- MPCalcs$FM_Pret
    Z_P <- MPCalcs$Z_P
    retA_P <- MPCalcs$retA_P
    retL_P <- MPCalcs$retL_P
    V_P <- MPCalcs$V_P
    SLarray_P <- MPCalcs$SLarray_P
    FMa[, mm, y] <- MPCalcs$Ftot
    LR5_P <- MPCalcs$LR5_P
    LFR_P <- MPCalcs$LFR_P
    Rmaxlen_P <- MPCalcs$Rmaxlen_P
    L5_P <- MPCalcs$L5_P
    LFS_P <- MPCalcs$LFS_P
    Vmaxlen_P <- MPCalcs$Vmaxlen_P
    Fdisc_P <- MPCalcs$Fdisc_P
    DR_P <- MPCalcs$DR_P
    RetainCatch <- apply(CB_Pret[, , y, ], 1, sum)
    RetainCatch[RetainCatch <= 0] <- tiny
    Cost_out[, mm, y] <- Effort[, mm, y] * CostCurr * (1 + 
                                                         CostInc/100)^y
    Rev_out[, mm, y] <- (RevPC * (1 + RevInc/100)^y * RetainCatch)
    PMargin <- 1 - Cost_out[, mm, y]/Rev_out[, mm, y]
    Profit <- Rev_out[, mm, y] - Cost_out[, mm, y]
    Effort_pot <- Effort_pot + Response * Profit
    Effort_pot[Effort_pot < 0] <- tiny
    LatEffort_out[, mm, y] <- LastTAE - Effort[, mm, y]
    TAE_out[, mm, y] <- LastTAE
    for (y in 2:proyears) {
      if (!silent) {
        setTxtProgressBar(pb, y)
      }
      SelectChanged <- FALSE
      if (any(range(retA_P[, , nyears + y] - FleetPars$retA_real[, 
                                                                 , nyears + y]) != 0)) 
        SelectChanged <- TRUE
      if (any(range(V_P[, , nyears + y] - FleetPars$V_real[, 
                                                           , nyears + y]) != 0)) 
        SelectChanged <- TRUE
      if (SelectChanged) {
        y1 <- nyears + y
        MSYrefsYr <- sapply(1:nsim, optMSY_eq, StockPars$M_ageArray, 
                            StockPars$Wt_age, StockPars$Mat_age, Fec_age = StockPars$Fec_Age, 
                            V_P, StockPars$maxage, StockPars$R0, StockPars$SRrel, 
                            StockPars$SSBpR, StockPars$hs, yr.ind = y1, 
                            plusgroup = StockPars$plusgroup)
        MSY_y[, mm, y1] <- MSYrefsYr[1, ]
        FMSY_y[, mm, y1] <- MSYrefsYr[2, ]
        SSBMSY_y[, mm, y1] <- MSYrefsYr[3, ]
        per_recruit_F <- lapply(1:nsim, per_recruit_F_calc, 
                                M_ageArray = StockPars$M_ageArray, Wt_age = StockPars$Wt_age, 
                                Mat_age = StockPars$Mat_age, Fec_age = StockPars$Fec_Age, 
                                V = FleetPars$V_real, maxage = StockPars$maxage, 
                                yr.ind = y1, plusgroup = StockPars$plusgroup, 
                                SPR_target = SPR_target, StockPars = StockPars)
        F_SPR_y[, mm, , y1] <- sapply(per_recruit_F, 
                                      getElement, 1) %>% t()
        F01_YPR_y[, mm, y1] <- sapply(per_recruit_F, 
                                      function(x) x[[2]][1])
        Fmax_YPR_y[, mm, y1] <- sapply(per_recruit_F, 
                                       function(x) x[[2]][2])
      }
      SAYRt <- as.matrix(expand.grid(1:nsim, 1:n_age, 
                                     y + nyears, 1:nareas))
      SAYt <- SAYRt[, 1:3]
      SAYtMP <- cbind(SAYt, mm)
      SYt <- SAYRt[, c(1, 3)]
      SAY1R <- as.matrix(expand.grid(1:nsim, 1:n_age, 
                                     y - 1, 1:nareas))
      SAYR <- as.matrix(expand.grid(1:nsim, 1:n_age, y, 
                                    1:nareas))
      SY <- SAYR[, c(1, 3)]
      SA <- SAYR[, 1:2]
      S1 <- SAYR[, 1]
      SAY <- SAYR[, 1:3]
      S <- SAYR[, 1]
      SR <- SAYR[, c(1, 4)]
      SA2YR <- as.matrix(expand.grid(1:nsim, 2:n_age, 
                                     y, 1:nareas))
      SA1YR <- as.matrix(expand.grid(1:nsim, 1:(n_age - 
                                                  1), y - 1, 1:nareas))
      NextYrN <- lapply(1:nsim, function(x) popdynOneTScpp(nareas, 
                                                           StockPars$maxage, Ncurr = N_P[x, , y - 1, ], 
                                                           Zcurr = Z_P[x, , y - 1, ], plusgroup = StockPars$plusgroup))
      N_P[, , y, ] <- aperm(array(unlist(NextYrN), dim = c(n_age, 
                                                           nareas, nsim, 1)), c(3, 1, 4, 2))
      Biomass_P[SAYR] <- N_P[SAYR] * StockPars$Wt_age[SAYt]
      SSN_P[SAYR] <- N_P[SAYR] * StockPars$Mat_age[SAYt]
      SSB_P[SAYR] <- N_P[SAYR] * StockPars$Fec_Age[SAYt]
      SSBcurr <- apply(SSB_P[, , y, ], c(1, 3), sum)
      recdev <- StockPars$Perr_y[, y + nyears + n_age - 
                                   1]
      rec_area <- sapply(1:nsim, calcRecruitment, SRrel = StockPars$SRrel, 
                         SSBcurr = SSBcurr, recdev = recdev, hs = StockPars$hs, 
                         aR = StockPars$aR, bR = StockPars$bR, R0a = StockPars$R0a, 
                         SSBpR = StockPars$SSBpR, SSB0 = StockPars$SSB0)
      N_P[, 1, y, ] <- t(rec_area)
      Ntemp <- lapply(1:nsim, function(x) movestockCPP(nareas, 
                                                       StockPars$maxage, StockPars$mov[x, , , , y + 
                                                                                         nyears], N_P[x, , y, ]))
      N_P[, , y, ] <- array(unlist(Ntemp), dim = c(n_age, 
                                                   nareas, nsim)) %>% aperm(c(3, 1, 2))
      Biomass_P[SAYR] <- N_P[SAYR] * StockPars$Wt_age[SAY1]
      VBiomass_P[SAYR] <- Biomass_P[SAYR] * V_P[SAYt]
      SSN_P[SAYR] <- N_P[SAYR] * StockPars$Mat_age[SAYt]
      SSB_P[SAYR] <- N_P[SAYR] * StockPars$Fec_Age[SAYt]
      StockPars$N_P <- N_P
      myData_MP1 <<- Data_MP
      if (y %in% upyrs) {
        Data_MP <- updateData(Data = Data_MP, OM, MPCalcs, 
                              Effort, Biomass = StockPars$Biomass, StockPars$N, 
                              Biomass_P, CB_Pret, N_P, SSB = StockPars$SSB, 
                              SSB_P, VBiomass = StockPars$VBiomass, VBiomass_P, 
                              RefPoints = ReferencePoints, retA_P, retL_P, 
                              StockPars, FleetPars, ObsPars, ImpPars, V_P, 
                              upyrs, interval, y, mm, Misc = Data_p@Misc, 
                              RealData, Sample_Area = ObsPars$Sample_Area)
        myData_MP2 <<- Data_MP
        runMP <- applyMP_NK(Data = Data_MP, MPs = MPs[mm], 
                         reps = reps, silent = TRUE)
        MPRecs <- runMP[[1]][[1]]
        Data_p <- runMP[[2]]
        Data_p@TAC <- MPRecs$TAC
      }
      TACused <- apply(Data_p@TAC, 2, quantile, p = pstar, 
                       na.rm = T)
      if (length(MPRecs$TAC) > 0) {
        checkNA[y] <- sum(is.na(TACused))
        TACused[is.na(TACused)] <- LastTAC[is.na(TACused)]
        TACused[TACused < tiny] <- tiny
        TACa[, mm, y] <- TACused
      }
      MPCalcs <- CalcMPDynamics(MPRecs, y, nyears, proyears, 
                                nsim, Biomass_P, VBiomass_P, LastTAE, histTAE, 
                                LastSpatial, LastAllocat, LastTAC, TACused, 
                                maxF, LR5_P, LFR_P, Rmaxlen_P, retL_P, retA_P, 
                                L5_P, LFS_P, Vmaxlen_P, SLarray_P, V_P, Fdisc_P, 
                                DR_P, FM_P, FM_Pret, Z_P, CB_P, CB_Pret, Effort_pot, 
                                StockPars, FleetPars, ImpPars, control = control)
      LastSpatial <- MPCalcs$Si
      LastAllocat <- MPCalcs$Ai
      LastTAE <- MPCalcs$TAE
      Effort[, mm, y] <- MPCalcs$Effort
      FMa[, mm, y] <- MPCalcs$Ftot
      CB_P <- MPCalcs$CB_P
      CB_Pret <- MPCalcs$CB_Pret
      LastTAC <- TACa[, mm, y]
      FM_P <- MPCalcs$FM_P
      FM_Pret <- MPCalcs$FM_Pret
      Z_P <- MPCalcs$Z_P
      retA_P <- MPCalcs$retA_P
      retL_P <- MPCalcs$retL_P
      V_P <- MPCalcs$V_P
      SLarray_P <- MPCalcs$SLarray_P
      LR5_P <- MPCalcs$LR5_P
      LFR_P <- MPCalcs$LFR_P
      Rmaxlen_P <- MPCalcs$Rmaxlen_P
      L5_P <- MPCalcs$L5_P
      LFS_P <- MPCalcs$LFS_P
      Vmaxlen_P <- MPCalcs$Vmaxlen_P
      Fdisc_P <- MPCalcs$Fdisc_P
      DR_P <- MPCalcs$DR_P
      RetainCatch <- apply(CB_Pret[, , y, ], 1, sum)
      RetainCatch[RetainCatch <= 0] <- tiny
      Cost_out[, mm, y] <- Effort[, mm, y] * CostCurr * 
        (1 + CostInc/100)^y
      Rev_out[, mm, y] <- (RevPC * (1 + RevInc/100)^y * 
                             RetainCatch)
      Profit <- Rev_out[, mm, y] - Cost_out[, mm, y]
      Effort_pot <- Effort_pot + Response * Profit
      Effort_pot[Effort_pot < 0] <- tiny
      LatEffort_out[, mm, y] <- LastTAE - Effort[, mm, 
                                                 y]
      TAE_out[, mm, y] <- LastTAE
      if ("progress" %in% names(control)) {
        if (control$progress) {
          if (requireNamespace("shiny", quietly = TRUE)) {
            shiny::setProgress((mm - 1 + y/proyears)/nMP, 
                               detail = paste0(round((mm - 1 + y/proyears)/nMP * 
                                                       100), "% \n Management procedure ", 
                                               mm, "/", nMP, " (", MPs[mm], ")"))
          }
          else {
            warning("package `shiny` needs to be installed for progress bar")
          }
        }
      }
    }
    if (!silent) 
      close(pb)
    if (max(upyrs) < proyears) {
      Data_MP <- updateData(Data = Data_MP, OM, MPCalcs, 
                            Effort, StockPars$Biomass, StockPars$N, Biomass_P, 
                            CB_Pret, N_P, StockPars$SSB, SSB_P, StockPars$VBiomass, 
                            VBiomass_P, RefPoints = ReferencePoints, retA_P, 
                            retL_P, StockPars, FleetPars, ObsPars, ImpPars, 
                            V_P, upyrs = c(upyrs, proyears), interval = rep(proyears - 
                                                                              max(upyrs), length(interval)), y, mm, Misc = Data_p@Misc, 
                            RealData, ObsPars$Sample_Area)
    }
    SB_SBMSY_a[, mm, ] <- apply(SSB_P, c(1, 3), sum, na.rm = TRUE)/SSBMSY_y[, 
                                                                            mm, (OM@nyears + 1):(OM@nyears + OM@proyears)]
    F_FMSYa[, mm, ] <- FMa[, mm, ]/FMSY_y[, mm, (OM@nyears + 
                                                   1):(OM@nyears + OM@proyears)]
    Ba[, mm, ] <- apply(Biomass_P, c(1, 3), sum, na.rm = TRUE)
    SSBa[, mm, ] <- apply(SSB_P, c(1, 3), sum, na.rm = TRUE)
    VBa[, mm, ] <- apply(VBiomass_P, c(1, 3), sum, na.rm = TRUE)
    Ca[, mm, ] <- apply(CB_P, c(1, 3), sum, na.rm = TRUE)
    CaRet[, mm, ] <- apply(CB_Pret, c(1, 3), sum, na.rm = TRUE)
    SPReqa[, mm, ] <- CalcSPReq(FM_P, StockPars, n_age, 
                                nareas, nyears, proyears, nsim)
    SPRdyna[, mm, ] <- CalcSPRdyn(abind::abind(StockPars$FM, 
                                               FM_P, along = 3), StockPars, n_age, nareas, nyears, 
                                  proyears, nsim)
    if (!silent) {
      cat("\n")
      if (all(checkNA[upyrs] != nsim) & !all(checkNA == 
                                             0)) {
        ntot <- sum(checkNA[upyrs])
        totyrs <- sum(checkNA[upyrs] > 0)
        nfrac <- round(ntot/(length(upyrs) * nsim), 
                       2) * 100
        message(totyrs, " years had TAC = NA for some simulations (", 
                nfrac, "% of total simulations)")
        message("Used TAC_y = TAC_y-1")
      }
    }
    N_P_mp[, , mm, , ] <- N_P
    B_P_mp[, , mm, , ] <- Biomass_P
    SB_P_mp[, , mm, , ] <- SSB_P
    VB_P_mp[, , mm, , ] <- VBiomass_P
    Catch_P_mp[, , mm, , ] <- CB_Pret
    Removals_P_mp[, , mm, , ] <- CB_P
    FM_P_mp[, , mm, , ] <- FM_P
    FMret_P_mp[, , mm, , ] <- FM_Pret
    Data_MP@Misc$StockPars <- Data_MP@Misc$FleetPars <- Data_MP@Misc$ReferencePoints <- NULL
    MSElist[[mm]] <- Data_MP
  }
  CB_hist <- apply(Hist@TSdata$Landings, 1:2, sum)
  FM_hist <- Hist@TSdata$Find * Hist@SampPars$Fleet$qs
  SSB_hist <- apply(Hist@TSdata$SBiomass, 1:2, sum)
  Misc <- list()
  Misc$extended <- list()
  if (extended) {
    histN <- replicate(nMP, StockPars$N) %>% aperm(c(1, 
                                                     2, 5, 3, 4))
    N_all <- abind::abind(histN, N_P_mp, along = 4)
    histB <- replicate(nMP, StockPars$Biomass) %>% aperm(c(1, 
                                                           2, 5, 3, 4))
    B_all <- abind::abind(histB, B_P_mp, along = 4)
    histSSB <- replicate(nMP, StockPars$SSB) %>% aperm(c(1, 
                                                         2, 5, 3, 4))
    SB_all <- abind::abind(histSSB, SB_P_mp, along = 4)
    histVB <- replicate(nMP, StockPars$VBiomass) %>% aperm(c(1, 
                                                             2, 5, 3, 4))
    VB_all <- abind::abind(histVB, VB_P_mp, along = 4)
    histCatch <- replicate(nMP, StockPars$CBret) %>% aperm(c(1, 
                                                             2, 5, 3, 4))
    Catch_all <- abind::abind(histCatch, Catch_P_mp, along = 4)
    histRemovals <- replicate(nMP, StockPars$CB) %>% aperm(c(1, 
                                                             2, 5, 3, 4))
    Removals_all <- abind::abind(histRemovals, Removals_P_mp, 
                                 along = 4)
    histF <- replicate(nMP, StockPars$FM) %>% aperm(c(1, 
                                                      2, 5, 3, 4))
    FM_all <- abind::abind(histF, FM_P_mp, along = 4)
    histFret <- replicate(nMP, StockPars$FMret) %>% aperm(c(1, 
                                                            2, 5, 3, 4))
    FMret_all <- abind::abind(histFret, FMret_P_mp, along = 4)
    Misc$extended <- list(N = N_all, B = B_all, SSB = SB_all, 
                          VB = VB_all, Catch = Catch_all, Removals = Removals_all, 
                          FM = FM_all, FMret = FMret_all)
    Hist_out <- Hist
  }
  else {
    Hist_out <- new("Hist")
  }
  MSEout <- new("MSE", Name = OM@Name, nyears = nyears, proyears = proyears, 
                nMPs = nMP, MPs = MPs, nsim = nsim, OM = Data@OM, Obs = Data@Obs, 
                SB_SBMSY = SB_SBMSY_a, F_FMSY = F_FMSYa, N = apply(N_P_mp, 
                                                                   c(1, 3, 4), sum), B = Ba, SSB = SSBa, VB = VBa, 
                FM = FMa, SPR = list(Equilibrium = SPReqa, Dynamic = SPRdyna), 
                Catch = CaRet, Removals = Ca, Effort = Effort, TAC = TACa, 
                TAE = TAE_out, BioEco = list(LatEffort = LatEffort_out, 
                                             Revenue = Rev_out, Cost = Cost_out), RefPoint = list(MSY = MSY_y, 
                                                                                                  FMSY = FMSY_y, SSBMSY = SSBMSY_y, F_SPR = F_SPR_y, 
                                                                                                  Dynamic_Unfished = Hist@Ref$Dynamic_Unfished, ByYear = Hist@Ref$ByYear), 
                CB_hist = CB_hist, FM_hist = FM_hist, SSB_hist = SSB_hist, 
                Hist = Hist_out, PPD = MSElist, Misc = Misc)
  attr(MSEout, "version") <- packageVersion("MSEtool")
  attr(MSEout, "date") <- date()
  attr(MSEout, "R.version") <- R.version
  MSEout
}
