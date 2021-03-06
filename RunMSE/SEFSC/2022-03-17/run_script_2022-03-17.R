library(abind)
library(bamExtras)
library(bamMSE)
library(openMSE)
library(magrittr)

rm(list=ls())
t_list <- Sys.time()
myseed <- 8675309

source("fn/make_interim_MP_lag.R")
source("fn/make_projection_MP_lag.R")
source("fn/SCA_lag.R")

nsim <- 48 #250
runScenarios <- TRUE # Run scenarios to do MSE or just generate historical data?
runMSE_args <- list("parallel"=TRUE,"extended"=TRUE,"silent"=FALSE)
lag_Assess_Init <- 2 # Number of years between terminal year of assessment and first year of management. May be modified in scenarios
MSY_frac_Init <- 0.75 # Fraction of MSY for setting TAC. May be modified in scenarios
AddInd_val_Init <- 1  # Index used in interim procedures and possibly assessments. May be modified in scenarios
AddInd_all_assess_Init <- TRUE # Should all available indices be used in assessments? Only works with SCA_lag() May be modified in scenarios

OMName_scen_complete <- NA
#   local({
#   a <- gsub(".rds","",list.files("MSE_obj"))
#   b <- gsub("MSE","OM",a)
#   b
# })

OMName_all <- gsub(".rds","",list.files("OM"))
OMName <- OMName_all#[!OMName_all%in%OMName_complete]
OMName <- c("OM_BlackSeaBass", # Runs
             "OM_RedPorgy"#, # Runs
             ,"OM_VermilionSnapper" # Runs
             #,
             #"OM_SnowyGrouper" # Runs 2022-1-20 batch 1 for base took only 25 min
             #"OM_RedGrouper", # 2022-1-19 This took 7 hours just to run the base scenario batch 1 so I interrupted it
             #,"OM_GagGrouper"#, # Runs
             #,"OM_GrayTriggerfish" # Problems with lightly fished scenario where "More than 5 % of simulations can't get to the specified level of depletion with these Operating Model parameters"
             #"OM_RedSnapper" # Seems to run but takes a long time
)

# Setup loops
scenario <- c(#"base", #"hs", "hd",
              #"nolag"#,
               "noaddind"
              #"lf",
              #"dep",
              #"ucvhi",
              #"ucvlo",
              #"ubias"#,
              #"epiM",
              #"rc",     # Regime change
              #"perfobs", # Perfect observations
              #"minerr",    # Minimize error and variation in operating model
              #"minrecdev",  # Minimize recruitment deviations
              #"minrecac"#, # Minimize recruitment autocorrelation
              #"tachi",
              #"vtiv"#,    # Vulnerability time-invariant
              # "mconst"
              ##"genfle"   # Generic fleet (kicked out error in RedPorgy "'Len_age' must be array with dimensions: nsim, maxage+1, nyears + proyears.")
)

# Depletion scenario args
dep_args <- list("scale"=0.5, "min"=0.05)

# Lightly fished scenario args
lf_args <- list("scale"=2,"max"=1)

## Episodic M scenario args
# yrprop:  Proportion of years to apply a multiplier on M
# M_mult_max:  Limit on how high the multiplier on M can be (Huynh set it to 4 for age-invariant M
#     but when I used that, the TACs were absurdly high like 1e+11 causing other values to be absurd)
#     Values as high as 0.5 resulted in some MSY values equal to zero for RedPorgy though it was fine for
#     BlackSeaBass and VermilionSnapper. I think it has to do with how depleted Red Porgy is. Perhaps it is crashing the population
#     when M is too high.
# M_lognorm_sd: lognormal sd on distribution of M_mult. Huynh used 2 but that results in a lot of values above
#     M_mult_max, so when  you apply pmin to limit the maximum value, then M_mult_max
#     ends up being one of the most common M_mult values
epiM_args <- list("yrprop"=0.1,
                  "M_mult_max"=0.2,
                  "M_lognorm_sd" <- 0.2
                  )

# Index cv high scenario args
ucvhi_args <- list("scale"=2)

# Index cv low scenario args
ucvlo_args <- list("scale"=0.5)

# Index bias scenario args
ubias_args <- list("int"=0,
                     "slope"=-0.01, # Slope of change in index error per year
                     "yr1diff"=10  # Number of years between the beginning of the projection period and start of change in errors
                    )

# Regime change scenario args
rc_args <-  list("yr1diff"=  10,   # Number of years between the beginning of the projection period and start of change in rec devs
                 "transdur"= 10,   # Duration (in years) of transition between regime 1 and 2
                 "r2_mult" =  0.75  # Multiplier on rec devs for regime 2. (a value of 1 would mean recruitment was not changing)
)

# seeds <- setNames(sample(1:10000,length(OMName),replace = FALSE),OMName)

MSEtool::setup(12) # Run in parallel over 12 cores
  sfLibrary("magrittr", character.only = TRUE, verbose = FALSE)

for(OMName_k in OMName) { ######### Loop over operating model

  MSEName_k <- gsub("OM","MSE",OMName_k)
  DataName_k <- gsub("OM","Data",OMName_k)

  OMInit_k <- readRDS(paste0("OM/", OMName_k, ".rds"))
  DataInit_k <- readRDS(paste0("Data/", DataName_k, ".rds"))
  Data_k <- DataInit_k

    for(scenario_i in scenario) { ######### Loop over scenario

      set.seed(myseed)
      # All scenarios
      OMName_scen <- paste0(OMName_k, "_", scenario_i)
      if(!OMName_scen%in%OMName_scen_complete){
        OM_k <- OMInit_k

        # Set lower limit for recruitment autocorrelation
        OM_k@AC[1] <- OM_k@AC[1]*.25

        # Set M to be age-invariant in the OM as it is in the stock assessment
        OM_k@cpars <- OM_k@cpars[names(OM_k@cpars)!="M_ageArray"]
        OM_k@M <- rep(Data_k@Mort,2)

        ## SCENARIOS

        # Hyperstable
        if(scenario_i == "hs") OM_k@beta <- c(1/3, 2/3) # beta values below 1 lead to hyperstability

        # Hyperdeplete
        if(scenario_i == "hd") OM_k@beta <- c(1.5, 3)

        # Depleted
        if(scenario_i == "dep") {
          OM_k@cpars$D <- pmax(dep_args$scale * OM_k@cpars$D, dep_args$min)
          # Remove qs so that runMSE will Optimize for user-specified depletion in last historical year
          OM_k@cpars <- OM_k@cpars[names(OM_k@cpars)[names(OM_k@cpars)!="qs"]]
        }

        # Lightly fished
        if(scenario_i == "lf")  {
          OM_k@cpars$D <- pmin(lf_args$scale * OM_k@cpars$D, lf_args$max)
          # Remove qs so that runMSE will Optimize for user-specified depletion in last historical year
          OM_k@cpars <- OM_k@cpars[names(OM_k@cpars)[names(OM_k@cpars)!="qs"]]
        }

        # Episodic M
        # NK modified this section to apply to M-at-age
        if(scenario_i == "epiM") {
          OM_k@cpars$M_ageArray <- with(epiM_args,{
          M_mult <- rbinom(OM_k@proyears * OM_k@nsim, 1, yrprop) * pmin(exp(rnorm(OM_k@proyears * OM_k@nsim, 0, 2)), M_mult_max)
          M_mult_age <- rep(M_mult,each=OM_k@maxage+1) # Vector of multipliers repeating for each age
          M_array_hist <- OM_k@cpars$M_ageArray[,,1:OM_k@nyears]
          M_array_proj1 <- OM_k@cpars$M_ageArray[,,-(1:OM_k@nyears)]
          a1 <- as.numeric(aperm(M_array_proj1, perm = c(2, 1, 3))) # vectorize array and rearrange dimensions
          M_y <- a1 * (1 + M_mult_age) # Note that one is added to the multiplier so that the observed M is actually the minimum
          M_array_proj <- aperm(array(M_y, dim = c(OM_k@maxage+1,OM_k@nsim, OM_k@proyears)), perm = c(2, 1, 3))
           return(abind::abind(M_array_hist, M_array_proj, along = 3))
          })
        }

        # Index CV high
        if(scenario_i=="ucvhi"){
          # OM_k@Iobs <- OM_k@Iobs*ucvhi_args$scale
          OM_k@cpars$Data@CV_AddInd[,1,] <- OM_k@cpars$Data@CV_AddInd[,1,]*ucvhi_args$scale
        }

        # Index CV low
        if(scenario_i=="ucvlo"){
          #OM_k@Iobs <- OM_k@Iobs*ucvlo_args$scale
          OM_k@cpars$Data@CV_AddInd[,1,] <- OM_k@cpars$Data@CV_AddInd[,1,]*ucvlo_args$scale
        }

        # Index bias trend
        if(scenario_i=="ubias"){
          args <- get(paste0(scenario_i,"_args"))
          proyears <- OM_k@proyears
          nyears <- OM_k@nyears
          years <- OM_k@nyears+proyears

          # Setup empty array
          AddInd <- OM_k@cpars$Data@AddInd
          CV_AddInd <- OM_k@cpars$Data@CV_AddInd
          AddIerr_hist <- AddInd*NA
          AddIerr_proj <- array(NA,
                                dim=c(dim(AddIerr_hist)[1:2],proyears),
                                dimnames = list(dimnames(AddIerr_hist)[[1]],
                                                dimnames(AddIerr_hist)[[2]],
                                                rev(as.numeric(dimnames(AddIerr_hist)[[3]]))[1]+1:proyears
                                                )
                                )

          for(i in  1:dim(CV_AddInd)[2]){
            CV_AddInd_i <- CV_AddInd[,i,]
            AddIerr_hist[,i,] <- t(apply(CV_AddInd_i,1,function(x){
              lnorm_vector_boot(x=x/x,cv=x)
            }))
            AddIerr_proj[,i,] <- t(apply(CV_AddInd_i,1,function(x){
              x_proj <- sample(as.numeric(x[!is.na(x)]),size=proyears,replace=TRUE)
              lnorm_vector_boot(x=x_proj/x_proj,cv=x_proj)
              }))
          }
          AddIerr <- abind::abind(AddIerr_hist,AddIerr_proj,along=3)
          yr1 <- OM_k@nyears+args$yr1diff+1
          yrsmod <- yr1:max(years) # Years to modify
          yrdiff <- yrsmod-yr1
          AddIerr[,1,yrsmod] <- t(t(AddIerr[,1,yrsmod])+(args$int+args$slope*yrdiff))

          OM_k@cpars$AddIerr <- AddIerr
        }

        # Regime change (change in average recruitment deviations)
        if(scenario_i=="rc"){
          args <- get(paste0(scenario_i,"_args"))
          Perr_y <- OM_k@cpars$Perr_y
          years <- dim(Perr_y)[2]
          y_mult <- local({
            x <- rep(1,years)
            yr1 <- (years-OM_k@proyears)+args$yr1diff+1
            # Compute multiplier for transitional period
            slope <- (args$r2_mult-x[yr1])/args$transdur
            yrdiff <- 0:args$transdur
            yrtrans <- yr1+yrdiff
            x[yrtrans] <- 1+(slope*yrdiff)
            # Fill in multipliers for years of regime 2
            yrr2 <- (yr1+args$transdur+1):years
            x[yrr2] <- args$r2_mult
            x
          })
          val <- t(t(Perr_y)*y_mult)
          OM_k@cpars$Perr_y <- val
        }

        # Vulnerability constant (during historic and projection years)
        if(scenario_i=="vtiv"){
          OM_k@cpars$V <- local({
            V <- OM_k@cpars$V
            V2 <- aperm(V, perm = c(2, 1, 3))
            Vc <- V[1,,OM_k@nyears+1]
            V3 <- array(Vc,dim=dim(V2))
            aperm(V3, perm = c(2, 1, 3))
          })
        }

        # Perfect observation
        if(scenario_i=="perfobs"){
          OM_k <- Replace(OM_k, Perfect_Info)
          }

        # Generic fleet
        if(scenario_i=="genfle"){
          OM_k <- Replace(OM_k, Generic_Fleet)
        }

        if(scenario_i=="minerr"){
          # Minimize observation error
          OM_k <- Replace(OM_k, Perfect_Info)
          ## Minimize other sources of error (but don't reduce to zero)
          OM_k@cpars <- OM_k@cpars[names(OM_k@cpars)!="Perr_y"]
          OM_k@Perr <- c(0,0.05)
          OM_k@cpars$Data@CV_AddInd[!is.na(OM_k@cpars$Data@CV_AddInd)] <- 0.05
          OM_k@cpars$LenCV[!is.na(OM_k@cpars$LenCV)] <- 0.05

          # Make vulnerability constant (during historic and projection years)
          OM_k@cpars$V <- local({
            V <- OM_k@cpars$V
            V2 <- aperm(V, perm = c(2, 1, 3))
            Vc <- V[1,,OM_k@nyears+1]
            V3 <- array(Vc,dim=dim(V2))
            aperm(V3, perm = c(2, 1, 3))
          })
        }

        if(scenario_i=="minrecdev"){
          ## Minimize recruitment deviation
          OM_k@cpars <- OM_k@cpars[names(OM_k@cpars)!="Perr_y"]
          OM_k@Perr <- c(0,0.05)
        }

        if(scenario_i=="minrecac"){
        # Minimize recruitment autocorrelation
          OM_k@AC <- c(0,0.05)
        }



        if(scenario_i=="tachi"){
        # Set TAC to high value
          MSY_frac <- 1.25
        }else{
          MSY_frac <- MSY_frac_Init
        }

        if(scenario_i=="mconst"){
          OM_k@cpars <- OM_k@cpars[names(OM_k@cpars)!="M_ageArray"]
          OM_k@M <- rep(Data_k@Mort,2)
        }

        if(scenario_i=="noaddind"){
          # Clear AddInd and related slots from OM_k$cpars$Data
          # Data_e <- Data_empty()
          # # OM_k@cpars$Data
          # for(slotName_m in c("AddInd","CV_AddInd","AddIndV","AddIndType","AddIunits")){
          # slot(OM_k@cpars$Data,slotName_m) <- slot(Data_e,slotName_m)
          # }
          # Just remove the Data which contains only
          #OM_k@cpars <- OM_k@cpars[names(OM_k@cpars)[names(OM_k@cpars)!="Data"]]

          # Actually, let AddInd needs to exist because of the way SCA_lag is coded, but in
          # case it will be generated by the operating model but not used by SCA or the interim approaches
          AddInd_val <- "VB"
          AddInd_all_assess <- FALSE
        }else{
          AddInd_val <- AddInd_val_Init
          AddInd_all_assess <- AddInd_all_assess_Init
        }

        # No lag between stock assessment and management
        if(scenario_i=="nolag"){
          lag_Assess <- 0
        }else{
          lag_Assess <- lag_Assess_Init
        }

        source('fn/iMP.R') # Define MPs


        OM_k <- SubCpars(OM_k, sims = 1:nsim) # Limit number of simulations

        # Save OM to object
        OMName_ki <- paste0(OMName_k, "_", scenario_i)
        assign(OMName_ki,OM_k)
        # Save OM to file
        saveRDS(get(OMName_ki),
                file=paste0("OM_modified/",paste0(OMName_ki, ".rds")))

        ######## All MPs
        OM_k@interval <- c(1, 1, 1, 1, 1,
                           1,5,10,
                           1,1,
                           1,1,
                           1,1)
        #set.seed(myseed)
        # myHist_Init <- Simulate(OMInit_k)
        set.seed(myseed)
        # myHist <- Simulate(OM_k)

        message(paste("Run all MPs for:", OMName_ki))
        t_list <- c(t_list,Sys.time())
        message(paste0("at: ",tail(t_list,1),".(",round(diff(tail(t_list,2)),2)," since start)"))

        # Run all MPs together so that the Hist objects are always identical
        set.seed(myseed)
        sfExport(list = "SCA_lag")
        MSE_batch_1 <- runMSE(OM_k,
                              MPs = c("AvC", "CC1", "DCAC", "DBSRA", "SPMSY" # simple MPs
                                      ,"SCA_1","SCA_5","SCA_10"       # assessment only MPs
                                      ,"iMP_avg_5", "iMP_avg_10"       # interim average MPs
                                      ,"iMP_buffer_5","iMP_buffer_10"   # interim buffered MPs
                                      ,"pMP_5","pMP_10"                # projection MPs
                                      ),
                              parallel = runMSE_args$parallel, extended=runMSE_args$extended, silent=runMSE_args$silent)

        t_list <- c(t_list,Sys.time())
        message(paste0("batch 1 finished at ",tail(t_list,1),".(",round(diff(tail(t_list,2)),2)," duration"))


        res <- MSE_batch_1
        names(res@PPD) <- res@MPs

        saveRDS(res,file = paste0("MSE_obj/", MSEName_k, "_", scenario_i, ".rds"))
      }
    }
}

save.image("run_script.RData")

sfStop()
