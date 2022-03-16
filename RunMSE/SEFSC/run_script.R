library(abind)
library(bamExtras)
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
lag_Assess_init <- 2 # Number of years between terminal year of assessment and first year of management. May be modified in scenarios
MSY_frac_init <- 0.75 # Fraction of MSY for setting TAC. May be modified in scenarios

OM_name_scen_complete <- NA
#   local({
#   a <- gsub(".rds","",list.files("MSE_obj"))
#   b <- gsub("MSE","OM",a)
#   b
# })

OM_name_all <- gsub(".rds","",list.files("OM"))
OM_name <- OM_name_all#[!OM_name_all%in%OM_name_complete]
OM_name <- c("OM_BlackSeaBass", # Runs
             "OM_RedPorgy"#, # Runs
             ,"OM_VermilionSnapper" # Runs
             ,"OM_SnowyGrouper" # Runs 2022-1-20 batch 1 for base took only 25 min
             #"OM_RedGrouper", # 2022-1-19 This took 7 hours just to run the base scenario batch 1 so I interrupted it
             #,"OM_GagGrouper"#, # Runs
             #,"OM_GrayTriggerfish" # Problems with lightly fished scenario where "More than 5 % of simulations can't get to the specified level of depletion with these Operating Model parameters"
             #"OM_RedSnapper" # Seems to run but takes a long time
)

# Setup loops
scenario <- c(#"base", #"hs", "hd",
              #"nolag",
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
              "minrecac", # Minimize recruitment autocorrelation
              #"tachi",
              "vtiv"#,    # Vulnerability time-invariant
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

# seeds <- setNames(sample(1:10000,length(OM_name),replace = FALSE),OM_name)

MSEtool::setup(12) # Run in parallel over 12 cores
  sfLibrary("magrittr", character.only = TRUE, verbose = FALSE)

for(OM_name_k in OM_name) { ######### Loop over operating model

  MSE_name_k <- gsub("OM","MSE",OM_name_k)
  Data_name_k <- gsub("OM","Data",OM_name_k)

  myOM_init <- readRDS(paste0("OM/", OM_name_k, ".rds"))
  # myData_init <- readRDS(paste0("Data/", Data_name_k, ".rds"))


    for(scenario_i in scenario) { ######### Loop over scenario

      set.seed(myseed)
      # All scenarios
      OM_name_scen <- paste0(OM_name_k, "_", scenario_i)
      if(!OM_name_scen%in%OM_name_scen_complete){
        myOM <- myOM_init

        # Add indices observed in the assessment, their CVs, and their selectivities
          # myData <- new("Data")
          # slot(myData,"AddInd") <- slot(myData_init,"AddInd")
          # slot(myData,"CV_AddInd") <- slot(myData_init,"CV_AddInd")
          # slot(myData,"AddIndV") <- slot(myData_init,"AddIndV")1
          # myOM@cpars$Data <- myData
          # AddInd_colnums <- 1:dim(slot(myData_init,"AddInd"))[2]

        ## SCENARIOS

        # Hyperstable
        if(scenario_i == "hs") myOM@beta <- c(1/3, 2/3) # beta values below 1 lead to hyperstability

        # Hyperdeplete
        if(scenario_i == "hd") myOM@beta <- c(1.5, 3)

        # Depleted
        if(scenario_i == "dep") {
          myOM@cpars$D <- pmax(dep_args$scale * myOM@cpars$D, dep_args$min)
          # Remove qs so that runMSE will Optimize for user-specified depletion in last historical year
          myOM@cpars <- myOM@cpars[names(myOM@cpars)[names(myOM@cpars)!="qs"]]
        }

        # Lightly fished
        if(scenario_i == "lf")  {
          myOM@cpars$D <- pmin(lf_args$scale * myOM@cpars$D, lf_args$max)
          # Remove qs so that runMSE will Optimize for user-specified depletion in last historical year
          myOM@cpars <- myOM@cpars[names(myOM@cpars)[names(myOM@cpars)!="qs"]]
        }

        # Episodic M
        # NK modified this section to apply to M-at-age
        if(scenario_i == "epiM") {
          myOM@cpars$M_ageArray <- with(epiM_args,{
          M_mult <- rbinom(myOM@proyears * myOM@nsim, 1, yrprop) * pmin(exp(rnorm(myOM@proyears * myOM@nsim, 0, 2)), M_mult_max)
          M_mult_age <- rep(M_mult,each=myOM@maxage+1) # Vector of multipliers repeating for each age
          M_array_hist <- myOM@cpars$M_ageArray[,,1:myOM@nyears]
          M_array_proj1 <- myOM@cpars$M_ageArray[,,-(1:myOM@nyears)]
          a1 <- as.numeric(aperm(M_array_proj1, perm = c(2, 1, 3))) # vectorize array and rearrange dimensions
          M_y <- a1 * (1 + M_mult_age) # Note that one is added to the multiplier so that the observed M is actually the minimum
          M_array_proj <- aperm(array(M_y, dim = c(myOM@maxage+1,myOM@nsim, myOM@proyears)), perm = c(2, 1, 3))
           return(abind::abind(M_array_hist, M_array_proj, along = 3))
          })
        }

        # Index CV high
        if(scenario_i=="ucvhi"){
          # myOM@Iobs <- myOM@Iobs*ucvhi_args$scale
          myOM@cpars$Data@CV_AddInd[,1,] <- myOM@cpars$Data@CV_AddInd[,1,]*ucvhi_args$scale
        }

        # Index CV low
        if(scenario_i=="ucvlo"){
          #myOM@Iobs <- myOM@Iobs*ucvlo_args$scale
          myOM@cpars$Data@CV_AddInd[,1,] <- myOM@cpars$Data@CV_AddInd[,1,]*ucvlo_args$scale
        }

        # Index bias trend
        if(scenario_i=="ubias"){
          args <- get(paste0(scenario_i,"_args"))
          proyears <- myOM@proyears
          nyears <- myOM@nyears
          years <- myOM@nyears+proyears

          # Setup empty array
          AddInd <- myOM@cpars$Data@AddInd
          CV_AddInd <- myOM@cpars$Data@CV_AddInd
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
          yr1 <- myOM@nyears+args$yr1diff+1
          yrsmod <- yr1:max(years) # Years to modify
          yrdiff <- yrsmod-yr1
          AddIerr[,1,yrsmod] <- t(t(AddIerr[,1,yrsmod])+(args$int+args$slope*yrdiff))

          myOM@cpars$AddIerr <- AddIerr
        }

        # Regime change (change in average recruitment deviations)
        if(scenario_i=="rc"){
          args <- get(paste0(scenario_i,"_args"))
          Perr_y <- myOM@cpars$Perr_y
          years <- dim(Perr_y)[2]
          y_mult <- local({
            x <- rep(1,years)
            yr1 <- (years-myOM@proyears)+args$yr1diff+1
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
          myOM@cpars$Perr_y <- val
        }

        # Vulnerability constant (during historic and projection years)
        if(scenario_i=="vtiv"){
          myOM@cpars$V <- local({
            V <- myOM@cpars$V
            V2 <- aperm(V, perm = c(2, 1, 3))
            Vc <- V[1,,myOM@nyears+1]
            V3 <- array(Vc,dim=dim(V2))
            aperm(V3, perm = c(2, 1, 3))
          })
        }

        # Perfect observation
        if(scenario_i=="perfobs"){
          myOM <- Replace(myOM, Perfect_Info)
          }

        # Generic fleet
        if(scenario_i=="genfle"){
          myOM <- Replace(myOM, Generic_Fleet)
        }

        if(scenario_i=="minerr"){
          # Minimize observation error
          myOM <- Replace(myOM, Perfect_Info)
          ## Minimize other sources of error (but don't reduce to zero)
          myOM@cpars <- myOM@cpars[names(myOM@cpars)!="Perr_y"]
          myOM@Perr <- c(0,0.05)
          myOM@cpars$Data@CV_AddInd[!is.na(myOM@cpars$Data@CV_AddInd)] <- 0.05
          myOM@cpars$LenCV[!is.na(myOM@cpars$LenCV)] <- 0.05

          # Make vulnerability constant (during historic and projection years)
          myOM@cpars$V <- local({
            V <- myOM@cpars$V
            V2 <- aperm(V, perm = c(2, 1, 3))
            Vc <- V[1,,myOM@nyears+1]
            V3 <- array(Vc,dim=dim(V2))
            aperm(V3, perm = c(2, 1, 3))
          })
        }

        if(scenario_i=="minrecdev"){
          ## Minimize recruitment deviation
          myOM@cpars <- myOM@cpars[names(myOM@cpars)!="Perr_y"]
          myOM@Perr <- c(0,0.05)
        }

        if(scenario_i=="minrecac"){
        # Minimize recruitment deviations
          myOM@AC <- c(0,0.05)
        }



        if(scenario_i=="tachi"){
        # Set TAC to high value
          MSY_frac <- 1.25
        }else{
          MSY_frac <- MSY_frac_init
        }

        # No lag between stock assessment and management
        if(scenario_i=="nolag"){
          lag_Assess <- 0
          source('fn/iMP_nolag.R') # Define MPs
          #source('fn/iMP.R') # This should work too
        }else{
          lag_Assess <- lag_Assess_init
          source('fn/iMP.R') # Define MPs
        }

        myOM <- SubCpars(myOM, sims = 1:nsim) # Limit number of simulations

        # Save OM to object
        OM_name_ki <- paste0(OM_name_k, "_", scenario_i)
        assign(OM_name_ki,myOM)
        # Save OM to file
        saveRDS(get(OM_name_ki),
                file=paste0("OM_modified/",paste0(OM_name_ki, ".rds")))

        ######## All MPs
        myOM@interval <- c(1, 1, 1, 1, 1,
                           1,5,10,
                           1,1,
                           1,1,
                           1,1)
        #set.seed(myseed)
        # myHist_init <- Simulate(myOM_init)
        set.seed(myseed)
        # myHist <- Simulate(myOM)

        message(paste("Run all MPs for:", OM_name_ki))
        t_list <- c(t_list,Sys.time())
        message(paste0("at: ",tail(t_list,1),".(",round(diff(tail(t_list,2)),2)," since start)"))

        # Run all MPs together so that the Hist objects are always identical
        set.seed(myseed)
        sfExport(list = "SCA_lag")
        MSE_batch_1 <- runMSE(myOM,
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

        saveRDS(res,file = paste0("MSE_obj/", MSE_name_k, "_", scenario_i, ".rds"))
      }
    }
}

save.image("run_script.RData")

sfStop()
