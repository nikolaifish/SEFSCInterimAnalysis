library(openMSE)
library(TMB)

rm(list=ls())
t_list <- Sys.time()
set.seed(1358)

nsim <- 48
runScenarios <- TRUE # Run scenarios to do MSE or just generate historical data?
runParallel <- TRUE
ObsUser <- "Perfect_Info" # Set to NULL to avoid replacing Obs
ImpUser <- "Perfect_Imp"  # Set to NULL to avoid replacing Imp

OM_name_scen_complete <- local({
  a <- gsub(".rds","",list.files("MSE_obj"))
  b <- gsub("MSE","OM",a)
  b
})

OM_name_all <- gsub(".rds","",list.files("OM"))
OM_name <- OM_name_all#[!OM_name_all%in%OM_name_complete]
OM_name <- c("OM_BlackSeaBass","OM_RedPorgy","OM_VermilionSnapper"
              #,"OM_RedGrouper", "OM_GagGrouper", "OM_GrayTriggerfish"
             )

# Setup loops
scenario <- c("base", "nr" #"hs", "hd",
              #"lf","dep",
              # "epiM"
)

seeds <- setNames(sample(1:10000,length(OM_name),replace = FALSE),OM_name)

MSEtool::setup(12) # Run in parallel over 12 cores

for(OM_name_k in OM_name) { ######### Loop over operating model
  source('fn/iMP_NK.R')

  MSE_name_k <- gsub("OM","MSE",OM_name_k)

  myOM_init <- readRDS(paste0("OM/", OM_name_k, ".rds"))

  if(runScenarios){
    for(scenario_i in scenario) { ######### Loop over scenario
      # All scenarios
      OM_name_scen <- paste0(OM_name_k, "_", scenario_i)
      if(!OM_name_scen%in%OM_name_scen_complete){
      myOM <- myOM_init
      myOM <- SubCpars(myOM, sims = 1:nsim) # Limit number of simulations
      if(!is.null(ObsUser)) myOM <- Replace(myOM, get(ObsUser))
      if(!is.null(ImpUser)) myOM <- Replace(myOM, get(ImpUser))

      # No recruitment process error scenario
      if(scenario_i == "nr") {
        myOM@cpars$Perr_y[] <- 1
      }

      # Hyperstable scenario
      if(scenario_i == "hs") myOM@beta <- c(1/3, 2/3) # beta values below 1 lead to hyperstability

      # Hyperdeplete scenario
      if(scenario_i == "hd") myOM@beta <- c(1.5, 3)

      # Depleted scenario
      if(scenario_i == "dep") {
        myOM@cpars$D <- 0.5 * myOM@cpars$D
        # Remove qs so that runMSE will Optimize for user-specified depletion in last historical year
        myOM@cpars <- myOM@cpars[names(myOM@cpars)[names(myOM@cpars)!="qs"]]
      }

      # Lightly fished scenario
      if(scenario_i == "lf")  {
        myOM@cpars$D <- 1.5 * myOM@cpars$D
        # Remove qs so that runMSE will Optimize for user-specified depletion in last historical year
        myOM@cpars <- myOM@cpars[names(myOM@cpars)[names(myOM@cpars)!="qs"]]
      }

      # Episodic M scenario
      # NK modified this section to apply to M-at-age
      if(scenario_i == "epiM") {
        set.seed(seeds[OM_name_k])
        M_mult <- rbinom(myOM@proyears * myOM@nsim, 1, 0.1) * pmin(exp(rnorm(myOM@proyears * myOM@nsim, 0, 2)), 4)
        M_mult_age <- rep(M_mult,each=myOM@maxage+1) # Vector of multipliers repeating for each age
        M_array_hist <- myOM@cpars$M_ageArray[,,1:myOM@nyears]
        M_array_future1 <- myOM@cpars$M_ageArray[,,-(1:myOM@nyears)]
        a1 <- as.numeric(aperm(M_array_future1, perm = c(2, 1, 3))) # vectorize array and rearrange dimensions
        M_array_future <- aperm(array(a1*(1+M_mult_age), dim = c(myOM@maxage+1,myOM@nsim, myOM@proyears)), perm = c(2, 1, 3))
        myOM@cpars$M_ageArray <- abind::abind(M_array_hist, M_array_future, along = 3)
      }

      # Save OM to object
      OM_name_ki <- paste0(OM_name_k, "_", scenario_i)
      assign(OM_name_ki,myOM)
      # Save OM to file
      saveRDS(get(OM_name_ki),
              file=paste0("OM_modified/",paste0(OM_name_ki, ".rds")))

      ######## Fixed TAC MPs and Averaged Index MPs
      myOM@interval <- c(10, 1)

      message(paste("Fixed TAC and Averaged MPs for:", OM_name_ki))
      t_list <- c(t_list,Sys.time())
      message(paste0("at: ",tail(t_list,1),".(",round(diff(tail(t_list,2)),2)," since start)"))

      MSE_batch_1 <- runMSE(myOM, MPs = c("SCA_10","iMP_avg_10"), parallel = FALSE)
      message(paste0("batch 1 finished at ",tail(t_list,1),".(",round(diff(tail(t_list,2)),2)," duration"))

      ####### Annual Assessment MP
      t_list <- c(t_list,Sys.time())
      myOM@interval <- 1

      MSE_batch_2 <- runMSE(myOM, MPs = "SCA_1", parallel = runParallel)
      t_list <- c(t_list,Sys.time())
      message(paste0("batch 2 finished at ",tail(t_list,1),".(",round(diff(tail(t_list,2)),2)," duration"))

      # ####### Perfect MP
      myOM@interval <- 1
       MSE_batch_4 <- runMSE(myOM, MPs = c("MP_PerfectMSY"), parallel = FALSE) # Doesn't work in parallel, but runs fast anyway (see ?Perfect)
      t_list <- c(t_list,Sys.time())
      message(paste0("batch 4 finished at ",tail(t_list,1),".(",round(diff(tail(t_list,2)),2)," duration"))

      ######## Merge MSE and save output
      res <- merge_MSE(MSE_batch_1,
                       MSE_batch_2,
                       #,MSE_batch_3
                       MSE_batch_4
      )

      saveRDS(res,file = paste0("MSE_obj/", MSE_name_k, "_", scenario_i, ".rds"))
    }
    }
  }
}

sfStop()
