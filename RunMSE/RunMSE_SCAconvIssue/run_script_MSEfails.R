library(openMSE)
library(TMB)

rm(list=ls())
t_list <- Sys.time()
set.seed(1358)

nsim <- 48
runParallel <- TRUE
ObsUser <- "Perfect_Info" # Set to NULL to avoid replacing Obs
ImpUser <- "Perfect_Imp"  # Set to NULL to avoid replacing Imp

OM_name <- "OM_RedGrouper"

# Setup loops
scenario <- "base"

seeds <- setNames(sample(1:10000,length(OM_name),replace = FALSE),OM_name)

MSEtool::setup(12) # Run in parallel over 12 cores

iMP_avg_10 <- make_interim_MP(SCA_Pope, .HCR = "HCR_MSY", MSY_frac = 1,
                              assessment_interval = 10,
                              type = "mean", type_par = 3)

for(OM_name_k in OM_name) { ######### Loop over operating model
  MSE_name_k <- gsub("OM","MSE",OM_name_k)

  myOM_init <- readRDS(paste0(OM_name_k, ".rds"))


    for(scenario_i in scenario) { ######### Loop over scenario
      # All scenarios
      OM_name_scen <- paste0(OM_name_k, "_", scenario_i)
      myOM <- myOM_init
      myOM <- SubCpars(myOM, sims = 1:nsim) # Limit number of simulations
      if(!is.null(ObsUser)) myOM <- Replace(myOM, get(ObsUser))
      if(!is.null(ImpUser)) myOM <- Replace(myOM, get(ImpUser))

      # # Save OM to object
      OM_name_ki <- paste0(OM_name_k, "_", scenario_i)

      ######## Fixed TAC MPs and Averaged Index MPs
      myOM@interval <- 1

      message(paste("Fixed TAC and Averaged MPs for:", OM_name_ki))
      t_list <- c(t_list,Sys.time())
      message(paste0("at: ",tail(t_list,1),".(",round(diff(tail(t_list,2)),2)," since start)"))

      MSE_batch_1 <- runMSE(myOM, MPs = "iMP_avg_10", parallel = TRUE)
      message(paste0("batch 1 finished at ",tail(t_list,1),".(",round(diff(tail(t_list,2)),2)," duration"))

      ######## Merge MSE and save output
      # res <- merge_MSE(MSE_batch_1
      # )

      #saveRDS(res,file = paste0("MSE_obj/", MSE_name_k, "_", scenario_i, ".rds"))
    }
}

sfStop()
