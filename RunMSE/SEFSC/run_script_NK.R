library(openMSE)

rm(list=ls())
t1 <- Sys.time()
set.seed(1358)

runScenarios <- TRUE # Run scenarios to do MSE or just generate historical data?

# assessment_type <-  "SCA_Pope"
# HCR_type <-  "HCR_MSY"
# MSY_frac <- 0.75 # Passed to HCR_MSY and multiplied by Fmsy

# Setup loops
scenario <- c("base"#, #"hs", "hd",
              #"lf","dep",
              # "epiM"
               )
 OM_name <- c(
   "OM_RedPorgy"#,
   #"OM_BlackSeaBass",
   #"OM_VermilionSnapper"#,
   #"OM_RedGrouper",
   # "OM_RedSnapper",
   # "OM_SnowyGrouper"
             )
seeds <- setNames(sample(1:10000,length(OM_name),replace = FALSE),OM_name)

MSEtool::setup(12) # Run in parallel over 12 cores

for(OM_name_k in OM_name) { ######### Loop over operating model
  source('fn/iMP_NK.R')
  #source('fn/iMP.R')

  MSE_name_k <- gsub("OM","MSE",OM_name_k)

  myOM_init <- readRDS(paste0("OM/", OM_name_k, ".rds"))


  if(runScenarios){
    for(scenario_i in scenario) { ######### Loop over scenario
      myOM <- myOM_init

      if(scenario_i == "hs") myOM@beta <- c(1/3, 2/3) # beta values below 1 lead to hyperstability
      if(scenario_i == "hd") myOM@beta <- c(1.5, 3)
      if(scenario_i == "dep") {
        myOM@cpars$D <- 0.5 * myOM@cpars$D
        # Remove qs so that runMSE will Optimize for user-specified depletion in last historical year
        myOM@cpars <- myOM@cpars[names(myOM@cpars)[names(myOM@cpars)!="qs"]]
        }
      if(scenario_i == "lf")  {
        myOM@cpars$D <- 1.5 * myOM@cpars$D
        # Remove qs so that runMSE will Optimize for user-specified depletion in last historical year
        myOM@cpars <- myOM@cpars[names(myOM@cpars)[names(myOM@cpars)!="qs"]]
      }

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
      myOM@interval <- c(5, 10, 1, 1)

      message(paste("Fixed TAC and Averaged MPs for:", OM_name_ki))
      t2 <- Sys.time()
      message(paste0("at: ",t2,".(",round(t2-t1,2)," since start)"))

      MSE_batch_1 <- runMSE(myOM, MPs = c("SCA_5", "SCA_10", "iMP_avg_5", "iMP_avg_10"),
                            parallel = TRUE)
      t3 <- Sys.time()
      message(paste0("batch 1 finished at ",t3,".(",round(t3-t2,2)," duration"))

      ####### Annual Assessment MP
      t4 <- Sys.time()
      myOM@interval <- 1

      MSE_batch_2 <- runMSE(myOM, MPs = "SCA_1", parallel = TRUE)
      t5 <- Sys.time()
      message(paste0("batch 2 finished at ",t5,".(",round(t5-t4,2)," duration"))

      ####### Buffered interim MPs
      myOM@interval <- 1
      MSE_batch_3 <- runMSE(myOM, MPs = c("iMP_buffer_5", "iMP_buffer_10"),
                            parallel = TRUE)
      t6 <- Sys.time()
      message(paste0("batch 3 finished at ",t6,".(",round(t6-t5,2)," duration"))

      ######## Merge MSE and save output
      res <- merge_MSE(MSE_batch_1,
                       MSE_batch_2,
                       MSE_batch_3
      )
      #saveRDS(MSE_batch_1,file = paste0("MSE_obj/", OM_name_k, "_", scenario_i,"_b1", ".rds"))
      #saveRDS(MSE_batch_2,file = paste0("MSE_obj/", OM_name_k, "_", scenario_i,"_b2", ".rds"))
      #saveRDS(MSE_batch_3,file = paste0("MSE_obj/", OM_name_k, "_", scenario_i,"_b3", ".rds"))
      saveRDS(res,file = paste0("MSE_obj/", MSE_name_k, "_", scenario_i, ".rds"))
    }
  }
}

sfStop()
