# Redefine MPs to use all additional indices in stock assessments

# Define custom stock assessment function
BAM_SCA_args <- list(SR="BH",
                     vulnerability = "logistic",
                     early_dev = "comp_onegen",
                     late_dev = "comp50",
                     CAA_dist = "multinomial",
                     AddInd = AddInd_colnums,
                     lag=2)

# Averaged Index and Buffered Index MPs with (numbers indicate assessment interval)
iMP_avg_5 <- do.call(make_interim_MP_lag,
                     c(list(assessment_interval = 5,type = "mean", type_par = 3),
                       MP_args, BAM_SCA_args)
                     )

iMP_avg_10 <- do.call(make_interim_MP_lag,
                      c(list(assessment_interval = 10,type = "mean", type_par = 3),
                        MP_args, BAM_SCA_args)
                      )

iMP_buffer_5 <-  do.call(make_interim_MP_lag,
                         c(list(assessment_interval = 5, type = "buffer", type_par = 1),
                           MP_args, BAM_SCA_args)
                         )

iMP_buffer_10 <-  do.call(make_interim_MP_lag,
                           c(list(assessment_interval = 10, type = "buffer", type_par = 1),
                             MP_args, BAM_SCA_args)
)

# Fixed TAC MPs (numbers indicate assessment interval) and annual assessment MPs
SCA_5 <- SCA_10 <- SCA_1 <- do.call(make_MP,
                                    c(list(diagnostic = "min"),
                                      MP_args, BAM_SCA_args)
                                    )
