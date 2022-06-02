# Averaged Index and Buffered Index MPs with (numbers indicate assessment interval)
iMP_avg_5 <- do.call(make_interim_MP_NK,
                     c(list(assessment_interval = 5,type = "mean", type_par = 3),
                       MP_args, BAM_SCA_args)
                     )

iMP_avg_10 <- do.call(make_interim_MP_NK,
                      c(list(assessment_interval = 10,type = "mean", type_par = 3),
                        MP_args, BAM_SCA_args)
                      )

iMP_buffer_5 <-  do.call(make_interim_MP_NK,
                         c(list(assessment_interval = 5, type = "buffer", type_par = 1),
                           MP_args, BAM_SCA_args)
                         )

iMP_buffer_10 <-  do.call(make_interim_MP_NK,
                           c(list(assessment_interval = 10, type = "buffer", type_par = 1),
                             MP_args, BAM_SCA_args)
                          )

# Fixed TAC MPs (numbers indicate assessment interval) and annual assessment MPs
SCA_5 <- SCA_10 <- SCA_1 <- do.call(make_MP_NK,
                                    c(MP_args, BAM_SCA_args)
                                    )
# Projection MPs
pMP_5 <- do.call(make_projection_MP_NK,
                 c(list(assessment_interval = 5,process_error = NULL, p_sim = 1,
                        Ftarget = expression(MSY_frac * Assessment@FMSY)),
                   MP_args, BAM_SCA_args)
                 )
pMP_10 <- do.call(make_projection_MP_NK,
                 c(list(assessment_interval = 10,process_error = NULL, p_sim = 1,
                        Ftarget = expression(MSY_frac * Assessment@FMSY)),
                   MP_args, BAM_SCA_args)
)
