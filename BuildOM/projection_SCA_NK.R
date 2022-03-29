# pmin(FMort, max_F) # Don't let me exceed max_F


Assessment=res
constrain = c("F", "Catch")
Ftarget=res@FMSY
#Catch
p_years = 1
p_sim = 1
obs_error = list(array(1, c(1, 1, res@info$data$nsurvey)), matrix(1, 1, 1))
process_error = matrix(1, 1, 1)
max_F = 3
seed = 499

projection_SCA_NK <- function (Assessment,
                               constrain = c("F", "Catch"),
                               Ftarget,
                               Catch,
                               p_years = 50,
                               p_sim = 200,
                               obs_error,
                               process_error,
                               max_F = 3,
                               seed = 499,
                               ...) 
{
  constrain <- "F"#match.arg(constrain)
  TMB_report <- Assessment@TMB_report
  TMB_data <- Assessment@obj$env$data
  Pope <- TMB_data$catch_eq == "Pope"
  set.seed(seed)
  if (is.matrix(process_error)) {
    if (nrow(process_error) != p_sim) 
      stop("Number of rows of process_error not equal to ", 
           p_sim, call. = FALSE)
    if (ncol(process_error) != p_years) 
      stop("Number of columns of process_error not equal to ", 
           p_years, call. = FALSE)
    p_log_rec_dev <- process_error
  } else {
    tau <- ifelse(!is.null(process_error), process_error[1], 
                  TMB_report$tau)
    if (!is.null(tau) && tau > 0) 
      p_log_rec_dev <- exp(matrix(rnorm(p_years * p_sim, 
                                        0, tau), p_sim, p_years) - 0.5 * tau^2)
  }
  if (!exists("p_log_rec_dev", inherits = FALSE)) 
    p_log_rec_dev <- matrix(1, p_sim, p_years)
  nsurvey <- ifelse(is.matrix(Assessment@Index), Assessment@Index %>% 
                      ncol(), 1)
  if (is.array(obs_error[[1]])) {
    Iobs_err <- obs_error[[1]]
  } else {
    samps <- rnorm(p_years * p_sim * nsurvey, 0, obs_error[[1]]) %>% 
      array(c(nsurvey, p_sim, p_years)) %>% aperm(c(2, 
                                                    3, 1))
    Iobs_err <- exp(samps)
  }
  if (is.array(obs_error[[2]])) {
    Cobs_err <- obs_error[[2]]
  } else {
    omega <- ifelse(!is.null(obs_error[[2]]), obs_error[[2]], 
                    0)
    Cobs_err <- exp(matrix(rnorm(p_years * p_sim, 0, omega), 
                           p_sim, p_years))
  }
  p_output <- lapply(1:p_sim, SAMtool:::projection_SCA_internal, p_log_rec_dev = p_log_rec_dev, 
                     Cerr = Cobs_err, Ierr = Iobs_err, FMort = Ftarget, Catch = Catch, 
                     constrain = constrain, TMB_report = TMB_report, TMB_data = TMB_data, 
                     Pope = Pope, max_F = max_F)
  CAApred <- lapply(p_output, getElement, "CAApred") %>% simplify2array()
  Ipred <- array(NA_real_, dim(Iobs_err)[c(2, 3, 1)])
  Ipred[] <- lapply(p_output, getElement, "Ipred") %>% simplify2array()
  output <- new("project",
                Catch = do.call(rbind, lapply(p_output, getElement, "Cpred")),
                C_at_age = aperm(CAApred, c(3,1, 2)),
                Index = Ipred %>% aperm(c(3, 1, 2)),
                SSB = do.call(rbind, lapply(p_output, getElement, "E")),
                R = do.call(rbind, lapply(p_output, function(x) x$N[, 1])),
                N = do.call(rbind, lapply(p_output, function(x) rowSums(x$N))),
                VB = do.call(rbind, lapply(p_output, getElement, "VB")),
                B = do.call(rbind, lapply(p_output, getElement, "B")),
                FMort = do.call(rbind, lapply(p_output, getElement, "Fout")))
  return(output)
}
