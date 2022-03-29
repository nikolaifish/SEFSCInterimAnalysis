SCA_NK <- function (x = 1, Data, AddInd = "B", SR = c("BH", "Ricker",
                                            "none"), vulnerability = c("logistic", "dome"), catch_eq = c("Baranov",
                                                                                                         "Pope"), CAA_dist = c("multinomial", "lognormal"), CAA_multiplier = 50,
          rescale = "mean1", max_age = Data@MaxAge, start = NULL,
          prior = list(), fix_h = TRUE, fix_F_equilibrium = TRUE,
          fix_omega = TRUE, fix_tau = TRUE, LWT = list(), early_dev = c("comp_onegen",
                                                                        "comp", "all"), late_dev = "comp50", integrate = FALSE,
          silent = TRUE, opt_hess = FALSE, n_restart = ifelse(opt_hess,
                                                              0, 1), control = list(iter.max = 2e+05, eval.max = 4e+05),
          inner.control = list(), ...)
{
  out <- SCA__NK(x = x, Data = Data, AddInd = AddInd, SR = SR,
              vulnerability = vulnerability, catch_eq = catch_eq,
              comp = "age", comp_dist = CAA_dist, comp_multiplier = CAA_multiplier,
              rescale = rescale, max_age = max_age, start = start,
              prior = prior, fix_h = fix_h, fix_F_equilibrium = fix_F_equilibrium,
              fix_omega = fix_omega, fix_tau = fix_tau, LWT = LWT,
              early_dev = early_dev, late_dev = late_dev, integrate = integrate,
              silent = silent, opt_hess = opt_hess, n_restart = n_restart,
              control = control, inner.control = inner.control, ...)
  return(out)
}
class(SCA_NK) <- "Assess"
