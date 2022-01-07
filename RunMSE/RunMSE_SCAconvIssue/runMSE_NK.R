runMSE_NK <- function (OM = MSEtool::testOM, MPs = NA, Hist = FALSE, silent = FALSE,
                       parallel = FALSE, extended = FALSE, checkMPs = TRUE)
{
  CheckMPs <- MSEtool:::CheckMPs
  if (class(OM) == "OM") {
    if (OM@nsim <= 1)
      stop("OM@nsim must be > 1", call. = FALSE)
  }
  else if (class(OM) == "Hist") {
    if (!silent)
      message("Using `Hist` object to reproduce historical dynamics")
  }
  else {
    stop("You must specify an operating model")
  }
  if (checkMPs & !Hist)
    MPs <- CheckMPs(MPs = MPs, silent = silent)
  if (class(OM) == "OM") {
    HistSims <- Simulate(OM, parallel, silent)
  }
  else {
    HistSims <- OM
  }
  if (Hist) {
    if (!silent)
      message("Returning historical simulations")
    return(HistSims)
  }
  myHistSims <<- HistSims
  if (!silent)
    message("Running forward projections")
  MSEout <- try(Project_NK(Hist = HistSims, MPs, parallel, silent,
                        extended = extended, checkMPs = FALSE), silent = TRUE)
  if (class(MSEout) == "try-error") {
    message("The following error occured when running the forward projections: ",
            crayon::red(attributes(MSEout)$condition))
    message("Returning the historical simulations (class `Hist`). To avoid re-running spool up, ",
            "the forward projections can be run with ", "`runMSE(Hist, MPs, ...)`")
    return(HistSims)
  }
  MSEout
}
