make_MP_NK <- function (.Assess, .HCR, diagnostic = c("min", "full", "none"),
          ...)
{
  diagnostic <- match.arg(diagnostic)
  if (is.character(.Assess)) {
    Assess_char <- .Assess
    .Assess <- as.symbol(.Assess)
  }
  else {
    .Assess <- substitute(.Assess)
    Assess_char <- as.character(.Assess)
  }
  if (is.character(.HCR)) {
    .HCR <- as.symbol(.HCR)
  }
  else {
    .HCR <- substitute(.HCR)
  }
  if (!inherits(eval(.Assess), "Assess")) {
    stop(paste(.Assess, "does not belong to class 'Assess'. Use: avail('Assess') to find eligible objects."))
  }
  if (!inherits(eval(.HCR), "HCR")) {
    stop(paste(.HCR, "does not belong to class 'HCR.' Use: avail('HCR') to find eligible objects."))
  }
  dots <- list(...)
  dots_in_Assess <- dots[match(names(formals(eval(.Assess))),
                               names(dots), nomatch = 0)]
  dots_in_HCR <- dots[match(names(formals(eval(.HCR))), names(dots),
                            nomatch = 0)]
  Assess_call <- as.call(c(.Assess, x = quote(x), Data = quote(Data),
                           dots_in_Assess))
  HCR_call <- as.call(c(.HCR, Assessment = quote(do_Assessment),
                        reps = quote(reps), dots_in_HCR))
  MP_body <- bquote({
    dependencies <- .(SAMtool:::get_dependencies(Assess_char, dots))
    do_Assessment <- .(Assess_call)
    Rec <- .(HCR_call)
    if (diagnostic != "none")
      Rec@Misc <- Assess_diagnostic_NK(x, Data, do_Assessment,
                                    include_assessment = .(diagnostic == "full"))
    # # Add other assessment outputs to Rec@Misc$diagnostic
    # Rec@Misc$diagnostic[[1]]$MSY <- do_Assessment@MSY
    if (!is.null(do_Assessment@info$Misc))
      Rec@Misc <- c(Rec@Misc, do_Assessment@info$Misc)
    return(Rec)
  })
  custom_MP <- eval(call("function", as.pairlist(alist(x = ,
                                                       Data = , reps = 1)), MP_body))
  formals(custom_MP)$diagnostic <- diagnostic
  return(structure(custom_MP, class = "MP"))
}
