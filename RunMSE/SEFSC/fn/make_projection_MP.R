# Author: Quang Huynh

make_projection_MP <- function(...) {
  fn <- projection_MP
  dots <- list(...)
  arg_ind <- pmatch(names(dots), names(formals(fn)))
  formals(fn)[arg_ind] <- dots
  class(fn) <- "MP"
  return(fn)
}