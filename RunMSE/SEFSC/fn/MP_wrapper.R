# Script copied 2022-02-25 from: https://github.com/Blue-Matter/WSKJ/blob/master/03a_MPs.R#L6
# Author: Quang Huynh

# MP wrapper that puts the BB fleet in Data@Ind
MP_wrapper <- function(MP, delay = 2, ...) {
  MP <- substitute(MP)
  MP_call <- as.call(c(MP, x = quote(x), Data = quote(Data), reps = quote(reps), list(...)))
  force(delay)

  MP_body <- bquote({
    Data@Ind <- Data@AddInd[, 1, ]
    Data@CV_Ind <- Data@CV_AddInd[, 1, ]
    if(delay > 0) {
      delay_ind <- seq(ncol(Data@Cat) - delay + 1, ncol(Data@Cat), 1)
      Data@Ind <- Data@Ind[, -delay_ind]
      Data@CV_Ind <- Data@CV_Ind[, -delay_ind]
      Data@Cat <- Data@Cat[, -delay_ind, drop = FALSE]
      Data@CV_Cat <- Data@CV_Cat[, -delay_ind, drop = FALSE]
      Data@Year <- Data@Year[-delay_ind]
    }
    Rec <- .(MP_call)
    return(Rec)
  })

  MP_out <- eval(call("function", as.pairlist(alist(x = 1, Data = , reps = 1)), MP_body))
  class(MP_out) <- "MP"
  return(MP_out)
}

