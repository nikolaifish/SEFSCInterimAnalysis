
MP=SCA; lag = 2; Data <- rdat_to_Data(rdat_BlackSeaBass)

# MP wrapper that puts the BB fleet in Data@Ind
MP_wrapper <- function(MP, lag = 2, ...) {
  MP <- substitute(MP)
  MP_call <- as.call(c(MP, x = quote(x), Data = quote(Data), reps = quote(reps), list(...)
                       )
                     )
  force(lag)

  MP_body <- bquote({
    Data@Ind <- Data@AddInd[, 1, ]
    Data@CV_Ind <- Data@CV_AddInd[, 1, ]
    if(lag > 0) {
      lag_ind <- seq(ncol(Data@Cat) - lag + 1, ncol(Data@Cat), 1)
      Data@Ind <- Data@Ind[, -lag_ind]
      Data@CV_Ind <- Data@CV_Ind[, -lag_ind]
      Data@Cat <- Data@Cat[, -lag_ind, drop = FALSE]
      Data@CV_Cat <- Data@CV_Cat[, -lag_ind, drop = FALSE]
      Data@Year <- Data@Year[-lag_ind]
    }
    Rec <- .(MP_call)
    return(Rec)
  })

  MP_out <- eval(call("function", as.pairlist(alist(x = 1, Data = , reps = 1)), MP_body))
  class(MP_out) <- "MP"
  return(MP_out)
}

MP_wrapper2 <- function(MP, lag = 3, ...) {
  MP <- substitute(MP)
  MP_call <- as.call(c(MP, x = quote(x), Data = quote(Data), reps = quote(reps), list(...)
  )
  )
  force(lag)

  MP_body <- bquote({
    Data <- Data_sub_dim(Data,nyr_rm=lag)
    Rec <- .(MP_call)
    return(Rec)
  })

  MP_out <- eval(call("function", as.pairlist(alist(x = 1, Data = , reps = 1)), MP_body))
  class(MP_out) <- "MP"
  return(MP_out)
}
