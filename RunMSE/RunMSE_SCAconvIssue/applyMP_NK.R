applyMP_NK <- function (Data, MPs = NA, reps = 100, nsims = NA, silent = FALSE) 
{
  print("applyMP_NK: running applyMP_NK")
  if (class(Data) != "Data") 
    stop("First argument must be object of class 'Data'", 
         call. = FALSE)
  Dataout <- Data
  if (is.na(nsims)) 
    nsims <- nrow(Data@Cat)
  nMPs <- length(MPs)
  if (.hasSlot(Data, "nareas")) {
    nareas <- Data@nareas
  }
  else {
    nareas <- 2
  }
  returnList <- list()
  recList <- list()
  TACout <- array(NA, dim = c(nMPs, reps, nsims))
  refMPs <- c("FMSYref", "FMSYref50", "FMSYref75", "NFref")
  runParallel <- snowfall::sfIsRunning()
  if (!silent) 
    message("Attempting to run ", length(MPs), " MPs:")
  for (mp in 1:nMPs) {
    if (!silent) 
      message(MPs[mp])
    mp_ns <- find(MPs[mp])
    dlmmp <- grepl("DLMtool", mp_ns)
    if (length(dlmmp) < 1) 
      dlmmp <- FALSE
    msemmp <- grepl("MSEtool", mp_ns)
    if (length(msemmp) < 1) 
      msemmp <- FALSE
    if (MPs[mp] %in% c("LBSPR", "LBSPR_MLL")) 
      dlmmp <- FALSE
    if (dlmmp | msemmp) 
      runParallel <- FALSE
    print("applyMP_NK: Made it here")
    if (runParallel) {
      temp <- try(snowfall::sfLapply(1:nsims, MPs[mp], 
                                     Data = Data, reps = reps), silent = TRUE)
    }
    else {
      temp <- try(lapply(1:nsims, MPs[mp], Data = Data, 
                         reps = reps), silent = TRUE)
    }
    mytemp <<-temp
    myrecList <<- recList
    mynsims <<- nsims
    myDataout <<- Dataout
    mynareas <<- nareas
    myTACout<<-TACout
    if (class(temp) == "try-error") {
      warning("Method ", MPs[mp], " failed with error: ", 
              temp)
    }
    else {
      slots <- slotNames(temp[[1]])
      for (X in slots) {
        if (X == "Misc") {
          rec <- lapply(temp, slot, name = X)
        }
        else {
          rec <- do.call("cbind", lapply(temp, slot, 
                                         name = X))
        }
        if (X == "Spatial") {
          rec <- matrix(rec, nareas, nsims, byrow = FALSE)
        }
        recList[[X]] <- rec
        for (x in 1:nsims) Dataout@Misc[[x]] <- recList$Misc[[x]]
        recList$Misc <- NULL
        if (MPs[mp] %in% refMPs) {
          recList$type <- "reference"
        }
        else {
          recList$type <- "mp"
        }
      }
      
      if (length(recList$TAC) > 0) 
        TACout[mp, , ] <- recList$TAC
      returnList[[mp]] <- recList
      if (!silent && any(apply(is.na(recList$TAC), 2, 
                               sum) > rep(0.5 * reps, nsims))) 
        message("Method ", MPs[mp], " produced greater than 50% NA values")
    }
  }
  Dataout@TAC <- TACout
  Dataout@MPs <- MPs
  nms <- names(Data@Misc)
  nms <- nms[nchar(nms) > 0]
  Dataout@Misc[nms] <- Data@Misc[nms]
  list(returnList, Dataout)
}
