Analyze <- function (index, indexQnt, dataList, NumDensBasis = 7, norder = 4, ncycle = 10, 
          itdisp = FALSE, verbose = FALSE) 
{
  nbin <- dataList$nbin
  markers <- dataList$PcntMarkers/100
  bdry0 <- seq(0, 2 * nbin, 1)/(2 * nbin)
  parmListvec <- vector("list", ncycle)
  pdffinemat <- matrix(0, 101, ncycle)
  Qvecmat <- matrix(0, 5, ncycle)
  logdensbasis <- fda::create.bspline.basis(c(0, 100), NumDensBasis, norder) # JL: BUG in the CRAN version 3.3.3
  SfdList <- dataList$SfdList
  n <- length(SfdList)
  Sdim <- 0
  for (i in 1:n) {
    SStri <- SfdList[[i]]
    Sdim <- Sdim + SStri$M
  }
  chcemat <- dataList$chcemat
  indfine <- seq(0, 100, len = 101)
  HALmat <- matrix(0, ncycle, 2)
  for (icycle in 1:ncycle) {
    print(paste("----------  Cycle ", icycle, "-----------"))
    SfdResult <- Sbinsmth(index, dataList)
    SfdList <- SfdResult$SfdList
    binctr <- SfdResult$binctr
    bdry <- SfdResult$bdry
    freq <- SfdResult$freq
    Fvec <- Ffun(index, SfdList, chcemat)
    meanF <- mean(Fvec)
    HALmat[icycle, 1] <- meanF
    if (verbose) 
      print(paste("Mean data fit = ", round(meanF, 3)))
    indexfunList <- index_fun(index, SfdList, chcemat, 20, 
                              0.001, itdisp = itdisp)
    index <- indexfunList$index_out
    Fval <- indexfunList$Fval
    DFval <- indexfunList$DFval
    D2Fval <- indexfunList$D2Fval
    active <- indexfunList$active
    indexdens <- index[0 < index & index < 100]
    index_distnList <- index_distn(indexdens, logdensbasis)
    pdffine <- index_distnList$pdffine
    denscdf <- as.numeric(index_distnList$denscdf)
    indcdf <- as.numeric(index_distnList$indcdf)
    Qvec <- pracma::interp1(denscdf, indcdf, markers)
    pdffinemat[, icycle] <- pdffine
    Qvecmat[, icycle] <- Qvec
    binctr <- indexQnt[seq(2, 2 * nbin, 2)]
    bdry <- indexQnt[seq(1, 2 * nbin + 1, 2)]
    DSfine <- matrix(0, 101, Sdim)
    m2 <- 0
    for (i in 1:n) {
      SListi <- SfdList[[i]]
      Mi <- SListi$M
      m1 <- m2 + 1
      m2 <- m2 + Mi
      DSfine[, m1:m2] <- SListi$DSmatfine
    }
    indfine <- seq(0, 100, len = 101)
    infoSurp <- max(pracma::cumtrapz(indfine, sqrt(apply(DSfine^2, 
                                                         1, sum))))
    HALmat[icycle, 2] <- infoSurp
    if (verbose) 
      print(paste("infoSurp in bits = ", round(infoSurp, 
                                               1)))
    Result <- index_search(SfdList, dataList$chcemat, index, 
                           Fval, DFval, D2Fval)
    index <- Result$index
    Fval <- Result$Fval
    DFval <- Result$DFval
    D2Fval <- Result$D2Fval
    parmList_icycle <- list(index = index, indexQnt = indexQnt, 
                            SfdList = SfdList, meanF = meanF, binctr = binctr, 
                            bdry = bdry, freq = freq, Qvec = Qvec, Fval = Fval, 
                            DFval = DFval, D2Fval = D2Fval, active = active, 
                            infoSurp = infoSurp)
    parmListvec[[icycle]] <- parmList_icycle
  }
  return(list(parmListvec = parmListvec, pdffinemat = pdffinemat, 
              Qvecmat = Qvecmat, HALmat = HALmat))
}
