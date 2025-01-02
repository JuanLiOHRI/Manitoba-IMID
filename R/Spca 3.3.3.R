Spca <- function(SfdList, nharm=2, Sdim=NULL, rotate=TRUE, plotindex = NULL) {
  
  #  Last modified 12 December 2023 by Jim Ramsay
  
  #  set up the dimension of the over-space containing the test info curve
  
  n <- length(SfdList) # JL
  if (is.null(plotindex)) plotindex <- 1:n # JL
  
  if (is.null(Sdim)) {
    Sdim <- 0
    for (i in plotindex) {
      SfdListi <- SfdList[[i]]
      Sdim <- Sdim + SfdListi$M
    }
  }
  
  #  compute grid indfine the probability and surprisal values,
  #  and the first derivative of the surprisal values for the total 
  #  of the curves
  
  nfine   <- 101
  indfine <- seq(0,100,len=nfine)
  
  Pmat_full  <- matrix(0,nfine, Sdim)
  Smat_full  <- matrix(0,nfine, Sdim)
  DSmat_full <- matrix(0,nfine, Sdim)
  
  m2 <- 0
  for (item in plotindex) { # JL
    SListi <- SfdList[[item]]
    Sfdi   <- SListi$Sfd
    Mi     <- SListi$M
    Zmati  <- SListi$Zmat
    m1 <- m2 +  1
    m2 <- m2 + Mi
    Smat_full[,m1:m2]  <- eval.surp(indfine, Sfdi, Zmati)
    DSmat_full[,m1:m2] <- eval.surp(indfine, Sfdi, Zmati, 1)
    Pmat_full[,m1:m2]  <- Mi^(Smat_full[,m1:m2])
  }
  
  #  set up basis for (smoothing over [0,100])
  
  Snbasis <- 7
  Snorder <- 4
  Sbasis  <- fda::create.bspline.basis(c(0,100), Snbasis, Snorder) 
  SfdPar  <- fda::fdPar(fd(matrix(0,Snbasis,1),Sbasis))
  
  #  smooth all Sdim surprisal curves to get best fitting fd curves
  
  Sfd <- fda::smooth.basis(indfine, Smat_full, SfdPar)$fd
  
  #  functional PCA of the fd versions of surprisal curves
  #  the output is nharm principal component functions
  
  pcaList <- fda::pca.fd(Sfd, nharm, SfdPar, FALSE)
  
  #  set up the unrotated harmonic functional data object
  
  harmfd <- pcaList$harmonics
  
  #  compute variance proportions for unrotated solution
  
  varprop  <- pcaList$varprop
  if (!rotate) {
    varmxList    <- pcaList
    harmvarmxfd  <- harmfd
    varpropvarmx <- varprop
  } else {
    varmxList    <- varmx.pca.fd(pcaList)
    harmvarmxfd  <- varmxList$harmonics
    varpropvarmx <- varmxList$varprop
  }
  
  # JL: Order the principal components 
  order <- order(varpropvarmx, decreasing = TRUE)
  varpropvarmx <- varpropvarmx[order]
  
  harmvarmxfd$coefs <- harmvarmxfd$coefs[,order]
  colnames(harmvarmxfd$coefs) <- paste0("PC", 1:3)
  
  varmxList$harmonics$coefs <- varmxList$harmonics$coefs[,order]
  colnames(varmxList$harmonics$coefs) <- paste0("PC", 1:3)
  varmxList$values <- varmxList$values[order]
  varmxList$scores <- varmxList$scores[,order]
  varmxList$varprop <- varmxList$varprop[order]
  varmxList$rotmat <- varmxList$rotmat[,order]
  
  return(list(varmxList=varmxList, harmvarmxfd=harmvarmxfd, varpropvarmx=varpropvarmx))
  
}