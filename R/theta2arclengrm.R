deriv_mirt_graded = function(x, Theta, surp = T, surp_base = 2){
  Theta <- as.matrix(Theta)
  a <- mirt:::ExtractLambdas(x)
  P <- mirt:::ProbTrace(x, Theta, itemexp = FALSE)
  grad <- vector('list', x@ncat)
  for(i in seq_len(x@ncat))
    grad[[i]] <- matrix(0, nrow(Theta), x@nfact)
  for(j in seq_len(x@nfact)){
    for(i in seq_len(ncol(P)-1L)){
      w1 <- P[,i] * (1-P[,i]) * a[j]
      w2 <- P[,i+1L] * (1-P[,i+1L]) * a[j]
      grad[[i]][ ,j] <- w1 - w2
      if (surp) {
        grad[[i]][, j] <- -grad[[i]][, j] / (P[, i] * log(surp_base))
      }
    }
  }
  return(do.call(cbind, grad))
}

theta2arclengrm <- function(mod, thetas, total_cat, surp_base=rep(2, extract.mirt(mod, "nitems")), max_surp = NULL, item_ids = NULL, mesh_accuracy = 1001) {
  if (is.vector(surp_base) && length(surp_base) == 1) {
    surp_base=rep(surp_base, extract.mirt(mod, "nitems"))
  }
  derivs <- matrix(0, mesh_accuracy, total_cat)
  m2 <- 0
  if (is.null(item_ids)) item_ids <- 1:extract.mirt(mod, "nitems")
  # We want item category derivatives in columns and theta values in rows
  theta_mesh <- seq(min(thetas), max(thetas), length.out = mesh_accuracy)
  for (i in item_ids) {
    item <- extract.item(mod, i)
    m1 <- m2 + 1
    m2 <- m2 + item@ncat
    deriv_resc <- (max(theta_mesh) - min(theta_mesh)) / (length(theta_mesh) - 1)
    if ("graded" %in% mod@Model$itemtype) {
      derivs[, m1:m2] <- deriv_mirt_graded(item, theta_mesh, surp = T, surp_base[i]) * deriv_resc
    }
    if ("gpcm" %in% mod@Model$itemtype) {
      derivs[, m1:m2] <- deriv_mirt_gpc(item, theta_mesh, surp = T, surp_base[i]) * deriv_resc
    }
    if (!is.null(max_surp)) {
      capped_surp <- -log(probtrace(item, theta_mesh), surp_base[i]) > max_surp
      derivs[, m1:m2] <- (1 - capped_surp) * derivs[, m1:m2]
    }
  }

  al_mesh <- pracma::cumtrapz(sqrt(apply(derivs^2, 1, sum)))
  theta_al <- pracma::interp1(x = theta_mesh, y = as.numeric(al_mesh), xi = thetas)
  return(cbind(al = theta_al, theta = thetas))
}

item_arclength <- function(mod, thetas, surp_base=2, max_surp = NULL, item_ids = NULL, mesh_accuracy = 1001) {
  no_thetas <- length(thetas)
  m2 <- 0
  if (is.null(item_ids)) item_ids <- 1:extract.mirt(mod, "nitems")
  # We want item category derivatives in columns and theta values in rows
  theta_mesh <- seq(min(thetas), max(thetas), length.out = mesh_accuracy)
  scopes <- list()
  for (i in item_ids) {
    item <- extract.item(mod, i)
    m1 <- m2 + 1
    m2 <- m2 + item@ncat
    deriv_resc <- (max(theta_mesh) - min(theta_mesh)) / (length(theta_mesh) - 1)
    if ("graded" %in% mod@Model$itemtype) {
      derivs <- deriv_mirt_graded(item, theta_mesh, surp = T, surp_base) * deriv_resc
    }
    if ("gpcm" %in% mod@Model$itemtype) {
      derivs <- deriv_mirt_gpc(item, theta_mesh, surp = T, surp_base) * deriv_resc
    }
    if (!is.null(max_surp)) {
      capped_surp <- -log(probtrace(item, theta_mesh), surp_base) > max_surp
      derivs <- (1 - capped_surp) * derivs
    }
    al_mesh <- pracma::cumtrapz(sqrt(apply(derivs^2, 1, sum)))
    scopes[[i]] <- pracma::interp1(x = theta_mesh, y = as.numeric(al_mesh), xi = thetas)
  }
  return(scopes)
}