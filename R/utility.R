get_tg_thetas <- function(dat, mod, sum_score=NULL, percntrnk=NULL) {
  if (is.null(percntrnk)) {
    if (is.null(sum_score)) {
      stop("Either sum_score or percntrnk must be provided.")
    }
    N <- nrow(dat)
    scrjit <- sum_score + runif(nrow(dat), -0.5, 0.5)
    scrrnk <- vector("numeric", N)
    for (j in 1:N) {
      scrrnk[j] <- sum(scrjit <= scrjit[j])
    }
    percntrnk <- 100 * scrrnk / N
  }
  
  if (packageVersion("TestGardener") == "3.3.3") {
    thetas <- index_fun(percntrnk, SfdList = mod$parmListvec[[length(mod$parmListvec)]]$SfdList, chcemat = dat, crit = 1e-6)
    Result <- index_search(SfdList = mod$parmListvec[[length(mod$parmListvec)]]$SfdList, dat, 
                           thetas$index_out, thetas$Fval, thetas$DFval, thetas$D2Fval)
    return(Result$index)
  } else {
    thetas <- thetafun(percntrnk, WfdList = mod$parList[[length(mod$parList)]]$WfdList, U = dat, crit = 1e-6)
    Result <- thetasearch(mod$parList[[length(mod$parList)]]$WfdList, dat, thetas$theta_out, thetas$Hval, thetas$DHval, thetas$D2Hval)
    return(Result$theta)
  }
}

# ------------------
get_tg_probs <- function(theta, mod) {
  if (packageVersion("TestGardener") == "3.3.3") {
    iter <- length(mod$parmListvec)
    no_items <- length(mod$parmListvec[[iter]]$SfdList)
    probs <- list()
    for (item in 1:no_items) {
      WListi <- mod$parmListvec[[iter]]$SfdList[[item]]
      surpisal <- eval.surp(theta, WListi$Sfd, WListi$Zmat)
      probs[[item]] <- WListi$M^(-surpisal)
    }
  } else {
    iter <- length(mod$parList)
    no_items <- length(mod$parList[[iter]]$WfdList)
    probs <- list()
    for (item in 1:no_items) {
      WListi <- mod$parList[[iter]]$WfdList[[item]]
      surpisal <- eval.surp(theta, WListi$Wfd)
      probs[[item]] <- WListi$M^(-surpisal)
    }
  }
  
  probs
}

# ------------------
get_mirt_probs <- function(theta, mod) {
  probs <- list()
  for (item in 1:extract.mirt(mod, "nitems")) {
    extr <- extract.item(mod, item)
    probs[[item]] <- probtrace(extr, theta)
  }
  probs
}

# ------------------
log_likelihood <- function(dat, prob_list, log_base = exp(1)) {
  loglikelihood <- 0
  for (item in 1:length(prob_list)) {
    item_probs <- prob_list[[item]]
    item_resp_mat <- matrix(0, nrow = nrow(item_probs), ncol = ncol(item_probs))
    score_pos <- dat[, item]
    for (person in 1:length(score_pos)) {
      item_resp_mat[person, score_pos[person]] <- 1
    }
    loglikelihood <- loglikelihood +
      sum(
        item_resp_mat * log(item_probs, base = log_base)
      )
  }
  loglikelihood
}

