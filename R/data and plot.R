data4Plot <- function(test_str, model = "TG", test = "PHQ-9", i_test = NULL) {
  if (model == "TG") {
    if (packageVersion("TestGardener") == "3.3.3") {
      filestr       <- paste0("data/saved results/", test_str, "_")
    } else {
      filestr       <- paste0("data/saved results/3.2.6/", test_str, "_")
    }
    dataList      <- readRDS(paste0(filestr, "dataList.rds"))
    analyzeResult <- readRDS(paste0(filestr, "analyzeResult.rds"))
    infoList      <- readRDS(paste0(filestr, "infoList.rds"))
    
    n <- dataList$n
    noption <- dataList$noption
    
    if (packageVersion("TestGardener") == "3.3.3") {
      ncycle        <- length(analyzeResult$parmListvec)
      parList       <- analyzeResult$parmListvec[[ncycle]]
      WfdList       <- parList$SfdList
      arclengthvec  <- infoList$infoSurpvec
      Qvec_al       <- infoList$Qinfovec
      binctr_al     <- infoList$bininfoctr
    } else {
      ncycle        <- length(analyzeResult$parList)
      parList       <- analyzeResult$parList[[ncycle]]
      WfdList       <- parList$WfdList
      arclengthvec  <- infoList$arclengthvec
      Qvec_al       <- infoList$Qvec_al
      binctr_al     <- infoList$binctr_al
    }
    Qvec          <- parList$Qvec
    binctr        <- parList$binctr
    itemArclength <- matrix(0,n,1)
    itemSurpmat   <- matrix(0,101, n)
    indfine  <- seq(0,100,len=101)
    for (item in 1:n) {
      if (packageVersion("TestGardener") == "3.3.3") {
        Resulti <- index2info(indfine, Qvec, WfdList, NULL, item, c(0,100), TRUE) 
        itemArclength[item] <- Resulti$infoSurp
        itemSurpmat[,item] <- Resulti$infoSurpvec
      } else {
        Resulti <- theta2arclen(indfine, Qvec, WfdList, NULL, item, c(0,100), TRUE)
        itemArclength[item] <- Resulti$arclength
        itemSurpmat[,item] <- Resulti$arclengthvec
      }
    }
    return(list(test          = test,
                model         = NULL,
                n             = n,
                noption       = noption,
                arclengthvec  = arclengthvec,
                WfdList       = WfdList,
                dataList      = dataList,
                Qvec_al       = Qvec_al,
                binctr_al     = binctr_al,
                itemArclength = itemArclength,
                itemSurpmat   = itemSurpmat))
  } else {
    filestr <- paste0("data/saved results/grm_", test_str)
    grm <- readRDS(paste0(filestr, ".rds"))
    df_theta <- read.csv("data/saved results/grm_df_theta.csv", header = T)
    df_theta_al <- read.csv("data/saved results/grm_df_theta_al.csv", header = T)
    item_al <- readRDS("data/saved results/grm_item_al.rds")
    
    theta <- df_theta[, i_test]
    Qvec   <- quantile(theta, c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = T)
    theta_interval <- seq(min(theta, na.rm = T), max(theta, na.rm = T), len = 2001)
    theta_al <- df_theta_al[, i_test]
    Qvec_al   <- quantile(theta_al, c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = T)
    itemArclength <- item_al[[i_test]]
    
    n <- grm@Data$nitems
    noption <- rep(0, n)
    for (i in 1:n) noption[i] <- length(unique(grm@Data$data[, i]))
    
    return(list(test            = test,
                model           = grm,
                n               = n,
                noption         = noption,
                theta_interval  = theta_interval,
                Qvec            = Qvec,
                Qvec_al         = Qvec_al,
                itemArclength = itemArclength))
  }
}

# ----------------------------
multi.plot <- function(data_2b, data_4b, data_grm, titlestr = "Item ", plotIndex = NULL) {
  plot_w_2b <- list()
  plot_p_2b <- list()
  plot_w_4b <- list()
  plot_p_4b <- list()
  plot_w_grm <- list()
  plot_p_grm <- list()
  
  # Plot
  n       <- data_2b$n
  noption <- data_2b$noption
  test    <- data_2b$test
  
  if (is.null(plotIndex)) {
    plotIndex <- 1:n
  }
  
  for (i in 1:length(plotIndex)) {
    plotTitle <- T
    
    if (i == length(plotIndex)) {
      ldgendPos <- "bottom"
      xlab = paste0("Test Information (", noption[i], "-bits)")
    } else {
      ldgendPos <- "none"
      xlab = ""
    }
    
    ploti <- plotIndex[i]
    # GRM
    plots_grm <- grm_4_icc_plots(data_grm$model, i, titlestr, data_grm$theta_interval, 
                                 data_grm$Qvec, data_grm$Qvec_al, labsize = 12,
                                 xlab = xlab, lgdpos = ldgendPos)
    plot_w_grm[[i]]   <- plots_grm$p3 +
      labs(title = paste0("Item ", ploti,", scope ",round(data_grm$itemArclength[ploti],1)," (", noption[i], "-bits)"))
    plot_p_grm[[i]]   <- plots_grm$p4 +
      labs(title = paste0("Item ", ploti,", scope ",round(data_grm$itemArclength[ploti],1)," (", noption[i], "-bits)"))
    
    # TG: 2b
    plot_w_2b[[i]] <- ICC.plot(data_2b$arclengthvec, data_2b$WfdList, data_2b$dataList, data_2b$Qvec_al, 
                               data_2b$binctr_al, data_point=TRUE, 
                               plotType=c("W"), Wrng=c(0,4), plotindex=i, plotMissing = T,
                               plotTitle = FALSE, xlab = xlab, lgdpos = ldgendPos) +
      labs(title = paste0("Item ", ploti,", scope ",round(data_2b$itemArclength[ploti],1)," (", noption[i], "-bits)"))
    
    plot_p_2b[[i]] <- ICC.plot(data_2b$arclengthvec, data_2b$WfdList, data_2b$dataList, data_2b$Qvec_al, 
                               data_2b$binctr_al, data_point=TRUE, 
                               plotType=c("P"), Wrng=c(0,4), plotindex=i, plotMissing = T,
                               plotTitle = FALSE, xlab = xlab, lgdpos = ldgendPos) +
      labs(title = paste0("Item ", ploti,", scope ",round(data_2b$itemArclength[ploti],1)," (", noption[i], "-bits)"))
    
    # TG: 4b
    plot_w_4b[[i]] <- ICC.plot(data_4b$arclengthvec, data_4b$WfdList, data_4b$dataList, data_4b$Qvec_al, 
                               data_4b$binctr_al, data_point=TRUE, 
                               plotType=c("W"), Wrng=c(0,4), plotindex=i, plotMissing = T,
                               plotTitle = FALSE, xlab = xlab, lgdpos = ldgendPos) +
      labs(title = paste0("Item ", ploti,", scope ",round(data_4b$itemArclength[ploti],1)," (", noption[i], "-bits)"))
    
    plot_p_4b[[i]] <- ICC.plot(data_4b$arclengthvec, data_4b$WfdList, data_4b$dataList, data_4b$Qvec_al, 
                               data_4b$binctr_al, data_point=TRUE, 
                               plotType=c("P"), Wrng=c(0,4), plotindex=i, plotMissing = T,
                               plotTitle = FALSE, xlab = xlab, lgdpos = ldgendPos) +
      labs(title = paste0("Item ", ploti,", scope ",round(data_4b$itemArclength[ploti],1)," (", noption[i], "-bits)"))
  }
  
  likelihoods <- readRDS("data/saved results/likelihoods.rds")
  
  ttlsize <- 16
  group_p_grm <- set_ggarrange(plot_p_grm)
  group_p_grm <- annotate_figure(group_p_grm, 
                                 top = text_grob(paste0("GRM (", 
                                                        round(likelihoods$value[which(likelihoods$test == test &
                                                                                        likelihoods$model == "GRM")],2), ")"),
                                                 face = "bold", size = ttlsize))
  
  group_w_grm <- set_ggarrange(plot_w_grm)
  group_w_grm <-annotate_figure(group_w_grm, 
                                top = text_grob(paste0("GRM (",
                                                       round(likelihoods$value[which(likelihoods$test == test &
                                                                                       likelihoods$model == "GRM")],2), ")"),
                                                face = "bold", size = ttlsize))
  
  group_p_2b <- set_ggarrange(plot_p_2b)
  group_p_2b <- annotate_figure(group_p_2b, 
                                top = text_grob(paste0("TG: 2-basis (",
                                                       round(likelihoods$value[which(likelihoods$test == test &
                                                                                       likelihoods$model == "2-basis")],2), ")"),
                                                face = "bold", size = ttlsize))
  
  group_w_2b <- set_ggarrange(plot_w_2b)
  group_w_2b <- annotate_figure(group_w_2b, 
                                top = text_grob(paste0("TG: 2-basis (",
                                                       round(likelihoods$value[which(likelihoods$test == test &
                                                                                       likelihoods$model == "2-basis")],2), ")"),
                                                face = "bold", size = ttlsize))
  
  group_p_4b <- set_ggarrange(plot_p_4b)
  group_p_4b <- annotate_figure(group_p_4b, 
                                top = text_grob(paste0("TG: 4-basis (", 
                                                       round(likelihoods$value[which(likelihoods$test == test &
                                                                                       likelihoods$model == "4-basis")],2), ")"),
                                                face = "bold", size = ttlsize))
  
  group_w_4b <- set_ggarrange(plot_w_4b)
  group_w_4b <- annotate_figure(group_w_4b, 
                                top = text_grob(paste0("TG: 4-basis (", 
                                                       round(likelihoods$value[which(likelihoods$test == test &
                                                                                       likelihoods$model == "4-basis")],2), ")"),
                                                face = "bold", size = ttlsize))
  
  print(ggarrange(group_p_grm, group_p_2b, group_p_4b, ncol = 3))
  print(ggarrange(group_w_grm, group_w_2b, group_w_4b, ncol = 3))
}

# ---------------
set_ggarrange <- function(plotList) {
  n <- length(plotList)
  if (n == 2) {
    return(ggarrange(plotList[[1]], plotList[[2]], nrow = n, heights = c(rep(0.9, n-1), 1)))
  } else if (n == 8) {
    return(ggarrange(plotList[[1]], plotList[[2]], plotList[[3]], plotList[[4]], plotList[[5]],
                     plotList[[6]], plotList[[7]], plotList[[8]], nrow = n, heights = c(rep(0.9, n-1), 1)))
  } else {
    return(ggarrange(plotList[[1]], plotList[[2]], plotList[[3]], plotList[[4]], plotList[[5]],
                     plotList[[6]], plotList[[7]], plotList[[8]], plotList[[9]], nrow = n, heights = c(rep(0.9, n-1), 1)))
  }
}

