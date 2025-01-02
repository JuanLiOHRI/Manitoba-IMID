

ICC.plot <- function(scrfine, SfdList, dataList, Qvec, binctr, plotType = "P", 
                     plotindex=1:n, plotrange=c(min(scrfine),max(scrfine)), 
                     shaderange = NULL, Wrng=c(0,5), DWrng=c(-0.2, 0.2), 
                     data_point = FALSE, ci = FALSE, 
                     titlestr = NULL, autoplot = FALSE,
                     ttlsz = NULL, axisttl = NULL, axistxt = NULL, 
                     lgdlab = NULL, lgdpos = "bottom",
                     plotMissing = TRUE, plotTitle = TRUE,
                     xlab = "Score Index") {
  
  # scrfine    Vector of length 101 containing plotting points
  # SfdList    List of Sbinsmth results
  # dataList   List of make.dataList
  # Qvec       vector of quantile positions
  # binctr     vector of bin centers
  # plotType   Type(s) of plot, default as "P" for probability, can also be "W" for surprisal, 
  #            "DW" for sensitivity, and any combination of the three
  
  # plotindex  a vector of item indices; set if uses want to focus on specific items 
  # plotrange  a vector of length 2; set if users want to focus on specific score range 
  # shaderange a list of length 2 vector(s); set if users want to gray out specific score range(s)
  
  # Wrng       make ylim of all suprisal plots the same
  # DWrng      make ylim of all sensitivity plots the same
  
  # titlestr   plot title
  # autoplot   in Vignette, plot all items in a batch
  
  # data_point default as FALSE; set to TRUE for plotting data points
  # ci         default as FALSE; set to TRUE for plotting confidence limits of probability and surprisal
  
  # ttlsz      font size of plot title
  # axisttl    font size of axis titles
  # axistxt    font size of axis labels
  # lgdlab     font size of legend labels
  
  # Last modified 14 February 2023 by Jim Ramsay
  
  n <- length(SfdList)
  
  optionList <- dataList$optList 

  if (!is.null(shaderange)) {
    if (is.list(shaderange))
    {
      shadeValid <- rep(TRUE, length(shaderange))
      for (ishade in 1:length(shaderange))
      {
        rng <- shaderange[[ishade]]
        if (rng[1] >= rng[2])
        {
          shadeValid[ishade] <- FALSE
        } else if (rng[1] < plotrange[1])
        {
          shaderange[[ishade]][1] = plotrange[1]
        } else if (rng[2] > plotrange[2])
        {
          shaderange[[ishade]][2] = plotrange[2]
        }
      }
      
      if (sum(shadeValid) > 0) {
        shaderangeV <- shaderange[shadeValid] 
        } else {
          shaderangeV <- NULL
        }
    } else {
      stop("shaderange is not a list.")
    }
  } else {
    shaderangeV <- NULL
  }
  
  #   Invoke plotCore that actually controls the plotting
  
  plotlist <- list() # 2022-06-17
  if (autoplot) {
    #  In Vignette, plot all items in a batch
    for (iplot in plotindex)
    {
      p <- plotCore(iplot, SfdList, optionList, dataList,
                    scrfine, titlestr, plotType,
                    plotrange, plotindex, shaderangeV, binctr, 
                    Wrng, DWrng, Qvec, data_point, ci, plotMissing, plotTitle, xlab, lgdpos)    
      #print(p)
      plotlist[[iplot]] <- p
    }
  } else {
    #  Usual case of user control of plotting
    nplot <- length(plotindex)
    jplot <- 1
    iplot <- plotindex[jplot]
    while (jplot >= 1 & jplot <= nplot) {
      p <- plotCore(iplot, SfdList, optionList, dataList,
                    scrfine, titlestr, plotType,
                    plotrange, shaderangeV, binctr, 
                    Wrng, DWrng, Qvec, data_point, ci, plotMissing, plotTitle, xlab, lgdpos)    
      #print(p)
      plotlist[[iplot]] <- p
      #  request next plot
      
      if (length(plotindex) > 1) {
        n1 <- readline(prompt=
                         paste("Press [N/n/enter] to next item;", 
                               "[P/p] to previous item;", 
                               "or position in plotindex: "))
        if (n1 == "N" | n1 == "n" | n1 == "")
        {
          jplot <- jplot + 1
          iplot <- plotindex[jplot]
        } else if (n1 == "P" | n1 == "p") {
          jplot <- jplot - 1
          iplot <- plotindex[jplot]
        } else {
          n1 <- suppressWarnings(as.integer(n1))
          if (!is.na(n1) & (n1 >= 1 & n1 <= nplot)){
            jplot <- n1
            iplot <- plotindex[jplot]
          } else {
            break
          }
        }
      } else {
        break
      }
      gc() # to  free up memory from the previous run    } # end iplot
    }
  } # end iplot
  
  return(p)
}

#  ----------------------------------------------------------------------------

plotCore   <- function(iplot, SfdList, optionList, dataList,
                       scrfine, titlestr, plotType,
                       plotrange, shaderangeV, binctr, 
                       Wrng, DWrng, Qvec, data_point, ci,
                       plotMissing, plotTitle, xlab, lgdpos) {
  SListi    <- SfdList[[iplot]]
  itemStri  <- optionList$itemLab[iplot]  # Juan Li 2021-02-18
  optStri   <- optionList$optLab[[iplot]] # Juan Li 2021-02-18
  
  Sfdi       <- SListi$Sfd
  Mi         <- SListi$M # 2020-12-14
  Pbini      <- SListi$Pbin
  Sbini      <- SListi$Sbin
  Sfitfinei  <- SListi$Smatfine
  DSfitfinei <- SListi$DSmatfine
  Pfitfinei  <- SListi$Pmatfine
  
  Mi_label <- Mi # 2024-08-16 plotMissing == F, only omit the corresponding curve, but the surprisal were still estimated with the orginal Mi
  if (!plotMissing)
  {
    Mi         <- Mi - 1
    Pbini      <- SListi$Pbin[,1:Mi]
    Sbini      <- SListi$Sbin[,1:Mi]
    Sfitfinei  <- SListi$Smatfine[,1:Mi]
    DSfitfinei <- SListi$DSmatfine[,1:Mi]
    Pfitfinei  <- SListi$Pmatfine[,1:Mi]
  }
  
  if (ci) # plot confidence interval
  {
    PStdErr    <- SListi$PStdErr   # Probabilities over fine mesh
    SStdErr    <- SListi$SStdErr
    
    if (!plotMissing)
    {
      PStdErr    <- SListi$PStdErr[,1:Mi]
      SStdErr    <- SListi$SStdErr[,1:Mi]
    }
  } else
  {
    PStdErr    <- NULL
    SStdErr    <- NULL
  }
  
  if (is.null(dataList$key))
  {
    keyi <- NULL
  } else
  {
    keyi <- dataList$key[iplot]
  }
  
  if (!is.null(itemStri))
  {
    ttllab <- paste(dataList$titlestr,' ',iplot,': ',itemStri,sep="") 
  } else
  {
    if (!is.null(titlestr))
    {
      ttllab <- paste(titlestr, iplot,sep="")
    } else
    {
      ttllab <- paste(dataList$titlestr, iplot,sep="")
    }
  }
  
  # ----------- Juan Li 2021-02-17 ---------
  if (!is.null(optStri)) 
  {
    optionVec <- optStri
  } else
  {
    optionVec <- NULL
  }
  # -----------------------------------------
  nplot  <- length(plotType)
  pList <- list()
  for (itype in 1:nplot) {
    if (plotType[itype] == "P") {
      # probability
      pList[[itype]] <- plotICC(Mi, scrfine, Pbini, Pfitfinei, 
                                Qvec, binctr, c(0,1), 0.5, "Proportion/Probability",
                                plotrange, shaderangeV, data_point, ci, 
                                PStdErr, keyi, optVec=optionVec, xlab = xlab, lgdpos = lgdpos)
    } else {
      if (plotType[itype] == "S" || plotType[itype] == "W") {
        ystr <- paste("Surprisal (",Mi_label,"-bits)",sep="")
        # surprisal
        pList[[itype]] <- plotICC(Mi, scrfine, Sbini, Sfitfinei, 
                                  Qvec, binctr, Wrng, 0, ystr, 
                                  plotrange, shaderangeV, data_point, ci,
                                  SStdErr, keyi, optVec=optionVec, xlab = xlab, lgdpos = lgdpos)
      } else {
        if (plotType[itype] == "DS" || plotType[itype] == "DW") {
          # sensitivity
          pList[[itype]] <- plotICC(Mi, scrfine, NULL, DSfitfinei, 
                                Qvec, binctr, DWrng, 0, "Sensitivity", 
                                plotrange, shaderangeV, data_point, ci,
                                NULL, keyi, optVec=optionVec, xlab = xlab, lgdpos = lgdpos)
        } else {
          stop("Can't recognize the plot type.")
        }
      }
    }
  }
  
  if(length(pList) > 1) {
    p <- ggpubr::ggarrange(plotlist = pList, ncol = 1) #, common.legend = TRUE,legend="bottom"
  } else {
    p <- pList[[1]]
  }
  
  
  if (plotTitle)
  {
    p <- ggpubr::annotate_figure(p, top = ggpubr::text_grob(ttllab, size = 16))
  }
  
  return(p)
}

#  ----------------------------------------------------------------------------

plotICC   <- function(Mi, scrfine, bin1, fitfinei, 
                      Qvec, binctr, range, intercept, ylabel, 
                      plotrange, shaderange, data_point, ci,
                      StdErr=NULL, keyi=NULL, 
                      ttllab=NULL, optVec=NULL,
                      ttlsz = NULL, axisttl = NULL, axistxt = NULL, lgdlab = NULL,
                      xlab = "Score Index", lgdpos = "bottom")
{
  # Last modified 19 January 2022 by Juan Li
  value    <- 0
  variable <- 0
  if (!is.null(binctr)) nbin <- length(binctr)
  linesize <- 1
  ptssize  <- 2
  ind  <- 1:Mi 
  if (is.null(keyi)) indW <- NULL  else indW  <- ind[ind!=keyi] 
  
  # ----------- Juan Li 2021-02-17 ---------
  if (is.null(optVec)) 
  {
    optVec <- rep("",Mi)
    for(m in 1:Mi)
      optVec[m] <- m #paste("Option: ",m, sep="")
  }
  # -----------------------------------------
  default_size =16
  default_size1=12
  default_size2=10
  My_Theme = ggplot2::theme(
    plot.title  = ggplot2::element_text(size = ifelse(is.null(ttlsz),  default_size,ttlsz)),
    axis.title  = ggplot2::element_text(size = ifelse(is.null(axisttl),default_size1,axisttl)),
    axis.text   = ggplot2::element_text(size = ifelse(is.null(axistxt),default_size2,axistxt)),
    legend.text = ggplot2::element_text(size = ifelse(is.null(lgdlab), default_size2,lgdlab)))
  
  ymin <- NULL
  ymax <- NULL
  ymin_r <- NULL
  ymax_r <- NULL
  
  if (!is.null(keyi)) {
    
    #  ---------------------------------------------------------------------
    #       Plot probability curve(s) for multiple choice test items
    #  ---------------------------------------------------------------------
    
    # plot the curve of right option
    dffine_r <- data.frame(scrfine=scrfine,
                           fitfinei_r= fitfinei[,keyi])
    pp <- ggplot2::ggplot(
      data=dffine_r, ggplot2::aes(scrfine,fitfinei[,keyi])) +
      ggplot2::geom_line(ggplot2::aes(lty=optVec[keyi]),
                         color="blue", linewidth=linesize*2, na.rm = TRUE) +
      ggplot2::scale_linetype('Right') +
      ggplot2::geom_hline(yintercept = intercept,  color="black", linetype = "dashed") +
      ggplot2::geom_vline(xintercept = Qvec, color="black", linetype = "dashed")
    
    # plot the curve of wrong options
    dffine_w <- data.frame(scrfine=scrfine, fitfinei_w= fitfinei[,indW])
    names(dffine_w)[2:ncol(dffine_w)] <- optVec[indW]
    dffine_w <- tidyr::gather(dffine_w, key = "variable", value = "value", -scrfine)
    
    pp <- pp +
      ggplot2::geom_line(data = dffine_w, ggplot2::aes(scrfine,value,color=variable), 
                         linewidth=linesize, na.rm = TRUE) +
      ggplot2::scale_colour_hue(name="Wrong")  #,palette=scales::hue_pal(direction = -1)
    
  } else {
    
    #  ---------------------------------------------------------------------
    #             Plot probability curves for scale items
    #  ---------------------------------------------------------------------
    
    dffine <- data.frame(scrfine=scrfine, fitfinei = fitfinei[,ind]) 
    names(dffine)[2:ncol(dffine)]<- optVec[ind]
    dffine <-   tidyr::gather(dffine, key = "variable", value = "value", -scrfine)
    
    pp <- ggplot2::ggplot(dffine, ggplot2::aes(scrfine,value,color=variable)) +
      ggplot2::geom_line(linewidth=linesize, na.rm = TRUE) +
      ggplot2::scale_linetype('Right') +
      ggplot2::geom_hline(yintercept = intercept,  color="black", linetype = "dashed") +
      ggplot2::geom_vline(xintercept = Qvec, color="black", linetype = "dashed") +
      ggplot2::scale_colour_hue()  # palette=scales::hue_pal(direction = -1)
    
  }
  
  # Plot data points (optional)
  
  if (data_point & !is.null(bin1)) # plot observed proportions
  {
    if (!is.null(keyi)) {
      
      #  ---------------------------------------------------------------------
      #       Plot probability points for multiple choice test items
      #  ---------------------------------------------------------------------
      
      # plot data points of right option
      dfpts_r <- data.frame(binctr=binctr,
                            bin1_r= bin1[,keyi])
      
      pp <- pp + ggplot2::geom_point(
        data=dfpts_r, ggplot2::aes(binctr,bin1[,keyi]), shape = 21, fill = "blue", size = ptssize, 
        na.rm = TRUE)
      
      # plot data points of wrong options
      dfpts_w <- data.frame(binctr=binctr,
                            bin1_w= bin1[,indW]) 
      names(dfpts_w)[2:ncol(dfpts_w)]<- optVec[indW]
      dfpts_w <-   tidyr::gather(dfpts_w, key = "variable", value = "value", -binctr)
      
      pp <- pp +
        ggplot2::geom_point(data=dfpts_w, ggplot2::aes(binctr,value,fill=variable), 
                            shape = 21, size = ptssize, na.rm = TRUE)+
        ggplot2::scale_fill_hue(guide = "none") # ,palette=scales::hue_pal(direction = -1)
    } else {
      
      #  ---------------------------------------------------------------------
      #       Plot probability points for scale items
      #  ---------------------------------------------------------------------
      
      dfpts <- data.frame(binctr=binctr,
                          bin1= bin1[,ind]) 
      names(dfpts)[2:ncol(dfpts)]<- optVec[ind]
      dfpts <-   tidyr::gather(dfpts, key = "variable", value = "value", -binctr)
      
      pp <- pp +
        ggplot2::geom_point(data=dfpts, ggplot2::aes(binctr,value,fill=variable),
                            shape = 21, size = ptssize, na.rm = TRUE)+
        ggplot2::scale_fill_hue(guide = "none")  #, palette=scales::hue_pal(direction = -1)
    }
  }
  
  # Plot confidence interval (optional)
  if (ci & !is.null(StdErr) & !is.null(binctr))  
  {
    binfit <- matrix(0,nbin,Mi)
    for (m in 1:Mi) {
      binfit[,m] <- pracma::interp1(as.numeric(scrfine), as.numeric(fitfinei[,m]), 
                                    as.numeric(binctr))
    }
    
    if (!is.null(keyi)) {
      
      #  ---------------------------------------------------------------------
      #       Plot confidence intervals for multiple choice test items
      #  ---------------------------------------------------------------------
      
      # ci of thr right option
      ci_r <- data.frame(binctr = binctr,
                         ymin_r = binfit[,keyi]-2*StdErr[,keyi],
                         ymax_r = binfit[,keyi]+2*StdErr[,keyi])
      
      pp <- pp + ggplot2::geom_ribbon(data = ci_r, ggplot2::aes(ymin = ymin_r, ymax = ymax_r, x = binctr), 
                                      inherit.aes = FALSE, fill = "blue", alpha = 0.5)
      
      # ci of wrong options
      cimin_w <- data.frame(binctr = binctr,
                            ymin = binfit[,indW]-2*StdErr[,indW]) 
      names(cimin_w)[2:ncol(cimin_w)]<- optVec[indW]
      cimin_w <-   tidyr::gather(cimin_w, key = "variable", value = "ymin", -binctr)
      
      cimax_w <- data.frame(binctr = binctr,
                            ymax   = binfit[,indW]+2*StdErr[,indW]) 
      names(cimax_w)[2:ncol(cimax_w)]<- optVec[indW]
      cimax_w <-   tidyr::gather(cimax_w, key = "variable", value = "ymax", -binctr)
      
      ci_w    <- cimin_w %>% dplyr::left_join(cimax_w, by = c("binctr", "variable"))
      
      pp <- pp +
        ggplot2::geom_ribbon(data=ci_w, ggplot2::aes(ymin = ymin, ymax = ymax, x = binctr, fill=variable), 
                             inherit.aes = FALSE, alpha = 0.5)+
        ggplot2::scale_fill_hue(guide = "none") #,palette=scales::hue_pal(direction = -1)
      
    } else {
      
      #  ---------------------------------------------------------------------
      #       Plot confidence intervals for scale items
      #  ---------------------------------------------------------------------
      
      cimin <- data.frame(binctr = binctr,
                          ymin   = binfit[,ind]-2*StdErr[,ind]) 
      names(cimin)[2:ncol(cimin)]<- optVec[ind]
      cimin <-   tidyr::gather(cimin, key = "variable", value = "ymin", -binctr)
      
      cimax <- data.frame(binctr = binctr,
                          ymax   = binfit[,ind]+2*StdErr[,ind]) 
      names(cimax)[2:ncol(cimax)]<- optVec[ind]
      cimax <-   tidyr::gather(cimax, key = "variable", value = "ymax", -binctr)
      
      ci    <- cimin %>% dplyr::left_join(cimax, by = c("binctr", "variable"))
      
      pp <- pp +
        ggplot2::geom_ribbon(data=ci, ggplot2::aes(ymin = ymin, ymax = ymax, x = binctr, fill=variable), 
                             inherit.aes = FALSE, alpha = 0.5)+
        ggplot2::scale_fill_hue(guide = "none") #,palette=scales::hue_pal(direction = -1)
    }
  }
  
  if (!is.null(shaderange))
  {
    rectAlpha <- 0.9
    for (ishade in 1:length(shaderange))
    {
      pp <- pp + ggplot2::annotate("rect", xmin=shaderange[[ishade]][1], xmax=shaderange[[ishade]][2], 
                                   ymin=range[1], ymax=range[2],
                                   alpha = rectAlpha)
    }
  }
  
  pp <- pp +
    ggplot2::ylim(range[1], range[2]) +
    ggplot2::xlim(plotrange[1], plotrange[2]) +
    ggplot2::ylab(ylabel)
  
  pp <- pp +
    ggplot2::theme_bw()+
    ggplot2::theme(axis.title=ggplot2::element_text(size=16,face="bold"))
  
  if (xlab == "") {
    pp <- pp + ggplot2::theme(axis.title.x = element_blank())
  } else {
    pp <- pp + ggplot2::xlab(xlab) 
  }
  
  pp <- pp + ggplot2::theme(legend.position = lgdpos,
                            legend.title=ggplot2::element_blank()) # remove legend title
  
  pp <- pp + My_Theme
  return(pp)
}

