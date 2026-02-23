###########################################################
##### Visualization of Simulation Results             #####
##### Linear interactive fixed effects model for      ##### 
##### Synthetic control methods                       #####
##### (1) Constrained OLS                             #####
##### (2) Unconstrained OLS                           #####
##### (3) Proximal inference                          #####
##### kendrick.li@stjude.org                          #####
##### 02/23/2026                                      #####
###########################################################


rm(list=ls())

library(dplyr)
library(RColorBrewer)

colp.OLS <- rep(0, 3)
colp.SC <- rep(1, 3)
colp.NC <- rep(2, 3)
colp.NC_cons <- rep(5, 3)

add.col <- NA 

##########################################################
## Generate input data for plotting for SC methods without 
## baseline covariate adjustment
## input:
## - rslt.all: data.frame combining results. Should contain
##            the following columns:
##   - t0: number of pre-treatment periods
##   - n.units: number of units
##   - SC_est, SC_OLS_est, NC_SC_est, NC_SC_constrained_est:
##         Point estimates for different methods
##########################################################

getplotdata_nocov <- function(rslt.all, removeoutlier = F){

  n.rslt <- length(rslt.all)
  n.t0 <- length(t0.all)
  plotdat <- SCNC <- NULL
  mycol <- NULL
  
  for(n.units in unique(rslt.all$n.units)){
    SC <- OLS <- NC <- NC_cons <- list()
    for(t0 in unique(rslt.all$t0)){
      rowind <- (rslt.all$n.units == n.units) & (rslt.all$t0 == t0) 
      SC <- c(SC, list(rslt.all[rowind, "SC_est"]))
      OLS <- c(OLS, list(rslt.all[rowind, "SC_OLS_est"]))
      NC <- c(NC, list(rslt.all[rowind, "NC_SC_est"]))
      NC_cons <- c(NC_cons, list(rslt.all[rowind, "NC_SC_constrained_est"]))
    }
    
    SCNC <- c(SCNC, OLS, SC, NC, NC_cons)
    space <- NA
    if(n.units != max(rslt.all$n.units)){
      plotdat <- c(plotdat,
                   OLS, space,
                   SC, space,
                   NC, space,
                   NC_cons, space,
                   space, space)
      mycol <- c(mycol,
                 colp.OLS, add.col,
                 colp.SC, add.col,
                 #colp.SC,add.col,
                 colp.NC, add.col,
                 colp.NC_cons, add.col,
                 add.col, add.col)
    } else {
      plotdat <- c(plotdat,
                   OLS, space,
                   SC, space,
                   NC, space,
                   NC_cons)
      mycol <- c(mycol,
                 colp.OLS, add.col,
                 colp.SC, add.col,
                 #colp.SC,add.col,
                 colp.NC, add.col,
                 colp.NC_cons)
    }
  }
  
  if(removeoutlier == T){
    ## remove outlier
    ind <- lapply(plotdat, FUN = function(x){
      if(!is.na(x[1])){
        which(x > quantile(x, 1 - xbound) | x < quantile(x, xbound))
      }else{
        NULL
      }
    })
    for (i in 1:length(plotdat)) {
      if (!is.null(ind[[i]])) {
        plotdat[[i]] <- plotdat[[i]][-ind[[i]]]
      }
    }
  }
  
  return(list(SCNC = SCNC, plotdat = plotdat, mycol = mycol))
}

################################################################
## Generate input data for plotting for SC methods with  
## baseline covariate adjustment
## input:
## - rslt.all: data.frame combining results. Should contain
##            the following columns:
##   - t0: number of pre-treatment periods
##   - n.units: number of units
##   - SC_est, SC_OLS_est, NC_SC_est, NC_SC_constrained_est:
##     point estimates for different methods, adjusted for
##     baseline covariates;
##   - SC_est2, SC_OLS_est2, NC_SC_est2, NC_SC_constrained_est2:
##     unadjusted point estimates for different methods
################################################################

getplotdata <- function(rslt.all, removeoutlier = F){
  n.rslt <- length(rslt.all)
  n.t0 <- length(t0.all)
  mycol <- NULL
  plotdat <- SCNC <- list()
  
  for(n.units in unique(rslt.all$n.units)){
    SC <- SC.ignoreC <- OLS <- OLS.ignoreC <- 
      NC <- NC.ignoreC <- NC_cons <- NC_cons.ignoreC <- NULL
    for(t0 in unique(rslt.all$t0)){
      rowind <- (rslt.all$n.units == n.units) & (rslt.all$t0 == t0) 
      SC <- append(SC, list(rslt.all[rowind, "SC_est"]))
      SC.ignoreC <- append(SC.ignoreC, list(rslt.all[rowind, "SC_est2"]))
      OLS <- append(OLS, list(rslt.all[rowind, "SC_OLS_est"]))
      OLS.ignoreC <- append(OLS.ignoreC, list(rslt.all[rowind, "SC_OLS_est2"]))
      NC <- append(NC, list(rslt.all[rowind, "NC_SC_est"]))
      NC.ignoreC <- append(NC.ignoreC, list(rslt.all[rowind, "NC_SC_est2"]))
      NC_cons <- append(NC_cons, list(rslt.all[rowind, "NC_SC_constrained_est"]))
      NC_cons.ignoreC <- append(NC_cons.ignoreC, 
                                list(rslt.all[rowind, "NC_SC_constrained_est2"]))
    }
    SCNC <- c(SCNC, OLS, OLS.ignoreC, SC, NC.ignoreC, NC)
    space <- NA
    if(n.units != max(rslt.all$n.units)){
      ## organizing columns
      plotdat <- c(plotdat,
                   OLS.ignoreC, space, space,
                   OLS, space, space,
                   SC.ignoreC, space, space,
                   SC, space, space,
                   NC.ignoreC, space, space,
                   NC, space, space,
                   NC_cons.ignoreC, space, space,
                   NC_cons, space,
                   space, space)
      ## organizing colors
      mycol <- c(mycol,
                 colp.OLS, add.col, add.col,
                 colp.OLS, add.col, add.col,
                 colp.SC, add.col, add.col,
                 colp.SC, add.col, add.col,
                 colp.NC, add.col, add.col,
                 colp.NC, add.col, add.col,
                 colp.NC_cons, add.col, add.col,
                 colp.NC_cons, add.col,
                 add.col, add.col)
    } else {
      ## organizing columns
      plotdat <- c(plotdat,
                   OLS.ignoreC, space, space,
                   OLS, space, space,
                   SC.ignoreC, space,  space,
                   SC, space,  space,
                   NC.ignoreC, space,  space,
                   NC, space,  space,
                   NC_cons.ignoreC, space,  space,
                   NC_cons)
      ## organizing colors
      mycol <- c(mycol,
                 colp.OLS, add.col, add.col,
                 colp.OLS, add.col, add.col,
                 colp.SC, add.col, add.col,
                 colp.SC, add.col, add.col,
                 colp.NC, add.col, add.col,
                 colp.NC, add.col, add.col,
                 colp.NC_cons, add.col, add.col,
                 colp.NC_cons)
    }
  }
  
  if(removeoutlier == T){
    ind <- lapply(plotdat, FUN = function(x){
      if(!is.na(x[1])){
        which(x > quantile(x, 1 - xbound) | x < quantile(x, xbound))
      }else{
        NULL
      }
    })
    for (i in 1:length(plotdat)) {
      if (!is.null(ind[[i]])) {
        plotdat[[i]] <- plotdat[[i]][-ind[[i]]]
      }
    }
  }
  
  return(list(SCNC = SCNC, plotdat = plotdat, mycol = mycol))
}

########################################################
## Making plots to show the mean and SD of SC estimators
## with baseline covariate adjustment
## Input:
## - plotdat: output from getplotdata
## - mycol: colors for each column, contained in the 
##          output data frame of
## - myrange: deprecated
########################################################
makeplot2 <- function(plotdat, mycol, myrange){
  mean.bias <- sapply(plotdat, function(x) mean(x - true.beta))
  sd <- sapply(plotdat, function(x) sd(x - true.beta))
  LU <- lapply(plotdat, function(x){
    if(is.na(x[1])){
      return(c(NA,NA))
    }else{
      return(quantile(x - true.beta, probs = c(0.025,0.975)))
    }
  })
  
  
  breakat <- 1000; scale <- 0.2
  
  U <<- U.save <<- sapply(LU, function(x) x[2])
  L <<- L.save <<- sapply(LU, function(x) x[1])
  ind.U <- which(U > breakat)
  if (length(ind.U) > 0){
    U[ind.U] <- mean.bias[ind.U] + breakat + abs(U[ind.U] - mean.bias[ind.U] - breakat) * scale
  }
  ind.L <- which(L < (-1) * breakat)
  if(length(ind.L) > 0){
    L[ind.L] <- mean.bias[ind.L] - breakat - abs(mean.bias[ind.L] - L[ind.L] - breakat) * scale
  }
  ymin <- max(min(L, na.rm = T), -8.5)
  ymax <- min(max(U, na.rm = T), 8.5)
  
  ind.OLS_ignoreC <- c(1:3 + 41 * 0, 1:3 + 41 * 1, 1:3 + 41 * 2)
  ind.OLS <- c(6:8 + 41 * 0, 6:8 + 41 * 1, 6:8 + 41 * 2)
  ind.SC_ignoreC <- c(11:13 + 41 * 0, 11:13 + 41 * 1, 11:13 + 41 * 2)
  ind.SC <- c(16:18 + 41 * 0, 16:18 + 41 * 1, 16:18 + 41 * 2)
  ind.NC_ignoreC <- c(21:23 + 41 * 0, 21:23 + 41 * 1, 21:23 + 41 * 2)
  ind.NC <- c(26:28 + 41 * 0, 26:28 + 41 * 1, 26:28 + 41 * 2)
  ind.NC_cons_ignoreC <- c(31:33 + 41 * 0, 31:33 + 41 * 1, 31:33 + 41 * 2)
  ind.NC_cons <- c(36:38 + 41 * 0, 36:38 + 41 * 1, 36:38 + 41 * 2)
  
  plot.loc <- 1:length(plotdat)
  
  save.loc.OLS_ignoreC <- plot.loc[ind.OLS_ignoreC]
  plot.loc[ind.OLS_ignoreC] <- save.loc.OLS_ignoreC + rep(c(0, 2, 4), 3)
  plot.loc[ind.OLS] <- save.loc.OLS_ignoreC + rep(0.9 + c(0, 2, 4), 3)
  
  save.loc.SC_ignoreC <- plot.loc[ind.SC_ignoreC]
  plot.loc[ind.SC_ignoreC] <- save.loc.SC_ignoreC + rep(c(0, 2, 4), 3)
  plot.loc[ind.SC] <- save.loc.SC_ignoreC + rep(0.9 + c(0, 2, 4), 3)
  
  save.loc.NC_ignoreC <- plot.loc[ind.NC_ignoreC]
  plot.loc[ind.NC_ignoreC] <- save.loc.NC_ignoreC + rep(c(0, 2, 4),3)
  plot.loc[ind.NC] <- save.loc.NC_ignoreC + rep(0.9 + c(0, 2, 4), 3)
  
  save.loc.NC_cons_ignoreC <- plot.loc[ind.NC_cons_ignoreC]
  plot.loc[ind.NC_cons_ignoreC] <- save.loc.NC_cons_ignoreC + rep(c(0, 2, 4), 3)
  plot.loc[ind.NC_cons] <- save.loc.NC_cons_ignoreC + rep(0.9 + c(0, 2, 4), 3)
  
  ## specify colors
  mypal <- brewer.pal(4, "Set1")
  mycolors <- mycol
  mycolors[mycolors == colp.OLS] <- mypal[1]
  mycolors[mycolors == colp.SC] <- mypal[2]
  mycolors[mycolors == colp.NC] <- mypal[3]
  mycolors[mycolors == colp.NC_cons] <- mypal[4]
  
  ## change point types for the estimators w/o adjustments
  mycol2 <- mycol
  mycol2[ind.OLS_ignoreC] <- 15
  mycol2[ind.SC_ignoreC] <- 16
  mycol2[ind.NC_ignoreC] <- 17
  mycol2[ind.NC_cons_ignoreC] <- 18
  
  plot.new()
  plot.window(xlim = c(1, length(plotdat)), ylim = c(ymin, ymax))
  
  axis(side = 2, at = 0,line = 1, labels = "Bias", tick = F)
  
  abline(v = c(39.95, 80.95), col = "grey")
  axis(side = 1, at = c(-3.76, 39.95, 80.95, 124.7), labels = F)
  abline(h=0)
  
  box()
  points(x = plot.loc, y = mean.bias, pch = mycol2,
         ylim = c(ymin,ymax),cex = 0.9,
         col = mycolors)
  
  abline(h = 0)
  xHigh <- plot.loc#1:ncol(plotdat)
  xLow <- plot.loc#1:ncol(plotdat)
  yHigh <- U
  yLow <- L
  plot.lty <- rep(1, length(plotdat))
  
  ## specify the line types
  plot.lty[c(ind.OLS_ignoreC[c(1, 4, 7)], ind.OLS[c(1, 4, 7)],
             ind.SC_ignoreC[c(1, 4, 7)], ind.SC[c(1, 4, 7)],
             ind.NC_ignoreC[c(1, 4, 7)], ind.NC[c(1, 4, 7)],
             ind.NC_cons_ignoreC[c(1, 4, 7)], ind.NC_cons[c(1, 4, 7)])] <- 3
  plot.lty[c(ind.OLS_ignoreC[c(2, 5, 8)], ind.OLS[c(2, 5, 8)],
             ind.SC_ignoreC[c(2, 5, 8)], ind.SC[c(2, 5, 8)],
             ind.NC_ignoreC[c(2, 5, 8)], ind.NC[c(2, 5, 8)],
             ind.NC_cons_ignoreC[c(2, 5, 8)], ind.NC_cons[c(2, 5, 8)])] <- 4
  plot.lty[c(ind.OLS_ignoreC[c(3, 6, 9)], ind.OLS[c(3, 6, 9)],
             ind.SC_ignoreC[c(3, 6, 9)], ind.SC[c(3, 6, 9)],
             ind.NC_ignoreC[c(3, 6, 9)], ind.NC[c(3, 6, 9)],
             ind.NC_cons_ignoreC[c(3, 6, 9)], ind.NC_cons[c(3, 6, 9)])] <- 1
  
  arrows(xHigh, yHigh, xLow, yLow, angle = 90, length = 0, code = 3, lty = plot.lty,
         col = mycolors, lwd = 1.7)
  arrows(xHigh, yHigh, xLow, yHigh - 0.01, angle = 90, length = 0.05, code = 1, lty = 1,
         col = mycolors, lwd = 1.7)
  arrows(xHigh, yLow, xLow, yLow + 0.01, angle = 90, length = 0.05, code = 1, lty = 1,
         col = mycolors, lwd = 1.7)
  
  ind.ymin <- which.min(L); ind.ymax <- which.max(U)
  (y.label <- c(max(round(L[ind.ymin], 1), -8.5), -breakat, 0, breakat, min(round(U[ind.ymax], 1), 8.5)))
  (y.label.print <- c(max(round(L.save[ind.ymin], 1), -8.5), -breakat, 0, breakat, min(round(U.save[ind.ymax], 1), 8.5)))
  

  axis(side = 2, las = 1, cex.axis = 0.9, lwd = 0, lwd.ticks = 1,
       at = y.label,labels = y.label.print)
  

  text(x = c(4.5, 14.5, 24.5, 34.5,
             45.5, 55.5, 65.5, 75.5,
             86.5, 96.5, 106.5, 116.5),
       y = ymin - (ymax - ymin) * 0.09, srt = 0, adj = 0.5,
       labels = rep(c("OLS", "SC", "PI", "cPI"), 3), 
       xpd = TRUE, cex = 0.9,lwd = 1)

  xlabels = sapply(c(
    bquote(No.~ctrl~units==.(n.units.all[1]-1)),
    bquote(No.~ctrl~units==.(n.units.all[2]-1)),
    bquote(No.~ctrl~units==.(n.units.all[3]-1))
  ),as.expression)
  text(x = c(19.5, 60.5, 101.5),
       y = ymin - (ymax - ymin) * 0.18, srt = 0, adj = 0.5,
       labels = xlabels, xpd = TRUE, cex = 0.9, lwd = 1)
  if(length(ind.U) > 0){
    plotrix::axis.break(axis = 2, breakpos = c(breakat), pos = -0.5, bgcol = "white")
  }

}

########################################################
## Making plots to show the mean and SD of SC estimators
## without baseline covariate adjustment
## Input:
## - plotdat: output from getplotdata_nocov
## - mycol: colors for each column, contained in the 
##          output data frame of
## - myrange: deprecated
########################################################
makeplot2_nocov <- function(plotdat, mycol, myrange){
  mean.bias <- sapply(plotdat, function(x) mean(x - true.beta))
  sd <- sapply(plotdat, sd)
  LU <- lapply(plotdat, FUN = function(x){
    if (is.na(x[1])){
      return(c(NA, NA))
    } else {
      return(quantile(x - true.beta, probs = c(0.025, 0.975)))
    }
  })
  
  breakat = 1000; scale = 0.2
  
  
  U <<- U.save <<- sapply(LU, function(x) x[2])
  L <<- L.save <<- sapply(LU, function(x) x[1])
  ind.U = which(U > breakat)
  if(length(ind.U) > 0){
    U[ind.U] <- mean.bias[ind.U] + breakat + abs(U[ind.U] - mean.bias[ind.U] - breakat) * scale
  }
  ind.L <- which(L < (-1) * breakat)
  if(length(ind.L) > 0){
    L[ind.L] <- mean.bias[ind.L] - breakat - abs(mean.bias[ind.L] - L[ind.L] - breakat) * scale
  }
  ymin <- max(min(L, na.rm = T), -6.5)
  ymax <- min(max(U, na.rm = T), 6.5)
  
  ## specify colors
  mypal <- brewer.pal(4, "Set1")
  mycolors <- mycol
  mycolors[mycolors == colp.OLS] <- mypal[1]
  mycolors[mycolors == colp.SC] <- mypal[2]
  mycolors[mycolors == colp.NC] <- mypal[3]
  mycolors[mycolors == colp.NC_cons] <- mypal[4]
  
  
  plot.new()
  plot.window(xlim = c(1, length(plotdat)), ylim = c(ymin, ymax))
  
  axis(side = 2, at = 0,line = 1, labels = "Bias", tick = F)
  
  abline(v = c(17, 34), col = "grey")
  axis(side = 1, at = c(-1, 17, 34, 53), labels = F)
  
  abline(h=0)
  box()
  points(x = 1:length(plotdat),y = mean.bias, pch = mycol,
         ylim = c(ymin,ymax),cex = 0.9,
         col = mycolors)
  
  
  xHigh = 1:length(plotdat)
  xLow = 1:length(plotdat)
  
  yHigh = U
  yLow = L
  
  ## specify line types
  plot.lty <- mean.bias
  plot.lty[!is.na(plot.lty)] <- c(3, 4, 1)
  plot.lty[is.na(plot.lty)] <- 0
  
  arrows(xHigh, yHigh, xLow, yLow, angle = 90, length = 0, code = 3,
         lty = plot.lty, col = mycolors, lwd = 1.7)
  
  ## plot arrow heads with solid lines
  arrows(xHigh, yHigh, xLow, yHigh - 0.01, angle = 90, length = 0.05, code = 1,
         lty = 1, col = mycolors, lwd = 1.7)
  arrows(xHigh, yLow, xLow, yLow + 0.01, angle = 90, length = 0.05, code = 1,
         lty = 1, col = mycolors, lwd = 1.7)
  
  
  ind.ymin = which.min(L); ind.ymax = which.max(U)
  (y.label = c(max(round(L[ind.ymin],1), -6.5), -breakat, 0, breakat, min(round(U[ind.ymax],1), 6.5)))

  (y.label.print = c(max(round(L.save[ind.ymin],1), -6.5),-breakat, 0, breakat, min(round(U.save[ind.ymax],1), 6.5)))
  

  axis(side=2, las=1, cex.axis=0.9,lwd=0,lwd.ticks=1,
       at = y.label, labels = y.label.print)
  
  
  text(x = c(2, 6, 10, 14,
             20, 24, 28, 32,
             38, 42, 46, 50),
       y = ymin - (ymax - ymin) * 0.09, srt = 0, adj = 0.5,
       labels = rep(c("OLS", "SC", "PI", "cPI"), 3), 
       xpd = TRUE, cex = 0.9,lwd = 1)

  xlabels = sapply(c(
    bquote(No.~ctrl~units==.(n.units.all[1]-1)),
    bquote(No.~ctrl~units==.(n.units.all[2]-1)),
    bquote(No.~ctrl~units==.(n.units.all[3]-1))
  ),as.expression)
  text(x = c(8, 26, 44),
       y = ymin - (ymax - ymin) * 0.18, srt = 0, adj = 0.5,
       labels = xlabels, xpd = TRUE, cex = 0.9, lwd = 1)
  # legend("bottomleft",
  #        legend=sapply(c(
  #          bquote(T[0]==.(t0.all[1])),
  #          bquote(T[0]==.(t0.all[2])),
  #          bquote(T[0]==.(t0.all[3]))
  #        ),as.expression),
  #        pch=colp.SC,
  #        #fill=colp,
  #        bty="n",cex=0.6)
  if(length(ind.U)>0){plotrix::axis.break(axis=2,breakpos=c(breakat),pos=-0.5,bgcol="white")}
  if(length(ind.L)>0){plotrix::axis.break(axis=2,breakpos=c(-breakat),pos=-0.5,bgcol="white")}
}

########################################################################
#################################### bias plot #########################
########################################################################
true.beta <- 2; n.rep <- 500
t0.all <- c(30, 100, 200)
n.units.all <- c(1 + 4, 1 + 6, 1 + 10);
xbound <- 0.005

myrange <- 200

param_grid <- rbind(expand.grid(dist.epsilon = c("iid"),
                                dist.lambda = c("stationary", "nonstationary"),
                                U.setting = c("constrained", "unconstrained"),
                                addcov = c(TRUE, FALSE)),
                    expand.grid(dist.epsilon = c("AR"),
                                dist.lambda = c("stationary", "nonstationary"),
                                U.setting = c("constrained", "unconstrained"),
                                addcov = c(FALSE))) %>% unique()


for (ii in 1:nrow(param_grid)) {      
  addcov <- param_grid[ii, "addcov"]
  dist.epsilon <- param_grid[ii, "dist.epsilon"]
  dist.lambda <- param_grid[ii, "dist.lambda"]
  U.setting <- param_grid[ii, "U.setting"]
  
  if (dist.lambda == "nonstationary") {
    mysd <- 1.5
  } else {
    mysd <- 1.5
  }
  pdf(file = paste0("plot_LM_est/LM_v5",
                    "_sd_", mysd,
                    "_cov_", addcov,
                    "_epsilon_", dist.epsilon,
                    "_lambda_", dist.lambda,
                    "_U_", U.setting,
                    "_t0_", paste(t0.all, collapse="_"),
                    "_nunits_", paste(n.units.all, collapse="_"),
                    "_nrep", n.rep,
                    ".pdf"), width = (9 + 3) * 0.8, height = 5 * 0.8)
  par(mar = c(3.5, 3, 1, 1))

  
  rslt.all.combine <- NULL
  
  ## Reading the simulation results and organizing into a single data frame
  for (n.units in n.units.all) {
    for (t0 in t0.all) {
      for(batch in 1:10){
        tmp_list <- try(readRDS(file = paste0("results_linear_est/LM_v5",
                                              "_epsilon_", dist.epsilon,
                                              "_U_", U.setting,
                                              "_lambda_", dist.lambda,
                                              "_cov_", addcov,
                                              "_mysd_", mysd,
                                              "_t0_", t0,
                                              "_nunits_", n.units,
                                              "_nrep_", n.rep,
                                              "_batch_", batch,
                                              ".rds")),
                        silent=T)
        
        if (class(tmp_list)[[1]] != "try-error"){
          tmp_list <- tmp_list[unlist(sapply(tmp_list, 
                                             function(x) class(x) != "try-error"))]
          tmp <- data.frame(n.units = n.units,
                            t0 = t0,
                            SC_est = sapply(tmp_list, function(x) x["SC_est"]),
                            SC_OLS_est = sapply(tmp_list, function(x) x["SC_OLS_est"]),
                            NC_SC_est = sapply(tmp_list, function(x) x["NC_SC_est"]),
                            NC_SC_constrained_est = sapply(tmp_list, function(x) x["NC_SC_constrained_est"]),
                            SC_OLS_se = sapply(tmp_list, function(x) x["SC_OLS_se"]),
                            SC_OLS_se_hac = sapply(tmp_list, function(x) x["SC_OLS_se_hac"]),
                            NC_SC_se = sapply(tmp_list, function(x) x["NC_SC_se"]),
                            NC_SC_se_hac = sapply(tmp_list, function(x) x["NC_SC_se_hac"]),
                            SC_est2 = sapply(tmp_list, function(x) x["SC_est2"]),
                            SC_OLS_est2 = sapply(tmp_list, function(x) x["SC_OLS_est2"]),
                            NC_SC_est2 = sapply(tmp_list, function(x) x["NC_SC_est2"]),
                            NC_SC_constrained_est2 = sapply(tmp_list, function(x) x["NC_SC_constrained_est2"]),
                            SC_OLS_se2 = sapply(tmp_list, function(x) x["SC_OLS_se2"]),
                            SC_OLS_se_hac2 = sapply(tmp_list, function(x) x["SC_OLS_se_hac2"]),
                            NC_SC_se2 = sapply(tmp_list, function(x) x["NC_SC_se2"]),
                            NC_SC_se_hac2 = sapply(tmp_list, function(x) x["NC_SC_se_hac2"]))
          
          rslt.all.combine <- rbind(rslt.all.combine, tmp)
        }
      } 
    }
  }
  
  #for(argseeds in 1:5){
  rslt.all <- rslt.all.combine
  
  if(addcov == T){
    plotdat2 <- getplotdata(rslt.all, removeoutlier = T)
    #plotdat2=getplotdata_08242021(rslt.all,removeoutlier=F)
  }else{
    plotdat2 <- getplotdata_nocov(rslt.all, removeoutlier = T)
  }
  
  plotdat <- plotdat2$plotdat; mycol <- plotdat2$mycol
  if (addcov == T){
    makeplot2(plotdat, mycol, myrange)
  } else {
    makeplot2_nocov(plotdat, mycol, myrange)
  }
  dev.off()
  
}



par(mar = c(0, 0, 0, 0))
plot.new()
legend("left",#x=-0,y=0.58,
       legend=sapply(c(
         bquote("OLS, "~T[0]==.(t0.all[1])),
         bquote("OLS, "~T[0]==.(t0.all[2])),
         bquote("OLS, "~T[0]==.(t0.all[3])),
         ###

         bquote("PI (w/o adj), "~T[0]==.(t0.all[1])),
         bquote("PI (w/ adj), "~~~T[0]==.(t0.all[1])),
         
         
         bquote("PI (w/o adj), "~T[0]==.(t0.all[2])),
         bquote("PI (w/ adj), "~~~T[0]==.(t0.all[2])),
         
         bquote("PI (w/o adj), "~T[0]==.(t0.all[3])),
         bquote("PI (w/ adj), "~~~T[0]==.(t0.all[3]))
       ),as.expression),
       pch=c(colp.SC,rep(colp.NC,each=2)
       ),
       lty=c(rep(1,3),rep(c(3,1),3)),
       #fill=colp,
       ncol=1, inset=.1,
       # bty="n",
       cex=1.2)
dev.off()

png("legend_LM_nocov.png", width = 720, height = 120)
par(mar = c(0, 0, 0, 0))
plot.new()
legend(x = 0.5, y = 0.66, xjust = 0.5, yjust = 0.5,
       legend = c("OLS", "SC", "PI", "cPI"),
       col = brewer.pal(4, "Set1"),
       pch = c(0, 1, 2, 5), cex = 2.4, 
       lty = 0, lwd = 2.5, ncol = 4, bty = "n")
legend(x = 0.5, y = 0.33, xjust = 0.5, yjust = 0.5,
       lty = c(3, 4, 1), lwd = 2.5,
       legend = c(bquote(T[0] == .(t0.all[1])),
                  bquote(T[0] == .(t0.all[2])),
                  bquote(T[0] == .(t0.all[3]))),
       ncol = 3, cex = 2.4, bty = "n")
dev.off()

png("legend_LM_cov.png", width = 720, height = 240)
par(mar = c(0, 0, 0, 0))
plot.new()
text(0.03, 0.90, adj = 0, labels = c("Without covariate adjustment"),
     cex = 2.4)
legend(x = 0, y = 0.75, xjust = 0, yjust = 0.5,
       legend = c("OLS", "SC", "PI", "cPI"),
       col = brewer.pal(4, "Set1"),
       pch = 15:18, cex = 2.4, 
       lty = 0, lwd = 2.5, ncol = 4, bty = "n")
text(0.03, 0.55, adj = 0, labels = c("With covariate adjustment"),
     cex = 2.4)
legend(x = 0, y = 0.40, xjust = 0, yjust = 0.5,
       legend = c("OLS", "SC", "PI", "cPI"),
       col = brewer.pal(4, "Set1"),
       pch = c(0, 1, 2, 5), cex = 2.4, 
       lty = 0, lwd = 2.5, ncol = 4, bty = "n")
legend(x = 0, y = 0.10, xjust = 0, yjust = 0.5,
       lty = c(3, 4, 1), lwd = 2.5,
       legend = c(bquote(T[0] == .(t0.all[1])),
                  bquote(T[0] == .(t0.all[2])),
                  bquote(T[0] == .(t0.all[3]))),
       ncol = 3, cex = 2.4, bty = "n")
dev.off()

rue.beta &
#                                           SC.ate + qnorm(qq) * SC.se > true.beta ),
#                   NC.cover = as.numeric(NC.ate2 - qnorm(qq) * NC.se2.HAC < true.beta &
#                                           NC.ate2 + qnorm(qq) * NC.se2.HAC > true.beta ),
#                   CORR = dist.epsilon,
#                   addcov = addcov,
#                   nunits = n.units,
#                   t0 = t0)
#       rslt.all.combine <- rbind(rslt.all.combine, rslt)
#     } 
#   }
# }
# 
# rslt.all.summary <- rslt.all.combine %>% 
#   group_by(CORR, nunits, addcov, t0) %>%
#   summarize(SC.coverage = mean(SC.cover), NC.coverage = mean(NC.cover))
# 
# write.csv(rslt.all.summary, file = "tab_LM_coverage.csv")
# 
# 
