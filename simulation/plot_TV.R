rm(list=ls())

library(dplyr)
library(RColorBrewer)

colp.OLS <- rep(0, 3)
colp.SC <- rep(1, 3)#c("gray70","gray90","gray100")
colp.NC <- rep(2, 3)#c("gray70","gray90","gray100")
colp.NC_cons <- rep(5, 3)

add.col <- NA #c("gray10")
getplotdata_tv1 <- function(rslt.all, removeoutlier = F){
  # c(SC.ate=SC.ate, ## unconstrained OLS
  #   SE.ate2=SC.ate.constrained, ## constrained OLS
  #   NC.ate=NC.ate,##NC ignore C
  #   NC.ate2=NC.ate2,##NC adjust C
  #   SC.se =SC.se, ## unconstrained OLS
  #   NC.se =NC.se,NC.se.HAC=NC.se.HAC, ##NC ignore C
  #   NC.se2=NC.se2,NC.se2.HAC=NC.se2.HAC##NC adjust C
  # )
  n.rslt <- length(rslt.all)
  n.t0 <- length(t0.all)
  plotdat <- SCNC <- NULL
  mycol <- NULL
  
  for(n.units in unique(rslt.all$n.units)){
    SC <- OLS <- NC <- NC_cons <- list()
    for(t0 in unique(rslt.all$t0)){
      rowind <- (rslt.all$n.units == n.units) & (rslt.all$t0 == t0) 
      SC <- c(SC, list(rslt.all[rowind, "SC_est1"]))
      OLS <- c(OLS, list(rslt.all[rowind, "SC_OLS_est1"]))
      NC <- c(NC, list(rslt.all[rowind, "NC_SC_est1"]))
      NC_cons <- c(NC_cons, list(rslt.all[rowind, "NC_SC_constrained_est1"]))
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

getplotdata_tv2 <- function(rslt.all, removeoutlier = F){
  # c(SC.ate=SC.ate, ## unconstrained OLS
  #   SE.ate2=SC.ate.constrained, ## constrained OLS
  #   NC.ate=NC.ate,##NC ignore C
  #   NC.ate2=NC.ate2,##NC adjust C
  #   SC.se =SC.se, ## unconstrained OLS
  #   NC.se =NC.se,NC.se.HAC=NC.se.HAC, ##NC ignore C
  #   NC.se2=NC.se2,NC.se2.HAC=NC.se2.HAC##NC adjust C
  # )
  n.rslt <- length(rslt.all)
  n.t0 <- length(t0.all)
  plotdat <- SCNC <- NULL
  mycol <- NULL
  
  for(n.units in unique(rslt.all$n.units)){
    SC <- OLS <- NC <- NC_cons <- list()
    for(t0 in unique(rslt.all$t0)){
      rowind <- (rslt.all$n.units == n.units) & (rslt.all$t0 == t0) 
      SC <- c(SC, list(rslt.all[rowind, "SC_est2"]))
      OLS <- c(OLS, list(rslt.all[rowind, "SC_OLS_est2"]))
      NC <- c(NC, list(rslt.all[rowind, "NC_SC_est2"]))
      NC_cons <- c(NC_cons, list(rslt.all[rowind, "NC_SC_constrained_est2"]))
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


makeplot_tv <- function(plotdat, cord, mycol, myrange, ymin = -6.5, ymax = 6.5){
  mean.bias <- sapply(plotdat, function(x) mean(x - true.beta[cord]))
  sd <- sapply(plotdat, sd)
  LU <- lapply(plotdat, FUN = function(x){
    if (is.na(x[1])){
      return(c(NA, NA))
    } else {
      return(c(mean(x) - sd(x) - true.beta[cord], mean(x) + sd(x) - true.beta[cord]))
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
  ymin <- max(min(L, na.rm = T), ymin)
  # ymin=-myrange
  ymax <- min(max(U, na.rm = T), ymax)
  # ymax=myrange
  
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
  # yHigh = mean.bias-qnorm(0.975)*sd
  # yLow = mean.bias+qnorm(0.975)*sd
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
  # y.label[which.min(y.label)]=round(ymin,2)
  # y.label[which.max(y.label)]=round(ymax,2)
  (y.label.print = c(max(round(L.save[ind.ymin],1), -6.5),-breakat, 0, breakat, min(round(U.save[ind.ymax],1), 6.5)))
  
  # (y.label=c(seq(round(ymin,1),breakat,round(ymax,1),length.out=3),0))
  # y.label=y.label[order(y.label)]
  # # y.label[which.min(y.label)]=round(ymin,2)
  # # y.label[which.max(y.label)]=round(ymax,2)
  # y.label.print=y.label
  # # y.label.print[1]=round(min(LU,na.rm=T),1)
  # # y.label.print[5]=round(max(LU,na.rm=T),1)
  axis(side=2, las=1, cex.axis=0.9,lwd=0,lwd.ticks=1,
       at = y.label, labels = y.label.print)
  
  
  text(x = c(2, 6, 10, 14,
             20, 24, 28, 32,
             38, 42, 46, 50),
       y = ymin - (ymax - ymin) * 0.09, srt = 0, adj = 0.5,
       labels = rep(c("OLS", "SC", "PI", "cPI"), 3), 
       xpd = TRUE, cex = 0.9,lwd = 1)
  # text(x=c(3.3,7.4,13.3,17.4,23.3,27.4)-0.5,
  #      y=ymin-(ymax-ymin)*0.09,srt=0,adj=1,
  #      labels=rep(c("OLS","PI"),3),xpd=TRUE,cex=0.9,lwd=1)
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
# fixed.lambda <- F; 
true.beta <- c(2, 0.4); n.rep <- 500
t0.all <- c(100, 150, 200)
#n.units.all <- c(1 + 4, 1 + 10, 1 + 20);
n.units.all <- c(1 + 2, 1 + 4, 1 + 6);
xbound <- 0.01
# n.units.all=c(1+10,1+20,1+30);


# U_is_fixed <- T; 
myrange <- 200

# if (mysd == 0.5) { myrange <- 100 }

param_grid <- expand.grid(dist.epsilon = c("iid"),
                          dist.lambda = c("stationary"),
                          U.setting = c("constrained", "unconstrained"))

# param_grid <- expand.grid(dist.epsilon = c("iid"),
#                           dist.lambda = c("stationary", "nonstationary"),
#                           U.setting = c("constrained"),
#                           addcov = c(FALSE)) 





for (ii in 1:nrow(param_grid)) {      
  dist.epsilon <- param_grid[ii, "dist.epsilon"]
  dist.lambda <- param_grid[ii, "dist.lambda"]
  U.setting <- param_grid[ii, "U.setting"]
  
  
  mysd <- 1
  
    # for(arg3 in c(1,2,4,6,7)[1]){
  
  
  rslt.all.combine <- NULL
  
  for (n.units in n.units.all) {
    for (t0 in t0.all) {
      for(batch in 1:10){
        tmp_list <- try(readRDS(file = paste0("results_tv/LM",
                                              "_epsilon_", dist.epsilon,
                                              "_U_", U.setting,
                                              "_lambda_", dist.lambda,
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
                            SC_est1 = sapply(tmp_list, function(x) x$SC_est[1]),
                            SC_OLS_est1 = sapply(tmp_list, function(x) x$SC_OLS_est[1]),
                            NC_SC_est1 = sapply(tmp_list, function(x) x$NC_SC_est[1]),
                            NC_SC_constrained_est1 = sapply(tmp_list, function(x) x$NC_SC_constrained_est[1]),
                            SC_est2 = sapply(tmp_list, function(x) x$SC_est[2]),
                            SC_OLS_est2 = sapply(tmp_list, function(x) x$SC_OLS_est[2]),
                            NC_SC_est2 = sapply(tmp_list, function(x) x$NC_SC_est[2]),
                            NC_SC_constrained_est2 = sapply(tmp_list, function(x) x$NC_SC_constrained_est[2]),
                            SC_OLS_se1 = sapply(tmp_list, function(x) x$SC_OLS_se[1]),
                            SC_OLS_se_hac1 = sapply(tmp_list, function(x) x$SC_OLS_se_hac[1]),
                            NC_SC_se1 = sapply(tmp_list, function(x) x$NC_SC_se[1]),
                            NC_SC_se_hac1 = sapply(tmp_list, function(x) x$NC_SC_se_hac[1]),
                            SC_OLS_se2 = sapply(tmp_list, function(x) x$SC_OLS_se[2]),
                            SC_OLS_se_hac2 = sapply(tmp_list, function(x) x$SC_OLS_se_hac[2]),
                            NC_SC_se2 = sapply(tmp_list, function(x) x$NC_SC_se[2]),
                            NC_SC_se_hac2 = sapply(tmp_list, function(x) x$NC_SC_se_hac[2]))
          
          rslt.all.combine <- rbind(rslt.all.combine, tmp)
        }
      } 
    }
  }
  
  #for(argseeds in 1:5){
  rslt.all <- rslt.all.combine
  
  plotdat_obj1 <- getplotdata_tv1(rslt.all, removeoutlier = T)
  plotdat1 <- plotdat_obj1$plotdat; mycol1 <- plotdat_obj1$mycol
  
  plotdat_obj2 <- getplotdata_tv2(rslt.all, removeoutlier = T)
  plotdat2 <- plotdat_obj2$plotdat; mycol2 <- plotdat_obj2$mycol
  
  
  pdf(file = paste0("plot_TV/beta1",
                    "_sd_", mysd,
                    "_lambda_", dist.lambda,
                    "_U_", U.setting,
                    "_t0_", paste(t0.all, collapse="_"),
                    "_nunits_", paste(n.units.all, collapse="_"),
                    "_nrep", n.rep,
                    # "_U_is_fixed", U_is_fixed,
                    ".pdf"), width = (9 + 3) * 0.8, height = 5 * 0.8)
  par(mar = c(3.5, 3, 1, 1))

  makeplot_tv(plotdat1, cord = 1, mycol1, myrange, ymin = -2.5, ymax = 2.5)
  dev.off()
  
  
  pdf(file = paste0("plot_TV/beta2",
                    "_sd_", mysd,
                    "_lambda_", dist.lambda,
                    "_U_", U.setting,
                    "_t0_", paste(t0.all, collapse="_"),
                    "_nunits_", paste(n.units.all, collapse="_"),
                    "_nrep", n.rep,
                    # "_U_is_fixed", U_is_fixed,
                    ".pdf"), width = (9 + 3) * 0.8, height = 5 * 0.8)
  par(mar = c(3.5, 3, 1, 1))
  
  makeplot_tv(plotdat2, cord = 2, mycol1, myrange, ymin = -2.5, ymax = 2.5)
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
         # bquote("OLS, "~T[0]==.(t0.all[1])~", "~T==.(2*t0.all[1])),
         # bquote("OLS, "~T[0]==.(t0.all[2])~", "~T==.(2*t0.all[2])),
         # bquote("OLS, "~T[0]==.(t0.all[3])~", "~T==.(2*t0.all[3]))
         bquote("PI (w/o adj), "~T[0]==.(t0.all[1])),
         bquote("PI (w/ adj), "~~~T[0]==.(t0.all[1])),
         
         
         bquote("PI (w/o adj), "~T[0]==.(t0.all[2])),
         bquote("PI (w/ adj), "~~~T[0]==.(t0.all[2])),
         
         bquote("PI (w/o adj), "~T[0]==.(t0.all[3])),
         bquote("PI (w/ adj), "~~~T[0]==.(t0.all[3]))
       ),as.expression),
       pch=c(colp.SC,rep(colp.NC,each=2)
             # colp.SC[1],colp.NC[1],
             # colp.SC[2],colp.NC[2],
             # colp.SC[3],colp.NC[3]
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




# ########################################################################
# #################################### 95%CI ####################################
# ########################################################################
# rm(list=ls())
# 
# qq <- 1 - (1 - 0.95) / 2 ##95%CI
# 
# 
# param_grid <- expand.grid(n.rep = 200,
#                           t0 = c(50, 100, 200),
#                           n.units = c(1 + 2, 1 + 10, 1 + 20),
#                           addcov = c(T, F),
#                           dist.epsilon = c("iid", "AR"))
# true.beta <- 2
# rslt.all.combine <- NULL
# for (ii in 1:nrow(param_grid)) {
#   n.rep <- param_grid[ii, "n.rep"]
#   t0 <- param_grid[ii, "t0"]
#   n.units <- param_grid[ii, "n.units"]
#   dist.epsilon <- param_grid[ii, "dist.epsilon"]
#   addcov <- param_grid[ii, "addcov"]
#   
#   for (batch in 1:10) {
#     if (dist.epsilon == "iid") {
#       rslt <- readRDS(file = paste0("results/LM_",
#                                     "epsilon_", dist.epsilon,
#                                     "cov_", addcov,
#                                     "t0_", t0,
#                                     "nunits_", n.units,
#                                     "nrep", n.rep,
#                                     "batch", batch,
#                                     ".rds")) %>% 
#         as.data.frame %>%
#         transmute(SC.cover = as.numeric(SC.ate - qnorm(qq) * SC.se < true.beta &
#                                           SC.ate + qnorm(qq) * SC.se > true.beta ),
#                   NC.cover = as.numeric(NC.ate2 - qnorm(qq) * NC.se2 < true.beta &
#                                           NC.ate2 + qnorm(qq) * NC.se2 > true.beta ),
#                   CORR = dist.epsilon,
#                   addcov = addcov,
#                   nunits = n.units,
#                   t0 = t0)
#       rslt.all.combine <- rbind(rslt.all.combine, rslt)
#     } else {
#       rslt <- readRDS(file = paste0("results/LM_",
#                                     "epsilon_", dist.epsilon,
#                                     "cov_", addcov,
#                                     "t0_", t0,
#                                     "nunits_", n.units,
#                                     "nrep", n.rep,
#                                     "batch", batch,
#                                     ".rds")) %>% 
#         as.data.frame %>%
#         transmute(SC.cover = as.numeric(SC.ate - qnorm(qq) * SC.se < true.beta &
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
