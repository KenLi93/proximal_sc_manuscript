rm(list=ls())
placebo.test = F
library(dplyr)
library(foreign)
library(tidyr)
library(scpi)
library(Rsolnp)
library(RColorBrewer)

mycols <- brewer.pal(8, "Set1")
my.filepath = "application/"
source(paste0(my.filepath, "german_data_functions.R"))
source(paste0(my.filepath, "functions_conformal.R"))
source(paste0(my.filepath, "scpi_ate.R"))
source(paste0(my.filepath, "functions_conformal_avgeff.R"))
#Load data
germany_data <- as.data.frame(scpi_germany)


first.time = 1960
if (placebo.test == F) {
  last.ctrl.time = 1990
  last.trt.time = 2003
} else {
  last.ctrl.time = 1975
  last.trt.time = 1990
}

gdat <- germany_data %>%
  filter(year %in% c(first.time:last.trt.time))

## Set parameters for data preparation
id.var <- "country" # ID variable
time.var <- "year" # Time variable
period.pre <- first.time:last.ctrl.time # Pre - treatment period
period.post <- (last.ctrl.time+1):last.trt.time # Post - treatment period
myt0 = length(period.pre)
myt1 = length(period.post)
unit.tr <- "West Germany" # Treated unit
unit.co <- unique(germany_data$country)[-7] # Donor pool
pi.donors <- c("Austria", "Italy", "Japan", "Netherlands", "Switzerland", "USA")
pi.proxies <- setdiff(unique(germany_data$country), c(unit.tr, pi.donors))
outcome.var <- "gdp" # Outcome variable
constant <- TRUE # Include constant term
cointegrated.data <- TRUE # Cointegrated data


outcomes = getWYZ(dat = gdat, unit.tr, pi.donors, pi.proxies, outcome.var, id.var, time.var)
X <- c(rep(0,myt0),rep(1,myt1))
Y <<- as.matrix(outcomes$Y)
W <<- as.matrix(outcomes$W)
Z <<- as.matrix(outcomes$Z)
V <- cbind(W, Z)
V_Z <- Z[1:length(period.pre)] %*% t(Z[1:length(period.pre)]) / length(period.pre)
########################
#### do estimation ####
########################
data.pre = list(Y = Y[1:myt0], W = W[1:myt0,], Z = Z[1:myt0,])
data.all = list(X = X, Y = Y, W = W, Z = Z)


# Data preparation
SC_est <- SC_nocov(data.all, myt0)
ols_est <- OLS_nocov(data.all, myt0)

# proximal inference
cPI_est <- NC_constrained_nocov(data.all, myt0)
PI_est <- NC_nocov(data.all, myt0)

## get predicted Y0 from PI
ols.coef.pre = ols_est
Y0.ols = matrix(V %*% ols_est, ncol = 1)
ATT.ols <- mean(Y[myt0+1:myt1]-Y0.ols[myt0+1:myt1])
ols.est <- c(ols.coef.pre, ATT.ols)

SC.coef.pre = SC_est
Y0.SC = matrix(V %*% SC_est, ncol = 1)
ATT.SC <- mean(Y[myt0+1:myt1]-Y0.SC[myt0+1:myt1])
SC.est <- c(SC.coef.pre, ATT.SC)

PI.coef.pre = PI_est
Y0.PI = matrix(W %*% PI_est, ncol = 1)
ATT.PI <- mean(Y[myt0+1:myt1]-Y0.PI[myt0+1:myt1])
PI.est <- c(PI.coef.pre, ATT.PI)

cPI.coef.pre = cPI_est
Y0.cPI = matrix(W %*% cPI_est,
                ncol = 1)
ATT.cPI <- mean(Y[myt0+1:myt1]-Y0.cPI[myt0+1:myt1])
cPI.est <- c(cPI.coef.pre, ATT.cPI)


########################
#### Prep plot      ####
########################
## (1) get avg of all 38 states
outcomes.tmp = getWYZ(dat = gdat, Y.name = NULL, 
                      W.name = setdiff(unique(gdat$country), unit.tr), 
                      Z.name = NULL, outcome.var, "country", "year")
All_controls = outcomes.tmp$W; 
Y0.avg = cbind(apply(All_controls, 1, mean))
## change row names to time
time = unique(gdat$year)
row.names(Y0.avg) <- row.names(Y0.ols) <- 
  row.names(Y0.SC) <- row.names(Y0.PI) <- row.names(Y0.cPI) <- time
myrows = time




########################
#### End Prep plot  ####
########################
mygrey = "grey90"
lty.orig = 1; col.orig = 1;
lty.controls = 3; col.controls = mygrey; 
lty.average = 4; col.average = "grey50"
lty.ols = 1; col.ols = mycols[1]
lty.SC = 3; col.SC = mycols[2]
lty.PI = 4; col.PI = mycols[3]
lty.cPI = 5; col.cPI = mycols[4]

if(placebo.test==F){
  pdf(file=paste0(my.filepath,"SC_NewApplication.pdf"),w=6,h=5)
  myrange=c(0,33000)#range(c(0,Y,Y.PI.Sync))
}else{
  pdf(file=paste0(my.filepath,"SC_NewApplication_placebo.pdf"),w=6,h=5)
  myrange=c(0,31000)#range(c(0,Y,Y.PI.Sync))
}

par(mfrow=c(1,1),mar=c(0,0,0,0),mai=c(0.7,0.7,0.1,0.1),omi=c(0,0,0,0))
plot(x=time,Y * 1000, ylim=myrange,type="l",ylab="",xlab="",col=col.orig,lty=lty.orig,lwd=2)
axis(side=2,line=1,tick=F,at=mean(myrange),labels="Per Capita GDP")
axis(side=1,line=1,tick=F,at=mean(time),labels="Year")
# for(i in 1:ncol(All_controls)){
#   lines(x=time,y=All_controls[,i],lty=lty.controls,col=col.controls)
# }
# lines(x=time,y=Y,col=col.orig,lty=lty.orig)
abline(v=last.ctrl.time,lty=3)
# axis(side=1,line=-2,tick=F,at=last.ctrl.time,labels=last.ctrl.time,cex.axis=0.8,font=4)
if(last.ctrl.time==1990){
  axis(side=1,line=-2,tick=F,at=1975,labels="Pre-treatment",cex.axis=0.8,font=4)
  axis(side=1,line=-2,tick=F,at=1999,labels="Post-treatment",cex.axis=0.8,font=4)
  arrows(x0=1993, y0=2500, x1 = 1990, y1=2500, length=0.1)
  text(x=1997,y=2500,labels="Reunification",cex=0.8)
}else{
  axis(side=1,line=-2,tick=F,at=1967,labels="Pre-treatment",cex.axis=0.8,font=4)
  axis(side=1,line=-2,tick=F,at=1985,labels="Post-treatment",cex.axis=0.8,font=4)
  arrows(x0=1977.2, y0=2500, x1 = 1975, y1=2500, length=0.1)
  text(x=1981.8,y=2500,labels="Placebo Reunification",cex=0.8)
}
lines(x=time,y=Y0.avg * 1000, lty=lty.average,col=col.average,lwd=2)
lines(x=time,y=Y0.ols * 1000, lty=lty.ols,col=col.ols,lwd=2)
lines(x=time,y=Y0.SC * 1000, lty=lty.SC,col=col.SC,lwd=2)
lines(x=time,y=Y0.PI * 1000, lty=lty.PI,col=col.PI,lwd=2)
lines(x=time,y=Y0.cPI * 1000, lty=lty.cPI,col=col.cPI,lwd=2)

legend("topleft",
       legend = c("West Germany (treated)",
                  # "Control countries (N=16)",
                  "Average of 16 control countries",
                  "Syn West Germany (OLS)",
                  "Syn West Germany (SC)",
                  "Syn West Germany (PI)",
                  "Syn West Germany (cPI)"),
       lty=c(lty.orig,#lty.controls,
             lty.average,
             lty.ols, lty.SC, 
             lty.PI, lty.cPI),
       col=c(col.orig,#col.controls,
             col.average, 
             mycols[1:4]),
       bty="n",lwd=2,cex=1)
dev.off()

if (placebo.test == T) {
  MAE_OLS <- mean(abs(Y[(myt0 + 1):(myt0 + myt1)] - Y0.ols[(myt0 + 1):(myt0 + myt1)])) ## 3.316175
  MAE_SC <- mean(abs(Y[(myt0 + 1):(myt0 + myt1)] - Y0.SC[(myt0 + 1):(myt0 + myt1)])) ## 1.350597
  MAE_PI <- mean(abs(Y[(myt0 + 1):(myt0 + myt1)] - Y0.PI[(myt0 + 1):(myt0 + myt1)])) ## 1.115383
  MAE_cPI <- mean(abs(Y[(myt0 + 1):(myt0 + myt1)] - Y0.cPI[(myt0 + 1):(myt0 + myt1)])) ## 0.163842
}

if (placebo.test == F) {
  crude_grid <- replicate(myt1, seq(-3.8, 0.9, 0.1), simplify = F)
  cPI_crude <- ConformalPointwiseCI(data.all, myt0, myt1, crude_grid, 
                                    output = "both", alpha = 0.1, 
                                    method = "NC-constrained")
  
  cPI_grid <- lapply(1:myt1, function(tt) {
    seq(cPI_crude$lb[tt] - 0.1, cPI_crude$ub[tt] + 0.1, 0.01)
  })
  
  cPI_pointwise <- ConformalPointwiseCI(data.all, myt0, myt1, cPI_grid, 
                                        output = "both", alpha = 0.1, 
                                        method = "NC-constrained")
  
  SC_crude <- ConformalPointwiseCI(data.all, myt0, myt1, crude_grid, 
                                       output = "both", alpha = 0.1, 
                                       method = "SC")
  
  SC_grid <- lapply(1:myt1, function(tt) {
    seq(cPI_crude$lb[tt] - 0.1, SC_crude$ub[tt] + 0.1, 0.01)
  })
  
  SC_pointwise <- ConformalPointwiseCI(data.all, myt0, myt1, SC_grid, 
                                        output = "both", alpha = 0.1, 
                                        method = "SC")
  
  
  
  #### plot differences
  pdf(file = paste0(my.filepath,"SC_conformal_PI.pdf"),w = 6,h = 5)
  myrange = c(-4000, 1000)
  
  
  plot(x = time, y = rep(0, length(time)), ylim = myrange,
       type = "l", ylab = "", xlab = "", col = col.orig, lty = lty.orig, lwd=2)
  abline(h = seq(-4000, 1000, 500), col = mygrey)
  axis(side=2,line=1,tick=F,at=mean(myrange),labels="Per Capita GDP")
  axis(side=1,line=1,tick=F,at=mean(time),labels="Year")
  # for(i in 1:ncol(All_controls)){
  #   lines(x=time,y=All_controls[,i],lty=lty.controls,col=col.controls)
  # }
  # lines(x=time,y=Y,col=col.orig,lty=lty.orig)
  abline(v=last.ctrl.time,lty=3)
  # axis(side=1,line=-2,tick=F,at=last.ctrl.time,labels=last.ctrl.time,cex.axis=0.8,font=4)
  
  axis(side=1,line=-2,tick=F,at=1975,labels="Pre-treatment",cex.axis=0.8,font=4)
  axis(side=1,line=-2,tick=F,at=1999,labels="Post-treatment",cex.axis=0.8,font=4)
  arrows(x0=1986, y0 = -3000, x1 = 1990, y1 = -3000, length=0.1)
  text(x=1982,y= -3000,labels="Reunification",cex=0.8)
  lines(x = time, y = (Y - Y0.SC) * 1000, lty = 1, col = col.SC, lwd = 2)
  # points(x = time[(myt0 + 1):(myt0 + myt1)], 
  #        y = (Y - Y0.SC)[(myt0 + 1):(myt0 + myt1)] * 1000, 
  #        pch = 4, col = col.SC, lwd = 2)
  arrows(x0 = time[(myt0 + 1):(myt0 + myt1)],
         y0 = SC_pointwise$lb * 1000,
         x1 =  time[(myt0 + 1):(myt0 + myt1)],
         y1 = SC_pointwise$ub * 1000,
         code = 3, angle = 90, length = 0.04,
         lwd = 2, col = col.SC)
  dev.off()
  
  pdf(file = paste0(my.filepath,"cPI_conformal_PI.pdf"),w = 6,h = 5)
  myrange = c(-4000, 1000)
  
  
  plot(x = time, y = rep(0, length(time)), ylim = myrange,
       type = "l", ylab = "", xlab = "", col = col.orig, lty = lty.orig, lwd=2)
  abline(h = seq(-4000, 1000, 500), col = mygrey)
  axis(side=2,line=1,tick=F,at=mean(myrange),labels="Per Capita GDP")
  axis(side=1,line=1,tick=F,at=mean(time),labels="Year")
  # for(i in 1:ncol(All_controls)){
  #   lines(x=time,y=All_controls[,i],lty=lty.controls,col=col.controls)
  # }
  # lines(x=time,y=Y,col=col.orig,lty=lty.orig)
  abline(v=last.ctrl.time,lty=3)
  # axis(side=1,line=-2,tick=F,at=last.ctrl.time,labels=last.ctrl.time,cex.axis=0.8,font=4)
  
  axis(side=1,line=-2,tick=F,at=1975,labels="Pre-treatment",cex.axis=0.8,font=4)
  axis(side=1,line=-2,tick=F,at=1999,labels="Post-treatment",cex.axis=0.8,font=4)
  arrows(x0=1986, y0 = -3000, x1 = 1990, y1 = -3000, length=0.1)
  text(x=1982,y= -3000,labels="Reunification",cex=0.8)
  lines(x = time, y = (Y - Y0.cPI) * 1000, lty = 1, col = col.cPI, lwd = 2)
  # points(x = time[(myt0 + 1):(myt0 + myt1)], 
  #        y = (Y - Y0.SC)[(myt0 + 1):(myt0 + myt1)] * 1000, 
  #        pch = 4, col = col.SC, lwd = 2)
  arrows(x0 = time[(myt0 + 1):(myt0 + myt1)],
         y0 = cPI_pointwise$lb * 1000,
         x1 =  time[(myt0 + 1):(myt0 + myt1)],
         y1 = cPI_pointwise$ub * 1000,
         code = 3, angle = 90, length = 0.04,
         lwd = 2, col = col.cPI)
  dev.off()
  
  
  
  scpi_ate_sc_int(data.all, myt0, method = "ols")
  scpi_ate_sc_int(data.all, myt0, method = "simplex")
  scpi_ate_nc_int(data.all, myt0, method = "ols")
  scpi_ate_nc_int(data.all, myt0, method = "simplex")
  
  OLS_nocov_gmm(data.all)
  NC_nocov_gmm(data.all)
}


