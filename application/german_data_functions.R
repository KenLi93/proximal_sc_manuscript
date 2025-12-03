predY0_each = function(para, data = data.all){
  ### get data and estimates
  X = cbind(data$X)
  Y = cbind(data$Y);C.Y=cbind(data$C.Y)
  Z = cbind(data$Z);C.Z=cbind(data$C.Z)
  W = cbind(data$W);C.W=cbind(data$C.W)
  theta.est = para[1:ncol(C.Y)]
  omega.est = para[ncol(C.Y) + 1:ncol(W)]
  ### predict Y0
  Womega = NULL
  for(i in predictor.name){
    Womega = cbind(
      Womega,
      C.W[,grep(i,colnames(C.W))] %*% omega.est
    )
  }
  # equivalent to below
  # Womega=NULL
  # for(t in 1:nrow(C.W)){
  #   ABC = matrix(C.W[t,],ncol=ncol(C.Y),byrow=T)
  #   Womega = rbind(
  #     Womega, t(omega.est)%*%ABC
  #   )
  # }
  Comega = C.Y - Womega
  if(ncol(C.Y) == 1){
    Comega = cbind(apply(Comega,1,mean))
  }
  Y0 = as.matrix(W) %*% c(omega.est) + Comega %*% c(theta.est)
  return(Y0)
}

getWYZ = function(dat=d, Y.name, W.name,Z.name,
                  outcome.name="gdp", unit.name="country", time.name="year"){
  Synth.Y = reshape(dat[,c(unit.name,outcome.name,time.name)],timevar=unit.name,direction="wide",idvar=time.name)
  Synth.Y = Synth.Y[,-grep(time.name,names(Synth.Y))]
  ind.Y = which(names(Synth.Y)==c(sapply(outcome.name,FUN=function(x){paste0(x,".",Y.name)})))
  ind.W = which(names(Synth.Y)%in%c(sapply(outcome.name,FUN=function(x){paste0(x,".",W.name)})))
  ind.Z = which(names(Synth.Y)%in%c(sapply(outcome.name,FUN=function(x){paste0(x,".",Z.name)})))
  Y = cbind(Synth.Y[,ind.Y])
  W = cbind(Synth.Y[,ind.W])
  Z = cbind(Synth.Y[,ind.Z])
  return(list(W=W,Y=Y,Z=Z))
}