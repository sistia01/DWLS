#find a dampening constant for the weights using cross-validation
findDampeningConstant<-function(S,B,goldStandard){
  solutionsSd<-NULL
  #goldStandard is used to define the weights
  sol<-goldStandard
  ws<-as.vector((1/(S%*%sol))^2)
  wsScaled<-ws/min(ws)
  wsScaledMinusInf<-wsScaled
  #ignore infinite weights
  if(max(wsScaled)=="Inf"){
    wsScaledMinusInf<-wsScaled[-which(wsScaled=="Inf")]
  }
  #try multiple values of the dampening constant (multiplier)
  #for each, calculate the variance of the dampened weighted solution for a subset of genes
  for (j in 1:ceiling(log2(max(wsScaledMinusInf)))){
    multiplier<-1*2^(j-1)
    wsDampened<-wsScaled
    wsDampened[which(wsScaled>multiplier)]<-multiplier
    solutions<-NULL
    seeds<-c(1:100)
    for (i in 1:100){
      set.seed(seeds[i]) #make nondeterministic
      subset<-sample(length(ws),size=length(ws)*0.5) #randomly select half of gene set
      #solve dampened weighted least squares for subset
      fit = lm (B[subset] ~ -1+S[subset,],weights=wsDampened[subset])
      sol<-fit$coef*sum(goldStandard)/sum(fit$coef)
      solutions<-cbind(solutions,sol)
    }
    solutionsSd<-cbind(solutionsSd,apply(solutions,1,sd))
  }
  #choose dampening constant that results in least cross-validation variance
  j<-which.min(colMeans(solutionsSd^2))
  return(j)
}
