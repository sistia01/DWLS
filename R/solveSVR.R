solveSVR<-function(S,B){
  #scaling
  ub=max(c(as.vector(S),B)) #upper bound
  lb=min(c(as.vector(S),B)) #lower bound
  Bs=((B-lb)/ub)*2-1
  Ss=((S-lb)/ub)*2-1

  #perform SVR
  model<-svm(Ss,Bs, nu=0.5,scale = TRUE, type = "nu-regression",kernel ="linear",cost = 1)
  coef <- t(model$coefs) %*% model$SV
  coef[which(coef<0)]<-0
  coef<-as.vector(coef)
  names(coef)<-colnames(S)
  print(round(coef/sum(coef),5))
  return(coef/sum(coef))
}
