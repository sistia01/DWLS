#' @title solveSVR
#'
#' @description Performs a support vector regression (SVR). First, the data is
#' scalled then it solves for the SVR. An svm model is used with the following
#' specifications nu=0.5,scale = TRUE, type = "nu-regression",
#' kernel ="linear",cost = 1.
#'
#' @param S List output from trimData$sig (S)
#' @param B List output from trimData$bulk (B)
#'
#' @return Value (SVR)
#'
#' @examples
#' trimData(Sig, dataBulk)
#  S <- test$sig
#' B <- test$bulk
#' solveSVR(S, B)
#'
#' @export solveSVR
#'
#' @importFrom dplyr "%>%"



solveSVR<-function(S,B){
  ub=max(c(as.vector(S),B)) #upper bound
  lb=min(c(as.vector(S),B)) #lower bound
  Bs=((B-lb)/ub)*2-1
  Ss=((S-lb)/ub)*2-1
  model<-svm(Ss,Bs, nu=0.5,scale = TRUE, type = "nu-regression",kernel ="linear",cost = 1)
  coef <- t(model$coefs) %*% model$SV
  coef[which(coef<0)]<-0
  coef<-as.vector(coef)
  names(coef)<-colnames(S)
  print(round(coef/sum(coef),5))
  return(coef/sum(coef))
}
