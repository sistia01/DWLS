#' @title solveDampenedWLSj
#'
#' @description Solve dampened weighted least squares given a dampening constant.
#'
#'
#' Note: The function uses solveDampenedWLSj() and findDampeningConstant().
#'
#' @param S List output from trimData$sig (S)
#' @param B List output from trimData$bulk (B)
#' @param goldStandard Starting point for the weights, this can be determined
#' using solveOLSInternal(S,B)
#' @param j The dampening constant, this can be determined using
#' findDampeningConstant(S,B,goldStandard)
#'
#' @return value (Dampened weighted least squares estimation values)
#'
#' @examples
#' trimData(Sig, dataBulk)
#  S <- test$sig
#' B <- test$bulk
#' solution <- solveOLSInternal(S,B)
#' j <- findDampeningConstant(S,B,solution)
#' solveDampenedWLSj(S,B,goldStandard,j)
#'
#' @export
#'
#' @importFrom dplyr "%>%"


solveDampenedWLSj<-function(S,B,goldStandard,j){
  multiplier<-1*2^(j-1)
  sol<-goldStandard
  ws<-as.vector((1/(S%*%sol))^2)
  wsScaled<-ws/min(ws)
  wsDampened<-wsScaled
  wsDampened[which(wsScaled>multiplier)]<-multiplier
  W<-diag(wsDampened)
  D<-t(S)%*%W%*%S
  d<- t(S)%*%W%*%B
  A<-cbind(diag(dim(S)[2]))
  bzero<-c(rep(0,dim(S)[2]))
  sc <- norm(D,"2")
  solution<-solve.QP(D/sc,d/sc,A,bzero)$solution
  names(solution)<-colnames(S)
  return(solution)
}
