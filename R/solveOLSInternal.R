#' @title solveOLSInternal
#'
#' @description This function solves or the unknown parameters
#' using ordinary least squares (OLS) without printing the output.
#' It returns the cell numbers, not the proportions (see solveOLS).
#'
#' @param S List output from trimData$sig (S)
#' @param B List output from trimData$bulk (B)
#'
#' @return Cell numbers
#'
#' @examples
#'
#' #load("data/Sig.RData")
#' #load("data/dataBulk.RData")
#' data('Sig', package = "DWLS")
#' data('dataBulk', package = "DWLS")
#'
#' trimmed <- trimData(Sig, dataBulk)
#' S <- trimmed$sig
#' B <- trimmed$bulk
#' solveOLSInternal(S, B)
#'
#' @export solveOLSInternal
#'
#' @importFrom dplyr "%>%"
#' @importFrom quadprog "solve.QP"
#'


solveOLSInternal<-function(S,B){
  D<-t(S)%*%S
  d<-t(S)%*%B
  A<-cbind(diag(dim(S)[2]))
  bzero<-c(rep(0,dim(S)[2]))
  solution<-solve.QP(D,d,A,bzero)$solution
  names(solution)<-colnames(S)
  return(solution)
}
