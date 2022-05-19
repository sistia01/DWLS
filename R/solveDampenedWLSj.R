#' @title solveDampenedWLSj
#'
#' @description Solve dampened weighted least squares given a dampening constant.
#'
#'
#' Note: The function uses solveDampenedWLS() and findDampeningConstant().
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
#' \donttest{
#' #Sig
#' #url <- "https://github.com/sistia01/DWLS/raw/main/inst/extdata/Sig.RData"
#' #dest <- "data/Sig.RData"
#' #download.file(url, dest)
#' #load("data/Sig.RData")
#' load(system.file("extdata", "Sig.RData", package = "DWLS"))
#'
#' #dataBulk
#' #url <- "https://github.com/sistia01/DWLS/raw/main/inst/extdata/dataBulk.RData"
#' #dest <- "data/dataBulk.RData"
#' #download.file(url, dest)
#' #load("data/dataBulk.RData")
#' load(system.file("extdata", "dataBulk.RData", package = "DWLS"))
#'
#' trimmed <- trimData(Sig, dataBulk)
#' S <- trimmed$sig
#' B <- trimmed$bulk
#' solution <- solveOLSInternal(S,B)
#' j <- findDampeningConstant(S,B,solution)
#' goldStandard <- solveOLSInternal(S,B)
#' solveDampenedWLSj(S,B,goldStandard,j)
#'}
#' @export solveDampenedWLSj
#'
#' @importFrom dplyr "%>%"
#' @importFrom quadprog "solve.QP"

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
