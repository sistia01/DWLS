#' @title solveOLS
#'
#' @description This function solves using OLS. It is constrained such that cell type numbers are greater than 0
#'
#' @param Signature_Matrix
#' @param bulkdata
#'
#' @return NULL
#'
#' @examples
#'
#' @export
#'
#' @importFrom dplyr "%>%"
#'

solveOLS<-function(S,B){
  D<-t(S)%*%S
  d<-t(S)%*%B
  A<-cbind(diag(dim(S)[2]))
  bzero<-c(rep(0,dim(S)[2]))
  solution<-solve.QP(D,d,A,bzero)$solution
  names(solution)<-colnames(S)
  print(round(solution/sum(solution),5))
  return(solution/sum(solution))
}
