#' @title solveDampenedWLS
#'
#' @description Dampened weighted least squares (DLWS) is an estimation
#' method for gene expression deconvolution, in which the cell-type
#' composition of a bulk RNA-seq data set is computationally inferred.
#' This method corrects common biases towards cell types that are
#' characterized by highly expressed genes and/or are highly prevalent,
#' to provide accurate detection across diverse cell types. To begin,
#' the user must input a bulk RNA-seq data set, along with a labeled
#' representative single-cell RNA-seq data set that will serve to generate
#' cell-type-specific gene expression profiles. Ideally, the single-cell
#' data set will contain cells from all cell types that may be found in the
#' bulk data. DWLS will return the cell-type composition of the bulk data.
#' First, solve OLS then use the solution to find a starting point for the
#' weights. Next, the dampened weighted least squares is performed. The weights
#' are iterated until convergence then the dampening constant for weights is
#' found using cross-validation (with decreasing step size for convergence).
#'
#' Note: The function uses solveDampenedWLSj() and findDampeningConstant().
#'
#' @param S List output from trimData$sig (S)
#' @param B List output from trimData$bulk (B)
#'
#' @return value (Dampened weighted least squares estimation values)
#'
#' @examples
#' trimData(Sig, dataBulk)
#  S <- test$sig
#' B <- test$bulk
#' solveDampenedWLS(S, B)
#'
#' @export solveDampenedWLS
#'
#' @importFrom dplyr "%>%"


solveDampenedWLS<-function(S,B){
  solution<-solveOLSInternal(S,B)
  iterations<-0
  changes<-c()
  j<-findDampeningConstant(S,B,solution)
  change<-1
  while(change>.01 & iterations<1000){
    newsolution<-solveDampenedWLSj(S,B,solution,j)
    solutionAverage<-rowMeans(cbind(newsolution,matrix(solution,nrow = length(solution),ncol = 4)))
    change<-norm(as.matrix(solutionAverage-solution))
    solution<-solutionAverage
    iterations<-iterations+1
    changes<-c(changes,change)
  }
  print(round(solution/sum(solution),5))
  return(solution/sum(solution))
}
