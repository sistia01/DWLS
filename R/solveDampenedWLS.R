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
#' @param S List output from trimData
#' @param B List output from trimData
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
#' solveDampenedWLS(S, B)
#'}
#' @export solveDampenedWLS
#'
#' @importFrom dplyr "%>%"
#' @importFrom quadprog "solve.QP"


solveDampenedWLS<-function(S,B){
  solution<-solveOLSInternal(S,B)
  #solution <- solveSVR(S,B)
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
