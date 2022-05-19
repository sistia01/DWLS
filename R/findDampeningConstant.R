#' @title findDampeningConstant
#'
#' @description Finds a dampening constant for the weights using cross-validation.
#' The goldStandard is used to define the weights. Multiple values of the
#' dampening constant (multiplier) are tried. For each attempt,
#' the variance of the dampened weighted solution for a subset of genes is
#' calculated (on a randomly selected half of the genes).
#' Note that infinite weights are ignored.The dampening constant that
#' results in least cross-validation variance is chosen. It functions
#' in a nondeterministic manner. The dampening constant defines the maximum
#' value that any weight can take on.
#'
#' @param S List output from trimData$Sig (S)
#' @param B List output from trimData$dataBulk (B)
#' @param goldStandard Starting point for the weights, determined by solving OLS
#'
#' @return value (dampening constant value)
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
#' findDampeningConstant(S, B, solution)
#'}
#'
#' @export findDampeningConstant
#'
#' @importFrom dplyr "%>%"
#' @importFrom stats "lm"
#' @importFrom stats "sd"

findDampeningConstant<-function(S,B,goldStandard){
  solutionsSd<-NULL
  sol<-goldStandard
  ws<-as.vector((1/(S%*%sol))^2)
  wsScaled<-ws/min(ws)
  wsScaledMinusInf<-wsScaled
  if(max(wsScaled)=="Inf"){
    wsScaledMinusInf<-wsScaled[-which(wsScaled=="Inf")]
  }
  for (j in 1:ceiling(log2(max(wsScaledMinusInf)))){
    multiplier<-1*2^(j-1)
    wsDampened<-wsScaled
    wsDampened[which(wsScaled>multiplier)]<-multiplier
    solutions<-NULL
    seeds<-c(1:100)
    for (i in 1:100){
      set.seed(seeds[i]) #make nondeterministic
      subset<-sample(length(ws),size=length(ws)*0.5)
      fit = lm (B[subset] ~ -1+S[subset,],weights=wsDampened[subset])
      sol<-fit$coef*sum(goldStandard)/sum(fit$coef)
      solutions<-cbind(solutions,sol)
    }
    solutionsSd<-cbind(solutionsSd,apply(solutions,1,sd))
  }
  j<-which.min(colMeans(solutionsSd^2))
  return(j)
}
