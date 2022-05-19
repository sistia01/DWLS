#'
#' @title solveSVR
#'
#' @description Performs a support vector regression (SVR). First, the data is
#' scaled then it solves for the SVR. An svm model is used with the following
#' specifications nu=0.5,scale = TRUE, type = "nu-regression",
#' kernel ="linear",cost = 1.
#'
#' Nu-support vector regression was performed using the svm function in the
#' e1071 package in R. Parameters were set to nu = 0.5, type = “nu-regression”,
#' kernel = “linear”, cost = 1, and all others to the default values.
#' Bulk data and signature matrices were scaled to -1, 1. These parameter
#' and scaling choices match those specified in Schelker et al. in their
#' MATLAB code, accessed through https://figshare.com/s/865e694ad06d5857db4b.
#' As in Newman et al., model coefficients are extracted from the svm model
#' using t(model$coefs) model$SV, and any negative coefficients are set
#' to zero. The coefficients are then scaled by the sum of the coefficients,
#' such that the scaled coefficients will sum to one.
#'
#' Citations:
#' Newman, A. M. et al. Robust enumeration of cell subsets from tissue
#' expression profiles. Nat. Methods 12, 453–457 (2015).
#'
#' Schelker, M. et al. Estimation of immune cell content in tumor
#' tissue using single-cell RNA-seq data. Nat. Commun. 8, 2032 (2017).
#'
#'
#' @param S List output from trimData$sig (S)
#' @param B List output from trimData$bulk (B)
#'
#' @return Value (SVR)
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
#' solveSVR(S, B)
#'}
#' @export solveSVR
#'
#' @importFrom dplyr "%>%"
#' @importFrom e1071 "svm"

solveSVR<-function(S,B){
  ub=max(c(as.vector(S),B)) #upper bound
  lb=min(c(as.vector(S),B)) #lower bound
  Bs=((B-lb)/ub)*2-1
  Ss=((S-lb)/ub)*2-1
  model<-svm(Ss,Bs, nu=0.5,scale = TRUE, type = "nu-regression",
             kernel ="linear",cost = 1)
  coef <- t(model$coefs) %*% model$SV
  coef[which(coef<0)]<-0
  coef<-as.vector(coef)
  names(coef)<-colnames(S)
  print(round(coef/sum(coef),5))
  return(coef/sum(coef))
  }
