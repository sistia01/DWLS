#' @title v.auc
#'
#' @description Uses the prediction() function in order to create standardized
#' output from the data in order to perform an AUC calculation. The calculation
#' results are rounded to the third decimal place. This function serves mainly
#' to support the DWLS function.
#'
#' @param data.v Data
#' @param group.v Data subdivision
#'
#' @return Matrix of standardized output of AUC calculation
#'
#' @examples
#'
#' \dontrun{
#' m.auc=function(data.m,group.v) {
#' AUC=apply(data.m, 1, function(x) v.auc(x,group.v))
#' AUC[is.na(AUC)]=0.5
#' return(AUC)}
#' }
#'
#' @export v.auc
#'
#' @importFrom dplyr "%>%"
#' @importFrom ROCR "prediction"
#' @importFrom ROCR "performance"


v.auc = function(data.v,group.v) {
  prediction.use=prediction(data.v, group.v, 0:1)
  perf.use=performance(prediction.use,"auc")
  auc.use=round(perf.use@y.values[[1]],3)
  return(auc.use)
}
