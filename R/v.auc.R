#' @title Signature Matrix Using MAST
#'
#' @description This function builds a signature matrix using genes identified
#'   by the DEAnalysisMAST() function.
#'
#' @param symbol
#'
#' @return Signature Matrix
#'
#' @examples
#'
#' @export
#'
#' @importFrom dplyr "%>%"


v.auc = function(data.v,group.v) {
  prediction.use=prediction(data.v, group.v, 0:1)
  perf.use=performance(prediction.use,"auc")
  auc.use=round(perf.use@y.values[[1]],3)
  return(auc.use)
}
