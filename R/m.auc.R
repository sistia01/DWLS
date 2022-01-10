#' @title m.auc
#' @description Calculates the AUC of a dataset. The function mainly serves
#' to support the DWLS function.
#' @param dataset Data
#' @param grouping Data subdivision
#' @return Matrix of standardized output of AUC calculation
#' @examples
#' \dontrun{
#' pseudo.count = 0.1
#' data.used.log2   <- log2(scdata+pseudo.count)
#' colnames(data.used.log2)<-make.unique(colnames(data.used.log2))
#' diff.cutoff=0.5
#' for (i in unique(id)){
#'   cells.symbol.list2 = colnames(data.used.log2)[which(id==i)]
#'   cells.coord.list2 = match(cells.symbol.list2, colnames(data.used.log2))
#'   cells.symbol.list1 = colnames(data.used.log2)[which(id != i)]
#'   cells.coord.list1= match(cells.symbol.list1, colnames(data.used.log2))
#'   data.used.log2.ordered = cbind(data.used.log2[,cells.coord.list1],
#'                                          data.used.log2[,cells.coord.list2])
#'   group.v <- c(rep(0,length(cells.coord.list1)),
#'                                rep(1, length(cells.coord.list2)))
#'   #ouput
#'   log2.stat.result <- stat.log2(data.used.log2.ordered,
#'                                     group.v, pseudo.count)
#'   Auc <- m.auc(data.used.log2.ordered, group.v)}
#' }
#'
#' @export m.auc
#' @importFrom dplyr "%>%"

m.auc = function(dataset, grouping)
  { AUC=apply(dataset, 1, function(x) v.auc(x,grouping))
  AUC[is.na(AUC)]=0.5
  return(AUC)
}
