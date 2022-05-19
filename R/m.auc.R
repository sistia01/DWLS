#' @title m.auc
#' @description Calculates the AUC of a dataset. The function mainly serves
#' to support the DWLS function.
#' @param dataset Data
#' @param grouping Data subdivision
#' @return Matrix of standardized output of AUC calculation
#' @examples
#'
#' \donttest{
#'
#' #dataSC
#' #url <- "https://github.com/sistia01/DWLS/raw/main/inst/extdata/dataSC.RData"
#' #dest <- "data/dataSC.RData"
#' #load(download.file(url, tempfile(data/dataSC.RData))
#' #load("dataSC.RData")
#' #SOLUTION
#' load(system.file("extdata", "dataSC.RData", package = "DWLS"))
#'
#' #dataBulk
#' #url <- "https://github.com/sistia01/DWLS/raw/main/inst/extdata/dataBulk.RData"
#' #dest <- "data/dataBulk.RData"
#' #load(download.file(url, tempfile(dest)))
#' #load("data/dataBulk.RData")
#' load(system.file("extdata", "dataBulk.RData", package = "DWLS"))
#'
#' #labels
#' #url <- "https://github.com/sistia01/DWLS/raw/main/inst/extdata/labels.RData"
#' #dest <- "data/labels.RData"
#' #download.file(url, dest)
#' #load("data/labels.RData")
#' load(system.file("extdata", "labels.RData", package = "DWLS"))
#'
#' #data('trueLabels', package = "DWLS")
#' #url <- "https://github.com/sistia01/DWLS/raw/main/inst/extdata/trueLabels.RData"
#' #dest <- "data/trueLabels.RData"
#' #download.file(url, dest)
#' #load("data/trueLabels.RData")
#' load(system.file("extdata", "trueLabels.RData", package = "DWLS"))
#'
#' pseudo.count = 0.1
#' data.used.log2 <- log2(dataSC+pseudo.count)
#' colnames(data.used.log2)<-make.unique(colnames(data.used.log2))
#' diff.cutoff=0.5
#' id = labels
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
#'}
#'
#' @export m.auc
#'
#' @importFrom dplyr "%>%"
#' @importFrom ROCR "prediction"
#' @importFrom ROCR "performance"

m.auc = function(dataset, grouping)
  { AUC=apply(dataset, 1, function(x) v.auc(x,grouping))
  AUC[is.na(AUC)]=0.5
  return(AUC)
}
