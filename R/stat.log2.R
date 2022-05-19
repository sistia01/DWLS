#' @title stat.log2
#'
#' @description One of the functions required for the  differential expression
#' analysis using MAST (DEAnalysisMast()) function.
#'
#' @param data.m Data
#' @param group.v Groupings
#' @param pseudo.count A pseudocount value
#'
#' @return A dataframe of the log2 applied results
#'
#' @examples
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
#'   cells.coord.list1 = match(cells.symbol.list1, colnames(data.used.log2))
#'   data.used.log2.ordered = cbind(data.used.log2[,cells.coord.list1],
#'                                         data.used.log2[,cells.coord.list2])}
#' group.v <- c(rep(0,length(cells.coord.list1)),
#'                                          rep(1, length(cells.coord.list2)))
#'}
#' @export stat.log2
#'
#' @importFrom dplyr "%>%"
#' @importFrom stats "aggregate"

stat.log2=function(data.m, group.v, pseudo.count)
  #data.m=data.used.log2
  { log2.mean.r <- aggregate(t(data.m), list(as.character(group.v)),
                           function(x) Mean.in.log2space(x,pseudo.count))
  log2.mean.r <- t(log2.mean.r)
  colnames(log2.mean.r) <- paste("mean.group",log2.mean.r[1,], sep="")
  log2.mean.r = log2.mean.r[-1,]
  log2.mean.r = as.data.frame(log2.mean.r)
  log2.mean.r = varhandle::unfactor(log2.mean.r)  #from varhandle
  log2.mean.r[,1] = as.numeric(log2.mean.r[,1])
  log2.mean.r[,2] = as.numeric(log2.mean.r[,2])
  log2_foldchange = log2.mean.r$mean.group1-log2.mean.r$mean.group0
  results = data.frame(cbind(log2.mean.r$mean.group0,log2.mean.r$mean.group1,
                             log2_foldchange))
  colnames(results) = c("log2.mean.group0","log2.mean.group1","log2_fc")
  rownames(results) = rownames(log2.mean.r)
  return(results)
}


