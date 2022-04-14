#' @title trimData
#'
#' @description This function trims bulk and single-cell data to
#' contain the same genes. The result is a list of the intersecting genes
#' within the two datasets.
#'
#' @param Signature_Matrix A single-cell signature matrix
#' @param bulkdata A bulk dataset
#'
#' @return A list of trimmed bulk and single-cell data.
#'
#' @examples
#'  \donttest{
#' download.file(
#' "https://github.com/sistia01/DWLS/raw/main/inst/extdata/Sig.RData",
#' "Sig.RData")
#' load("Sig.RData")
#'
#' #data('dataBulk', package = "DWLS")
#' download.file(
#' "https://github.com/sistia01/DWLS/raw/main/inst/extdata/dataBulk.RData",
#' "dataBulk.RData")
#' load("dataBulk.RData")
#'
#' trimData(Signature_Matrix = Sig, bulkdata = dataBulk)
#'}
#' @export trimData
#'
#' @importFrom dplyr "%>%"

trimData<-function(Signature_Matrix,bulkdata){
  Genes<-intersect(rownames(Signature_Matrix),names(bulkdata))
  B<-bulkdata[Genes]
  S<-Signature_Matrix[Genes,]
  return(list("sig"=S,"bulk"=B))
}
