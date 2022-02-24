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
#' trimData(Signature_Matrix = Sig, bulkdata = dataBulk)
#'
#' @export trimData
#'
#' @importFrom dplyr "%>%"
trimData<-function(Signature_Matrix,bulkdata){
  Genes<-intersect(rownames(Signature_Matrix),names(bulkdata))
  B<-bulkdata[Genes]
  S<-Signature_Matrix[Genes,]
  return(list("sig"=S,"bulk"=B))
}


