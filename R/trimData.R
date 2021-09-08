#' @title trimData
#'
#' @description This function trims bulk and single-cell data to contain the same genes
#'
#' @param Signature_Matrix  A single-cell signature matrix
#' @param bulkdata A bulk dataset
#'
#' @return A list of trimmed bulk and single-cell data.
#'
#' @examples
#' trimData(Signature_Matrix = Sig, bulkdata = dataBulk)
#'
#' @export
#'
#' @importFrom dplyr "%>%"

trimData<-function(Signature_Matrix,bulkdata){
  Genes<-intersect(rownames(Signature),names(bulkData))
  B<-bulkData[Genes]
  S<-Signature[Genes,]
  return(list("sig"=S,"bulk"=B))
}


