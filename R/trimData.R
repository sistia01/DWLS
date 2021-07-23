#' @title trimData
#'
#' @description This function trims bulk and single-cell data to contain the same genes
#'
#' @param Signature matrix and bulkdata
#'
#' @return NULL
#'
#' @examples trimData(Signature, bulkData)
#'
#' @export

trimData<-function(Signature,bulkData){
  Genes<-intersect(rownames(Signature),names(bulkData))
  B<-bulkData[Genes]
  S<-Signature[Genes,]
  return(list("sig"=S,"bulk"=B))
}
