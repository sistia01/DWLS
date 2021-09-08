#' @title Differential Expression with Idents (Seurat)
#'
#' @description This function calculates the differential expression values
#' along with identifying the Idents (through Seurat). The output is
#' saved in an RDS file.
#'
#' @param symbol
#'
#' @return Differential Expression Values
#'
#' @examples
#'
#' @export DEAnalysis_Ident
#'
#'
DEAnalysis_Ident<-function(scdata,id,path){
  exprObj<-CreateSeuratObject(counts=as.data.frame(scdata), project = "DE")
  Idents(object = exprObj) <- as.vector(labels)
  print("Calculating differentially expressed genes:")
  for (i in unique(id)){
    de_group <- FindMarkers(object=exprObj, ident.1 = i, ident.2 = NULL,
                            only.pos = TRUE, test.use = "bimod")
    saveRDS(de_group,file=paste(path,"/de_",i,".rds",sep=""))
    save(de_group,file=paste(path,"/de_",i,".RData",sep=""))
  }
}
