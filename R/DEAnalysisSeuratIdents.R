#' @title Differential Expression with Idents (Seurat)
#'
#' @description This function calculates the differential expression values
#' along with identifying the Idents (through Seurat). The output is
#' saved in an RData file for each unique identity (id).
#'
#' @param scdata The gene expression datafarme
#' @param id The unique identities within the data
#' @param path The path for the RData results
#'
#' @return An RData file with the differential expression analysis results
#' for each unique id.
#'
#' @examples
#' \donttest{
#'
#' #dataSC
#' load(system.file("extdata", "dataSC.RData", package = "DWLS"))
#'
#' #dataBulk
#' load(system.file("extdata", "dataBulk.RData", package = "DWLS"))
#'
#' #labels
#' load(system.file("extdata", "labels.RData", package = "DWLS"))
#'
#' #trueLabels
#' load(system.file("extdata", "trueLabels.RData", package = "DWLS"))
#'
#' labels<-trueLabels
#' #Change to real labels
#' newcat<-c("NonCycISC","CycISC","TA","Ent","PreEnt","Goblet","Paneth",
#' "Tuft","EE")
#' for (i in 1:length(newcat)){
#'   labels[which(labels==(i-1))]<-newcat[i]
#'   }
#' #Run deconvolution
#' Seurat_test2 <- DEAnalysisSeuratIdents(dataSC, labels, "results")
#'}
#'
#' @export DEAnalysisSeuratIdents
#'
#' @importFrom Seurat "CreateSeuratObject"
#' @importFrom Seurat "FindMarkers"
#' @importFrom Seurat "Idents"

DEAnalysisSeuratIdents<-function(scdata,id,path)
  { exprObj<- CreateSeuratObject(counts=as.data.frame(scdata), project = "DE")
  Idents(object = exprObj) <- as.vector(id)
  #print("Calculating differentially expressed genes:")
  for (i in unique(id)){
    #de_group <- FindMarkers(object=exprObj, ident.1 = i, ident.2 = NULL,
                            #only.pos = TRUE, test.use = "bimod")
    de_group <- FindMarkers(object=exprObj, ident.1 = i, ident.2 = NULL,
                            only.pos = TRUE, test.use = "MAST")
    saveRDS(de_group,file=paste(path,"/de_",i,".rds",sep=""))
    save(de_group,file=paste(path,"/de_",i,".RData",sep=""))
  }
}
