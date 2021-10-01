#' @title Differential Expression without Idents (Seurat)
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
#' \dontrun{
#' load("data/dataSC_1.RData")
#' load("data/dataSC_2.RData")
#' dataSC <- cbind(dataSC_1, dataSC_2)
#' load("data/trueLabels.RData")
#' load("data/dataBulk.RData") #read in bulk data for WT1 (control condition #1)
#' load("data/labels.RData") #read in single-cell labels from clustering
#' labels<-trueLabels
# #Change to real labels
#' newcat<-c("NonCycISC","CycISC","TA","Ent","PreEnt","Goblet","Paneth","Tuft","EE")
#' for (i in 1:length(newcat)){
#'   labels[which(labels==(i-1))]<-newcat[i]
#'   }
#' #Run deconvolution
#' Seurat_test2 <- DEAnalysisSeurat(dataSC, labels, "results")
#' }
#'
#' @export DEAnalysisSeurat
#'
#' @importFrom dplyr "%>%"
#' @importFrom Seurat "FindMarkers"
#'
DEAnalysisSeurat<-function(scdata,id,path){
  exprObj<-CreateSeuratObject(counts=as.data.frame(scdata), project = "DE")
  print("Calculating differentially expressed genes:")
  for (i in unique(id)){
    de_group <- FindMarkers(object=exprObj, ident.1 = i, ident.2 = NULL,
                            only.pos = TRUE, test.use = "bimod", group.by = as.vector(id))
    saveRDS(de_group,file=paste(path,"/de_",i,".rds",sep=""))
    save(de_group,file=paste(path,"/de_",i,".RData",sep=""))
    print("The RData differential expression results are in the 'results' folder")
  }
}
