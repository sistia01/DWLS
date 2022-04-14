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
#' @return An RData and rds file with the differential expression analysis
#' results for each unique id.
#'
#' @examples
#' \donttest{
#' download.file(
#' "https://github.com/sistia01/DWLS/raw/main/inst/extdata/dataSC.RData",
#' "dataSC.RData")
#' load("dataSC.RData")
#'
#' data('dataBulk', package = "DWLS")
#' download.file(
#' "https://github.com/sistia01/DWLS/raw/main/inst/extdata/dataBulk.RData",
#' "dataBulk.RData")
#' load("dataBulk.RData")
#'
#' #data('labels', package = "DWLS")
#' download.file(
#' "https://github.com/sistia01/DWLS/raw/main/inst/extdata/labels.RData",
#' "labels.RData")
#' load("labels.RData")
#'
#' #data('trueLabels', package = "DWLS")
#' download.file(
#' "https://github.com/sistia01/DWLS/raw/main/inst/extdata/trueLabels.RData",
#' "trueLabels.RData")
#' load("trueLabels.RData")
#'
#' #Old Method
#' #load("data/dataBulk.RData") #read in bulk data for WT1 (control condition #1)
#' #load("data/labels.RData") #read in single-cell labels from clustering
#'
#' labels<-trueLabels
# #Change to real labels
#' newcat<-c("NonCycISC","CycISC","TA","Ent","PreEnt","Goblet","Paneth","Tuft","EE")
#' for (i in 1:length(newcat)){
#'   labels[which(labels==(i-1))]<-newcat[i]
#'   }
#' #Run deconvolution
#' Seurat_DE <- DEAnalysisSeurat(dataSC, labels, "inst/extdata/results")
#'}
#' @export DEAnalysisSeurat
#'
#' @importFrom dplyr "%>%"
#' @importFrom Seurat "FindMarkers"
#'
DEAnalysisSeurat<-function(scdata,id,path)
  { exprObj<-CreateSeuratObject(counts=as.data.frame(scdata), project = "DE")
  #print("Calculating differentially expressed genes:")
  for (i in unique(id)){
    de_group <- FindMarkers(object=exprObj, ident.1 = i, ident.2 = NULL,
                            only.pos = TRUE,
                            test.use = "bimod", group.by = as.vector(id))
    saveRDS(de_group,file=paste(path,"/de_",i,".rds",sep=""))
    save(de_group,file=paste(path,"/de_",i,".RData",sep=""))
    #print("RData differential expression results are in the'results' folder")
    print(i)
  }
}

