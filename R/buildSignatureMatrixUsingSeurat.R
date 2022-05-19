#' @title Signature Matrix Using Seurat
#'
#' @description This function builds a signature matrix using genes identified
#'   by the DEAnalysis() function.
#'
#' @param scdata The data
#' @param id The identities of the genes
#' @param path The path to the file results
#' @param diff.cutoff This is automatically set to 0.5
#' @param pval.cutoff The p-value cutoff. This is automatically set to 0.01
#' @param f The maximum number of genes (when creating the signature matrix,
#' need to reduce number of genes, between 50:f number of significant genes are
#' chosen). If not set, this number is automatically set to 200.
#'
#' @return Signature Matrix built using the Seurat algorithm
#'
#' @examples
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
#' #trueLabels
#' #data('trueLabels', package = "DWLS")
#' #url <- "https://github.com/sistia01/DWLS/raw/main/inst/extdata/trueLabels.RData"
#' #dest <- "data/trueLabels.RData"
#' #download.file(url, dest)
#' #load("data/trueLabels.RData")
#' load(system.file("extdata", "trueLabels.RData", package = "DWLS"))
#'
#' #Old Method
#' #load("data/dataSC_3.RData")
#' #load("data/trueLabels.RData")
#' #load("data/dataBulk.RData") #read in bulk data for WT1 (control condition#1)
#' #load("data/labels.RData") #read in single-cell labels from clustering
#'
#' labels<-trueLabels
#'
#' #Change to real labels
#' newcat<-c("NonCycISC","CycISC","TA","Ent","PreEnt","Goblet","Paneth","Tuft",
#' "EE")
#' for (i in 1:length(newcat)){
#'   labels[which(labels==(i-1))]<-newcat[i]
#'   }
#'
#' #Run on local with inst/extdata/results folder
#' #Signature<-buildSignatureMatrixUsingSeurat(dataSC,labels,
#' # "inst/extdata/results", diff.cutoff=0.5,pval.cutoff=0.01)
#'}
#'
#' @export buildSignatureMatrixUsingSeurat

buildSignatureMatrixUsingSeurat<-function(scdata,id,path,diff.cutoff=0.5,
                                          pval.cutoff=0.01, f=200)
  { DEAnalysisSeuratIdents(scdata,id,path)

  numberofGenes<-c()
  for (i in unique(id)){
    #load(file=paste(path,"/de_",i,".RData",sep=""))
    de_group <- readRDS(file=paste(path,"/de_",i,".rds", sep =""))
    #load(file=paste(path,"/de_",i,".rds", sep =""))
    DEGenes<-rownames(de_group)[intersect(which(de_group$p_val_adj<pval.cutoff),
                                        which(de_group$avg_log2FC>diff.cutoff))]
    nonMir = grep("MIR|Mir", DEGenes, invert = T)
    assign(paste("cluster_lrTest.table.",i,sep=""),
           de_group[which(rownames(de_group)%in%DEGenes[nonMir]),])
    numberofGenes<-c(numberofGenes,length(DEGenes[nonMir]))
  }

  #need to reduce number of genes
  #for each subset, order significant genes by decreasing fold change,
  #choose between 50 and 200 genes
  #choose matrix with lowest condition number
  conditionNumbers<-c()
  for(G in 50:f){
    #changed it from for(G in 50:100) -- should find a way to do "best"
    Genes<-c()
    j=1
    for (i in unique(id)){
      if(numberofGenes[j]>0){
        temp<-paste("cluster_lrTest.table.",i,sep="")
        temp<-as.name(temp)
        temp<-eval(parse(text = temp))
        temp<-temp[order(temp$p_val_adj,decreasing=TRUE),]
        Genes<-c(Genes,(rownames(temp)[1:min(G,numberofGenes[j])]))
      }
      j=j+1
    }
    Genes<-unique(Genes)
    #make signature matrix
    ExprSubset<-scdata[Genes,]
    Sig<-NULL
    for (i in unique(id)){
      Sig<-cbind(Sig,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
    }
    colnames(Sig)<-unique(id)
    conditionNumbers<-c(conditionNumbers,kappa(Sig))
  }
  G<-which.min(conditionNumbers)+min(49,numberofGenes-1)
  #G is optimal gene number
  Genes<-c()
  j=1
  for (i in unique(id)){
    if(numberofGenes[j]>0){
      temp<-paste("cluster_lrTest.table.",i,sep="")
      temp<-as.name(temp)
      temp<-eval(parse(text = temp))
      temp<-temp[order(temp$p_val_adj,decreasing=TRUE),]
      Genes<-c(Genes,(rownames(temp)[1:min(G,numberofGenes[j])]))
    }
    j=j+1
  }
  Genes<-unique(Genes)
  ExprSubset<-scdata[Genes,]
  Sig<-NULL
  for (i in unique(id)){
    Sig<-cbind(Sig,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
  }
  colnames(Sig)<-unique(id)
  save(Sig,file=paste(path,"/Sig.RData",sep=""))
  return(Sig)
}
