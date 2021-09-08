#' @title Signature Matrix Using Seurat
#'
#' @description This function builds a signature matrix using genes identified
#'   by the DEAnalysis() function.
#'
#' @param symbol
#'
#' @return Signature Matrix
#'
#' @examples
#'
#' @export buildSignatureMatrixUsingSeurat

buildSignatureMatrixUsingSeurat<-function(scdata,id,path,diff.cutoff=0.5,pval.cutoff=0.01){

  #perform differential expression analysis
  DEAnalysis(scdata,id,path)

  numberofGenes<-c()
  for (i in unique(id)){
    #load(file=paste(path,"/de_",i,".RData",sep=""))
    readRDS(file=paste(path,"/de_",i,".rds", sep =""))
    #load(file=paste("results","/de_",i,".rds", sep =""))
    DEGenes<-rownames(de_group)[intersect(which(de_group$p_val_adj<pval.cutoff),which(de_group$avg_log2FC>diff.cutoff))]
    nonMir = grep("MIR|Mir", DEGenes, invert = T)
    assign(paste("cluster_lrTest.table.",i,sep=""),de_group[which(rownames(de_group)%in%DEGenes[nonMir]),])
    numberofGenes<-c(numberofGenes,length(DEGenes[nonMir]))
  }

  #need to reduce number of genes
  #for each subset, order significant genes by decreasing fold change, choose between 50 and 200 genes
  #choose matrix with lowest condition number
  conditionNumbers<-c()
  for(G in 50:200){
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
  G<-which.min(conditionNumbers)+min(49,numberofGenes-1) #G is optimal gene number
  #
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
