#' @title Signature Matrix Using MAST
#'
#' @description This function builds a signature matrix using genes identified
#'   by the DEAnalysisMAST() function.
#'
#' @param symbol
#'
#' @return Signature Matrix
#'
#' @examples
#'
#' @export buildSignatureMatrixMAST


buildSignatureMatrixMAST<-function(scdata,id,path,diff.cutoff=0.5,pval.cutoff=0.01){
  #compute differentially expressed genes for each cell type
  DEAnalysisMAST(scdata,id,path)

  #for each cell type, choose genes in which FDR adjusted p-value is less than 0.01 and the estimated fold-change
  #is greater than 0.5
  numberofGenes<-c()
  for (i in unique(id)){
    if(file.exists(paste(path,"/",i,"_MIST.RData", sep=""))){
      load(file=paste(path,"/",i,"_MIST.RData", sep=""))
      pvalue_adjusted<-p.adjust(cluster_lrTest.table[,3], method = "fdr", n = length(cluster_lrTest.table[,3]))
      cluster_lrTest.table<-cbind(cluster_lrTest.table,pvalue_adjusted)
      DEGenes<-cluster_lrTest.table$Gene[intersect(which(pvalue_adjusted<pval.cutoff),which(cluster_lrTest.table$log2fold_change>diff.cutoff))]
      nonMir = grep("MIR|Mir", DEGenes, invert = T)  # because Mir gene is usually not accurate
      assign(paste("cluster_lrTest.table.",i,sep=""),cluster_lrTest.table[which(cluster_lrTest.table$Gene%in%DEGenes[nonMir]),])
      numberofGenes<-c(numberofGenes,length(DEGenes[nonMir]))
    }
  }

  #need to reduce number of genes
  #for each subset, order significant genes by decreasing fold change, choose between 50 and 200 genes
  #for each, iterate and choose matrix with lowest condition number
  conditionNumbers<-c()
  for(G in 50:200){
    Genes<-c()
    j=1
    for (i in unique(id)){
      if(numberofGenes[j]>0){
        temp<-paste("cluster_lrTest.table.",i,sep="")
        temp<-as.name(temp)
        temp<-eval(parse(text = temp))
        temp<-temp[order(temp$log2fold_change,decreasing=TRUE),]
        Genes<-c(Genes,varhandle::unfactor(temp$Gene[1:min(G,numberofGenes[j])]))
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
  Genes<-c()
  j=1
  for (i in unique(id)){
    if(numberofGenes[j]>0){
      temp<-paste("cluster_lrTest.table.",i,sep="")
      temp<-as.name(temp)
      temp<-eval(parse(text = temp))
      temp<-temp[order(temp$log2fold_change,decreasing=TRUE),]
      Genes<-c(Genes,varhandle::unfactor(temp$Gene[1:min(G,numberofGenes[j])]))
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
  #save(Sig,file=paste(path,"/Sig.RData",sep=""))
  saveRDS(Sig,file=paste(path,"/Sig.rds",sep=""))
  return(Sig)
}
