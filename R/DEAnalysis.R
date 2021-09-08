DEAnalysis<-function(scdata,id,path){
  exprObj<-CreateSeuratObject(counts=as.data.frame(scdata), project = "DE")
  print("Calculating differentially expressed genes:")
  for (i in unique(id)){
    de_group <- FindMarkers(object=exprObj, ident.1 = i, ident.2 = NULL,
                            only.pos = TRUE, test.use = "bimod", group.by = as.vector(id))
    saveRDS(de_group,file=paste(path,"/de_",i,".rds",sep=""))
    save(de_group,file=paste(path,"/de_",i,".RData",sep=""))
  }
}
