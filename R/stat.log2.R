stat.log2=function(data.m, group.v, pseudo.count){
  #data.m=data.used.log2
  log2.mean.r <- aggregate(t(data.m), list(as.character(group.v)), function(x) Mean.in.log2space(x,pseudo.count))
  log2.mean.r <- t(log2.mean.r)
  colnames(log2.mean.r) <- paste("mean.group",log2.mean.r[1,], sep="")
  log2.mean.r = log2.mean.r[-1,]
  log2.mean.r = as.data.frame(log2.mean.r)
  log2.mean.r = varhandle::unfactor(log2.mean.r)  #from varhandle
  log2.mean.r[,1] = as.numeric(log2.mean.r[,1])
  log2.mean.r[,2] = as.numeric(log2.mean.r[,2])
  log2_foldchange = log2.mean.r$mean.group1-log2.mean.r$mean.group0
  results = data.frame(cbind(log2.mean.r$mean.group0,log2.mean.r$mean.group1,log2_foldchange))
  colnames(results) = c("log2.mean.group0","log2.mean.group1","log2_fc")
  rownames(results) = rownames(log2.mean.r)
  return(results)
}
