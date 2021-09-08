Mean.in.log2space=function(x,pseudo.count) {
  return(log2(mean(2^(x)-pseudo.count)+pseudo.count))
}
