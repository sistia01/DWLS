#' @title DEAnalysisMAST
#'
#' @description Perform DE analysis using MAST.
#' Dampened weighted least squares (DLWS) is an estimation
#' method for gene expression deconvolution, in which the cell-type
#' composition of a bulk RNA-seq data set is computationally inferred.
#' This method corrects common biases towards cell types that are
#' characterized by highly expressed genes and/or are highly prevalent,
#' to provide accurate detection across diverse cell types. To begin,
#' the user must input a bulk RNA-seq data set, along with a labeled
#' representative single-cell RNA-seq data set that will serve to generate
#' cell-type-specific gene expression profiles. Ideally, the single-cell
#' data set will contain cells from all cell types that may be found in the
#' bulk data. DWLS will return the cell-type composition of the bulk data.
#' First, solve OLS then use the solution to find a starting point for the
#' weights. Next, the dampened weighted least squares is performed. The weights
#' are iterated until convergence then the dampening constant for weights is
#' found using cross-validation (with decreasing step size for convergence).
#'
#' DWLS captures ISC composition changes across conditions.
#' One of the most important applications of deconvolution methods is in
#' the identification of cell-type composition variations across conditions.
#'
#' Note: The function uses solveDampenedWLSj() and findDampeningConstant().
#'
#' @param scdata The gene expression datafarme
#' @param id The unique identities within the data
#' @param path The path for the RData results
#'
#' @return matrix
#' The resulting matrix is a gene by cell-type signature matrix.
#' The cell-type signature matrix is constructed using a representative
#' single-cell data set, such that all cell types expected in the bulk
#' data are also represented in the single-cell data (the converse need not be
#' true). The single-cell data is first clustered to reveal
#' its constituent cell types.The function will return 3 different files:
#' an RData file, an rds file, and a csv file.
#'
#' @examples
#' \donttest{
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
#' url <- "https://github.com/sistia01/DWLS/raw/main/inst/extdata/dataBulk.RData"
#' dest <- "data/dataBulk.RData"
#' load(download.file(url, tempfile(dest)))
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
#' #data('trueLabels', package = "DWLS")
#' #url <- "https://github.com/sistia01/DWLS/raw/main/inst/extdata/trueLabels.RData"
#' #dest <- "data/trueLabels.RData"
#' #download.file(url, dest)
#' #load("data/trueLabels.RData")
#' load(system.file("extdata", "trueLabels.RData", package = "DWLS"))
#'
#' labels<-trueLabels
#'
#' #Old Method
#' #load("data/dataBulk.RData") #read in bulk data for WT1 (control condition #1)
#' #load("data/labels.RData") #read in single-cell labels from clustering
#' #data('dataSC_3', package = "DWLS")
#' #dataSC <- dataSC_3
#'
#' labels<-trueLabels
#'
#' #Old Method
#' #load("data/dataBulk.RData") #read in bulk data for WT1 (control condition #1)
#' #load("data/labels.RData") #read in single-cell labels from clustering
#' #data('dataSC_3', package = "DWLS")
#' #dataSC <- dataSC_3
#'
#' labels<-trueLabels
#'
#' #Change to real labels
#' newcat<-c("NonCycISC","CycISC","TA","Ent","PreEnt","Goblet","Paneth","Tuft",
#'  "EE")
#' for (i in 1:length(newcat)){
#'   labels[which(labels==(i-1))]<-newcat[i]
#'   }
#'
#' #Run deconvolution
#' #Results are in inst/extdata/results folder -- run on local
#' #Example code below
#' Mast_test <- DEAnalysisMAST(dataSC, labels, "inst/extdata/results")
#'}
#'
#' @export DEAnalysisMAST
#'
#' @importFrom dplyr "%>%"
#' @importFrom MAST "FromFlatDF"
#' @importFrom MAST "zlm"
#' @importFrom MAST "show"
#' @importFrom MAST "lrTest"
#' @importFrom SummarizedExperiment "colData"
#' @importFrom reshape "melt"
#' @importFrom utils "head"
#' @importFrom utils "write.csv"

DEAnalysisMAST<-function(scdata,id,path)
  { Number.of.Cells <- NULL
  pseudo.count = 0.1
  data.used.log2  <- log2(scdata+pseudo.count)
  colnames(data.used.log2)<-make.unique(colnames(data.used.log2))
  diff.cutoff=0.5
  for (i in unique(id)){
    cells.symbol.list2     = colnames(data.used.log2)[which(id==i)]
    cells.coord.list2      = match(cells.symbol.list2, colnames(data.used.log2))
    cells.symbol.list1     = colnames(data.used.log2)[which(id != i)]
    cells.coord.list1      = match(cells.symbol.list1, colnames(data.used.log2))
    data.used.log2.ordered  = cbind(data.used.log2[,cells.coord.list1],
                                    data.used.log2[,cells.coord.list2])
    group.v <- c(rep(0,length(cells.coord.list1)),
                 rep(1, length(cells.coord.list2)))
    #ouput
    log2.stat.result <- stat.log2(data.used.log2.ordered, group.v, pseudo.count)
    Auc <- m.auc(data.used.log2.ordered, group.v)
    bigtable <- data.frame(cbind(log2.stat.result, Auc))

    DE <- bigtable[bigtable$log2_fc >diff.cutoff,]
    dim(DE)
    if(dim(DE)[1]>1){
      data.1 = data.used.log2[,cells.coord.list1]
      data.2 = data.used.log2[,cells.coord.list2]
      genes.list = rownames(DE)
      log2fold_change        = cbind(genes.list, DE$log2_fc)
      colnames(log2fold_change) = c("gene.name", "log2fold_change")
      counts  = as.data.frame(cbind( data.1[genes.list,], data.2[genes.list,] ))
      groups  = c(rep("Cluster_Other", length(cells.coord.list1) ), rep(i, length(cells.coord.list2) ) )
      groups  = as.character(groups)
      data_for_MIST <- as.data.frame(cbind(rep(rownames(counts), dim(counts)[2]), melt(counts),rep(groups, each = dim(counts)[1]), rep(1, dim(counts)[1] * dim(counts)[2]) ))
      colnames(data_for_MIST) = c("Gene", "Subject.ID", "Et", "Population", "Number.of.Cells")
      vbeta = data_for_MIST
      #vbeta.1 <- vbeta.fa[vbeta.fa[["Number.of.Cells"]] == 1, , drop = FALSE]
      vbeta.fa <-FromFlatDF(vbeta, idvars=c("Subject.ID"),
                            primerid='Gene', measurement='Et', ncells='Number.of.Cells',
                            geneid="Gene",  cellvars=c('Number.of.Cells', 'Population'),
                            phenovars=c('Population'), id='vbeta all')
      vbeta.1 <- subset(vbeta.fa,Number.of.Cells==1)
      # .3 MAST
      head(colData(vbeta.1))
      zlm.output <- zlm(~ Population, vbeta.1, method='bayesglm', ebayes=TRUE, parallel = TRUE)
      show(zlm.output)
      coefAndCI <- summary(zlm.output, logFC=TRUE)
      zlm.lr <- lrTest(zlm.output, 'Population')
      zlm.lr_pvalue <- melt(zlm.lr[,,'Pr(>Chisq)'])
      zlm.lr_pvalue <- zlm.lr_pvalue[which(zlm.lr_pvalue$test.type == 'hurdle'),]

      lrTest.table <- merge(zlm.lr_pvalue, DE, by.x = "primerid", by.y = "row.names")
      colnames(lrTest.table) <- c("Gene", "test.type", "p_value", paste("log2.mean.", "Cluster_Other", sep=""), paste("log2.mean.",i,sep=""), "log2fold_change", "Auc")
      cluster_lrTest.table <- lrTest.table[rev(order(lrTest.table$Auc)),]
      #. 4 save results
      write.csv(cluster_lrTest.table, file=paste(path,"/",i,"_lrTest.csv", sep=""))
      save(cluster_lrTest.table, file=paste(path,"/",i,"_MIST.RData", sep=""))
      saveRDS(cluster_lrTest.table, file=paste(path,"/",i,"_MIST.rds", sep=""))
      #print("The RData differential expression results are in the 'results' folder for:")
      print(i)
    }
  }
}
