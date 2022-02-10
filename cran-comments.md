## R CMD check results
0 errors ✓ | 0 warnings ✓ | 2 notes x
There were 0 ERRORs and 0 WARNINGS

Comments on the notes: 

Note 1: 
> checking installed package size ... NOTE
    installed size is 80.2Mb
    sub-directories of 1Mb or more:
      data  80.0Mb

There was one note about the size of the files in the /data. I know there 
are (at times) exceptions that are made. I am hopeful that one can be made in 
this case as the data is single-cell gene expression data along with labels/
bulk data in order for the user to see a best real-world example. 

---- 

Note 2: 
> checking R code for possible problems ... NOTE
  DEAnalysisMAST: no visible binding for global variable
    ‘Number.of.Cells’
  Undefined global functions or variables:
    Number.of.Cells
    
There was one note about there not being visible binding for global variables 
'Number.of.cells'. This come from this line of code in the DEAnalysisMAST.R
file.

''' 
colnames(data_for_MIST) = c("Gene", "Subject.ID", "Et", "Population", "Number.of.Cells")
      vbeta = data_for_MIST
      vbeta.fa <-FromFlatDF(vbeta, idvars=c("Subject.ID"),
                  primerid='Gene', measurement='Et', ncells='Number.of.Cells',
                  geneid="Gene",  cellvars=c('Number.of.Cells', 'Population'),
                  phenovars=c('Population'), id='vbeta all')
      vbeta.1 <- subset(vbeta.fa,Number.of.Cells==1)
      
'''
Number.of.cells is not a global variable (rather it comes from the FromFlatDF -  
https://rdrr.io/bioc/MAST/man/FromFlatDF.html within the MAST function) 
therefore I believe this note to be spurious.


The notice ' * checking installed package size ... NOTE ',
CRAN has a limit of 100 MB for a *.tar.gz file (Package is under size requirements)
