# DWLS: Gene Expression Deconvolution Using Dampened Weighted Least Squares

Dampened weighted least squares (DWLS) is an estimation method for gene 
expression deconvolution, in which the cell-type composition of a bulk RNA-seq 
data set is computationally inferred. This method corrects common biases 
towards cell types that are characterized by highly expressed genes and/or 
are highly prevalent, to provide accurate detection across diverse cell types. 
To begin, the user must input a bulk RNA-seq data set, along with a labeled 
representative single-cell RNA-seq data set that will serve to generate 
cell-type-specific gene expression profiles. Ideally, the single-cell data 
set will contain cells from all cell types that may be found in the bulk data. 
DWLS will return the cell-type composition of the bulk data.

# Data Sources

[1] Schelker M, Feau S, Du J, Ranu N, Klipp E, MacBeath G, Schoeberl B, Raue 
A: Estimation of immune cell content in tumour tissue using single-cell 
RNA-seq data. Nat Commun 2017, 8:2032. 

[2] Yan KS, Janda CY, Chang J, Zheng GXY, Larkin KA, Luca VC, Chia LA, 
Mah AT, Han A, Terry JM, et al: Non-equivalence of Wnt and R-spondin ligands 
during Lgr5. Nature 2017, 545:238-242. 

[3] Han X, Wang R, Zhou Y, Fei L, Sun H, Lai S, Saadatpour A, Zhou Z, Chen H, 
Ye F, et al: Mapping the Mouse Cell Atlas by Microwell-Seq. Cell 2018, 
172:1091-1107.e1017.

# References

[Tsoucas D, Dong R, Chen H, Zhu Q, Guo G, Yuan GC. Accurate estimation of 
cell-type composition from gene expression data. Nat Commun. 
2019 Jul 5;10(1):2975.](https://www.nature.com/articles/s41467-019-10802-z)


-- 
<!-- badges: start -->
[![R-CMD-check](https://github.com/sistia01/DWLS/workflows/R-CMD-check/badge.svg)](https://github.com/sistia01/DWLS/actions)
<!-- badges: end -->
