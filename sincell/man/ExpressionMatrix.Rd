\docType{data}
\name{ExpressionMatrix}
\alias{ExpressionMatrix}
\title{Single-cell expression data for genes differentially expressed in differentiating human skeletal muscle myoblasts cells}
\format{Numeric matrix}
\source{
sincell
}
\usage{
data(ExpressionMatrix)
ExpressionMatrix
}
\description{
This dataset contains expression profiles from a time-series study of differentiating human skeletal muscle myoblasts (object HSMM in Bioconductor package monocle). Expression values are in FPKM units. Data is part of a publicly available single-cell RNA-seq dataset from Trapnell et al 2014. In this work, authors generated single-cell RNA-seq libraries for differentiating myoblasts at 0, 24, 48 and 72 hours. Original data can be accessed at GEO database accession number GSE52529. Following Trapnell et al 2014 and the vignette of its associated Bioconductor package Monocle, the expression matrix is restricted to genes differentially expressed between cells from times 0 and the ensemble of cells of times 24, 28 and 72 hours of differentiation. Steps to achieve this are reported in monocle's  vignette. Those steps produce the matrix ExpressionMatrix represent¡ng the expression profiles of those differentially expressed genes. ExpressionMatrix is provided as part of Sincell package in order to keep the running time of its vignette short.
}

