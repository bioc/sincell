\name{sc_StatisticalSupportByGeneSubsampling}
\alias{sc_StatisticalSupportByGeneSubsampling}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Statistical support of cell-state hierarchies by gene subsampling
}
\description{
Function to provide statistical support to the connected graph in SincellObject[["cellstateHierarchy"]] assessed by function sc_GraphBuilderObj() representing a cell-state hierarchy. 
sc_StatisticalSupportByGeneSubsampling() performs "num_it" times a random subsampling of a given number "num_genes" of genes on the original gene expression matrix data in SincellObject[["expressionmatrix"]]. Then, for each resampling, a new connected graph of cells is assessed by calling sc_GraphBuilderObj() with same parameters as for the original SincellObject[["cellstateHierarchy"]]. In each subsampling, the similarity between the resulting connected graph and the original one is assessed as the spearman rank correlation between the two graphs of the shortest distance for all pairs of cells. The distribution of spearman rank correlation values of all iterations is stored as a vector in SincellObject[["StatisticalSupportbyGeneSubsampling"]] and a summary is printed in the standard output.

}
\usage{
sc_StatisticalSupportByGeneSubsampling(SincellObject, num_it=100, 
  num_genes=as.integer(nrow(SincellObject[["expressionmatrix"]])*0.5), 
  cores=ifelse(detectCores()>=4, 4, detectCores()))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SincellObject}{
A SincellObject named list as created by function sc_GraphBuilderObj(), containing in member "cellstateHierarchy" a connected graph representing a cell-state hierarchy.
}
  \item{num_it}{
Number of subsamplings to perform on the original gene expression matrix data contained in SincellObject[["expressionmatrix"]]
}
  \item{num_genes}{
Number of genes to sample in each subsampling. Default is fifty percent of the genes in the original gene expression matrix.
}
  \item{cores}{
Number of threads used to paralyze the computation. Under Unix platforms, by default the function uses all cores up to 4 (to avoid possilbe issues while running on a cluster with the default parameter) detected by the operating system. Under non Unix based platforms, this parameter will be automatically set to 1.
}
}
%\details{
%The correlation between trees is computed as a spearman correlation.
%}
\value{
The SincellObject named list provided as input where following list members are added: SincellObject[["StatisticalSupportbyGeneSubsampling"]] representing the vector of spearman rank correlation values of all "num_it" iterations. Each element of SincellObject[["StatisticalSupportbyGeneSubsampling"]] represents the similarity between the connected graph resulting from one subsampling and the original graph, and it is assessed as the spearman rank correlation between the two graphs of the shortest distance for all pairs of cells.
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
## Generate some random data
Data <- matrix(abs(rnorm(3000, sd=2)),ncol=10,nrow=30)

## Initializing SincellObject named list
mySincellObject <- sc_InitializingSincellObject(Data)

## Assessmet of cell-to-cell distance matrix after dimensionality reduction 
## with Principal Component Analysis (PCA) 
mySincellObject <- sc_DimensionalityReductionObj(mySincellObject, method="PCA",dim=2)

## Cluster
mySincellObject <- sc_clusterObj (mySincellObject, clust.method="max.distance", 
  max.distance=0.5)

## Assessment of cell-state hierarchy
mySincellObject<- sc_GraphBuilderObj(mySincellObject, graph.algorithm="SST", 
  graph.using.cells.clustering=TRUE)

## Assessment statistical support by gene subsampling
mySincellObject<- sc_StatisticalSupportByGeneSubsampling(mySincellObject, 
  num_it=1000)

## To access the distribution of Spearman rank correlations:
StatisticalSupportbyGeneSubsampling<-
  mySincellObject[["StatisticalSupportbyGeneSubsampling"]]
summary(StatisticalSupportbyGeneSubsampling)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
