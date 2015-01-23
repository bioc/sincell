\name{sc_AssociationOfCellsHierarchyWithAGeneSet}
\alias{sc_AssociationOfCellsHierarchyWithAGeneSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Association of a cell-state hierarchy with a functional gene set
}
\description{
First, this function assesses a cell-state hierarchy where only the expression levels of the genes in a given functional gene set are considered. Second, it calculates the similarity of that hierarchy with the one assessed by function sc_GraphBuilderObj() on the initial gene expression matrix. Third it provides an empirical p-value of the observed similarity between the two hierarchies. The hierarchy resulting when considering only the genes in the gene set is assessed with exactly the same parameters used to obtain the reference hierarchy. The similarity between the two hierarchies is computed as the spearman rank correlation between the two graphs of the shortest distance for all pairs of cells. The empirical p-value is calculated from a distribution of similarities resulting from random samplings of gene sets of the same size.
}
\usage{
sc_AssociationOfCellsHierarchyWithAGeneSet(SincellObject,GeneSet, 
  minimum.geneset.size=50,p.value.assessment=TRUE,
  spearman.rank.threshold=0.5,num_it=1000, 
  cores=ifelse(detectCores()>=4, 4, detectCores()))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{SincellObject}{
A SincellObject named list as created by function sc_GraphBuilderObj(), containing in member "cellstateHierarchy" a connected graph representing a cell-state hierarchy.
}
\item{GeneSet}{
A character vector containing the gene names of a functional gene set. Gene names should be of the same type as those used in the gene expression matrix.
}
\item{minimum.geneset.size}{
Minimum number of genes from the gene set that should be present in the original gene expression matrix and that have a non-zero variance across cells. If that overlap is is lower than this parameter, the association will not be computed.
}
\item{p.value.assessment}{
A logical value indicating whether an empirical p-value of the similarity should be calculated.
}
\item{spearman.rank.threshold}{
The minimum value of the spearman rank correlation that the two hierarchies should have to allow computation of an empirical p-value. This limit is set in order to avoid an extra computation time invested in getting an empirical p-value for a low correlation not worthy of consideration.
}
\item{num_it}{
Number of subsamplings to perform on the original gene expression matrix data contained in SincellObject[["expressionmatrix"]] to obtain the empirical p-value
}
  \item{cores}{
Number of threads used to paralyze the computation. Under Unix platforms, by default the function uses all cores up to 4 (to avoid possilbe issues while running on a cluster with the default parameter) detected by the operating system. Under non Unix based platforms, this parameter will be automatically set to 1.
}
}
%\details{
%The correlation between trees is computed as a spearman correlation.
%}
\value{
The SincellObject named list provided as input where following list members are added: The similarity between the reference hierarchy and the hierarchy obtained from the gene set, stored in SincellObject[["AssociationOfCellsHierarchyWithAGeneSet"]]; and its empirical p-value, stored in SincellObject[["AssociationOfCellsHierarchyWithAGeneSet.pvalue"]]
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
Data <- matrix(abs(rnorm(3000, sd=2)),ncol=10,nrow=50)
rownames(Data)<-character(dim(Data)[1])

## Generate gene names from index
for (i in 1:dim(Data)[1]){rownames(Data)[i]<-as.character(i)}

## Generate a hypothetical gene list from the first 10 gene names
myGeneSet<-rownames(Data)[1:10]

## Initializing SincellObject named list
mySincellObject <- sc_InitializingSincellObject(Data)

## Assessmet of cell-to-cell distance matrix after dimensionality reduction with 
## Principal Component Analysis (PCA) 
mySincellObject <- sc_DimensionalityReductionObj(mySincellObject, method="PCA",dim=2)

## Cluster
mySincellObject <- sc_clusterObj (mySincellObject, clust.method="max.distance", 
  max.distance=0.5)

## Assessment of cell-state hierarchy
mySincellObject<- sc_GraphBuilderObj(mySincellObject, graph.algorithm="SST", 
  graph.using.cells.clustering=TRUE)

## Assessment of association of the hierarchy with a gene set
mySincellObject<-sc_AssociationOfCellsHierarchyWithAGeneSet(mySincellObject,
  myGeneSet, minimum.geneset.size=9,p.value.assessment=TRUE,
  spearman.rank.threshold=0.5,num_it=1000)

## To access the similarity between the reference hierarchy and the hierarchy obtained 
## from the gene set
myAssociationOfCellsHierarchyWithGeneSet<-
  mySincellObject[["AssociationOfCellsHierarchyWithAGeneSet"]]
myAssociationOfCellsHierarchyWithGeneSet.pvalue<-
  mySincellObject[["AssociationOfCellsHierarchyWithAGeneSet.pvalue"]]

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
