\name{sc_GraphBuilderObj}
\alias{sc_GraphBuilderObj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Graph building function for assessment of cell-state hierarchies
}
\description{
Function to build a connected graph from a cell-to-cell distance matrix that will be regarded as a cell-state hierarchy. Three algorithms are available: the Minimum Spanning Tree (MST), the Maximum Similarity Spanning Tree (SST) and the Iterative Mutual Clustering Graph (IMC). Optionally, algorithms in sc_GraphBuilderObj can use a precalculated clustering of cells. In the case of MST, this is used to overlay connections between pairs of cells belonging to the same cluster. In the case of SST, clusters of cells are treated as atomic elements in the graph-building process together with non-clustered cells. By definition, IMC builds a connected graph through iterations on the clustering results produced the K-Mutual Nearest Neighbour (K-MNN) algorithm.
}
\usage{
sc_GraphBuilderObj(SincellObject, graph.algorithm="MST", 
  graph.using.cells.clustering=FALSE,k=3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SincellObject}{
A SincellObject named list as created by function sc_distanceObj() or 
sc_DimensionalityReductionObj(), containing in member "cell2celldist" a distance matrix representing a cell-to-cell distance matrix assessed on a gene expression matrix with a metric of choice.
}
  \item{graph.algorithm}{
Graph building algorithm to be used: the Minimum Spanning Tree 
(graph.algorithm="MST"), the Maximum Similarity Spanning Tree 
(graph.algorithm="SST") and the Iterative Mutual Clustering Graph 
(graph.algorithm="IMC").
}
  \item{graph.using.cells.clustering}{
If graph.using.cells.clustering=TRUE and graph.algorithm="MST" or graph.algorithm="MST", a precalculated clustering of cells is used. The clustering of cells is taken from SincellObject[["cellsClustering"]] as calculated by function sc_clusterObj().
}
  \item{k}{
If IMC algorithm is selected, the number of nearest neighbors used in the underlying K-Mutual Nearest Neighbour (K-MNN) algorithm is set to k.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The SincellObject named list provided as input where following list members are added: 
"cellstateHierarchy"=cellstateHierarchy,
"graph.algorithm"=graph.algorithm,
"graph.using.cells.clustering"=graph.using.cells.clustering. 
The member "cellstateHierarchy" contains an igraph graph object (see "igraph" R package documentation) representing a totally connected graph built with the indicated parameters.
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
## Generate some data
## Generate some random data
Data <- matrix(abs(rnorm(3000, sd=2)),ncol=10,nrow=300)

## Initializing SincellObject named list
mySincellObject <- sc_InitializingSincellObject(Data)

## Assessmet of cell-to-cell distance matrix after dimensionality reduction with 
## Principal Component Analysis (PCA) 
mySincellObject <- sc_DimensionalityReductionObj(mySincellObject, method="PCA",dim=2)

## Cluster
mySincellObject <- sc_clusterObj (mySincellObject, clust.method="max.distance", 
  max.distance=0.5)

## Assessment of cell-state hierarchy
mySincellObject<- sc_GraphBuilderObj(mySincellObject, graph.algorithm="MST", 
  graph.using.cells.clustering=FALSE)
mySincellObject<- sc_GraphBuilderObj(mySincellObject, graph.algorithm="SST", 
  graph.using.cells.clustering=TRUE)
mySincellObject<- sc_GraphBuilderObj(mySincellObject, graph.algorithm="IMC")

## To access the totally connected graph (igraph object)
cellstateHierarchy<-mySincellObject[["cellstateHierarchy"]]

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
