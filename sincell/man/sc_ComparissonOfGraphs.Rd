\name{sc_ComparissonOfGraphs}
\alias{sc_ComparissonOfGraphs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Comparisson of graphs
}
\description{
Function to assess a distance matrix comparing the graphs from Sincell objects that were generated with function sc_GraphBuilderObj(). The distance between two graphs is assessed as 1 minus their similarity, which is calculated as the spearman rank correlation between the two graphs of the shortest distance for all pairs of cells. Cell-state hierarchies are igraph graph objects (see "igraph" R package documentation) representing a totally connected graph.
}
\usage{
sc_ComparissonOfGraphs(cellstateHierarchy1,cellstateHierarchy2, ...,
  graph.names=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cellstateHierarchy1}{
A first cell-state hierarchy as created by function sc_GraphBuilderObj() on a SincellObject.
}
  \item{cellstateHierarchy2}{
A second cell-state hierarchy as created by function sc_GraphBuilderObj() on a SincellObject.
}
  \item{...}{
Further cell-state hierarchies
}
  \item{graph.names}{
A vector of characters indicating the names of the cell-state hierarchies provided as arguments.
}
}
% \details{
%%  ~~ If necessary, more details than the description above ~~
% }
\value{
A distance matrix comparing the graphs.
}
% \references{
%% ~put references to the literature/web site here ~
% }
% \author{
%%  ~~who you are~~
% }
% \note{
%%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
% }
\examples{

## Generate some random data
Data <- matrix(abs(rnorm(3000, sd=2)),ncol=10,nrow=30)

## Initializing SincellObject
mySincellObject <- sc_InitializingSincellObject(Data)

## Assessmet of cell-to-cell distance matrix after dimensionality reduction 
## with Principal Component Analysis (PCA), with Independent Component 
## Analysis (ICA), or with non-metric Multidimensional Scaling (nonmetric-MDS)
mySincellObject_PCA <- sc_DimensionalityReductionObj(mySincellObject, 
  method="PCA",dim=2)
mySincellObject_ICA <- sc_DimensionalityReductionObj(mySincellObject, 
  method="ICA",dim=2)
mySincellObject_classicalMDS <- sc_DimensionalityReductionObj(mySincellObject, 
  method="classical-MDS",dim=2)
mySincellObject_nonmetricMDS <- sc_DimensionalityReductionObj(mySincellObject, 
  method="nonmetric-MDS",dim=2)

## Assessment of cell-state hierarchy
mySincellObject_PCA<- sc_GraphBuilderObj(mySincellObject_PCA, 
  graph.algorithm="SST")
mySincellObject_ICA<- sc_GraphBuilderObj(mySincellObject_ICA, 
  graph.algorithm="SST")
mySincellObject_classicalMDS<- sc_GraphBuilderObj(mySincellObject_classicalMDS, 
  graph.algorithm="SST")
mySincellObject_nonmetricMDS<- sc_GraphBuilderObj(mySincellObject_nonmetricMDS, 
  graph.algorithm="SST")

## Comparisson of hierarchies obtained from different methods
myComparissonOfGraphs<-sc_ComparissonOfGraphs(
  mySincellObject_PCA[["cellstateHierarchy"]],
  mySincellObject_ICA[["cellstateHierarchy"]],
  mySincellObject_classicalMDS[["cellstateHierarchy"]],
  mySincellObject_nonmetricMDS[["cellstateHierarchy"]],
  graph.names=c("PCA","ICA","classicalMDS","nonmetricMDS")
  )
plot(hclust(myComparissonOfGraphs))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
