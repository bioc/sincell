\name{sc_DimensionalityReductionObj}
\alias{sc_DimensionalityReductionObj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dimensionality reduction of an expression matrix
}
\description{
Function to perform a dimensionality reduction upon the original gene expression matrix data through a method of choice, either linear or no-linear, among the following: Principal Component Analysis (PCA), Independent Component Analysis (ICA), t-Distributed Stochastic Neighbor Embedding (tSNE), classical Multidimensional Scaling (MDS) and non-metric Multidimensional Scaling.
}
\usage{
sc_DimensionalityReductionObj(SincellObject, method="PCA", dim=2, 
MDS.distance="spearman", bins=c(-Inf,0,1,2,Inf),tsne.perplexity=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SincellObject}{
A SincellObject named list as created by function sc_InitializingSincellObject with a named member "expressionmatrix" containing a numeric matrix that represents a gene expression matrix gathering the expression levels of each single-cell in the experiment (displayed by columns) for each detected gene (displayed by rows).
}
  \item{method}{
Dimensionality reduction algorithm to be used.  Options are: Principal Component Analysis (method="PCA"), Independent Component Analysis (method="ICA"; using fastICA() function in fastICA package), t-Distributed Stochastic Neighbor Embedding (method="tSNE"; using Rtsne() function in Rtsne package with parameters tsne.perplexity=1 and tsne.theta=0.25), classical Multidimensional Scaling (method="classical-MDS"; using the cmdscale() function) and  non-metric Multidimensional Scaling (method="nonmetric-MDS";using the isoMDS() function in MASS package). if method="PCA" is chosen, the proportion of variance explained by each of the principal axes is plotted.
}
  \item{dim}{
Number of dimensions in low-dimensional space to be retained. Default is dim=2.
}
  \item{MDS.distance}{
Distance method to be used if method="classical-MDS" or method="nonmetric-MDS" is selected. The available distances are the Euclidean distance (method="euclidean"), Manhattan distance  (also called L1 distance, method="L1"), distance based on Pearson (method="pearson") or Spearman (method="spearman") correlation coefficients, and distance based on Mutual Information (method="MI"). Intervals used to assess Mutual Information are indicated in the parameter “bins” (see below).
}
  \item{bins}{
Intervals used to discretize the data in the case that Mutual Information distance (MDS.distance="MI") is selected.
}
  \item{tsne.perplexity}{
perplexity parameter for tSNE algorithm. We refer the reader to the Frequently Asked Questions in http://homepage.tudelft.nl/19j49/t-SNE.html
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A SincellObject named list whose members are: 
expressionmatrix=SincellObject[["expressionmatrix"]],
cellsLowDimensionalSpace=cellsLowDimensionalSpace,
cell2celldist=distance,method=method, dim=dim, 
MDS.distance=MDS.distance, bins=bins, 
where cellsLowDimensionalSpace contains the coordinates of each cell (by columns) in each low dimensional axis (by rows), and "cell2celldist" contains the numeric matrix representing the cell-to-cell distance matrix assessed in low dimensional space
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
Data <- matrix(abs(rnorm(3000, sd=2)),ncol=10,nrow=300)

## Initializing SincellObject named list
mySincellObject <- sc_InitializingSincellObject(Data)

## To access the gene expression matrix
expressionmatrix<-mySincellObject[["expressionmatrix"]]

## Dimensionality reduction
# Principal Component Analysis (PCA)
mySincellObject <- sc_DimensionalityReductionObj(mySincellObject, method="PCA",dim=2)
# Independent Component Analysis (ICA)
mySincellObject <- sc_DimensionalityReductionObj(mySincellObject, method="ICA",dim=2)
# t-Distributed Stochastic Neighbor Embedding (t-SNE)
mySincellObject <- sc_DimensionalityReductionObj(mySincellObject, method="tSNE",dim=2)
# Classic Multidimensional Scaling (classic-MDS).
mySincellObject <- sc_DimensionalityReductionObj(mySincellObject, method="classical-MDS",dim=2)
# Non-metric Multidimensional Scaling (nonmetric-MDS).
mySincellObject <- sc_DimensionalityReductionObj(mySincellObject, method="nonmetric-MDS",dim=2)

## To access the coordinates of cells (by columns) in low dimensional space (axes by rows)
cellsLowDimensionalSpace<-mySincellObject[["cellsLowDimensionalSpace"]]

## To access the cell-to-cell distance matrix assessed on low dimensional space
cell2celldist<-mySincellObject[["cell2celldist"]]

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
