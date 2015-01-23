\name{sc_distanceObj}
\alias{sc_distanceObj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Assessment of a cell-to-cell distance matrix with a metric of choice
}
\description{
Function to assess a cell-to-cell distance matrix from a gene expression matrix with a metric of choice among the following: Euclidean distance, Mutual Information, L1 distance (Manhattan distance), Pearson correlation or Spearman correlation.
}
\usage{
sc_distanceObj(SincellObject, method="euclidean", bins=c(-Inf,0,1,2,Inf))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SincellObject}{
A SincellObject named list as created by function sc_InitializingSincellObject with a named member "expressionmatrix" containing a numeric matrix that represents a gene expression matrix gathering the expression levels of each single-cell in the experiment (displayed by columns) for each detected gene (displayed by rows).
}
  \item{method}{
Distance method to be used. The available distances are the Euclidean distance (method="euclidean"), Manhattan distance  (also called L1 distance, method="L1"), cosine distance (method="cosine") , distance based on Pearson (method="pearson") or Spearman (method="spearman") correlation coefficients, and distance based on Mutual Information (method="MI"). Intervals used to assess Mutual Information are indicated in the parameter “bins” (see below).
}
  \item{bins}{
Intervals used to discretize the data in the case that Mutual Information distance (method="MI") is selected.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A SincellObject named list whose members are: 
expressionmatrix=SincellObject[["expressionmatrix"]],
cell2celldist=cell2celldist,method=method,bins=bins, 
where "cell2celldist" contains the numeric matrix representing the cell-to-cell distance matrix assessed by sc_distanceObj with the indicated parameters
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

## Distance
mySincellObject<-sc_distanceObj(mySincellObject)
mySincellObject<-sc_distanceObj(mySincellObject, method="MI",bins=c(-Inf,0,2,4,6,8,Inf))
mySincellObject<-sc_distanceObj(mySincellObject, method="spearman")

## To access the cell-to-cell distance matrix
cell2celldist<-mySincellObject[["cell2celldist"]]

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
