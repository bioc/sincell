\name{pseudoreplicatesbynoise}
\alias{pseudoreplicatesbynoise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Auxiliary function of sc_InSilicoCellsReplicatesObj function used when its parameter method="variance.deciles"
}
\description{
Auxiliary function implemented in C++ making part of the sc_InSilicoCellsReplicatesObj function when its parameter method="variance.deciles".
}
\usage{
pseudoreplicatesbynoise(originaldata, rows, colums, deciles, 
  lengthdeciles, coorsorted, vargenessorted, positive, seed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{originaldata}{
"expressionmatrix" within Sincell object: numeric matrix containing a gene expression matrix gathering the expression levels of each single-cell in the experiment (displayed by columns) for each detected gene (displayed by rows)
}
  \item{rows}{
number of rows in list "expressionmatrix" within Sincell object
}
  \item{colums}{
number of colums in list "expressionmatrix" within Sincell object
}
  \item{deciles}{
array containing the indexes indicating the limits of the deciles based on mean of gene expression
}
  \item{lengthdeciles}{
length(deciles)
}
  \item{coorsorted}{
order of permutated indexes
}
  \item{vargenessorted}{
Vector containing for each gene in "expressionmatrix" the variance of the expression levels. Order of genes corresponds to mean expression levels (increasing order) 
}
  \item{positive}{
Force the new matrix to be positive. 1 for TRUE, 0 for FALSE
}
  \item{seed}{
seed integer for random generation
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A numeric matrix is returned as described in sc_InSilicoCellsReplicatesObj 
  when method="variance.deciles"
}
%\references{
%
%}
\seealso{
sc_InSilicoCellsReplicatesObj()
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~
%\examples{

%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
