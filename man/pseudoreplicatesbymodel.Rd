\name{pseudoreplicatesbymodel}
\alias{pseudoreplicatesbymodel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Auxiliary function of sc_InSilicoCellsReplicatesObj function used when its parameter method="lognormal-3parameters"
}
\description{
Auxiliary function implemented in C++ making part of the sc_InSilicoCellsReplicatesObj function
}
\usage{
pseudoreplicatesbymodel(rows, colums, alpha, vargenes, meangenes, positive, f, seed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rows}{
number of rows in list "expressionmatrix" within Sincell object
}
  \item{colums}{
number of colums in list "expressionmatrix" within Sincell object
}
  \item{alpha}{
Vector containing for each gene in "expressionmatrix" the proportion of cells where expression was detected above a given threshold level (parameter "no_expr" in function sc_InSilicoCellsReplicatesObj() )
}
  \item{vargenes}{
Vector containing for each gene in "expressionmatrix" the variance of the expression levels in those cells where expression was detected above a given threshold level (parameter "no_expr" in function sc_InSilicoCellsReplicatesObj() )
}
  \item{meangenes}{
Vector containing for each gene in "expressionmatrix" the average of the expression levels in those cells where expression was detected above a given threshold level (parameter "no_expr" in function sc_InSilicoCellsReplicatesObj() )
}
  \item{positive}{
Force the new matrix to be positive. 1 for TRUE, 0 for FALSE
}
  \item{f}{
R function rnorm
}
  \item{seed}{
seed integer for random generation
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A numeric matrix is returned as described in sc_InSilicoCellsReplicatesObj when method="lognormal-3parameters"
}
%\references{
%
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
sc_InSilicoCellsReplicatesObj()
}
%\examples{

%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
