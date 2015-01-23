\name{knnalgorithm}
\alias{knnalgorithm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Auxiliary function for KNN and IMC algorithms
}
\description{
Auxiliary function
}
\usage{
knnalgorithm(distance, mutual, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distance}{
distance matrix
}
  \item{mutual}{
logical specifying if the connection between neighbors must be mutual
}
  \item{k}{
maximum order of neighbors
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
An adjacency matrix is returned.
}
%\references{
%See sc_GraphBuilderObj()
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
sc_GraphBuilderObj(), sc_clusterObj()
}
%\examples{

%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
