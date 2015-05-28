\name{sstalgorithm}
\alias{sstalgorithm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Auxiliary function for SST algorithm within function sc_GraphBuilderObj()
}
\description{
Auxiliary function implemented in C++ making part of the SST algorithm in function 
sc_GraphBuilderObj()
}
\usage{
sstalgorithm(membership, num_cells, distance)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{membership}{
a numeric array
}
  \item{num_cells}{
total number of cells in the sample
}
  \item{distance}{
a distance matrix
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A numeric array of length 3 is returned. The first element of the array is the minimum distance, and the second and third ones are the coordinates.
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
sc_GraphBuilderObj()
}
%\examples{

%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
