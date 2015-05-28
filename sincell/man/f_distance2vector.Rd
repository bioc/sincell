\name{f_distance2vector}
\alias{f_distance2vector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conversion of the lower triangular matrix of a distance matrix into an array
}
\description{
Auxiliary function to convert of the lower triangular matrix of a distance matrix into an array
}
\usage{
f_distance2vector(distance)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distance}{
A distance matrix or a symmetric matrix
}
}
% \details{
%%  ~~ If necessary, more details than the description above ~~
% }
\value{
Array resulting from the concatenation of the rows of the lower triangular matrix of the input symetric matrix. Array length is n*(n-1)/2, where n is the number of rows of the symetric matrix.
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
## Generate some data
Data <- matrix(rnorm(300),ncol=10,nrow=30)

## Calculate distance matrix and transform its lower triangular matrix into a one 
## dimensional array
d <- f_distance2vector(as.matrix(dist(Data)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
