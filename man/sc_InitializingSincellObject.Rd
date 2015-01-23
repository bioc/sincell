\name{sc_InitializingSincellObject}
\alias{sc_InitializingSincellObject}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to initialize a sincell object
}
\description{
Function initializes a named list with a unique member so-called "expressionmatrix" containing the input gene expression matrix. Genes with a variance equal to zero are filtered out from the gene expression matrix at this step.
}
\usage{
sc_InitializingSincellObject(BaseData)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BaseData}{
A numeric matrix representing a gene expression matrix gathering the normalized expression levels of each single-cell in the experiment (displayed by columns) for each detected gene (displayed by rows).
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
a named list: list(expressionmatrix=BaseData)
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
