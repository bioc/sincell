\name{sc_marker2color}
\alias{sc_marker2color}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Palette of colors from the expression values of a marker gene
}
\description{
Function that transforms the expression values of a marker gene into a vector of colors that can be used as a color code for the intensity of expression. First, the function extracts the vector of values form the expression matrix row in SincellObject[["expressionmatrix"]] whose name equals the indicated marker. Then those values are transformed into a color scale in which the minimum value is assigned the color "minimum" and the maximum value the color "maximum". If relative.to.marker=TRUE, the minimum and maximum values are taken from the expression values of the marker. If relative.to.marker=FALSE, the minimum and maximum values are taken from the expression values of the entire expression matrix.
}
\usage{
sc_marker2color(SincellObject, marker, color.minimum="green", 
  color.maximum="red", relative.to.marker=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SincellObject}{
A SincellObject named list as created by function sc_InitializingSincellObject with a named member "expressionmatrix" containing a numeric matrix that represents a gene expression matrix gathering the expression levels of each single-cell in the experiment (displayed by columns) for each detected gene (displayed by rows).
}
  \item{marker}{
Name of the gene marker. It should correspond to a row name in the expression matrix in SincellObject[["expressionmatrix"]]
}
  \item{color.minimum}{
Color that will be assigned to the minimum expression value
}
  \item{color.maximum}{
Color that will be assigned to the maximum expression value
}
  \item{relative.to.marker}{
Logic indicating whether the minimum and maximum values are taken from the expression values of the marker (relative.to.marker=TRUE) or from the entire expression matrix (relative.to.marker=FALSE)
}
}
% \details{
%%  ~~ If necessary, more details than the description above ~~
% }
\value{
The function returns an array of colors in hexadecimal format
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

## Initializing SincellObject named list
mySincellObject <- sc_InitializingSincellObject(Data)

## Adding gene names to expression matrix
rownames(mySincellObject[["expressionmatrix"]]) <- 1:30

## Getting the color vector coding for the expression values of a marker gene
mymarkerColorCodeA<-sc_marker2color(mySincellObject, marker="3", 
  color.minimum="green", color.maximum="red", relative.to.marker=TRUE)
mymarkerColorCodeB<-sc_marker2color(mySincellObject, marker="5", 
  color.minimum="yellow", color.maximum="blue", relative.to.marker=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
