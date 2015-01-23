\name{sc_InSilicoCellsReplicatesObj}
\alias{sc_InSilicoCellsReplicatesObj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
In silico generation of replicates of individual cells
}
\description{
Function to generate in silico replicates of individual cells under different models of noise. These in silico replicates  will be used by function sc_StatisticalSupportByReplacementWithInSilicoCellReplicates() in order to provide statistical support to the connected graph in SincellObject[["cellstateHierarchy"]] assessed by function sc_GraphBuilderObj() representing a cell-state hierarchy.
}
\usage{
sc_InSilicoCellsReplicatesObj(SincellObject, method="variance.deciles", 
  multiplier=100, no_expr=0.5, cores=ifelse(detectCores()>=4, 4, detectCores()))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SincellObject}{
A SincellObject named list as created by function sc_GraphBuilderObj(), containing i) in member "cellstateHierarchy" a connected graph representing a cell-state hierarchy; and ii) in member "expressionmatrix" a numeric matrix that represents a gene expression matrix gathering the expression levels of each single-cell in the experiment (displayed by columns) for each detected gene (displayed by rows).
}
  \item{method}{
Method to generate in silico replicates of individual cells. Options are: 

i) method="variance.deciles": the mean and variance of all genes in the original gene expression matrix is assessed. Genes are assigned to classes according to the deciles of mean they belong to. Next, for a given gene g, a variance v is randomly chosen from the set of variances within the class of the gene. Then, a random value drawn from a uniform distribution U(0,v) of mean zero and variance v is added to the expression value of a gene g in a cell c. By perturbing in this way all genes in a reference cell c we obtain an in silico replicate c'. Redoing the process N times, N stochastic replicates are generated for each original cell.

ii) method= "cv2.deciles": Same as i) but a squared coefficient of variation cv2 is randomly chosen from the set of coefficient of variation values within the class of the gene (defined by deciles of mean). Then, the parameter v for the uniform distribution is assessed by v=cv2*(mean**2).

iii) method="lognormal-3parameters": random perturbations of gene expression levels are drawn from a log normal distribution log(x)~N(m,v) (where m is the mean and v the variance of the gene levels across all samples) with a third parameter alpha describing the proportion of cells where transcript expression was detected above a given threshold level (parameter "no_expr"; see Shalek et al. Nature 2014). NOTICE: This option assumes that the expression data has been log-transformed. You may want to check whether your Sincell object contains a gene expression matrix transformed that way.

}
  \item{multiplier}{
Number of in silico replicates of individual cells to generate for each cell in the original data
}
  \item{no_expr}{
Threshold value in gene expression levels of SincellObject[["expressionmatrix"]] under which a gene will be considered as non-expressed. In the case that log-transformed RPKM are used, a recomended value is 0,5.
}
  \item{cores}{
Number of threads used to paralyze the computation. Under Unix platforms, by default the function uses all cores up to 4 (to avoid possilbe issues while running on a cluster with the default parameter) detected by the operating system. Under non Unix based platforms, this parameter will be automatically set to 1.
}
}
% \details{
%%  ~~ If necessary, more details than the description above ~~
% }
\value{
The SincellObject named list provided as input where list member "InSilicoCellsReplicates" is added. SincellObject[["InSilicoCellsReplicates"]] contains the concatenation by columns of the original expression matrix together with the matrix containing the expression values per gene (by rows) of the in silico generated cells replicates (by columns).
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

## Assessmet of cell-to-cell distance matrix after dimensionality reduction 
## with Principal Component Analysis (PCA) 
mySincellObject <- sc_DimensionalityReductionObj(mySincellObject, method="PCA",dim=2)

## Cluster
mySincellObject <- sc_clusterObj (mySincellObject, clust.method="max.distance", 
  max.distance=0.5)

## Assessment of cell-state hierarchy
mySincellObject<- sc_GraphBuilderObj(mySincellObject, graph.algorithm="SST", 
  graph.using.cells.clustering=TRUE)

## In silico generation of replicates of individual cells
mySincellObject <- sc_InSilicoCellsReplicatesObj(mySincellObject, 
  method="variance.deciles", multiplier=100, no_expr=0.5)

# To access the in silico generated cells replicates
InSilicoCellsReplicates<-mySincellObject[["InSilicoCellsReplicates"]]

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
