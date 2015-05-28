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
  dispersion.statistic = NULL, multiplier=100, no_expr=0.5, 
  LogTransformedData = T, baseLogTransformation=exp(1),
  pseudocounts.added.before.log.transformation=1, 
  cores=ifelse(detectCores()>=4, 4, detectCores()))
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

iv) method="negative.binomial": random perturbations of gene expression levels are drawn from a negative binomial (NB) distribution NB(m,r), where m is the mean and r is the size (i.e. the dispersion parameter). Under this parameterization, the variance is v=m+(m^2/r), therefore r=m^2/(v-m). For each gene, its mean m is estimated from the expression levels of the expression matrix. There are three alternative ways of defining the variance v for each gene, which are indicated in parameter dispersion.statistic. Some works has found that, for most genes, the variability observed among their expression levels across individual cells was better described by a negative binomial (NB) distribution rather than a lognormal distribution (Grün et al., 2014). Grün and colleagues used NB distribution to model not only technical noise but also true biological gene expression noise. Their assumption was that endogenous mRNA abundance follows a NB as supported by a physical model of bursting expression (Raj et al., 2006). A negative binomial noise model was also adopted in (Zeisel et al., 2015). As pointed out in these works, NB is frequently used to model overdispersed count data and has been previously used for bulk RNA-seq data (Anders and Huber, 2010; Robinson et al., 2010). We recommend this approach only if normalized count data is used (i.e. not length-normalized RPKM/FPKM). Sincell can follow an NB distribution parameterized on the observed gene expression levels to generate random perturbations and produce in silico cell replicates accordingly. If log-transformed normalized counts are used, Sincell would unlog the perturbed data through a NB and afterwards will redo the log trans-formation. Parameters "LogTransformedData", "baseLogTransformation", "pseudocounts.added.before.log.transformation", should be indicated to help Sincell perform de unlog and log in a consistent way with user's transformations.

}
  \item{dispersion.statistic}{
  if parameter method=="negative.binomial", there are three alternative ways of definining the variance v that will be used to parameterize the negative binomial distribution
a) dispersion.statistic==NULL ; variance is estimated from the input expression levels of the expression matrix
b) is.numeric(dispersion.statistic) ; vector provided by the user of length equal to the number of genes in the input expression matrix. This vector should contain cv2 estimates reflecting e.g. estimated technical noise. Estimates of technical noise for each gene can be obtained by modeling the dependence of the coefficient of variation (cv2) of spike-in molecules as a function of their average expression. For instance, in Brennecke et al. 2013, for each technical gene i (e.g. the spike-ins), the sample mean (m) and sample variance of its normalized counts are estimated. Then, the observed squared coefficients of variation (cv2) are fitted against the sample mean (m) with a generalized linear model of the gamma family with identity link and parameterization cv2=a1/m+ alpha0. Applying the fitted formula to the sample mean expression levels of a gene provides an estimate of cv2 arising from technical noise. Sincell permits the incorporation of a technical cv2 estimate per gene in the assessment of in silico cell replicates based on normalized counts (i.e. following the previously described negative binomial distribution whose dispersion is parameterized using the estimated technical cv2). 

c) dispersion.statistic!="cv2.fitted.to.data" ; alternatively, in the absence of spike-in molecules, Sincell implements the fit described in Brennecke et al. 2013 using the cv2 and m values of all genes in the input expression matrix to provide a surrogate of technical noise estimates. However, this alternative should not be used if the user has previously followed our recommendation in Section 1 of using such an approach to identify highly variable genes in order to decrease the size of the input matrix (http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html; Section "Identifying highly variable genes").

}
  \item{multiplier}{
Number of in silico replicates of individual cells to generate for each cell in the original data
}
  \item{no_expr}{
Threshold value in gene expression levels of SincellObject[["expressionmatrix"]] under which a gene will be considered as non-expressed. In the case that log-transformed RPKM are used, a recomended value is 0,5.
}
  \item{LogTransformedData}{
T (TRUE) or F (FALSE). Indicating whether the input expression matrix used to assessed hierarchies was previously logtransformed
}
  \item{baseLogTransformation}{
if LogTransformedData==T, the base used for the logtransformation
}
  \item{pseudocounts.added.before.log.transformation}{
if LogTransformedData==T, the number of pseudocounts added to the normalized count data before performing logtransformation
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
