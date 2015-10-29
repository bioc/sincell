## Dependences
# install.packages(c("igraph", "entropy", "scatterplot3d", "MASS", "TSP", "ggplot2", "reshape2", "fields", "spam", "grid", "maps", "proxy","parallel", "Rtsne", "fastICA"))
## Optional
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("biomaRt", "GO.db", "KEGGgraph", "KEGG.db"))

# library(igraph)
# library(entropy)
# library(scatterplot3d)
# library(MASS)
# library(TSP)
# library(ggplot2)
# library(reshape2)
# library(fields)
# library(proxy)
# library(parallel)
# library(Rtsne)
# library(fastICA)
# library(cluster)
# library(biomaRt)
# library(statmod)

# Loading required packages: spam, grid, maps
#===================================================
## Creating sincell object
sc_InitializingSincellObject<-function(BaseData){
  if(is.matrix(BaseData)==FALSE){
    errormesseage<-paste(deparse(substitute(BaseData))," is not a matrix")
    stop(errormesseage)
  }else{
    SincellObject<-list(
      "expressionmatrix"=as.matrix(unique(BaseData[apply(BaseData,1,var)>0,]))
    )
    return(SincellObject)
  }
}
#===================================================
## Auxiliary function for knn algorithm
knnalgorithm <- function(distance, mutual=TRUE, k=3){
  n <- nrow(distance)
  A <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    d <- sort(distance[i, ])
    A[i, distance[i, ] <= d[k + 1]] <- 1
  }
  diag(A) <- 0
  if (mutual) {
    for (i in 1:n) {
      A[i, ] <- A[i, ] & A[, i]
      A[, i] <- A[i, ]
    }
  }
  A <- A*distance
  g <- graph.adjacency(A, weighted=TRUE, mode = "undirected")
  #g <- as.undirected(g, mode="collapse")
  return(g)
}
#===================================================
## Distance
sc_distanceObj <- function(SincellObject, method="euclidean", bins=c(-Inf,0,1,2,Inf)){
  if(is.list(SincellObject)==FALSE){
    errormesseage<-paste(deparse(substitute(SincellObject)),"is not a list")
    stop(errormesseage)
  }else if(is.null(SincellObject[["expressionmatrix"]])==TRUE){
    errormesseage<-paste(deparse(substitute(SincellObject)),"does not contain a valid expressionmatrix member. Please initialize using function sc_InitializingSincellObject()")
    stop(errormesseage)
  }else if(!method %in% c("euclidean", "L1", "MI", "pearson", "spearman", "cosine")){
    errormesseage<-paste(method,"is not a valid option for parameter method. Please read documentation of sc_distanceObj() in Sincell manual")
    stop(errormesseage)
  }else{
    SincellObject<-SincellObject["expressionmatrix"]
    BaseData<-SincellObject[["expressionmatrix"]]
    if(method=="euclidean"){
      BaseData <- t(BaseData)
      distance <- rdist(BaseData)
      rownames(distance)<-rownames(BaseData)
      colnames(distance)<-rownames(BaseData)
    }
    if(method=="L1"){
      BaseData <- t(BaseData)
      distance <- as.matrix(stats::dist(BaseData, method = "manhattan"))
    }
    if(method=="MI"){
      bins.labels <- c() 
      for(i in 1:(length(bins)-1)){bins.labels <- append(bins.labels, paste(bins[i],bins[i+1],sep="-"))}
      distance <- matrix(as.numeric(cut(array(unlist(BaseData)),breaks = bins,labels = bins.labels)),ncol=ncol(BaseData),nrow=nrow(BaseData))
      dist_bins <- matrix(0, ncol=ncol(BaseData), nrow=ncol(BaseData))
      for(i in 1:(ncol(distance)-1)){
        for(j in (i+1):ncol(distance)){
          dist_bins[i,j] <- mi.empirical(rbind(distance[,i],distance[,j]),unit="log2")
          dist_bins[j,i] <- dist_bins[i,j]
        }
      }
      distance <- dist_bins
      rownames(distance)<-colnames(BaseData)
      colnames(distance)<-colnames(BaseData)
    }
    if(method=="pearson"){
      distance <- 1-cor(BaseData, method="pearson")
    }
    if(method=="spearman"){
      distance <- 1-cor(BaseData, method="spearman")
    }
    if(method=="cosine"){
      y <- t(BaseData) %*% BaseData
      distance <- 1 - y / (sqrt(diag(y)) %*% t(sqrt(diag(y))))
    }
    SincellObject[["cell2celldist"]]=distance
    SincellObject[["method"]]=method
    SincellObject[["bins"]]=bins
    return(SincellObject)
  }
}

#===================================================
## DImensional reduction
sc_DimensionalityReductionObj <- function(SincellObject, method="PCA", dim=2, MDS.distance="spearman", bins=c(-Inf,0,1,2,Inf),tsne.perplexity=1,tsne.theta=0.25){
  if(is.list(SincellObject)==FALSE){
    errormesseage<-paste(deparse(substitute(SincellObject)),"is not a list")
    stop(errormesseage)
  }else if(is.null(SincellObject[["expressionmatrix"]])==TRUE){
    errormesseage<-paste(deparse(substitute(SincellObject)),"does not contain a valid expressionmatrix member. Please initialize using function sc_InitializingSincellObject()")
    stop(errormesseage)
  }else if(!method %in% c("PCA", "classical-MDS", "nonmetric-MDS", "tSNE", "ICA")){
    errormesseage<-paste(method,"is not a valid option for parameter method. Please read documentation of sc_DimensionalityReductionObj() in Sincell manual")
    stop(errormesseage)
  }else if(!MDS.distance %in% c("euclidean", "L1", "MI", "pearson", "spearman")){
    errormesseage<-paste(method,"is not a valid option for parameter MDS.distance . Please read documentation of sc_DimensionalityReductionObj() in Sincell manual")
    stop(errormesseage)
  }else{ 
    SincellObject<-SincellObject["expressionmatrix"]
    BaseData<-SincellObject[["expressionmatrix"]]
    if(method=="PCA"){
      newData <- prcomp(t(BaseData), center = TRUE, scale. = TRUE)
      # plot(summary(newData)$importance[2,],las=1, main="Proportion of variance explained by\neach principal axis", ylab="Proportion of variance",xlab="Principal axes",pch=16,ylim=c(0,0.25))
      SincellObject[["EigenValuesPCA"]]<-summary(newData)$importance[2,]
      newData <- newData$x[,1:dim]
    }
    if(method=="classical-MDS"){
      SincellObject<-sc_distanceObj(SincellObject, method=MDS.distance, bins=bins)
      
      newData <- cmdscale((stats::as.dist(SincellObject[["cell2celldist"]])),eig=TRUE, k=dim)$points
    }
    if(method=="nonmetric-MDS"){
      SincellObject<-sc_distanceObj(SincellObject, method=MDS.distance, bins=bins)
      newData <- isoMDS((stats::as.dist(SincellObject[["cell2celldist"]])), k=dim)$points
    }
    if(method=="tSNE"){
      newData <- Rtsne(t(BaseData),dims=dim, initial_dims=ncol(BaseData), perplexity=tsne.perplexity, theta=tsne.theta, check_duplicates=FALSE)$Y
      newData <- newData[,1:dim]
    }
    if(method=="ICA"){
      newData <- fastICA(t(BaseData), n.comp=dim, alg.typ = "parallel", fun = "logcosh", alpha  = 1,
                         method = "C", row.norm = FALSE, maxit = 200,
                         tol = 0.0001, verbose = FALSE)$S
      newData <- newData[,1:dim]
    }
    rownames(newData)<-colnames(BaseData)
    cellsLowDimensionalSpace<-t(as.matrix(newData)); 
    cell2celldist <- rdist(t(cellsLowDimensionalSpace))
    rownames(cell2celldist)<-colnames(cellsLowDimensionalSpace)
    colnames(cell2celldist)<-colnames(cellsLowDimensionalSpace)

    SincellObject[["cellsLowDimensionalSpace"]]=cellsLowDimensionalSpace
    SincellObject[["cell2celldist"]]=cell2celldist
    SincellObject[["method"]]=method
    SincellObject[["dim"]]=dim
    SincellObject[["MDS.distance"]]=MDS.distance
    SincellObject[["bins"]]=bins
    return(SincellObject)
  }    
}
#===================================================
## Extract color from marker
sc_marker2color <- function(SincellObject, marker, color.minimum="green", color.maximum="red", relative.to.marker=TRUE){
  if(is.list(SincellObject)==FALSE){
    errormesseage<-paste(deparse(substitute(SincellObject)),"is not a list")
    stop(errormesseage)
  }else if(is.null(SincellObject[["expressionmatrix"]])==TRUE){
    errormesseage<-paste(deparse(substitute(SincellObject)),"does not contain a valid expressionmatrix member. Please initialize using function sc_InitializingSincellObject()")
    stop(errormesseage)
  }else if(!is.character(marker)|!is.vector(marker)|!(length(marker)==1)){
    errormesseage<-paste("marker does not contain a valid value. Please provide a string with a single gene name")
    stop(errormesseage)
  }else if(!(marker %in% rownames(SincellObject[["expressionmatrix"]]))){
    errormesseage<-paste(deparse(substitute(SincellObject)),"does not contain any gene named",marker,". Please check that you are using the same type of gene identifiers")
    stop(errormesseage)
  }else{
    breaks=1000
    BaseData<-SincellObject[["expressionmatrix"]]
    colPal <- colorRampPalette(c(color.minimum,color.maximum))
    if(is.logical(relative.to.marker) && relative.to.marker==TRUE){
      marker.color <- colPal(breaks)[as.numeric(cut(BaseData[rownames(BaseData)==marker,],breaks = breaks))]
    }else{
      maxVal<-max(BaseData)
      minVal<-min(BaseData)
      markerVals<-BaseData[rownames(BaseData)==marker,]
      markerValsScaled<-floor((markerVals-minVal)*(1/(maxVal-minVal))*(breaks-1)+1)
      marker.color <- colPal(breaks)[markerValsScaled]
    }
    return(marker.color)
  }
}
#===================================================
## Clusters methods
sc_clusterObj <- function(SincellObject, clust.method="knn", mutual=TRUE, k=3, max.distance=0, shortest.rank.percent=10){
  if(is.list(SincellObject)==FALSE){
    errormesseage<-paste(deparse(substitute(SincellObject))," is not a list")
    stop(errormesseage)
  }else if(is.null(SincellObject[["expressionmatrix"]])==TRUE){
    errormesseage<-paste(deparse(substitute(SincellObject)),"does not contain a valid expressionmatrix member. Please initialize using function sc_InitializingSincellObject()")
    stop(errormesseage)
  }else if(is.null(SincellObject[["cell2celldist"]])==TRUE){
    errormesseage<-paste(deparse(substitute(SincellObject)),"does not contain a valid cell2celldist member. Please assess a cell-2-cell distance matrix using either function sc_distanceObj() or function sc_DimensionalityReductionObj()")
    stop(errormesseage)
  }else{ 
    SincellObject<-SincellObject[c("expressionmatrix","cellsLowDimensionalSpace","cell2celldist","method","dim","MDS.distance","bins")]
    distance<-SincellObject[["cell2celldist"]]
    if(clust.method=="knn"){
      g <- knnalgorithm(distance, mutual=TRUE, k=k)
    }
    if(clust.method=="max.distance"){
      distance[ row(distance) == col(distance) ] <- 0
      g <- graph.adjacency(distance ,weighted=TRUE,mode="undirected")
      g <- delete.edges(g, which(E(g)$weight > max.distance))
      g <- delete.edges(g, which(E(g)$weight == 0))
    }
    if(clust.method=="percent"){
      distance[ row(distance) == col(distance) ] <- 0
      g <- graph.adjacency(distance ,weighted=TRUE,mode="undirected")
      g <- delete.edges(g, which(E(g)$weight > sort(E(g)$weight)[shortest.rank.percent/100*length(E(g)$weight)]))
      g <- delete.edges(g, which(E(g)$weight == 0))
    }
#     # NOTICE: KNNG is implemented here but commented. It would require adding to the function header the parameters: train, test, cl
#     if(clust.method=="knng"){
#       classification <- class::knn(train=train, test=test, cl=cl, k=k, l=(k%/%2)+1)
#       g <- graph.adjacency(distance, weighted=TRUE, mode="undirected")
#       g <- delete.edges(g, E(g))
#       for(i in unique(classification)){
#         if(!is.na(i)){
#           tmp <- which(classification==i)
#           tmp <- rbind(tmp, classification[tmp])
#           edge_list <- array(unlist(tmp))
#           edge_weights <- c()
#           for(j in 1:ncol(tmp)){
#             edge_weights[j] <- distance[tmp[1,j],tmp[2,j]]
#           }
#           g <- add.edges(g, edge_list, attr=list(weight=edge_weights))
#         }
#       }
#     }
    if(
      # "hclust.method" can take one of "ward.D", "ward.D2", "single", "complete","average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
      (clust.method=="ward.D")|(clust.method=="ward.D2")|
      (clust.method=="single")|(clust.method=="complete")|
      (clust.method=="average")|(clust.method=="mcquitty")|
      (clust.method=="median")|(clust.method=="centroid")
    ){
      d <- as.dist(distance)
      fit <- hclust(d, method=clust.method) 
      classification <- cutree(fit, k=k)
      g <- graph.adjacency(distance, weighted=TRUE, mode="undirected")
      g <- delete.edges(g, E(g))
      for(i in unique(classification)){
        if(!is.na(i)){
          tmp <- which(classification==i)
          tmp <- rbind(tmp, classification[tmp])
          edge_list <- array(unlist(tmp))
          edge_weights <- c()
          for(j in 1:ncol(tmp)){
            edge_weights[j] <- distance[tmp[1,j],tmp[2,j]]
          }
          g <- add.edges(g, edge_list, attr=list(weight=edge_weights))
        }
      }
    }
    if(clust.method=="k-medoids"){
      classification <-pam(as.dist(distance), k=k,diss=TRUE,cluster.only = TRUE)
      g <- graph.adjacency(distance, weighted=TRUE, mode="undirected")
      g <- delete.edges(g, E(g))
      for(i in unique(classification)){
        if(!is.na(i)){
          tmp <- which(classification==i)
          tmp <- rbind(tmp, classification[tmp])
          edge_list <- array(unlist(tmp))
          edge_weights <- c()
          for(j in 1:ncol(tmp)){
            edge_weights[j] <- distance[tmp[1,j],tmp[2,j]]
          }
          g <- add.edges(g, edge_list, attr=list(weight=edge_weights))
        }
      }
    }
    SincellObject[["cellsClustering"]]=g
    SincellObject[["clust.method"]]=clust.method
    SincellObject[["mutual"]]=mutual
    SincellObject[["k"]]=k
    SincellObject[["max.distance"]]=max.distance
    SincellObject[["shortest.rank.percent"]]=shortest.rank.percent
    return(SincellObject)
  }
}
#===================================================
sc_GraphBuilderObj <- function(SincellObject, graph.algorithm="MST", graph.using.cells.clustering=FALSE,k=3){
  if(is.list(SincellObject)==FALSE){
    errormesseage<-paste(deparse(substitute(SincellObject))," is not a list")
    stop(errormesseage)
  }else if(is.null(SincellObject[["expressionmatrix"]])==TRUE){
    errormesseage<-paste(deparse(substitute(SincellObject)),"does not contain a valid expressionmatrix member. Please initialize using function sc_InitializingSincellObject()")
    stop(errormesseage)
  }else if(is.null(SincellObject[["cell2celldist"]])==TRUE){
    errormesseage<-paste(deparse(substitute(SincellObject)),"does not contain a valid cell2celldist member. Please assess a cell-2-cell distance matrix using either function sc_distanceObj() or function sc_DimensionalityReductionObj()")
    stop(errormesseage)
  }else if(is.logical(graph.using.cells.clustering)&(graph.using.cells.clustering==TRUE)&((graph.algorithm=="MST")|(graph.algorithm=="SST"))&(is.null(SincellObject[["cellsClustering"]]))){
    errormesseage<-paste("graph.using.cells.clustering=TRUE is indicated, however",deparse(substitute(SincellObject)),"does not contain a valid cellsClustering member. Please assess a cells clustering using either function sc_clusterObj() or switch to graph.using.cells.clustering=FALSE")
    stop(errormesseage)
  }else if(!graph.algorithm %in% c("MST", "SST", "IMC")){
    errormesseage<-paste(graph.algorithm,"is not a valid option for parameter graph.algorithm . Please read documentation of sc_GraphBuilderObj() in Sincell manual")
    stop(errormesseage)
  }else{ 
    SincellObject<-SincellObject[c("expressionmatrix","cellsLowDimensionalSpace","cell2celldist","method","dim","MDS.distance","bins","cellsClustering","clust.method","mutual","k","max.distance","shortest.rank.percent")]
    cell2celldist<-SincellObject[["cell2celldist"]]
    if(is.logical(graph.using.cells.clustering)&(graph.using.cells.clustering==TRUE)){cellsClustering=SincellObject[["cellsClustering"]]}
    
    if(graph.algorithm=="MST"){
      g <- minimum.spanning.tree(graph.adjacency(cell2celldist, weighted=TRUE, mode="undirected"))
      if(!is.logical(graph.using.cells.clustering)){
        g <- get.adjacency(g,type="both",attr="weight",sparse=FALSE)
        gt <- get.adjacency(cellsClustering,type="both",attr="weight",sparse=FALSE)
        for(i in which(gt > g)){
          g[i] <- gt[i]
        }
        g <- graph.adjacency(g, weighted=TRUE, mode="undirected")
      }
    }
    if(graph.algorithm=="SST"){
      if(is.logical(graph.using.cells.clustering) && graph.using.cells.clustering==FALSE){
        g <- graph.adjacency(cell2celldist ,weighted=TRUE,mode="undirected")
        g <- delete.edges(g, E(g))
      } else{
        g <- cellsClustering
      }
      g.clusters <- clusters(g)
      while(g.clusters$no > 1){
        tmp <- sstalgorithm(membership=g.clusters$membership, num_cells=length(g.clusters$membership), distance=cell2celldist)
        dist.min <- tmp[1]
        dist.min.coor.x <- tmp[2]
        dist.min.coor.y <- tmp[3]
        g <- add.edges(g, c(dist.min.coor.x,dist.min.coor.y), attr=list(weight=dist.min))
        g.clusters <- clusters(g)
      }
    }
    if(graph.algorithm=="IMC"){
      SincellObject[["clust.method"]]="knn"
      SincellObject[["mutual"]]=TRUE
      SincellObject[["k"]]=k
      m <- max(cell2celldist)
      g <- knnalgorithm(cell2celldist, mutual=TRUE, k=k)
      SincellObject[["cellsClustering"]]<-g
      graph.using.cells.clustering=TRUE
      g.clusters <- clusters(g)
      while(g.clusters$no > 1){
        for(i in which(g.clusters$csize > 1)){
          tmp <- which(g.clusters$membership==i)
          tmp <- t(combn(tmp,2))
          cell2celldist[tmp] <- m+1
          cell2celldist[cbind(tmp[,2], tmp[,1])] <- m+1
        }
        gt <- knnalgorithm(cell2celldist, mutual=TRUE, k=k)
        g <- get.adjacency(g,type="both",attr="weight",sparse=FALSE)
        gt <- get.adjacency(gt,type="both",attr="weight",sparse=FALSE)
        gt[gt>m] <- 0
        g <- graph.adjacency(g+gt, weighted=TRUE, mode="undirected")
        g.clusters <- clusters(g)    
      }
    }
    SincellObject[["cellstateHierarchy"]]=g
    SincellObject[["graph.algorithm"]]=graph.algorithm
    SincellObject[["graph.using.cells.clustering"]]=graph.using.cells.clustering
    return(SincellObject)
  }
}
#===================================================
## Statistical Support By Gene Subsampling
sc_StatisticalSupportByGeneSubsampling <- function(SincellObject, num_it=100, num_genes=as.integer(nrow(SincellObject[["expressionmatrix"]])*0.5), cores=ifelse(detectCores()>=4, 4, detectCores())){
  if(is.list(SincellObject)==FALSE){
    errormesseage<-paste(deparse(substitute(SincellObject))," is not a list")
    stop(errormesseage)
  }else if(is.null(SincellObject[["expressionmatrix"]])==TRUE){
    errormesseage<-paste(deparse(substitute(SincellObject)),"does not contain a valid expressionmatrix member. Please initialize using function sc_InitializingSincellObject()")
    stop(errormesseage)
  }else if(is.null(SincellObject[["cell2celldist"]])==TRUE){
    errormesseage<-paste(deparse(substitute(SincellObject)),"does not contain a valid cell2celldist member. Please assess a cell-2-cell distance matrix using either function sc_distanceObj() or function sc_DimensionalityReductionObj()")
    stop(errormesseage)
  }else if(is.null(SincellObject[["cellstateHierarchy"]])==TRUE){
    errormesseage<-paste(deparse(substitute(SincellObject)),"does not contain a valid cellstateHierarchy member. Please assess a cellstateHierarchy using function sc_GraphBuilderObj()")
    stop(errormesseage)
  }else if(num_genes>nrow(SincellObject[["expressionmatrix"]])){
    errormesseage<-paste("The number of genes indicated for resampling is higher than the number of genes (rows) in member expressionmatrix of",deparse(substitute(SincellObject)),". Please set num_genes to a lower number")
    stop(errormesseage)
  }else{ 
    if(!(.Platform$OS.type=="unix")){cores=1;}
    if(is.null(SincellObject[["StatisticalSupportbyGeneSubsampling"]])==FALSE){
      SincellObject[["StatisticalSupportbyGeneSubsampling"]]<-NULL
    }
    ## Trees from subsets
    f_mclapply <- function(x){
      gene_subset <- sample(1:nrow(SincellObject[["expressionmatrix"]]), num_genes)
      tmp <- SincellObject[["expressionmatrix"]][gene_subset,]
      mySincellObjecttmp <- sc_InitializingSincellObject(tmp)
      if(
        (SincellObject[["method"]]=="euclidean") | (SincellObject[["method"]]=="L1") | (SincellObject[["method"]]=="pearson") | (SincellObject[["method"]]=="MI")
      ){
        mySincellObjecttmp<- sc_distanceObj(mySincellObjecttmp, method=SincellObject[["method"]], bins=SincellObject[["bins"]])
      }
      if(
        (SincellObject[["method"]]=="PCA") | (SincellObject[["method"]]=="ICA") | (SincellObject[["method"]]=="tSNE") | (SincellObject[["method"]]=="classical-MDS") | (SincellObject[["method"]]=="nonmetric-MDS")
      ){
        mySincellObjecttmp<-  sc_DimensionalityReductionObj(mySincellObjecttmp, method=SincellObject[["method"]], dim=SincellObject[["dim"]], MDS.distance=SincellObject[["MDS.distance"]], bins=SincellObject[["bins"]])
      }
      if(((SincellObject[["graph.algorithm"]]=="MST")|(SincellObject[["graph.algorithm"]]=="SST"))&(SincellObject[["graph.using.cells.clustering"]]==TRUE)){
        mySincellObjecttmp<- sc_clusterObj (mySincellObjecttmp, clust.method=SincellObject[["clust.method"]], mutual=SincellObject[["mutual"]], k=SincellObject[["k"]], max.distance=SincellObject[["max.distance"]], shortest.rank.percent=SincellObject[["shortest.rank.percent"]])
        mySincellObjecttmp<- sc_GraphBuilderObj(mySincellObjecttmp, graph.algorithm=SincellObject[["graph.algorithm"]], graph.using.cells.clustering=SincellObject[["graph.using.cells.clustering"]],k=SincellObject[["k"]])
      }
      if(((SincellObject[["graph.algorithm"]]=="MST")|(SincellObject[["graph.algorithm"]]=="SST"))&(SincellObject[["graph.using.cells.clustering"]]==FALSE)){
        mySincellObjecttmp<- sc_GraphBuilderObj(mySincellObjecttmp, graph.algorithm=SincellObject[["graph.algorithm"]], graph.using.cells.clustering=SincellObject[["graph.using.cells.clustering"]],k=SincellObject[["k"]])
      }
      if(SincellObject[["graph.algorithm"]]=="IMC"){
        mySincellObjecttmp<- sc_GraphBuilderObj(mySincellObjecttmp, graph.algorithm=SincellObject[["graph.algorithm"]], graph.using.cells.clustering=SincellObject[["graph.using.cells.clustering"]],k=SincellObject[["k"]])
      }
      node.dist <- shortest.paths(mySincellObjecttmp[["cellstateHierarchy"]])
      return(f_distance2vector(node.dist))
    }
    tree.list <-  mclapply(1:num_it, f_mclapply, mc.cores=cores, mc.preschedule = TRUE)
    #if(sum(unlist(lapply(tree.list, is.null)))){
      #while(sum(unlist(lapply(tree.list, is.null)))){
      while(!all(unlist(lapply(tree.list, is.vector))) || !all(unlist(lapply(tree.list, is.numeric))) || any(unlist(lapply(tree.list, is.na))) || any(unlist(lapply(tree.list, is.null))) || length(tree.list)!=num_it || !all(unlist(lapply(tree.list, length))==nrow(SincellObject[["cell2celldist"]])*(nrow(SincellObject[["cell2celldist"]])-1)/2) ){
        # print("WARNING: memory error. Re-executing.")
        tree.list <-  mclapply(1:num_it, f_mclapply, mc.cores=cores, mc.preschedule = TRUE)
      }
    #}
    node.dist <- shortest.paths(SincellObject[["cellstateHierarchy"]])
    # SincellObject[["tree.list"]]<-tree.list
    correlation<-cor(f_distance2vector(node.dist),matrix(unlist(tree.list),ncol=num_it,nrow=(ncol(SincellObject[["expressionmatrix"]])*(ncol(SincellObject[["expressionmatrix"]])-1)/2)), method="spearman")
    SincellObject[["StatisticalSupportbyGeneSubsampling"]]<-as.vector(correlation)
    messeage<-paste("The summary of the distribution of Spearman rank correlations\nbetween the original hierarchy and the hierarchies obtained\nfrom",num_it,"resamplings of",num_genes,"genes in the initial expression matrix is:\n")
    cat(messeage)
    print(summary(as.vector(correlation)))
    # plot(density(correlation),main="Distribution of correlations",xlim=c(0,1),col="red",ylab="Density",xlab="Spearman rank correlation")
    return(SincellObject)
  }
}
# f_GraphBuilderObjfromTemplate<- function(expressionmatrix,ReferenceSincellObject=SincellObject){  
#   mySincellObjecttmp <- sc_InitializingSincellObject(expressionmatrix)
#   if(
#     (ReferenceSincellObject[["method"]]=="euclidean") | (ReferenceSincellObject[["method"]]=="L1") | (ReferenceSincellObject[["method"]]=="pearson") | (ReferenceSincellObject=="MI")
#   ){
#     mySincellObjecttmp<- sc_distanceObj(mySincellObjecttmp, method=ReferenceSincellObject[["method"]], bins=ReferenceSincellObject[["bins"]])
#   }
#   if(
#     (ReferenceSincellObject[["method"]]=="PCA") | (ReferenceSincellObject[["method"]]=="ICA") | (ReferenceSincellObject[["method"]]=="tSNE") | (ReferenceSincellObject[["method"]]=="classical-MDS") | (ReferenceSincellObject[["method"]]=="nonmetric-MDS")
#   ){
#     mySincellObjecttmp<-  sc_DimensionalityReductionObj(mySincellObjecttmp, method=ReferenceSincellObject[["method"]], dim=ReferenceSincellObject[["dim"]], MDS.distance=ReferenceSincellObject[["MDS.distance"]], bins=ReferenceSincellObject[["bins"]])
#   }
#   if(((ReferenceSincellObject[["graph.algorithm"]]=="MST")|(ReferenceSincellObject[["graph.algorithm"]]=="SST"))&(ReferenceSincellObject[["graph.using.cells.clustering"]]==TRUE)){
#     mySincellObjecttmp<- sc_clusterObj (mySincellObjecttmp, clust.method=ReferenceSincellObject[["clust.method"]], mutual=ReferenceSincellObject[["mutual"]], k=ReferenceSincellObject[["k"]], max.distance=ReferenceSincellObject[["max.distance"]], shortest.rank.percent=ReferenceSincellObject[["shortest.rank.percent"]])
#     mySincellObjecttmp<- sc_GraphBuilderObj(mySincellObjecttmp, graph.algorithm=ReferenceSincellObject[["graph.algorithm"]], graph.using.cells.clustering=ReferenceSincellObject[["graph.using.cells.clustering"]],k=ReferenceSincellObject[["k"]])
#   }
#   if(((ReferenceSincellObject[["graph.algorithm"]]=="MST")|(ReferenceSincellObject[["graph.algorithm"]]=="SST"))&(ReferenceSincellObject[["graph.using.cells.clustering"]]==FALSE)){
#     mySincellObjecttmp<- sc_GraphBuilderObj(mySincellObjecttmp, graph.algorithm=ReferenceSincellObject[["graph.algorithm"]], graph.using.cells.clustering=ReferenceSincellObject[["graph.using.cells.clustering"]],k=ReferenceSincellObject[["k"]])
#   }
#   if(ReferenceSincellObject[["graph.algorithm"]]=="IMC"){
#     mySincellObjecttmp<- sc_GraphBuilderObj(mySincellObjecttmp, graph.algorithm=ReferenceSincellObject[["graph.algorithm"]], graph.using.cells.clustering=ReferenceSincellObject[["graph.using.cells.clustering"]],k=ReferenceSincellObject[["k"]])
#   }
#   return(mySincellObjecttmp)
# }
#===================================================
## Convert distance matrix to vector of dim*(dim-1)/2 elements
f_distance2vector <- function(distance){
  return(as.vector(as.matrix(distance)[upper.tri(as.matrix(distance), diag = FALSE)]))
}
#===================================================
## Generation of In Silico Cells Replicates
sc_InSilicoCellsReplicatesObj <- function(SincellObject, method="variance.deciles",dispersion.statistic = NULL, multiplier=100, no_expr=0.5, LogTransformedData = T, baseLogTransformation=exp(1),pseudocounts.added.before.log.transformation=1, cores=ifelse(detectCores()>=4, 4, detectCores())){
  positive=TRUE
  if(is.list(SincellObject)==FALSE){
    errormesseage<-paste(deparse(substitute(SincellObject))," is not a list")
    stop(errormesseage)
  }else if(is.null(SincellObject[["expressionmatrix"]])==TRUE){
    errormesseage<-paste(deparse(substitute(SincellObject)),"does not contain a valid expressionmatrix member. Please initialize using function sc_InitializingSincellObject()")
    stop(errormesseage)
  }else if(method!="variance.deciles" && method!="cv2.deciles" && method!="lognormal-3parameters" && method!="negative.binomial"){
    errormesseage<-paste(deparse(substitute(method)),"is not a valid option for parameter method . Please read documentation of sc_StatisticalSupportByReplacementWithInSilicoCellsReplicates() in Sincell manual")
    stop(errormesseage)
  }else if(!is.null(dispersion.statistic) && dispersion.statistic!="cv2.fitted.to.data" && !is.numeric(dispersion.statistic)){
    errormesseage<-paste(deparse(substitute(dispersion.statistic)),"is not a valid option for parameter dispersion.statistic . Please read documentation of sc_InSilicoCellsReplicatesObj() in Sincell manual")
    stop(errormesseage)
  }else if(method=="lognormal-3parameters" && LogTransformedData==F){
    errormesseage<-paste(deparse(substitute(method)),"is designed for Log-transformed expression data, however parameter LogTransformedData was set to FALSE. Please read documentation of sc_InSilicoCellsReplicatesObj() in Sincell manual")
    stop(errormesseage)
  }else{
    undo.unlog <- F
    if(!(.Platform$OS.type=="unix")){cores=1;}
    if(is.null(SincellObject[["InSilicoCellsReplicates"]])==FALSE){
      SincellObject[["InSilicoCellsReplicates"]]<-NULL
      SincellObject[["multiplier"]]<-NULL
      SincellObject[["methodInSilicoCellsReplicates"]]<-NULL
    }
    BaseData<-SincellObject[["expressionmatrix"]]
    nr <- nrow(BaseData)
    nc <- ncol(BaseData)
    
    if(method=="variance.deciles" | method=="cv2.deciles"){
      var_genes <- apply(BaseData,1,var)
      mean_genes <- apply(BaseData,1,mean)
      cv2_genes <- var_genes/mean_genes/mean_genes
      index <- 1:nr
      values_df <- data.frame(index, mean_genes, var_genes, cv2_genes)
      values_df <- values_df[order(values_df$mean_genes),]
      
      deciles <- c(1,length(var_genes)/10, 2*length(var_genes)/10, 3*length(var_genes)/10,
                   4*length(var_genes)/10, 5*length(var_genes)/10, 6*length(var_genes)/10,
                   7*length(var_genes)/10, 8*length(var_genes)/10, 9*length(var_genes)/10,
                   length(var_genes))
      if(method=="variance.deciles"){
        var_genes_sorted <- as.array(values_df$var_genes)
        indexes <- (1:nr)[order(as.numeric(values_df$index))]
        f_mclapply <- function(x){
          newData <- matrix(0,nrow=nr,ncol=nc)
          newData <- pseudoreplicatesbynoise(originaldata=as.matrix(BaseData), rows=nr, colums=nc, deciles=as.integer(deciles), lengthdeciles=length(deciles), coorsorted=indexes, vargenessorted=var_genes_sorted, positive=as.numeric(positive), seed=as.integer(runif(1,-2141716335,2146848561)))
          return(newData)
        }
      }
      if(method=="cv2.deciles"){
        var_genes_sorted <- as.array(values_df$cv2_genes)
        indexes <- (1:nr)[order(as.numeric(values_df$index))]
        f_mclapply <- function(x){
          newData <- matrix(0,nrow=nr,ncol=nc)
          newData <- pseudoreplicatesbynoise_cv2(originaldata=as.matrix(BaseData), rows=nr, colums=nc, deciles=as.integer(deciles), lengthdeciles=length(deciles), coorsorted=indexes, vargenessorted=var_genes_sorted, means=mean_genes, positive=as.numeric(positive), seed=as.integer(runif(1,-2141716335,2146848561)))
          return(newData)
        }
      }
    }
    if(method=="lognormal-3parameters"){
      alpha <- apply(BaseData>no_expr,1,sum)/nc
      expressedvar <- function(x, no.expr=no_expr){
        expressedcells <- x[x>no.expr]
        return(var(expressedcells))
      }
      var.genes <- apply(BaseData,1,expressedvar)
      Data <- BaseData
      Data[Data<=no_expr] <- 0
      mean.genes <- apply(Data,1,sum)/apply(BaseData>no_expr,1,sum)
      f_mclapply <- function(x){
        newData <- matrix(0,nrow=nr,ncol=nc)
        newData <- pseudoreplicatesbymodel(rows=nr, colums=nc, alpha=alpha, vargenes=var.genes, meangenes=mean.genes, positive=as.numeric(positive), f=rnorm, seed=as.integer(runif(1,-2141716335,2146848561)))
        return(newData)
      }
    }
    if(method=="negative.binomial"){
      if(LogTransformedData==T){
        BaseData <- baseLogTransformation**BaseData-pseudocounts.added.before.log.transformation
        undo.unlog <- T
      }
      var_genes <- apply(BaseData,1,var)
      mean_genes <- apply(BaseData,1,mean)
      cv2_genes <- var_genes/mean_genes/mean_genes
      
      if(is.null(dispersion.statistic)){
        r <- mean_genes*mean_genes/(var_genes-mean_genes)
      }else if(dispersion.statistic=="cv2.fitted.to.data"){
        minMeanForFit <- unname( quantile( mean_genes[ which( cv2_genes > .3 ) ], .95 ) )
        useForFit <- mean_genes >= minMeanForFit
        fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/mean_genes[useForFit] ),cv2_genes[useForFit] )
        a0 <- unname( fit$coefficients["a0"] )
        a1 <- unname( fit$coefficients["a1tilde"])
        cv2.estimated <- a1/mean_genes + a0
        var.genes <- cv2.estimated*mean_genes*mean_genes
        r <- mean_genes*mean_genes/(var.genes-mean_genes)
      }else if(is.numeric(dispersion.statistic)){
        var.genes <- cv2.estimated*mean_genes*mean_genes
        r <- mean_genes*mean_genes/(var.genes-mean_genes)
      }
      
      if(as.numeric(summary(r<=0)[3])>0){
        messeage<-paste("WARNING: A total of ",as.numeric(summary(r<0)[3]),"genes have dispersion value r<=0 and will not be perturbed with a random value drawn from the negative binomial distribution\n")
        cat(messeage)
      }
      f_mclapply <- function(x){
        newData <- matrix(0,nrow=nr,ncol=nc)
        for(i in 1:nr){
          if(r[i]>0){
            newData[i,] <- rnbinom(nc,mu=mean_genes[i],size=r[i])
          }else{
            newData[i,] <- BaseData[i,]
          }
        }
        return(newData)
      }
    }
    replicates <- mclapply(1:multiplier, f_mclapply, mc.cores=cores, mc.preschedule = TRUE)
    #if(sum(unlist(lapply(replicates, is.null)))>0){
    #while(sum(unlist(lapply(tree.list, is.null)))){
    while(!all(unlist(lapply(replicates, is.matrix))) || !all(unlist(lapply(replicates, is.numeric))) || any(unlist(lapply(replicates, is.na))) || any(unlist(lapply(replicates, is.null))) || length(replicates)!=multiplier || !all(unlist(lapply(replicates, ncol))==ncol(SincellObject[["expressionmatrix"]])) || !all(unlist(lapply(replicates, nrow))==nrow(SincellObject[["expressionmatrix"]])) ){
      # print("WARNING: memory error. Re-executing.")
      replicates <- mclapply(1:multiplier, f_mclapply, mc.cores=cores, mc.preschedule = TRUE)
    }
    #}
    if(undo.unlog==F){
      SincellObject[["InSilicoCellsReplicates"]]<-cbind(BaseData, matrix(unlist(replicates),nrow=nr))
    }else{
      SincellObject[["InSilicoCellsReplicates"]]<-cbind(SincellObject[["expressionmatrix"]],log(matrix(unlist(replicates),nrow=nr)+pseudocounts.added.before.log.transformation),base=baseLogTransformation)
    }
    TotalNumNA<-(sum(unlist(lapply(replicates, is.na))))
    if(TotalNumNA>0){
      messeage<-paste("WARNING: Total NAs in Replicates",TotalNumNA,"\n")
      cat(messeage)
    }
    SincellObject[["multiplier"]]<-multiplier
    SincellObject[["methodInSilicoCellsReplicates"]]<-method
    return(SincellObject)
  }
}
#=============================================
## Statistical Support By Replacement With InSilico Cell Replicates
sc_StatisticalSupportByReplacementWithInSilicoCellsReplicates <- function(SincellObject, method="own", num_it=100, fraction.cells.to.replace=0.15, cores=ifelse(detectCores()>=4, 4, detectCores())){
  if(is.list(SincellObject)==FALSE){
    errormesseage<-paste(deparse(substitute(SincellObject))," is not a list")
    stop(errormesseage)
  }else if(is.null(SincellObject[["expressionmatrix"]])==TRUE){
    errormesseage<-paste(deparse(substitute(SincellObject)),"does not contain a valid expressionmatrix member. Please initialize using function sc_InitializingSincellObject()")
    stop(errormesseage)
  }else if(is.null(SincellObject[["cell2celldist"]])==TRUE){
    errormesseage<-paste(deparse(substitute(SincellObject)),"does not contain a valid cell2celldist member. Please assess a cell-2-cell distance matrix using either function sc_distanceObj() or function sc_DimensionalityReductionObj()")
    stop(errormesseage)
  }else if(is.null(SincellObject[["cellstateHierarchy"]])==TRUE){
    errormesseage<-paste(deparse(substitute(SincellObject)),"does not contain a valid cellstateHierarchy member. Please assess a cellstateHierarchy using function sc_GraphBuilderObj()")
    stop(errormesseage)
  }else if(is.null(SincellObject[["InSilicoCellsReplicates"]])==TRUE){
    errormesseage<-paste(deparse(substitute(SincellObject)),"does not contain in silico generated cells replicates (list member 'InSilicoCellsReplicates' is NULL). Please generate cell replicates using function sc_InSilicoCellsReplicatesObj()")
    stop(errormesseage)
  }else if(is.character(method)&!(method=="own")&!(method=="all")){
    errormesseage<-paste(deparse(substitute(method)),"is not a valid option for parameter method . Please read documentation of sc_StatisticalSupportByReplacementWithInSilicoCellsReplicates() in Sincell manual")
    stop(errormesseage)
  }else{ 
    if(!(.Platform$OS.type=="unix")){cores=1;}
    if(is.null(SincellObject[["StatisticalSupportByReplacementWithInSilicoCellReplicates"]])==FALSE){
      SincellObject[["StatisticalSupportByReplacementWithInSilicoCellReplicates"]]<-NULL
    }
    if(is.numeric(method)){
      max.order <- max(shortest.paths(SincellObject[["cellstateHierarchy"]], weights=NA))
      if(method > max.order){
        errormesseage<-paste("WARNING: The highest neighbor order in current member 'cellstateHierarchy' of",deparse(substitute(SincellObject)),"is", max.order, ". Please set parameter 'method' to a number equal or lower than", max.order,". For more information please read documentation of sc_StatisticalSupportByReplacementWithInSilicoCellsReplicates() in Sincell manual.\n")
        stop(errormesseage)
      } else {
         if(method == 0){method <- "own"}
         if(method == max.order){method <- "all"}
      }
    }
    multiplier<-SincellObject[["multiplier"]]
    num.original.cells<-ncol(SincellObject[["expressionmatrix"]])
    Pseudoreplicates<- SincellObject[["InSilicoCellsReplicates"]]
    tree.original <- SincellObject[["cellstateHierarchy"]]
    f_mclapply <- function(x){
      ## Changes for its own pseudoreplicate
      if(method=="own"){
        tmp <- Pseudoreplicates[,1:num.original.cells]
        changes <- sample(1:num.original.cells,fraction.cells.to.replace*num.original.cells)
        for(i in changes){
          tmp[,i] <- Pseudoreplicates[,i+sample(1:multiplier,1)*num.original.cells]
        }
      }
      ## Changes for any pseudoreplicate
      if(method=="all"){
        tmp <- Pseudoreplicates[,1:num.original.cells]
        changes <- sample(1:num.original.cells,fraction.cells.to.replace*num.original.cells)
        for(i in changes){
          tmp[,i] <- Pseudoreplicates[,sample((num.original.cells+1):ncol(Pseudoreplicates),1)]
        }
      }
      ## Changes for its order n neighbors' pseudoreplicates
      if(is.numeric(method)){
        tmp <- Pseudoreplicates[,1:num.original.cells]
        changes <- sample(1:num.original.cells,fraction.cells.to.replace*num.original.cells)
        for(i in changes){
          tmp[,i] <- Pseudoreplicates[,sample(igraph::neighborhood(tree.original,method,i)[[1]],1)+sample(1:multiplier,1)*num.original.cells]
        }
      }
      mySincellObjecttmp <- sc_InitializingSincellObject(tmp)
      if(
        (SincellObject[["method"]]=="euclidean") | (SincellObject[["method"]]=="L1") | (SincellObject[["method"]]=="pearson") | (SincellObject[["method"]]=="MI")
      ){
        mySincellObjecttmp<- sc_distanceObj(mySincellObjecttmp, method=SincellObject[["method"]], bins=SincellObject[["bins"]])
      }
      if(
        (SincellObject[["method"]]=="PCA") | (SincellObject[["method"]]=="ICA") | (SincellObject[["method"]]=="tSNE") | (SincellObject[["method"]]=="classical-MDS") | (SincellObject[["method"]]=="nonmetric-MDS")
      ){
        mySincellObjecttmp<-  sc_DimensionalityReductionObj(mySincellObjecttmp, method=SincellObject[["method"]], dim=SincellObject[["dim"]], MDS.distance=SincellObject[["MDS.distance"]], bins=SincellObject[["bins"]])
      }
      if(((SincellObject[["graph.algorithm"]]=="MST")|(SincellObject[["graph.algorithm"]]=="SST"))&(SincellObject[["graph.using.cells.clustering"]]==TRUE)){
        mySincellObjecttmp<- sc_clusterObj (mySincellObjecttmp, clust.method=SincellObject[["clust.method"]], mutual=SincellObject[["mutual"]], k=SincellObject[["k"]], max.distance=SincellObject[["max.distance"]], shortest.rank.percent=SincellObject[["shortest.rank.percent"]])
        mySincellObjecttmp<- sc_GraphBuilderObj(mySincellObjecttmp, graph.algorithm=SincellObject[["graph.algorithm"]], graph.using.cells.clustering=SincellObject[["graph.using.cells.clustering"]],k=SincellObject[["k"]])
      }
      if(((SincellObject[["graph.algorithm"]]=="MST")|(SincellObject[["graph.algorithm"]]=="SST"))&(SincellObject[["graph.using.cells.clustering"]]==FALSE)){
        mySincellObjecttmp<- sc_GraphBuilderObj(mySincellObjecttmp, graph.algorithm=SincellObject[["graph.algorithm"]], graph.using.cells.clustering=SincellObject[["graph.using.cells.clustering"]],k=SincellObject[["k"]])
      }
      if(SincellObject[["graph.algorithm"]]=="IMC"){
        mySincellObjecttmp<- sc_GraphBuilderObj(mySincellObjecttmp, graph.algorithm=SincellObject[["graph.algorithm"]], graph.using.cells.clustering=SincellObject[["graph.using.cells.clustering"]],k=SincellObject[["k"]])
      }
      
      node.dist.matrix <- matrix(0,0,(num.original.cells*(num.original.cells-1)/2))
      node.dist <- shortest.paths(SincellObject[["cellstateHierarchy"]])
      node.dist.matrix <- rbind(node.dist.matrix,f_distance2vector(node.dist))    
      
      node.dist <- shortest.paths(mySincellObjecttmp[["cellstateHierarchy"]])
      return(f_distance2vector(node.dist))
    }
    tree.list <-  mclapply(1:num_it, f_mclapply, mc.cores=cores, mc.preschedule = TRUE)
    #if(sum(unlist(lapply(tree.list, is.null)))){
      #while(sum(unlist(lapply(tree.list, is.null)))){
       while(!all(unlist(lapply(tree.list, is.vector))) || !all(unlist(lapply(tree.list, is.numeric))) || any(unlist(lapply(tree.list, is.na))) || any(unlist(lapply(tree.list, is.null))) || length(tree.list)!=num_it || !all(unlist(lapply(tree.list, length))==nrow(SincellObject[["cell2celldist"]])*(nrow(SincellObject[["cell2celldist"]])-1)/2) ){
        # print("WARNING: memory error. Re-executing.")
        tree.list <-  mclapply(1:num_it, f_mclapply, mc.cores=cores, mc.preschedule = TRUE)
      }
    #}
    node.dist <- shortest.paths(SincellObject[["cellstateHierarchy"]])
    correlation<-cor(f_distance2vector(node.dist),matrix(unlist(tree.list),ncol=num_it,nrow=(ncol(SincellObject[["expressionmatrix"]])*(ncol(SincellObject[["expressionmatrix"]])-1)/2)), method="spearman")
    SincellObject[["StatisticalSupportByReplacementWithInSilicoCellReplicates"]]<-as.vector(correlation)
    messeage<-paste("The summary of the distribution of spearman rank correlations\n between the original hierarchy and the hierarchies obtained from\n substitution of in silico generated replicates of order",method,"\nin the initial expression matrix is:\n")
    cat(messeage)
    print(summary(as.vector(correlation)))
    # plot(density(correlation),main="Distribution of correlations",xlim=c(0,1),col="red",ylab="Density",xlab="Spearman rank correlation")
    return(SincellObject)
  }
}
# mySincellObject <- sc_InitializingSincellObject(ExpressionMatrix_StdLog2)
# mySincellObject <- sc_DimensionalityReductionObj(mySincellObject, method="ICA",dim=2)
# mySincellObject <- sc_clusterObj (mySincellObject, clust.method="max.distance", max.distance=0.5)
# mySincellObject <- sc_clusterObj (mySincellObject, clust.method="knn", mutual=TRUE)
# mySincellObject<- sc_GraphBuilderObj(mySincellObject, graph.algorithm="SST", graph.using.cells.clustering=TRUE)
# mySincellObject<- sc_GraphBuilderObj(mySincellObject, graph.algorithm="MST", graph.using.cells.clustering=FALSE)
# system.time(mySincellObjectTest<-sc_StatisticalSupportByGeneSubsampling(mySincellObject,num_it=100))
# mySincellObject<- sc_InSilicoCellsReplicatesObj(mySincellObject, method="noise", multiplier=100, no_expr=0.5)
# system.time(mySincellObjectTest<-sc_StatisticalSupportByReplacementWithInSilicoCellsReplicates(mySincellObject, method="own", multiplier=100, num_it=100, fraction.cells.to.replace=0.15))
# system.time(mySincellObjectTest2<-sc_StatisticalSupportByReplacementWithInSilicoCellsReplicates(mySincellObject, method=2, multiplier=100, num_it=100, fraction.cells.to.replace=0.15))
# system.time(mySincellObjectTest2<-sc_StatisticalSupportByReplacementWithInSilicoCellsReplicates(mySincellObject, method=16, multiplier=100, num_it=100, fraction.cells.to.replace=0.15))
# system.time(mySincellObjectTest2<-sc_StatisticalSupportByReplacementWithInSilicoCellsReplicates(mySincellObject, method="all", multiplier=100, num_it=100, fraction.cells.to.replace=0.15))
#===================================================
# Statistical Test for Association Of Cells Hierarchy With A Gene Set
sc_AssociationOfCellsHierarchyWithAGeneSet <- function(SincellObject,GeneSet, minimum.geneset.size=50,p.value.assessment=TRUE,spearman.rank.threshold=0.5,num_it=1000, cores=ifelse(detectCores()>=4, 4, detectCores())){
  if(is.list(SincellObject)==FALSE){
    errormesseage<-paste(deparse(substitute(SincellObject))," is not a list")
    stop(errormesseage)
  }else if(is.null(SincellObject[["expressionmatrix"]])==TRUE){
    errormesseage<-paste(deparse(substitute(SincellObject)),"does not contain a valid expressionmatrix member. Please initialize using function sc_InitializingSincellObject()")
    stop(errormesseage)
  }else if(is.null(SincellObject[["cell2celldist"]])==TRUE){
    errormesseage<-paste(deparse(substitute(SincellObject)),"does not contain a valid cell2celldist member. Please assess a cell-2-cell distance matrix using either function sc_distanceObj() or function sc_DimensionalityReductionObj()")
    stop(errormesseage)
  }else if(is.null(SincellObject[["cellstateHierarchy"]])==TRUE){
    errormesseage<-paste(deparse(substitute(SincellObject)),"does not contain a valid cellstateHierarchy member. Please assess a cellstateHierarchy using function sc_GraphBuilderObj()")
    stop(errormesseage)
  }else if((is.vector(GeneSet)==FALSE)|(is.character(GeneSet)==FALSE)){
    errormesseage<-paste(deparse(substitute(GeneSet))," is not a vector of characters")
    stop(errormesseage)
  }else if(sum(rownames(SincellObject[["expressionmatrix"]]) %in% GeneSet)==0){
    errormesseage<-paste("The expression matrix in",deparse(substitute(SincellObject)), "does not contain any gene indicated in",deparse(substitute(GeneSet)),". Please check that gene names in the expression matrix and gene names in your gene set are using the same type of identifier (e.g. HGNC, Ensembl, etc.)")
    stop(errormesseage)
  }else if(sum(rownames(SincellObject[["expressionmatrix"]]) %in% GeneSet)<minimum.geneset.size){
    EffectiveGeneSetSizetmp<-sum(rownames(SincellObject[["expressionmatrix"]]) %in% GeneSet)
    errormesseage<-paste("The expression matrix in",deparse(substitute(SincellObject)), "only contains",EffectiveGeneSetSizetmp,"genes, which is lower than the value indicated for parameter minimum.geneset.size=",minimum.geneset.size,". Please try with a larger gene set")
    stop(errormesseage)
  }else{ 
    if(!(.Platform$OS.type=="unix")){cores=1;}
    if(is.null(SincellObject[["AssociationOfCellsHierarchyWithAGeneSet"]])==FALSE){
      SincellObject[["AssociationOfCellsHierarchyWithAGeneSet"]]<-NULL
      SincellObject[["AssociationOfCellsHierarchyWithAGeneSet.pvalue"]]<-NULL
    }
    EffectiveGeneSetSize<-sum(rownames(SincellObject[["expressionmatrix"]]) %in% GeneSet)
    # Assessing the tree for the geneset    
    tmp <- SincellObject[["expressionmatrix"]][(rownames(SincellObject[["expressionmatrix"]]) %in% GeneSet),]
    mySincellObjecttmp <- sc_InitializingSincellObject(tmp)
    if(
      (SincellObject[["method"]]=="euclidean") | (SincellObject[["method"]]=="L1") | (SincellObject[["method"]]=="pearson") | (SincellObject[["method"]]=="MI")
    ){
      mySincellObjecttmp<- sc_distanceObj(mySincellObjecttmp, method=SincellObject[["method"]], bins=SincellObject[["bins"]])
    }
    if(
      (SincellObject[["method"]]=="PCA") | (SincellObject[["method"]]=="ICA") | (SincellObject[["method"]]=="tSNE") | (SincellObject[["method"]]=="classical-MDS") | (SincellObject[["method"]]=="nonmetric-MDS")
    ){
      mySincellObjecttmp<-  sc_DimensionalityReductionObj(mySincellObjecttmp, method=SincellObject[["method"]], dim=SincellObject[["dim"]], MDS.distance=SincellObject[["MDS.distance"]], bins=SincellObject[["bins"]])
    }
    if(((SincellObject[["graph.algorithm"]]=="MST")|(SincellObject[["graph.algorithm"]]=="SST"))&(SincellObject[["graph.using.cells.clustering"]]==TRUE)){
      mySincellObjecttmp<- sc_clusterObj (mySincellObjecttmp, clust.method=SincellObject[["clust.method"]], mutual=SincellObject[["mutual"]], k=SincellObject[["k"]], max.distance=SincellObject[["max.distance"]], shortest.rank.percent=SincellObject[["shortest.rank.percent"]])
      mySincellObjecttmp<- sc_GraphBuilderObj(mySincellObjecttmp, graph.algorithm=SincellObject[["graph.algorithm"]], graph.using.cells.clustering=SincellObject[["graph.using.cells.clustering"]],k=SincellObject[["k"]])
    }
    if(((SincellObject[["graph.algorithm"]]=="MST")|(SincellObject[["graph.algorithm"]]=="SST"))&(SincellObject[["graph.using.cells.clustering"]]==FALSE)){
      mySincellObjecttmp<- sc_GraphBuilderObj(mySincellObjecttmp, graph.algorithm=SincellObject[["graph.algorithm"]], graph.using.cells.clustering=SincellObject[["graph.using.cells.clustering"]],k=SincellObject[["k"]])
    }
    if(SincellObject[["graph.algorithm"]]=="IMC"){
      mySincellObjecttmp<- sc_GraphBuilderObj(mySincellObjecttmp, graph.algorithm=SincellObject[["graph.algorithm"]], graph.using.cells.clustering=SincellObject[["graph.using.cells.clustering"]],k=SincellObject[["k"]])
    }
    node.dist.geneset <- shortest.paths(mySincellObjecttmp[["cellstateHierarchy"]])
    node.dist.originaltree <- shortest.paths(SincellObject[["cellstateHierarchy"]])
    correlation.geneset<-cor(f_distance2vector(node.dist.originaltree),f_distance2vector(node.dist.geneset), method="spearman")
    SincellObject[["AssociationOfCellsHierarchyWithAGeneSet"]]<-correlation.geneset
    if(spearman.rank.threshold<correlation.geneset){
      messeage<-paste("The spearman rank correlation between the original hierarchy \nand the hierarchy obtained when using only the",EffectiveGeneSetSize,"genes \ncommon with gene list",as.name(GeneSet),"is\nr=",correlation.geneset,"\n")
      # Assessing subsampings for p-value
      if(p.value.assessment==TRUE){
        ## Trees from subsets
        f_mclapply <- function(x){
          gene_subset <- sample(1:nrow(SincellObject[["expressionmatrix"]]), EffectiveGeneSetSize)
          tmp <- SincellObject[["expressionmatrix"]][gene_subset,]
          mySincellObjecttmp <- sc_InitializingSincellObject(tmp)
          if(
            (SincellObject[["method"]]=="euclidean") | (SincellObject[["method"]]=="L1") | (SincellObject[["method"]]=="pearson") | (SincellObject[["method"]]=="MI")
          ){
            mySincellObjecttmp<- sc_distanceObj(mySincellObjecttmp, method=SincellObject[["method"]], bins=SincellObject[["bins"]])
          }
          if(
            (SincellObject[["method"]]=="PCA") | (SincellObject[["method"]]=="ICA") | (SincellObject[["method"]]=="tSNE") | (SincellObject[["method"]]=="classical-MDS") | (SincellObject[["method"]]=="nonmetric-MDS")
          ){
            mySincellObjecttmp<-  sc_DimensionalityReductionObj(mySincellObjecttmp, method=SincellObject[["method"]], dim=SincellObject[["dim"]], MDS.distance=SincellObject[["MDS.distance"]], bins=SincellObject[["bins"]])
          }
          if(((SincellObject[["graph.algorithm"]]=="MST")|(SincellObject[["graph.algorithm"]]=="SST"))&(SincellObject[["graph.using.cells.clustering"]]==TRUE)){
            mySincellObjecttmp<- sc_clusterObj (mySincellObjecttmp, clust.method=SincellObject[["clust.method"]], mutual=SincellObject[["mutual"]], k=SincellObject[["k"]], max.distance=SincellObject[["max.distance"]], shortest.rank.percent=SincellObject[["shortest.rank.percent"]])
            mySincellObjecttmp<- sc_GraphBuilderObj(mySincellObjecttmp, graph.algorithm=SincellObject[["graph.algorithm"]], graph.using.cells.clustering=SincellObject[["graph.using.cells.clustering"]],k=SincellObject[["k"]])
          }
          if(((SincellObject[["graph.algorithm"]]=="MST")|(SincellObject[["graph.algorithm"]]=="SST"))&(SincellObject[["graph.using.cells.clustering"]]==FALSE)){
            mySincellObjecttmp<- sc_GraphBuilderObj(mySincellObjecttmp, graph.algorithm=SincellObject[["graph.algorithm"]], graph.using.cells.clustering=SincellObject[["graph.using.cells.clustering"]],k=SincellObject[["k"]])
          }
          if(SincellObject[["graph.algorithm"]]=="IMC"){
            mySincellObjecttmp<- sc_GraphBuilderObj(mySincellObjecttmp, graph.algorithm=SincellObject[["graph.algorithm"]], graph.using.cells.clustering=SincellObject[["graph.using.cells.clustering"]],k=SincellObject[["k"]])
          }
          node.dist <- shortest.paths(mySincellObjecttmp[["cellstateHierarchy"]])
          return(f_distance2vector(node.dist))
        }
        tree.list <-  mclapply(1:num_it, f_mclapply, mc.cores=cores, mc.preschedule = TRUE)
        #if(sum(unlist(lapply(tree.list, is.null)))){
          #while(sum(unlist(lapply(tree.list, is.null)))){
          while(!all(unlist(lapply(tree.list, is.vector))) || !all(unlist(lapply(tree.list, is.numeric))) || any(unlist(lapply(tree.list, is.na))) || any(unlist(lapply(tree.list, is.null))) || length(tree.list)!=num_it || !all(unlist(lapply(tree.list, length))==nrow(SincellObject[["cell2celldist"]])*(nrow(SincellObject[["cell2celldist"]])-1)/2) ){
            # print("WARNING: memory error. Re-executing.")
            tree.list <-  mclapply(1:num_it, f_mclapply, mc.cores=cores, mc.preschedule = TRUE)
          }
        #}
        node.dist <- shortest.paths(SincellObject[["cellstateHierarchy"]])
        correlation<-cor(f_distance2vector(node.dist),matrix(unlist(tree.list),ncol=num_it,nrow=(ncol(SincellObject[["expressionmatrix"]])*(ncol(SincellObject[["expressionmatrix"]])-1)/2)), method="spearman")
        correlation<-sort(append(as.vector(correlation),correlation.geneset),decreasing=TRUE)
        p.value=match(correlation.geneset,correlation)/length(correlation)
        SincellObject[["AssociationOfCellsHierarchyWithAGeneSet.pvalue"]]<-p.value
        messeage<-paste(messeage,"with an empirical p-value=",format(p.value, scientific=TRUE),"\ndrawn from",num_it,"random subsamplings of equal gene set size=",EffectiveGeneSetSize,"\n")
      }
      cat(messeage)
    }
    return(SincellObject)
  }
}
# ## Compare graphs
# 
sc_ComparissonOfGraphs <- function(cellstateHierarchy1,cellstateHierarchy2, ...,graph.names=NULL){
  input_list <- list(...)
  node.dist <- shortest.paths(cellstateHierarchy1)
  node.dist.matrix <- f_distance2vector(node.dist)
  node.dist <- shortest.paths(cellstateHierarchy2)
  node.dist.matrix <- rbind(node.dist.matrix, f_distance2vector(node.dist))
  if(length(input_list)>0){
    for(i in 1:length(input_list)){
      node.dist <- shortest.paths(input_list[[i]])
      node.dist.matrix <- rbind(node.dist.matrix, f_distance2vector(node.dist))
    }
  }
  Comparisson.distance <- as.matrix(1-cor(t(as.matrix(node.dist.matrix)), method="spearman"))
  
  if(is.character(graph.names)&is.vector(graph.names)&(length(graph.names)==(2+length(input_list)))){
    for(i in 1:length(graph.names)){
      rownames(Comparisson.distance)[i]<-graph.names[i]
    }
  }else{
    rownames(Comparisson.distance)[1]<-as.character(1)
    rownames(Comparisson.distance)[2]<-as.character(2)
    if(length(input_list)>0){
      for(i in 1:length(input_list)){
        j=i+2
        rownames(Comparisson.distance)[j]<-as.character(j)
      }
    }
  }
  colnames(Comparisson.distance)<-rownames(Comparisson.distance)

  if(length(input_list)>0){
    plot(hclust(as.dist(Comparisson.distance)))
  }
  return(as.dist(Comparisson.distance))
}



# sc_ComparissonOfGraphs <- function(SincellObject1,SincellObject2, ...,graph.names=NULL){
#   if(is.list(SincellObject1)==FALSE){
#     errormesseage<-paste(deparse(substitute(SincellObject1))," is not a list")
#     stop(errormesseage)
#   }else if(is.null(SincellObject1[["expressionmatrix"]])==TRUE){
#     errormesseage<-paste(deparse(substitute(SincellObject1)),"does not contain a valid expressionmatrix member. Please initialize using function sc_InitializingSincellObject()")
#     stop(errormesseage)
#   }else if(is.null(SincellObject1[["cell2celldist"]])==TRUE){
#     errormesseage<-paste(deparse(substitute(SincellObject1)),"does not contain a valid cell2celldist member. Please assess a cell-2-cell distance matrix using either function sc_distanceObj() or function sc_DimensionalityReductionObj()")
#     stop(errormesseage)
#   }else if(is.null(SincellObject1[["cellstateHierarchy"]])==TRUE){
#     errormesseage<-paste(deparse(substitute(SincellObject1)),"does not contain a valid cellstateHierarchy member. Please assess a cellstateHierarchy using function sc_GraphBuilderObj()")
#     stop(errormesseage)
#   }else if(is.list(SincellObject2)==FALSE){
#     errormesseage<-paste(deparse(substitute(SincellObject2))," is not a list")
#     stop(errormesseage)
#   }else if(is.null(SincellObject2[["expressionmatrix"]])==TRUE){
#     errormesseage<-paste(deparse(substitute(SincellObject2)),"does not contain a valid expressionmatrix member. Please initialize using function sc_InitializingSincellObject()")
#     stop(errormesseage)
#   }else if(is.null(SincellObject2[["cell2celldist"]])==TRUE){
#     errormesseage<-paste(deparse(substitute(SincellObject2)),"does not contain a valid cell2celldist member. Please assess a cell-2-cell distance matrix using either function sc_distanceObj() or function sc_DimensionalityReductionObj()")
#     stop(errormesseage)
#   }else if(is.null(SincellObject2[["cellstateHierarchy"]])==TRUE){
#     errormesseage<-paste(deparse(substitute(SincellObject2)),"does not contain a valid cellstateHierarchy member. Please assess a cellstateHierarchy using function sc_GraphBuilderObj()")
#     stop(errormesseage)
#   }else{
#     input_list <- list(...)
#     for(i in 1:length(input_list)){
#       if(is.list(input_list[[i]])==FALSE){
#         errormesseage<-paste(deparse(substitute(input_list[[i]]))," is not a list")
#         stop(errormesseage)
#       }else if(is.null(input_list[[i]][["expressionmatrix"]])==TRUE){
#         errormesseage<-paste(deparse(substitute(input_list[[i]])),"does not contain a valid expressionmatrix member. Please initialize using function sc_InitializingSincellObject()")
#         stop(errormesseage)
#       }else if(is.null(input_list[[i]][["cell2celldist"]])==TRUE){
#         errormesseage<-paste(deparse(substitute(input_list[[i]])),"does not contain a valid cell2celldist member. Please assess a cell-2-cell distance matrix using either function sc_distanceObj() or function sc_DimensionalityReductionObj()")
#         stop(errormesseage)
#       }else if(is.null(input_list[[i]][["cellstateHierarchy"]])==TRUE){
#         errormesseage<-paste(deparse(substitute(input_list[[i]])),"does not contain a valid cellstateHierarchy member. Please assess a cellstateHierarchy using function sc_GraphBuilderObj()")
#         stop(errormesseage)
#       }
#     }
#     node.dist <- shortest.paths(SincellObject1[["cellstateHierarchy"]])
#     node.dist.matrix <- f_distance2vector(node.dist)
#     node.dist <- shortest.paths(SincellObject2[["cellstateHierarchy"]])
#     node.dist.matrix <- rbind(node.dist.matrix, f_distance2vector(node.dist))
#     for(i in 1:length(input_list)){
#       node.dist <- shortest.paths(input_list[[i]][["cellstateHierarchy"]])
#       node.dist.matrix <- rbind(node.dist.matrix, f_distance2vector(node.dist))
#     }
#     Comparisson.distance <- as.matrix(1-cor(t(as.matrix(node.dist.matrix)), method="pearson"))
#     
#     if(is.character(graph.names)&is.vector(graph.names)&(length(graph.names)==(2+length(input_list)))){
#       for(i in 1:length(graph.names)){
#         rownames(Comparisson.distance)[i]<-graph.names[i]
#       }
#     }else{
#       rownames(Comparisson.distance)[1]<-as.character(1)
#       rownames(Comparisson.distance)[2]<-as.character(2)
#       for(i in 1:length(input_list)){
#         j=i+2
#         rownames(Comparisson.distance)[j]<-as.character(j)
#       }
#     }
#     colnames(Comparisson.distance)<-rownames(Comparisson.distance)
#     plot(hclust(as.dist(Comparisson.distance)))
#     return(as.dist(Comparisson.distance))
#   }
# }


