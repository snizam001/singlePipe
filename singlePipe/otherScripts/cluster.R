#___________________  
#!/usr/bin/env Rscript
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
if(!require(SummarizedExperiment)){
       BiocManager::install("SummarizedExperiment")}
library(SummarizedExperiment)

if(!require(optparse)){
        install.packages("optparse")}
library(optparse)

if(!require(Matrix)){
        install.packages("Matrix")}
library(Matrix)

if(!require(Seurat)){
        install.packages("Seurat")}
library(Seurat)

if(!require(data.table)){
       install.packages("data.table")}
library(data.table)

if(!require(tidyverse)){
       install.packages("tidyverse")}
library(tidyverse)

if(!require(HGNChelper)){
       install.packages("HGNChelper")}
library(HGNChelper)

if(!require(uwot)){
        install.packages("uwot")}
library(uwot)

if(!require(edgeR)){
        install.packages("edgeR")}
library(edgeR)

if(!require(matrixStats)){
        install.packages("matrixStats")}
library(matrixStats)

if(!require(Rcpp)){
        install.packages("Rcpp")}
library(Rcpp)

#--- 
option_list = list(

make_option(c("-i", "--inputFile"), 
    type="character", 
    default=NA,
    help="Input file", 
    metavar="character"),

make_option(c("-o", "--outputDir"), 
    type="character", 
    default=NA,
    help="output directory",
    metavar="character"
           ),

make_option(c("-t", "--threads"), 
    type="integer", 
    default=1,
    help="How many threads? Default is 1"),
                      
make_option(c("-d", "--dimPCA"), 
    type="integer", 
    default=40,
    help="How many eigenvectors or dimension should be used to find neighbors. Default is 40"),
        
make_option(c("-r", "--resPCA"), 
    type="character", 
    default='0.5',
    help="Value to find clusters. Default is 0.5"
           ),

make_option(c("-n", "--clusterPlotN1"), 
    type="integer", 
    default=1,
    help="Dimension number of UMAP/tsne/PCA to generate graph. Default is 1"),

make_option(c("-m", "--clusterPlotN2"), 
    type="integer", 
    default=2,
    help="Dimension number of UMAP/tsne/PCA to generate graph. Default is 2"),

make_option(c("-c", "--clustMethod"), 
    type="character", 
    default="Seurat",
    help="Should this pipeline use default method (Seurat) or method of Granja et al (2019). Default is Seurat."),
             ########################### you can add new method of clustering here - part0    

make_option(c("-s", "--selMethod"), 
    type="character", 
    default='vst',
    help="selection method for variables. Options are vst, dispersion and mean.var.plot. For detail see: https://satijalab.org/seurat/reference/findvariablefeatures. Default is vst.",
    metavar="character"
           ),
           
make_option(c("-v", "--VariableNumberFeatures"), 
    type="integer", 
    default=2000,
    help="Number of the features which are variable and used for dimension reduction. Default is 2000."
           )
)

parser<- OptionParser(option_list=option_list)
opt <- parse_args(parser)
inputFile=opt$inputFile
outputDir=opt$outputDir
dimPCA=opt$dimPCA
resPCA=opt$resPCA
clusterPlotN1=opt$clusterPlotN1
clusterPlotN2=opt$clusterPlotN2
clustMethod=opt$clustMethod
threads=opt$threads
VariableNumberFeatures=opt$VariableNumberFeatures
selMethod=opt$selMethod
set.seed(1)

#Binarize Sparse Matrix
binarizeMat <- function(mat){
    mat@x[mat@x > 0] <- 1
    mat
}

#LSI Adapted from fly-atac with information for re-projection analyses
calcLSI <- function(mat, nComponents = 50, binarize = TRUE, nFeatures = NULL){

    set.seed(1)

    #TF IDF LSI adapted from flyATAC
    if(binarize){
        message(paste0("Binarizing matrix..."))
        mat@x[mat@x > 0] <- 1 
    }

    if(!is.null(nFeatures)){
        message(paste0("Getting top ", nFeatures, " features..."))
        idx <- head(order(Matrix::rowSums(mat), decreasing = TRUE), nFeatures)
        mat <- mat[idx,] 
    }else{
        idx <- which(Matrix::rowSums(mat) > 0)
        mat <- mat[idx,]
    }

    #Calc RowSums and ColSums
    colSm <- Matrix::colSums(mat)
    rowSm <- Matrix::rowSums(mat)

    #Calc TF IDF
    message("Computing Term Frequency IDF...")
    freqs <- t(t(mat)/colSm)
    idf   <- as(log(1 + ncol(mat) / rowSm), "sparseVector")
    tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs

    #Calc SVD then LSI
    message("Computing SVD using irlba...")
    svd <- irlba::irlba(tfidf, nComponents, nComponents)
    svdDiag <- matrix(0, nrow=nComponents, ncol=nComponents)
    diag(svdDiag) <- svd$d
    matSVD <- t(svdDiag %*% t(svd$v))
    rownames(matSVD) <- colnames(mat)
    colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))

    #Return Object
    out <- list(
        matSVD = matSVD, 
        rowSm = rowSm, 
        colSm = colSm, 
        idx = idx, 
        svd = svd, 
        binarize = binarize, 
        nComponents = nComponents,
        date = Sys.Date(),
        seed = 1)

    out

}

#Clustering function using seurat SNN (Seurat v2.3.4)
seuratSNN <- function(matSVD, dims.use = 1:50, ...){
  set.seed(1)
  message("Making Seurat Object...")
  mat <- matrix(rnorm(nrow(matSVD) * 3, 1000), ncol = nrow(matSVD), nrow = 3)
  colnames(mat) <- rownames(matSVD)
  rownames(mat) <- colnames(matSVD)[1:3]
  obj <- Seurat::CreateSeuratObject(mat, project='scATAC', min.cells=0, min.genes=0)
  obj [['pca']] = CreateDimReducObject(embeddings  = matSVD, key = "PC")
  obj <- FindNeighbors(obj, dims = 1:25)
  obj <- FindClusters(obj, resolution = as.numeric(0.2))
#  obj <- Seurat::SetDimReduction(object = obj, reduction.type = "pca", slot = "cell.embeddings", new.data = matSVD)
#  obj <- Seurat::SetDimReduction(object = obj, reduction.type = "pca", slot = "key", new.data = "PC")
#  obj <- Seurat::FindClusters(object = obj, reduction.type = "pca", dims.use = dims.use, print.output = TRUE, ...)
  clust <- obj@meta.data[,ncol(obj@meta.data)]
  paste0("Cluster",match(clust, unique(clust)))
}

#Sparse Variances Rcpp
sourceCpp(code='
  #include <Rcpp.h>

  using namespace Rcpp;
  using namespace std;

  // [[Rcpp::export]]
  Rcpp::NumericVector computeSparseRowVariances(IntegerVector j, NumericVector val, NumericVector rm, int n) {
    const int nv = j.size();
    const int nm = rm.size();
    Rcpp::NumericVector rv(nm);
    Rcpp::NumericVector rit(nm);
    int current;
    // Calculate RowVars Initial
    for (int i = 0; i < nv; ++i) {
      current = j(i) - 1;
      rv(current) = rv(current) + (val(i) - rm(current)) * (val(i) - rm(current));
      rit(current) = rit(current) + 1;
    }
    // Calculate Remainder Variance
    for (int i = 0; i < nm; ++i) {
      rv(i) = rv(i) + (n - rit(i))*rm(i)*rm(i);
    }
    rv = rv / (n - 1);
    return(rv);
  }'
)

#Compute Fast Sparse Row Variances
sparseRowVariances <- function (m){
    rM <- Matrix::rowMeans(m)
    rV <- computeSparseRowVariances(m@i + 1, m@x, rM, ncol(m))
    return(rV)
}

#Helper function for summing sparse matrix groups
groupSums <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
    stopifnot(!is.null(groups))
    stopifnot(length(groups) == ncol(mat))
    gm <- lapply(unique(groups), function(x) {
        if (sparse) {
            Matrix::rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
        else {
            rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
    }) %>% Reduce("cbind", .)
    colnames(gm) <- unique(groups)
    return(gm)
}

#Optimized LSI for scRNA-seq analysis
optimizeLSI <- function(sobj, scaleTo = 10000, priorCount = 3, pcsUse = 1:25, 
    resolution = c(0.2, 0.4, 0.8), varFeatures = c(2500, 2500, 2500), seed = 1){
    mat = sobj[['RNA']]@counts
    set.seed(seed)
   stopifnot(length(resolution) > 1)

    #Initialize List
    lsiOut <- list()

    #Initial LSI uses variances that are across all single cells and will have larger batch relationships
    i <- 1
    message("Initial LSI...")
    matNorm <- t(t(mat)/Matrix::colSums(mat)) * scaleTo
    matNorm@x <- log2(matNorm@x + 1)
#    matNorm@x <- sobj[['RNA']]@data
    idVarFeatures <- head(order(sparseRowVariances(matNorm),decreasing=TRUE), varFeatures[i])
    lsiObj <- calcLSI(mat[idVarFeatures,], binarize = FALSE, nComponents = max(pcsUse))
    clusters <- seuratSNN(lsiObj$matSVD, dims.use = pcsUse, resolution = resolution[i], n.start = 10, print.output = FALSE)

    #Store
    lsiOut[[paste0("iter", i)]] <- list(
        lsiMat = lsiObj$matSVD, 
        varFeatures = idVarFeatures, 
        clusters = clusters
        )

    for(i in seq(2, length(varFeatures))){

       message(sprintf("Additional LSI %s...", i))

        #Run LSI
        clusterMat <- edgeR::cpm(groupSums(mat, clusters, sparse = TRUE), log=TRUE, prior.count = priorCount)
        idVarFeatures <- head(order(rowVars(clusterMat), decreasing=TRUE), varFeatures[i])
        lsiObj <- calcLSI(mat[idVarFeatures,], binarize = FALSE, nComponents = max(pcsUse))
        clusters <- seuratSNN(lsiObj$matSVD, dims.use = pcsUse, resolution = resolution[i], n.start = 10, print.output = FALSE)

        if(i == length(varFeatures)){
            #Save All Information from LSI Attempt
            lsiOut[[paste0("iter", i)]] <- list(
                lsiObj = lsiObj, 
                varFeatures = idVarFeatures, 
                clusters = clusters,
                matNorm = matNorm
                )
        }else{
            lsiOut[[paste0("iter", i)]] <- list(
                lsiMat = lsiObj$matSVD, 
                varFeatures = idVarFeatures, 
                clusters = clusters
                )
        }

    }

    return(lsiOut)

}

sobj = readRDS(inputFile)

if (clustMethod == "Seurat") {

	sobj <- FindNeighbors(sobj, dims = 1:dimPCA)
	sobj <- FindClusters(sobj, resolution = as.numeric(resPCA))
	sobj <- RunUMAP(sobj, dims = 1:dimPCA)

	print("createing")
	z=paste(outputDir, '/', clustMethod, 'DimPlot-UMAP.jpeg', sep = "")
	print(z)
	jpeg(paste(outputDir, '/', clustMethod, 'DimPlot-UMAP.jpeg', sep = ""),
		unit = "in",
		res = 300,
		height = 6,
		width =  6)
		
	DimPlot(sobj, reduction = "umap", dims = c(clusterPlotN1,clusterPlotN2))

	dev.off()

	sobj <- RunTSNE(sobj, dims = 1:dimPCA)

	jpeg(paste(outputDir, '/', clustMethod, 'DimPlot-tSNE.jpeg', sep = ""),
		unit = "in",
		res = 300,
		height = 6,
		width =  6)
		
	DimPlot(sobj, reduction = "tsne", dims = c(clusterPlotN1,clusterPlotN2))

	dev.off()

	sobj <- RunPCA(sobj, dims = 1:dimPCA)

	jpeg(paste(outputDir, '/', clustMethod, 'DimPlot-PCA.jpeg', sep = ""),
		unit = "in",
		res = 300,
		height = 6,
		width =  6)
		
	DimPlot(sobj, reduction = "pca", dims = c(clusterPlotN1,clusterPlotN2))

	dev.off()
		
} else if (clustMethod == "GranjaEtal")  {
	
	nTop = rep(VariableNumberFeatures,3)
	resolution = as.numeric(unlist(
		strsplit(
		    resPCA, 
		    ",")))	
	#c(0.2,0.5,1.0)
	nPCs=1:dimPCA
	lsiObj <- optimizeLSI(sobj, 
	  resolution = resolution, 
	  pcsUse = nPCs,
	  varFeatures = nTop)
	sobj@misc$optimizeLSI <- lsiObj
	sobj@misc$matSVD <- lsiObj[[length(lsiObj)]][[1]][[1]] #Last one
	sobj@misc$variableGenes <- rownames(sobj)[lsiObj[[length(lsiObj)]]$varFeatures] #Variable genes
	sobj@misc$clusters <- sobj@misc$optimizeLSI[[length(lsiObj)]]$clusters
	sobj@meta.data$seurat_clusters  <- sobj@misc$optimizeLSI[[length(lsiObj)]]$clusters
	sobj@meta.data$seurat_clusters  <- as.numeric(sub("Cluster","",sobj@meta.data$seurat_clusters))


	matSVD <- sobj@misc$matSVD
	clusters <- sobj@misc$clusters

	#Set Seed and perform UMAP on LSI-SVD Matrix
	set.seed(1)
	uwotUmap <- uwot::umap(
	    matSVD, 
	    n_neighbors = 35, 
	    min_dist = 0.45, 
	    metric = "euclidean", 
	    n_threads = threads, 
	    verbose = TRUE, 
	    ret_nn = TRUE,
	    ret_model = TRUE
	    )

	jpeg(paste(outputDir, '/', clustMethod, 'DimPlot-UMAP.jpeg', sep = ""),
		unit = "in",
		res = 300,
		height = 6,
		width =  6)
		
	df <- data.frame(
	    x = uwotUmap[[1]][,1],
	    y = -uwotUmap[[1]][,2], 
	    color = clusters
	    )
	ggplot(df,aes(x,y,color=color)) + 
	    geom_point() + 
	    theme_bw() + 
	    xlab("UMAP Dimension 1") + 
	    ylab("UMAP Dimension 2")
	dev.off()

	#Add UMAP coordinates to column data in summarized experiment
	sobj@meta.data$UMAP1 <- uwotUmap[[1]][,1]
	sobj@meta.data$UMAP2 <- uwotUmap[[1]][,2]

	#Save UMAP embedding
	saveRDS(uwotUmap, file = paste(outputDir, '/', clustMethod, 'Umap-model.rds', sep = ""))
} ########################### you can add new method of clustering here - part1
saveRDS(sobj, file = paste(outputDir, '/', 'sobj.rds', sep = ""))


