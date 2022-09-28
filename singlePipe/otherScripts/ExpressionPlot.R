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
       BiocManager::install("data.table")}
library(data.table)

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

make_option(c("-e", "--ExpressionPlotGenes"), 
    type="character", 
    default=NA,
    help="list of genes",
    metavar="character"
           )

)

parser<- OptionParser(option_list=option_list)
opt <- parse_args(parser)
inputFile=opt$inputFile
outputDir=opt$outputDir
ExpressionPlotGenes=opt$ExpressionPlotGenes

ExpressionPlotGenes = unlist(
	strsplit(
	    ExpressionPlotGenes, 
	    ","))
		    
set.seed(1)
sobj = readRDS(inputFile)

jpeg(paste(outputDir, '/', 'ExpressionPlot.jpeg', sep = ""),
	unit = "in",
	res = 300,
	height = 8,
	width =  8)
	
FeaturePlot(sobj, features = ExpressionPlotGenes)

dev.off()

