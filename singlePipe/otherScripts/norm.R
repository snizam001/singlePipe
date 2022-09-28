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
           
make_option(c("-m", "--normMethod"), 
    type="character", 
    default='LogNormalize',
    help="Normalization method: LogNormalize, CLR, RC. Default is LogNormalize.", 
    metavar="character"),
        
make_option(c("-s", "--normScalingFactor"), 
    type="integer", 
    default=10000,
    help="Scaling factor for the normalization. Default is 10000."
           )
)

parser<- OptionParser(option_list=option_list)
opt <- parse_args(parser)
inputFile=opt$inputFile
outputDir=opt$outputDir
normMethod=opt$normMethod
normScalingFactor=opt$normScalingFactor

set.seed(1)
sobj = readRDS(inputFile)
sobj <- NormalizeData(sobj, normalization.method = normMethod, scale.factor = normScalingFactor)
saveRDS(sobj, file = paste(outputDir, '/', 'sobj.rds', sep = ""))


