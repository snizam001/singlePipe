#!/usr/bin/env python3

import pkg_resources
import sys
import subprocess
import argparse
import os
from itertools import combinations
from termcolor import colored
from datetime import datetime
from itertools import product
import itertools
import pkg_resources as psource
from multiprocessing import Pool
from itertools import repeat
from singlePipe import universal

#_____________________________________________________________________________________________________

__description__ = """

singleCell analysis pipeline for medical University of Vienna (version 1.0).
------------------------------------------------------------------------------------
Use this command after installation for usage instruction: 

singlePipe --help

This pipeline has different modes:
1) qc: 
2) norm:
3) featureSelectDimentionReduction
4) clustering
5) diffExprGenes
6) annoT

"""

epilogue = """
Authors: 
        (a) Sheikh Nizamuddin:  
        snizam001@gmail.com, 
        (b) ..
"""
#_____________________________________________________________________________________________________

parser = argparse.ArgumentParser('singlePipe',
                                 usage = __description__,
                                 epilog = colored(epilogue, 'yellow', attrs = ['bold']))

reqNamed = parser.add_argument_group('required arguments\n   _____________________')

reqNamed.add_argument("--modes",
                      help=colored("run mode.", 'green', attrs = ['bold']) + " Multiple modes can be provided as "+
                      "comma seperated values e.g. --modes mode1,mode2. Options are "+
                      "qc,norm,featureSelectDimentionReduction,clustering,diffExprGenes,annoT,.....",
                      type=str,
                      required=True)

reqNamed.add_argument("--outputdir",
                      help="output directory (give full path)",
                      type=str,
                      required=True)

comNamed = parser.add_argument_group('common arguments\n   _____________________')

comNamed.add_argument("--inputData",
                      help="input data (give full path)",
                      default="NA")

comNamed.add_argument("--inputDataType",
                      help="input data type. Is it tsv, csv, rds, fastq or bamfiles. Default is rds.",
                      choices=['tsv',
                               'csv',
                               'rds',
                               '..'
                               ],
                      default="rds")
                      
comNamed.add_argument("--threads", 
                    help="number of threads (default is: total CPU - 2)",
                    type=int,
                    default=(os.cpu_count()-2))

parser.add_argument("--nFeatureMin", 
                    help=colored("qc mode: ", 'green', attrs = ['bold']) + 
                    "Cells less than this feature will be filtered. Default value is 200",
                    default='200')

parser.add_argument("--nFeatureMax", 
                    help=colored("qc mode: ", 'green', attrs = ['bold']) + 
                    "Cells more than this feature will be filtered. Default value is 2500",
                    default='2500')
                    
parser.add_argument("--percMT", 
                    help=colored("qc mode: ", 'green', attrs = ['bold']) + 
                    "Cells more than this percentage of mitochondria will be filtered. Default value is 5 percent",
                    default='5')

parser.add_argument("--normMethod", 
                    help=colored("norm (normalization) mode: ", 'green', attrs = ['bold']) + 
                    "Normalization method: LogNormalize, CLR, RC. Default is LogNormalize.",
                    default='LogNormalize')

parser.add_argument("--normScalingFactor", 
                    help=colored("norm (normalization) mode: ", 'green', attrs = ['bold']) + 
                    "Scaling factor for the normalization. Default is 10000.",
                    default='10000')

parser.add_argument("--selMethod", 
                    help=colored("featureSelectDimentionReduction and cluster mode: ", 'green', attrs = ['bold']) + 
                    "selection method for variables. Options are vst, dispersion and mean.var.plot. "+
                    "For detail see: https://satijalab.org/seurat/reference/findvariablefeatures. Default is vst.",
                    default='vst')

parser.add_argument("--VariableNumberFeatures", 
                    help=colored("featureSelectDimentionReduction and clustering mode: ", 'green', attrs = ['bold']) + 
                    "Number of the features/genes which are highly variable and used for dimension reduction. Default is 2000.",
                    default='2000')

parser.add_argument("--dimPCA", 
                    help=colored("clustering mode: ", 'green', attrs = ['bold']) + 
                    "How many eigenvectors or dimension should be used to find neighbors. Default is 40."+
                    " This number should be identified from the elbow plot which is generated in the featureSelectDimentionReduction mode",
                    default='40')

parser.add_argument("--resPCA", 
                    help=colored("clustering mode: ", 'green', attrs = ['bold']) + 
                    "Value to find clusters. Default is 0.5. If --clustMethod is "+
                    "GranjaEtal, give three comma seperated values e.g. 0.2,0.6,1.0",
                    default='0.5')

parser.add_argument("--clusterPlotN1", 
                    help=colored("clustering and annoT mode: ", 'green', attrs = ['bold']) + 
                    "Dimension number of UMAP/tsne/PCA to generate graph. Default is 1",
                    default='1')

parser.add_argument("--clusterPlotN2", 
                    help=colored("clustering and annoT mode: ", 'green', attrs = ['bold']) + 
                    "Dimension number of UMAP/tsne/PCA to generate graph. Default is 2",
                    default='2')

parser.add_argument("--clustMethod", 
                    help=colored("clustering and diffExprGenes mode: ", 'green', attrs = ['bold']) + 
                    "Should this pipeline use default method (Seurat) or method of Granja et al (2019)."+
                    " Default is Seurat.",
                    choices = ['Seurat','GranjaEtal','...'],
                    default='Seurat')

parser.add_argument("--findMarkersLogFC", 
                    help=colored("diffExprGenes mode: ", 'green', attrs = ['bold']) + 
                    "Limit testing to genes which show, on average,"+
                    " at least X-fold difference (log-scale) between the two groups of cells. Default is 0.25.",
                    default='0.25')

parser.add_argument("--findMarkersMinPCT", 
                    help=colored("diffExprGenes mode: ", 'green', attrs = ['bold']) + 
                    "Only test genes that are detected in a minimum fraction of " +
                    "min.pct cells in either of the two populations. Default is 0.25.",
                    default='0.25')

parser.add_argument("--ScTypeDB", 
                    help=colored("annoT mode: ", 'green', attrs = ['bold']) + 
                    "Input excel file containing information of markers for cells of a particular tissue"+
                    ". By default it uses internal datasset of package sc-type. But if you are using your own"+
                    " marker then just prepare an input XLSX file in the same format as our DB file. "+
                    "DB file should contain four columns (tissueType - tissue type, cellName - cell type,"+ 
                    "geneSymbolmore1 - positive marker genes, geneSymbolmore2 - marker genes not expected to "+
                    "be expressed by a cell type)",
                    default='NA')

parser.add_argument("--cellAnnotateTissue", 
                    help=colored("annoT mode: ", 'green', attrs = ['bold']) + 
                    "cells of this tissue should be used for annotation. Default is Immune system."+
                    " By default choices are Adrenal, Brain,Eye,Heart,Immune system,Intestine,"+
                    "Kidney,Liver,Lung,Muscle,Pancreas,Placenta,Spleen,Stomach and Thymus"+
                    ". If you are using your own list with --ScTypeDB, then define tissue type here.",
                    default='Immune system')

parser.add_argument("--ExpressionPlotGenes", 
                    help=colored("ExpressionPlot mode: ", 'green', attrs = ['bold']) + 
                    "Expression of gene on UMAP plot. Give the list of genes as comma seperated values",
                    default='NA')


args = parser.parse_args()
#------
outputdir=args.outputdir
modes=args.modes
#------
inputData=args.inputData
inputDataType=args.inputDataType
threads=args.threads
#----
nFeatureMin=args.nFeatureMin
nFeatureMax=args.nFeatureMax
percMT=args.percMT
#----
normMethod=args.normMethod
normScalingFactor=args.normScalingFactor
#----
selMethod=args.selMethod
VariableNumberFeatures=args.VariableNumberFeatures
#----
dimPCA=args.dimPCA
resPCA=args.resPCA
clusterPlotN1=args.clusterPlotN1
clusterPlotN2=args.clusterPlotN2
clustMethod=args.clustMethod
#----
findMarkersLogFC=args.findMarkersLogFC
findMarkersMinPCT=args.findMarkersMinPCT
#----
cellAnnotateTissue=args.cellAnnotateTissue
ScTypeDB=args.ScTypeDB
ExpressionPlotGenes=args.ExpressionPlotGenes
#--- check python packages installed or not

#-- checking if other neccessary packages installed or not 
      
#-- all about required arguments

#-- checking if modes are written correctly

for mode in modes.split(','):
    if mode == 'qc' or mode == 'norm' or mode == 'featureSelectDimentionReduction' or mode == 'clustering' or mode == 'diffExprGenes' or mode == 'annoT':
        print('')
    else:
        print(colored(mode + ': check spelling. Is this correct mode?',
                      'red',
                      attrs=['bold'])
             )
        exit()
#-- checking files/folders present or not?
if not os.path.exists(inputData):
    print(colored(inputData + ': file does not exist.',
                  'red', 
                  attrs=['bold'])
         )
    
if not os.path.exists(outputdir):
    print(colored(outputdir + ': directory does not exist, creating ..',
                  'red', 
                  attrs=['bold'])
         )
    os.makedirs(outputdir)
    
if not os.path.exists(outputdir+'/'+'QC'):
    print(colored(outputdir +'/'+'QC' + ': directory does not exist, creating ..',
                  'red', 
                  attrs=['bold'])
         )
    os.makedirs(outputdir+'/'+'QC')

if not os.path.exists(outputdir+'/'+'norm'):
    print(colored(outputdir +'/'+'norm' + ': directory does not exist, creating ..',
                  'red', 
                  attrs=['bold'])
         )
    os.makedirs(outputdir+'/'+'norm')

if not os.path.exists(outputdir+'/'+'featureSelectDimentionReduction'):
    print(colored(outputdir +'/'+'featureSelectDimentionReduction' + ': directory does not exist, creating ..',
                  'red', 
                  attrs=['bold'])
         )
    os.makedirs(outputdir+'/'+'featureSelectDimentionReduction')

if not os.path.exists(outputdir+'/'+'Cluster'):
    print(colored(outputdir +'/'+'Cluster' + ': directory does not exist, creating ..',
                  'red', 
                  attrs=['bold'])
         )
    os.makedirs(outputdir+'/'+'Cluster')

if not os.path.exists(outputdir+'/'+'diffExprGenes'):
    print(colored(outputdir +'/'+'diffExprGenes' + ': directory does not exist, creating ..',
                  'red', 
                  attrs=['bold'])
         )
    os.makedirs(outputdir+'/'+'diffExprGenes')    

if not os.path.exists(outputdir+'/'+'annoT'):
    print(colored(outputdir +'/'+'annoT' + ': directory does not exist, creating ..',
                  'red', 
                  attrs=['bold'])
         )
    os.makedirs(outputdir+'/'+'annoT')   

if not os.path.exists(outputdir+'/'+'ExpressionPlot'):
    print(colored(outputdir +'/'+'ExpressionPlot' + ': directory does not exist, creating ..',
                  'red', 
                  attrs=['bold'])
         )
    os.makedirs(outputdir+'/'+'ExpressionPlot') 
    
#-- logfile 

with open(outputdir+'/'+'log.txt','w') as logfile:
    logfile.write(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))

#---
def main ():
    for mode in modes.split(','):
        #-- qc pass
        #______________________________________________________________________________________________________
        if mode == 'qc':
            print(colored('==============> '+
                          'QC mode',
                          'green',
                          attrs=['bold'])
                 )
            myqc = psource.resource_filename(__name__, "otherScripts/qc.R")
            cmd=['Rscript', myqc,
                 '--inputFile', inputData,
                 '--inType', inputDataType,
                 '--outputDir', outputdir+'/'+'QC',
                 '--nFeatureMin', nFeatureMin,
                 '--nFeatureMax', nFeatureMax,
                 '--percMT', percMT
                ]
            universal.run_cmd(cmd,outputdir)
            
        #-- normalization
        #______________________________________________________________________________________________________
        if mode == 'norm':
            print(colored('==============> '+
                          'norm mode',
                          'green',
                          attrs=['bold'])
                 )
            mynorm = psource.resource_filename(__name__, "otherScripts/norm.R")
            if 'qc' in modes.split(','):
                inNorm = outputdir+'/'+'QC'+'/'+'sobj.rds'
            else:
                inNorm = inputData

            cmd=['Rscript', mynorm,
                 '--inputFile', inNorm,
                 '--outputDir', outputdir+'/'+'norm',
                 '--normMethod', normMethod,
                 '--normScalingFactor', normScalingFactor
                ]
            universal.run_cmd(cmd,outputdir)

        #-- featureSelectDimentionReduction
        #______________________________________________________________________________________________________
        if mode == 'featureSelectDimentionReduction':
            print(colored('==============> '+
                          'featureSelectDimentionReduction mode. ' + 
                          'Go and drink coffee, this step will take time :-)',
                          'green',
                          attrs=['bold'])
                 )
            myFeature = psource.resource_filename(__name__, "otherScripts/featureSelectDimentionReduction.R")
            if 'norm' in modes.split(','):
                inFeature = outputdir+'/'+'norm'+'/'+'sobj.rds'
            else:
                inFeature = inputData

            cmd=['Rscript', myFeature,
                 '--inputFile', inFeature,
                 '--outputDir', outputdir+'/'+'featureSelectDimentionReduction',
                 '--selMethod', selMethod,
                 '--VariableNumberFeatures', VariableNumberFeatures
                ]
            universal.run_cmd(cmd,outputdir)

            print(colored('==============> '+
                          'Open ElbowPlot.jpeg in the your output folder and'+
                          ' find number of eigenvectors required for the '+
                          'downstream analysis. The right number is on x-axis after which plateau reached.',
                          'green',
                          attrs=['bold'])
                 )
            
        #-- clustering
        #______________________________________________________________________________________________________
        if mode == 'clustering':
            print(colored('==============> '+
                          'clustering mode',
                          'green',
                          attrs=['bold'])
                 )
            myCluster = psource.resource_filename(__name__, "otherScripts/cluster.R")
            if 'featureSelectDimentionReduction' in modes.split(','):
                inCluster = outputdir+'/'+'featureSelectDimentionReduction'+'/'+'sobj.rds'
            else:
                inCluster = inputData
            
            if clustMethod == "GranjaEtal" and len(resPCA.split(','))!=3:
                print(colored('Method '+clustMethod+" needs three --resPCA values. Since, it is not given"+
                              " automatically setting 0.2,0.6,1.0.",
                              'green')
                     )
                cmd=['Rscript', myCluster,
                     '--inputFile', inCluster,
                     '--outputDir', outputdir+'/'+'Cluster',
                     '--dimPCA', dimPCA,
                     '--resPCA', "0.2,0.6,1.0",
                     '--clusterPlotN1', clusterPlotN1,
                     '--clusterPlotN2', clusterPlotN2,
                     '--clustMethod', clustMethod,
                     '--threads', str(threads),
                     '--VariableNumberFeatures', VariableNumberFeatures,
                     '--selMethod', selMethod
                    ]
            else:
                cmd=['Rscript', myCluster,
                     '--inputFile', inCluster,
                     '--outputDir', outputdir+'/'+'Cluster',
                     '--dimPCA', dimPCA,
                     '--resPCA', resPCA,
                     '--clusterPlotN1', clusterPlotN1,
                     '--clusterPlotN2', clusterPlotN2,
                     '--clustMethod', clustMethod,
                     '--threads', str(threads),
                     '--VariableNumberFeatures', VariableNumberFeatures,
                     '--selMethod', selMethod
                    ]
            universal.run_cmd(cmd,outputdir)

        #-- differentially expressed genes among clusters
        #______________________________________________________________________________________________________
        if mode == 'diffExprGenes':
            print(colored('==============> '+
                          'diffExprGenes mode',
                          'green',
                          attrs=['bold'])
                 )
            mydiffExprGenes = psource.resource_filename(__name__, "otherScripts/diffExprGenes.R")
            if 'clustering' in modes.split(','):
                indiffExprGenes = outputdir+'/'+'Cluster'+'/'+'sobj.rds'
            else:
                indiffExprGenes = inputData

            cmd=['Rscript', mydiffExprGenes,
                 '--inputFile', indiffExprGenes,
                 '--outputDir', outputdir+'/'+'diffExprGenes',
                 '--findMarkersLogFC', findMarkersLogFC,
                 '--findMarkersMinPCT', findMarkersMinPCT ,
                 '--clustMethod', clustMethod
                ]
            universal.run_cmd(cmd,outputdir)

        #-- anotating cells with known markers
        #______________________________________________________________________________________________________
        if mode == 'annoT':
            print(colored('==============> '+
                          'annoT mode',
                          'green',
                          attrs=['bold'])
                 )
            myannoT = psource.resource_filename(__name__, "otherScripts/annoT.R")
            gene_sets_prepare = psource.resource_filename(__name__, "otherScripts/gene_sets_prepare.R")
            sctype_score =  psource.resource_filename(__name__, "otherScripts/sctype_score_.R")
            if ScTypeDB == "NA":
                ScTypeDB_full = psource.resource_filename(__name__, "data/ScTypeDB_full.xlsx")
            else:
                ScTypeDB_full = ScTypeDB
            if 'clustering' or diffExprGenes in modes.split(','):
                inmyannoT = outputdir+'/'+'Cluster'+'/'+'sobj.rds'
            else:
                inmyannoT = inputData

            cmd=['Rscript', myannoT,
                 '--inputFile', inmyannoT,
                 '--outputDir', outputdir+'/'+'annoT',
                 '--gene_sets_prepare', gene_sets_prepare,
                 '--sctype_score', sctype_score,
                 '--ScTypeDB_full', ScTypeDB_full,
                 '--cellAnnotateTissue', cellAnnotateTissue,
                 '--clusterPlotN1', clusterPlotN1,
                 '--clusterPlotN2', clusterPlotN2,
                 '--clustMethod', clustMethod
                ]
            universal.run_cmd(cmd,outputdir)
            
        #-- To generate plot: expression of genes on UMAP plot
        #______________________________________________________________________________________________________
        if mode == 'ExpressionPlot':
            print(colored('==============> '+
                          'ExpressionPlot mode',
                          'green',
                          attrs=['bold'])
                 )
            myExpressionPlot = psource.resource_filename(__name__, "otherScripts/ExpressionPlot.R")
            
            inExpressionPlot = outputdir+'/'+'Cluster'+'/'+'sobj.rds'

            if not os.path.exists(inmyannoT):
                print(colored(inmyannoT + ': file does not exist',
                              'red', 
                              attrs=['bold'])
                     )
                exit()
            else:
                cmd=['Rscript', myExpressionPlot,
                     '--inputFile', inExpressionPlot,
                     '--outputDir', outputdir+'/'+'annoT',
                     '--ExpressionPlotGenes', ExpressionPlotGenes                     
                    ]
                universal.run_cmd(cmd,outputdir)
            
    
if __name__ == "__main__":
    main()

