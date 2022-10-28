# singlePipe
singleCell analysis pipeline for medical University of Vienna (version 1.0). Unpublished work.

Installation:

        git clone https://github.com/snizam001/singlePipe.git
        cd singlePipe
        pip install ./
    
Usage:

        Use this command after installation for usage instruction: 

        singlePipe --help

        optional arguments:
          -h, --help            show this help message and exit
          --nFeatureMin NFEATUREMIN
                                qc mode: Cells less than this feature will be filtered. 
                                Default value is 200
          --nFeatureMax NFEATUREMAX
                                qc mode: Cells more than this feature will be filtered. 
                                Default value is 2500
          --percMT PERCMT       qc mode: Cells more than this percentage of mitochondria will be filtered. 
                                Default value is 5 percent
          --normMethod NORMMETHOD
                                norm (normalization) mode: Normalization method: 
                                LogNormalize, CLR, RC. Default is LogNormalize.
          --normScalingFactor NORMSCALINGFACTOR
                                norm (normalization) mode: Scaling factor for the normalization. 
                                Default is 10000.
          --selMethod SELMETHOD
                                featureSelectDimentionReduction and cluster mode: selection method 
                                for variables. Options are vst, dispersion and
                                mean.var.plot. For detail see: 
                                https://satijalab.org/seurat/reference/findvariablefeatures. 
                                Default is vst.
          --VariableNumberFeatures VARIABLENUMBERFEATURES
                                featureSelectDimentionReduction and clustering mode: 
                                Number of the features/genes which are highly variable and used for
                                dimension reduction. Default is 2000.
          --dimPCA DIMPCA       clustering mode: How many eigenvectors or dimension should be used 
                                to find neighbors. Default is 40. This number should be
                                identified from the elbow plot which is generated in the 
                                featureSelectDimentionReduction mode
          --resPCA RESPCA       clustering mode: Value to find clusters. Default is 0.5. 
                                If --clustMethod is GranjaEtal, give three comma seperated values 
                                e.g. 0.2,0.6,1.0
          --clusterPlotN1 CLUSTERPLOTN1
                                clustering and annoT mode: Dimension number of 
                                UMAP/tsne/PCA to generate graph. Default is 1
          --clusterPlotN2 CLUSTERPLOTN2
                                clustering and annoT mode: Dimension number of 
                                UMAP/tsne/PCA to generate graph. Default is 2
          --clustMethod {Seurat,GranjaEtal,...}
                                clustering and diffExprGenes mode: Should this pipeline 
                                use default method (Seurat) or method of Granja et al (2019). Default
                                is Seurat.
          --findMarkersLogFC FINDMARKERSLOGFC
                                diffExprGenes mode: Limit testing to genes which show, 
                                on average, at least X-fold difference (log-scale) between the two
                                groups of cells. Default is 0.25.
          --findMarkersMinPCT FINDMARKERSMINPCT
                                diffExprGenes mode: Only test genes that are detected in 
                                a minimum fraction of min.pct cells in either of the two populations.
                                Default is 0.25.
          --ScTypeDB SCTYPEDB   annoT mode: Input excel file containing information of markers 
                                for cells of a particular tissue. By default it uses internal
                                datasset of package sc-type. But if you are using your own 
                                marker then just prepare an input XLSX file in the same format 
                                as our DB file. DB file should contain four columns (
                                tissueType - tissue type, cellName - cell type,
                                geneSymbolmore1 - positive marker genes, 
                                geneSymbolmore2 -marker genes not expected to be expressed by a 
                                cell type)
          --cellAnnotateTissue CELLANNOTATETISSUE
                                annoT mode: cells of this tissue should be used for annotation. 
                                Default is Immune system. By default choices are Adrenal,
                                Brain,Eye,Heart,Immune system,Intestine,Kidney,Liver,Lung,
                                Muscle,Pancreas,Placenta,Spleen,Stomach and Thymus. 
                                If you are using your own list
                                with --ScTypeDB, then define tissue type here.
          --ExpressionPlotGenes EXPRESSIONPLOTGENES
                                ExpressionPlot mode: Expression of gene on UMAP plot. 
                                Give the list of genes as comma seperated values

        required arguments
           _____________________:
          --modes MODES         run mode. Multiple modes can be provided as comma seperated values 
                                e.g. --modes mode1,mode2. Options are
                                qc,norm,featureSelectDimentionReduction,
                                clustering,diffExprGenes,annoT,ExpressionPlot,.....
          --outputdir OUTPUTDIR
                                output directory (give full path)

        common arguments
           _____________________:
          --inputData INPUTDATA
                                input data (give full path)
          --inputDataType {tsv,csv,rds,..}
                                input data type. Is it tsv, csv, rds, fastq or bamfiles. 
                                Default is rds.
          --threads THREADS     number of threads (default is: total CPU - 2)



Example:

        Download the example dataset from 
        https://jeffgranja.s3.amazonaws.com/MPAL-10x/Supplementary_Data/Healthy-Data/scRNA-Healthy-Hematopoiesis-191120.rds of 
        Granja et al. (Nat. Biotech 2019; Pubmed central ID = PMC7258684) and use following command:

        singlePipe 
        --modes qc,norm,featureSelectDimentionReduction,clustering,diffExprGenes,annoT,ExpressionPlot \
        --outputdir ~/Desktop/ViennaAssignment/out \
        --inputDataType rds \
        --threads 20 \
        --inputData /home/sheikh/Desktop/ViennaAssignment/data/scRNA-Healthy-Hematopoiesis-191120.rds \
        --ExpressionPlotGenes CD3D,CD14,CD19,CD8A,CEBPB,GATA1,TBX21,PAX5
    
    
    

