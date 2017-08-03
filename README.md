-   [Intalling rROMA](#intalling-rroma)
-   [Using rROMA](#using-rroma)
-   [An example](#an-example)
    -   [Getting the dataset](#getting-the-dataset)
    -   [Performing ROMA](#performing-roma)
    -   [Module activity](#module-activity)
    -   [Statistical comparison across
        samples](#statistical-comparison-across-samples)
    -   [Top contributing genes](#top-contributing-genes)
    -   [Gene weights](#gene-weights)
    -   [Sample projections](#sample-projections)
    -   [Recurrent genes](#recurrent-genes)
    -   [Exploring the dynamics of a selected
        gene](#exploring-the-dynamics-of-a-selected-gene)
    -   [Looking at the details](#looking-at-the-details)
    -   [Using the interactive
        dashboard](#using-the-interactive-dashboard)
    -   [Converting gene names](#converting-gene-names)
    -   [Visualising on ACSN](#visualising-on-acsn)
    -   [Session information](#session-information)
-   [Using the docker image](#using-the-docker-image)
-   [Using the docker container](#using-the-docker-container)

This package provides an R implementation of
[ROMA](http://journal.frontiersin.org/article/10.3389/fgene.2016.00018/full).
The package is under active development and is currently being tested.

A Java implementation developed by Andrei Zynovyev and is also
[available](https://github.com/sysbio-curie/Roma).

Intalling rROMA
===============

The rRoma package relies on the `scater` and `biomaRt` packages, which
are available only on BioConductor. These packagee can be installed with
the following command

    source("https://bioconductor.org/biocLite.R")

    ## Bioconductor version 3.5 (BiocInstaller 1.26.0), ?biocLite for help

    if(!requireNamespace("scater")){
      biocLite("scater")
    }

    ## Loading required namespace: scater

    if(!requireNamespace("biomaRt")){
      biocLite("biomaRt")
    }

rRoma can then be installed using `devtools`

    if(!requireNamespace("devtools")){
      install.packages("devtools")
    }

    if(!requireNamespace("rRoma")){
      devtools::install_github("Albluca/rROMA")
    }

To fill missing values rRoma uses the mice package. This package needs
to be installed manually if datasets with missing values need to be
analysed:

    if(!requireNamespace("tictoc")){
      install.packages("tictoc")
    }

Finally, rRoma allow projecting the results of the analysis on [ACSN
maps](https://acsn.curie.fr). To use this functionality it is necessary
to install the `RNaviCell` package:

    if(!requireNamespace("devtools")){
      install.packages("devtools")
    }

    if(!requireNamespace("RNaviCell")){
      devtools::install_github("sysbio-curie/RNaviCell")
    }

The packages `GEOquery`, `tictoc`, and `readr` are not required to run
rRoma, but are used in the following example and need to be installed to
reproduce the analysis

    if(!requireNamespace("GEOquery")){
      source("https://bioconductor.org/biocLite.R")
      biocLite("GEOquery")
    }

    ## Loading required namespace: GEOquery

    ## Setting options('download.file.method.GEOquery'='auto')

    ## Setting options('GEOquery.inmemory.gpl'=FALSE)

    if(!requireNamespace("readr")){
      install.packages("readr")
    }

    ## Loading required namespace: readr

    if(!requireNamespace("tictoc")){
      install.packages("tictoc")
    }

    ## Loading required namespace: tictoc

Using rROMA
===========

The package can be loaded with the usual syntax, i.e., by typing

    library(rRoma)

rRoma requires a gene expression matrix - with column names indicating
samples and row names indicating gene names - and a module file
containing information on the genesets that need to be evaluated. The
module file can be loaded from a GMT file. Functions to automate the
generation of then GMT file are also available.

Various functions are then available to explore the analysis, including
plotting and statistical cross sample analysis.

An example
==========

To show a concrete example of rROMA we will use a dataset available on
GEO.

Getting the dataset
-------------------

Let us begin by getting the description of the dataset

    library(GEOquery)

    ## Loading required package: Biobase

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, cbind, colMeans,
    ##     colnames, colSums, do.call, duplicated, eval, evalq, Filter,
    ##     Find, get, grep, grepl, intersect, is.unsorted, lapply,
    ##     lengths, Map, mapply, match, mget, order, paste, pmax,
    ##     pmax.int, pmin, pmin.int, Position, rank, rbind, Reduce,
    ##     rowMeans, rownames, rowSums, sapply, setdiff, sort, table,
    ##     tapply, union, unique, unsplit, which, which.max, which.min

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    gse <- getGEO("GSE50760", GSEMatrix = TRUE, destdir = getwd())

    ## https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50760/matrix/

    ## OK

    ## Found 1 file(s)

    ## GSE50760_series_matrix.txt.gz

    ## Using locally cached version: /bioinfo/users/lalberga/R/rRoma/GSE50760_series_matrix.txt.gz

    ## Using locally cached version of GPL11154 found here:
    ## /bioinfo/users/lalberga/R/rRoma/GPL11154.soft

Then we get the actual expression expression files

    filePaths = getGEOSuppFiles("GSE50760")

    ## https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50760/suppl/

    ## OK

Now we can construct the expression matrix. Note that the code below is
designed to work on a Unix-like environment (e.g. MacOS). The execution
on a Windows environment may require replacing `"/"` with `"\"`.

    library(readr)
    Content <- untar(row.names(filePaths)[1], list = TRUE)
    untar(row.names(filePaths)[1], list = FALSE)

    MatData <- NULL

    for(i in 1:length(Content)){
      Exp <- read_delim(Content[i], "\t", escape_double = FALSE, trim_ws = TRUE)
      
      if(is.null(MatData)){
        MatData <- cbind(unlist(Exp[,1]), unlist(Exp[,2]))
      } else {
        if(any(MatData[,1] != unlist(Exp[,1]))){
          stop("Incompatible samples")
        }
        MatData <- cbind(MatData, unlist(Exp[,2]))
      }
      
      file.remove(Content[i])
      
    }

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_2.1_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_3.1_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_5.1_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_6.1_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_7.1_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_8.1_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_9.1_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_10.1_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_12.1_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_13.1_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_17.1_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_18.1_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_19.1_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_20.1_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_21.1_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_22.1_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_23.1_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_24.1_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_2.2_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_3.2_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_5.2_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_6.2_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_7.2_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_8.2_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_9.2_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_10.2_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_12.2_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_13.2_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_17.2_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_18.2_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_19.2_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_20.2_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_21.2_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_22.2_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_23.2_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_24.2_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_2.3_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_3.3_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_5.3_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_6.3_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_7.3_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_8.3_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_9.3_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_10.3_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_12.3_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_13.3_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_17.3_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_18.3_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_19.3_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_20.3_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_21.3_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_22.3_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_23.3_FPKM = col_double()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character(),
    ##   AMC_24.3_FPKM = col_double()
    ## )

    SplitPath <- unlist(strsplit(rownames(filePaths)[1], "/"))
    SplitPath <- SplitPath[-length(SplitPath)]
    unlink(x = paste(SplitPath, collapse = "/"), recursive = TRUE)

    Genes <- MatData[,1]
    MatData <- data.matrix(data.frame(MatData[,-1]))
    rownames(MatData) <- Genes
    colnames(MatData) <- unlist(lapply(strsplit(Content, "_"), "[[", 1))

And look at the different groups of cells present

    Type <- as.character(gse$GSE50760_series_matrix.txt.gz[[1]])
    Type <- unlist(lapply(strsplit(Type, " "), "[[", 1))
    names(Type) = as.character(gse$GSE50760_series_matrix.txt.gz[[2]])

    table(Type)

    ## Type
    ## metastasized       normal      primary 
    ##           18           18           18

For convenience, we will transform `Type` into a factor:

    Type <- factor(Type, levels = c("normal", "primary", "metastasized"))

This will allow plotting functions, such as `Plot.Genesets`, to use more
meaningful ordering when reporting information.

### Selecting the module list

At this point we can create the metagene files. We will extract all the
"HALLMARK" geneset from MSig. We will also remove "HALLMARK\_" from the
names to simplify the graphical representation.

    AllHall <- SelectFromMSIGdb("HALLMARK")

    ## [1] "Searching in MsigDB v5.2"
    ## [1] "The following genesets have been selected:"
    ##  [1] "HALLMARK_TNFA_SIGNALING_VIA_NFKB (200 genes)"          
    ##  [2] "HALLMARK_HYPOXIA (200 genes)"                          
    ##  [3] "HALLMARK_CHOLESTEROL_HOMEOSTASIS (74 genes)"           
    ##  [4] "HALLMARK_MITOTIC_SPINDLE (200 genes)"                  
    ##  [5] "HALLMARK_WNT_BETA_CATENIN_SIGNALING (42 genes)"        
    ##  [6] "HALLMARK_TGF_BETA_SIGNALING (54 genes)"                
    ##  [7] "HALLMARK_IL6_JAK_STAT3_SIGNALING (87 genes)"           
    ##  [8] "HALLMARK_DNA_REPAIR (150 genes)"                       
    ##  [9] "HALLMARK_G2M_CHECKPOINT (200 genes)"                   
    ## [10] "HALLMARK_APOPTOSIS (161 genes)"                        
    ## [11] "HALLMARK_NOTCH_SIGNALING (32 genes)"                   
    ## [12] "HALLMARK_ADIPOGENESIS (200 genes)"                     
    ## [13] "HALLMARK_ESTROGEN_RESPONSE_EARLY (200 genes)"          
    ## [14] "HALLMARK_ESTROGEN_RESPONSE_LATE (200 genes)"           
    ## [15] "HALLMARK_ANDROGEN_RESPONSE (101 genes)"                
    ## [16] "HALLMARK_MYOGENESIS (200 genes)"                       
    ## [17] "HALLMARK_PROTEIN_SECRETION (96 genes)"                 
    ## [18] "HALLMARK_INTERFERON_ALPHA_RESPONSE (97 genes)"         
    ## [19] "HALLMARK_INTERFERON_GAMMA_RESPONSE (200 genes)"        
    ## [20] "HALLMARK_APICAL_JUNCTION (200 genes)"                  
    ## [21] "HALLMARK_APICAL_SURFACE (44 genes)"                    
    ## [22] "HALLMARK_HEDGEHOG_SIGNALING (36 genes)"                
    ## [23] "HALLMARK_COMPLEMENT (200 genes)"                       
    ## [24] "HALLMARK_UNFOLDED_PROTEIN_RESPONSE (113 genes)"        
    ## [25] "HALLMARK_PI3K_AKT_MTOR_SIGNALING (105 genes)"          
    ## [26] "HALLMARK_MTORC1_SIGNALING (200 genes)"                 
    ## [27] "HALLMARK_E2F_TARGETS (200 genes)"                      
    ## [28] "HALLMARK_MYC_TARGETS_V1 (200 genes)"                   
    ## [29] "HALLMARK_MYC_TARGETS_V2 (58 genes)"                    
    ## [30] "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION (200 genes)"
    ## [31] "HALLMARK_INFLAMMATORY_RESPONSE (200 genes)"            
    ## [32] "HALLMARK_XENOBIOTIC_METABOLISM (200 genes)"            
    ## [33] "HALLMARK_FATTY_ACID_METABOLISM (158 genes)"            
    ## [34] "HALLMARK_OXIDATIVE_PHOSPHORYLATION (200 genes)"        
    ## [35] "HALLMARK_GLYCOLYSIS (200 genes)"                       
    ## [36] "HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY (49 genes)"   
    ## [37] "HALLMARK_P53_PATHWAY (200 genes)"                      
    ## [38] "HALLMARK_UV_RESPONSE_UP (158 genes)"                   
    ## [39] "HALLMARK_UV_RESPONSE_DN (144 genes)"                   
    ## [40] "HALLMARK_ANGIOGENESIS (36 genes)"                      
    ## [41] "HALLMARK_HEME_METABOLISM (200 genes)"                  
    ## [42] "HALLMARK_COAGULATION (138 genes)"                      
    ## [43] "HALLMARK_IL2_STAT5_SIGNALING (200 genes)"              
    ## [44] "HALLMARK_BILE_ACID_METABOLISM (112 genes)"             
    ## [45] "HALLMARK_PEROXISOME (104 genes)"                       
    ## [46] "HALLMARK_ALLOGRAFT_REJECTION (200 genes)"              
    ## [47] "HALLMARK_SPERMATOGENESIS (135 genes)"                  
    ## [48] "HALLMARK_KRAS_SIGNALING_UP (200 genes)"                
    ## [49] "HALLMARK_KRAS_SIGNALING_DN (200 genes)"                
    ## [50] "HALLMARK_PANCREAS_BETA_CELLS (40 genes)"

    AllHall <- lapply(AllHall, function(x){
      x$Name <- sub("HALLMARK_", "", x$Name)
      x
    })

Performing ROMA
---------------

To reduce potential problems we will remove genes with a duplicated name

    if(any(duplicated(rownames(MatData)))){
      MatData <- MatData[!(rownames(MatData) %in% rownames(MatData)[duplicated(rownames(MatData))]), ]
    }

Moreover, we will use pseudocount transformation to boost the normality
of the data

    MatData <- log2(MatData + 1)

And now we are ready to perform ROMA by using the `rRoma.R` function.
The function has a number of parameters that can be used to control its
behaviour, but can be run with a minimal set of input: the expression
matrix (ExpressionMatrix) and the module list (ModuleList). In the
following example we will also used `FixedCenter` (which control if PCA
with fuixed center should be used), `UseParallel` (which control if the
computation should be done in parallel, analyzsis with fuixed center
should be used), `nCores` (the number of cores to use), and `ClusType`
(the type of parallel environment to use). We will also use a PC
orientation mode that will try to correlates gene expression level and
module score (`PCSignMode="CorrelateAllWeightsByGene"`). Other
parameters are described in the fucntion manual.

To perform ROMA fixed center in parallel, by using 8 cores, it is
sufficient write:

    tictoc::tic()
    Data.FC <- rRoma.R(ExpressionMatrix = MatData, ModuleList = AllHall, FixedCenter = TRUE,
                       UseParallel = TRUE, nCores = 8, ClusType = "FORK", MaxGenes = 100,
                       PCSignMode="CorrelateAllWeightsByGene", PCAType = "DimensionsAreSamples")

    ## [1] "Centering gene expression over samples (genes will have 0 mean)"
    ## [1] "Using global center (centering over genes)"
    ## [1] "The following geneset(s) will be ignored due to the number of genes being outside the specified range"
    ##  [1] "TNFA_SIGNALING_VIA_NFKB"          
    ##  [2] "HYPOXIA"                          
    ##  [3] "MITOTIC_SPINDLE"                  
    ##  [4] "DNA_REPAIR"                       
    ##  [5] "G2M_CHECKPOINT"                   
    ##  [6] "APOPTOSIS"                        
    ##  [7] "ADIPOGENESIS"                     
    ##  [8] "ESTROGEN_RESPONSE_EARLY"          
    ##  [9] "ESTROGEN_RESPONSE_LATE"           
    ## [10] "ANDROGEN_RESPONSE"                
    ## [11] "MYOGENESIS"                       
    ## [12] "INTERFERON_GAMMA_RESPONSE"        
    ## [13] "APICAL_JUNCTION"                  
    ## [14] "COMPLEMENT"                       
    ## [15] "UNFOLDED_PROTEIN_RESPONSE"        
    ## [16] "PI3K_AKT_MTOR_SIGNALING"          
    ## [17] "MTORC1_SIGNALING"                 
    ## [18] "E2F_TARGETS"                      
    ## [19] "MYC_TARGETS_V1"                   
    ## [20] "EPITHELIAL_MESENCHYMAL_TRANSITION"
    ## [21] "INFLAMMATORY_RESPONSE"            
    ## [22] "XENOBIOTIC_METABOLISM"            
    ## [23] "FATTY_ACID_METABOLISM"            
    ## [24] "OXIDATIVE_PHOSPHORYLATION"        
    ## [25] "GLYCOLYSIS"                       
    ## [26] "P53_PATHWAY"                      
    ## [27] "UV_RESPONSE_UP"                   
    ## [28] "UV_RESPONSE_DN"                   
    ## [29] "HEME_METABOLISM"                  
    ## [30] "COAGULATION"                      
    ## [31] "IL2_STAT5_SIGNALING"              
    ## [32] "BILE_ACID_METABOLISM"             
    ## [33] "PEROXISOME"                       
    ## [34] "ALLOGRAFT_REJECTION"              
    ## [35] "SPERMATOGENESIS"                  
    ## [36] "KRAS_SIGNALING_UP"                
    ## [37] "KRAS_SIGNALING_DN"                
    ## [1] "Too many cores selected! 7 will be used"
    ## [1] "2017-08-03 14:53:56 CEST"
    ## [1] "[1/13] Working on NOTCH_SIGNALING - http://www.broadinstitute.org/gsea/msigdb/cards/HALLMARK_NOTCH_SIGNALING"
    ## [1] "32 genes available for analysis"
    ## [1] "Detecting outliers using leave one out and median-absolute-deviations away from median (scater package)"
    ## [1] "1 gene(s) will be filtered:"
    ## [1] "WNT2"
    ## [1] "Not using weigths for PCA computation"
    ## [1] "Pre-filter data"
    ## [1] "L1 = 0.785374600359913 L1/L2 = 23.3872431178885"
    ## [1] "Median expression (uncentered): 13.5150221827314"
    ## [1] "Median expression (centered/weighted): 0.0288003550423865"
    ## [1] "Post-filter data"
    ## [1] "L1 = 0.217954829386733 L1/L2 = 3.96754668957679"
    ## [1] "Median expression (uncentered): 13.5291252684908"
    ## [1] "Median expression (centered/weighted): 0.0230072380009256"
    ## [1] "Previous sample size: 0"
    ## [1] "Next sample size: 32"
    ## [1] "Computing samples"
    ##    user  system elapsed 
    ##   0.094   0.170   2.509 
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "2017-08-03 14:54:00 CEST"
    ## [1] "[2/13] Working on HEDGEHOG_SIGNALING - http://www.broadinstitute.org/gsea/msigdb/cards/HALLMARK_HEDGEHOG_SIGNALING"
    ## [1] "36 genes available for analysis"
    ## [1] "Detecting outliers using leave one out and median-absolute-deviations away from median (scater package)"
    ## [1] "5 gene(s) will be filtered:"
    ## [1] "SCG2"   "SLIT1"  "PLG"    "NKX6-1" "CNTFR" 
    ## [1] "Not using weigths for PCA computation"
    ## [1] "Pre-filter data"
    ## [1] "L1 = 0.372953723336721 L1/L2 = 1.80271296604357"
    ## [1] "Median expression (uncentered): 13.1529183587471"
    ## [1] "Median expression (centered/weighted): 0.00223034662800391"
    ## [1] "Post-filter data"
    ## [1] "L1 = 0.259331522163147 L1/L2 = 2.57440442646472"
    ## [1] "Median expression (uncentered): 13.2823644791553"
    ## [1] "Median expression (centered/weighted): 0.0064586698739066"
    ## [1] "Previous sample size: 32"
    ## [1] "Next sample size: 36"
    ## [1] "Computing samples"
    ##    user  system elapsed 
    ##   0.089   0.173   2.889 
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "2017-08-03 14:54:04 CEST"
    ## [1] "[3/13] Working on ANGIOGENESIS - http://www.broadinstitute.org/gsea/msigdb/cards/HALLMARK_ANGIOGENESIS"
    ## [1] "36 genes available for analysis"
    ## [1] "Detecting outliers using leave one out and median-absolute-deviations away from median (scater package)"
    ## [1] "6 gene(s) will be filtered:"
    ## [1] "PF4"      "SERPINA5" "OLR1"     "PGLYRP1"  "PRG2"     "CXCL6"   
    ## [1] "Not using weigths for PCA computation"
    ## [1] "Pre-filter data"
    ## [1] "L1 = 0.349739414901111 L1/L2 = 1.37054382146703"
    ## [1] "Median expression (uncentered): 13.211812247159"
    ## [1] "Median expression (centered/weighted): 0.0245384923834588"
    ## [1] "Post-filter data"
    ## [1] "L1 = 0.231551014863786 L1/L2 = 1.34914572571376"
    ## [1] "Median expression (uncentered): 13.429864256795"
    ## [1] "Median expression (centered/weighted): -0.00680503589949048"
    ## [1] "Previous sample size: 36"
    ## [1] "Next sample size: 36"
    ## [1] "Reusing previous sampling (Same metagene size)"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "2017-08-03 14:54:04 CEST"
    ## [1] "[4/13] Working on PANCREAS_BETA_CELLS - http://www.broadinstitute.org/gsea/msigdb/cards/HALLMARK_PANCREAS_BETA_CELLS"
    ## [1] "40 genes available for analysis"
    ## [1] "Detecting outliers using leave one out and median-absolute-deviations away from median (scater package)"
    ## [1] "2 gene(s) will be filtered:"
    ## [1] "NEUROD1" "G6PC2"  
    ## [1] "Not using weigths for PCA computation"
    ## [1] "Pre-filter data"
    ## [1] "L1 = 0.350681958682222 L1/L2 = 3.08437236883014"
    ## [1] "Median expression (uncentered): 12.5489421629889"
    ## [1] "Median expression (centered/weighted): 0.109141826550011"
    ## [1] "Post-filter data"
    ## [1] "L1 = 0.346306300186471 L1/L2 = 2.77556448944163"
    ## [1] "Median expression (uncentered): 12.6889056658956"
    ## [1] "Median expression (centered/weighted): 0.113287760321446"
    ## [1] "Previous sample size: 36"
    ## [1] "Next sample size: 40"
    ## [1] "Computing samples"
    ##    user  system elapsed 
    ##   0.111   0.152   3.969 
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "2017-08-03 14:54:08 CEST"
    ## [1] "[5/13] Working on WNT_BETA_CATENIN_SIGNALING - http://www.broadinstitute.org/gsea/msigdb/cards/HALLMARK_WNT_BETA_CATENIN_SIGNALING"
    ## [1] "42 genes available for analysis"
    ## [1] "Detecting outliers using leave one out and median-absolute-deviations away from median (scater package)"
    ## [1] "4 gene(s) will be filtered:"
    ## [1] "WNT6" "DKK4" "DKK1" "WNT1"
    ## [1] "Not using weigths for PCA computation"
    ## [1] "Pre-filter data"
    ## [1] "L1 = 0.526845433876468 L1/L2 = 2.54306235780074"
    ## [1] "Median expression (uncentered): 13.3219279229047"
    ## [1] "Median expression (centered/weighted): 0.0252941445542472"
    ## [1] "Post-filter data"
    ## [1] "L1 = 0.256348671897308 L1/L2 = 7.07449466502326"
    ## [1] "Median expression (uncentered): 13.4160717745688"
    ## [1] "Median expression (centered/weighted): -0.00505501664700372"
    ## [1] "Previous sample size: 40"
    ## [1] "Next sample size: 42"
    ## [1] "Reusing previous sampling (Comparable metagene size)"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "2017-08-03 14:54:08 CEST"
    ## [1] "[6/13] Working on APICAL_SURFACE - http://www.broadinstitute.org/gsea/msigdb/cards/HALLMARK_APICAL_SURFACE"
    ## [1] "44 genes available for analysis"
    ## [1] "Detecting outliers using leave one out and median-absolute-deviations away from median (scater package)"
    ## [1] "10 gene(s) will be filtered:"
    ##  [1] "RHCG"     "MAL"      "PKHD1"    "ATP6V0A4" "GHRL"     "SLC34A3" 
    ##  [7] "RTN4RL1"  "CD160"    "SLC22A12" "NTNG1"   
    ## [1] "Not using weigths for PCA computation"
    ## [1] "Pre-filter data"
    ## [1] "L1 = 0.299561449815429 L1/L2 = 1.64452196865508"
    ## [1] "Median expression (uncentered): 13.078317612573"
    ## [1] "Median expression (centered/weighted): 0.0615200797577556"
    ## [1] "Post-filter data"
    ## [1] "L1 = 0.332285744487647 L1/L2 = 3.28815805594463"
    ## [1] "Median expression (uncentered): 13.4014128725526"
    ## [1] "Median expression (centered/weighted): 0.0245102665025946"
    ## [1] "Previous sample size: 40"
    ## [1] "Next sample size: 44"
    ## [1] "Computing samples"
    ##    user  system elapsed 
    ##   0.106   0.156   3.308 
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "2017-08-03 14:54:12 CEST"
    ## [1] "[7/13] Working on REACTIVE_OXIGEN_SPECIES_PATHWAY - http://www.broadinstitute.org/gsea/msigdb/cards/HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY"
    ## [1] "48 genes available for analysis"
    ## [1] "Detecting outliers using leave one out and median-absolute-deviations away from median (scater package)"
    ## [1] "1 gene(s) will be filtered:"
    ## [1] "MPO"
    ## [1] "Not using weigths for PCA computation"
    ## [1] "Pre-filter data"
    ## [1] "L1 = 0.573232544690332 L1/L2 = 15.0864407341184"
    ## [1] "Median expression (uncentered): 13.6307790036278"
    ## [1] "Median expression (centered/weighted): 0.0204939834579211"
    ## [1] "Post-filter data"
    ## [1] "L1 = 0.149879220496982 L1/L2 = 2.70307339280201"
    ## [1] "Median expression (uncentered): 13.6419952649137"
    ## [1] "Median expression (centered/weighted): 0.0139898251723049"
    ## [1] "Previous sample size: 44"
    ## [1] "Next sample size: 48"
    ## [1] "Computing samples"
    ##    user  system elapsed 
    ##   0.111   0.158   3.495 
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "2017-08-03 14:54:15 CEST"
    ## [1] "[8/13] Working on TGF_BETA_SIGNALING - http://www.broadinstitute.org/gsea/msigdb/cards/HALLMARK_TGF_BETA_SIGNALING"
    ## [1] "54 genes available for analysis"
    ## [1] "Detecting outliers using leave one out and median-absolute-deviations away from median (scater package)"
    ## [1] "2 gene(s) will be filtered:"
    ## [1] "LEFTY2" "NOG"   
    ## [1] "Not using weigths for PCA computation"
    ## [1] "Pre-filter data"
    ## [1] "L1 = 0.56523572781357 L1/L2 = 6.28863174085824"
    ## [1] "Median expression (uncentered): 13.6414872140507"
    ## [1] "Median expression (centered/weighted): 0.0335821943623831"
    ## [1] "Post-filter data"
    ## [1] "L1 = 0.147813066834024 L1/L2 = 2.44916461749659"
    ## [1] "Median expression (uncentered): 13.6720935948085"
    ## [1] "Median expression (centered/weighted): 0.0286294392766931"
    ## [1] "Previous sample size: 48"
    ## [1] "Next sample size: 54"
    ## [1] "Computing samples"
    ##    user  system elapsed 
    ##   0.125   0.163   4.052 
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "2017-08-03 14:54:20 CEST"
    ## [1] "[9/13] Working on MYC_TARGETS_V2 - http://www.broadinstitute.org/gsea/msigdb/cards/HALLMARK_MYC_TARGETS_V2"
    ## [1] "58 genes available for analysis"
    ## [1] "Detecting outliers using leave one out and median-absolute-deviations away from median (scater package)"
    ## [1] "No gene will be filtered"
    ## [1] "Not using weigths for PCA computation"
    ## [1] "Pre-filter data"
    ## [1] "L1 = 0.30711252440909 L1/L2 = 5.75080346266076"
    ## [1] "Median expression (uncentered): 13.5235005936785"
    ## [1] "Median expression (centered/weighted): 0.0355963366138567"
    ## [1] "Post-filter data"
    ## [1] "L1 = 0.30711252440909 L1/L2 = 5.75080346266791"
    ## [1] "Median expression (uncentered): 13.5235005936785"
    ## [1] "Median expression (centered/weighted): 0.0355963366138567"
    ## [1] "Previous sample size: 54"
    ## [1] "Next sample size: 58"
    ## [1] "Computing samples"
    ##    user  system elapsed 
    ##   0.134   0.157   4.934 
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "2017-08-03 14:54:25 CEST"
    ## [1] "[10/13] Working on CHOLESTEROL_HOMEOSTASIS - http://www.broadinstitute.org/gsea/msigdb/cards/HALLMARK_CHOLESTEROL_HOMEOSTASIS"
    ## [1] "74 genes available for analysis"
    ## [1] "Detecting outliers using leave one out and median-absolute-deviations away from median (scater package)"
    ## [1] "2 gene(s) will be filtered:"
    ## [1] "ADH4"   "AVPR1A"
    ## [1] "Not using weigths for PCA computation"
    ## [1] "Pre-filter data"
    ## [1] "L1 = 0.184778656708759 L1/L2 = 1.72337310774601"
    ## [1] "Median expression (uncentered): 13.5770747151966"
    ## [1] "Median expression (centered/weighted): 0.0261952514023445"
    ## [1] "Post-filter data"
    ## [1] "L1 = 0.0827272783459245 L1/L2 = 0.715041263902501"
    ## [1] "Median expression (uncentered): 13.6116012020924"
    ## [1] "Median expression (centered/weighted): 0.0278583738048012"
    ## [1] "Previous sample size: 58"
    ## [1] "Next sample size: 74"
    ## [1] "Computing samples"
    ##    user  system elapsed 
    ##   0.135   0.159   5.511 
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "2017-08-03 14:54:31 CEST"
    ## [1] "[11/13] Working on IL6_JAK_STAT3_SIGNALING - http://www.broadinstitute.org/gsea/msigdb/cards/HALLMARK_IL6_JAK_STAT3_SIGNALING"
    ## [1] "87 genes available for analysis"
    ## [1] "Detecting outliers using leave one out and median-absolute-deviations away from median (scater package)"
    ## [1] "9 gene(s) will be filtered:"
    ## [1] "IL6"   "REG1A" "INHBE" "CRLF2" "PF4"   "DNTT"  "CSF2"  "CNTFR" "CCL7" 
    ## [1] "Not using weigths for PCA computation"
    ## [1] "Pre-filter data"
    ## [1] "L1 = 0.260911836642868 L1/L2 = 1.20988202150609"
    ## [1] "Median expression (uncentered): 13.4361909954814"
    ## [1] "Median expression (centered/weighted): 0.0218640317615592"
    ## [1] "Post-filter data"
    ## [1] "L1 = 0.134356840056698 L1/L2 = 0.937336409812835"
    ## [1] "Median expression (uncentered): 13.4997835811167"
    ## [1] "Median expression (centered/weighted): 0.011905293892675"
    ## [1] "Previous sample size: 74"
    ## [1] "Next sample size: 87"
    ## [1] "Computing samples"
    ##    user  system elapsed 
    ##   0.141   0.170   6.375 
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "2017-08-03 14:54:37 CEST"
    ## [1] "[12/13] Working on PROTEIN_SECRETION - http://www.broadinstitute.org/gsea/msigdb/cards/HALLMARK_PROTEIN_SECRETION"
    ## [1] "96 genes available for analysis"
    ## [1] "Detecting outliers using leave one out and median-absolute-deviations away from median (scater package)"
    ## [1] "2 gene(s) will be filtered:"
    ## [1] "SH3GL2"   "ATP6V1B1"
    ## [1] "Not using weigths for PCA computation"
    ## [1] "Pre-filter data"
    ## [1] "L1 = 0.473571242899042 L1/L2 = 2.52388526667304"
    ## [1] "Median expression (uncentered): 13.52851515655"
    ## [1] "Median expression (centered/weighted): -0.0081114967003805"
    ## [1] "Post-filter data"
    ## [1] "L1 = 0.0713155467495194 L1/L2 = 0.565493800424161"
    ## [1] "Median expression (uncentered): 13.5401280182475"
    ## [1] "Median expression (centered/weighted): -0.0114954339367517"
    ## [1] "Previous sample size: 87"
    ## [1] "Next sample size: 96"
    ## [1] "Computing samples"
    ##    user  system elapsed 
    ##   0.145   0.175   7.207 
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "2017-08-03 14:54:45 CEST"
    ## [1] "[13/13] Working on INTERFERON_ALPHA_RESPONSE - http://www.broadinstitute.org/gsea/msigdb/cards/HALLMARK_INTERFERON_ALPHA_RESPONSE"
    ## [1] "97 genes available for analysis"
    ## [1] "Detecting outliers using leave one out and median-absolute-deviations away from median (scater package)"
    ## [1] "3 gene(s) will be filtered:"
    ## [1] "SAMD9"  "SAMD9L" "IL7"   
    ## [1] "Not using weigths for PCA computation"
    ## [1] "Pre-filter data"
    ## [1] "L1 = 0.1322894563459 L1/L2 = 1.08589917483378"
    ## [1] "Median expression (uncentered): 13.5266826196238"
    ## [1] "Median expression (centered/weighted): 0.0263565592174932"
    ## [1] "Post-filter data"
    ## [1] "L1 = 0.125731998097361 L1/L2 = 1.5321298006571"
    ## [1] "Median expression (uncentered): 13.52851515655"
    ## [1] "Median expression (centered/weighted): 0.0248646010644798"
    ## [1] "Previous sample size: 96"
    ## [1] "Next sample size: 97"
    ## [1] "Reusing previous sampling (Comparable metagene size)"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"

    tictoc::toc()

    ## 54.631 sec elapsed

Module activity
---------------

We can now explore the overdispersed genesets in ROMA with fixed center.
Since we have information on the groups, we can also aggregate the
module weigths by groups and consider the mean and standard deviation.
Note how we used the `SelectGeneSets` function to determind the
overdispersed genesets.

    AggData.FC <- Plot.Genesets(RomaData = Data.FC,
                  Selected = SelectGeneSets(RomaData = Data.FC, VarThr = 5e-2,
                                            VarMode = "Wil", VarType = "Over"),
                  GenesetMargin = 20, SampleMargin = 14, cluster_cols = FALSE,
                  GroupInfo = Type, AggByGroupsFL = c("mean", "sd"),
                  HMTite = "Overdispersed genesets (Fixed center)")

    ## [1] "Using genestes overdispersed according to Wilcoxon test. VarThr = 0.05"
    ## [1] "No coordinatedness filter selected"
    ## [1] "No expression filter selected"
    ## [1] "13 geneset selected"

![](README_files/figure-markdown_strict/unnamed-chunk-16-1.png)![](README_files/figure-markdown_strict/unnamed-chunk-16-2.png)

    ## [1] "54 samples have an associated group"

![](README_files/figure-markdown_strict/unnamed-chunk-16-3.png)![](README_files/figure-markdown_strict/unnamed-chunk-16-4.png)![](README_files/figure-markdown_strict/unnamed-chunk-16-5.png)

The fucntion returns a list of matrices with the aggregated data, which
we saved in the `AggData.FC` variable.

    AggData.FC

    ## $mean
    ##                                      normal      primary metastasized
    ## CHOLESTEROL_HOMEOSTASIS          0.08424088  0.015763817  -0.10000469
    ## WNT_BETA_CATENIN_SIGNALING       0.14079850 -0.042992332  -0.09780617
    ## TGF_BETA_SIGNALING              -0.13305532  0.035013898   0.09804142
    ## IL6_JAK_STAT3_SIGNALING         -0.10286116  0.023109146   0.07975202
    ## NOTCH_SIGNALING                  0.12418462 -0.002713360  -0.12147126
    ## PROTEIN_SECRETION                0.05055215  0.004004691  -0.05455684
    ## INTERFERON_ALPHA_RESPONSE       -0.09395742  0.054325536   0.03963188
    ## APICAL_SURFACE                   0.13922203 -0.015341109  -0.12388093
    ## HEDGEHOG_SIGNALING               0.14012739 -0.042143382  -0.09798401
    ## MYC_TARGETS_V2                  -0.14113970  0.045970991   0.09516871
    ## REACTIVE_OXIGEN_SPECIES_PATHWAY  0.10889324 -0.033756284  -0.07513696
    ## ANGIOGENESIS                    -0.13936076  0.011379056   0.12798170
    ## PANCREAS_BETA_CELLS              0.12323313  0.013303027  -0.13653615
    ## 
    ## $sd
    ##                                     normal    primary metastasized
    ## CHOLESTEROL_HOMEOSTASIS         0.07033183 0.11182778   0.15167465
    ## WNT_BETA_CATENIN_SIGNALING      0.06902014 0.10448563   0.10033111
    ## TGF_BETA_SIGNALING              0.07054766 0.10630552   0.11101780
    ## IL6_JAK_STAT3_SIGNALING         0.10705038 0.12262768   0.11757308
    ## NOTCH_SIGNALING                 0.09765782 0.07189569   0.11026098
    ## PROTEIN_SECRETION               0.07908624 0.10602382   0.18829096
    ## INTERFERON_ALPHA_RESPONSE       0.12606849 0.10546781   0.13293445
    ## APICAL_SURFACE                  0.05105022 0.09204174   0.10355855
    ## HEDGEHOG_SIGNALING              0.04728758 0.11403235   0.10366845
    ## MYC_TARGETS_V2                  0.04035685 0.13110017   0.08418963
    ## REACTIVE_OXIGEN_SPECIES_PATHWAY 0.09654433 0.12260551   0.12137211
    ## ANGIOGENESIS                    0.06220049 0.08366833   0.09955138
    ## PANCREAS_BETA_CELLS             0.02589543 0.10425537   0.10619772

We can also look at underdispersed datasets. In this case we will also
look at the expression level and select only the genesets that are
significantly overexpressed (`MedType = "Over"`).

    AggData2.FC <- Plot.Genesets(RomaData = Data.FC,
                  Selected = SelectGeneSets(RomaData = Data.FC, VarThr = 5e-2,
                                            VarMode = "Wil", VarType = "Under",
                                            MedThr = 5e-3, MedMode = "Wil", MedType = "Over"),
                  GenesetMargin = 20, SampleMargin = 14, cluster_cols = FALSE,
                  GroupInfo = Type, AggByGroupsFL = c("mean", "sd"),
                  HMTite = "Underdispersed genesets (Fixed center)")

    ## [1] "Using genestes underdispersed according to Wilcoxon test. VarThr = 0.05"
    ## [1] "No coordinatedness filter selected"
    ## [1] "Using genestes overexpressed according to Wilcoxon test. MedThr = 0.005"
    ## [1] "10 geneset selected"

![](README_files/figure-markdown_strict/unnamed-chunk-18-1.png)![](README_files/figure-markdown_strict/unnamed-chunk-18-2.png)

    ## [1] "54 samples have an associated group"

![](README_files/figure-markdown_strict/unnamed-chunk-18-3.png)![](README_files/figure-markdown_strict/unnamed-chunk-18-4.png)![](README_files/figure-markdown_strict/unnamed-chunk-18-5.png)

Statistical comparison across samples
-------------------------------------

Visual inspection already stresses the difference between the groups. It
is possible to quantify theses differnces performing statistical
testing:

    CompareAcrossSamples(RomaData = Data.FC,
                         Selected = SelectGeneSets(RomaData = Data.FC, VarThr = 5e-2,
                                                   VarMode = "Wil", VarType = "Over"),
                         Groups = Type)

    ## [1] "Using genestes overdispersed according to Wilcoxon test. VarThr = 0.05"
    ## [1] "No coordinatedness filter selected"
    ## [1] "No expression filter selected"
    ## [1] "13 geneset selected"
    ## [1] "Performing Type III AOV (R default)"
    ##              Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Group         2  0.318 0.15883   8.754 0.000176 ***
    ## Residuals   699 12.682 0.01814                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

![](README_files/figure-markdown_strict/unnamed-chunk-19-1.png)![](README_files/figure-markdown_strict/unnamed-chunk-19-2.png)

    ## [1] "A significant difference is observed across groups"
    ## [1] "Calculating Tukey Honest Significant Differences"
    ## [1] "2 significant differences found"

![](README_files/figure-markdown_strict/unnamed-chunk-19-3.png)

    ## [1] "Performing Type III AOV (R default)"
    ##                Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Group           2  0.318 0.15883   14.89 4.73e-07 ***
    ## Group:GeneSet  36  5.610 0.15584   14.61  < 2e-16 ***
    ## Residuals     663  7.072 0.01067                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

![](README_files/figure-markdown_strict/unnamed-chunk-19-4.png)![](README_files/figure-markdown_strict/unnamed-chunk-19-5.png)![](README_files/figure-markdown_strict/unnamed-chunk-19-6.png)![](README_files/figure-markdown_strict/unnamed-chunk-19-7.png)

    ## [1] "A significant difference is observed across groups and metagenes"
    ## [1] "Calculating Tukey Honest Significant Differences"
    ## [1] "260 significant differences found"

![](README_files/figure-markdown_strict/unnamed-chunk-19-8.png)![](README_files/figure-markdown_strict/unnamed-chunk-19-9.png)![](README_files/figure-markdown_strict/unnamed-chunk-19-10.png)![](README_files/figure-markdown_strict/unnamed-chunk-19-11.png)![](README_files/figure-markdown_strict/unnamed-chunk-19-12.png)![](README_files/figure-markdown_strict/unnamed-chunk-19-13.png)![](README_files/figure-markdown_strict/unnamed-chunk-19-14.png)![](README_files/figure-markdown_strict/unnamed-chunk-19-15.png)![](README_files/figure-markdown_strict/unnamed-chunk-19-16.png)![](README_files/figure-markdown_strict/unnamed-chunk-19-17.png)![](README_files/figure-markdown_strict/unnamed-chunk-19-18.png)![](README_files/figure-markdown_strict/unnamed-chunk-19-19.png)

Top contributing genes
----------------------

We can also explore the top contributing genes, i.e., the genes with the
largers weights in the selected modules:

    GeneMat <- GetTopContrib(Data.FC, Selected = SelectGeneSets(RomaData = Data.FC, VarThr = 1e-5,
                                            VarMode = "Wil", VarType = "Over"), nGenes = 6, OrderType = "Abs")

    ## [1] "Using genestes overdispersed according to Wilcoxon test. VarThr = 1e-05"
    ## [1] "No coordinatedness filter selected"
    ## [1] "No expression filter selected"

![](README_files/figure-markdown_strict/unnamed-chunk-20-1.png)

This allow us to quicky determine if a set of genes play a strong
(potentially opposite) role across different genesets. The function also
returns the matrix used to plot the heatmap.

Gene weights
------------

We can also explore the gene weigths with the fixed center

    PlotGeneWeight(RomaData = Data.FC, PlotGenes = 30,
                   ExpressionMatrix = MatData, LogExpression = FALSE,
                   Selected = SelectGeneSets(RomaData = Data.FC, VarThr = 5e-2,
                                             VarMode = "Wil", VarType = "Over"),
                   PlotWeigthSign = TRUE)

    ## [1] "Using genestes overdispersed according to Wilcoxon test. VarThr = 0.05"
    ## [1] "No coordinatedness filter selected"
    ## [1] "No expression filter selected"
    ## [1] "13 geneset selected"

![](README_files/figure-markdown_strict/unnamed-chunk-21-1.png)![](README_files/figure-markdown_strict/unnamed-chunk-21-2.png)![](README_files/figure-markdown_strict/unnamed-chunk-21-3.png)![](README_files/figure-markdown_strict/unnamed-chunk-21-4.png)![](README_files/figure-markdown_strict/unnamed-chunk-21-5.png)![](README_files/figure-markdown_strict/unnamed-chunk-21-6.png)![](README_files/figure-markdown_strict/unnamed-chunk-21-7.png)![](README_files/figure-markdown_strict/unnamed-chunk-21-8.png)![](README_files/figure-markdown_strict/unnamed-chunk-21-9.png)![](README_files/figure-markdown_strict/unnamed-chunk-21-10.png)![](README_files/figure-markdown_strict/unnamed-chunk-21-11.png)![](README_files/figure-markdown_strict/unnamed-chunk-21-12.png)![](README_files/figure-markdown_strict/unnamed-chunk-21-13.png)

Sample projections
------------------

Moreover, we can look at the projections of the samples with the fixed
center

    PlotSampleProjections(RomaData = Data.FC, PlotSamples = 30,
                          ExpressionMatrix = MatData, LogExpression = FALSE,
                          Selected = SelectGeneSets(RomaData = Data.FC, VarThr = 5e-3,
                                                    VarMode = "Wil", VarType = "Over"))

    ## [1] "Using genestes overdispersed according to Wilcoxon test. VarThr = 0.005"
    ## [1] "No coordinatedness filter selected"
    ## [1] "No expression filter selected"
    ## [1] "13 geneset selected"

![](README_files/figure-markdown_strict/unnamed-chunk-22-1.png)![](README_files/figure-markdown_strict/unnamed-chunk-22-2.png)![](README_files/figure-markdown_strict/unnamed-chunk-22-3.png)![](README_files/figure-markdown_strict/unnamed-chunk-22-4.png)![](README_files/figure-markdown_strict/unnamed-chunk-22-5.png)![](README_files/figure-markdown_strict/unnamed-chunk-22-6.png)![](README_files/figure-markdown_strict/unnamed-chunk-22-7.png)![](README_files/figure-markdown_strict/unnamed-chunk-22-8.png)![](README_files/figure-markdown_strict/unnamed-chunk-22-9.png)![](README_files/figure-markdown_strict/unnamed-chunk-22-10.png)![](README_files/figure-markdown_strict/unnamed-chunk-22-11.png)![](README_files/figure-markdown_strict/unnamed-chunk-22-12.png)![](README_files/figure-markdown_strict/unnamed-chunk-22-13.png)

Recurrent genes
---------------

Finally, we can look at genes which appear across different genesests
and explore thier weithgs

    PlotRecurringGenes(RomaData = Data.FC,
                       Selected = SelectGeneSets(RomaData = Data.FC, VarThr = 5e-3,
                                                 VarMode = "Wil", VarType = "Over"),
                       GenesByGroup = 25, MinMult = 2)

    ## [1] "Using genestes overdispersed according to Wilcoxon test. VarThr = 0.005"
    ## [1] "No coordinatedness filter selected"
    ## [1] "No expression filter selected"
    ## [1] "13 geneset selected"

![](README_files/figure-markdown_strict/unnamed-chunk-23-1.png)![](README_files/figure-markdown_strict/unnamed-chunk-23-2.png)![](README_files/figure-markdown_strict/unnamed-chunk-23-3.png)![](README_files/figure-markdown_strict/unnamed-chunk-23-4.png)![](README_files/figure-markdown_strict/unnamed-chunk-23-5.png)![](README_files/figure-markdown_strict/unnamed-chunk-23-6.png)![](README_files/figure-markdown_strict/unnamed-chunk-23-7.png)

Exploring the dynamics of a selected gene
-----------------------------------------

If we are interested in the dynamics of a specific gene, it is possible
to look at it using the function `ExploreGeneProperties`. This function
will repost a number of plots that provide information on the working of
a gene. Note the `DO_tSNE`, which will produce a tSNE plot when set to
`TRUE`, and the `initial_dims` and `perplexity` parameters that control
the tSNE plot.

    ExploreGeneProperties(RomaData = Data.FC, 
                          Selected = SelectGeneSets(RomaData = Data.FC, VarThr = 5e-3,
                                                 VarMode = "Wil", VarType = "Over"),
                          GeneName = "JAG1", ExpressionMatrix = MatData, GroupInfo = Type,
                          DO_tSNE = TRUE, initial_dims = 50, perplexity = 15)

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning: Removed 11466 rows containing non-finite values (stat_density).

![](README_files/figure-markdown_strict/unnamed-chunk-24-1.png)![](README_files/figure-markdown_strict/unnamed-chunk-24-2.png)![](README_files/figure-markdown_strict/unnamed-chunk-24-3.png)![](README_files/figure-markdown_strict/unnamed-chunk-24-4.png)![](README_files/figure-markdown_strict/unnamed-chunk-24-5.png)

    ## [1] "Using genestes overdispersed according to Wilcoxon test. VarThr = 0.005"
    ## [1] "No coordinatedness filter selected"
    ## [1] "No expression filter selected"
    ## [1] "Gene found in 4 modules"

![](README_files/figure-markdown_strict/unnamed-chunk-24-6.png)

    ## `geom_smooth()` using method = 'loess'

![](README_files/figure-markdown_strict/unnamed-chunk-24-7.png)

    ExploreGeneProperties(RomaData = Data.FC, 
                          Selected = SelectGeneSets(RomaData = Data.FC, VarThr = 5e-3,
                                                 VarMode = "Wil", VarType = "Over"),
                          GeneName = "NEUROG3", ExpressionMatrix = MatData, GroupInfo = Type,
                          DO_tSNE = TRUE, initial_dims = 50, perplexity = 15)

    ## Warning: Transformation introduced infinite values in continuous x-axis

    ## Warning: Removed 11466 rows containing non-finite values (stat_density).

![](README_files/figure-markdown_strict/unnamed-chunk-25-1.png)

    ## Warning in wilcox.test.default(xi, xj, paired = paired, ...): cannot
    ## compute exact p-value with ties

    ## Warning in wilcox.test.default(xi, xj, paired = paired, ...): cannot
    ## compute exact p-value with ties

    ## Warning in wilcox.test.default(xi, xj, paired = paired, ...): cannot
    ## compute exact p-value with ties

![](README_files/figure-markdown_strict/unnamed-chunk-25-2.png)

    ## Warning in wilcox.test.default(c(12.0610336016811, 12.5423064225826,
    ## 11.6302671250268, : cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(10.0251395622785, 13.9943534368589,
    ## 11.2871351860867, : cannot compute exact p-value with ties

![](README_files/figure-markdown_strict/unnamed-chunk-25-3.png)![](README_files/figure-markdown_strict/unnamed-chunk-25-4.png)![](README_files/figure-markdown_strict/unnamed-chunk-25-5.png)

    ## [1] "Using genestes overdispersed according to Wilcoxon test. VarThr = 0.005"
    ## [1] "No coordinatedness filter selected"
    ## [1] "No expression filter selected"
    ## [1] "Gene found in 1 modules"

![](README_files/figure-markdown_strict/unnamed-chunk-25-6.png)

    ## `geom_smooth()` using method = 'loess'

![](README_files/figure-markdown_strict/unnamed-chunk-25-7.png)

Looking at the details
----------------------

By default, rRoma will perform the appropiate anaysis without showing
all the details to the users. However, it is possible to obtain
additioanl graphical and textual information. This will results in a
large amount of inromation being derived an plotted. Hence, it is not
advisable to do that for large analysis. To show an example of the
extended information that can be produced, we will first obtain a
smaller module list.

    RedGMT <- SelectFromMSIGdb(SearchString = c("xenobiotic", "keg"), Mode = "ALL")

    ## [1] "Searching in MsigDB v5.2"
    ## [1] "The following genesets have been selected:"
    ## [1] "KEGG_METABOLISM_OF_XENOBIOTICS_BY_CYTOCHROME_P450 (70 genes)"

Wen perfoming ROMA we will now set `PlotData = TRUE` (to produce
diagnostic plots), `MoreInfo = TRUE` (to print additional diagnostic
information), and `FullSampleInfo = TRUE` (to compute and save
additional information on the samples). Note that rRoma will ask to
confim cetain choiches if R is run interactivelly.

    tictoc::tic()
    RedData.NFC <- rRoma.R(ExpressionMatrix = MatData, ModuleList = RedGMT, FixedCenter = FALSE,
                        UseParallel = TRUE, nCores = 8, ClusType = "FORK",
                        PCSignMode="CorrelateAllWeightsByGene", 
                        PlotData = TRUE, MoreInfo = TRUE, FullSampleInfo = TRUE,
                        Grouping = Type, GroupPCSign = FALSE, PCAType = "DimensionsAreGenes")

    ## [1] "Centering gene expression over samples (genes will have 0 mean)"
    ## [1] "Using local center (NOT centering over genes)"
    ## [1] "All the genesets will be used"
    ## [1] "Too many cores selected! 7 will be used"
    ## [1] "2017-08-03 14:56:12 CEST"
    ## [1] "[1/1] Working on KEGG_METABOLISM_OF_XENOBIOTICS_BY_CYTOCHROME_P450 - http://www.broadinstitute.org/gsea/msigdb/cards/KEGG_METABOLISM_OF_XENOBIOTICS_BY_CYTOCHROME_P450"
    ## [1] "70 genes available for analysis"
    ## [1] "The following genes will be used:"
    ##  [1] "CYP2F1"  "CYP2E1"  "AKR1C4"  "EPHX1"   "CYP3A5"  "UGT2B28" "CYP3A4" 
    ##  [8] "GSTP1"   "GSTT2"   "GSTT1"   "AKR1C3"  "GSTZ1"   "CYP2C18" "ADH1B"  
    ## [15] "CYP2S1"  "ADH1C"   "ADH4"    "ADH5"    "GSTO1"   "ADH1A"   "GSTA5"  
    ## [22] "MGST2"   "MGST1"   "MGST3"   "GSTA3"   "GSTM1"   "UGT2B11" "CYP2C9" 
    ## [29] "GSTA4"   "CYP2C19" "GSTM4"   "CYP2C8"  "GSTM3"   "CYP2B6"  "GSTM2"  
    ## [36] "UGT2A3"  "UGT1A4"  "UGT1A1"  "CYP3A7"  "UGT1A3"  "GSTM5"   "UGT1A10"
    ## [43] "UGT1A8"  "UGT1A7"  "GSTA1"   "UGT1A6"  "GSTA2"   "ALDH3A1" "UGT1A5" 
    ## [50] "GSTK1"   "AKR1C2"  "AKR1C1"  "CYP1A1"  "ADH7"    "CYP1A2"  "CYP1B1" 
    ## [57] "ADH6"    "UGT2A1"  "ALDH1A3" "ALDH3B1" "ALDH3B2" "CYP3A43" "UGT1A9" 
    ## [64] "DHDH"    "GSTO2"   "UGT2B10" "UGT2B7"  "UGT2B4"  "UGT2B17" "UGT2B15"
    ## [1] "Detecting outliers using leave one out and median-absolute-deviations away from median (scater package)"

![](README_files/figure-markdown_strict/unnamed-chunk-27-1.png)

    ## [1] "19 gene(s) will be filtered:"
    ##  [1] "CYP2F1"  "AKR1C4"  "UGT2B28" "GSTT1"   "UGT2B11" "UGT1A4"  "UGT1A3" 
    ##  [8] "UGT1A8"  "UGT1A7"  "UGT1A5"  "CYP1A1"  "CYP1A2"  "UGT2A1"  "ALDH3B2"
    ## [15] "CYP3A43" "UGT1A9"  "DHDH"    "UGT2B4"  "UGT2B17"
    ## [1] "Not using weigths for PCA computation"
    ## [1] "Pre-filter data"
    ## [1] "L1 = 0.32484762172228 L1/L2 = 3.36077904254904"
    ## [1] "Median expression (uncentered): 12.7811547706959"
    ## [1] "Median expression (centered/weighted): 0.10923220260671"
    ## [1] "Post-filter data"
    ## [1] "L1 = 0.391959147712557 L1/L2 = 3.06919422461758"
    ## [1] "Median expression (uncentered): 13.2177305375661"
    ## [1] "Median expression (centered/weighted): 0.0819972157791584"
    ## [1] "Previous sample size: 0"
    ## [1] "Next sample size: 70"
    ## [1] "Computing samples"
    ##    user  system elapsed 
    ##   0.115   0.160   8.764

![](README_files/figure-markdown_strict/unnamed-chunk-27-2.png)![](README_files/figure-markdown_strict/unnamed-chunk-27-3.png)![](README_files/figure-markdown_strict/unnamed-chunk-27-4.png)

    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"

![](README_files/figure-markdown_strict/unnamed-chunk-27-5.png)

    ## [1] "Plotting expression VS sample score"

![](README_files/figure-markdown_strict/unnamed-chunk-27-6.png)![](README_files/figure-markdown_strict/unnamed-chunk-27-7.png)![](README_files/figure-markdown_strict/unnamed-chunk-27-8.png)

    ## [1] "Plotting correlations of expression VS sample score"

    ## Warning in cor(x, y): l'écart type est nulle

![](README_files/figure-markdown_strict/unnamed-chunk-27-9.png)![](README_files/figure-markdown_strict/unnamed-chunk-27-10.png)

    ## Warning: Removed 1 rows containing missing values (geom_errorbar).

    ## Warning: Removed 1 rows containing missing values (geom_point).

![](README_files/figure-markdown_strict/unnamed-chunk-27-11.png)![](README_files/figure-markdown_strict/unnamed-chunk-27-12.png)![](README_files/figure-markdown_strict/unnamed-chunk-27-13.png)

    ## [1] "Plotting expression VS gene weight"

![](README_files/figure-markdown_strict/unnamed-chunk-27-14.png)![](README_files/figure-markdown_strict/unnamed-chunk-27-15.png)![](README_files/figure-markdown_strict/unnamed-chunk-27-16.png)

    ## [1] "Plotting correlation of expression VS gene weight"

![](README_files/figure-markdown_strict/unnamed-chunk-27-17.png)![](README_files/figure-markdown_strict/unnamed-chunk-27-18.png)![](README_files/figure-markdown_strict/unnamed-chunk-27-19.png)

    tictoc::toc()

    ## 25.688 sec elapsed

    PlotPCProjections(RomaData = RedData.NFC, Selected = NULL,
                      PlotPCProj = c('Points', 'Density',  'Bins'))

    ## [1] "1 geneset selected"

    ## Warning: Removed 10 rows containing missing values (geom_point).

![](README_files/figure-markdown_strict/unnamed-chunk-28-1.png)

    ## Warning: Removed 10 rows containing non-finite values (stat_density2d).

![](README_files/figure-markdown_strict/unnamed-chunk-28-2.png)

    ## Warning: Removed 10 rows containing non-finite values (stat_bin2d).

![](README_files/figure-markdown_strict/unnamed-chunk-28-3.png)

Using the interactive dashboard
-------------------------------

An interactive dasboard that can be used to perform ROMA and visualize
the results is available as a separate r package on
[GitHub](https://github.com/Albluca/rRomaDash)

Converting gene names
---------------------

Sometimes it may be necessary to convert gene names between different
conventions and possibly organisms. The rRoma package provides the
fucntion `ConvertNames` that can be used to conver the genes specified
by the Module list. The fucntion relies on the `biomaRt` package and
therefore a working internet connection is required. The `HOST` and
`PATH` argument are optional and can be used to specify an alternative
host when ensembl.org is having problems.

To derive and convert the MyoGS module list into the Ensembl format it
is sufficient to write:

    MyoGS <- SelectFromInternalDB("myogenesis")

    ## [1] "Searching in MsigDB v6.0"
    ## [1] "The following genesets have been selected:"
    ## [1] "REACTOME_MYOGENESIS (28 genes)"                 
    ## [2] "VANOEVELEN_MYOGENESIS_SIN3A_TARGETS (220 genes)"
    ## [3] "HALLMARK_MYOGENESIS (200 genes)"

    MyoGS.ENS <- ConvertModuleNames(ModuleList = MyoGS,
                                    SourceOrganism = "hsapiens", TargetOrganism = "hsapiens",
                       SourceTypes = "Names", TargetTypes = "Ensembl",
                       HOST = "grch37.ensembl.org", PATH = "/biomart/martservice")

    MyoGS.ENS[[1]]$Genes

    ##  [1] "ENSG00000064309" "ENSG00000044115" "ENSG00000067141"
    ##  [4] "ENSG00000188130" "ENSG00000129910" "ENSG00000196628"
    ##  [7] "ENSG00000185386" "ENSG00000111046" "ENSG00000140299"
    ## [10] "ENSG00000111049" "ENSG00000108984" "ENSG00000168036"
    ## [13] "ENSG00000081189" "ENSG00000144857" "ENSG00000280641"
    ## [16] "ENSG00000140262" "ENSG00000179242" "ENSG00000116604"
    ## [19] "ENSG00000122180" "ENSG00000112062" "ENSG00000071564"
    ## [22] "ENSG00000162068" "ENSG00000070831" "ENSG00000129152"
    ## [25] "ENSG00000097007" "ENSG00000170558" "ENSG00000068305"

It is also possible to convert gene names across organisms. To minimize
potential problematic results this is done only on genes with an
homology score:

    MyoGS.Mouse <- ConvertModuleNames(ModuleList = MyoGS, SourceOrganism = "hsapiens", TargetOrganism = "mmusculus",
                       SourceTypes = "Names", TargetTypes = "Names", HomologyLevel = 1,
                       HOST = "grch37.ensembl.org", PATH = "/biomart/martservice")

    ## [1] "Homology-supported gene conversion. Gene weigths will be deleted"

    MyoGS.Mouse[[1]]$Genes

    ##      CDON    CTNNA1      NEO1    MAPK12     CDH15      TCF4    MAPK11 
    ##    "Cdon"  "Ctnna1"    "Neo1"  "Mapk12"   "Cdh15"    "Tcf4"  "Mapk11" 
    ##      MYF6     BNIP2      MYF5    MAP2K6    CTNNB1     MEF2C       BOC 
    ##    "Myf6"   "Bnip2"    "Myf5"  "Map2k6"  "Ctnnb1"   "Mef2c"     "Boc" 
    ##     TCF12      CDH4     MEF2D      MYOG    MAPK14      TCF3      NTN3 
    ##   "Tcf12"    "Cdh4"   "Mef2d"    "Myog"  "Mapk14"    "Tcf3" "Tbc1d24" 
    ##     CDC42     MYOD1      ABL1      CDH2     MEF2A 
    ##   "Cdc42"   "Myod1"    "Abl1"    "Cdh2"   "Mef2a"

Converting across species, will suppress the weigth information.

Visualising on ACSN
-------------------

It is possible to project the results of rRoma analysis on a map. This
functionality is currently experimental. ACSN maps needs gene lists.
Therefore it is necessary to map the module activity to a gene specific
value. Since a gene can be associated with multiple modules, it is
possible to filter genes with a large variation in the gene weigths, as
this is associated with a larger variation, and hence uncertainlty, on
the contribution of the gene. The function takes several parameters:
`SampleName`, which indicates the name of the sample(s) to consider,
`AggScoreFun`, which describes the function used to group the module
scores when different samples are present, `FilterByWei`, which is a
parameters used to filter genes based on the variance of the gene weight
across modules, and `AggGeneFun`, which describes the function used to
obtain the value associated to a gene when the multiples weigths and
scores are present. The parameters `QTop`, `QBottom`, and `Steps`
control the heatmap of the plot. Finally, `DispMode` indicates which
quantities should be projected on the maps, the module score ("Module"),
or the gene weigths ("Gene").

In this example we will use the cell survival map and therefore we will
recompute ROMA using the associated GMT.

    tictoc::tic()
    Data.NFC.CS <- rRoma.R(ExpressionMatrix = MatData,
                           ModuleList = ReadGMTFile("http://acsn.curie.fr/files/survival_v1.1.gmt"),
                           FixedCenter = FALSE, MaxGenes = 500,
                           UseParallel = TRUE, nCores = 8, ClusType = "FORK",
                           PCSignMode="CorrelateAllWeightsByGene")

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_character()
    ## )

    ## See spec(...) for full column specifications.

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_character()
    ## )

    ## See spec(...) for full column specifications.

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_character()
    ## )

    ## See spec(...) for full column specifications.

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_character()
    ## )

    ## See spec(...) for full column specifications.

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_character()
    ## )

    ## See spec(...) for full column specifications.

    ## Parsed with column specification:
    ## cols(
    ##   X = col_character()
    ## )

    ## [1] "The following genesets have been loaded and selected:"
    ## [1] "HEDGEHOG (281 genes)"          "MAPK (208 genes)"             
    ## [3] "PI3K_AKT_MTOR (296 genes)"     "WNT_CANONICAL (431 genes)"    
    ## [5] "WNT_NON_CANONICAL (425 genes)"
    ## [1] "Centering gene expression over samples (genes will have 0 mean)"
    ## [1] "Using local center (NOT centering over genes)"
    ## [1] "All the genesets will be used"
    ## [1] "Too many cores selected! 7 will be used"
    ## [1] "2017-08-03 14:57:04 CEST"
    ## [1] "[1/5] Working on MAPK - na"
    ## [1] "207 genes available for analysis"
    ## [1] "Detecting outliers using leave one out and median-absolute-deviations away from median (scater package)"
    ## [1] "17 gene(s) will be filtered:"
    ##  [1] "GRK1"  "GRK7"  "TCF15" "TCF23" "TCF24" "BCL2"  "GPR1"  "GPR6" 
    ##  [9] "GPR12" "GPR15" "GPR17" "GPR20" "PRKCB" "PRKCG" "ERBB4" "FLT3" 
    ## [17] "MMP13"
    ## [1] "Not using weigths for PCA computation"
    ## [1] "Pre-filter data"
    ## [1] "L1 = 0.272356623388721 L1/L2 = 2.12933935438333"
    ## [1] "Median expression (uncentered): 13.4653750856794"
    ## [1] "Median expression (centered/weighted): 0.0364940949992434"
    ## [1] "Post-filter data"
    ## [1] "L1 = 0.308532038175362 L1/L2 = 2.44906266390702"
    ## [1] "Median expression (uncentered): 13.5243580596037"
    ## [1] "Median expression (centered/weighted): 0.0341023905754732"
    ## [1] "Previous sample size: 0"
    ## [1] "Next sample size: 207"
    ## [1] "Computing samples"
    ##    user  system elapsed 
    ##   0.117   0.171  48.487 
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "2017-08-03 14:57:56 CEST"
    ## [1] "[2/5] Working on HEDGEHOG - na"
    ## [1] "275 genes available for analysis"
    ## [1] "Detecting outliers using leave one out and median-absolute-deviations away from median (scater package)"
    ## [1] "28 gene(s) will be filtered:"
    ##  [1] "AKAP4"  "AKAP5"  "AKAP14" "GAS2"   "GPC5"   "GRK1"   "GRK7"  
    ##  [8] "GNG3"   "GNG4"   "GNG7"   "GNG8"   "GNG13"  "TCF15"  "TCF23" 
    ## [15] "TCF24"  "GNGT1"  "ARC"    "BCL2"   "FOXC2"  "FOXE1"  "FOXQ1" 
    ## [22] "NKX2-2" "PAX9"   "S100A7" "SFRP1"  "SIX1"   "SNAI1"  "HRK"   
    ## [1] "Not using weigths for PCA computation"
    ## [1] "Pre-filter data"
    ## [1] "L1 = 0.239363333828642 L1/L2 = 2.31207905299507"
    ## [1] "Median expression (uncentered): 13.382353827338"
    ## [1] "Median expression (centered/weighted): 0.0202799308952928"
    ## [1] "Post-filter data"
    ## [1] "L1 = 0.160857660690806 L1/L2 = 1.75047781401604"
    ## [1] "Median expression (uncentered): 13.4594315674689"
    ## [1] "Median expression (centered/weighted): 0.0225028382059111"
    ## [1] "Previous sample size: 207"
    ## [1] "Next sample size: 275"
    ## [1] "Computing samples"
    ##    user  system elapsed 
    ##   0.141   0.154  78.147 
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "2017-08-03 14:59:19 CEST"
    ## [1] "[3/5] Working on PI3K_AKT_MTOR - na"
    ## [1] "293 genes available for analysis"
    ## [1] "Detecting outliers using leave one out and median-absolute-deviations away from median (scater package)"
    ## [1] "7 gene(s) will be filtered:"
    ## [1] "PRKCG"   "ERBB4"   "FLT3"    "EIF4E1B" "PPP3R2"  "PPP2R2C" "FASLG"  
    ## [1] "Not using weigths for PCA computation"
    ## [1] "Pre-filter data"
    ## [1] "L1 = 0.232218378608347 L1/L2 = 1.95407750611921"
    ## [1] "Median expression (uncentered): 13.5306501745991"
    ## [1] "Median expression (centered/weighted): 0.0212445487028692"
    ## [1] "Post-filter data"
    ## [1] "L1 = 0.204674888557722 L1/L2 = 2.5716859618024"
    ## [1] "Median expression (uncentered): 13.5487015183355"
    ## [1] "Median expression (centered/weighted): 0.020895886038816"
    ## [1] "Previous sample size: 275"
    ## [1] "Next sample size: 293"
    ## [1] "Computing samples"
    ##    user  system elapsed 
    ##   0.167   0.162  88.649 
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "2017-08-03 15:00:53 CEST"
    ## [1] "[4/5] Working on WNT_NON_CANONICAL - na"
    ## [1] "412 genes available for analysis"
    ## [1] "Detecting outliers using leave one out and median-absolute-deviations away from median (scater package)"
    ## [1] "68 gene(s) will be filtered:"
    ##  [1] "CDH2"     "GNG3"     "GNG4"     "GNG7"     "GNG8"     "GNG13"   
    ##  [7] "GNGT1"    "MAPK4"    "MAPK10"   "PRKCB"    "PRKCG"    "ERBB4"   
    ## [13] "FLT3"     "PPP2R2C"  "PPP2R3A"  "ATP6V0D2" "FZD9"     "WNT1"    
    ## [19] "WNT2"     "WNT2B"    "WNT3"     "WNT3A"    "WNT5B"    "WNT6"    
    ## [25] "WNT7A"    "WNT7B"    "WNT8A"    "WNT8B"    "WNT9A"    "WNT9B"   
    ## [31] "WNT10B"   "WNT16"    "CTNND2"   "CTNNA2"   "CTNNA3"   "PRKAG3"  
    ## [37] "PRKAR2B"  "PRKAA2"   "PRKACG"   "ADCY2"    "ADCY5"    "ADCY8"   
    ## [43] "ADCY10"   "CAMK2A"   "CAMK2B"   "CAPN8"    "CAPN11"   "CAPN14"  
    ## [49] "CAPNS2"   "RSPO1"    "RSPO2"    "RSPO3"    "RSPO4"    "F2"      
    ## [55] "CASR"     "MEF2B"    "PPP3R2"   "TUBA3C"   "TUBA3D"   "TUBA3E"  
    ## [61] "ESR2"     "GATA1"    "GATA2"    "GATA4"    "GATA5"    "ACTA1"   
    ## [67] "SNAI1"    "IL2"     
    ## [1] "Not using weigths for PCA computation"
    ## [1] "Pre-filter data"
    ## [1] "L1 = 0.186860153257785 L1/L2 = 2.29958601053513"
    ## [1] "Median expression (uncentered): 13.3679607031876"
    ## [1] "Median expression (centered/weighted): 0.0365646440368517"
    ## [1] "Post-filter data"
    ## [1] "L1 = 0.154944841090472 L1/L2 = 1.77294425622764"
    ## [1] "Median expression (uncentered): 13.4977271172587"
    ## [1] "Median expression (centered/weighted): 0.0288109011789315"
    ## [1] "Previous sample size: 293"
    ## [1] "Next sample size: 412"
    ## [1] "Computing samples"
    ##    user  system elapsed 
    ##   0.190   0.174 171.266 
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "2017-08-03 15:03:54 CEST"
    ## [1] "[5/5] Working on WNT_CANONICAL - na"
    ## [1] "420 genes available for analysis"
    ## [1] "Detecting outliers using leave one out and median-absolute-deviations away from median (scater package)"
    ## [1] "59 gene(s) will be filtered:"
    ##  [1] "TCF15"    "TCF23"    "TCF24"    "GNG3"     "GNG4"     "GNG7"    
    ##  [7] "GNG8"     "GNG13"    "GNGT1"    "ATP6V0D2" "DAB1"     "FZD9"    
    ## [13] "KLHL1"    "KLHL4"    "KLHL10"   "KLHL30"   "KLHL33"   "KLHL34"  
    ## [19] "KLHL35"   "KLHL38"   "HECW1"    "HECW2"    "PPP1R1A"  "PPP1R1C" 
    ## [25] "PPP1R17"  "PPP1R27"  "PPP1R36"  "PPP1R42"  "WNT1"     "WNT2"    
    ## [31] "WNT2B"    "WNT3"     "WNT3A"    "WNT5B"    "WNT6"     "WNT7A"   
    ## [37] "WNT7B"    "WNT8A"    "WNT8B"    "WNT9A"    "WNT9B"    "WNT10B"  
    ## [43] "WNT16"    "CTNND2"   "DKK1"     "KREMEN2"  "NKD1"     "RSPO1"   
    ## [49] "CTNNA2"   "CTNNA3"   "SIX3"     "PYGO1"    "PRKAG3"   "PRKAR2B" 
    ## [55] "PRKAA2"   "PRKACG"   "ACTA1"    "MIR34A"   "WIF1"    
    ## [1] "Not using weigths for PCA computation"
    ## [1] "Pre-filter data"
    ## [1] "L1 = 0.180229465059627 L1/L2 = 2.16790854203311"
    ## [1] "Median expression (uncentered): 13.3356694801966"
    ## [1] "Median expression (centered/weighted): 0.00857729655253436"
    ## [1] "Post-filter data"
    ## [1] "L1 = 0.170229862732015 L1/L2 = 2.17315215905343"
    ## [1] "Median expression (uncentered): 13.477252323938"
    ## [1] "Median expression (centered/weighted): 0.00396167735379649"
    ## [1] "Previous sample size: 412"
    ## [1] "Next sample size: 420"
    ## [1] "Reusing previous sampling (Comparable metagene size)"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"
    ## [1] "Missing gene weights will be replaced by 1"
    ## [1] "Orienting PC by correlating gene expression and sample score (pearson)"
    ## [1] "Not using groups"
    ## [1] "Computing correlations"
    ## [1] "Correcting using weights"

    tictoc::toc()

    ## 422.034 sec elapsed

We can now project the information obtained by using the following
commands, which will open a windows in the default browser to visualize
the map. We will start by displaying information relative to normal
samples

    PlotOnACSN(RomaData = Data.NFC.CS, SampleName = names(Type[Type == "normal"]),
               AggScoreFun = "median", FilterByWei = 30, 
               DispMode = c("Module", "Gene"), DataName = "Normal", 
               QTop = .99, QBottom = .1,
               Selected = NULL,
               MapURL = 'https://acsn.curie.fr/navicell/maps/survival/master/index.php',
               AggGeneFun = "median")

    ## [1] "18 sample(s) selected"
    ## [1] "5 geneset(s) selected"

![](README_files/figure-markdown_strict/unnamed-chunk-34-1.png)

    ## [1] "The following genes will be ignored:"
    ## [1] "APC"   "CDH2"  "GSK3B" "MSH6"  "RUNX2" "UBE3D"

![](README_files/figure-markdown_strict/unnamed-chunk-34-2.png)![](README_files/figure-markdown_strict/unnamed-chunk-34-3.png)

    ## [1] "User-defined url"

    ## waiting for data to be imported...

    ## data imported.

![](README_files/figure-markdown_strict/unnamed-chunk-34-4.png)

    ## [1] "User-defined url"

    ## waiting for data to be imported...
    ## data imported.

![](README_files/figure-markdown_strict/unnamed-chunk-34-5.png)

This command will result in the following maps to be drawn:

![Roma Scores](README_files/Map1.png)

![Gene weights](README_files/Map2.png)

We can then compare the information with metastatic samples.

    PlotOnACSN(RomaData = Data.NFC.CS, SampleName = names(Type[Type == "metastasized"]),
               AggScoreFun = "median", FilterByWei = 30, 
               DispMode = "Module", DataName = "Metastasized", 
               QTop = .99, QBottom = .1,
               Selected = NULL,
               MapURL = 'https://acsn.curie.fr/navicell/maps/survival/master/index.php',
               AggGeneFun = "median")

    ## [1] "18 sample(s) selected"
    ## [1] "5 geneset(s) selected"

![](README_files/figure-markdown_strict/unnamed-chunk-35-1.png)

    ## [1] "The following genes will be ignored:"
    ## [1] "APC"   "CDH2"  "GSK3B" "MSH6"  "RUNX2" "UBE3D"

![](README_files/figure-markdown_strict/unnamed-chunk-35-2.png)![](README_files/figure-markdown_strict/unnamed-chunk-35-3.png)

    ## [1] "User-defined url"

    ## waiting for data to be imported...

    ## data imported.

![](README_files/figure-markdown_strict/unnamed-chunk-35-4.png)

![Roma Scores](README_files/Map3.png)

Session information
-------------------

    sessionInfo()

    ## R version 3.4.0 (2017-04-21)
    ## Platform: x86_64-redhat-linux-gnu (64-bit)
    ## Running under: CentOS Linux 7 (Core)
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /usr/lib64/R/lib/libRblas.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=fr_FR.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=fr_FR.UTF-8        LC_COLLATE=fr_FR.UTF-8    
    ##  [5] LC_MONETARY=fr_FR.UTF-8    LC_MESSAGES=fr_FR.UTF-8   
    ##  [7] LC_PAPER=fr_FR.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ## [1] readr_1.1.1          GEOquery_2.42.0      Biobase_2.36.2      
    ## [4] BiocGenerics_0.22.0  rRoma_0.0.4.2000     BiocInstaller_1.26.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] viridis_0.4.0        httr_1.2.1           edgeR_3.18.1        
    ##  [4] bit64_0.9-7          viridisLite_0.2.0    shiny_1.0.3         
    ##  [7] assertthat_0.2.0     stats4_3.4.0         blob_1.1.0          
    ## [10] vipor_0.4.5          yaml_2.1.14          RSQLite_2.0         
    ## [13] backports_1.1.0      lattice_0.20-35      glue_1.1.1          
    ## [16] limma_3.32.2         digest_0.6.12        ggsignif_0.3.0      
    ## [19] RColorBrewer_1.1-2   colorspace_1.3-2     htmltools_0.3.6     
    ## [22] httpuv_1.3.3         Matrix_1.2-9         plyr_1.8.4          
    ## [25] XML_3.98-1.9         pkgconfig_2.0.1      pheatmap_1.0.8      
    ## [28] biomaRt_2.32.1       zlibbioc_1.22.0      xtable_1.8-2        
    ## [31] scales_0.4.1         Rtsne_0.13           tibble_1.3.3        
    ## [34] IRanges_2.10.2       tictoc_1.0           ggplot2_2.2.1       
    ## [37] lazyeval_0.2.0       RJSONIO_1.3-0        magrittr_1.5        
    ## [40] mime_0.5             memoise_1.1.0        evaluate_0.10.1     
    ## [43] MASS_7.3-47          beeswarm_0.2.3       shinydashboard_0.6.1
    ## [46] tools_3.4.0          scater_1.4.0         data.table_1.10.4   
    ## [49] hms_0.3              matrixStats_0.52.2   stringr_1.2.0       
    ## [52] S4Vectors_0.14.3     munsell_0.4.3        locfit_1.5-9.1      
    ## [55] irlba_2.2.1          AnnotationDbi_1.38.1 bindrcpp_0.2        
    ## [58] colorRamps_2.3       compiler_3.4.0       rlang_0.1.1         
    ## [61] rhdf5_2.20.0         grid_3.4.0           RCurl_1.95-4.8      
    ## [64] tximport_1.4.0       rjson_0.2.15         labeling_0.3        
    ## [67] bitops_1.0-6         rmarkdown_1.6        gtable_0.2.0        
    ## [70] curl_2.7             reshape_0.8.6        DBI_0.7             
    ## [73] reshape2_1.4.2       R6_2.2.2             gridExtra_2.2.1     
    ## [76] knitr_1.16           dplyr_0.7.1          bit_1.1-12          
    ## [79] bindr_0.1            rprojroot_1.2        RNaviCell_0.2.1     
    ## [82] stringi_1.1.5        ggbeeswarm_0.5.3     Rcpp_0.12.11

Using the docker image
======================

A docker image contains rRoma, the [rRomaDash
package](https://github.com/Albluca/rRomaDash), and all the necessary
dependencies is available. It can be installed by typing

    docker pull albluca/rroma

After installation, the image can be started with the command

    docker run -p 8787:8787 albluca/rromadock

At this point the RStudio web interface containing rRoma will be
availabe, in Unix-like systems, at the address <http://localhost:8787>.
It is then possible to access the interface by using "rstudio" as the
username and password.

On windows, it may be necessary to replace localhost with the address of
the machine. If any problem is encountered see the instruction
[here](https://github.com/rocker-org/rocker/wiki/Using-the-RStudio-image).

Using the docker container
==========================

A docker image containing the RStusio web server and all the packages
needed to run rRoma is also available. After installing
[docker](http://www.docker.com), the image can be obtained using

    docker pull albluca/rroma

After installation, the image can be started with the command

    docker run -p 8787:8787 albluca/rroma

At this point, the RStudio web interface containing rRoma and the
[rRomaDash interface](https://github.com/Albluca/rRomaDash) will be
availabe.

In Unix-like systems, the interface will be available at the address
<http://localhost:8787>. It is then possible to access the interface by
using "rstudio" as the username and password.

On windows systems, it may be necessary to replace localhost with a
different address. On windows 10 professional, the docker image will be
accessuible at the address <http://127.0.0.1:8787>.

Further instructions to troubleshoot any difficulty related to accessing
the rstudio interface can be found
[here](https://github.com/rocker-org/rocker/wiki/Using-the-RStudio-image).

After starting the interface, it is possible to write r code as in the
desktop version of RStudio.
