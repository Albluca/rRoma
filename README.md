The rROMA package is being developed to simplify the use of [ROMA](https://github.com/sysbio-curie/Roma) from R. It is currently under development.

Intalling rROMA
---------------

The package is currently under development, non very well documented, and only available on GitHub. A functional Java Virtual machine for your system is necessary. Before installing the package, it is advisable to install rJava from sources by using

``` r
install.packages(pkgs = "rJava", repos="http://rforge.net", type = 'source')
```

Compiling from source requires the appropriate development tools, e.g., C/C++ compiler. The installation of the package requires the devtools package, which is available from CRAN. The rpgraph package can be installed using

``` r
install.packages("devtools")
library(devtools)
install_github("Albluca/rROMA")
```

Using rROMA
-----------

rROMA exports the single function `rRoma` which provides a subset of the functionalites of ROMA. The function requires a matrix with gene expression data across multiple samples and a module list. In the following example artificial gene expression matrix and module list will be created and ROMA is executed. Just as the ROMA executable, only `ExpressionMatrix` and `ModuleList` are mandatory. All of the parameter of the function have the same function of the Java executable. Parameters that are not present in the argument list have not been implemented yet.

``` r
library(rRoma)
```

    ## Loading required package: rJava

``` r
nGenes <- 300
nSamples <- 100
nModules <- 5
GenesPerModule <- 100

GeneNames <- paste("Gene_", 1:nGenes, sep='')
SampleNames <- paste("Sample_", 1:nSamples, sep='')

NumericalExpressionMatrix <- matrix(rlnorm(n = nGenes*nSamples), nrow = nGenes, nSamples)
ExpressionMatrix <- cbind(GeneNames, NumericalExpressionMatrix)
colnames(ExpressionMatrix) <- c("Gene", SampleNames)

ModuleList <- list()

for(i in 1:nModules){
  ModuleList[[i]] <- list(Name=paste("Module", i), Desc="Descripion", Genes=sample(GeneNames, GenesPerModule))
}


Tables <- rRoma(ExpressionMatrix = ExpressionMatrix, ModuleList = ModuleList, fillMissingValues = 0, centerData = TRUE,
                doubleCenterData = FALSE, outlierThreshold = 2, typeOfPCAUsage = 0, robustPCAcalculation = TRUE,
                robustPCAcalculationForSampling = TRUE, numberOfGeneSetSizesToSample = 5, mostContributingGenesZthreshold = 1,
                diffSpotGenesZthreshold = 1, correlationThreshold = 0.6, graphicalOutputThreshold = 0.05,
                minimalNumberOfGenesInModule = 10, maximalNumberOfGenesInModule = 1000, minimalNumberOfGenesInModuleFound = 8,
                numberOfPermutations = 0, typeOfModuleFile = 0, saveDecomposedFiles = TRUE)
```

    ## [1] "Constructing basic classes"
    ## [1] "Setting up working parameters"
    ## [1] "Loading gene expression information"
    ## [1] "Cleaning up memory"
    ## [1] "Filling up missing values"
    ## [1] "Centering data (zero mean for any gene)..."
    ## [1] "The operation is performed by R"
    ## [1] "Processing ModuleList"
    ## [1] "Decomposing datasets into modules"
    ## [1] "Compute module activity"
    ## [1] "Produce graphical output"
    ## [1] "Top contributing genes determination"
    ## [1] "Creating modile activity table"
    ## [1] "Cleaning up"

After running the function, `Tables` will be a list containing the tables that the Java command would have produced, with the only exception of tables encoded in both ".txt" and ".dat" files, in which case only ".txt" files are considered. The names of the elements of the list indicate the name of the file produced by the Java executable. For example, the module scores can be obtained by typing

``` r
Tables$module_scores.xls
```

    ##      [,1]       [,2]          [,3]        [,4]             
    ## [1,] "MODULE"   "L1"          "L1/L2"     "NUMBER_OF_GENES"
    ## [2,] "Module 1" "0.03751716"  "1.0511351" "94"             
    ## [3,] "Module 2" "0.038088135" "0.9692078" "94"             
    ## [4,] "Module 3" "0.036027104" "1.0404766" "95"             
    ## [5,] "Module 4" "0.038958933" "1.0217744" "93"             
    ## [6,] "Module 5" "0.03880043"  "1.0502253" "96"

Note that `ModuleList` must possess exactly three fields: `Name`, which is the name of the module, `Desc`, which is a description string, and `Genes`, which is an array of genes. Gene weights can be specified as in GMT files by appending the value encloded in a square parenthesis to the gene name, e.g., "Gene\_1\[+5\]".

`ExpressionMatrix` must be a numeric, with missing values encoded by `NA`. The first column must indicate the gene namas and the names of the columns must indicate the names of the samples (except for the first one, which marks the gene name column).
