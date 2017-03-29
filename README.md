This package provides an R implementation of [ROMA](http://journal.frontiersin.org/article/10.3389/fgene.2016.00018/full). The package is under active development and is currently being tested.

A Java implementation developed by Andrei Zynovyev and is also [available](https://github.com/sysbio-curie/Roma).

Intalling rROMA
---------------

The rRoma package relies on the `scater` package, which is available only on BioConductor. This package can be installed with the following command

``` r
source("https://bioconductor.org/biocLite.R")
biocLite("scater")
```

rRoma can then be installed using `devtools`

``` r
if(!require("devtools")){
  install.packages("devtools")
}
devtools::install_github("Albluca/rROMA")
```

The packages `GEOquery`, `tictoc`, and `readr` are not required to run rROMA, but are used in the following example and need to be installed to reproduce the analysis

``` r
if(!require("GEOquery")){
  source("https://bioconductor.org/biocLite.R")
  biocLite("GEOquery")
}

if(!require("readr")){
  install.packages("readr")
}

if(!require("tictoc")){
  install.packages("tictoc")
}
```

Using rROMA
-----------

The package can be loaded with the usual syntax, i.e., by typing

``` r
library(rRoma)
```

rRoma requires a gene expression matrix - with column names indicating samples and row names indicating gene names - and a module file containing information on the genesets that need to be evaluated. The module file can be loaded from a GMT file. Functions to automate the generation of then GMT file are also available.

Various functions are then available to explore the analysis, including plotting and statistical cross sample analysis.

An example
----------

To show a concrete example of rROMA we will use a dataset publicly available on GEO. The output of the analysis is not reported due to limitation in the markdown files used by GitHub. Let us begin by getting the description of the dataset

``` r
library(GEOquery)
gse <- getGEO("GSE50760", GSEMatrix = TRUE)
```

Then we get the actual expression expression files

``` r
filePaths = getGEOSuppFiles("GSE50760")
```

Now we can construct the expression matrix. Note that the code below is designed to work on a Unix-like environment (e.g. MacOS). The execution on a Windows environment may require replacing `"/"` with `"\"`.

``` r
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

SplitPath <- unlist(strsplit(rownames(filePaths)[1], "/"))
SplitPath <- SplitPath[-length(SplitPath)]
unlink(x = paste(SplitPath, collapse = "/"), recursive = TRUE)

Genes <- MatData[,1]
MatData <- data.matrix(data.frame(MatData[,-1]))
rownames(MatData) <- Genes
colnames(MatData) <- unlist(lapply(strsplit(Content, "_"), "[[", 1))
```

And look at the different groups of cells present

``` r
Type <- as.character(gse$GSE50760_series_matrix.txt.gz[[1]])
Type <- unlist(lapply(strsplit(Type, " "), "[[", 1))
names(Type) = as.character(gse$GSE50760_series_matrix.txt.gz[[2]])

table(Type)
```

For convenience, we will modify the sample names to accoun for the tissue of origin

``` r
colnames(MatData) <- paste(colnames(MatData),
                           paste("(", substr(Type[colnames(MatData)], 1, 1), ")",
                                 sep=""))
```

At this point we can create the metagene files. We will extract all the "HALLPARK" geneset from MSig.

``` r
AllHall <- SelectFromMSIGdb("HALLMARK")

AllHall <- lapply(AllHall, function(x){
  x$Name <- sub("HALLMARK_", "", x$Name)
  x
})
```

To reorient the PCs in such a way that gene with high expression are predominanty positive we will need to create two additional module list

``` r
tModule <- list(list(Name = "AllGenes", Desc = NA, Genes = unique(rownames(MatData)),
                     Weigths = rep(NA, length(unique(rownames(MatData))))))

ModuleListWU <- InferBinaryWeigth(ExpressionMatrix = MatData, ModuleList = AllHall, FillAllNA = TRUE)

tModuleWU <- InferBinaryWeigth(ExpressionMatrix = MatData, ModuleList = tModule, FillAllNA = TRUE)

SamplingGeneWeights <- tModuleWU[[1]]$Weigths
names(SamplingGeneWeights) <- tModuleWU[[1]]$Genes
```

To reduce potential problmes we will also rename genes with a duplicated name

``` r
if(any(duplicated(rownames(MatData)))){
  for(i in rownames(MatData)[duplicated(rownames(MatData))]){
    rownames(MatData)[rownames(MatData) %in% i] <-
      paste(rownames(MatData)[rownames(MatData) %in% i], 1:sum(rownames(MatData) %in% i))
  }
}
```

And now we are ready to perform ROMA without fixed center

``` r
tictoc::tic()
DataWU.NFC <- rRoma.R(ExpressionMatrix = MatData, centerData = TRUE, ExpFilter = FALSE,
                      ApproxSamples = 5, ModuleList = ModuleListWU, MinGenes = 5,
                      MaxGenes = 1000, nSamples = 100, UseWeigths = FALSE,
                      DefaultWeight = 1, FixedCenter = FALSE,
                      GeneOutDetection = 'L1OutVarPerc', GeneOutThr = 1,
                      GeneSelMode = "All", SampleFilter = TRUE, MoreInfo = FALSE,
                      PlotData = FALSE, PCSignMode = "UseAllWeigths", OutGeneNumber = 5,
                      Ncomp = 100, OutGeneSpace = 5, PCADims = 5, PCSignThr = NULL,
                      UseParallel = TRUE, nCores = 3, ClusType = "FORK",
                      SamplingGeneWeights = SamplingGeneWeights)
tictoc::toc()
```

and with fixed center

``` r
tictoc::tic()
DataWU.FC <- rRoma.R(ExpressionMatrix = MatData, centerData = TRUE, ExpFilter = FALSE,
                     ApproxSamples = 5, ModuleList = ModuleListWU, MinGenes = 5,
                     MaxGenes = 1000, nSamples = 100, UseWeigths = FALSE,
                     DefaultWeight = 1, FixedCenter = TRUE,
                     GeneOutDetection = 'L1OutVarPerc', GeneOutThr = 1,
                     GeneSelMode = "All", SampleFilter = TRUE, MoreInfo = FALSE,
                     PlotData = FALSE, PCSignMode = "UseAllWeigths", OutGeneNumber = 5,
                     Ncomp = 100, OutGeneSpace = 5, PCADims = 5, PCSignThr = NULL,
                     UseParallel = TRUE, nCores = 3, ClusType = "FORK",
                     SamplingGeneWeights = SamplingGeneWeights)
tictoc::toc()
```

Let us have a look at the overdispersed genesets for using the fixed center

``` r
Plot.Genesets(RomaData = DataWU.FC,
              Selected = SelectGeneSets(RomaData = DataWU.FC, VarThr = 1e-3,
                                        VarMode = "Wil", VarType = "Over"),
              GenesetMargin = 20, SampleMargin = 14, cluster_cols = TRUE)
```

and without the fixed center

``` r
Plot.Genesets(RomaData = DataWU.NFC,
              Selected = SelectGeneSets(RomaData = DataWU.NFC, VarThr = 1e-3,
                                        VarMode = "Wil", VarType = "Over"),
              GenesetMargin = 20, SampleMargin = 14, cluster_cols = TRUE)
```

Now we can look at the underdispersed geneset with fixed center

``` r
Plot.Genesets(RomaData = DataWU.FC,
              Selected = SelectGeneSets(RomaData = DataWU.FC, VarThr = 1e-3,
                                        VarMode = "Wil", VarType = "Under",
                                        MedThr = 1e-3, MedMode = "Wil", MedType = "Over"),
              GenesetMargin = 20, SampleMargin = 14, cluster_cols = TRUE)
```

and without the fixed center

``` r
Plot.Genesets(RomaData = DataWU.NFC,
              Selected = SelectGeneSets(RomaData = DataWU.NFC, VarThr = 1e-3,
                                        VarMode = "Wil", VarType = "Under",
                                        MedThr = 1e-3, MedMode = "Wil", MedType = "Over"),
              GenesetMargin = 20, SampleMargin = 14, cluster_cols = TRUE)
```

We can also explore the gene weigths with the fixed center

``` r
PlotGeneWeight(RomaData = DataWU.FC, PlotGenes = 30,
               ExpressionMatrix = MatData, LogExpression = FALSE,
               Selected = SelectGeneSets(RomaData = DataWU.FC, VarThr = 1e-3,
                                         VarMode = "Wil", VarType = "Over"),
               PlotWeigthSign = TRUE)
```

and without the fixed center

``` r
PlotGeneWeight(RomaData = DataWU.NFC, PlotGenes = 30,
               ExpressionMatrix = MatData, LogExpression = FALSE,
               Selected = SelectGeneSets(RomaData = DataWU.NFC, VarThr = 1e-3,
                                         VarMode = "Wil", VarType = "Over"),
               PlotWeigthSign = TRUE)
```

Moreover, we can look at the projections of the samples with the fixed center

``` r
PlotSampleProjections(RomaData = DataWU.FC, PlotSamples = 30,
                      ExpressionMatrix = MatData, LogExpression = FALSE,
                      Selected = SelectGeneSets(RomaData = DataWU.FC, VarThr = 1e-6,
                                                VarMode = "Wil", VarType = "Over"),
                      PlotPCProj = "Density")
```

and without the fixed center

``` r
PlotSampleProjections(RomaData = DataWU.NFC, PlotSamples = 30,
                      ExpressionMatrix = MatData, LogExpression = FALSE,
                      Selected = SelectGeneSets(RomaData = DataWU.NFC, VarThr = 1e-6,
                                                VarMode = "Wil", VarType = "Over"),
                      PlotPCProj = "Density")
```

Additionaly, we can compare across samples

``` r
CompareAcrossSamples(RomaData = DataWU.FC,
                     Selected = SelectGeneSets(RomaData = DataWU.FC, VarThr = 1e-3,
                                               VarMode = "Wil", VarType = "Over"),
                     Groups = Type)
```

Finally, we can look at genes which appear across different genesests and explore thier weithgs

``` r
PlotRecurringGenes(RomaData = DataWU.NFC,
                   Selected = SelectGeneSets(RomaData = DataWU.NFC, VarThr = 1e-3,
                                             VarMode = "Wil", VarType = "Over"),
                   GenesByGroup = 25, MinMult = 3)
```
