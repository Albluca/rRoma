

#' Perform ROMA on a datasets
#'
#' @param ExpressionMatrix matrix, a numeric matrix containing the gene expression information. Columns indicate samples and rows indicated genes.
#' @param ModuleList list, gene module list
#' @param UseWeigths logical, should the weigths be used for PCA calculation?
#' @param ExpFilter logical, should the samples be filtered?
#' @param MinGenes integer, the minimum number of genes reported by a module available in the expression matrix to process the module
#' @param MaxGenes integer, the maximum number of genes reported by a module available in the expression matrix to process the module
#' @param nSamples integer, the number of randomized gene sampled (per module)
#' @param ApproxSamples integer between 0 and 100 the approximation parameter to reuse samples. This is the minimal percentage variation to reuse samples.
#' For example 5, means that samples re recalculated only if the number of genes in the geneset has increased by at least 5\%.
#' @param OutGeneNumber scalar, number of median-absolute-deviations away from median required for the total number of genes expressed in a sample to be called an outlier
#' @param Ncomp iteger, number of principal components used to filter samples in the gene expression space
#' @param OutGeneSpace scalar, number of median-absolute-deviations away from median required for in a sample to be called
#' an outlier in the gene expression space. If set to NULL, the gene space filtering will not be performed.
#' @param FixedCenter logical, should PCA with fixed center be used?
#' @param GeneOutDetection character scalar, the algorithm used to filter genes in a module. Possible values are
#' \itemize{
#' \item 'L1OutVarPerc': percentage variation relative to the median variance explained supporgted by a leave one out approach 
#' \item 'L1OutVarDC': dendrogram clustering statistics on variance explained supported by a leave one out approach
#' \item 'L1OutExpOut': number of median-absolute-deviations away from median explined variance
#' \item 'L1OutSdMean': Number of standard deviations away from the mean
#' }
#' The option "L1OutExpOut" requires the scater package to be installed.
#' @param GeneOutThr scalar, threshold used by gene filtering algorithm in the modules. It can represent maximum size of filtered cluster ("L1OutVarDC"), 
#' minimal percentage variation (L1OutVarPerc) or the number of median-absolute-deviations away from median ("L1OutExpOut")
#' @param GeneSelMode character scalar, mode used to sample genes: all available genes ("All") or genes not present in the module ("Others")
#' @param centerData logical, should the gene expression values be centered over the samples?
#' @param MoreInfo logical, shuold detailed information on the computation by printed?
#' @param PlotData logical, shuold debugging plots by produced ?
#' @param SampleFilter logical, should outlier detection be applied to sampled data as well?
#' @param PCADims integer, the number of PCA dimensions to compute. Should be >= 2. Note that, the value 1 is allowed,
#' but is not advisable under normal circumstances.
#' Larger values decrease the error in the estimation of the explained variance but increase the computation time.
#' @param DefaultWeight integer scalar, the default weigth to us if no weith is specified by the modile file and an algorithm requiring weigths is used
#' @param PCSignMode characrter scalar, the modality to use to determine the direction of the principal components. The following options are currentlhy available:
#' \itemize{
#' \item 'none' (The direction is chosen at random)
#' \item 'PreferActivation': the direction is chosen in such a way that the sum of the projection is positive
#' \item 'UseAllWeights': as 'PreferActivation', but the projections are multiplied by the weigths, missing weights are set to DefaultWeight
#' \item 'UseKnownWeights': as 'UseAllWeights', but missing weigth are set to 0
#' \item 'CorrelateAllWeightsByGene': the direction is chosen in such a way to maximise the positive correlation between the expression of genes with a positive (negative) weights
#' and the (reversed) PC projections, missing weights are set to DefaultWeight
#' \item 'CorrelateKnownWeightsByGene': as 'CorrelateAllWeights', but missing weights are set to 0
#' \item 'CorrelateAllWeightsBySample': the direction is chosen in such a way to maximise the positive correlation between the expression of genes and the PC corrected weigth
#' (i.e., PC weigths are multiplied by gene weigths), missing weights are set to DefaultWeight
#' \item 'CorrelateKnownWeightsBySample': as 'CorrelateAllWeightsBySample', but missing weights are set to 0
#' }
#' If 'CorrelateAllWeights', 'CorrelateKnownWeights', 'CorrelateAllWeightsBySample' or 'CorrelateKnownWeightsBySample' are used
#' and GroupPCSign is TRUE, the correltions will be computed on the groups defined by Grouping.
#' @param PCSignThr numeric scalar, a quantile threshold to limit the projections (or weights) to use, e.g., if equal to .9
#' only the 10\% of genes with the largest projection (or weights) in absolute value will be considered.
#' @param UseParallel boolean, shuold a parallel environment be used? Note that using a parallel environment will increase the memorey usage as a
#' copy of the gene expression matrix is needed for each core
#' @param nCores integer, the number of cores to use if UseParallel is TRUE. Set to NULL for auto-detection
#' @param ClusType string, the cluster type to use. The default value ("PSOCK") should be available on most systems, unix-like environments also support the "FORK",
#' which should be faster.
#' @param SamplingGeneWeights named vector, numeric. Weigth so use when correcting the sign of the PC for sampled data.
#' @param FillNAMethod names list, additional parameters to pass to the mice function
#' @param Grouping named vector, the groups associated with the sample.
#' @param FullSampleInfo boolean, should full PC information be computed and saved for all the randomised genesets?
#' @param GroupPCSign boolean, should grouping information to be used to orient PCs?
#' @param CorMethod character string indicating which correlation coefficient is to be used
#' for orienting the principal components. Can be "pearson", "kendall", or "spearman".
#' 
#' @return
#' @export
#'
#' @examples
PCrRoma.R <- function(ExpressionMatrix, centerData = TRUE, ExpFilter=FALSE, ModuleList, UseWeigths = FALSE,
                    DefaultWeight = 1, MinGenes = 10, MaxGenes = 1000, ApproxSamples = 5,
                    nSamples = 100, OutGeneNumber = 5, Ncomp = 10, OutGeneSpace = NULL, FixedCenter = TRUE,
                    GeneOutDetection = "L1OutExpOut", GeneOutThr = 5, GeneSelMode = "All", SampleFilter = TRUE,
                    MoreInfo = FALSE, PlotData = FALSE, PCADims = 2, PCSignMode ='none', PCSignThr = NULL,
                    UseParallel = FALSE, nCores = NULL, ClusType = "PSOCK", SamplingGeneWeights = NULL,
                    FillNAMethod = list(), Grouping = NULL, FullSampleInfo = FALSE, GroupPCSign = FALSE,
                    CorMethod = "pearson") {
  
  if(PCADims < 1){
    stop("PCADims should be >= 1")
  }
  
  if(PCADims == 1){
    print("PCADims = 1. Be aware that this is not advisable under normal circumstances")
    readline("Press any key")
  }
  
  if(is.null(SamplingGeneWeights)){
    SamplingGeneWeights = rep(DefaultWeight, nrow(ExpressionMatrix))
    names(SamplingGeneWeights) <- rownames(ExpressionMatrix)
  }
  
  if(any(is.na(ExpressionMatrix))){
    
    if(!require(mice)){
      stop("Unable to load mice. Impossible to proceed")
    } else {
      print("Filling NA with mice")
    }
    
    imp <- do.call(what = mice, args = append(list(data = ExpressionMatrix), FillNAMethod))
    ExpressionMatrix <- complete(imp)
    
  }
  
  ExpressionMatrix <- data.matrix(ExpressionMatrix)
  
  SAMPLE_WARNING <- 3
  
  AllGenesModule <- unique(unlist(lapply(ModuleList, "[[", "Genes")))
  AllGenesMatrix <- rownames(ExpressionMatrix)
  
  if(sum(AllGenesMatrix %in% AllGenesModule) < 3){
    stop("Not enough module genes found in the matrix. Impossible to proceed")
  }
  
  if(ncol(ExpressionMatrix) <= SAMPLE_WARNING & interactive()){
    print(paste("Only", ncol(ExpressionMatrix), "sample found"))
    print("The number of samples may be too small to guarantee a reliable analysis")
    Ans <- readline("Do you want to continue anyway? (y/n)")
    if(Ans != "y" & Ans != "Y"){
      return(NULL)
    }
  }
  
  if(any(AllGenesMatrix[duplicated(AllGenesMatrix)] %in% AllGenesModule)){
    stop("Module gene are not unique in the matrix. Impossible to proceed.")
  }
  
  if(any(duplicated(AllGenesMatrix)) & interactive()){
    print("Duplicated gene names detected. This may create inconsistencies in the analysis. Consider fixing this problem.")
    readline("Press any key")
  }
  
  if(FullSampleInfo & interactive()){
    print("PC projections and weigths will be computed and reoriented for sampled genesets. This is potentially very time consuming.")
    Ans <- readline("Are you sure you want to do that? (y/n)")
    if(Ans != "y" & Ans != "Y"){
      FullSampleInfo <- FALSE
      print("PC projections and weigths will NOT be computed and reoriented for sampled genesets.")
    } else {
      print("PC projections and weigths will be computed and reoriented for sampled genesets.")
    }
  }
  
  
  
  if(PlotData & interactive()){
    print("Diagnostic plots will be produced. This is time consuming and can produce errors expecially if done interactivelly.")
    print("It is advisable to only use this option if a relatively small number of genesets is analyzed and/or to redirect the graphic out (e.g. with pdf())")
    Ans <- readline("Are you sure you want to do that? (y/n)")
    if(Ans != "y" & Ans != "Y"){
      FullSampleInfo <- FALSE
      print("Diagnostic plots will NOT be produced.")
    } else {
      print("Diagnostic plots will be produced.")
    }
  }
  
  # Cleanup data
  
  if(ExpFilter){
    
    print("Filtering samples with abnormal numbers of genes detected")
    
    GeneDetectesOUT <- scater::isOutlier(colSums(ExpressionMatrix==0), nmads = OutGeneNumber)
    
    print(paste(sum(GeneDetectesOUT), "sample(s) will be filtered:"))
    if(sum(GeneDetectesOUT)>0){
      print(names(which(GeneDetectesOUT)))
    }
    
    print("Filtering samples in the gene space")
    
    if(is.null(Ncomp) | Ncomp > ncol(ExpressionMatrix)){
      Ncomp <- ncol(ExpressionMatrix)
    }
    
    if(!is.null(OutGeneSpace)){
      GeneSpaceOUT <- scater::isOutlier(rowSums(prcomp(t(ExpressionMatrix), center = TRUE, retx = TRUE)$x[,1:Ncomp]^2),
                                        nmads = OutGeneSpace)
      
      print(paste(sum(GeneSpaceOUT), "sample(s) will be filtered:"))
      if(sum(GeneSpaceOUT)>0){
        print(names(which(GeneSpaceOUT)))
      }
      
      ExpressionMatrix <- ExpressionMatrix[, !GeneDetectesOUT & !GeneSpaceOUT]
    }
    
  }
  
  # Look at groups
  
  if(!is.null(Grouping)){
    GrpNames <- unique(Grouping)
    GrpCol <- rainbow(length(GrpNames))
    names(GrpCol) <- GrpNames
    
    NotFoundSampNames <- which(!(colnames(ExpressionMatrix) %in% names(Grouping)))
    FoundSampNames <- which((colnames(ExpressionMatrix) %in% names(Grouping)))
    
    if(length(NotFoundSampNames)>1){
      print(paste("The following samples don't have an associated group:"))
      print(colnames(ExpressionMatrix)[NotFoundSampNames])
      print("They will NOT be considered for group associated analysis")
    }
    
  } else {
    FoundSampNames <- NULL
    Grouping <- rep("N/A", ncol(ExpressionMatrix))
    names(Grouping) <- colnames(ExpressionMatrix)
  }
  
  OrgExpMatrix = ExpressionMatrix
  
  if(centerData){
    print("Centering gene expression over samples")
    ExpressionMatrix <- t(scale(t(ExpressionMatrix), center = TRUE, scale = FALSE))
    GeneCenters <- attr(ExpressionMatrix, "scaled:center")
    attr(ExpressionMatrix, "scaled:center") <- NULL
  } else {
    warning("Skipping centering of gene expression over samples. This will generate inconsistencies if the data are not centered.")
    GeneCenters = rep(0, nrow(ExpressionMatrix))
  }
  
  if(FixedCenter){
    print("Using global center (centering over genes)")
    ExpressionMatrix <- scale(ExpressionMatrix, center = TRUE, scale = FALSE)
    SampleCenters <- attr(ExpressionMatrix, "scaled:center")
    attr(ExpressionMatrix, "scaled:center") <- NULL
    ModulePCACenter = FALSE
  } else {
    print("Using local center (NOT centering over genes)")
    ModulePCACenter = TRUE
    SampleCenters = rep(0, ncol(ExpressionMatrix))
  }
  
  ModuleSummary <- list()
  ModuleMatrix <- NULL
  
  ProjMatrix <- NULL
  WeigthList <- list()
  
  PVVectMat <- NULL
  
  OutLiersList <- list()
  UsedModules <- NULL
  
  # Filter genes for compatibility with the expression matrix
  
  for(i in 1:length(ModuleList)){
    Preserve <- ModuleList[[i]]$Genes %in% rownames(ExpressionMatrix)
    ModuleList[[i]]$Genes <- ModuleList[[i]]$Genes[Preserve]
    ModuleList[[i]]$Weigths <- ModuleList[[i]]$Weigths[Preserve]
  }
  
  # Filter genesets depending on the number of genes
  
  nGenes <- unlist(lapply(lapply(ModuleList, "[[", "Genes"), length))
  ToFilter <- (nGenes > MaxGenes | nGenes < MinGenes)
  ToUse <- !ToFilter
  
  if(sum(ToFilter)>1){
    print("The following geneset(s) will be ignored due to the number of genes being outside the specified range")
    print(unlist(lapply(ModuleList[ToFilter], "[[", "Name")))
    ModuleList <- ModuleList[ToUse]
  } else {
    print("All the genesets will be used")
  }
  
  
  
  # if(UseParallel){
  #   
  #   no_cores <- parallel::detectCores()
  #   
  #   if(is.null(nCores)){
  #     # Calculate the number of cores
  #     nCores = no_cores - 1
  #   }
  #   
  #   if(nCores >= no_cores){
  #     nCores = no_cores - 1
  #   }
  #   
  #   # Initiate cluster
  #   cl <- parallel::makeCluster(nCores, type = ClusType)
  #   
  #   parallel::clusterExport(cl=cl, varlist=c("SampleFilter", "GeneOutDetection", "GeneOutThr",
  #                                            "ModulePCACenter", "ExpressionMatrix", "DetectOutliers",
  #                                            "PCADims", "OrgExpMatrix", "FullSampleInfo", "FoundSampNames"), 
  #                           envir = environment())
  # }
  
  ModuleOrder <- order(unlist(lapply(lapply(ModuleList, "[[", "Genes"), length)))
  ModuleList <- ModuleList[ModuleOrder]
  
  OldSamplesLen <- 0
  
  for(i in 1:length(ModuleList)){
    
    print(Sys.time())
    
    print(paste("[", i, "/", length(ModuleList), "] Working on ", ModuleList[[i]]$Name, " - ", ModuleList[[i]]$Desc, sep = ""))
    
    CompatibleGenes <- ModuleList[[i]]$Genes
    
    print(paste(length(CompatibleGenes), "genes available for analysis"))
    
    if(MoreInfo){
      print("The following genes will be used:")
      print(CompatibleGenes)
    }
    
    # Filtering genes
    # SelGenes <- DetectOutliers(GeneOutDetection = GeneOutDetection, GeneOutThr = GeneOutThr, ModulePCACenter = ModulePCACenter,
    #                            CompatibleGenes = CompatibleGenes, ExpressionData = ExpressionMatrix[CompatibleGenes, ],
    #                            PlotData = PlotData, ModuleName = ModuleList[[i]]$Name)
    # 
    # if(length(SelGenes) > MaxGenes | length(SelGenes) < MinGenes){
    #   print("Number of selected genes outside the specified range")
    #   print("Skipping module")
    #   next()
    # } else {
    #   UsedModules <- c(UsedModules, i)
    # }
    # 
    # Keep track of outliers
    # OutLiersList[[i]] <- setdiff(CompatibleGenes, SelGenes)
    
    # Computing PC on the unfiltered data (only for reference)
    
    if(UseWeigths){
      print("Using weigths")
      Correction <- ModuleList[[i]]$Weigths
      names(Correction) <- CompatibleGenes
      Correction[!is.finite(Correction)] <- DefaultWeight
    } else {
      print("Not using weigths for PCA computation")
      Correction <- rep(1, length(CompatibleGenes))
      names(Correction) <- CompatibleGenes
    }
    
    BaseMatrix <- t(Correction[CompatibleGenes]*ExpressionMatrix[CompatibleGenes, ])
    
    BasePC <- analogue::prcurve(BaseMatrix, method = "ca", trace = TRUE, maxit = 50)
    
    ExpVar <- 1 - BasePC$dist/BasePC$totalDist
    
    PCBaseUnf <- PCBase
    ExpVarUnf <- ExpVar
    
    MedianExp <- median(OrgExpMatrix[CompatibleGenes, ])
    
    print("Pre-filter data")
    print(paste("L1 =", ExpVar[1], "L1/L2 =", ExpVar[1]/ExpVar[2]))
    print(paste("Median expression (uncentered):", MedianExp))
    print(paste("Median expression (centered/weighted):", median(BaseMatrix)))
    
    # Computing PC on the filtered data
    
    BaseMatrix <- t(Correction[SelGenes]*ExpressionMatrix[SelGenes, ])
    
    if(PCADims > min(dim(ExpressionMatrix[SelGenes, ])/3)){
      PCBase <- prcomp(x = BaseMatrix, center = ModulePCACenter, scale. = FALSE)
    } else {
      PCBase <- irlba::prcomp_irlba(x = BaseMatrix, n = PCADims, center = ModulePCACenter, scale. = FALSE, maxit = 10000)
    }
    
    ExpVar <- (PCBase$sdev^2)/sum(apply(scale(BaseMatrix, center = ModulePCACenter, scale = FALSE), 2, var))
    
    # ExpVar <- c(
    #   1/sum(apply(scale(Correction[CompatibleGenes]*ExpressionMatrix[CompatibleGenes, ], center = ModulePCACenter, scale = FALSE), 2, var)/PCBase$sdev[1]^2),
    #   1/sum(apply(scale(Correction[CompatibleGenes]*ExpressionMatrix[CompatibleGenes, ], center = ModulePCACenter, scale = FALSE), 2, var)/PCBase$sdev[2]^2)
    # )
    
    MedianExp <- median(OrgExpMatrix[SelGenes, ])
    
    print("Post-filter data")
    print(paste("L1 =", ExpVar[1], "L1/L2 =", ExpVar[1]/ExpVar[2]))
    print(paste("Median expression (uncentered):", MedianExp))
    print(paste("Median expression (centered/weighted):", median(BaseMatrix)))
    
    print(paste("Previous sample size:", OldSamplesLen))
    print(paste("Next sample size:", length(CompatibleGenes)))
    
    # Comparison with sample genesets
    if(nSamples > 0){
      
      if((length(CompatibleGenes)*(100-ApproxSamples)/100 <= OldSamplesLen) & (OldSamplesLen > 0)){
        
        if(length(CompatibleGenes) == OldSamplesLen){
          print("Reusing previous sampling (Same metagene size)")
        } else {
          print("Reusing previous sampling (Comparable metagene size)")
        }
        
      } else {
        
        print("Computing samples")
        
        # Define base analysis function
        
        TestGenes <- function(Gl){
          library(irlba)
          if(SampleFilter){
            SampleSelGenes <- DetectOutliers(GeneOutDetection = GeneOutDetection, GeneOutThr = GeneOutThr, ModulePCACenter = ModulePCACenter,
                                             CompatibleGenes = Gl, ExpressionData = ExpressionMatrix[Gl, ], PlotData = FALSE,
                                             ModuleName = '', PrintInfo = FALSE)
            if(length(SampleSelGenes)<PCADims){
              warning(paste("Size of filtered sample geneset too small (",  length(SampleSelGenes), "). This may cause inconsitencies. Increase MinGenes or GeneOutThr to prevent the problem"))
            }
          } else {
            SampleSelGenes <- Gl
          }
          
          if(length(SampleSelGenes) <= 1){
            SampMedian <- median(ExpressionMatrix[SampleSelGenes])
            warning(paste("Size of filtered sample geneset extremely small (",  length(SampleSelGenes), "). This may cause inconsitencies. Increase MinGenes or GeneOutThr to prevent the problem"))
            return(list("ExpVar"=c(1, rep(0, PCADims - 1)), "MedianExp"= SampMedian))
          }
          
          BaseMatrix <- t(ExpressionMatrix[SampleSelGenes, ])
          SampMedian <- median(OrgExpMatrix[SampleSelGenes, ])
          
          if(FullSampleInfo){
            
            if(length(SampleSelGenes) >= 3*PCADims){
              PCSamp <- irlba::prcomp_irlba(x = BaseMatrix, n = PCADims, center = ModulePCACenter, scale. = FALSE, retx = TRUE)
            } else {
              PCSamp <- prcomp(x = BaseMatrix, center = ModulePCACenter, scale. = FALSE, retx = TRUE)
            }
            
            VarVect <- PCSamp$sdev^2
            
            if(length(VarVect)>PCADims){
              VarVect <- VarVect[1:PCADims]
            }
            
            if(length(VarVect)<PCADims){
              VarVect <- c(VarVect, rep(0, PCADims - length(VarVect)))
            }
            
            ExpMat <- NULL
            if(PCSignMode %in% c("CorrelateAllWeightsByGene", "CorrelateKnownWeightsByGene",
                                 "CorrelateAllWeightsBySample", "CorrelateKnownWeightsBySample")){
              ExpMat <- OrgExpMatrix[SampleSelGenes, ]
            }
            
            CorrectSign1 <- FixPCSign(PCWeigth = PCSamp$rotation[,1], PCProj = PCSamp$x[,1],
                                      Wei = SamplingGeneWeights[SampleSelGenes],
                                      Mode = PCSignMode, DefWei = DefaultWeight, Thr = PCSignThr,
                                      Grouping = Grouping, ExpMat = ExpMat, CorMethod = CorMethod)
            if(PCADims > 1){
              CorrectSign2 <- FixPCSign(PCWeigth = PCSamp$rotation[,2], PCProj = PCSamp$x[,2],
                                        Wei = SamplingGeneWeights[SampleSelGenes],
                                        Mode = PCSignMode, DefWei = DefaultWeight, Thr = PCSignThr,
                                        Grouping = Grouping, ExpMat = ExpMat, CorMethod = CorMethod)
            } else {
              CorrectSign2 = NULL
            }
            
            
            return(list("ExpVar"=VarVect/sum(apply(scale(BaseMatrix, center = ModulePCACenter, scale = FALSE), 2, var)),
                        "MedianExp"= SampMedian,
                        "PCProj"=cbind(CorrectSign1*PCSamp$x[,1], CorrectSign2*PCSamp$x[,2]))
            )
            
          } else {
            
            if(length(SampleSelGenes) >= 3*PCADims){
              PCSamp <- irlba::prcomp_irlba(x = BaseMatrix, n = PCADims, center = ModulePCACenter, scale. = FALSE, retx = FALSE)
            } else {
              PCSamp <- prcomp(x = BaseMatrix, center = ModulePCACenter, scale. = FALSE, retx = FALSE)
            }
            
            VarVect <- PCSamp$sdev^2
            
            if(length(VarVect)>PCADims){
              VarVect <- VarVect[1:PCADims]
            }
            
            if(length(VarVect)<PCADims){
              VarVect <- c(VarVect, rep(0, PCADims - length(VarVect)))
            }
            
            return(list("ExpVar"=VarVect/sum(apply(scale(BaseMatrix, center = ModulePCACenter, scale = FALSE), 2, var)),
                        "MedianExp"= SampMedian, "PCProj"=NULL)
            )
            
          }
          
        }
        
        
        if(SampleFilter & !UseParallel){
          pb <- txtProgressBar(min = 0, max = nSamples, initial = 0, style = 3)
          GeneToSample <- length(SelGenes)
        } else {
          GeneToSample <- length(CompatibleGenes)
        }
        
        if(GeneSelMode == "All"){
          SampledsGeneList <- lapply(as.list(1:nSamples), function(i){sample(x = rownames(ExpressionMatrix), size = GeneToSample, replace = FALSE)})
        }
        
        if(GeneSelMode == "Others"){
          SampledsGeneList <- lapply(as.list(1:nSamples), function(i){sample(x = setdiff(rownames(ExpressionMatrix), SelGenes), size = GeneToSample, replace = FALSE)})
        }
        
        if(!exists("SampledsGeneList")){
          stop("Incorrect sampling mode")
        }
        
        
        
        if(!UseParallel){
          Time <- system.time(SampledExp <- lapply(SampledsGeneList, TestGenes), gcFirst = FALSE)
        } else {
          Time <- system.time(SampledExp <- parallel::parLapply(cl, SampledsGeneList, TestGenes), gcFirst = FALSE)
        }
        
        print(Time)
        
        SampleExpVar <- sapply(SampledExp, "[[", "ExpVar")
        SampleMedianExp <- sapply(SampledExp, "[[", "MedianExp")
        
        if(PCADims >= 2){
          SampleExpVar <- rbind(SampleExpVar[1,], SampleExpVar[1,]/SampleExpVar[2,], SampleExpVar[2,])
        } else {
          SampleExpVar <- rbind(SampleExpVar, rep(NA, length(SampleExpVar)), rep(NA, length(SampleExpVar)))
        }
        
        rownames(SampleExpVar) <- c("Sampled L1", "Sampled L1/L2", "Sampled L2")
        
        OldSamplesLen <- length(CompatibleGenes)
        
      }
      
    }
    
    
    if(PlotData){
      boxplot(SampleExpVar[1,], at = 1, ylab = "Explained variance",
              main = ModuleList[[i]]$Name, ylim = range(c(SampleExpVar[1,], SampleExpVar[1])))
      points(x=1, y=SampleExpVar[1], pch = 20, col="red", cex= 2)
      
      if(PCADims >= 2){
        boxplot(SampleExpVar[2,], at = 1, ylab = "Explained variance (PC1) / Explained variance (PC2)",
                log = "y", main = ModuleList[[i]]$Name, ylim=range(c(SampleExpVar[2,], SampleExpVar[1]/SampleExpVar[2])))
        points(x=1, y=SampleExpVar[1]/SampleExpVar[2], pch = 20, col="red", cex= 2)
      }
      
      boxplot(SampleMedianExp, at = 1, ylab = "Median expression",
              main = ModuleList[[i]]$Name, ylim=range(c(SampleMedianExp, median(BaseMatrix))))
      points(x=1, y=median(BaseMatrix), pch = 20, col="red", cex= 2)
      
    }
    
    L1Vect <- SampleExpVar[1,] - ExpVar[1]
    L1Vect <- L1Vect[is.finite(L1Vect)]
    
    L1L2Vect <- SampleExpVar[2,] - ExpVar[1]/ExpVar[2]
    L1L2Vect <- L1L2Vect[is.finite(L1L2Vect)]
    
    MedianVect <- SampleMedianExp - MedianExp
    MedianVect <- MedianVect[is.finite(MedianVect)]
    
    PVVect <- rep(NA, 6)
    
    if(length(L1Vect) > 5){
      PVVect[1] <- wilcox.test(L1Vect, alternative = "less")$p.value
      PVVect[2] <- wilcox.test(L1Vect, alternative = "greater")$p.value
    }
    
    if(length(L1L2Vect) > 5){
      PVVect[3] <- wilcox.test(L1L2Vect, alternative = "less")$p.value
      PVVect[4] <- wilcox.test(L1L2Vect, alternative = "greater")$p.value
    }
    
    if(length(MedianVect) > 5){
      PVVect[5] <- wilcox.test(MedianVect, alternative = "less")$p.value
      PVVect[6] <- wilcox.test(MedianVect, alternative = "greater")$p.value
    }
    
    PVVectMat <- rbind(PVVectMat, PVVect)
    
    ModuleMatrix <- rbind(ModuleMatrix,
                          c(ExpVar[1], sum(sign(SampleExpVar[1,] - ExpVar[1])==1)/nSamples,
                            ExpVar[1]/ExpVar[2], sum(sign(SampleExpVar[2,] - ExpVar[1]/ExpVar[2])==1)/nSamples,
                            MedianExp, sum(sign(MedianVect - MedianExp)==1)/nSamples))
    
    # Compute the sign correction
    
    if(GroupPCSign){
      GroupPCsVect <- Grouping
    } else {
      GroupPCsVect <- NULL
    }
    
    
    
    ExpMat <- NULL
    if(PCSignMode %in% c("CorrelateAllWeightsByGene", "CorrelateKnownWeightsByGene",
                         "CorrelateAllWeightsBySample", "CorrelateKnownWeightsBySample")){
      ExpMat <- OrgExpMatrix[CompatibleGenes, ]
    }
    
    CorrectSignUnf <- FixPCSign(PCWeigth = PCBaseUnf$rotation[,1], PCProj = PCBaseUnf$x[,1],
                                Wei = ModuleList[[i]]$Weigths[ModuleList[[i]]$Genes %in% CompatibleGenes],
                                Mode = PCSignMode, DefWei = DefaultWeight, Thr = PCSignThr,
                                Grouping = GroupPCsVect, ExpMat = ExpMat, CorMethod = CorMethod)
    
    
    
    ExpMat <- NULL
    if(PCSignMode %in% c("CorrelateAllWeightsByGene", "CorrelateKnownWeightsByGene",
                         "CorrelateAllWeightsBySample", "CorrelateKnownWeightsBySample")){
      ExpMat <- OrgExpMatrix[SelGenes, ]
    }
    
    CorrectSign1 <- FixPCSign(PCWeigth = PCBase$rotation[,1], PCProj = PCBase$x[,1],
                              Wei = ModuleList[[i]]$Weigths[ModuleList[[i]]$Genes %in% SelGenes],
                              Mode = PCSignMode, DefWei = DefaultWeight, Thr = PCSignThr,
                              Grouping = GroupPCsVect, ExpMat = ExpMat, CorMethod = CorMethod)
    
    
    if(PCADims >= 2){
      CorrectSign2 <- FixPCSign(PCWeigth = PCBase$rotation[,2], PCProj = PCBase$x[,2],
                                Wei = ModuleList[[i]]$Weigths[ModuleList[[i]]$Genes %in% SelGenes],
                                Mode = PCSignMode, DefWei = DefaultWeight, Thr = PCSignThr,
                                Grouping = GroupPCsVect, ExpMat = ExpMat, CorMethod = CorMethod)
      
    } else {
      CorrectSign2 <- NULL
    }
    
    if(PlotData){
      
      if(PCADims >= 2){
        DF <- data.frame(PC1 = CorrectSign1*PCBase$x[,1], PC2 = CorrectSign2*PCBase$x[,2], Group = Grouping[colnames(OrgExpMatrix[SelGenes, ])])
        
        p <- ggplot2::ggplot(DF, ggplot2::aes(x=PC1, y=PC2, color = Group)) + ggplot2::geom_point() +
          ggplot2::labs(title = ModuleList[[i]]$Name)
        print(p)
      }
      
      WeiVect <- ModuleList[[i]]$Weigths[ModuleList[[i]]$Genes %in% SelGenes]
      names(WeiVect) <- SelGenes
      
      if(PCSignMode %in% c('UseAllWeights', 'CorrelateAllWeightsBySample', 'CorrelateAllWeightsByGene')){
        WeiVect[is.na(WeiVect)] <- DefaultWeight
      }
      
      # print(WeiVect)
      
      if(any(!is.na(WeiVect))){
        
        LocMat <- OrgExpMatrix[SelGenes[!is.na(WeiVect)], ]
        
        MeltData <- reshape::melt(LocMat)
        colnames(MeltData) <- c("Gene", "Sample", "Exp")
        
        MeltData <- cbind(MeltData, Grouping[as.character(MeltData$Sample)])
        colnames(MeltData)[4] <- c("Group")
        
        MeltData <- cbind(MeltData, WeiVect[as.character(MeltData$Gene)])
        colnames(MeltData)[5] <- c("Wei")
        
        CorrProj <- CorrectSign1*PCBase$x[,1]
        names(CorrProj) <- colnames(OrgExpMatrix[SelGenes,])
        
        MeltData <- cbind(MeltData, CorrProj[as.character(MeltData$Sample)])
        colnames(MeltData)[6] <- c("Proj")
        
        CorrLoading <- CorrectSign1*PCBase$rotation[,1]*WeiVect
        names(CorrLoading) <- rownames(OrgExpMatrix[SelGenes,])
        
        MeltData <- cbind(MeltData, CorrLoading[as.character(MeltData$Gene)])
        colnames(MeltData)[7] <- c("Load")
        
        MeltData$Wei <- factor(MeltData$Wei)
        
        SplitGroups <- cut(seq(from=1, by=1, to = length(unique(MeltData$Gene))),
                           breaks = seq(from=0, to=length(unique(MeltData$Gene))+16, by=16))
        
        names(SplitGroups) <- unique(MeltData$Gene)
        
        print("Plotting expression VS projections")
        
        for(GeneGroup in levels(SplitGroups)){
          
          GenesToUse <- names(SplitGroups[as.character(SplitGroups) == GeneGroup])
          
          if(length(GenesToUse)>0){
            p <- ggplot2::ggplot(MeltData[as.character(MeltData$Gene) %in% GenesToUse,], ggplot2::aes(y=Exp, x=Proj, shape = Wei, color = Group)) + ggplot2::geom_point() +
              ggplot2::facet_wrap( ~ Gene) + ggplot2::labs(title = ModuleList[[i]]$Name, x = "PC1 projections", y = "Expression") +
              ggplot2::scale_shape_discrete(name = "Weight") + ggplot2::scale_color_discrete(name = "Group")
            print(p)
          }
          
          
        }
        
        print("Plotting correlations of expression VS projections")
        
        CorData <- apply(LocMat, 1, function(x){
          CT <- cor.test(x, CorrProj)
          return(c(CT$estimate, CT$conf.int))
        })
        
        CorData <- t(rbind(CorData, sign(WeiVect[colnames(CorData)])*sign(CorData[1,])))
        
        CorData <- cbind(rownames(CorData), CorData)
        colnames(CorData) <- c("Gene",  "Est", "Low", "High", "Conc")
        
        CorData <- data.frame(CorData)
        
        CorData$Est <- as.numeric(as.character(CorData$Est))
        CorData$Low <- as.numeric(as.character(CorData$Low))
        CorData$High <- as.numeric(as.character(CorData$High))
        
        # SplitGroups <- cut(seq(from=1, by=1, to = length(unique(CorData$Gene))),
        #                    breaks = seq(from=0, to=length(unique(CorData$Gene))+16, by=16))
        # 
        # names(SplitGroups) <- unique(MeltData$Gene)
        
        # print(CorData)
        
        for(GeneGroup in levels(SplitGroups)){
          
          GenesToUse <- names(SplitGroups[as.character(SplitGroups) == GeneGroup])
          
          # print(GenesToUse)
          
          if(length(GenesToUse)>0){
            p <- ggplot2::ggplot(CorData[as.character(CorData$Gene) %in% GenesToUse,], ggplot2::aes(x =  Gene, y = Est, ymin = Low, ymax = High, color = Conc)) +
              ggplot2::geom_hline(yintercept = 0, linetype = 2) + ggplot2::geom_errorbar() +
              ggplot2::geom_point() + ggplot2::coord_flip() + 
              ggplot2::labs(title = ModuleList[[i]]$Name, y = "Estimated correlation (Exp VS PC1 Proj - 95% CI)", x = "")
            
            print(p)
          }
          
        }
        
        
        
        print("Plotting expression VS PC weigths")
        
        for(GroupID in levels(MeltData$Group)){
          
          if(sum(as.character(MeltData$Group) == GroupID, na.rm = TRUE)>0){
            p <- ggplot2::ggplot(MeltData[as.character(MeltData$Group) == GroupID & !is.na(MeltData$Group),], ggplot2::aes(y=Exp, x=Load, shape = Wei, color = Group)) + ggplot2::geom_point() +
              ggplot2::facet_wrap( ~ Sample) + ggplot2::labs(title = ModuleList[[i]]$Name, x = "PC1 weights", y = "Expression") +
              ggplot2::scale_shape_discrete(name = "Gene weight") + ggplot2::scale_color_discrete(name = "Group")
            print(p)
          }
          
        }
        
        if(any(is.na(MeltData$Group))){
          tData <- MeltData[is.na(MeltData$Group),]
          tData$Group <- 'N/A'
          p <- ggplot2::ggplot(tData, ggplot2::aes(y=Exp, x=Load, shape = Wei, color = Group)) + ggplot2::geom_point() +
            ggplot2::facet_wrap( ~ Sample) + ggplot2::labs(title = ModuleList[[i]]$Name, x = "PC1 weights", y = "Expression") +
            ggplot2::scale_shape_discrete(name = "Gene weight") + ggplot2::scale_color_discrete(name = "Group")
          print(p)
        }
        
        
        
        
        print("Plotting correlation of expression VS PC weigths")
        
        CorData <- apply(LocMat, 2, function(x){
          CT <- cor.test(x[!is.na(CorrLoading)], CorrLoading[!is.na(CorrLoading)])
          return(c(CT$estimate, CT$conf.int))
        })
        
        CorData <- t(rbind(CorData, Grouping[colnames(CorData)]))
        
        CorData <- cbind(rownames(CorData), CorData)
        colnames(CorData) <- c("Gene",  "Est", "Low", "High", "Group")
        
        CorData <- data.frame(CorData)
        
        CorData$Est <- as.numeric(as.character(CorData$Est))
        CorData$Low <- as.numeric(as.character(CorData$Low))
        CorData$High <- as.numeric(as.character(CorData$High))
        
        
        # print(CorData)
        
        for(GroupID in levels(CorData$Group)){
          
          if(sum(as.character(CorData$Group) == GroupID, na.rm = TRUE)>0){
            p <- ggplot2::ggplot(CorData[as.character(CorData$Group) == GroupID & !is.na(CorData$Group),], ggplot2::aes(x =  Gene, y = Est, ymin = Low, ymax = High, color = Group)) +
              ggplot2::geom_hline(yintercept = 0, linetype = 2) + ggplot2::geom_errorbar() +
              ggplot2::geom_point() + ggplot2::coord_flip() + 
              ggplot2::labs(title = ModuleList[[i]]$Name, y = "Estimated correlation (Exp VS PC1 Wei - 95% CI)", x = "")
            
            print(p)
          }
          
        }
        
        if(any(is.na(CorData$Group))){
          tData <- CorData[is.na(CorData$Group),]
          tData$Group <- 'N/A'
          p <- ggplot2::ggplot(tData, ggplot2::aes(x =  Gene, y = Est, ymin = Low, ymax = High, color = Group)) +
            ggplot2::geom_hline(yintercept = 0, linetype = 2) + ggplot2::geom_errorbar() +
            ggplot2::geom_point() + ggplot2::coord_flip() + 
            ggplot2::labs(title = ModuleList[[i]]$Name, y = "Estimated correlation (Exp VS PC1 Wei - 95% CI)", x = "")
          
          print(p)
        }
        
      }
      
    }
    
    # Correct the sign of the first PC projections
    ModProjSamples <- CorrectSign1*PCBase$x[,1]
    names(ModProjSamples) <- colnames(ExpressionMatrix)
    
    ProjMatrix <- rbind(ProjMatrix, ModProjSamples)
    
    # Correct the sign of the first PC weigths (Filtered and unfiltered)
    tWeigths <- CorrectSign1*PCBase$rotation[,1]
    names(tWeigths) <- SelGenes
    
    tWeigthsUnf <- CorrectSignUnf*PCBaseUnf$rotation[,1]
    names(tWeigthsUnf) <- CompatibleGenes
    
    WeigthList[[length(WeigthList)+1]] <- tWeigths
    
    ModuleSummary[[length(ModuleSummary)+1]] <- list(ModuleName = ModuleList[[i]]$Name, ModuleDesc = ModuleList[[i]]$Desc,
                                                     OriginalGenes = CompatibleGenes, UsedGenes = SelGenes, SampledGenes = SampledsGeneList,
                                                     PCABase = PCBase, PCBaseUnf = PCBaseUnf,
                                                     CorrectSign1 = CorrectSign1, CorrectSign2 = CorrectSign2, ExpVarBase = ExpVar, ExpVarBaseUnf = ExpVarUnf, SampledExp = SampledExp,
                                                     PC1Weight.SignFixed = tWeigths, PC1WeightUnf.SignFixed = tWeigthsUnf,
                                                     GMTWei = ModuleList[[i]]$Weigths[ModuleList[[i]]$Genes %in% SelGenes])
    
  }
  
  if(UseParallel){
    
    # Stop cluster
    parallel::stopCluster(cl)
    
  }
  
  # Makes sure ModuleMatrix is treated as a matrix
  if(length(ModuleMatrix) == 6){
    dim(ModuleMatrix) <- c(1, 6)
  }
  
  colnames(ModuleMatrix) <- c("L1", "ppv L1", "L1/L2", "ppv L1/L2", "Median Exp", "ppv Median Exp")
  rownames(ModuleMatrix) <- unlist(lapply(ModuleList, "[[", "Name"))[UsedModules]
  
  # Makes sure PVVectMat is treated as a matrix
  if(length(PVVectMat) == 6){
    dim(PVVectMat) <- c(1, 6)
  }
  
  colnames(PVVectMat) <- c("L1 WT less pv", "L1 WT greater pv", "L1/L2 WT less pv", "L1/2 WT greater pv",
                           "Median Exp WT less pv", "Median Exp WT greater pv")
  rownames(PVVectMat) <- unlist(lapply(ModuleList, "[[", "Name"))[UsedModules]
  
  # Makes sure ProjMatrix is treated as a matrix
  if(length(ProjMatrix) == ncol(ExpressionMatrix)){
    dim(ProjMatrix) <- c(1, ncol(ExpressionMatrix))
  }
  
  colnames(ProjMatrix) <- colnames(ExpressionMatrix)
  rownames(ProjMatrix) <- unlist(lapply(ModuleList, "[[", "Name"))[UsedModules]
  
  if(length(UsedModules)>1){
    ReorderIdxs <- order(ModuleOrder[UsedModules])
    
    return(list(ModuleMatrix = ModuleMatrix[ReorderIdxs,], ProjMatrix = ProjMatrix[ReorderIdxs,], ModuleSummary = ModuleSummary[ReorderIdxs],
                WeigthList = WeigthList[ReorderIdxs], PVVectMat = PVVectMat[ReorderIdxs,], OutLiersList = OutLiersList[ReorderIdxs],
                GeneCenters = GeneCenters, SampleCenters = SampleCenters))
  } else {
    return(list(ModuleMatrix = ModuleMatrix, ProjMatrix = ProjMatrix, ModuleSummary = ModuleSummary,
                WeigthList = WeigthList, PVVectMat = PVVectMat, OutLiersList = OutLiersList,
                GeneCenters = GeneCenters, SampleCenters = SampleCenters))
  }
  
  
  
}




