
#' Detect outliers
#'
#' @param GeneOutDetection 
#' @param GeneOutThr 
#' @param ModulePCACenter 
#' @param Genes 
#' @param ExpressionData 
#' @param PlotData 
#' @param ModuleName 
#'
#' @return
#' @export
#'
#' @examples
DetectOutliers <- function(GeneOutDetection, GeneOutThr, ModulePCACenter, CompatibleGenes, ExpressionData, PlotData = FALSE, ModuleName = '', PrintInfo = TRUE) {
  
  SelGenes <- CompatibleGenes
  
  if(GeneOutDetection == "PC1IQR"){
    if(PrintInfo){
      print("Detecting outliers using order statistics on PC1")
    }
    
    
    # Computing PC1 over samples
    PCA1Genes <- irlba::prcomp_irlba(x = ExpressionData, n = 1, center = ModulePCACenter, scale. = FALSE, retx = TRUE)$x
    
    
    if(PlotData){
      # Plotting rotated genes over PC1
      B <- boxplot(x = PCA1Genes, at = 1, horizontal = FALSE, ylab = "PC1 projection", main = ModuleName)
    }
    
    # Selecting outliers using mean and sd
    # Use Rank statistics
    GenesOut <- scater::isOutlier(PCA1Genes, nmads = GeneOutThr)
    
    if(PrintInfo){
      print(paste(sum(GenesOut), "genes will be filtered:"))
      if(sum(GenesOut)>0){
        print(CompatibleGenes[GenesOut])
      }
    }
    
    if(PlotData){
      # Highliting outliers
      points(y= PCA1Genes[GenesOut], x=rep(1, sum(GenesOut)), col='red', pch=20)
      legend("center", pch = 20, col='red', legend = "Outlier(s)")
    }
    
    # Updating gene list
    SelGenes <- CompatibleGenes[!GenesOut]
  }
  
  if(GeneOutDetection == "PC1DC"){
    if(PrintInfo){
      print("Detecting outliers using dendrogram clusterization on distanc from the center of PC1")
    }
    
    # Computing PC1 over samples
    PCA1Genes <- as.vector(irlba::prcomp_irlba(x = ExpressionData, n = 1, center = ModulePCACenter, scale. = FALSE, retx = TRUE)$x)
    PCA1GenesABS <- abs(PCA1Genes)
    names(PCA1GenesABS) <- CompatibleGenes
    
    HC <- hclust(dist(PCA1GenesABS))
    
    if(PlotData){
      # Plotting rotated genes over PC1
      plot(HC, xlab = "Genes", main = ModuleName)
    }
    
    DGroup <- cutree(HC, k = 2)
    G1 <- which(DGroup == 1)
    G2 <- which(DGroup == 2)
    
    if(length(G1) <= GeneOutThr | length(G2) <= GeneOutThr){
      
      if(length(G1) < length(G2)){
        GenesOut <- DGroup == 1
      } else {
        GenesOut <- DGroup == 2
      }
      
    }
    
    if(PrintInfo){
      print(paste(sum(GenesOut), "genes will be filtered:"))
      if(sum(GenesOut)>0){
        print(CompatibleGenes[GenesOut])
      }
    }
    
    # Updating gene list
    SelGenes <- CompatibleGenes[!GenesOut]
  }
  
  if(GeneOutDetection == "L1OutVarPerc"){
    
    if(PrintInfo){
      print("Detecting outliers using leave one out and percentage variation on variance explained by PC1")
    }
    
    # Computing all the PC1
    AllPCA1 <- sapply(as.list(1:length(CompatibleGenes)), function(i){
      PC1Var <- irlba::prcomp_irlba(x = ExpressionData[-i, ], n = 1, center = ModulePCACenter, scale. = FALSE)$sdev^2
      return(PC1Var/sum(apply(scale(ExpressionData[-i, ], center = ModulePCACenter, scale = FALSE), 2, var)))
    })
    
    
    # Getting the distance from the median PC1
    GenesOut <- abs(AllPCA1-median(AllPCA1)) > GeneOutThr/100
    
    if(PlotData){
      # Plotting the distances from the median PC1
      B <- boxplot(x = AllPCA1, at = 1, horizontal = FALSE, ylab = "Variance explained by PC1", main = ModuleName)
    }
    
    if(PrintInfo){
      print(paste(sum(GenesOut), "genes will be filtered:"))
      if(sum(GenesOut)>0){
        print(CompatibleGenes[GenesOut])
      }
    }

    
    if(PlotData){
      # Highliting outliers
      points(y= AllPCA1[GenesOut], x=rep(1, sum(GenesOut)), col='red', pch=20)
      legend("left", pch = 20, col='red', legend = "Outlier(s)")
    }
    
    # Updating gene list
    SelGenes <- CompatibleGenes[!GenesOut]
    
  }
  
  if(GeneOutDetection == "L1OutVarDC"){
    print("Detecting outliers using leave one out and percentage variation on variance explained by PC1")
    
    # Computing all the PC1
    AllPCA1 <- sapply(as.list(1:length(CompatibleGenes)), function(i){
      PC1Var <- irlba::prcomp_irlba(x = ExpressionData[-i, ], n = 1, center = ModulePCACenter, scale. = FALSE)$sdev^2
      return(PC1Var/sum(apply(scale(ExpressionData[-i, ], center = ModulePCACenter, scale = FALSE), 2, var)))
    })
    
    
    AllPCA1ABS <- abs(AllPCA1 - median(AllPCA1))
    names(AllPCA1ABS) <- CompatibleGenes
    
    HC <- hclust(dist(AllPCA1ABS))
    
    if(PlotData){
      # Plotting rotated genes over PC1
      plot(HC, xlab = "Genes", main = ModuleName)
    }
    
    DGroup <- cutree(HC, k = 2)
    G1 <- which(DGroup == 1)
    G2 <- which(DGroup == 2)
    
    if(length(G1) <= GeneOutThr | length(G2) <= GeneOutThr){
      
      if(length(G1) < length(G2)){
        GenesOut <- DGroup == 1
      } else {
        GenesOut <- DGroup == 2
      }
      
    }
    
    
    print(paste(sum(GenesOut), "genes will be filtered:"))
    if(sum(GenesOut)>0){
      print(CompatibleGenes[GenesOut])
    }
    
    # Updating gene list
    SelGenes <- CompatibleGenes[!GenesOut]
    
  }
  
  return(SelGenes)
  
}














#' Correct the sign of the principal component 
#'
#' @param PC1Rotation vector, numeric Principal component values
#' @param Wei vector, numeric optional vector of weigths
#' @param Mode scalar, character. Mode to correct the sign
#' @param Thr scalar, numeric quantile threshold 
#'
#' @return
#' @export
#'
#' @examples
FixPCSign <- function(PC1Projections, Wei = NULL, Mode ='none', DefWei = 1, Thr = NULL) {

  if(Mode == 'none'){
    return(1)
  }
  
  if(Mode == 'PreferActivation'){
    
    ToUse <- rep(TRUE, length(PC1Projections))
    if(!is.null(Thr)){
      ToUse <- abs(PC1Projections) > quantile(abs(PC1Projections), Thr) 
    }
      
    if(sum(PC1Projections[ToUse])<0){
      return(-1)
    } else {
      return(+1)
    }
  }
  
  if(Mode == 'UseAllWeigths'){
    Wei[is.na(Wei)] <- DeafaultWei
    
    ToUse <- rep(TRUE, length(PC1Projections))
    if(!is.null(Thr)){
      ToUse <- abs(PC1Projections) > quantile(abs(PC1Projections), Thr) 
    }
    
    if(sum(Wei[ToUse]*PC1Projections[ToUse])<0){
      return(-1)
    } else {
      return(+1)
    }
  }
  
  if(Mode == 'UseKnownWeigths'){
    
    ToUse <- rep(TRUE, length(PC1Projections))
    if(!is.null(Thr)){
      ToUse <- abs(PC1Projections) > quantile(abs(PC1Projections), Thr) 
    }
    
    if(sum(Wei[is.na(Wei) & ToUse]*PC1Projections[is.na(Wei) & ToUse])<0){
      return(-1)
    } else {
      return(+1)
    }
  }
  
}










#' Perform ROMA on a datasets
#'
#' @param ExpressionMatrix matrix, a numeric matrix containing the gene expression information. Columns indicate samples and rows indicated genes.
#' @param ModuleList list, gene module list
#' @param UseWeigths logical, should the weigths be used
#' @param DefaultWeigt integer, the value to use when gene weigth is not present
#' @param ExpFilter logical, should the samples be filtered?
#' @param MinGenes integer, the minimum number of genes reported by a module available in the expression matrix to process the module
#' @param MaxGenes integer, the maximum number of genes reported by a module available in the expression matrix to process the module
#' @param nSamples integer, the number of randomized gene sampled (per module)
#' @param ApproxSamples integer, the approximation parameter to reuse samples
#' @param OutGeneNumber scalar, number of median-absolute-deviations away from median required for the total number of genes expressed in a sample to be called an outlier
#' @param Ncomp iteger, number of principal components used to filter samples in the gene expression space
#' @param OutGeneSpace scalar, number of median-absolute-deviations away from median required for in a sample to be called an outlier in the gene expression space
#' @param FixedCenter logical, should PCA with fixed center be used?
#' @param GeneOutDetection character scalar, the algorithm used to filter genes in a module. Possible values are "PC1IQR" (Projection on PC1, order statistics),
#' "PC1DC" (Projection on PC1, dendrogram clustering), "L1OutVarPerc" (Percentage variation relative to the median variance explained supporgted by a leave one out approach), and
#' "L1OutVarDC" (Dendrogram clustering statistics on variance explained supported by a leave one out approach)
#' @param GeneOutThr scalar, threshold used by gene filtering algorithm in the modules. It can represent maximum size of filtered cluster ("PC1DC" and "L1OutVarDC"),
#' minimal percentage variation (L1OutVarPerc), or minimal median-absolute-deviations.
#' @param GeneSelMode character scalar, mode used to sample genes: all available genes ("All") or genes not present in the module ("Others")
#' @param centerData logical, should the gene expression values be centered over the samples?
#' @param MoreInfo logical, shuold detailed information on the computation by printed?
#' @param PlotData logical, shuold debugging plots by produced ?
#' @param SampleFilter logical, should outlier detection be applied to sampled to sample data as well?
#' @param PCADims integer, the number of PCA dimensions to compute. Should be >= 2.
#' Larger values decrease the error in the estimation of the explained variance but increase the computation time.
#' @param Mode scalar, string. The mode to be used to correct the sign of the 1st PC (Not implemented yet)
#' @param Thr scalar, numeric, the threshould to be used whn correcting the 1st PC (Not implemented yet) 
#'
#' @return
#' @export
#'
#' @examples
rRoma.R <- function(ExpressionMatrix, centerData = TRUE, ExpFilter=FALSE, ModuleList, UseWeigths = FALSE,
                    DefaultWeight = 1, MinGenes = 10, MaxGenes = 1000, ApproxSamples = 2,
                    nSamples = 100, OutGeneNumber = 5, Ncomp = 10, OutGeneSpace = 5, FixedCenter = TRUE,
                    GeneOutDetection = "PC1IQR", GeneOutThr = 5, GeneSelMode = "All", SampleFilter = FALSE,
                    MoreInfo = FALSE, PlotData = FALSE, PCADims = 2, PC1SignMode ='none', PC1SignThr = NULL) {

  AllGenes <- NULL
  for(i in 1:length(ModuleList)){
    AllGenes <- c(AllGenes, ModuleList[[i]]$Genes)
  }
  
  
  AllGenesModule <- unique(AllGenes)
  AllGenesMatrix <- rownames(ExpressionMatrix)
  
  if(any(AllGenesMatrix[duplicated(AllGenesMatrix)] %in% AllGenesModule)){
    stop("Module gene are not unique in the matrix. Impossible to proceed.")
  }
  
  if(any(duplicated(AllGenesMatrix))){
    warning("Duplicated gene names detected. This may create inconsistencies in the analysis. Consider fixing this problem.")
    readline("Press any key")
  }
  
  # Cleanup data
  
  if(ExpFilter){
    
    print("Filtering samples with abnormal numbers of genes detected")
    
    GeneDetectesOUT <- scater::isOutlier(colSums(ExpressionMatrix==0), nmads = OutGeneNumber)
    
    print(paste(sum(GeneDetectesOUT), "samples will be filtered:"))
    if(sum(GeneDetectesOUT)>0){
      print(names(which(GeneDetectesOUT)))
    }
    
    print("Filtering samples in the gene space")
    
    if(is.null(Ncomp) | Ncomp > ncol(ExpressionMatrix)){
      Ncomp <- ncol(ExpressionMatrix)
    }
    
    GeneSpaceOUT <- scater::isOutlier(rowSums(prcomp(t(ExpressionMatrix), center = TRUE, retx = TRUE)$x[,1:Ncomp]^2),
                                      nmads = OutGeneSpace)
    
    print(paste(sum(GeneSpaceOUT), "samples will be filtered:"))
    if(sum(GeneSpaceOUT)>0){
      print(names(which(GeneSpaceOUT)))
    }
    
    ExpressionMatrix <- ExpressionMatrix[, !GeneDetectesOUT & !GeneSpaceOUT]
    
  }
  
  if(centerData){
    print("Centering gene expression over samples")
    ExpressionMatrix <- t(scale(t(ExpressionMatrix), center = centerData, scale = FALSE))
  }
  
  if(FixedCenter){
    print("Using global center")
    ExpressionMatrix <- scale(ExpressionMatrix, center = TRUE, scale = FALSE)
    ModulePCACenter = FALSE
  } else {
    ModulePCACenter = TRUE
  }
  
  ModuleSummary <- list()
  ModuleMatrix <- NULL
  
  ProjLists <- list()
  PC1Matrix <- NULL
  
  PVVectMat <- NULL
  
  OutLiersList <- list()
  
  
  # Filter genes for compatibility with the expression matrix
  
  for(i in 1:length(ModuleList)){
    Preserve <- ModuleList[[i]]$Genes %in% rownames(ExpressionMatrix)
    ModuleList[[i]]$Genes <- ModuleList[[i]]$Genes[Preserve]
    ModuleList[[i]]$Weigths <- ModuleList[[i]]$Weigths[Preserve]
  }
  
  ModuleList <- ModuleList[order(unlist(lapply(lapply(ModuleList, "[[", "Genes"), length)))]
  
  OldSamplesLen <- 0
  
  for(i in 1:length(ModuleList)){
    
    print(paste("Working on", ModuleList[[i]]$Name, "-", ModuleList[[i]]$Desc))
    
    CompatibleGenes <- ModuleList[[i]]$Genes
    
    print(paste(length(CompatibleGenes), "genes available for analysis"))
    
    if(length(CompatibleGenes) > MaxGenes | length(CompatibleGenes) < MinGenes){
      print("Number of genes outside the specified range")
      print("Skipping module")
      next()
    }
    
    if(MoreInfo){
      print("The following genes will be used:")
      print(CompatibleGenes)
    }
    
    # Filtering genes
    SelGenes <- DetectOutliers(GeneOutDetection = GeneOutDetection, GeneOutThr = GeneOutThr, ModulePCACenter = ModulePCACenter,
                               CompatibleGenes = CompatibleGenes, ExpressionData = ExpressionMatrix[CompatibleGenes, ], PlotData = PlotData,
                               ModuleName = ModuleList[[i]]$Name)
    
    # Keep track of outliers
    OutLiersList[[i]] <- setdiff(CompatibleGenes, SelGenes)
    
    
    # Computing PC on the unfiltered data (only for reference)
    
    if(UseWeigths){
      print("Using weigths")
      Correction <- ModuleList[[i]]$Weigths
      names(Correction) <- CompatibleGenes
      Correction[!is.finite(Correction)] <- DefaultWeight
    } else {
      print("Not using weigths")
      Correction <- rep(1, length(CompatibleGenes))
      names(Correction) <- CompatibleGenes
    }
    
    PCBase <- irlba::prcomp_irlba(x = Correction[CompatibleGenes]*ExpressionMatrix[CompatibleGenes, ], n = PCADims, center = ModulePCACenter, scale. = FALSE)
    ExpVar <- (PCBase$sdev^2)/sum(apply(scale(Correction[CompatibleGenes]*ExpressionMatrix[CompatibleGenes, ], center = ModulePCACenter, scale = FALSE), 2, var))
    
    print("Pre-filter data")
    print(paste("L1 =", ExpVar[1], "L1/L2 =", ExpVar[1]/ExpVar[2]))
    
    # Computing PC on the filtered data
    PCBase <- irlba::prcomp_irlba(x = Correction[SelGenes]*ExpressionMatrix[SelGenes, ], n = PCADims, center = ModulePCACenter, scale. = FALSE)
    ExpVar <- (PCBase$sdev^2)/sum(apply(scale(Correction[SelGenes]*ExpressionMatrix[SelGenes, ], center = ModulePCACenter, scale = FALSE), 2, var))
    
    print("Post-filter data")
    print(paste("L1 =", ExpVar[1], "L1/L2 =", ExpVar[1]/ExpVar[2]))
    
    print(paste("Previous sample size:", OldSamplesLen))
    print(paste("Nex sample size:", length(CompatibleGenes)))
    
    # Comparison with sample genesets
    if(SampleFilter){
      
      if((length(CompatibleGenes) <= OldSamplesLen + ApproxSamples) & (OldSamplesLen > 0)){
        
        if(length(CompatibleGenes) == OldSamplesLen){
          print("Reusing previous sampling (Same metagene size)")
        } else {
          print("Reusing previous sampling (Comparable metagene size)")
        }
        
      } else {
        
        print("Computing samples")
        
        if(SampleFilter){
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
        
        SampledExp <- sapply(as.list(1:length(SampledsGeneList)), function(i){
          if(SampleFilter){
            SampleSelGenes <- DetectOutliers(GeneOutDetection = GeneOutDetection, GeneOutThr = GeneOutThr, ModulePCACenter = ModulePCACenter,
                                             CompatibleGenes = SampledsGeneList[[i]], ExpressionData = ExpressionMatrix[SampledsGeneList[[i]], ], PlotData = FALSE,
                                             ModuleName = '', PrintInfo = FALSE)
            setTxtProgressBar(pb, i)
          } else {
            SampleSelGenes <- SampledsGeneList[[i]]
          }
          PCSamp <- irlba::prcomp_irlba(x = ExpressionMatrix[SampleSelGenes, ], n = PCADims, center = ModulePCACenter, scale. = FALSE)
          return((PCSamp$sdev^2)/sum(apply(scale(ExpressionMatrix[SampleSelGenes, ], center = ModulePCACenter, scale = FALSE), 2, var)))
        })
        
        SampledExp <- rbind(SampledExp[1,], SampledExp[1,]/SampledExp[2,])
        
      }
      
      OldSamplesLen <- length(CompatibleGenes)
      
    }
    
    if(PlotData){
      boxplot(SampledExp[1,], at = 1, ylab = "Explained Variance",
              main = ModuleList[[i]]$Name, ylim = c(0,1))
      points(x=1, y=ExpVar[1], pch = 20, col="red", cex= 2)
      
      
      boxplot(SampledExp[2,], at = 1, ylab = "Explained variance (PC1) / Explained variance (PC2)", log = "y")
      points(x=1, y=ExpVar[1]/ExpVar[2], pch = 20, col="red", cex= 2)
    }
    
    PVVect <- c(
      wilcox.test(SampledExp[1,] - ExpVar[1], alternative = "less")$p.value,
      wilcox.test(SampledExp[1,] - ExpVar[1], alternative = "greater")$p.value,
      wilcox.test(SampledExp[2,] - ExpVar[1]/ExpVar[2], alternative = "less")$p.value,
      wilcox.test(SampledExp[2,] - ExpVar[1]/ExpVar[2], alternative = "greater")$p.value
    )
    
    PVVectMat <- rbind(PVVectMat, PVVect)
    
    ModuleMatrix <- rbind(ModuleMatrix,
                          c(ExpVar[1], sum(sign(SampledExp[1,] - ExpVar[1])==1)/nSamples,
                            ExpVar[1]/ExpVar[2], sum(sign(SampledExp[2,] - ExpVar[1]/ExpVar[2])==1)/nSamples))
    
    ###### WORK HERE
    
    CorrectSign <- FixPCSign(PCBase$x[,1], Wei = ModuleList[[i]]$Weigths[ModuleList[[i]]$Genes %in% SelGenes],
                             Mode = PC1SignMode, DefWei = DefaultWeight, Thr = PC1SignThr)
    
    ModProjGenes <- CorrectSign*PCBase$x[,1]
    names(ModProjGenes) <- SelGenes
    
    ProjLists[[i]] <- ModProjGenes
    
    PC1Matrix <- rbind(PC1Matrix, CorrectSign*PCBase$rotation[,1])
    
    ModuleSummary[[i]] <- list(ModuleName = ModuleList[[i]]$Name, ModuleDesc = ModuleList[[i]]$Desc,
                               UsedGenes = SelGenes, SampledGenes = SampledsGeneList, PCABase = PCBase,
                               ExpVarBase = ExpVar, PVVect = PVVect, SampledExp = SampledExp)
  
    
  }
  
  colnames(ModuleMatrix) <- c("L1", "ppv L1", "L1/L2", "ppv L1/L2")
  colnames(PVVectMat) <- c("L1 WT less pv", "L1 WT greater pv", "L1/L2 WT less pv", "L1/2 WT greater pv")

  colnames(PC1Matrix) <- colnames(ExpressionMatrix)
  
  return(list(ModuleMatrix = ModuleMatrix, PC1Matrix = PC1Matrix, ModuleSummary = ModuleSummary,
              ProjLists = ProjLists, PVVectMat = PVVectMat, OutLiersList = OutLiersList))
  
}











#' Produce a list of genesets with their overdispersion
#'
#' @param GeneExpressionMatrix 
#' @param nSamples 
#' @param nGenes 
#' @param GeneOutDetection 
#' @param GeneOutThr 
#' @param CenterPCA 
#' @param PlotData 
#' @param ModuleName 
#' @param PrintInfo 
#'
#' @return
#' @export
#'
#' @examples
SelectOverDispersedGeneSet <- function(GeneExpressionMatrix, nSamples, nGenes, GeneOutDetection = 'none', GeneOutThr = 1,
                                       CenterPCA = FALSE, PlotData = FALSE, ModuleName = '', PrintInfo = FALSE, PCADims = 2) {
  
  if(PCADims > 2){
    PCADims = 2
  }
  
  SampledsGeneList <- lapply(as.list(1:nSamples), function(i){sample(x = rownames(GeneExpressionMatrix), size = nGenes, replace = FALSE)})
  
  pb <- txtProgressBar(min = 0, max = nSamples, initial = 0, style = 3)

  SampledExp <- sapply(as.list(1:length(SampledsGeneList)), function(i){
    if(SampleFilter){
      SelGenes <- DetectOutliers(GeneOutDetection = GeneOutDetection, GeneOutThr = GeneOutThr, ModulePCACenter = ModulePCACenter,
                                 CompatibleGenes = SampledsGeneList[[i]], ExpressionData = ExpressionMatrix[SampledsGeneList[[i]], ],
                                 PlotData = FALSE, ModuleName = '', PrintInfo = FALSE)
      setTxtProgressBar(pb, i)
    } else {
      SelGenes <- SampledsGeneList[[i]]
    }
    PCSamp <- irlba::prcomp_irlba(x = ExpressionMatrix[SelGenes, ], n = PCADims, center = CenterPCA, scale. = FALSE)
    return((PCSamp$sdev^2)/sum(apply(scale(ExpressionData[SelGenes, ], center = ModulePCACenter, scale = FALSE), 2, var)))
  })
  
  SampledExp <- rbind(SampledExp[1,], SampledExp[1,]/SampledExp[2,])
  
  return(list(SampledsGenes = SampledsGeneList[order(SampledExp[1,], decreasing = TRUE)],
              SampledExp = SampledExp[, order(SampledExp[1,], decreasing = TRUE)]))
  
}


