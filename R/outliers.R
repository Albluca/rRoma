#' Detect outliers
#'
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
#' @param ModulePCACenter 
#' @param Genes 
#' @param ExpressionData 
#' @param PlotData 
#' @param ModuleName
#' @param PCAType character string, the type of PCA to perform. It can be "DimensionsAreGenes" or "DimensionsAreSamples"
#'
#' @return
#' @export
#'
#' @examples
DetectOutliers <- function(GeneOutDetection, GeneOutThr, ModulePCACenter,
                           CompatibleGenes, ExpressionData, PCAType = PCAType,
                           PlotData = FALSE, ModuleName = '', PrintInfo = TRUE,
                           Mode = 1) {
  
  if(!(PCAType %in% c("DimensionsAreGenes", "DimensionsAreSamples"))){
    print("Incompatible PCAType, no outlier filtering will be performed")
    return(CompatibleGenes)
  }
  
  SelGenes <- CompatibleGenes
  
  GetAllPC1Var <- function(i){
    
    if(PCAType == "DimensionsAreGenes"){
      tData <- t(ExpressionData[-i, ])
    }
    
    if(PCAType == "DimensionsAreSamples"){
      tData <- ExpressionData[-i, ]
    }
    
    PC1Var <- var(
      irlba::prcomp_irlba(x = tData, n = 1, center = ModulePCACenter,
                          scale. = FALSE, retx = TRUE)$x
    )
    return(PC1Var/sum(apply(scale(tData, center = ModulePCACenter, scale = FALSE), 2, var)))
    
  }
  
  # Computing all the PC1
  if(Mode == 1){
    AllPCA1 <- sapply(as.list(1:length(CompatibleGenes)), GetAllPC1Var)
  } else {
    AllPCA1 <- rep(0, length(CompatibleGenes))
    for(i in 1:length(CompatibleGenes)){
      AllPCA1[i] <- GetAllPC1Var(i)
    }
  }
  
  
  
  if(GeneOutDetection == "L1OutVarPerc"){
    
    if(PrintInfo){
      print("Detecting outliers using leave one out and percentage variation on variance explained by PC1")
    }
    
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
    
    if(PrintInfo){
      print("Detecting outliers using leave one out and dendrogram analysis on variance explained by PC1")
    }
    
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
        GenesOut <- (DGroup == 1)
      } else {
        GenesOut <- (DGroup == 2)
      }
      
    }
    
    if(PrintInfo){
      if(sum(GenesOut)>0){
        print(paste(sum(GenesOut), "gene(s) will be filtered:"))
        print(CompatibleGenes[GenesOut])
      } else {
        print("No gene will be filtered")
      }
    }
    
    # Updating gene list
    SelGenes <- CompatibleGenes[!GenesOut]
    
  }
  
  if(GeneOutDetection == "L1OutExpOut"){
    
    if(PrintInfo){
      print("Detecting outliers using leave one out and median-absolute-deviations away from median (scater package)")
    }
    
    # Getting the distance from the median PC1
    GenesOut <- scater::isOutlier(AllPCA1, nmads = GeneOutThr)
    
    if(PlotData){
      # Plotting the distances from the median PC1
      B <- boxplot(x = AllPCA1, at = 1, horizontal = FALSE, ylab = "Variance explained by PC1", main = ModuleName)
    }
    
    if(PrintInfo){
      if(sum(GenesOut)>0){
        print(paste(sum(GenesOut), "gene(s) will be filtered:"))
        print(CompatibleGenes[GenesOut])
      } else {
        print("No gene will be filtered")
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
  
  if(GeneOutDetection == "L1OutSdMean"){
    
    if(PrintInfo){
      print("Detecting outliers using leave one out and standards deviation away from the mean")
    }
    
    # Getting the distance from the median PC1
    ZScore <- (AllPCA1 - mean(AllPCA1))/sd(AllPCA1)
    GenesOut <- abs(ZScore) > GeneOutThr
    
    if(PlotData){
      # Plotting the distances from the median PC1
      B <- boxplot(x = AllPCA1, at = 1, horizontal = FALSE, ylab = "Variance explained by PC1", main = ModuleName)
    }
    
    if(PrintInfo){
      if(sum(GenesOut)>0){
        print(paste(sum(GenesOut), "gene(s) will be filtered:"))
        print(CompatibleGenes[GenesOut])
      } else {
        print("No gene will be filtered")
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
  
  return(SelGenes)
  
}
