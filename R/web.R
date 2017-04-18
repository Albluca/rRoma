#' Project ROMA results onto ACSN maps
#'
#' @param RomaData list, the analysis returned by rRoma
#' @param Selected vector, integer. The position of the genesets to plot
#' @param SampleName 
#' @param FilterByWei 
#' @param AggGeneFun 
#' @param DispMode 
#' @param DataName 
#' @param AggScoreFun 
#' @param QTop 
#' @param QBottom 
#' @param Steps 
#' @param MapURL 
#' @param DefDispVal 
#'
#' @return
#' @export
#'
#' @examples
PlotOnACSN <- function(RomaData, SampleName, AggScoreFun = "mean",
                       QTop = NULL, QBottom = NULL, Steps = 2,
                       MapURL = "", Selected = NULL, FilterByWei = NULL,
                       AggGeneFun = "mean", DispMode = "Module",
                       DataName = "Sample", DefDispVal = 0) {
  
  if(is.null(Selected)){
    Selected <- 1:nrow(RomaData$ProjMatrix)
  }
  
  if(sum(SampleName %in% colnames(RomaData$ProjMatrix))<0){
    print("Sample(s) not found")
    return()
  } else {
    SelSampleNames <- intersect(SampleName, colnames(RomaData$ProjMatrix))
  }
  
  if(length(intersect(Selected, 1:nrow(RomaData$ProjMatrix)))<1){
    print("No Genset selected")
    return(NULL)
  } else {
    print(paste(length(intersect(Selected, 1:nrow(RomaData$ProjMatrix))), "geneset selected"))
  }
  
  
  AllGenesWei <- lapply(RomaData$ModuleSummary, "[[", "PC1Weight.SignFixed")
  AllGenesWei <- unlist(AllGenesWei, use.names = TRUE)
  AllGenesWei.Split <- split(AllGenesWei, f = names(AllGenesWei))
  
  AllGenesWei.Var <- sapply(AllGenesWei.Split, var)
  
  ToPlot <- unique(names(AllGenesWei))
  
  par(mfcol=c(1,3))
  
  B <- boxplot(x = AllGenesWei.Var, at = 1, ylab = "Variance of gene weight")
  
  if(!is.null(FilterByWei)){
    Outliers <- scater::isOutlier(AllGenesWei.Var, nmads = FilterByWei)
    print("The following genes will be ignored:")
    print(names(which(Outliers)))
    ToPlot <- setdiff(ToPlot, names(which(Outliers)))
    points(x = rep(1, length(which(Outliers))), y = AllGenesWei.Var[names(which(Outliers))], col="red", pch=20)
    legend("topright", legend = "Outliers", pch=20, col="red")
  }

  if(any(DispMode == "Module")){
    ModuleScore <- RomaData$ProjMatrix[Selected, SelSampleNames]
    
    if(length(SelSampleNames)>1){
      FilteredScores <- apply(ModuleScore, 1, get(AggScoreFun))
    }
    
    AllGeneScores <- lapply(as.list(Selected), function(i){
      GenesNames <- RomaData$ModuleSummary[[i]]$UsedGenes
      GeneScores <- rep(FilteredScores[RomaData$ModuleSummary[[i]]$ModuleName], length(GenesNames))
      names(GeneScores) <- GenesNames
      GeneScores
    })
   
    names(AllGeneScores) <- NULL
    
    AllGeneScores <- unlist(AllGeneScores, use.names = TRUE)
    AllGeneScores.Split <- split(AllGeneScores, names(AllGeneScores))
    AllGeneScores.Var <- sapply(AllGeneScores.Split, var)
    
    barplot(table(is.na(AllGeneScores.Var)), names.arg = c("Multi", "Mono"),
            ylab = "Number of genes")
    boxplot(AllGeneScores.Var, ylab = "Variance of module score (per gene)")
    
    DataToPlot <- sapply(AllGeneScores.Split, get(AggGeneFun))
    
    if(any(is.na(DataToPlot))){
      DataToPlot[is.na(DataToPlot)] <- DefDispVal
    }
   
    DataToPlot <- data.matrix(DataToPlot)
    colnames(DataToPlot) <- "data"
    
    hist(DataToPlot)
    
    Session.Navicell <- RNaviCell::NaviCell(SetUrl = MapURL)
    
    Session.Navicell$launchBrowser()
    
    Session.Navicell$importDatatable("Continuous copy number data", DataName, DataToPlot)
    
    Session.Navicell$continuousConfigSwitchSampleTab(DataName, "color")
    Session.Navicell$continuousConfigSetStepCount("sample", 'color', DataName, Steps)
    Session.Navicell$continuousConfigSetColorAt(DataName, "sample", 0, 'FFFFFF')
    
    if(is.null(QTop)){
      Top <- quantile(DataToPlot, .95)
    } else {
      Top <- quantile(DataToPlot, QTop)
    }
    
    if(is.null(QBottom)){
      Bottom <- quantile(DataToPlot, .5)
    } else {
      Bottom <- quantile(DataToPlot, QBottom)
    }
    
    Session.Navicell$continuousConfigSetValueAt(DataName, "color", "sample", Bottom, -1)
    Session.Navicell$continuousConfigSetValueAt(DataName, "color", "sample", Top, 1)
    Session.Navicell$continuousConfigApply(DataName, "color")
    
    Session.Navicell$mapStainingEditorSelectDatatable(DataName)
    Session.Navicell$mapStainingEditorSelectSample('data')
    Session.Navicell$mapStainingEditorApply()
    
  }
  
  par(mfcol=c(1,1))
  
  if(any(DispMode == "Gene")){
    
    DataToPlot <- sapply(AllGenesWei.Split[ToPlot], get(AggGeneFun))
    
    DataToPlot <- data.matrix(DataToPlot)
    colnames(DataToPlot) <- "data"
    
    hist(DataToPlot)
    
    Session.Navicell <- RNaviCell::NaviCell(SetUrl = MapURL)
    
    Session.Navicell$launchBrowser()
    
    Session.Navicell$importDatatable("Continuous copy number data", DataName, DataToPlot)
    
    Session.Navicell$continuousConfigSwitchSampleTab(DataName, "color")
    Session.Navicell$continuousConfigSetStepCount("sample", 'color', DataName, Steps)
    Session.Navicell$continuousConfigSetColorAt(DataName, "sample", 0, 'FFFFFF')
    
    if(is.null(QTop)){
      Top <- quantile(DataToPlot, .95)
    } else {
      Top <- quantile(DataToPlot, QTop)
    }
    
    if(is.null(QBottom)){
      Bottom <- quantile(DataToPlot, .5)
    } else {
      Bottom <- quantile(DataToPlot, QBottom)
    }
    
    Session.Navicell$continuousConfigSetValueAt(DataName, "color", "sample", Bottom, -1)
    Session.Navicell$continuousConfigSetValueAt(DataName, "color", "sample", Top, 1)
    Session.Navicell$continuousConfigApply(DataName, "color")
    
    Session.Navicell$mapStainingEditorSelectDatatable(DataName)
    Session.Navicell$mapStainingEditorSelectSample('data')
    Session.Navicell$mapStainingEditorApply()
    
  }
  
}
