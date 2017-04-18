#' Title
#'
#' @param RomaData 
#' @param SampleName 
#' @param AggSampFun 
#' @param Selected 
#' @param FilterByWei 
#' @param AggGeneFun 
#' @param DispMode 
#' @param DataName 
#'
#' @return
#' @export
#'
#' @examples
PlotOnACSN <- function(RomaData, SampleName, AggScoreFun = "mean",
                       QTop = NULL, QBottom = NULL, Steps = 2,
                       MapURL = "", Selected = NULL, FilterByWei = NULL,
                       AggGeneFun = "mean", DispMode = "Module",
                       DataName = "Sample") {
  
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
  
  if(!is.null(FilterByWei)){
    boxplot(AllGenesWei.Var, ylab = "Variance of gene weight")
    Outliers <- scater::isOutlier(AllGenesWei.Var, nmads = FilterByWei, na.rm = TRUE)
    print("The following genes will be ignored:")
    print(names(which(Outliers)))
    ToPlot <- setdiff(ToPlot, names(which(Outliers)))
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
   
    AllGeneScores <- unlist(AllGeneScores, use.names = TRUE)
    AllGeneScores.Split <- split(AllGeneScores, names(AllGeneScores))
    DataToPlot <- sapply(AllGeneScores.Split, get(AggGeneFun))
    
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
