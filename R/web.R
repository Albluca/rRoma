#' Project ROMA results onto ACSN maps
#'
#' @param RomaData list, the analysis returned by rRoma
#' @param Selected vector, integer. The position of the genesets to plot
#' @param SampleName vector, string. The name(s) of the sample to consider
#' @param FilterByWei scalar, number of median-absolute-deviations away from median variace of gene weights
#' required for a value to be called an outlier and filtered. If NULL genes are not filtered. 
#' @param AggGeneFun string, scalar. The name of the function used to combine scored if multiple samples are considered 
#' @param DispMode string, vector. The quantity to be projected on the maps (in the form of an heatmap).
#' It can be any combination of the following:
#' \itemize{
#' \item 'Module': gene score is obtained by aggregating the score of selected modules containing that gene. The aggregating function is specified by the AggGeneFun parameter. 
#' \item 'Gene': gene score is obtained by aggregating the weigth of the genes across the selected modules. The aggregating function is specified by the AggGeneFun parameter. 
#' }
#' @param DataName string, scalar. The name of the dataset that will be reported in the web interface
#' @param AggScoreFun string, scalar. The name of the aggregating function used to compute the valus of the gene score.
#' Note that this function will be applied even to genes appearin in one geneset. To consider only genes found in a single module use the 'GetSingle' function and set DefDispVal to 0
#' @param DefDispVal numeric, scalar. The value to be used when the fucntion defined by AggScoreFun produce NA (e.g., sd applied to genes that appear only in one geneset)
#' @param QTop numeric, scalar. The quantile threshold associted with the top saturation color. If NULL, the max will be used.
#' @param QBottom numeric, scalar. The quantile threshold associted with the lower saturation color.  If NULL, the min will be used.
#' @param Steps numeric, scalar. The number of steps of the heatmap?? 
#' @param MapURL string, scalar. The url of the map to use. The available maps can be seen on acsn.curie.fr and navicell.curie.fr. Possible values include
#' \itemize{
#' \item 'https://navicell.curie.fr/navicell/maps/cellcycle/master/index.php'
#' \item 'https://acsn.curie.fr/navicell/maps/dnarepair/master/index.php'
#' \item 'https://acsn.curie.fr/navicell/maps/emtcellmotility/master/index.php'
#' \item 'https://navicell.curie.fr/navicell/maps/ewing/master/index.php'
#' \item 'https://navicell.curie.fr/navicell/maps/dendcells/master/index.php'
#' } 
#' @param PlotInfo boolean scalar. Should information plot be produced?
#' @param ReturnInfo boolean scalar. Should information be returned?
#' @param LocalRange boolean scalar. Should color bounds be relative to the considered samples only?
#'
#'
#' @return
#' @export
#'
#' @examples
PlotOnACSN <- function(RomaData, SampleName, AggScoreFun = "mean",
                       QTop = NULL, QBottom = NULL, Steps = 2,
                       MapURL = "", Selected = NULL, FilterByWei = NULL,
                       AggGeneFun = "mean", DispMode = "Module",
                       DataName = "Sample", DefDispVal = 0, PlotInfo = TRUE,
                       ReturnInfo = FALSE, LocalRange = FALSE) {
  
  if(!requireNamespace("RNaviCell", quietly = TRUE)){
    stop("Unable to load RNaviCell Impossible to proceed")
  }
  
  if(is.null(Selected)){
    Selected <- 1:nrow(RomaData$ProjMatrix)
  }
  
  if(sum(SampleName %in% colnames(RomaData$ProjMatrix))<0){
    print("Sample(s) not found")
    return()
  } else {
    SelSampleNames <- intersect(SampleName, colnames(RomaData$ProjMatrix))
    print(paste(length(SelSampleNames), "sample(s) selected"))
  }
  
  if(length(intersect(Selected, 1:nrow(RomaData$ProjMatrix)))<1){
    print("No Genset selected")
    return(NULL)
  } else {
    print(paste(length(intersect(Selected, 1:nrow(RomaData$ProjMatrix))), "geneset(s) selected"))
  }
  
  
  AllGenesWei <- lapply(RomaData$ModuleSummary[Selected], "[[", "PC1Weight.SignFixed")
  AllGenesWei <- unlist(AllGenesWei, use.names = TRUE)
  AllGenesWei.Split <- split(AllGenesWei, f = names(AllGenesWei))
  
  AllGenesWei.Var <- sapply(AllGenesWei.Split, var)
  
  ToPlot <- unique(names(AllGenesWei))
  
  if(any(!is.na(AllGenesWei.Var))){
    
    if(PlotInfo){
      B <- boxplot(x = AllGenesWei.Var, at = 1, ylab = "Variance of gene weight")
    }
    
    Outliers <- NULL
    
    if(!is.null(FilterByWei)){
      Outliers <- scater::isOutlier(AllGenesWei.Var[!is.na(AllGenesWei.Var)], nmads = FilterByWei)
      print("The following genes will be ignored:")
      print(names(which(Outliers)))
      ToPlot <- setdiff(ToPlot, names(which(Outliers)))
      if(PlotInfo){
        points(x = rep(1, length(which(Outliers))), y = AllGenesWei.Var[names(which(Outliers))], col="red", pch=20)
        legend("topright", legend = "Outliers", pch=20, col="red")
      }
    }
    
  }

  AllGeneScores.Var <- NULL
  AllGeneScores.Agg <- NULL
  AllGene.Mult <- table(names(AllGenesWei))
  
  if(any(DispMode == "Module")){
    
    ModuleScore <- RomaData$ProjMatrix[Selected, SelSampleNames]
    
    if(length(SelSampleNames)>1){
      if(length(Selected)>1){
        FilteredScores <- apply(ModuleScore, 1, get(AggScoreFun))
      } else {
        FilteredScores <- do.call(AggScoreFun, list(ModuleScore))
        names(FilteredScores) <- RomaData$ModuleSummary[[Selected]]$ModuleName
      }
    } else {
      print("Single sampled selected. Aggregating function will not be applied.")
      FilteredScores <- ModuleScore
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
    
    if(PlotInfo){
      barplot(c(sum(!is.na(AllGeneScores.Var)), sum(is.na(AllGeneScores.Var))),
              names.arg = c("Multi", "Mono"), ylab = "Number of genes")
    }
    
    if(any(!is.na(AllGeneScores.Var))){
      boxplot(AllGeneScores.Var, ylab = "Variance of module score (per gene)")
    }
    
    DataToPlot <- sapply(AllGeneScores.Split, get(AggGeneFun))
    
    if(any(is.na(DataToPlot))){
      DataToPlot[is.na(DataToPlot)] <- DefDispVal
    }
   
    DataToPlot <- data.matrix(DataToPlot)
    colnames(DataToPlot) <- "data"
    
    AllGeneScores.Agg <- DataToPlot
    
    if(PlotInfo){
      hist(DataToPlot, xlab="Score (per gene)", main = "Module score distribution across genes")
    }
    
    Session.Navicell <- RNaviCell::NaviCell(SetUrl = MapURL)
    
    Session.Navicell$launchBrowser()
    
    Session.Navicell$importDatatable("Continuous copy number data", DataName, DataToPlot)
    
    Session.Navicell$continuousConfigSwitchSampleTab(DataName, "color")
    Session.Navicell$continuousConfigSetStepCount("sample", 'color', DataName, Steps)
    Session.Navicell$continuousConfigSetColorAt(DataName, "sample", 0, 'FFFFFF')
    
    if(LocalRange){
      
      if(is.null(QTop)){
        Top <- max(DataToPlot)
      } else {
        Top <- quantile(DataToPlot, QTop)
      }
      
      if(is.null(QBottom)){
        Bottom <- min(DataToPlot)
      } else {
        Bottom <- quantile(DataToPlot, QBottom)
      }
      
    } else {
      
      if(is.null(QTop)){
        Top <- max(RomaData$ProjMatrix[Selected,])
      } else {
        Top <- quantile(RomaData$ProjMatrix[Selected,], QTop)
      }
      
      if(is.null(QBottom)){
        Bottom <- min(RomaData$ProjMatrix[Selected,])
      } else {
        Bottom <- quantile(RomaData$ProjMatrix[Selected,], QBottom)
      }
      
    }
    
    
    Session.Navicell$continuousConfigSetValueAt(DataName, "color", "sample", Bottom, -1)
    Session.Navicell$continuousConfigSetValueAt(DataName, "color", "sample", Top, 1)
    Session.Navicell$continuousConfigApply(DataName, "color")
    
    Session.Navicell$mapStainingEditorSelectDatatable(DataName)
    Session.Navicell$mapStainingEditorSelectSample('data')
    Session.Navicell$mapStainingEditorApply()
    
  }
  
  AllGeneWei.Agg <- NULL
  
  par(mfcol=c(1,1))
  
  if(any(DispMode == "Gene")){
    
    DataToPlot <- sapply(AllGenesWei.Split[ToPlot], get(AggGeneFun))
    
    DataToPlot <- data.matrix(DataToPlot)
    colnames(DataToPlot) <- "data"
    
    AllGeneWei.Agg <- DataToPlot
    
    if(PlotInfo){
      hist(DataToPlot, xlab="Weight (per gene)", main = "Weight distribution across genes")
    }
    
    Session.Navicell <- RNaviCell::NaviCell(SetUrl = MapURL)
    
    Session.Navicell$launchBrowser()
    
    Session.Navicell$importDatatable("Continuous copy number data", DataName, DataToPlot)
    
    Session.Navicell$continuousConfigSwitchSampleTab(DataName, "color")
    Session.Navicell$continuousConfigSetStepCount("sample", 'color', DataName, Steps)
    Session.Navicell$continuousConfigSetColorAt(DataName, "sample", 0, 'FFFFFF')
    
    
    if(is.null(QTop)){
      Top <- max(DataToPlot)
    } else {
      Top <- quantile(DataToPlot, QTop)
    }
    
    if(is.null(QBottom)){
      Bottom <- min(DataToPlot)
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
  
  if(ReturnInfo){
    return(list(GenesVar = AllGenesWei.Var, GeneOut = Outliers,
                ScoreVar = AllGeneScores.Var, GeneMult = AllGene.Mult,
                WeiDist = AllGeneWei.Agg, ScoreDist = AllGeneScores.Agg)
           )
  }
  
}






#' Return its input if it is a single value, else return NA
#'
#' @param x input vector
#'
#' @return
#' @export
#'
#' @examples
GetSingle <- function(x) {
  ifelse(length(x)==1, x, NA)
}
