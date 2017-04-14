#' Title
#'
#' @param RomaData list, the analysis returned by rRoma
#' @param Groups string vector, a vector of group identifier described by strings
#' @param Selected The genesets used to perform the analysis 
#' @param TestMode string, the type of statistical methodology to assess sample difference. Currentyl only ANOVA + Tukey is implemented ("Aov+Tuk")
#' @param TestPV1 numeric between 0 and 1, the threshold PV for the initial test (ANOVA)
#' @param TestPV2 numeric between 0 and 1, the threshold PV for difference test (Tukey)
#' @param PlotDiag boolean, should diagnostic plot be displyed?
#' @param PlotXGSDiff boolean, should all the significant differences by considered?
#'
#' @return
#' @export
#'
#' @examples
CompareAcrossSamples <- function(RomaData, Groups, Selected = NULL,
                                 TestMode = "Aov+Tuk", TestPV1 = 5e-2, TestPV2 = 5e-2,
                                 PlotDiag = FALSE, PlotXGSDiff = FALSE) {
  
  if(is.null(Selected)){
    Selected <- 1:nrow(RomaData$ProjMatrix)
  }
  
  if(length(intersect(Selected, 1:nrow(RomaData$ProjMatrix)))<1){
    print("No Genset selected")
    return(NULL)
  } else {
    print(paste(length(intersect(Selected, 1:nrow(RomaData$ProjMatrix))), "geneset selected"))
  }
  
  tMat <- RomaData$ProjMatrix[Selected,]
  names(Groups) <- colnames(tMat)
  
  MeltData <- reshape::melt(tMat)
  MeltData <- cbind(MeltData, Groups[as.character(MeltData$X2)])
  
  
  # MeltData <- data.frame(MeltData)
  
  colnames(MeltData) <- c("GeneSet", "Sample", "Value", "Group") 
  
  if(TestMode == "Aov+Tuk"){
    
    print("Performing Type III AOV (R default)")
    
    AOVFitTypeI <- aov(formula = Value ~ Group, data = MeltData)
    
    if(PlotDiag){
      plot(AOVFitTypeI)
    }
    
    SumData <- summary(AOVFitTypeI)
    print(SumData)
    
    
    p <- ggplot2::ggplot(MeltData, ggplot2::aes(y=Value, x=Group, fill=Group)) +
      ggplot2::geom_boxplot() + ggplot2::guides(fill = "none") +
      ggplot2::labs(y="PC1 projection", x="Groups", title = "Groups")
    
    print(p)
    
    p <- ggplot2::ggplot(MeltData, ggplot2::aes(y=Value, x=Sample, fill=Group)) +
      ggplot2::geom_boxplot() + ggplot2::coord_flip() +
      ggplot2::labs(y="PC1 projection", x="Samples", title = "Groups") +
      ggplot2::theme(axis.text.y = ggplot2::element_blank())
        
    print(p)
    
    
    if(SumData[[1]][1,5] < TestPV1){
      print("A significant difference is observed across groups")
      
      print("Calculating Tukey Honest Significant Differences")
      
      TukTest <- TukeyHSD(AOVFitTypeI, conf.level = 1-TestPV2)
      
      if(PlotDiag){
        plot(TukTest)
      }
      
      Diffs <- TukTest$Group
      Diffs <- data.frame(Diffs[Diffs[,4] < TestPV2,])
      Diffs <- cbind(rownames(Diffs), Diffs)
      colnames(Diffs)[1] <- "GDiff"
      
      print(paste(nrow(Diffs), "significant differences found"))
      
      if(nrow(Diffs)>1){
        
        Sep <- seq(from = 1, to = nrow(Diffs), by = 10)
        if(max(Sep) < nrow(Diffs)) Sep <- c(Sep, nrow(Diffs))
        
        for(i in 2:length(Sep)){
          p <- ggplot2::ggplot(Diffs[Sep[i-1]:Sep[i],],
                               ggplot2::aes(x = GDiff, y = diff, ymin = lwr, ymax = upr)) +
            ggplot2::geom_point() + ggplot2::geom_errorbar(width=0.25) + ggplot2::coord_flip() +
            ggplot2::facet_wrap( ~ GDiff, scales = "free_y", ncol = 1) +
            ggplot2::labs(y="Mean difference", x="Categories", title = paste("Groups - Part", i-1)) +
            ggplot2::theme(
              # axis.text.x = element_blank(),
              axis.text.y = ggplot2::element_blank(),
              axis.ticks = ggplot2::element_blank(),
              strip.text.x = ggplot2::element_text(size=6, face = "bold"))
          
          print(p)
          
        }
        
      }
      
      
      
    }
    
    print("Performing Type III AOV (R default)")
    
    AOVFitTypeI <- aov(formula = Value ~ Group/GeneSet, data = MeltData)
    
    if(PlotDiag){
      plot(AOVFitTypeI)
    }
    
    SumData <- summary(AOVFitTypeI)
    print(SumData)
    
    Sep <- seq(from = 1, to = length(unique(MeltData$GeneSet)), by = 4)
    if(max(Sep) < length(unique(MeltData$GeneSet))+1) Sep <- c(Sep, length(unique(MeltData$GeneSet))+1)
    
    for(i in 2:length(Sep)){
      
      p <- ggplot2::ggplot(MeltData[as.integer(MeltData$GeneSet) %in% Sep[i-1]:(Sep[i]-1),],
                           ggplot2::aes(y=Value, x=Group, fill=Group)) + ggplot2::geom_boxplot() +
        ggplot2::labs(y="PC1 projection", x="Groups", title = paste("Geneset VS Groups - Part", i-1)) +
        ggplot2::facet_wrap( ~ GeneSet, ncol = 2) + ggplot2::theme(strip.text.x = ggplot2::element_text(size=6, face = "bold")) +
        ggplot2::guides(fill = "none")
      
      print(p)
    }
    
    if(SumData[[1]][2,5] < TestPV1){
      print("A significant difference is observed across groups and metagenes")
      
      print("Calculating Tukey Honest Significant Differences")
      
      TukTest <- TukeyHSD(AOVFitTypeI, conf.level = 1-TestPV2)
      
      if(PlotDiag){
        plot(TukTest, las=2)
      }
      
      Diffs <- TukTest$`Group:GeneSet`
      Diffs <- data.frame(Diffs[Diffs[,4] < TestPV2,])
      Diffs <- cbind(rownames(Diffs), Diffs)
      colnames(Diffs)[1] <- "GGDiff"
      
      GSPairs <- lapply(strsplit(gsub(":", "-", as.character(Diffs$GGDiff)), c("-")),"[", c(2,4))
      SameGS <- unlist(lapply(lapply(GSPairs,duplicated),any))
      GSVect <- unlist(lapply(GSPairs[SameGS], "[[", 1), use.names = FALSE)
      
      print(paste(nrow(Diffs), "significant differences found"))
      
      if(nrow(Diffs)>0){
        
        if(PlotXGSDiff){
          
          Sep <- seq(from = 1, to = nrow(Diffs), by = 10)
          if(max(Sep) < nrow(Diffs)) Sep <- c(Sep, nrow(Diffs))
          
          
          for(i in 2:length(Sep)){
            p <- ggplot2::ggplot(Diffs[Sep[i-1]:Sep[i],],
                                 ggplot2::aes(x = GGDiff, y = diff, ymin = lwr, ymax = upr)) +
              ggplot2::geom_point() + ggplot2::geom_errorbar(width=0.25) + ggplot2::coord_flip() +
              ggplot2::facet_wrap( ~ GGDiff, scales = "free_y", ncol = 1) +
              ggplot2::labs(y="Mean difference", x="Categories", title = paste("Groups VS Genesets - Part", i-1)) +
              ggplot2::theme(
                # axis.text.x = element_blank(),
                axis.text.y = ggplot2::element_blank(),
                axis.ticks = ggplot2::element_blank(),
                strip.text.x = ggplot2::element_text(size=6, face = "bold"))
            
            print(p)
          }
          
        }
        
        SplitDiff <- split(Diffs[SameGS,], GSVect)
        
        for(i in 1:length(SplitDiff)){
          
          if(nrow(SplitDiff[[i]])<1){
            next
          }
          
          p <- ggplot2::ggplot(SplitDiff[[i]], ggplot2::aes(x = GGDiff, y = diff, ymin = lwr, ymax = upr)) +
            ggplot2::geom_point() + ggplot2::geom_errorbar(width=0.25) + ggplot2::coord_flip() +
            ggplot2::facet_wrap( ~ GGDiff, scales = "free_y", ncol = 1) +
            ggplot2::labs(y="Mean difference", x="Categories",
                          title = paste("Groups within", names(SplitDiff)[i])) +
            ggplot2::theme(
              # axis.text.x = element_blank(),
              axis.text.y = ggplot2::element_blank(),
              axis.ticks = ggplot2::element_blank(),
              strip.text.x = ggplot2::element_text(size=6, face = "bold"))
          
          print(p)
        }
        
        
        
      }
      
      
      
      # GSAOV <- lapply(split(MeltData, MeltData$GeneSet), function(DF){
      #   aov(formula = Value ~ Group, data = DF)
      # })
      
      # GSAOPVs <- lapply(lapply(lapply(lapply(GSAOV, summary), "[[", 1), "[[", 5), "[", 1)

      # barplot(unlist(GSAOPVs), log = "y", las=2)
      
      # TukList <- lapply(GSAOV, TukeyHSD, conf.level = 1-TestPV2)
      
      # PVPairsData <- cbind(t(matrix(unlist(strsplit(rownames(TukTest$Type), "-")), nrow = 2)), TukTest$Type[,4])
      # PVPairsMatrix <- matrix(rep(NA, length(unique(as.vector(PVPairsData[,1:2])))^2), length(unique(as.vector(PVPairsData[,1:2])))) 
      # colnames(PVPairsMatrix) <- unique(as.vector(PVPairsData[,1:2]))
      # rownames(PVPairsMatrix) <- unique(as.vector(PVPairsData[,1:2]))
      # 
      # apply(PVPairsData, 1, function(Vect) {
      #   PVPairsMatrix[Vect[1], Vect[2]] <<- as.numeric(Vect[3])
      #   PVPairsMatrix[Vect[2], Vect[1]] <<- as.numeric(Vect[3])
      # })
      # 
      
    }
    
  }

}










#' Select genesets accoding to specific conditions
#'
#' @param RomaData list, the analysis returned by rRoma 
#' @param VarThr numeric between 0 and 1, the threshold PV to select significantly over- or under-disperdes genesets
#' @param VarMode string, the test to use to select over- or under-underdisperdes genesets.
#' Currently it can be either 'Wil' (Wilcoxon test) or 'PPV' (permutation base p-value)
#' @param VarType string, the type of statistical difference to select. Currently it can be either 'Over' (overdispersed) or 'Under' (underdispersed)
#' @param MedThr numeric between 0 and 1, the threshold PV to select significantly over- or under-expressed genesets
#' @param MedMode string, the test to use to select over- or under-expressed genesets.
#' Currently it can be either 'Wil' (Wilcoxon test) or 'PPV' (permutation base p-value)
#' @param MedType string, the type of statistical difference to select. Currently it can be either 'Over' (overexpressed) or 'Under' (underexpressed)
#'
#' @return
#' @export
#'
#' @examples
SelectGeneSets <- function(RomaData, VarThr = 1e-3, VarMode = "Wil", VarType = "Over",
                           MedThr = NULL, MedMode = "Wil", MedType = "Over") {
  
  while(1){ 
    
    if(VarMode == 'Wil' & VarType == "Over"){
      print(paste("Using genestes overdispersed according to Wilcoxon test. VarThr =", VarThr))
      Selected <- which(RomaData$PVVectMat[,1]<VarThr)
      break
    }
    
    if(VarMode == 'Wil' & VarType == "Under"){
      print(paste("Using genestes underdispersed according to Wilcoxon test. VarThr =", VarThr))
      Selected <- which(RomaData$PVVectMat[,2]<VarThr)
      break
    }
    
    if(VarMode == 'PPV' & VarType == "Over"){
      print(paste("Using genestes overdispersed according to pseudo pv. VarThr =", VarThr))
      Selected <- which(RomaData$ModuleMatrix[,2]<VarThr)
      break
    }
    
    if(VarMode == 'PPV' & VarType == "Under"){
      print(paste("Using genestes underdispersed according to pseudo pv. VarThr =", VarThr))
      Selected <- which(RomaData$ModuleMatrix[,2]<1-VarThr)
      break
    }
    
    Selected <- NULL
    break
  }
  
  while(!is.null(Selected) & !is.null(MedThr)){
    
    if(MedMode == 'Wil' & MedType == "Over"){
      print(paste("Using genestes overexpressed according to Wilcoxon test. MedThr =", MedThr))
      Selected <- intersect(which(RomaData$PVVectMat[,5]<MedThr), Selected)
      break
    }
    
    if(MedMode == 'Wil' & MedType == "Under"){
      print(paste("Using genestes underexpressed according to Wilcoxon test. MedThr =", MedThr))
      Selected <- intersect(which(RomaData$PVVectMat[,6]<MedThr), Selected)
      break
    }
    
    if(MedMode == 'PPV' & MedType == "Over"){
      print(paste("Using genestes overexpressed according to pseudo pv. MedThr =", MedThr))
      Selected <- intersect(which(RomaData$ModuleMatrix[,6]<MedThr), Selected)
      break
    }
    
    if(MedMode == 'PPV' & MedType == "Under"){
      print(paste("Using genestes underexpressed according to pseudo pv. MedThr =", MedThr))
      Selected <- intersect(which(RomaData$ModuleMatrix[,6]<1-MedThr), Selected)
      break
    }
    
    Selected <- NULL
    break
    
  }
  
  return(Selected)
  
}
