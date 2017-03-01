#' Title
#'
#' @param RomaData 
#' @param Groups 
#' @param GSThr 
#' @param Mode 
#' @param Type 
#' @param TestMode 
#' @param TestPV1 
#' @param TestPV2 
#' @param PlotDiag 
#'
#' @return
#' @export
#'
#' @examples
CompareAcrossSamples <- function(RomaData, Groups, GSThr = 1e-3, Mode = "Wil", Type = "Over",
                                 TestMode = "Aov+Tuk", TestPV1 = 5e-2, TestPV2 = 5e-2, PlotDiag = FALSE) {
  
  while(1){
    
    if(Mode == 'Wil' & Type == "Over"){
      print(paste("Using genestes overdispersed according to Wilcoxon test. GSThr =", GSThr))
      Selected <- RomaData$PVVectMat[,1]<Thr
      break
    }
    
    if(Mode == 'Wil' & Type == "Under"){
      print(paste("Using genestes underdispersed according to Wilcoxon test. GSThr =", GSThr))
      Selected <- RomaData$PVVectMat[,2]<Thr
      break
    }
    
    if(Mode == 'PPV' & Type == "Over"){
      print(paste("Using genestes underdispersed according to pseudo pv. GSThr =", GSThr))
      Selected <- RomaData$ModuleMatrix[,2]<Thr
      break
    }
    
    if(Mode == 'PPV' & Type == "Over"){
      print(paste("Using genestes underdispersed according to pseudo pv. GSThr =", GSThr))
      Selected <- RomaData$ModuleMatrix[,2]<1-Thr
      break
    }
    
    return(NULL)
  }
  
  MeltData <- reshape::melt(RomaData$PC1Matrix[Selected,])
  
  MeltData <- cbind(MeltData, Samples)
  
  colnames(MeltData) <- c("GeneSet", "Sample", "Value", "Group") 
  
  if(TestMode == "Aov+Tuk"){
    
    print("Performing Type III AOV (R default)")
    
    AOVFitTypeI <- aov(formula = Value ~ Group/GeneSet, data = MeltData)
    
    if(PlotDiag){
      plot(AOVFitTypeI)
    }
    
    SumData <- summary(AOVFitTypeI)
    print(SumData)
    
    print("Calculating Tukey Honest Significant Differences")
    
    TukTest <- TukeyHSD(AOVFitTypeI, conf.level = 1-TestPV2)
    
    if(PlotDiag){
      plot(TukTest)
    }
    
    if(SumData[[1]][1,5] < TestPV1){
      print("A significant difference is observed across groups")
      
      p <- ggplot2::ggplot(MeltData, ggplot2::aes(y=Value, x=Group, fill=Group)) +
        ggplot2::geom_boxplot() + ggplot2::guides(fill = "none") +
        ggplot2::labs(y="PC1 weigth", x="Groups", title = "Groups")
      
      print(p)
      
      
      Diffs <- TukTest$Group
      Diffs <- data.frame(Diffs[Diffs[,4] < TestPV2,])
      Diffs <- cbind(rownames(Diffs), Diffs)
      colnames(Diffs)[1] <- "GDiff"
      
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
    
    print("Performing Type III AOV (R default)")
    
    AOVFitTypeI <- aov(formula = Value ~ Group/GeneSet, data = MeltData)
    
    if(PlotDiag){
      plot(AOVFitTypeI)
    }
    
    SumData <- summary(AOVFitTypeI)
    print(SumData)
    
    print("Calculating Tukey Honest Significant Differences")
    
    TukTest <- TukeyHSD(AOVFitTypeI, conf.level = 1-TestPV2)
    
    if(PlotDiag){
      plot(TukTest, las=2)
    }
    
    if(SumData[[1]][2,5] < TestPV1){
      print("A significant difference is observed across groups and metagenes")
      
      Sep <- seq(from = 1, to = length(unique(MeltData$GeneSet)), by = 4)
      if(max(Sep) < length(unique(MeltData$GeneSet))+1) Sep <- c(Sep, length(unique(MeltData$GeneSet))+1)
      
      for(i in 2:length(Sep)){
        
        p <- ggplot2::ggplot(MeltData[as.integer(MeltData$GeneSet) %in% Sep[i-1]:(Sep[i]-1),],
                    ggplot2::aes(y=Value, x=Group, fill=Group)) + ggplot2::geom_boxplot() +
          ggplot2::labs(y="PC1 weight", x="Groups", title = paste("Geneset VS Groups - Part", i-1)) +
          ggplot2::facet_wrap( ~ GeneSet, ncol = 2) + ggplot2::theme(strip.text.x = ggplot2::element_text(size=6, face = "bold")) +
          ggplot2::guides(fill = "none")
        
        print(p)
      }
      
      Diffs <- TukTest$`Group:GeneSet`
      Diffs <- data.frame(Diffs[Diffs[,4] < TestPV2,])
      Diffs <- cbind(rownames(Diffs), Diffs)
      colnames(Diffs)[1] <- "GGDiff"
      
      SameGS <- unlist(
        lapply(
          lapply(
            lapply(
              strsplit(gsub(":", "-", as.character(Diffs$GGDiff)), c("-")),
              "[", c(2,4)),
            duplicated),
          any))
      
      
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
      
      
      p <- ggplot2::ggplot(Diffs[SameGS,], ggplot2::aes(x = GGDiff, y = diff, ymin = lwr, ymax = upr)) +
        ggplot2::geom_point() + ggplot2::geom_errorbar(width=0.25) + ggplot2::coord_flip() +
        ggplot2::facet_wrap( ~ GGDiff, scales = "free_y", ncol = 1) +
        ggplot2::labs(y="Mean difference", x="Categories", title = "Groups within Genesets") +
        ggplot2::theme(
          # axis.text.x = element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          strip.text.x = ggplot2::element_text(size=6, face = "bold"))
      
      print(p)
      
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

