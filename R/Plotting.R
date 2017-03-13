#' Plot overdispersed genesets
#'
#' @param RomaData 
#' @param Thr 
#' @param Mode 
#' @param GenesetMargin 
#' @param SampleMargin 
#'
#' @return
#' @export
#'
#' @examples
Plot.Overdispersed <- function(RomaData, Thr = 1, Mode = "Wil",
                               GenesetMargin = 4, SampleMargin = 4,
                               ColorGradient = colorRamps::blue2red(50),
                               cluster_cols = FALSE){
  
  while(1){
    
    if(Mode == 'Wil'){
      Selected <- RomaData$PVVectMat[,1]<Thr
      break
    }
    
    if(Mode == 'PPV'){
      Selected <- RomaData$ModuleMatrix[,2]<Thr
      break
    }
    
    return(NULL)
  }
  
  op <- par(mar=c(GenesetMargin, 5, 4, 2))
  
  B <- boxplot(lapply(lapply(RomaData$ModuleSummary[Selected], "[[", "SampledExp"), function(x) x[1,]),
               at = 1:sum(Selected), las = 2, ylab = "Explained variance", main = "Overdispersed genesets",
               names = unlist(lapply(RomaData$ModuleSummary[Selected], "[[", "ModuleName")),
               ylim = c(0,1))
  points(x = 1:sum(Selected), y = RomaData$ModuleMatrix[Selected,1], pch = 20, col="red", cex = 2)
  
  par(op)
  
  
  # PC1ProjList <- RomaData$ProjLists
  # 
  # PC1ProjList.Scaled <- lapply(PC1ProjList, function(x){
  #   x[x>0] <- log10(x[x>0])
  #   x[x<0] <- -log10(-x[x<0])
  #   return(x)
  # }) 
  
  # for(i in 1:sum(Selected)){
  #   
  #   SortIds <- sort(PC1ProjList.Scaled[[i]], index.return=TRUE)
  #   
  #   B <- boxplot(list("Activators"=SortIds$x[SortIds$x>0], "InHIbitors"=-SortIds$x[SortIds$x<0]),
  #                main = RomaData$ModuleSummary[[which(Selected)[i]]]$ModuleName, ylab = "Controbution (Log10)")
  #   
  #   # B <- barplot(SortIds$x[SortIds$x>0], main = paste(RomaData$ModuleSummary[[which(Selected)[i]]]$ModuleName, "(Activators)"), horiz = TRUE)
  #   # text(y = B, x = 0, labels = (RomaData$ModuleSummary[[which(Selected)[i]]]$UsedGenes)[SortIds$ix[SortIds$x>0]],
  #   #      srt=0, cex=.5, col='red')
  #   # 
  #   # B <- barplot(-SortIds$x[SortIds$x<0], main = paste(RomaData$ModuleSummary[[which(Selected)[i]]]$ModuleName, "(Inhibitors)"), horiz = TRUE)
  #   # text(y = B, x = 0, labels = (RomaData$ModuleSummary[[which(Selected)[i]]]$UsedGenes)[SortIds$ix[SortIds$x<0]],
  #   #      srt=0, cex=.5, col='red')
  # }
  
  pheatmap::pheatmap(RomaData$ProjMatrix[Selected,], color = ColorGradient,
                     main = "Overdispersed genesets", cluster_cols = cluster_cols)

}














#' Plot underdispersed
#'
#' @param RomaData 
#' @param Thr 
#' @param Mode 
#' @param GenesetMargin 
#' @param SampleMargin 
#'
#' @return
#' @export
#'
#' @examples
Plot.Underdispersed <- function(RomaData, Thr = 1, Mode = "Wil",
                                GenesetMargin = 4, SampleMargin = 4,
                                ColorGradient = colorRamps::blue2red(50),
                                # ExpressionMatrix = NULL,
                                cluster_cols = FALSE){
  
  while(1){
    
    if(Mode == 'Wil'){
      Selected <- RomaData$PVVectMat[,2]<Thr
      break
    }
    
    if(Mode == 'PPV'){
      Selected <- RomaData$ModuleMatrix[,2]<1-Thr
      break
    }
    
    return(NULL)
  }
  
  op <- par(mar=c(GenesetMargin, 5, 4, 2))
  
  B <- boxplot(lapply(lapply(RomaData$ModuleSummary[Selected], "[[", "SampledExp"), function(x) x[1,]),
               at = 1:sum(Selected), las = 2, ylab = "Explained variance", main = "Underdispersed genesets",
               names = unlist(lapply(RomaData$ModuleSummary[Selected], "[[", "ModuleName")),
               ylim = c(0,1))
  points(x = 1:sum(Selected), y = RomaData$ModuleMatrix[Selected,1], pch = 20, col="red", cex = 2)
  
  par(op)
  
  # if(!is.null(ExpressionMatrix)){
  #   
  #   MedianByGS <- sapply(lapply(RomaData$ModuleSummary[Selected], "[[", "UsedGenes"), function(gl){
  #     median(ExpressionMatrix[gl,])
  #   })
  #   
  #   MedianByGSAll <- sapply(lapply(RomaData$ModuleSummary, "[[", "UsedGenes"), function(gl){
  #     median(ExpressionMatrix[gl,])
  #   })
  #   
  #   
  #   boxplot(MedianByGSAll)
  #   points(MedianByGS, col='red', pch=20)
  #   
  #   
  # }
  # 
  # PC1ProjList <- RomaData$ProjLists
  # 
  # PC1ProjList.Scaled <- lapply(PC1ProjList, function(x){
  #   x[x>0] <- log10(x[x>0])
  #   x[x<0] <- -log10(-x[x<0])
  #   return(x)
  # }) 
  # 
  # for(i in 1:sum(Selected)){
  #   SortIds <- sort(PC1ProjList.Scaled[[i]], index.return=TRUE)
  #   
  #   B <- boxplot(list("Activators"=SortIds$x[SortIds$x>0], "Inhibitors"=-SortIds$x[SortIds$x<0]),
  #                main = RomaData$ModuleSummary[[which(Selected)[i]]]$ModuleName, ylab = "Controbution (Log10)")
  #   
  #   
  #   # B <- barplot(SortIds$x[SortIds$x>0], main = paste(RomaData$ModuleSummary[[which(Selected)[i]]]$ModuleName, "(Activators)"), horiz = TRUE)
  #   # text(y = B, x = 0, labels = (RomaData$ModuleSummary[[which(Selected)[i]]]$UsedGenes)[SortIds$ix[SortIds$x>0]],
  #   #      srt=0, cex=.5, col='red')
  #   # 
  #   # B <- barplot(-SortIds$x[SortIds$x<0], main = paste(RomaData$ModuleSummary[[which(Selected)[i]]]$ModuleName, "(Inhibitors)"), horiz = TRUE)
  #   # text(y = B, x = 0, labels = (RomaData$ModuleSummary[[which(Selected)[i]]]$UsedGenes)[SortIds$ix[SortIds$x<0]],
  #   #      srt=0, cex=.5, col='red')
  # }
  
  pheatmap::pheatmap(RomaData$ProjMatrix[Selected,], color = ColorGradient,
                     main = "Underdispersed genesets", cluster_cols = cluster_cols)
  
}





#' Title
#'
#' @param RomaData 
#'
#' @return
#' @export
#'
#' @examples
PlotGeneWeight <- function(RomaData, PlotGenes = 40,
                                ExpressionMatrix = NULL, LogExpression = TRUE,
                                ModuleToPlot = NULL){
  
  if(is.null(ModuleToPlot)){
    ModuleToPlot <- 1:length(RomaData$ModuleSummary)
  } 
  
  for(i in ModuleToPlot){
    
    FiltGenes <- RomaData$ModuleSummary[[i]]$PC1Weight.SignFixed
    names(FiltGenes) <- RomaData$ModuleSummary[[i]]$UsedGenes
    
    GeneNames <- RomaData$ModuleSummary[[i]]$UsedGenes
    
    nGenes <- min(PlotGenes, length(GeneNames))
    
    DF <- data.frame(cbind(names(FiltGenes),
                           FiltGenes),
                     row.names = NULL)
    colnames(DF) <- c("Gene", "Value")
    
    DF$Value <- as.numeric(as.character(DF$Value))
    DF <- DF[order(abs(DF$Value), decreasing = TRUE)[1:nGenes],]
    DF <- DF[order(DF$Value), ]
    
    DF$Gene <- factor(as.character(DF$Gene), levels = DF$Gene)
    
    p <- ggplot2::ggplot(data = DF, ggplot2::aes(y = Gene, x = Value)) + ggplot2::geom_point() +
      # ggplot2::scale_y_continuous(breaks = DF$Position[1:nGenes], labels = DF$Gene[1:nGenes]) +
      # ggplot2::scale_y_discrete(breaks=levels(DF$Gene)) +
      ggplot2::labs(x = "PC1 weight", y = "Gene") +
      ggplot2::theme(panel.grid.minor = ggplot2::element_line(colour = NA))
    
    if(is.null(ExpressionMatrix)){
      
      print(p + ggplot2::ggtitle(paste(RomaData$ModuleSummary[[i]]$ModuleName,
                                       "L1=", signif(RomaData$ModuleMatrix[i,1], 3),
                                       "L1/L2=", signif(RomaData$ModuleMatrix[i,3], 3))))
      
    } else {
      
      ReshapedData <- reshape::melt(ExpressionMatrix[as.character(DF$Gene), ])
      
      ReshapedData$X1 <- factor(as.character(ReshapedData$X1), levels = levels(DF$Gene))
      
      p1 <- ggplot2::ggplot(ReshapedData, ggplot2::aes(x=X1, y=value)) + 
        ggplot2::geom_boxplot() + ggplot2::coord_flip() +
        ggplot2::labs(x = "Gene", y = "Expression")
      
      if(LogExpression){
        p1 <- p1 + ggplot2::scale_y_log10()
      }
      
      gridExtra::grid.arrange(p, p1, ncol=2, top=paste(RomaData$ModuleSummary[[i]]$ModuleName,
                                                               "L1=", signif(RomaData$ModuleMatrix[i,1], 3),
                                                               "L1/L2=", signif(RomaData$ModuleMatrix[i,3], 3)))
      
      # boxplot(t(ExpressionMatrix[as.character(DF$Gene[1:nGenes])[order(DF$Position[1:nGenes])], ]))
      
    }
    
  }
  
}

















#' Title
#'
#' @param RomaData 
#'
#' @return
#' @export
#'
#' @examples
PlotSampleProjections <- function(RomaData, PlotSamples = 40,
                                ExpressionMatrix = NULL, LogExpression = TRUE,
                                ModuleToPlot = NULL){
  
  if(is.null(ModuleToPlot)){
    ModuleToPlot <- 1:length(RomaData$ModuleSummary)
  } 
  
  for(i in ModuleToPlot){
    
    FiltSamp <- RomaData$ModuleSummary[[i]]$PCABase$x[,"PC1"]
    names(FiltSamp) <- colnames(RomaData$ProjMatrix)
    
    GeneNames <- RomaData$ModuleSummary[[i]]$UsedGenes
    
    nSamples <- min(PlotSamples, length(FiltSamp))
    
    DF <- data.frame(cbind(names(FiltSamp),
                           FiltSamp),
                     row.names = NULL)
    colnames(DF) <- c("Samples", "Value")
    
    DF$Value <- as.numeric(as.character(DF$Value))
    DF <- DF[order(abs(DF$Value), decreasing = TRUE)[1:nSamples],]
    DF <- DF[order(DF$Value), ]
    
    DF$Samples <- factor(as.character(DF$Samples), levels = DF$Samples)
    
    p <- ggplot2::ggplot(data = DF, ggplot2::aes(y = Samples, x = Value)) + ggplot2::geom_point() +
      # ggplot2::scale_y_continuous(breaks = DF$Position[1:nGenes], labels = DF$Gene[1:nGenes]) +
      # ggplot2::scale_y_discrete(breaks=levels(DF$Gene)) +
      ggplot2::labs(x = "PC1 projection", y = "Sample") +
      ggplot2::theme(panel.grid.minor = ggplot2::element_line(colour = NA))
    
    if(is.null(ExpressionMatrix)){
      
      print(p + ggplot2::ggtitle(paste(RomaData$ModuleSummary[[i]]$ModuleName,
                                       "L1=", signif(RomaData$ModuleMatrix[i,1], 3),
                                       "L1/L2=", signif(RomaData$ModuleMatrix[i,3], 3))))
      
    } else {
      
      ReshapedData <- reshape::melt(ExpressionMatrix[GeneNames, as.character(DF$Samples)])
      
      ReshapedData$X2 <- factor(as.character(ReshapedData$X2), levels = levels(DF$Samples))
      
      p1 <- ggplot2::ggplot(ReshapedData, ggplot2::aes(x=X2, y=value)) + 
        ggplot2::geom_boxplot() + ggplot2::coord_flip() +
        ggplot2::labs(x = "Sample", y = "Expression")
      
      if(LogExpression){
        p1 <- p1 + ggplot2::scale_y_log10()
      }
      
      gridExtra::grid.arrange(p, p1, ncol=2, top=paste(RomaData$ModuleSummary[[i]]$ModuleName,
                                                               "L1=", signif(RomaData$ModuleMatrix[i,1], 3),
                                                               "L1/L2=", signif(RomaData$ModuleMatrix[i,3], 3)))
      
      # boxplot(t(ExpressionMatrix[as.character(DF$Gene[1:nGenes])[order(DF$Position[1:nGenes])], ]))
      
    }
    
  }
  
}





