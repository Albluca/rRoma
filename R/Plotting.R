#' Plot genesets information
#'
#' @param RomaData list, the analysis returned by rRoma
#' @param Selected vector, integer. The position of the genesets to plot
#' @param GenesetMargin scalat, integer. The number of rows used to draw the geneset names.
#' @param SampleMargin scalat, integer. The number of rows used to draw the sample names.
#' @param ColorGradient vector, string. The colors used for the heatmap.
#' @param cluster_cols boolean, should the samp^le be reordered according to the dendrogram?
#'
#' @return
#' @export
#'
#' @examples
Plot.Genesets <- function(RomaData, Selected = NULL,
                               GenesetMargin = 4, SampleMargin = 4,
                               ColorGradient = colorRamps::blue2red(50),
                               cluster_cols = FALSE){
  
  if(is.null(Selected)){
    Selected <- 1:nrow(RomaData$ProjMatrix)
  }
  
  if(length(intersect(Selected, 1:nrow(RomaData$ProjMatrix)))<1){
    print("No Genset selected")
    return(NULL)
  } else {
    print(paste(length(intersect(Selected, 1:nrow(RomaData$ProjMatrix))), "geneset selected"))
  }
  
  op <- par(mar=c(GenesetMargin, 5, 4, 2))
  
  B <- boxplot(lapply(RomaData$ModuleSummary[Selected], function(x){
    sapply(x$SampledExp, "[[", "ExpVar")[1,]
  }), at = 1:length(Selected), las = 2, ylab = "Explained variance", main = "Selected genesets",
  names = unlist(lapply(RomaData$ModuleSummary[Selected], "[[", "ModuleName")),
  ylim = c(0,1))
  points(x = 1:length(Selected), y = RomaData$ModuleMatrix[Selected,1], pch = 20, col="red", cex = 2)
  
  PlotData <- lapply(RomaData$ModuleSummary[Selected], function(x){
    sapply(x$SampledExp, "[[", "MedianExp")
    })
  
  B <- boxplot(PlotData, at = 1:length(Selected), las = 2, ylab = "Median expression", main = "Selected genesets",
               names = unlist(lapply(RomaData$ModuleSummary[Selected], "[[", "ModuleName")),
               ylim = range(c(unlist(PlotData), RomaData$ModuleMatrix[Selected,5]), na.rm = TRUE))
  points(x = 1:length(Selected), y = RomaData$ModuleMatrix[Selected,5], pch = 20, col="red", cex = 2)
  
  par(op)
  
  pheatmap::pheatmap(RomaData$ProjMatrix[Selected,], color = ColorGradient,
                     main = "Selected genesets", cluster_cols = cluster_cols)

}










#' Plot gene weigth across selected samples
#'
#' @param RomaData list, the analysis returned by rRoma
#' @param Selected vector, integer. The position of the genesets to plot
#' @param PlotGenes scalar, numeric. The number of genes to plot
#' @param ExpressionMatrix matrix, numeric. The expression matrix used to produce gene expression boxplot. If NULL (default), no gene expression information is reported
#' @param LogExpression boolean, should gene expression be logtransformed?
#' @param PlotWeigthSign boolean, should the sign of the genes weigth be used to color the plots?
#'
#' @return
#' @export
#'
#' @examples
PlotGeneWeight <- function(RomaData, PlotGenes = 40,
                           ExpressionMatrix = NULL, LogExpression = TRUE,
                           Selected = NULL, PlotWeigthSign = FALSE){
  
  if(is.null(Selected)){
    Selected <- 1:nrow(RomaData$ProjMatrix)
  }
  
  if(length(intersect(Selected, 1:nrow(RomaData$ProjMatrix)))<1){
    print("No Genset selected")
    return(NULL)
  } else {
    print(paste(length(intersect(Selected, 1:nrow(RomaData$ProjMatrix))), "geneset selected"))
  }
  
  
  for(i in Selected){
    
    FiltGenes <- RomaData$ModuleSummary[[i]]$PC1Weight.SignFixed
    names(FiltGenes) <- RomaData$ModuleSummary[[i]]$UsedGenes
    
    GeneNames <- RomaData$ModuleSummary[[i]]$UsedGenes
    
    nGenes <- min(PlotGenes, length(GeneNames))
    
    GMTWei = sign(RomaData$ModuleSummary[[i]]$GMTWei)
    names(GMTWei) = names(FiltGenes)
    GMTWei[is.na(GMTWei)] <- 0
    
    DF <- data.frame(cbind(names(FiltGenes), FiltGenes,
                           GMTWei),
                     row.names = NULL)
    colnames(DF) <- c("Gene", "Value", "Wei")
    
    DF$Value <- as.numeric(as.character(DF$Value))
    DF <- DF[order(abs(DF$Value), decreasing = TRUE)[1:nGenes],]
    DF <- DF[order(DF$Value), ]
    
    DF$Gene <- factor(as.character(DF$Gene), levels = DF$Gene)
    DF$Wei <- factor(as.character(DF$Wei), levels = c("-1", "0", "1"))
    
    if(PlotWeigthSign){
      if(any(as.character(DF$Wei)=="0")){
        p <- ggplot2::ggplot(data = DF, ggplot2::aes(y = Gene, x = Value, color = Wei)) +
          ggplot2::scale_color_manual(values = c("red", "black", "blue"), guide=FALSE)
      } else {
        p <- ggplot2::ggplot(data = DF, ggplot2::aes(y = Gene, x = Value, color = Wei)) +
          ggplot2::scale_color_manual(values = c("red", "blue"), guide=FALSE)
      }
      
      
    } else {
      p <- ggplot2::ggplot(data = DF, ggplot2::aes(y = Gene, x = Value))
    }
    
     p <- p + ggplot2::geom_point() +
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
      ReshapedData <- cbind(ReshapedData, GMTWei[as.character(ReshapedData$X1)])
      colnames(ReshapedData)[4] <- "Wei"
      ReshapedData$Wei <- factor(as.character(ReshapedData$Wei), levels = c("-1", "0", "1"))
      
      if(PlotWeigthSign){
        if(any(as.character(ReshapedData$Wei)=="0")){
          p1 <- ggplot2::ggplot(ReshapedData, ggplot2::aes(x=X1, y=value, fill=Wei)) +
            ggplot2::scale_fill_manual(values = c("red", "white", "blue"), guide=FALSE)
        } else {
          p1 <- ggplot2::ggplot(ReshapedData, ggplot2::aes(x=X1, y=value, fill=Wei)) +
            ggplot2::scale_fill_manual(values = c("red", "blue"), guide=FALSE)
        }
      } else {
        p1 <- ggplot2::ggplot(ReshapedData, ggplot2::aes(x=X1, y=value))
      }
      
      p1 <- p1 + ggplot2::geom_boxplot() + ggplot2::coord_flip() + 
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

















#' Plot gene weigth across selected samples
#'
#' @param RomaData list, the analysis returned by rRoma
#' @param PlotSamples scalar, numeric. The number of samples to plot
#' @param ExpressionMatrix matrix, numeric. The expression matrix used to produce gene expression boxplot. If NULL (default), no gene expression information is reported
#' @param LogExpression boolean, should gene expression be logtransformed?
#' @param Selected vector, integer. The position of the genesets to plot
#' @param PlotPCProj vector, string. The plotting modality of projections. It can containing any combination of the following strings: 'Points', 'Density', or 'Bins'.
#' Any other value will result in projections not being plotted
#'
#' @return
#' @export
#'
#' @examples
PlotSampleProjections <- function(RomaData, PlotSamples = 40,
                                ExpressionMatrix = NULL, LogExpression = TRUE,
                                Selected = NULL, PlotPCProj = 'none'){
  
  if(is.null(Selected)){
    Selected <- 1:nrow(RomaData$ProjMatrix)
  }
  
  if(length(intersect(Selected, 1:nrow(RomaData$ProjMatrix)))<1){
    print("No Genset selected")
    return(NULL)
  } else {
    print(paste(length(intersect(Selected, 1:nrow(RomaData$ProjMatrix))), "geneset selected"))
  } 
  
  for(i in Selected){
    
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
    
    
    PrjList <- lapply(RomaData$ModuleSummary[[i]]$SampledExp, "[[", "PCProj")
    PC1Sample <- unlist(lapply(PrjList, function(x){if(all(!is.na(x))) x[,1]}), use.names = FALSE)
    PC2Sample <- unlist(lapply(PrjList, function(x){if(all(!is.na(x))) x[,2]}), use.names = FALSE)
    
    PC1Data <- RomaData$ModuleSummary[[i]]$PCABase$x[,1]*RomaData$ModuleSummary[[i]]$CorrectSign1
    PC2Data <- RomaData$ModuleSummary[[i]]$PCABase$x[,2]*RomaData$ModuleSummary[[i]]$CorrectSign2
    
    DF <- data.frame(PC1 = c(PC1Sample, PC1Data), PC2 = c(PC2Sample, PC2Data),
                     Source = c(rep("Sampling", length(PC1Sample)),
                                    rep("Data", length(PC1Data))
                                )
                     )
    
    XLims <- quantile(PC1Sample, c(.01, .99))
    XLims[1] <- min(XLims[1], min(PC1Data))
    XLims[2] <- max(XLims[2], max(PC1Data))
    
    YLims <- quantile(PC2Sample, c(.01, .99))
    YLims[1] <- min(YLims[1], min(PC2Data))
    YLims[2] <- max(YLims[2], max(PC2Data))
    
    if(any(PlotPCProj == 'Points')){
      p <- ggplot2::ggplot(DF, ggplot2::aes(x = PC1, y = PC2, colour = Source, alpha=Source)) + ggplot2::geom_point() +
        ggplot2::scale_alpha_manual(values=c(Data=1, Sampling=.2), guide=FALSE) +
        ggplot2::scale_x_continuous(limits = XLims) +
        ggplot2::scale_y_continuous(limits = YLims) +
        ggplot2::ggtitle(paste(RomaData$ModuleSummary[[i]]$ModuleName,
                               "L1=", signif(RomaData$ModuleMatrix[i,1], 3),
                               "L1/L2=", signif(RomaData$ModuleMatrix[i,3], 3))) +
        ggplot2::geom_hline(yintercept=0) + ggplot2::geom_vline(xintercept=0)
      
      print(p)
    }
    
    
    if(any(PlotPCProj == 'Density')){
      p <- ggplot2::ggplot(DF[DF$Source=="Sampling",], ggplot2::aes(x = PC1, y = PC2)) +
        ggplot2::stat_density_2d(ggplot2::aes(fill = ..density..), geom="raster", contour = FALSE, n = 250) +
        ggplot2::scale_x_continuous(limits = XLims) +
        ggplot2::scale_y_continuous(limits = YLims) +
        ggplot2::ggtitle(paste(RomaData$ModuleSummary[[i]]$ModuleName,
                               "L1=", signif(RomaData$ModuleMatrix[i,1], 3),
                               "L1/L2=", signif(RomaData$ModuleMatrix[i,3], 3))) +
        ggplot2::geom_hline(yintercept=0) + ggplot2::geom_vline(xintercept=0) +
        ggplot2::geom_point(data = DF[DF$Source=="Data",], mapping = ggplot2::aes(x = PC1, y = PC2), color="red")
      
      print(p)
    }
    
    
    
    if(any(PlotPCProj == 'Bins')){
      p <- ggplot2::ggplot(DF[DF$Source=="Sampling",], ggplot2::aes(x = PC1, y = PC2)) +
        ggplot2::geom_bin2d(bins=75) +
        ggplot2::scale_x_continuous(limits = XLims) +
        ggplot2::scale_y_continuous(limits = YLims) +
        ggplot2::ggtitle(paste(RomaData$ModuleSummary[[i]]$ModuleName,
                               "L1=", signif(RomaData$ModuleMatrix[i,1], 3),
                               "L1/L2=", signif(RomaData$ModuleMatrix[i,3], 3))) +
        ggplot2::geom_hline(yintercept=0) + ggplot2::geom_vline(xintercept=0) +
        ggplot2::geom_point(data = DF[DF$Source=="Data",], mapping = ggplot2::aes(x = PC1, y = PC2), color="red")
      
      print(p)
    }
    
  }
  
}





