#' Plot genesets information
#'
#' @param RomaData list, the analysis returned by rRoma
#' @param Selected vector, integer. The position of the genesets to plot
#' @param GenesetMargin scalat, integer. The number of rows used to draw the geneset names.
#' @param SampleMargin scalat, integer. The number of rows used to draw the sample names.
#' @param ColorGradient vector, string. The colors used for the heatmap.
#' @param cluster_cols boolean, should the samp^le be reordered according to the dendrogram?
#' @param GroupInfo vector, character. A vector describing the group association of each sample.
#' @param HMTite scalar, string. The title of the heatmap
#' @param AggByGroupsFL list, string. A list of function names (as strings) that will be used to aggregate the
#' geneset weigths and produce additional heatmaps.
#' @param Normalize boolean, shuold weights be normalized to c(-1, 1) for each geneset
#' @param Transpose boolean, should the samples by plotted on the rows instead of the columns?
#'
#' @return
#' @export
#'
#' @examples
Plot.Genesets <- function(RomaData, Selected = NULL,
                          GenesetMargin = 4, SampleMargin = 4,
                          ColorGradient = colorRamps::blue2red(50),
                          cluster_cols = FALSE, GroupInfo = NULL,
                          HMTite = "Selected Genesets", AggByGroupsFL = list(),
                          Normalize = FALSE, Transpose = FALSE){
  
  if(is.null(Selected)){
    Selected <- 1:nrow(RomaData$SampleMatrix)
  }
  
  if(length(intersect(Selected, 1:nrow(RomaData$SampleMatrix)))<1){
    print("No Genset selected")
    return(NULL)
  } else {
    print(paste(length(intersect(Selected, 1:nrow(RomaData$SampleMatrix))), "geneset selected"))
  }
  
  op <- par(mar=c(GenesetMargin, 5, 4, 2))
  
  B <- boxplot(lapply(RomaData$ModuleSummary[Selected], function(x){
    tSampleExp <- sapply(x$SampledExp, "[[", "ExpVar")
    if(!is.null(ncol(tSampleExp))){
      tSampleExp[1,]
    } else {
      tSampleExp
    }
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
  
  if(!is.null(GroupInfo)){
    
    FoundSamp <- intersect(colnames(RomaData$SampleMatrix), names(GroupInfo))
    print(paste(length(FoundSamp), "samples have an associated group"))
    
    GroupInfo <- GroupInfo[FoundSamp]
    
    if(length(FoundSamp) > 2){
      if(!is.factor(GroupInfo)){
        GroupInfo <- as.factor(GroupInfo)
      }
      AddInfo <- data.frame(Groups = GroupInfo)
    } else {
      AddInfo = NULL
    }
    
  } else {
    AddInfo = NULL
  }
  
  PlotMat <- RomaData$SampleMatrix[Selected,]
  
  if(Normalize){
    
    PlotMat <- t(apply(PlotMat, 1, function(x){
      x[x>0] <- x[x>0]/max(x[x>0])
      x[x<0] <- -x[x<0]/min(x[x<0])
      x
    }))
    
  }
  
  if(is.null(dim(PlotMat))){
    dim(PlotMat) <- c(1, length(PlotMat))
    colnames(PlotMat) <- names(RomaData$SampleMatrix[Selected,])
    rownames(PlotMat) <- rownames(RomaData$SampleMatrix)[Selected]
    tClusCol <- FALSE
    tClusRow <- FALSE
  } else {
    tClusCol <- cluster_cols
    tClusRow <- TRUE
  }
  
  if(Transpose){
    pheatmap::pheatmap(t(PlotMat), color = ColorGradient,
                       main = HMTite, cluster_rows = tClusCol,
                       cluster_cols = tClusRow,
                       annotation_row = AddInfo)
  } else {
    pheatmap::pheatmap(PlotMat, color = ColorGradient,
                       main = HMTite, cluster_cols = tClusCol,
                       cluster_rows = tClusRow,
                       annotation_col = AddInfo)
  }
  
  
  
  if(length(AggByGroupsFL)>0 & !is.null(AddInfo)){

    if(length(Selected) == 1){
      tDF <- data.frame(PlotMat[,FoundSamp])
      colnames(tDF) <- rownames(RomaData$SampleMatrix)[Selected]
      SplitData <- split(tDF, f=AddInfo$Groups)
    } else {
      SplitData <- split(data.frame(t(PlotMat[,FoundSamp])), f=AddInfo$Groups)
    }
    
    RetData <- list()
    
    for(i in 1:length(AggByGroupsFL)){
      
      Aggmat <- sapply(SplitData, function(x) {
        apply(x, 2, get(AggByGroupsFL[[i]]))
        })
      
      if(is.null(dim(Aggmat))){
        dim(Aggmat) <- c(1, length(Aggmat))
        colnames(Aggmat) <- names(SplitData)
        rownames(Aggmat) <- rownames(RomaData$SampleMatrix)[Selected]
        tClusCol <- FALSE
        tClusRow <- FALSE
      } else {
        tClusCol <- cluster_cols
        tClusRow <- TRUE
      }
      
      if(Transpose){
        pheatmap::pheatmap(t(Aggmat), color = ColorGradient,
                           main = paste(HMTite, "/", AggByGroupsFL[[i]]),
                           cluster_rows = tClusCol, cluster_cols = tClusRow)
      } else {
        pheatmap::pheatmap(Aggmat, color = ColorGradient,
                           main = paste(HMTite, "/", AggByGroupsFL[[i]]),
                           cluster_cols = tClusCol, cluster_rows = tClusRow)
      }
      
      RetData[[length(RetData)+1]] <- Aggmat

    }
    
    names(RetData) <- AggByGroupsFL

    return(RetData)
  }
  
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
    Selected <- 1:nrow(RomaData$SampleMatrix)
  }
  
  if(length(intersect(Selected, 1:nrow(RomaData$SampleMatrix)))<1){
    print("No Genset selected")
    return(NULL)
  } else {
    print(paste(length(intersect(Selected, 1:nrow(RomaData$SampleMatrix))), "geneset selected"))
  }
  
  
  for(i in Selected){
    
    FiltGenes <- RomaData$ModuleSummary[[i]]$GeneWeight * RomaData$ModuleSummary[[i]]$CorrectSign1
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
      ggplot2::labs(x = "Gene weight", y = "Gene") +
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

















#' Plot sample score across selected samples
#'
#' @param RomaData list, the analysis returned by rRoma
#' @param PlotSamples scalar, numeric. The number of samples to plot
#' @param ExpressionMatrix matrix, numeric. The expression matrix used to produce gene expression boxplot. If NULL (default), no gene expression information is reported
#' @param LogExpression boolean, should gene expression be logtransformed?
#' @param Selected vector, integer. The position of the genesets to plot
#'
#' @return
#' @export
#'
#' @examples
PlotSampleProjections <- function(RomaData, PlotSamples = 40,
                                ExpressionMatrix = NULL, LogExpression = TRUE,
                                Selected = NULL){
  
  if(is.null(Selected)){
    Selected <- 1:nrow(RomaData$SampleMatrix)
  }
  
  if(length(intersect(Selected, 1:nrow(RomaData$SampleMatrix)))<1){
    print("No Genset selected")
    return(NULL)
  } else {
    print(paste(length(intersect(Selected, 1:nrow(RomaData$SampleMatrix))), "geneset selected"))
  } 
  
  for(i in Selected){
    
    FiltSamp <- RomaData$ModuleSummary[[i]]$SampleScore * RomaData$ModuleSummary[[i]]$CorrectSign1
    # names(FiltSamp) <- colnames(RomaData$SampleMatrix)
    
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
      ggplot2::labs(x = "Sample Score", y = "Sample") +
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



  
#' Plot PC projections score across selected samples
#'
#' @param RomaData list, the analysis returned by rRoma
#' @param Selected vector, integer. The position of the genesets to plot
#' @param PlotPCProj vector, string. The plotting modality of projections. It can containing any combination of the following strings: 'Points', 'Density', or 'Bins'.
#' Any other value will result in projections not being plotted
#'
#' @return
#' @export
#'
#' @examples
PlotPCProjections <- function(RomaData, Selected = NULL, PlotPCProj = 'none'){
    
    if(is.null(Selected)){
      Selected <- 1:nrow(RomaData$SampleMatrix)
    }
    
    if(length(intersect(Selected, 1:nrow(RomaData$SampleMatrix)))<1){
      print("No Genset selected")
      return(NULL)
    } else {
      print(paste(length(intersect(Selected, 1:nrow(RomaData$SampleMatrix))), "geneset selected"))
    } 
    
    for(i in Selected){
      
      if(RomaData$InputPars$PCAType == "DimensionsAreSamples"){
        PrjList <- lapply(RomaData$ModuleSummary[[i]]$SampledExp, "[[", "GenesWei")
      }
      
      if(RomaData$InputPars$PCAType == "DimensionsAreGenes"){
        PrjList <- lapply(RomaData$ModuleSummary[[i]]$SampledExp, "[[", "SampleScore")
      }
      
      
      if(any(sapply(PrjList, is.null)) | !all(PlotPCProj != "none")){
        next()
      }
      
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



#' Plot the weigth of genes appearing across multiple genesets
#'
#' @param RomaData list, the analysis returned by rRoma
#' @param Selected vector, integer. The indices of the genesets to plot 
#' @param GenesByGroup scalar, integer. The number of genes per plot
#' @param MinMult scalar, integer. The minimal multiplicity to plot 
#'
#' @return
#' @export
#'
#' @examples
PlotRecurringGenes  <- function(RomaData, Selected = NULL,
                                GenesByGroup = 20, MinMult = 2){
  
  if(is.null(Selected)){
    Selected <- 1:nrow(RomaData$SampleMatrix)
  }
  
  if(length(intersect(Selected, 1:nrow(RomaData$SampleMatrix)))<1){
    print("No Genset selected")
    return(NULL)
  } else {
    print(paste(length(intersect(Selected, 1:nrow(RomaData$SampleMatrix))), "geneset selected"))
  }
  
  AllGenes <- lapply(RomaData$ModuleSummary, "[[", "UsedGenes")
  GeneMult <- table(unlist(AllGenes))
  GeneMult <- sort(GeneMult[GeneMult >= MinMult])

  if(length(GeneMult) == 0){
    return(NULL)
  }
  
  AllGenesWei <- lapply(RomaData$ModuleSummary, function(x){
    x$GeneWeight * x$CorrectSign1
  })
  
  GenesModuleMat <- sapply(AllGenesWei, function(x){x[names(GeneMult)]})

  rownames(GenesModuleMat) <- names(GeneMult)
  colnames(GenesModuleMat) <- unlist(lapply(RomaData$ModuleSummary, "[[", "ModuleName"))

  GenesByMult <- split(data.frame(GenesModuleMat), GeneMult)
  
  for(i in 1:length(GenesByMult)){
    ToPlotData <- t(data.matrix(GenesByMult[[i]]))
    # hist(apply(ToPlotData, 2, var, na.rm=TRUE),
    #      main = paste("Genes with multiplicity", names(GenesByMult)[i]))
    
    if(nrow(GenesByMult[[i]]) > 1){
      ToPlotData <- ToPlotData[,order(apply(ToPlotData, 2, var, na.rm=TRUE))]
    }
    
    
    
    par(mfcol=c(1,2))
    
    
    
    boxplot(apply(ToPlotData, 2, var, na.rm=TRUE), log='y', las=2,
            main = paste("Genes with multiplicity", names(GenesByMult)[i]),
            ylab = "Variance")
    
    boxplot(apply(ToPlotData, 2, var, na.rm=TRUE)/
              abs(apply(ToPlotData, 2, mean, na.rm=TRUE)), las=2, log='y',
            main = paste("Genes with multiplicity", names(GenesByMult)[i]),
            ylab = "Coefficient of variation")
    
    par(mfcol=c(1,1))
    
    Seps <- seq(from=0, to=ncol(ToPlotData), by=GenesByGroup)
    if(max(Seps) != ncol(ToPlotData)){
      Seps <- c(Seps, ncol(ToPlotData))
    }
    
    for(j in 2:length(Seps)){
      if(Seps[j-1]+1 != Seps[j]){
        BoxPlotData <- ToPlotData[,(Seps[j-1]+1):Seps[j]]
        YLim <- range(BoxPlotData, na.rm = TRUE)
        YLim <- YLim + c(-diff(YLim)*.05, diff(YLim)*.05)
        
        # CVData <- apply(ToPlotData[,(Seps[j-1]+1):Seps[j]], 2, var, na.rm=TRUE)/
        #   abs(apply(ToPlotData[,(Seps[j-1]+1):Seps[j]], 2, mean, na.rm=TRUE))
        # YLim2 <- range(CVData)
        
        boxplot(BoxPlotData, horizontal=FALSE, las = 2,
                main = paste("Genes with multiplicity", names(GenesByMult)[i]),
                ylab = "Weight (across genesets)", ylim = YLim)
        # points( ((CVData-min(CVData))/max(CVData-min(CVData)))*(YLim[2]-YLim[1])+YLim[1],
        #         pch = 18, col="blue", cex = 2)
        
      } else {
        boxplot(ToPlotData[,Seps[j]], horizontal=FALSE, las = 2,
                main = paste("Genes with multiplicity", names(GenesByMult)[i]),
                ylab = "Weight", xlab=colnames(ToPlotData)[Seps[j]])
        # barplot(var(ToPlotData[,Seps[j]], na.rm=TRUE)/
        #           abs(mean(ToPlotData[,Seps[j]], na.rm=TRUE)),
        #         las= 2, log = 'y', xlab=colnames(ToPlotData)[Seps[j]])
      }
    }
  }
}






#' Title
#'
#' @param RomaData 
#' @param GeneName 
#' @param ExpressionMatrix 
#' @param Selected 
#' @param GroupInfo 
#' @param DO_tSNE 
#' @param initial_dims 
#' @param perplexity 
#'
#' @return
#' @export
#'
#' @examples
#' 
ExploreGeneProperties <- function(
  RomaData,
  GeneName,
  ExpressionMatrix,
  Selected = NULL,
  GroupInfo = NULL,
  DO_tSNE = FALSE,
  initial_dims = 50,
  perplexity = 30
){
  
  
  VarVect <- apply(ExpressionMatrix, 1, var)
  
  if(!is.null(GroupInfo)){
    ExpMatByGroup <- split(data.frame(t(ExpressionMatrix)),
                           f=GroupInfo[colnames(ExpressionMatrix)])
    
    VarByGroup <- sapply(ExpMatByGroup, function(x) {
      apply(x, 2, var)
    })
    
    rownames(VarByGroup) <- names(VarVect)
    
    VarByGroup <- cbind(VarByGroup, Total = VarVect)
    
  } else {
    
    VarByGroup <- data.frame(X1 = names(VarVect), Total = VarVect)
    
  }
  
 
  VarByGroup.Melt <- reshape::melt(VarByGroup)
  colnames(VarByGroup.Melt) <- c("Gene", "Condition", "value")
  
  VarByGroup.Melt.Selected <- VarByGroup.Melt[VarByGroup.Melt$Gene == GeneName,]
  
  p <- ggplot2::ggplot(data = VarByGroup.Melt, mapping = ggplot2::aes(x=value, fill = Condition)) +
    ggplot2::geom_density() +
    ggplot2::geom_vline(data = VarByGroup.Melt.Selected,
                        mapping = ggplot2::aes(xintercept = value), color="black") +
    ggplot2::scale_x_log10() +
    ggplot2::facet_wrap(~Condition) +
    ggplot2::labs(x = "Variance of expression", title = GeneName) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  print(p)
  
  
  
  
  
  
  
  
  
  
  
  
  MeanVect <- apply(ExpressionMatrix, 1, mean)
  
  if(!is.null(GroupInfo)){
    
    MeanByGroup <- sapply(ExpMatByGroup, function(x) {
      apply(x, 2, mean)
    })
    
    rownames(MeanByGroup) <- names(MeanVect)
    
    MeanByGroup <- cbind(MeanByGroup, Total = MeanVect)
    
  } else {
    
    MeanByGroup <- data.frame(X1 = names(MeanVect), Total = MeanVect)
    
  }
  
  
  MeanByGroup.Melt <- reshape::melt(MeanByGroup)
  colnames(MeanByGroup.Melt) <- c("Gene", "Condition", "value")
  
  MeanByGroup.Melt.Selected <- MeanByGroup.Melt[MeanByGroup.Melt$Gene == GeneName,]
  
  p <- ggplot2::ggplot(data = MeanByGroup.Melt, mapping = ggplot2::aes(x=value, fill = Condition)) +
    ggplot2::geom_density() +
    ggplot2::geom_vline(data = MeanByGroup.Melt.Selected,
                        mapping = ggplot2::aes(xintercept = value), color="black") +
    ggplot2::facet_wrap(~Condition) +
    ggplot2::labs(x = "Mean expression", title = GeneName) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  print(p)
  
  
  
  
  
  
  
  
  
  
  ExpData <- data.frame(
    Sample = colnames(ExpressionMatrix),
    Exp = ExpressionMatrix[GeneName,],
    Condition = "Total",
    stringsAsFactors = FALSE
    )
  
  ToDisplay <- list()
  
  if(!is.null(GroupInfo)){
    ExpData <- rbind(
      ExpData,
      data.frame(
        Sample = colnames(ExpressionMatrix),
        Exp = ExpressionMatrix[GeneName,],
        Condition = GroupInfo[colnames(ExpressionMatrix)],
        stringsAsFactors = FALSE
      )
    )
    
    PWComp <- pairwise.wilcox.test(ExpData$Exp, ExpData$Condition)$p.value
    ExtPWComp <- matrix(rep(NA, length(unique(ExpData$Condition))^2), nrow = length(unique(ExpData$Condition)))
    colnames(ExtPWComp) <- sort(unique(ExpData$Condition))
    rownames(ExtPWComp) <- colnames(ExtPWComp)
    
    ExtPWComp[rownames(PWComp), colnames(PWComp)] <- PWComp
    
    ToDisplay <- apply(which(ExtPWComp <= .05, arr.ind = TRUE), 1, list)
    ToDisplay <- lapply(ToDisplay, function(x){as.integer(unlist(x))})
  }
  

  
  
  
  
  
  p <- ggplot2::ggplot(data = ExpData, mapping = ggplot2::aes(y = Exp, x = Condition, fill=Condition)) + 
    ggplot2::geom_boxplot() +
    ggplot2::labs(y = "Expression", title = GeneName) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggsignif::geom_signif(comparisons = ToDisplay,
                          map_signif_level=TRUE, test = "wilcox.test", step_increase = .1)
  
  print(p)
  
  
  
  
  
  
  
  
  
  PCA_TotData <- irlba::prcomp_irlba(x = t(ExpressionMatrix), n = 2, center = TRUE, retx = TRUE)
  
  p <- ggplot2::ggplot(data = data.frame(PCA_TotData$x,
                                    Exp = ExpressionMatrix[GeneName,]),
                  mapping = ggplot2::aes(x=PC1, y=PC2, color=Exp)) +
    ggplot2::geom_point() +
    ggplot2::labs(title = GeneName) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  print(p)
  
  
  
  
  
  
  if(DO_tSNE){
    Data <- Rtsne::Rtsne(X = t(ExpressionMatrix), initial_dims = initial_dims, perplexity = perplexity)
    
    p <- ggplot2::ggplot(data = data.frame(Dim1 = Data$Y[,1], Dim2 = Data$Y[,2],
                                           Exp = ExpressionMatrix[GeneName,]),
                         mapping = ggplot2::aes(x=Dim1, y=Dim2, color=Exp)) +
      ggplot2::geom_point() +
      ggplot2::labs(title = GeneName, y = "tSNE dimension 2", x = "tSNE dimension 1") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    
    print(p)
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  if(is.null(Selected)){
    Selected <- 1:nrow(RomaData$SampleMatrix)
  }
  
  ModulesGenes <- lapply(RomaData$ModuleSummary[Selected], "[[", "UsedGenes")
  names(ModulesGenes) <- sapply(RomaData$ModuleSummary[Selected], "[[", "ModuleName")
  
  WhichModules <- sapply(ModulesGenes, function(x){
    GeneName %in% x
  })
  
  Found <- sapply(RomaData$ModuleSummary, "[[", "ModuleName") %in%
    names(which(WhichModules))
  
  print(paste("Gene found in", sum(Found), "modules"))
  
  WeiList <- lapply(RomaData$ModuleSummary[Found], function(x){
    x$GeneWeight * x$CorrectSign1
  })
  
  names(WeiList) <- sapply(RomaData$ModuleSummary[Found], "[[", "ModuleName")
  
  WeiList <- lapply(1:length(WeiList), function(i) {
    x <- WeiList[[i]]
    return(
      data.frame(Wei = c(x, abs(x)),
                 Name = rep(names(x), 2),
                 Order = c(rank(x, ties.method = "average"),
                           rank(abs(x), ties.method = "average"))/length(x),
                 Type = rep(c("Signed", "Absolute"), each = length(x)),
                 Module = rep(names(WeiList)[i], 2*length(x)),
                 stringsAsFactors = FALSE)
    )
  })
  
  
  Wei.DF <- do.call(rbind, WeiList)
  Wei.DF.Selected <- Wei.DF[Wei.DF$Name == GeneName,]
  
  
  p <- ggplot2::ggplot(data = Wei.DF,
                       mapping = ggplot2::aes(y=Wei, x="")) +
    ggplot2::geom_boxplot() +
    ggplot2::facet_grid(Module ~ Type, scales = "free") + 
    ggplot2::geom_point(data = Wei.DF.Selected,
                        mapping = ggplot2::aes(y=Wei, x=""),
                        inherit.aes = FALSE,
                        color = "red",
                        size = 3) +
    ggplot2::labs(y = "Weight", x = "", title = GeneName) + 
    ggplot2::geom_text(data = Wei.DF.Selected,
                       mapping = ggplot2::aes(x = "", y = Wei,
                                              label = paste(signif(100*Order, 2), "th")
                       ),
                       hjust = 0, nudge_x = 0.05, size = 7,
                       inherit.aes = FALSE) + 
    ggplot2::geom_text(data = Wei.DF.Selected,
                       mapping = ggplot2::aes(x = "", y = Wei,
                                              label = signif(Wei, 4)
                       ),
                       hjust = 1, nudge_x = -0.05, size = 7,
                       inherit.aes = FALSE) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  print(p)
  
  
  
  GeneExp <- ExpressionMatrix[GeneName,]
  
  ModScoreList <- lapply(RomaData$ModuleSummary[Found], function(x){
    x$SampleScore * x$CorrectSign1
  })
  names(ModScoreList) <- sapply(RomaData$ModuleSummary[Found], "[[", "ModuleName")
  
  ModMat <- do.call(cbind, ModScoreList)
  Mod.DF <- data.frame(ModMat,
                       Sample = rownames(ModMat),
                       stringsAsFactors = FALSE)
  
  Mod.DF <- reshape::melt(Mod.DF, id.vars = "Sample")
  
  if(!is.null(GroupInfo)){
    
    Mod.DF <- data.frame(Mod.DF, Exp = GeneExp[Mod.DF$Sample], Groups = GroupInfo[Mod.DF$Sample])
    
    p <- ggplot2::ggplot(data = Mod.DF, mapping = ggplot2::aes(x = value, y = Exp)) +
      ggplot2::geom_smooth(color = "black") +
      ggplot2::geom_point(mapping = ggplot2::aes(color = Groups)) + 
      ggplot2::facet_wrap(~variable) +
      ggplot2::labs(x = "Sample Score", y = "Gene expression", title = GeneName) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    
   
    
  } else {
    
    Mod.DF <- data.frame(Mod.DF, Exp = GeneExp[Mod.DF$Sample])
    
    p <- ggplot2::ggplot(data = Mod.DF, mapping = ggplot2::aes(x = value, y = Exp)) +
      ggplot2::geom_smooth() + 
      ggplot2::geom_point() + 
      ggplot2::facet_wrap(~variable) +
      ggplot2::labs(x = "Sample Score", y = "Gene expression", title = GeneName) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    
  }
  
  print(p)
  
}
                      









