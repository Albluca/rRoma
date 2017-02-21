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
Plot.Overdispersed <- function(RomaData, Thr = 1, Mode = "Wil", GenesetMargin = 4, SampleMargin = 4){
  
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
  
  PC1ProjList <- RomaData$ProjLists
  
  PC1ProjList.Scaled <- lapply(PC1ProjList, function(x){
    x[x>0] <- log10(x[x>0])
    x[x<0] <- -log10(-x[x<0])
    return(x)
  }) 
  
  for(i in 1:sum(Selected)){
    
    SortIds <- sort(PC1ProjList.Scaled[[i]], index.return=TRUE)
    
    B <- boxplot(list("Activators"=SortIds$x[SortIds$x>0], "InHIbitors"=-SortIds$x[SortIds$x<0]),
                 main = RomaData$ModuleSummary[[which(Selected)[i]]]$ModuleName, ylab = "Controbution (Log10)")
    
    # B <- barplot(SortIds$x[SortIds$x>0], main = paste(RomaData$ModuleSummary[[which(Selected)[i]]]$ModuleName, "(Activators)"), horiz = TRUE)
    # text(y = B, x = 0, labels = (RomaData$ModuleSummary[[which(Selected)[i]]]$UsedGenes)[SortIds$ix[SortIds$x>0]],
    #      srt=0, cex=.5, col='red')
    # 
    # B <- barplot(-SortIds$x[SortIds$x<0], main = paste(RomaData$ModuleSummary[[which(Selected)[i]]]$ModuleName, "(Inhibitors)"), horiz = TRUE)
    # text(y = B, x = 0, labels = (RomaData$ModuleSummary[[which(Selected)[i]]]$UsedGenes)[SortIds$ix[SortIds$x<0]],
    #      srt=0, cex=.5, col='red')
  }
  
  gplots::heatmap.2(RomaData$PC1Matrix[Selected,],
                    margins = c(SampleMargin, GenesetMargin), trace = 'none',
                    col = colorRamps::blue2red)
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
Plot.Underdispersed <- function(RomaData, Thr = 1, Mode = "Wil", GenesetMargin = 4, SampleMargin = 4){
  
  while(1){
    
    if(Mode == 'Wil'){
      Selected <- RomaData$PVVectMat[,2]<Thr
      break
    }
    
    if(Mode == 'PPV'){
      Selected <- RomaData$ModuleMatrix[,2]>Thr
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
  
  PC1ProjList <- RomaData$ProjLists
  
  PC1ProjList.Scaled <- lapply(PC1ProjList, function(x){
    x[x>0] <- log10(x[x>0])
    x[x<0] <- -log10(-x[x<0])
    return(x)
  }) 
  
  for(i in 1:sum(Selected)){
    SortIds <- sort(PC1ProjList.Scaled[[i]], index.return=TRUE)
    
    B <- boxplot(list("Activators"=SortIds$x[SortIds$x>0], "InHIbitors"=-SortIds$x[SortIds$x<0]),
                 main = RomaData$ModuleSummary[[which(Selected)[i]]]$ModuleName, ylab = "Controbution (Log10)")
    
    
    # B <- barplot(SortIds$x[SortIds$x>0], main = paste(RomaData$ModuleSummary[[which(Selected)[i]]]$ModuleName, "(Activators)"), horiz = TRUE)
    # text(y = B, x = 0, labels = (RomaData$ModuleSummary[[which(Selected)[i]]]$UsedGenes)[SortIds$ix[SortIds$x>0]],
    #      srt=0, cex=.5, col='red')
    # 
    # B <- barplot(-SortIds$x[SortIds$x<0], main = paste(RomaData$ModuleSummary[[which(Selected)[i]]]$ModuleName, "(Inhibitors)"), horiz = TRUE)
    # text(y = B, x = 0, labels = (RomaData$ModuleSummary[[which(Selected)[i]]]$UsedGenes)[SortIds$ix[SortIds$x<0]],
    #      srt=0, cex=.5, col='red')
  }
  
  gplots::heatmap.2(RomaData$PC1Matrix[Selected,],
                    margins = c(SampleMargin, GenesetMargin), trace = 'none',
                    col = colorRamps::blue2red)
}





#' Title
#'
#' @param RomaData 
#'
#' @return
#' @export
#'
#' @examples
PlotGeneProjections <- function(RomaData, PlotGenes = 40){
  
  for(i in 1:length(RomaData$ModuleSummary)){
    
    FiltGenes <- RomaData$ModuleSummary[[i]]$PC1Projections.SignFixed
    names(FiltGenes) <- RomaData$ModuleSummary[[i]]$UsedGenes
    
    GeneNames <- RomaData$ModuleSummary[[i]]$UsedGenes
    
    nGenes <- min(PlotGenes, length(GeneNames))
    
    PC1FilStr <- rep(NA, length(GeneNames))
    names(PC1FilStr) <- GeneNames
    PC1FilStr[names(FiltGenes)] <- FiltGenes
    
    DF <- data.frame(cbind(1:length(PC1FilStr), GeneNames[order(PC1FilStr)],
                           PC1FilStr[order(PC1FilStr)],
                           rep("Original", length(PC1FilStr))
                           ),
                     row.names = NULL)
    
    DF <- DF[order(abs(PC1FilStr))[1:nGenes],]
    
    colnames(DF) <- c("Position", "Gene", "Value", "Filter")
    
    DF$Position <- as.numeric(as.character(DF$Position))
    DF$Value <- as.numeric(as.character(DF$Value))
    
    DF$Position <- rank(DF$Position, ties.method = "min")
    
    p <- ggplot2::ggplot(data = DF, ggplot2::aes(y = Position, x = Value)) + ggplot2::geom_point() +
      ggplot2::scale_y_continuous(breaks = DF$Position[1:nGenes], labels = DF$Gene[1:nGenes]) +
      ggplot2::ggtitle(RomaData$ModuleSummary[[i]]$ModuleName) +
      ggplot2::labs(x = "PC1 projection", y = "Gene") +
      ggplot2::theme(panel.grid.minor = ggplot2::element_line(colour = NA))
    
    print(p)
    
  }
  
}




