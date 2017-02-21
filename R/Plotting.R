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
PlotGeneProjections <- function(RomaData){
  
  for(i in 1:length(RomaData$ModuleSummary)){
    names(RomaData$ModuleSummary[[i]]$PC1Projections.SignFixed) <- RomaData$ModuleSummary[[i]]$UsedGenes
    SortedPC1 <- sort(RomaData$ModuleSummary[[i]]$PC1Projections.SignFixed)
    nGenes <- length(RomaData$ModuleSummary[[i]]$UsedGenes)
    plot(1:nGenes, SortedPC1, type = 'n', xlab = '', xaxt = 'n', ylab = "PC1 projection",
         main = RomaData$ModuleSummary[[i]]$ModuleName)
    abline(h=0)
    Cols = rep('red', nGenes)
    Cols[SortedPC1>0] <- "blue"
    text(x = 1:nGenes, y = SortedPC1,
         labels = names(SortedPC1), srt= 90, col = Cols)
  }
  
}




