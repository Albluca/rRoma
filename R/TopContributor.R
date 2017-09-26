
#' Get the top contributing genes per geneset
#'
#' @param RomaData list, the analysis returned by rRoma
#' @param Selected vector, integer. The indices of the genesets to plot 
#' @param nGenes either a value between 0 and 1 or a positive integer larger than 1. If the parameter
#' is an integer larger than 1, it will be interpreted as the number of genes to return per geneset.
#' If the parameter is a value between 1 and 0 it will be interpreted as the ratio or genes to return per
#' geneset (e.g., .1 indicates 10\% of the genes of the geneset)).
#' @param Mode string, the mode used to determine the top contributing genes. It can be "Wei"
#' (gene weigths will be used) or "Cor" (pearson correlation between gene expression and module score will be used)
#' @param Plot boolean, shoul the summary heatmap be plotted?
#' @param ExpressionMatrix numerical matrix, the expression matrix in the same format used by the rRoma.R function
#' @param OrderType string scalar. The mode of selection for the top contributing genes.
#' The current implementation allow "Abs" (genes with the largest weight in absolute value),
#' "Pos" (genes with the largest weight), and "Neg" (genes with the largest negative weight)
#' @param ColorGradient vector, string. The colors used for the heatmap.
#' @param cluster_cols boolean, should the genesets be reordered according to the dendrogram?
#' @param HMTite scalar, string. The title of the heatmap.
#' @param Transpose boolean, should the genes by plotted on the rows instead of the columns?
#' @param CorMethod string, the method to use to compute the correlation ("pearson", "kendall",
#' or "spearman"). Methods other than "pearson" will produce potenital problems in the derivation
#' of statistics accociated with the correlation
#'
#' @return
#' @export
#'
#' @examples
GetTopContrib <- function(RomaData, Selected = NULL, nGenes = .1,
                          Mode = "Wei", Plot = FALSE, ExpressionMatrix = NULL,
                          ColorGradient = colorRamps::blue2red(50),
                          cluster_cols = FALSE, HMTite = "Top contributing genes",
                          OrderType = "Abs", Transpose = FALSE, CorMethod = "pearson") {
  
  if(is.null(Selected)){
    Selected <- 1:nrow(RomaData$SampleMatrix)
  }
  
  if(is.null(ExpressionMatrix) & Mode == "Cor"){
    print("Impossible to compute cerrelations if expression matrix is not provided. Using weigths.")
    Mode <- "Wei"
  }
  
  if(Mode == "Cor"){
    
    AllData <- lapply(Selected, function(i){
      Scores <- RomaData$SampleMatrix[i,]
      Genes <- RomaData$ModuleSummary[[i]]$UsedGenes
      AllTests <- apply(ExpressionMatrix[Genes, names(Scores)], 1, cor.test, y = Scores, method = CorMethod)
      PVs <- sapply(AllTests, "[[", "p.value")
      Est <- sapply(AllTests, "[[", "estimate")
      CI <- sapply(AllTests, "[[", "conf.int")
      
      if(nGenes < 1){
        nGenes = round(length(Genes)*nGenes)
      }
      
      if(nGenes > length(Genes)){
        nGenes = length(Genes)
      }
      
      if(OrderType == "Abs"){
        Sel <- sort(abs(Est), decreasing = TRUE, index.return=TRUE)$ix[1:nGenes]
      }
      
      if(OrderType == "Pos"){
        Sel <- sort(Est, decreasing = TRUE, index.return=TRUE)$ix[1:nGenes]
        Sel <- Sel[Sel %in% which(Est > 0)]
      }
      
      if(OrderType == "Neg"){
        Sel <- sort(Est, decreasing = FALSE, index.return=TRUE)$ix[1:nGenes]
        Sel <- Sel[Sel %in% which(Est < 0)]
      }
      
      DF <- data.frame(
        Gene = names(PVs)[Sel], Module= RomaData$ModuleSummary[[i]]$ModuleName,
        p.value = PVs[Sel], estimate = Est[Sel], CI.low = CI[1,Sel], CI.high = CI[2,Sel]
      )
      rownames(DF) <- NULL
      
      return(DF)
      
    })
    
    RetTable <- do.call(rbind, AllData)
    
    VizMat <- matrix(0, nrow = length(unique(RetTable$Gene)), ncol = length(unique(RetTable$Module)))
    rownames(VizMat) <- sort(unique(as.character(RetTable$Gene)))
    colnames(VizMat) <- sort(unique(as.character(RetTable$Module)))
    
    tt <- apply(RetTable, 1, function(x) {
      VizMat[x["Gene"], x["Module"]] <<- as.numeric(x["estimate"])
      # VizMat["estimate"]
    })
    
    HMTite <- paste(HMTite, "(Cor)")
    
  }
  
  
  
  
  
  
  if(Mode == "Wei"){
    
    if(OrderType == "Abs"){
      GenesWei <- 
        lapply(
          lapply(RomaData$WeigthList[Selected], function(x){
            x[order(abs(x), decreasing = TRUE)]
          }),
          function(x){
            
            if(nGenes<1){
              nGenes <- round(length(x)*nGenes)
            }
            
            if(nGenes > length(x)){
              nGenes <- length(x)
            }
            
            x[1:nGenes]
          })

    }
    
    if(OrderType == "Pos"){
      GenesWei <- lapply(
        lapply(RomaData$WeigthList[Selected], sort, decreasing = TRUE), 
        function(x){
          
          if(nGenes<1){
            nGenes <- round(length(x)*nGenes)
          }
          
          if(nGenes > length(x)){
            nGenes <- length(x)
          }
          
          Ret <- x[1:nGenes]
          Ret[ (RomaData$WeigthList[Selected])[names(Ret)]>0 ]
          
        })
    }
    
    if(OrderType == "Neg"){
      GenesWei <- lapply(
        lapply(RomaData$WeigthList[Selected], sort, decreasing = FALSE), 
        function(x){
          
          if(nGenes<1){
            nGenes <- round(length(x)*nGenes)
          }
          
          if(nGenes > length(x)){
            nGenes <- length(x)
          }
          
          Ret <- x[1:nGenes]
          Ret[ (RomaData$WeigthList[Selected])[names(Ret)]<0 ]
          
        })
    }
    
    ModuleNameList <- lapply(RomaData$ModuleSummary, "[[", "ModuleName")
    GeneWeiList <- lapply(GenesWei, length)
    
    RetTable <- data.frame(
      Gene = unlist(lapply(GenesWei, names)),
      Module = unlist(
        lapply(1:length(ModuleNameList), function(i){rep(ModuleNameList[[i]], GeneWeiList[[i]])})
      ),
      Weigth = unlist(GenesWei)
    )
    
    VizMat <- matrix(0, nrow = length(unique(RetTable$Gene)), ncol = length(unique(RetTable$Module)))
    rownames(VizMat) <- sort(unique(as.character(RetTable$Gene)))
    colnames(VizMat) <- sort(unique(as.character(RetTable$Module)))
    
    tt <- apply(RetTable, 1, function(x) {
      VizMat[x["Gene"], x["Module"]] <<- as.numeric(x["Weigth"])
      # VizMat["estimate"]
    })
    
    HMTite <- paste(HMTite, "(Wei)")
    
  }
  
  if(Plot){
    if(Transpose){
      pheatmap::pheatmap(t(VizMat), color = ColorGradient, main = HMTite, cluster_rows = cluster_cols)
    } else {
      pheatmap::pheatmap(VizMat, color = ColorGradient, main = HMTite, cluster_cols = cluster_cols)
    }
  }
  
  return(list(Table = RetTable, VizMat = VizMat))
}
