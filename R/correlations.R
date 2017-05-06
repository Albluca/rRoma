GetCorrelations <- function(RomaData, Selected = NULL, MatData, Methods = "pearson", ConfLevel = .95) {
  
  if(is.null()){
    Selected <- 1:length(RomaData$ModuleSummary)
  }
  
  AllCor <- NULL
  
  SampleIds <- intersect(colnames(MatData), colnames(RomaData$ProjMatrix))
  
  for(i in Selected){
    
    CorVal <- apply(MatData[RomaData$ModuleSummary[[i]]$UsedGenes, SampleIds], 1, function(x){
      CT <- cor.test(x = x, y = RomaData$ProjMatrix[i, SampleIds], conf.level = ConfLevel, method = Methods)
      c(CT$estimate, CT$conf.int, CT$p.value)
    })
    
    CorVal <- rbind(colnames(CorVal), CorVal, rep(rownames(RomaData$ProjMatrix)[i], ncol(CorVal)))
    rownames(CorVal) <- c("gene", "cor", "ci.low", "ci.high", "p.val", "GS")
    
    AllCor <- cbind(AllCor, CorVal)
    
  }
  
  rownames(AllCor) <- NULL
  
  AllCorGenes = t(AllCor)
  
  
  AllCor <- NULL
  
  for(i in Selected){
    
    CorVal <- apply(MatData[RomaData$ModuleSummary[[i]]$UsedGenes, SampleIds], 2, function(x){
      CT <- cor.test(x = x, y = (RomaData$WeigthList[[i]])[names(x)], conf.level = ConfLevel, method = Methods)
      c(CT$estimate, CT$conf.int, CT$p.value)
    })
    
    CorVal <- rbind(colnames(CorVal), CorVal, rep(rownames(RomaData$ProjMatrix)[i], ncol(CorVal)))
    rownames(CorVal) <- c("sample", "cor", "ci.low", "ci.high", "p.val", "GS")
    
    AllCor <- cbind(AllCor, CorVal)
    
  }
  
  AllCorSamples = t(AllCor)
  
  return(Genes = AllCorGenes, Samples = AllCorSamples)
  
}