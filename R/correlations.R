#' Title
#'
#' @param RomaData
#' @param Selected
#' @param MatData
#' @param Methods
#' @param ConfLevel
#'
#' @return
#' @export
#'
#' @examples
GetCorrelations <- function(RomaData, Selected = NULL, MatData, Methods = "pearson",
                            ConfLevel = .95, Construct.Net = FALSE, Construct.Net.PVThr = 1e-2) {

  if(is.null(Selected)){
    Selected <- 1:length(RomaData$ModuleSummary)
  }

  GeneNets <- NULL

  AllCor <- NULL

  SampleIds <- intersect(colnames(MatData), colnames(RomaData$SampleMatrix))

  for(i in Selected){

    CorVal <- apply(MatData[RomaData$ModuleSummary[[i]]$UsedGenes, SampleIds], 1, function(x){
      CT <- cor.test(x = x, y = RomaData$SampleMatrix[i, SampleIds], conf.level = ConfLevel, method = Methods)
      c(CT$estimate, CT$conf.int, CT$p.value)
    })

    CorVal <- rbind(colnames(CorVal), CorVal, rep(rownames(RomaData$SampleMatrix)[i], ncol(CorVal)))

    if(Methods == "pearson"){
      rownames(CorVal) <- c("gene", "cor", "ci.low", "ci.high", "p.val", "GS")
    } else {
      rownames(CorVal) <- c("gene", "cor", "p.val", "GS")
    }

    AllCor <- cbind(AllCor, CorVal)

  }

  AllCorGenes = t(AllCor)


  AllCor <- NULL

  for(i in Selected){

    CorVal <- apply(MatData[RomaData$ModuleSummary[[i]]$UsedGenes, SampleIds], 2, function(x){
      CT <- cor.test(x = x, y = (RomaData$WeigthList[[i]])[names(x)], conf.level = ConfLevel, method = Methods)
      c(CT$estimate, CT$conf.int, CT$p.value)
    })

    CorVal <- rbind(colnames(CorVal), CorVal, rep(rownames(RomaData$SampleMatrix)[i], ncol(CorVal)))

    if(Methods == "pearson"){
      rownames(CorVal) <- c("sample", "cor", "ci.low", "ci.high", "p.val", "GS")
    } else {
      rownames(CorVal) <- c("sample", "cor", "p.val", "GS")
    }



    AllCor <- cbind(AllCor, CorVal)

  }

  AllCorSamples = t(AllCor)

  return(list(Genes = AllCorGenes, Samples = AllCorSamples))

}
