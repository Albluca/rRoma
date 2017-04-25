#' Correct the sign of the principal component 
#'
#' @param Wei vector, numeric optional vector of weigths
#' @param Mode scalar, character. Mode to correct the sign
#' @param Thr scalar, numeric quantile threshold 
#' @param PCWeigth 
#' @param PCProj 
#' @param DefWei 
#' @param Grouping 
#' @param ExpMat 
#' @param CorMethod c("pearson", "kendall", "spearman")
#'
#' @return
#' @export
#'
#' @examples
FixPCSign <-
  function(PCWeigth, PCProj, Wei = NULL, Mode = 'none', DefWei = 1,
           Thr = NULL, Grouping = NULL, ExpMat = NULL, CorMethod = "pearson") {
    
    if (Mode == 'none') {
      
      print("Orienting PC using a random direction")
      
      return(1)
    }
    
    if (Mode == 'PreferActivation') {
      
      print("Orienting PC by preferential activation")
      
      ToUse <- rep(TRUE, length(PC1Weigth))
      if (!is.null(Thr)) {
        ToUse <- abs(PCWeigth) >= quantile(abs(PCWeigth), Thr)
      }
      
      if (sum(PC1Weigth[ToUse]) < 0) {
        return(-1)
      } else {
        return(+1)
      }
    }
    
    if (Mode == 'UseAllWeights') {
      
      print(paste("Missing gene weights will be replaced by", DefWei))
      
      Wei[is.na(Wei)] <- DefWei
      Mode <-  'UseKnownWeights'
      
    }
    
    if (Mode == 'UseKnownWeights') {
      
      print("Orienting PC by combining PC weights and gene weights")
      print(paste("Missing gene weights will be replaced by", 0))
      
      ToUse <- rep(TRUE, length(PCWeigth))
      if (!is.null(Thr)) {
        ToUse <- abs(PCWeigth) >= quantile(abs(PCWeigth), Thr)
      }
      
      if (sum(Wei[!is.na(Wei) & ToUse] * PCWeigth[!is.na(Wei) & ToUse]) < 0) {
        return(-1)
      } else {
        return(+1)
      }
    }
    
    if (Mode == 'CorrelateAllWeightsByGene') {
      
      print(paste("Missing gene weights will be replaced by", DefWei))
      Wei[is.na(Wei)] <- DefWei
      Mode <-  'CorrelateKnownWeightsByGene'
      
    }
    
    if (Mode == 'CorrelateKnownWeightsByGene') {
      
      print(paste("Orienting PC by correlating gene expression and PC projections (", CorMethod, ")", sep = ''))
      
      if(sum(!is.na(Wei))<1){
        print("Not enough weights, PC will be oriented randomly")
        return(1)
      }
      
      
      if(!is.null(Grouping)){
        print("Using groups")
        
        AssocitedGroups <- Grouping[colnames(ExpMat)]
        TB <- table(AssocitedGroups, useNA = "no")
        
        if(sum(TB>0)<2){
          print("Not enough groups, PC will be oriented randomly")
          return(1)
        }
        
        GroupMedians <- apply(ExpMat[!is.na(Wei),], 1, function(x) {
          aggregate(x = x, by=list(AssocitedGroups), FUN = median)
        })
        
        MedianProj <- aggregate(PCProj, by=list(AssocitedGroups), FUN = median)
        
        print("Computing correlations")
        Cor.Test.Vect <- sapply(GroupMedians, function(x){
          if(length(unique(x[,2])) > 2 & length(unique(MedianProj[,2])) > 2){
            CT <- cor.test(x[,2], MedianProj[,2], method = CorMethod)
            c(CT$p.value, CT$estimate)
          } else {
            c(NA, NA)
          }
        })
        
        SelGenesWei <- Wei[!is.na(Wei)]
        names(SelGenesWei) <- names(GroupMedians)
        
        ToUse <- !is.na(Cor.Test.Vect[2,])
        if (!is.null(Thr)) {
          ToUse <- (Cor.Test.Vect[1,] < Thr) & ToUse
        }
        
        if(sum(Cor.Test.Vect[2,ToUse]*SelGenesWei[ToUse])>0){
          return(1)
        } else {
          return(-1)
        }
        
      } else {
        print("Not using groups")
        
        names(PCProj) <- colnames(ExpMat)
        
        print("Computing correlations")
        Cor.Test.Vect <- apply(ExpMat, 1, function(x){
          if(length(unique(x)) > 2 & length(unique(PCProj)) > 2){
            CT <- cor.test(x, PCProj, method = CorMethod)
            c(CT$p.value, CT$estimate)
          } else {
            c(NA, NA)
          }
        })
        
        print("Correcting using weights")
        if(sum(!is.na(Wei))>1){
          Cor.Test.Vect[2,!is.na(Wei)] <- Cor.Test.Vect[2,!is.na(Wei)]*Wei[!is.na(Wei)]
        }
        
        ToUse <- !is.na(Cor.Test.Vect[2,])
        if (!is.null(Thr)) {
          ToUse <- (Cor.Test.Vect[1,] < Thr) & ToUse
        }
        
        if(sum(Cor.Test.Vect[2,ToUse])>0){
          return(1)
        } else {
          return(-1)
        }
        
      }
      
    }
    
    if (Mode == 'CorrelateAllWeightsBySample') {
      
      print(paste("Missing gene weights will be replaced by", DefWei))
      Wei[is.na(Wei)] <- DefWei
      Mode <-  'CorrelateKnownWeightsBySample'
      
    }
    
    if (Mode == 'CorrelateKnownWeightsBySample') {
      
      print(paste("Orienting PC by correlating gene expression and PC weights (", CorMethod, ")", sep = ''))
      
      if(sum(!is.na(Wei))<1){
        print("Not enough weights, PC will be oriented randomly")
        return(1)
      }
      
      if(!is.null(Grouping)){
        print("Using groups")
        
        AssocitedGroups <- Grouping[colnames(ExpMat)]
        TB <- table(AssocitedGroups, useNA = "no")
        
        if(sum(TB>0)<2){
          print("Not enough groups, PC will be oriented randomly")
          return(1)
        }
        
        GroupMedians <- apply(ExpMat, 1, function(x) {
          aggregate(x = x, by=list(AssocitedGroups), FUN = median, na.rm=TRUE)
        })
        
        MediansByGroups <- sapply(GroupMedians, function(x) {
          x[,2]
        })
        
        print("Computing correlations")
        Cor.Test.Vect <- apply(MediansByGroups, 1, function(x){
          if(length(unique(x)) > 2 & length(unique(PCWeigth*Wei)) > 2){
            CT <- cor.test(x, PCWeigth*Wei, method = CorMethod)
            c(CT$p.value, CT$estimate)
          } else {
            c(NA, NA)
          }
        })
        
        ToUse <- !is.na(Cor.Test.Vect[2,])
        if (!is.null(Thr)) {
          ToUse <- (Cor.Test.Vect[1,] < Thr) & ToUse
        }
        
        if(sum(Cor.Test.Vect[2,ToUse])>0){
          return(1)
        } else {
          return(-1)
        }
        
      } else {
        print("Not using groups")
        
        names(PCWeigth) <- rownames(ExpMat)
        PCWeigth <- PCWeigth*Wei
        
        print("Computing correlations")
        Cor.Test.Vect <- apply(ExpMat, 2, function(x){
          if(length(unique(x)) > 2 & length(unique(PCWeigth)) > 2){
            CT <- cor.test(x, PCWeigth, method = CorMethod)
            c(CT$p.value, CT$estimate)
          } else {
            c(NA, NA)
          }
        })
        
        
        ToUse <- !is.na(Cor.Test.Vect[2,])
        if (!is.null(Thr)) {
          ToUse <- (Cor.Test.Vect[1,] < Thr) & ToUse
        }
        
        if(sum(Cor.Test.Vect[2,ToUse])>0){
          return(1)
        } else {
          return(-1)
        }
        
      }
      
    }
    
  }








