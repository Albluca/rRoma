#' Load a GMT file from disc
#'
#' @param FileLocation 
#'
#' @return
#' @export
#'
#' @examples
ReadGMTFile <- function(FileLocation) {
  
  Done <- FALSE
  GeneList <- list()
  Idx <- 0
  
  while(!Done){
    
    Idx <- Idx + 1
    
     tLine <- unlist(readr::read_delim(FileLocation, "\t", escape_double = FALSE,
                                         col_names = FALSE, trim_ws = TRUE, n_max = 1, skip = Idx - 1),
                              use.names = FALSE)
    
    if(length(tLine) < 3){
      Done <- TRUE
    } else {
      
      tGenes <- tLine[-c(1:2)]
      tGenes <- tGenes[!is.na(tGenes)]
      tGenes <- strsplit(tGenes, split = "]", fixed = TRUE)
      tGenes <- lapply(tGenes, strsplit, split = "[", fixed = TRUE)
      tGenes <- lapply(tGenes, "[[", 1)
      
      GeneNames <- unlist(lapply(tGenes, "[[", 1))
      
      Weigths <- rep(NA, length(GeneNames))

      if(any(lapply(tGenes, length) > 1)){
        Weigths[lapply(tGenes, length) > 1] <- unlist(lapply(tGenes[lapply(tGenes, length) > 1], "[[", 2))
      }
      
      Weigths <- as.numeric(Weigths)

      GeneList[[Idx]] <- list(Name = tLine[1], Desc = tLine[2], Genes = GeneNames, Weigths = Weigths)
    }
    
  }
  
  return(GeneList)
}





# Load data
# Msigdb.v6.0.symbols.all <- rRoma::ReadGMTFile("~/Downloads/msigdb.v6.0.symbols.gmt.txt")
# devtools::use_data(Msigdb.v6.0.symbols.all)



#' Return a subset of Msigdb containing the provided search string
#'
#' @param SearchString string scalar, the keyword to search 
#'
#' @return
#' @export
#'
#' @examples
SelectFromMSIGdb <- function(SearchString, Version = "6.0") {
  
  if(Version == "6.0"){
    print("Searching in MsigDB v6.0")
    
    return(rRoma::Msigdb.v6.0.symbols.all[grep(pattern = tolower(SearchString),
                                               x = tolower(unlist(lapply(rRoma::Msigdb.v6.0.symbols.all, "[[", "Name"))),
                                               fixed = TRUE)])
  }
  
  if(Version == "5.2"){
    print("Searching in MsigDB v5.2")
    
    return(rRoma::Msigdb.v5.2.symbols.all[grep(pattern = tolower(SearchString),
                                               x = tolower(unlist(lapply(rRoma::Msigdb.v5.2.symbols.all, "[[", "Name"))),
                                               fixed = TRUE)])
  }

  print("Versione unavailable")
}






