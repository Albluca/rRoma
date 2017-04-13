#' Load a GMT file from a file
#'
#' @param FileLocation the location of the GMT file to load. It can be a local or remote location.
#' @param SearchString string vector, the keywords to search in the loaded GMT file
#' @param Mode string scalar, the search mode for the keywords.
#' It can be "ALL" (genesets with all the keywords will be returned) or
#' "ANY" (genesets with at at least one keyword will be returned)
#'
#' @return
#' @export
#'
#' @examples
ReadGMTFile <- function(FileLocation, SearchString = NULL, Mode = "ANY") {
  
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
  
  if(!is.null(SearchString)){
    SelGeneSets <- NULL
    for(i in 1:length(SearchString)){
      if(Mode == "ANY"){
        SelGeneSets <- union(SelGeneSets,
                             grep(pattern = tolower(SearchString[i]),
                                  x = tolower(unlist(lapply(GeneList, "[[", "Name"))),
                                  fixed = TRUE)
        )
      }
      if(Mode == "ALL"){
        if(!is.null(SelGeneSets)){
          SelGeneSets <- intersect(SelGeneSets,
                                   grep(pattern = tolower(SearchString[i]),
                                        x = tolower(unlist(lapply(GeneList, "[[", "Name"))),
                                        fixed = TRUE)
          )
        } else {
          SelGeneSets <- grep(pattern = tolower(SearchString[i]),
                              x = tolower(unlist(lapply(GeneList, "[[", "Name"))),
                              fixed = TRUE)
        }
      }
    }
    
    SelGeneSets <- unique(SelGeneSets)
  } else {
    SelGeneSets <- 1:length(GeneList)
  }
  
  
  print("The following genesets have been loaded and selected:")
  print(paste(unlist(lapply(GeneList[SelGeneSets], "[[", "Name")),
              " (",
              unlist(lapply(lapply(GeneList[SelGeneSets], "[[", "Genes"), length)),
              " genes)", sep= ''))
  
  return(GeneList[SelGeneSets])
  
}





# Load data
# Msigdb.v6.0.symbols.all <- rRoma::ReadGMTFile("~/Downloads/msigdb.v6.0.symbols.gmt.txt")
# devtools::use_data(Msigdb.v6.0.symbols.all)



#' Return a subset of Msigdb containing the provided search string
#'
#' @param SearchString string vector, the keywords to search 
#' @param Version string scalar, the version of Msigdb to use. Only "5.2" and "6.0" are currently available.
#' @param Mode string scalar, the search mode for the keywords.
#' It can be "ALL" (genesets with all the keywords will be returned) or
#' "ANY" (genesets with at at least one keyword will be returned)
#'
#' @return
#' @export
#'
#' @examples
SelectFromMSIGdb <- function(SearchString, Version = "6.0", Mode = "ANY") {
  
  InternalDB <- list()
  
  if(Version == "6.0"){
    
    print("Searching in MsigDB v6.0")
    InternalDB <- rRoma::Msigdb.v6.0.symbols.all
    
  }
  
  if(Version == "5.2"){
    
    print("Searching in MsigDB v5.2")
    InternalDB <- rRoma::Msigdb.v5.2.symbols.all
    
  }

  SelGeneSets <- NULL
  for(i in 1:length(SearchString)){
    if(Mode == "ANY"){
      SelGeneSets <- union(SelGeneSets,
                           grep(pattern = tolower(SearchString[i]),
                                x = tolower(unlist(lapply(InternalDB, "[[", "Name"))),
                                fixed = TRUE)
      )
    }
    if(Mode == "ALL"){
      if(!is.null(SelGeneSets)){
        SelGeneSets <- intersect(SelGeneSets,
                                 grep(pattern = tolower(SearchString[i]),
                                      x = tolower(unlist(lapply(InternalDB, "[[", "Name"))),
                                      fixed = TRUE)
        )
      } else {
        SelGeneSets <- grep(pattern = tolower(SearchString[i]),
                            x = tolower(unlist(lapply(InternalDB, "[[", "Name"))),
                            fixed = TRUE)
      }
    }
  }
  
  SelGeneSets <- unique(SelGeneSets)
  
  print("The following genesets have been selected:")
  print(paste(unlist(lapply(InternalDB[SelGeneSets], "[[", "Name")),
        " (",
        unlist(lapply(lapply(InternalDB[SelGeneSets], "[[", "Genes"), length)),
        " genes)",sep=''))
  
  return(InternalDB[SelGeneSets])
  
}






