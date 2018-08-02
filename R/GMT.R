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
                                       col_names = FALSE, trim_ws = TRUE, n_max = 1, skip = Idx - 1,
                                       col_types = readr::cols(.default = readr::col_character())
                                       ),
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
SelectFromMSIGdb <- function(SearchString, Version = "6.2", Mode = "ANY") {
  
  return(SelectFromInternalDB(SearchString, BDName = "MsigDB", Version = Version, Mode = Mode))
  
}



















#' Return a an internal DB using the provided search string
#'
#' @param SearchString string vector, the keywords to search 
#' @param Version string scalar, the version of the internal DB to use. If NULL, the most recent version availble in the package will be used 
#' @param Mode string scalar, the search mode for the keywords.
#' It can be "ALL" (genesets with all the keywords will be returned) or
#' "ANY" (genesets with at at least one keyword will be returned)
#' @param BDName string, scalar. The internal DB to use, can be one of the following
#' \itemize{
#' \item 'MsigDB': Molecular signature database (http://software.broadinstitute.org/gsea/msigdb)
#' \item 'ACSN_Global': ACSN global map (https://acsn.curie.fr)
#' \item 'ACSN_Apoptosis': ACSN Apoptosis map (https://acsn.curie.fr)
#' \item 'ACSN_CellCycle': ACSN Cell Cycle map (https://acsn.curie.fr)
#' \item 'ACSN_DNARepair': ACSN DNA repair map (https://acsn.curie.fr)
#' \item 'ACSN_EMT': ACSN EMT map (https://acsn.curie.fr)
#' \item 'ACSN_Survival': ACSN Survival map (https://acsn.curie.fr)
#' }
#' 
#' @return
#' @export
#'
#'
#' @examples
SelectFromInternalDB <- function(SearchString, BDName = "MsigDB", Version = NULL, Mode = "ANY") {
  
  InternalDB <- list()
  
  if(BDName == "MsigDB"){
    
    if(is.null(Version)){
      
      print("Searching in MsigDB v6.2")
      InternalDB <- rRoma::Msigdb.v6.2.symbols.all
      
    } else {
      
      if(Version == "5.2"){
        print("Searching in MsigDB v5.2")
        InternalDB <- rRoma::Msigdb.v5.2.symbols.all
      }
      
      if(Version == "6.0"){
        print("Searching in MsigDB v6.0")
        InternalDB <- rRoma::Msigdb.v6.0.symbols.all
      }
      
      if(Version == "6.2"){
        print("Searching in MsigDB v6.2")
        InternalDB <- rRoma::Msigdb.v6.2.symbols.all
      }
      
      if(Version == "6.2E"){
        print("Searching in MsigDB v6.2 (Entrez)")
        InternalDB <- rRoma::Msigdb.v6.2.entrez.all
      }
      
    }
    
  }
  
  if(BDName == "ACSN_Global"){
    if(is.null(Version)){
      print("Searching in ACSN Global map v1.1")
      InternalDB <- rRoma::ACSN_Global.v1.1
    } else {
      if(Version == "1.1"){
        print("Searching in ACSN Global map v1.1")
        InternalDB <- rRoma::ACSN_Global.v1.1
      }
    }
  }
  
  if(BDName == "ACSN_Apoptosis"){
    if(is.null(Version)){
      print("Searching in ACSN Apoptosis map v1.1")
      InternalDB <- rRoma::ACSN_Apoptosis.v1.1
    } else {
      if(Version == "1.1"){
        print("Searching in ACSN Apoptosis map v1.1")
        InternalDB <- rRoma::ACSN_Apoptosis.v1.1
      }
    }
  }
  
  if(BDName == "ACSN_CellCycle"){
    if(is.null(Version)){
      print("Searching in ACSN Cell Cycle map v1.1")
      InternalDB <- rRoma::ACSN_CellCycle.v1.1
    } else {
      if(Version == "1.1"){
        print("Searching in ACSN Cell Cycle map v1.1")
        InternalDB <- rRoma::ACSN_CellCycle.v1.1
      }
    }
  }
  
  if(BDName == "ACSN_DNARepair"){
    if(is.null(Version)){
      print("Searching in ACSN DNA repair map v1.1")
      InternalDB <- rRoma::ACSN_DNARepair.v1.1
    } else {
      if(Version == "1.1"){
        print("Searching in ACSN DNA repair map v1.1")
        InternalDB <- rRoma::ACSN_DNARepair.v1.1
      }
    }
  }
  
  if(BDName == "ACSN_EMT"){
    if(is.null(Version)){
      print("Searching in ACSN EMT map v1.1")
      InternalDB <- rRoma::ACSN_EMT.v1.1
    } else {
      if(Version == "1.1"){
        print("Searching in ACSN EMT map v1.1")
        InternalDB <- rRoma::ACSN_EMT.v1.1
      }
    }
  }
  
  if(BDName == "ACSN_Survival"){
    if(is.null(Version)){
      print("Searching in ACSN Survival map v1.1")
      InternalDB <- rRoma::ACSN_Survival.v1.1
    } else {
      if(Version == "1.1"){
        print("Searching in ACSN Survival map v1.1")
        InternalDB <- rRoma::ACSN_Survival.v1.1
      }
    }
  }
  
  if(BDName == "InfoSig_Conserved"){
    if(is.null(Version)){
      print("Searching in InfoSig conserved")
      InternalDB <- rRoma::InfoSig_Conserved
    }
  }
  
  if(BDName == "InfoSig_Informative"){
    if(is.null(Version)){
      print("Searching in InfoSig")
      InternalDB <- rRoma::InfoSig_Informative
    }
  }
  
  if(length(InternalDB) == 0){
    print("No database selected")
    return(NULL)
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




















#' Infer weigths from expression data
#'
#' @param ExpressionMatrix 
#' @param ModuleList 
#' @param FillAllNA 
#'
#' @return
#' @export
#'
#' @examples
InferBinaryWeigth <- function(ExpressionMatrix, ModuleList, FillAllNA = TRUE) {
  
  MedianExpr <- apply(ExpressionMatrix, 1, median, na.rm=TRUE)
  GeneWei <- as.integer(MedianExpr >= median(ExpressionMatrix))
  GeneWei[GeneWei == 0] <- -1
  names(GeneWei) <- names(MedianExpr)
  
  for(i in 1:length(ModuleList)){
    
    if(!FillAllNA & any(!is.na(ModuleList[[i]]$Weigths))){
      next
    }
    
    ModuleList[[i]]$Weigths[is.na(ModuleList[[i]]$Weigths)] <-
      GeneWei[ModuleList[[i]]$Genes[is.na(ModuleList[[i]]$Weigths)]]
    
    if(any(is.na(ModuleList[[i]]$Weigths))){
      print(paste("Warning:", sum(is.na(ModuleList[[i]]$Weigths)), "gene(s) not found in the expression matrix"))
    }
    
  }
  
  return(ModuleList)
  
}

