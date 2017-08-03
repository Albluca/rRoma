#' Conversion between gene naming convenctions and organisms
#'
#' @param SourceOrganism Source organism (e.g., "hsapens" or "mmusculus"). Default is "hsapens".
#' @param TargetOrganism Target organism (e.g., "hsapens" or "mmusculus"). Default is "mmusculus".
#' @param Genes Vector of strings containing gene names. NAs can be present and will be removed
#' @param SourceTypes Type of gene names used as input (currently either "Names" or "Ensembl"). Default is "Names".
#' @param TargetTypes Type of gene names to be returned (currently either "Names" or "Ensembl"). Default is "Names".
#' @param HomologyLevel minimal level of homology (0 by default)
#' @param HOST alterntive host is the default one doe snot work
#' @param PATH path of the alterntive host
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#'
ConvertNames <- function(SourceOrganism = "hsapiens",
                         TargetOrganism = "mmusculus",
                         Genes,
                         SourceTypes = "Names",
                         TargetTypes = "Names",
                         HomologyLevel = 0,
                         HOST = NULL,
                         PATH = NULL) {

  if(is.null(PATH) & !is.null(HOST)){
    stop("PATH need to be specified if HOST is not NULL")
  } else {
    if(is.null(HOST)){
      Mart <- biomaRt::useMart("ensembl", dataset = paste(SourceOrganism, "_gene_ensembl", sep = ''))
    } else {
      Mart <- biomaRt::useMart("ensembl", dataset = paste(SourceOrganism, "_gene_ensembl", sep = ''),
                               host = HOST, path = PATH)
    }
  }
  
  
  
  

  ToList <- FALSE

  if(!is.list(Genes)){
    ToList <- TRUE
    Genes <- list(Genes)
  }

  if(SourceTypes == "Names"){
    FilterName = "external_gene_name"
  }

  if(SourceTypes == "Ensembl"){
    FilterName = "ensembl_gene_id"
  }
  
  if(SourceOrganism != TargetOrganism){
    
    CovGene <- lapply(Genes, function(x){
      x <- x[!is.na(x)]
      biomaRt::getBM(attributes = c(paste(TargetOrganism, "_homolog_orthology_confidence", sep=''),
                                    paste(TargetOrganism, "_homolog_ensembl_gene", sep=''),
                                    paste(TargetOrganism, "_homolog_associated_gene_name", sep=''),
                                    FilterName),
                     filters = FilterName,
                     values = x,
                     mart = Mart)
    })
    
    if(TargetTypes == "Names"){
      ColName = paste(TargetOrganism, "_homolog_associated_gene_name", sep='')
    }
    
    if(TargetTypes == "Ensembl"){
      ColName = paste(TargetOrganism, "_homolog_ensembl_gene", sep='')
    }
    
    RetGene <- lapply(CovGene, function(x){
      if(nrow(x) == 0){
        return(NULL)
      }
      DetectedHomology = x[,paste(TargetOrganism, "_homolog_orthology_confidence", sep='')]
      NewNames <- x[(DetectedHomology  >= HomologyLevel) & !is.na(DetectedHomology), ColName]
      names(NewNames) <- x[(DetectedHomology  >= HomologyLevel) & !is.na(DetectedHomology), FilterName]
      return(NewNames)
    })
    
  } else {
    
    CovGene <- lapply(Genes, function(x){
      x <- x[!is.na(x)]
      biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", FilterName),
                     filters = FilterName,
                     values = x,
                     mart = Mart)
    })
    
    if(TargetTypes == "Names"){
      ColName = "external_gene_name"
    }
    
    if(TargetTypes == "Ensembl"){
      ColName = "ensembl_gene_id"
    }
    
    RetGene <- lapply(CovGene, function(x){
      if(nrow(x) == 0){
        return(NULL)
      }
      
      NewNames <- x[, ColName]
      names(NewNames) <- x[, FilterName]

      return(NewNames)
    })
    
  }
  
  

  if(ToList){
    RetGene <- RetGene[[1]]
  }

  return(RetGene)

}












#' Conversion between gene naming convenctions and organisms of Modulelist
#'
#' @param SourceOrganism Source organism (e.g., "hsapens" or "mmusculus"). Default is "hsapens".
#' @param TargetOrganism Target organism (e.g., "hsapens" or "mmusculus"). Default is "mmusculus".
#' @param SourceTypes Type of gene names used as input (currently either "Names" or "Ensembl"). Default is "Names".
#' @param TargetTypes Type of gene names to be returned (currently either "Names" or "Ensembl"). Default is "Names".
#' @param HomologyLevel minimal level of homology (0 by default)
#' @param ModuleList The module list to be convert.
#' @param HOST alterntive host is the default one doe snot work
#' @param PATH path of the alterntive host
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#'
ConvertModuleNames <- function(
  ModuleList,
  SourceOrganism = "hsapiens",
  TargetOrganism = "mmusculus",
  SourceTypes = "Names",
  TargetTypes = "Names",
  HomologyLevel = 0,
  HOST = NULL,
  PATH = NULL) {
  
  
  NewGenes <- ConvertNames(SourceOrganism = SourceOrganism,
                           TargetOrganism = TargetOrganism,
                           Genes = lapply(ModuleList, "[[", "Genes"),
                           SourceTypes = SourceTypes,
                           TargetTypes = TargetTypes,
                           HomologyLevel = HomologyLevel) 
  
  if(SourceOrganism != TargetOrganism){
    print("Homology-supported gene conversion. Gene weigths will be deleted")
    RetGene <- lapply(as.list(1:length(ModuleList)), function(i){
      ModuleList[[i]]$Genes <- NewGenes[[i]]
      ModuleList[[i]]$Weigths <- rep(NA, length(NewGenes[[i]]))
      return(ModuleList[[i]])
    })
  } else {
    RetGene <- lapply(as.list(1:length(ModuleList)), function(i){
      
      OldWei <- ModuleList[[i]]$Weigths
      names(OldWei) <- ModuleList[[i]]$Genes
      
      OldWei <- OldWei[names(NewGenes[[i]])]
      names(OldWei) <- NULL
      ModuleList[[i]]$Weigths <- OldWei
      
      names(NewGenes[[i]]) <- NULL
      ModuleList[[i]]$Genes <- NewGenes[[i]]
      
      return(ModuleList[[i]])
    })
  }
  
  return(RetGene)
  
}



