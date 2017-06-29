#' Conversion between gene naming convenctions and organisms
#'
#' @param SourceOrganism Source organism (e.g., "hsapens" or "mmusculus"). Default is "hsapens".
#' @param TargetOrganism Target organism (e.g., "hsapens" or "mmusculus"). Default is "mmusculus".
#' @param Genes Vector of strings containing gene names. NAs can be present and will be removed
#' @param SourceTypes Type of gene names used as input (currently either "Names" or "Ensembl"). Default is "Names".
#' @param TargetTypes Type of gene names to be returned (currently either "Names" or "Ensembl"). Default is "Names".
#' @param HomologyLevel minimal level of homology (0 by default)
#'
#' @return
#' @export
#'
#' @examples
#'
#' Genes <- c("CCND3", "CD151", "CD2BP2", "CD81", "CDC34", "CDC37", "CDC42BPB", "CDC42EP1","CDC42EP4")
#'
#' ConvertNames("human", "mouse", SourceTypes = "Names", TargetTypes = "Names", GenNames)
#'
#' MoN <- ConvertNames("human", "mouse", SourceTypes = "Names", TargetTypes = "Names", GenNames)
#' MoEn <- ConvertNames("human", "mouse", SourceTypes = "Names", TargetTypes = "Ensembl", GenNames)
#'
#' ConvertNames("mouse", "human", SourceTypes = "Names", TargetTypes = "Names", MoN)
#' ConvertNames("mouse", "human", SourceTypes = "Ensembl", TargetTypes = "Ensembl", MoEn)
#'
#'
ConvertNames <- function(SourceOrganism = "hsapiens",
                         TargetOrganism = "mmusculus",
                         Genes,
                         SourceTypes = "Names",
                         TargetTypes = "Names",
                         HomologyLevel = 0) {

  Mart <- biomaRt::useMart("ensembl", dataset = paste(SourceOrganism, "_gene_ensembl", sep = ''))

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

  CovGene <- lapply(Genes, function(x){
    x <- x[!is.na(x)]
    biomaRt::getBM(attributes = c(paste(TargetOrganism, "_homolog_orthology_confidence", sep=''),
                                  paste(TargetOrganism, "_homolog_ensembl_gene", sep=''),
                                  paste(TargetOrganism, "_homolog_associated_gene_name", sep='')),
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
    if(length(x) <= 3){
      return(NULL)
    }
    DetectedHomology = x[,paste(TargetOrganism, "_homolog_orthology_confidence", sep='')]
    return(x[ColName, (DetectedHomology  >= HomologyLevel) & !is.na(DetectedHomology)])
  })

  if(ToList){
    RetGene <- unlist(RetGene)
  }

  return(RetGene)

}




MList <- SelectFromInternalDB("hall")
ConvNames <- ConvertNames(SourceOrganism = 'hsapiens', TargetOrganism = 'mmusculus',
                          Genes = lapply(MList, "[[", "Genes"))

