#' Returns all groups
#'
#' @param GrpLevs 
#'
#' @return
#' @export
#'
#' @examples
GetComb <- function(GrpLevs) {
  RetList <- list()
  for(i in 1:length(GrpLevs)){
    for(j in 1:length(GrpLevs)){
      if(i<j){
        RetList[[length(RetList)+1]] <- c(i, j)
      }
    }
  }
  return(RetList)
}