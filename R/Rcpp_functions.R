#' Calculate Spearman correlation given two matrices
#'
#' @param x matrix of samples. (Row: gene; column: samples)
#' @param y matrix of samples. (Row: gene; column: samples)
#'
#' @return A matrix.
#'
#' @useDynLib IDGenomic
#' @export
corrSpearman <- function( x, y){
  result <- corr_Spearman(x,y)
  rownames(result) <- colnames(x)
  colnames(result) <- colnames(y)
  return(result)
}

