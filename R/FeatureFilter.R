#' Select features/genes/probes based on a pre-define vector.
#'
#' This function will subset expression data based on a pre-define vector of genes or probes.
#'
#' @param GEOdata either a matrix or SummerizedExperiment object contains expression.
#' @param filter a vector of genes/probes identities one wish to keep
#' @param varFilter a character for the variable that is being subsetted
#' @param unique a boolean, after filter, whether to keep only those
#'
#' @return A SummerizedExperiment object contains both expression and metadata related to each sample.
#' @return Or a matrix of expression data.
#'
#' @import SummarizedExperiment
#'
#' @examples
#'
#' GeoExtract(GEOAccession = "GSE48762", getGPL = TRUE,  plots = TRUE, NSub = 20)
#'
