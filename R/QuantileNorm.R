#' Quantile normalization on a given eSet data.
#'
#' This function will quantile normalize expression from a SummerizedExperiment object.
#' Then produce a boxplots of NSub number of randomly selected samples.
#'
#' @param GEOdata either a matrix or SummerizedExperiment object contains expression.
#' @param plots boolean, whether plot boxplots of random chosen sample or not
#' @param NSub integer, number of random samples should be plotted.
#'
#' @return Either a matrix or SummerizedExperiment object contains both quantile normalized expression and metadata related to each sample.
#'
#' @import limma
#' @import SummarizedExperiment
#'
#' @examples
#'
#' QuantileNorm(GEOData = GSE48762, plots = TRUE, NSub = 10)
#'
#' @export
QuantileNorm <- function(GEOdata,
                         plots = TRUE,
                         NSub = 10){

  # Check GEOdata class
  if (class(GEOdata)[1] != "SummarizedExperiment" & !is.matrix(GEOdata)){
    stop(paste0("Wrong format of data. Must provide SummarizedExperiment object or Matrix.", GEOdata, " is ", class(GEOdata), "."))
  }

  # Check number of samples
  if (NSub > ncol(GEOdata) & plots) {
    NSub <- ncol(GEOdata)
    print(paste0("Only ",ncol(GEOdata), " sample(s) available for boxplots."))}

  # Quantile normalization based on type of data
  if(class(GEOdata)[1] == "SummarizedExperiment"){
    assays(GEOdata)$quantileNorm <- limma::normalizeQuantiles(assays(GEOdata)$exprs)

    boxplot(assays(GEOdata)$quantileNorm[, sample(1:ncol(assays(GEOdata)$quantileNorm), NSub, replace = FALSE)],
            main = "Distribution of 10 random samples after quantile normalization.")
  }

  else if(is.matrix(GEOdata)) {
    GEOdata <- limma::normalizeQuantiles(GEOdata)

    boxplot(GEOdata[, sample(1:ncol(GEOdata), NSub, replace = FALSE)],
            main = "Distribution of 10 random samples after quantile normalization.")
  }

  return(GEOdata)

}
