#' Quantile normalization on a given eSet data.
#'
#' This function will quantile normalize expression from a SummerizedExperiment object.
#' Then produce a boxplots of NSub number of randomly selected samples.
#'
#' @param GEOdata is either a matrix or SummerizedExperiment object.
#' If a SummarizedExperiment object is provided, it already contains the expression data and its corresponding metadata.
#' If a matrix is provided, rows are features and columns are samples (this structure is maintained throughout the entire package).
#' @param plots boolean, an option thats generates a boxplot of 10 (default) randomly selected samples, this option is set to TRUE by default.
#' @param NSub integer, an option that specifies the number of samples to be randomly selected for the boxplot, with a default value of 10.
#'
#' @return A SummerizedExperiment object contains both quantile-normalized expression and metadata related to each sample, or a matrix of quantile-normalized expression only.
#' For a SummerizedExperiment object, quantile-normalized data is stored in `assays(GEOdata)$quantileNorm`, while the original expression data is preserved.
#'
#' @import limma
#' @import SummarizedExperiment
#'
#' @examples
#'
#' QuantileNorm(GEOdata = GSE48762, plots = TRUE, NSub = 10)
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
