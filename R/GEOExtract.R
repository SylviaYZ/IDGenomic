#' Extract data from GEO.
#'
#' This function will extract data from GEO,
#' then list out available information given in GEO,
#' such as featureData and phenoData.
#' We will look into featureData and phenoData before proceeding for further processing.
#'
#' @param GEOAccession character contains GEO accession
#' @param getGPL boolean, whether extract GPL information or not
#' @param plots boolean, whether plot boxplots of random chosen sample or not
#' @param NSub integer, number of random samples should be plotted.
#'
#' @return One SummerizedExperiment object contains both expression and metadata related to each sample.
#'
#' @import GEOquery
#' @import Biobase
#' @import SummarizedExperiment
#'
#' @examples
#'
#' GeoExtract(GEOAccession = "GSE48762", getGPL = TRUE,  plots = TRUE, NSub = 20)
#'
#' @export
GeoExtract <- function(GEOAccession,
                       getGPL = TRUE,
                       plots = TRUE,
                       NSub = 10){

  if (!is.character(GEOAccession)) stop("GEO accession is invalid. Must provide character.")
  if (NSub <= 0 & plots) stop("NSub is invalid. Must provide positive integer.")

  # Retrieve first data from GEO
  GEOdata <- GEOquery::getGEO(GEO = GEOAccession,
                    GSEMatrix =TRUE,
                    getGPL=getGPL)[[1]]

  GEOdata.se <- SummarizedExperiment::SummarizedExperiment(assays=list(exprs=Biobase::assayData(GEOdata)$exprs),
                                                           rowData=Biobase::fData(GEOdata),
                                                           colData=Biobase::pData(GEOdata))

  # Display information
  print(paste0("Data retrived from GEO for ", GEOAccession))
  show(GEOdata)

  # Display phenotype data
  print(paste0("Available feature data from ", GEOAccession))
  print(fvarMetadata(GEOdata))

  # Display phenotype data
  print(paste0("Available phenotype data from ", GEOAccession))
  print(varMetadata(GEOdata))

  if(plots){

    if (NSub > ncol(GEOdata.se)) {
      NSub <- ncol(GEOdata.se)
      print(paste0("Only ", ncol(GEOdata.se), " sample available for boxplots."))}

    boxplot(SummarizedExperiment::assays(GEOdata.se)$exprs[, sample(1:ncol(SummarizedExperiment::assays(GEOdata.se)$exprs), NSub, replace = FALSE)],
            main = "Distribution of 10 random samples")
  }

  # Return SummarizedExperiment object
  return(GEOdata.se)

}
