#' Select features/genes/probes based on a pre-define vector.
#'
#' This function will subset expression data based on a pre-define vector of genes or probes.
#'
#' @param GEOdata either a matrix or SummerizedExperiment object contains expression.
#' @param filterVar a character for the variable that is being subsetted
#' @param filterValue a vector of genes/probes identities one wish to keep
#' @param unique a boolean, after filter, whether to keep only those
#'
#' @return A SummerizedExperiment object contains both expression and metadata related to each sample.
#' @return Or a matrix of expression data.
#'
#' @import SummarizedExperiment Biobase
#'
#' @examples
#'
#' FeatureFilter(GEOdata = GSE48762, filterVar = "Symbol", filterValue = gene_overlap, unique = TRUE)
#'
#' @export
FeatureFilter <- function(GEOdata,
                         filterVar = NA,
                         filterValue,
                         unique = FALSE){

  if (class(GEOdata)[1] != "SummarizedExperiment" & !is.matrix(GEOdata)){
    stop(paste0("Wrong format of data. Must provide SummarizedExperiment object or Matrix.", GEOdata, " is ", class(GEOdata), "."))
  }

  if (class(GEOdata)[1] == "SummarizedExperiment"){

    if(is.na(filterVar)) stop("Invalid filterVar option for SummarizedExperiment obbject.")

    if(!any(filterValue %in% rowData(GEOdata)[, filterVar]))  stop("Make sure filterVar corresponds to filterValue. Or none of the value in filterValue is in varFilter.")

    GEOdata <- GEOdata[rowData(GEOdata)[,filterVar] %in% filterValue, ]

    if (unique){
      genes.dup <-  unique(rowData(GEOdata)[,filterVar][duplicated(rowData(GEOdata)[,filterVar])])

      if ("quantileNorm" %in% names(assays(GEOdata))){ data.use = "quantileNorm" }
      else{data.use = "exprs"}

      for(g in genes.dup){
        index <- which(rowData(GEOdata)[,filterVar] == g)
        dup.rm <- index[ -which.max(Biobase::rowMedians(assays(GEOdata)[[data.use]][index,])) ]
        GEOdata <- GEOdata[-dup.rm,]
        }

      }
    }


  if (is.matrix(GEOdata) ){

    if(!is.na(filterVar)) stop("For matrix object, do not specify filterVar. Set rownames(GEOdata) = values of filterVar.")

    if(!any(filterValue %in% rownames(GEOdata))) stop("Make sure rownames(GEOdata) corresponds to filterValue. Or none of the value in filterValue is in rownames(GEOdata).")

    GEOdata <- GEOdata[which(rownames(GEOdata) %in% filterValue), ]

    if(unique){
      genes.dup <-  unique(rownames(GEOdata)[duplicated(rownames(GEOdata))])

      for(g in genes.dup){
        index <- which(rownames(GEOdata) == g)
        dup.rm <- index[ -which.max(Biobase::rowMedians(GEOdata[index,])) ]
        GEOdata <- GEOdata[-dup.rm,]
      }
    }
  }

  return(GEOdata)
}
