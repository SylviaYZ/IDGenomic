#' Select features/genes/probes based on a pre-define vector.
#'
#' This function will subset expression data based on a pre-define vector of genes or probes.
#'
#' @param GEOdata is either a matrix or SummerizedExperiment object. If a SummerizedExperiment object is provided, GEOdata contains both expression and metadata related to each sample and feature. If a matrix is provided, GEOdata contains expression data, and the `rownames(GEOdata)` must correspond to the `filterValue`.
#' @param filterVar should be a character indicating the variable that is being subsetted related to the feature of expression data for SummerizedExperiment object. The default value for filterVar is NA. Note that if a matrix is provided as GEOdata, there is no need to provide `filterVar`.
#' @param filterValue is a vector of gene/probe identities that one wishes to keep in the dataset after subsetting. If a SummerizedExperiment object is provided, `filterValue` is a subset of values in `filterVar`. If a matrix is provided as GEOdata, `filterValue` is a subset of `rownames(GEOdata)`.
#' @param unique is a boolean that determines whether to keep only unique genes/probes after filtering. Its default value is FALSE. If `unique=TRUE`, only unique genes/probes are retained by selecting the gene/probe with the highest median expression level.
#'
#' @return A SummerizedExperiment object or matrix.
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
