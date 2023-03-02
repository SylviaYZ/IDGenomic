#' Correlation plot
#'
#' This plot will generate correlation between samples, or correlation between genes.
#'
#' @param GEOdata either a matrix or SummerizedExperiment object.
#' If a SummarizedExperiment object is provided, it already contains the expression data and its corresponding metadata.
#' If a matrix is provided, rows are features and columns are samples (this structure is maintained throughout the entire package).
#' @param slot.use specify asasy data in GEOdata to calculate pair-wise correlation, default is "exprs"
#' @param corr.method method to calculate correlation c("pearson","spearman"), default is "spearman"
#' @param plot.method type of plot for output, c("heatmap", "histogram"). Heatmap is only recommended when sample size is less than 30.
#' @param corr.cutoff cutoff for exporting sample ids with correlation less than this cutoff
#'
#' @return Plot(s) of correlation.
#' Heatmap is drawn using [corrplot](https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html) package.
#'
#' @import SummarizedExperiment
#' @import corrplot
#'
#' @importFrom Hmisc rcorr
#'
#' @examples
#'
#' CorrPlot(GEOdata = GSE48762[,301:330], slot.use = "exprs", corr.method = "spearman", plot.method = c("heatmap", "hist"))
#'
#' @export
CorrPlot <- function(GEOdata,
                     slot.use = "exprs",
                     corr.method = "spearman",
                     plot.method = "histogram",
                     corr.cutoff = 0.5){

  # Check GEOdata class
  if (class(GEOdata)[1] != "SummarizedExperiment" & !is.matrix(GEOdata)){
    stop(paste0("Wrong format of data. Must provide SummarizedExperiment object or Matrix.", GEOdata, " is ", class(GEOdata), "."))
  }

  # Calculate correlation for SummarizedExperiment object
  if(class(GEOdata)[1] == "SummarizedExperiment"){
    sample.corr <- Hmisc::rcorr(as.matrix(assays(GEOdata)[[slot.use]]),  type=corr.method)
  }
  else if(is.matrix(GEOdata)) {
    sample.corr <- Hmisc::rcorr(GEOdata, type=corr.method)
  }

  # Set plot grid
  par(mfrow = c(1, length(plot.method)))
  # plot
  if ("heatmap" %in% plot.method){

    if (nrow(sample.corr$r) > 30) warning("Heatmap may be distorted with more than 30 samples.")

    corrplot::corrplot(as.matrix(sample.corr$r), order =  "hclust",
                       p.mat = sample.corr$p, sig.level = 0.05,
                       method = 'square',
                       tl.col = "black", tl.srt = 60,
                       cl.cex = 1, col = corrplot::COL2('PRGn'),
                       type = 'lower',diag = FALSE, mar=c(0,0,3,0),
                       main =  paste0("Pair-wise ", corr.method, " correlation"))
  }
  if ("histogram" %in% plot.method){
    hist(sample.corr$r[upper.tri(sample.corr$r, diag = FALSE)],
         xlab = "Pair-wise correlation",
         main = paste0("Pair-wise ", corr.method, " correlation"))

  }

  # colmean( , diag = FALSE)
  # RNA-seq: boxplot(logTransformed data)

}
