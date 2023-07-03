#' Assign label based on correlation result
#'
#' Assign label based on correlation result by given a quantile of correlation.
#'
#' @param cor.matrix matrix of correlation where each column contains correlations for each samples. Row corresponds to sample in label.
#' @param sample.label specify asasy data in GEOdata to calculate pair-wise correlation, default is "exprs"
#' @param quantile.cutoff lower bound for quantile of Spearman correlation to be averaged by each label
#'
#' @return A list containing (1) label result (2) corAVG.Top.pct: average correlation of each sample by each label (3) corSD.Top.pct: sd of correlation of each sample by each label
#'
#' @examples
#'
#' CorrAssign1(cor.matrix = spear.cor, sample.label = stage, quantile.cutoff = 0.5)
#'
#' @export
CorrAssign1 <- function(cor.matrix,
                       sample.label,
                       quantile.cutoff){

  # Classification method 1: average top x% correlation from each label

  #Step 1: Top X% correlation values from each label will selected

  #Step 2: Average correlation corresponding to each label will be calculated

  #Step 3: Label is assigned based on highest average correlation

  cor.matrix <- as.matrix(cor.matrix)

  mean.result <- matrix(NA, ncol = length(unique(sample.label)), nrow = ncol(cor.matrix))
  colnames(mean.result) <- unique(sample.label)
  rownames(mean.result) <- colnames(cor.matrix)

  sd.result <- mean.result

  # Compute by column.
  for (i in colnames(cor.matrix)){

    # Create data.frame holds correlation and label of top results

    label_i <- data.frame( cor_i = cor.matrix[,i],
                           label_sample_i = sample.label) %>%
      dplyr::group_by(label_sample_i) %>% dplyr::summarise(meanTop = mean(cor_i[cor_i > quantile(cor_i, quantile.cutoff)]),
                                                           sdTop = sd(cor_i[cor_i > quantile(cor_i, quantile.cutoff)]))

    match_values <- match(label_i$label_sample_i, colnames(mean.result))
    mean.result[i, match_values] <- label_i$meanTop
    sd.result[i, match_values] <- label_i$sdTop
  }

  cor.label <- colnames(mean.result)[apply(mean.result, 1, which.max)]

  return(list(label.result = cor.label,
              corAVG.Top.pct = mean.result,
              corSD.Top.pct = sd.result))

}

#' Assign label based on correlation result
#'
#' Assign label based on correlation result by given a quantile of correlation.
#'
#' @param cor.matrix matrix of correlation where each column contains correlations for each samples. Row corresponds to sample in label.
#' @param sample.label specify asasy data in GEOdata to calculate pair-wise correlation, default is "exprs"
#' @param quantile.cutoff lower bound for quantile of Spearman correlation to be averaged by each label
#'
#' @return A list containing (1) label result (2) avgCor.Top: average correlation of each sample by each label
#'
#' @examples
#'
#' CorrAssign2(cor.matrix = spear.cor, sample.label = stage, quantile.cutoff = 0.5)
#'
#' @export
CorrAssign2 <- function(cor.matrix,
                       sample.label,
                       quantile.cutoff){

  # Classification method 2: top X% correlation selected, then average by label

  # Step 1: Top X% correlation values will selected

  #Step 2: Average correlation corresponding to each label will be calculated

  # Step 3: Label is assigned based on highest average correlation

  cor.matrix <- as.matrix(cor.matrix)

  label.result <- matrix(NA, ncol = length(unique(sample.label)), nrow = ncol(cor.matrix))
  colnames(label.result) <- unique(sample.label)
  rownames(label.result) <- colnames(cor.matrix)

  samples <- colnames(cor.matrix)
  # Compute by column.
  for (i in colnames(cor.matrix)){

    # Generate index to top percentiles
    idx_top <- cor.matrix[,i] > quantile(cor.matrix[,i] , quantile.cutoff)

    # Create data.frame holds correlation and label of top results
    label_i <- data.frame( cor_i = cor.matrix[idx_top, i],
                           label_sample_i = sample.label[idx_top]) %>%
      group_by(label_sample_i) %>% summarise(meanTop = mean(cor_i))

    match_values <- match(label_i$label_sample_i, colnames(label.result))
    label.result[i, match_values] <- label_i$meanTop

  }

  cor.label <- colnames(label.result)[apply(label.result, 1, which.max)]

  return(list(label.result = cor.label,
              cor.Top.pct = label.result))

}
