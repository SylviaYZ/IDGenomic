#' Classification using Gaussian mixture model with EM algorithm.
#'
#' This function classify samples/cells into different pre-specified groups,
#' where mean and sd for Gaussian distribution are given as reference/initial input.
#' This is a supervised clustering method.
#'
#' @param expression a gene expression matrix, row = G genes, column = N samples/cells
#' @param ref_mean mean each marker for each cluster, row = L markers (L<G), column = K clusters
#' @param ref_sd sd each marker for each cluster, row = L markers (L<G), column = K clusters
#' @param subtype a vector contains names of clusters
#' @param true_label a vector contains actual cluster label of each sample/cells
#'
#' @return A list of three items:
#' @return 1: a matrix of probability of each sample belongs to each cluster
#' @return 2: assignment based on maximum likelihood
#' @return 3: EM algorithm iteration
#'
EMGMM <- function(expression,
                  ref_mean,
                  ref_sd,
                  subtype,
                  true_label){

  ######################################################
  # Consider get rid of subtype parameter, we can get it from colnames(ref_mean)
  #######################################################

  ref_mean <-  ref_mean[match(rownames(expression), rownames(ref_mean)), match(subtype, colnames(ref_mean))]
  ref_sd <-  ref_sd[match(rownames(expression), rownames(ref_sd)),match(subtype, colnames(ref_sd))]
  mean_diff <- matrix(1, nrow = nrow(ref_mean), ncol = ncol(ref_mean))
  i <- 0

  #####################################################
  # MODIFY: make change_iter to data.frame, then update it with index.
  change_iter <- matrix(NA, ncol=3)
  colnames(change_iter) <- c("num_iter", "max_mean_diff", "incorrect_label")
  ####################################################



  while( max(mean_diff) > 0.001 && i < 50 ){

    # E-step
    ref_sd <- ifelse(as.matrix(ref_sd) == 0, 0.0001, as.matrix(ref_sd))
    w_try <- apply(expression, 2, function(exp_i) (sapply(c(1:length(subtype)), function(i) sum(log10(dnorm(exp_i, ref_mean[,i], ref_sd[,i]))))))
    rownames(w_try) <- subtype
    w_try <- apply(w_try, 2, function(prob_i) prob_i + ifelse(max(prob_i) < 0 , abs(max(prob_i)),  -max(prob_i)))
    prob_matrix <- t(apply(w_try, 2, function(i) (10^i)/ sum((10^i))))
    colnames(prob_matrix) <- subtype
    prob_matrix <- ifelse(is.na(prob_matrix) == TRUE, 0, prob_matrix)
    prob_max <- apply(prob_matrix, 1,  function(cell) max(cell))
    prob_max_assign <- apply(prob_matrix, 1,  function(cell) ifelse((all(cell == 0) | all(is.na(cell))), "Unknown", subtype[which.max(cell)]))

    colnames(prob_matrix) <- subtype

    # M-step
    mean_update <- matrix(nrow = nrow(expression), ncol = length(subtype))
    sd_update <- matrix(nrow = nrow(expression), ncol = length(subtype))
    for (type_i in 1:length(subtype)){
      mean_update[,type_i] <- as.vector(apply(expression, 1, function(gene_i)sum(gene_i*prob_matrix[,type_i ])/ sum(prob_matrix[,type_i ])))
      sd_update[,type_i] <- sapply( c(1: nrow(expression)), function(gene_i) sqrt(sum( prob_matrix[,type_i ] *(expression[gene_i,] -  ref_mean[gene_i,type_i])^2  )/sum(prob_matrix[,type_i ]) ))
    }
    colnames(mean_update) <- subtype
    rownames(mean_update) <- rownames(expression)
    colnames(sd_update) <- subtype
    rownames(sd_update) <- rownames(expression)

    # update parameter
    mean_diff <- abs(mean_update - ref_mean)
    ref_mean <- mean_update
    ref_sd <- sd_update
    i <- i+1

    #####################################################
    # MODIFY: make change_iter to data.frame, then update it with index.
    change_iter <-rbind(change_iter,c(i,max(mean_diff),sum(prob_max_assign != true_label)))
    #####################################################
  }

  return_list <- list(prob_matrix,prob_max_assign,change_iter)
  return(return_list)
}

