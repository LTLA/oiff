#' Find the optimal filter threshold
#'
#' Find the threshold on an independent filter statistic that maximizes the number of discoveries
#' after applying the Benjamini-Hochberg method to control the false discovery rate.
#'
#' @param pvalues Numeric vector of p-values, one for each hypothesis.
#' @param filter Numeric vector of length equal to \code{pvalues}, 
#' containing the filter statistic for each hypothesis.
#' @param above Logical scalar indicating whether to retain hypotheses with \code{filter} above the threshold.
#' By default, lower values of \code{filter} are retained when choosing a filter statistic.
#' @param threshold Numeric scalar specifying the FDR threshold to use.
#' @param subsample Numeric scalar between 0 and 1 specifying the proportion of hypotheses to use for subsampling.
#' The optimal filter threshold is then estimated from the sampled subset instead of the full dataset.
#' This mitigates loss of FDR control from p-value-dependent adjustment of the filter threshold.
#' If \code{NULL}, no subsampling is performed.
#' @param iterations Integer scalar specifying the number of subsampling iterations to perform.
#' The returned filter threshold is defined as the mean of the optimal thresholds across all iterations.
#' Larger values improve the stability of the subsampled estimate at the cost of time.
#' Only used if \code{subsample} is not \code{NULL}.
#' @param num.threads Integer scalar specifying the number of threads to use for subsampling.
#' Only used if \code{subsample} is not \code{NULL}.
#'
#' @return A list containing:
#' \itemize{
#' \item \code{first}, the first (i.e., most stringent) filter threshold that achieves the maximum number of discoveries.
#' \item \code{last}, the last (i.e., most relaxed) filter threshold that achieves the maximum number of discoveries.
#' \item \code{middle}, an intermediate threshold between \code{first} and \code{last},
#' approximately the median of all thresholds that achieve the maximal number of discoveries.
#' In the presence of many tied thresholds, this is probably the best choice as it is more robust to the number of hypotheses.
#' \item \code{number}, the maximum number of discoveries.
#' }
#'
#' When subsampling, \code{first}, \code{last} and \code{middle} are the means of the corresponding values across iterations.
#' \code{number} is not reported as it is not guaranteed to be the same for the different means.
#' 
#' @author Aaron Lun
#'
#' @examples
#' pvalues <- c(runif(9900), rbeta(100, 1, 500))
#' filter <- c(rnorm(9900), rnorm(100) + 2)
#' findOptimalFilter(pvalues, filter, above=TRUE)
#' findOptimalFilter(pvalues, filter, above=TRUE, subsample=0.1)
#'
#' @export
#' @importFrom Rcpp sourceCpp
#' @useDynLib oiff
findOptimalFilter <- function(pvalues, filter, above=FALSE, threshold=0.05, subsample=NULL, iterations=100, num.threads=1) {
    args <- list(pvalues=pvalues, covariates=filter, threshold=threshold, larger=above)

    if (is.null(subsample)) {
        res <- do.call(find_optimal_filter, args)
    } else {
        args$subsample_proportion <- subsample
        args$num_iterations <- iterations
        args$random_seed <- sample(.Machine$integer.max, 1)
        args$num_threads <- num.threads
        sub <- do.call(find_optimal_filter_subsample, args)
        res <- lapply(sub[c("first", "last", "middle")], mean)
    }

    res
}
