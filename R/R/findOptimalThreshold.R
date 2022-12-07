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
#' @return A list containing \code{threshold}, the threshold to apply to \code{filter};
#' and \code{number}, the number of discoveries after applying the BH method at this filter threshold.
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
#' @importFrom stats p.adjust
#' @useDynLib oiff
findOptimalFilter <- function(pvalues, filter, above=FALSE, threshold=0.05, subsample=NULL, iterations=100, num.threads=1) {
    if (above) {
        filter <- -filter;
    }

    if (is.null(subsample)) {
        res <- find_optimal_filter(pvalues, filter, threshold)
    } else {
        sub <- find_optimal_filter_subsample(pvalues, filter, threshold, 
            subsample_proportion = subsample, 
            num_iterations = iterations, 
            random_seed = sample(.Machine$integer.max, 1),
            num_threads = num.threads)

        res <- list(threshold = mean(sub$threshold))
        res$number <- sum(p.adjust(pvalues[filter <= res$threshold], method="BH") <= threshold)
    }

    if (above) {
        res$threshold <- -res$threshold
    }

    res
}
