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
#'
#' @return A list containing \code{threshold}, the threshold to apply to \code{filter};
#' and \code{number}, the number of discoveries after applying the BH method at this filter threshold.
#' 
#' @author Aaron Lun
#'
#' @examples
#' pvalues <- c(runif(9900), rbeta(100, 1, 50))
#' filter <- c(rnorm(9900), rnorm(100) + 2)
#' findOptimalFilter(pvalues, filter, above=TRUE)
#'
#' @export
#' @importFrom Rcpp sourceCpp
#' @useDynLib oiff
findOptimalFilter <- function(pvalues, filter, above=FALSE, threshold=0.05) {
    if (above) {
        filter <- -filter;
    }

    res <- find_optimal_filter(pvalues, filter, threshold)

    if (above) {
        res$threshold <- -res$threshold
    }

    res
}
