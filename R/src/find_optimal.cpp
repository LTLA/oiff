#include "Rcpp.h"
#include "oiff/oiff.hpp"

// [[Rcpp::export(rng=false)]]
Rcpp::List find_optimal_filter(Rcpp::NumericVector pvalues, Rcpp::NumericVector covariates, double threshold) {
    size_t n = pvalues.size();
    if (n != static_cast<size_t>(covariates.size())) {
        throw std::runtime_error("'pvalues' and 'covariates' should have the same length");
    }

    auto res = oiff::find_optimal_filter(
        n, 
        static_cast<const double*>(pvalues.begin()), 
        static_cast<const double*>(covariates.begin()),
        threshold
    );

    return Rcpp::List::create(
        Rcpp::Named("threshold") = Rcpp::NumericVector::create(res.first),
        Rcpp::Named("number") = Rcpp::IntegerVector::create(res.second)
    );
}
