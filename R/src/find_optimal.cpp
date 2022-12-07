#include "Rcpp.h"
#include "oiff/oiff.hpp"

// [[Rcpp::export(rng=false)]]
Rcpp::List find_optimal_filter(Rcpp::NumericVector pvalues, Rcpp::NumericVector covariates, double threshold) {
    size_t n = pvalues.size();
    if (n != static_cast<size_t>(covariates.size())) {
        throw std::runtime_error("'pvalues' and 'covariates' should have the same length");
    }

    oiff::OptimizeFilter runner;
    runner.fdr_threshold = threshold;
    auto res = runner.run(n, static_cast<const double*>(pvalues.begin()), static_cast<const double*>(covariates.begin()));

    return Rcpp::List::create(
        Rcpp::Named("threshold") = Rcpp::NumericVector::create(res.first),
        Rcpp::Named("number") = Rcpp::IntegerVector::create(res.second)
    );
}

// [[Rcpp::export(rng=false)]]
Rcpp::List find_optimal_filter_subsample(Rcpp::NumericVector pvalues, Rcpp::NumericVector covariates, double threshold, double subsample_proportion, int num_iterations, int random_seed, int num_threads) {
    size_t n = pvalues.size();
    if (n != static_cast<size_t>(covariates.size())) {
        throw std::runtime_error("'pvalues' and 'covariates' should have the same length");
    }

    oiff::OptimizeFilter runner;
    runner.fdr_threshold = threshold;
    runner.num_iterations = num_iterations;
    runner.subsample_proportion = subsample_proportion;
    runner.random_seed = random_seed;
    runner.num_threads = num_threads;
    auto res = runner.run_subsample(n, static_cast<const double*>(pvalues.begin()), static_cast<const double*>(covariates.begin()));

    Rcpp::NumericVector thresholds(num_iterations);
    Rcpp::IntegerVector number(num_iterations);
    for (int i = 0; i < num_iterations; ++i) {
        thresholds[i] = res[i].first;
        number[i] = res[i].second;
    }

    return Rcpp::List::create(
        Rcpp::Named("threshold") = thresholds,
        Rcpp::Named("number") = number
    );
}

