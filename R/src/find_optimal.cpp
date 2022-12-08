#include "Rcpp.h"
#include "oiff/oiff.hpp"

// [[Rcpp::export(rng=false)]]
Rcpp::List find_optimal_filter(Rcpp::NumericVector pvalues, Rcpp::NumericVector covariates, double threshold, bool larger) {
    size_t n = pvalues.size();
    if (n != static_cast<size_t>(covariates.size())) {
        throw std::runtime_error("'pvalues' and 'covariates' should have the same length");
    }

    oiff::OptimizeFilter runner;
    runner.fdr_threshold = threshold;
    runner.retain_larger = larger;
    auto res = runner.run(n, static_cast<const double*>(pvalues.begin()), static_cast<const double*>(covariates.begin()));

    return Rcpp::List::create(
        Rcpp::Named("first") = Rcpp::NumericVector::create(res.first),
        Rcpp::Named("last") = Rcpp::NumericVector::create(res.last),
        Rcpp::Named("middle") = Rcpp::NumericVector::create(res.middle),
        Rcpp::Named("number") = Rcpp::IntegerVector::create(res.number)
    );
}

// [[Rcpp::export(rng=false)]]
Rcpp::List find_optimal_filter_subsample(Rcpp::NumericVector pvalues, Rcpp::NumericVector covariates, double threshold, bool larger, double subsample_proportion, int num_iterations, int random_seed, int num_threads) {
    size_t n = pvalues.size();
    if (n != static_cast<size_t>(covariates.size())) {
        throw std::runtime_error("'pvalues' and 'covariates' should have the same length");
    }

    oiff::OptimizeFilter runner;
    runner.fdr_threshold = threshold;
    runner.retain_larger = larger;
    runner.subsample_proportion = subsample_proportion;
    runner.num_iterations = num_iterations;
    runner.random_seed = random_seed;
    runner.num_threads = num_threads;
    auto res = runner.run_subsample(n, static_cast<const double*>(pvalues.begin()), static_cast<const double*>(covariates.begin()));

    Rcpp::NumericVector first_thresholds(num_iterations), last_thresholds(num_iterations), middle_thresholds(num_iterations);
    Rcpp::IntegerVector number(num_iterations);
    for (int i = 0; i < num_iterations; ++i) {
        first_thresholds[i] = res[i].first;
        last_thresholds[i] = res[i].last;
        middle_thresholds[i] = res[i].middle;
        number[i] = res[i].number;
    }

    return Rcpp::List::create(
        Rcpp::Named("first") = first_thresholds,
        Rcpp::Named("last") = last_thresholds,
        Rcpp::Named("middle") = middle_thresholds,
        Rcpp::Named("number") = number
    );
}

