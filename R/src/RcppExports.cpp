// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// find_optimal_filter
Rcpp::List find_optimal_filter(Rcpp::NumericVector pvalues, Rcpp::NumericVector covariates, double threshold);
RcppExport SEXP _oiff_find_optimal_filter(SEXP pvaluesSEXP, SEXP covariatesSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pvalues(pvaluesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type covariates(covariatesSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(find_optimal_filter(pvalues, covariates, threshold));
    return rcpp_result_gen;
END_RCPP
}
// find_optimal_filter_subsample
Rcpp::List find_optimal_filter_subsample(Rcpp::NumericVector pvalues, Rcpp::NumericVector covariates, double threshold, double subsample_proportion, int num_iterations, int random_seed, int num_threads);
RcppExport SEXP _oiff_find_optimal_filter_subsample(SEXP pvaluesSEXP, SEXP covariatesSEXP, SEXP thresholdSEXP, SEXP subsample_proportionSEXP, SEXP num_iterationsSEXP, SEXP random_seedSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pvalues(pvaluesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type covariates(covariatesSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< double >::type subsample_proportion(subsample_proportionSEXP);
    Rcpp::traits::input_parameter< int >::type num_iterations(num_iterationsSEXP);
    Rcpp::traits::input_parameter< int >::type random_seed(random_seedSEXP);
    Rcpp::traits::input_parameter< int >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(find_optimal_filter_subsample(pvalues, covariates, threshold, subsample_proportion, num_iterations, random_seed, num_threads));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_oiff_find_optimal_filter", (DL_FUNC) &_oiff_find_optimal_filter, 3},
    {"_oiff_find_optimal_filter_subsample", (DL_FUNC) &_oiff_find_optimal_filter_subsample, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_oiff(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
