#include <gtest/gtest.h>
#include "oiff/oiff.hpp"

#include <random>
#include <vector>

TEST(OiffTest, Basic) {
    size_t nobs = 100;
    std::vector<double> pvalues, covariates;
    std::mt19937_64 rng;
    std::uniform_real_distribution udist;
    std::normal_distribution ndist;

    for (size_t i = 0; i < nobs; ++i) {
        pvalues.push_back(udist(rng));
        covariates.push_back(ndist(rng));
    }

    auto best = oiff::find_optimal_filter(nobs, pvalues.data(), covariates.data(), 0.05);
    std::cout << best << std::endl;
}
