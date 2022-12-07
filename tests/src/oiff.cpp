#include <gtest/gtest.h>
#include "oiff/oiff.hpp"

#include <random>
#include <vector>
#include <algorithm>

class ScenarioTest : public ::testing::TestWithParam<std::tuple<int, int> > {
protected:
    static std::pair<double, int> reference(const std::vector<double>& pvalues, const std::vector<double>& covariates, double threshold) {
        size_t nobs = pvalues.size();
        std::vector<std::pair<double, double> > combined;
        combined.reserve(nobs);
        for (size_t i = 0; i < nobs; ++i) {
            combined.emplace_back(pvalues[i], covariates[i]);
        }
        std::sort(combined.begin(), combined.end());

        int max_hits = -1;
        double filter = 0;

        std::vector<int> collected;
        for (size_t i = 0; i < nobs; ++i) {
            auto current = covariates[i];

            collected.clear();
            for (size_t j = 0; j < nobs; ++j) {
                if (combined[j].second <= current) {
                    collected.push_back(j);
                }
            }

            int hits = 0;
            int counter = 0;
            auto adjusted_threshold = threshold / collected.size();
            for (auto j : collected) {
                ++counter;
                if (combined[j].first / counter <= adjusted_threshold) {
                    hits = counter;
                }
            }

            if (hits > max_hits || (hits == max_hits && filter < current)) {
                max_hits = hits;
                filter = current;
            }
        }

        return std::make_pair(filter, max_hits);
    }
};

TEST_P(ScenarioTest, Uniform) {
    auto param = GetParam();
    size_t nobs = std::get<0>(param);
    int seed = std::get<1>(param);

    std::mt19937_64 rng(seed);
    std::uniform_real_distribution udist;
    std::normal_distribution ndist;

    // Simulating a simple scenario with independence
    // between the p-values and the covariates.
    std::vector<double> pvalues, covariates;
    for (size_t i = 0; i < nobs; ++i) {
        pvalues.push_back(udist(rng));
        covariates.push_back(ndist(rng));
    }

    auto best = oiff::find_optimal_filter(nobs, pvalues.data(), covariates.data(), 0.05);
    auto ref = reference(pvalues, covariates, 0.05);
    EXPECT_EQ(best, ref);
}

TEST_P(ScenarioTest, PerfectCorrelation) {
    auto param = GetParam();
    size_t nobs = std::get<0>(param);
    int seed = std::get<1>(param);

    std::mt19937_64 rng(seed);
    std::uniform_real_distribution udist;
    std::normal_distribution ndist;

    std::vector<double> pvalues, covariates;
    for (size_t i = 0; i < nobs; ++i) {
        pvalues.push_back(udist(rng));
        covariates.push_back(ndist(rng));
    }

    // Now the pvalues and covariates are correlated.
    std::sort(pvalues.begin(), pvalues.end());
    std::sort(covariates.begin(), covariates.end());

    auto best = oiff::find_optimal_filter(nobs, pvalues.data(), covariates.data(), 0.05);
    auto ref = reference(pvalues, covariates, 0.05);
    EXPECT_EQ(best, ref);
}

TEST_P(ScenarioTest, ImperfectCorrelation) {
    auto param = GetParam();
    size_t nobs = std::get<0>(param);
    int seed = std::get<1>(param);

    std::mt19937_64 rng(seed);
    std::normal_distribution ndist;

    // Imperfect correlation between p-values and the covariates.
    std::vector<double> pvalues, covariates;
    for (size_t i = 0; i < nobs; ++i) {
        pvalues.push_back(static_cast<double>(i + 1) / (nobs + 1));
        covariates.push_back(ndist(rng) * 0.01 + static_cast<double>(i) / nobs);
    }

    auto best = oiff::find_optimal_filter(nobs, pvalues.data(), covariates.data(), 0.05);
    auto ref = reference(pvalues, covariates, 0.05);
    EXPECT_EQ(best, ref);
}

TEST_P(ScenarioTest, LowerCovariates) {
    auto param = GetParam();
    size_t nobs = std::get<0>(param);
    int seed = std::get<1>(param);

    std::mt19937_64 rng(seed);
    std::uniform_real_distribution udist;
    std::normal_distribution ndist;

    std::vector<double> pvalues, covariates;
    size_t most = nobs * 0.9;
    for (size_t i = 0; i < most; ++i) {
        pvalues.push_back(udist(rng));
        covariates.push_back(ndist(rng));
    }

    // Smaller p-values associated with lower covariates.
    for (size_t i = 0; i < nobs - most; ++i) {
        pvalues.push_back(udist(rng)/50);
        covariates.push_back(ndist(rng) - 1);
    }

    auto best = oiff::find_optimal_filter(nobs, pvalues.data(), covariates.data(), 0.05);
    auto ref = reference(pvalues, covariates, 0.05);
    EXPECT_EQ(best, ref);
}

TEST_P(ScenarioTest, HigherCovariates) {
    auto param = GetParam();
    size_t nobs = std::get<0>(param);
    int seed = std::get<1>(param);

    std::mt19937_64 rng(seed);
    std::uniform_real_distribution udist;
    std::normal_distribution ndist;

    std::vector<double> pvalues, covariates;
    size_t most = nobs * 0.5;
    for (size_t i = 0; i < most; ++i) {
        pvalues.push_back(udist(rng));
        covariates.push_back(ndist(rng));
    }

    // Smaller p-values associated with higher covariates.
    for (size_t i = 0; i < nobs - most; ++i) {
        pvalues.push_back(udist(rng)/50);
        covariates.push_back(ndist(rng) + 1);
    }

    auto best = oiff::find_optimal_filter(nobs, pvalues.data(), covariates.data(), 0.05);
    auto ref = reference(pvalues, covariates, 0.05);
    EXPECT_EQ(best, ref);
}

INSTANTIATE_TEST_SUITE_P(
    Scenarios,
    ScenarioTest,
    ::testing::Combine(
        ::testing::Values(100, 500, 1000), // number of observations.
        ::testing::Values(123456, 789, 0) // seed
    )
);
