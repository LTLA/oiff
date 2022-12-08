#include <gtest/gtest.h>
#include "oiff/oiff.hpp"

#include <random>
#include <vector>
#include <algorithm>

class ScenarioTest : public ::testing::TestWithParam<std::tuple<int, int> > {
protected:
    typedef oiff::OptimizeFilter::Thresholds<double> Thresholds;

    std::tuple<double, double, double, int> collapse(const Thresholds& x) {
        return std::make_tuple(x.first, x.last, x.last, x.number);
    }

    static Thresholds reference(const std::vector<double>& pvalues, const std::vector<double>& covariates, double threshold) {
        size_t nobs = pvalues.size();
        std::vector<std::pair<double, double> > combined;
        combined.reserve(nobs);
        for (size_t i = 0; i < nobs; ++i) {
            combined.emplace_back(pvalues[i], covariates[i]);
        }
        std::sort(combined.begin(), combined.end());

        int max_hits = -1;
        std::vector<double> best;

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

            if (hits > max_hits) {
                best.clear();
                max_hits = hits;
                best.push_back(current);
            } else if (hits == max_hits) {
                best.push_back(current);
            }
        }

        std::sort(best.begin(), best.end());

        Thresholds output;
        output.first = best.front();
        output.last = best.back();
        output.middle = best[best.size() / 2];
        output.number = max_hits;
        return output;
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

    auto best = oiff::OptimizeFilter().run(nobs, pvalues.data(), covariates.data());
    auto ref = reference(pvalues, covariates, 0.05);
    EXPECT_EQ(collapse(best), collapse(ref));
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

    auto best = oiff::OptimizeFilter().run(nobs, pvalues.data(), covariates.data());
    auto ref = reference(pvalues, covariates, 0.05);
    EXPECT_EQ(collapse(best), collapse(ref));
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

    auto best = oiff::OptimizeFilter().run(nobs, pvalues.data(), covariates.data());
    auto ref = reference(pvalues, covariates, 0.05);
    EXPECT_EQ(collapse(best), collapse(ref));
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

    auto best = oiff::OptimizeFilter().run(nobs, pvalues.data(), covariates.data());
    auto ref = reference(pvalues, covariates, 0.05);
    EXPECT_EQ(collapse(best), collapse(ref));
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

    oiff::OptimizeFilter runner;
    auto best = runner.run(nobs, pvalues.data(), covariates.data());
    auto ref = reference(pvalues, covariates, 0.05);
    EXPECT_EQ(collapse(best), collapse(ref));

    // What happpens if we retain the larger values instead?
    runner.retain_larger = true;
    auto revbest = runner.run(nobs, pvalues.data(), covariates.data());

    for (auto& x : covariates) {
        x *= -1;
    }
    auto revref = reference(pvalues, covariates, 0.05);
    revref.first *= -1;
    revref.middle *= -1;
    revref.last *= -1;

    EXPECT_EQ(collapse(revbest), collapse(revref));
}

TEST_P(ScenarioTest, ManyZeros) {
    auto param = GetParam();
    size_t nobs = std::get<0>(param);
    int seed = std::get<1>(param);

    std::mt19937_64 rng(seed);
    std::uniform_real_distribution udist;
    std::normal_distribution ndist;

    // Spike in with loads of zeros as forced true positives.
    size_t nulls = nobs * 0.8;
    std::vector<double> pvalues(nobs - nulls), covariates(nobs - nulls, 10);
    for (size_t i = 0; i < nulls; ++i) {
        pvalues.push_back(udist(rng));
        covariates.push_back(ndist(rng));
    }

    oiff::OptimizeFilter runner;
    auto best = runner.run(nobs, pvalues.data(), covariates.data());
    auto ref = reference(pvalues, covariates, 0.05);
    EXPECT_EQ(collapse(best), collapse(ref));
}

TEST_P(ScenarioTest, Ties) {
    auto param = GetParam();
    size_t nobs = std::get<0>(param);
    int seed = std::get<1>(param);

    std::mt19937_64 rng(seed);
    std::uniform_real_distribution udist;
    std::normal_distribution ndist;

    std::vector<double> pvalues, covariates;
    size_t nulls = nobs * 0.75;

    // Smaller p-values associated with higher covariates, but injecting ties.
    for (size_t i = 0; i < nobs - nulls; ++i) {
        if (i && udist(rng) <= 0.1) {
            pvalues.push_back(pvalues.back());
            covariates.push_back(covariates.back());
        } else {
            pvalues.push_back(udist(rng)/50);
            covariates.push_back(ndist(rng) + 1);
        }
    }

    // Adding lots of tied values.
    for (size_t i = 0; i < nulls; ++i) {
        if (i && udist(rng) <= 0.1) {
            pvalues.push_back(pvalues.back());
            covariates.push_back(covariates.back());
        } else {
            pvalues.push_back(udist(rng));
            covariates.push_back(ndist(rng));
        }
    }

    oiff::OptimizeFilter runner;
    auto best = runner.run(nobs, pvalues.data(), covariates.data());
    auto ref = reference(pvalues, covariates, 0.05);
    EXPECT_EQ(collapse(best), collapse(ref));
}

TEST_P(ScenarioTest, Subsampling) {
    auto param = GetParam();
    size_t nobs = std::get<0>(param);
    int seed = std::get<1>(param);

    // Just re-using a left-skewed distribution.
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution udist;
    std::normal_distribution ndist;

    std::vector<double> pvalues, covariates;
    size_t most = nobs * 0.5;
    for (size_t i = 0; i < most; ++i) {
        pvalues.push_back(udist(rng));
        covariates.push_back(ndist(rng));
    }

    for (size_t i = 0; i < nobs - most; ++i) {
        pvalues.push_back(udist(rng)/50);
        covariates.push_back(ndist(rng) - 1);
    }

    // Checking against the reference just in case.
    auto best = oiff::OptimizeFilter().run(nobs, pvalues.data(), covariates.data());
    auto ref = reference(pvalues, covariates, 0.05);
    EXPECT_EQ(collapse(best), collapse(ref));

    // Checking against subsampling with 100% retention.
    {
        oiff::OptimizeFilter runner;
        runner.subsample_proportion = 1;
        runner.num_iterations = 1;
        auto res = runner.run_subsample(nobs, pvalues.data(), covariates.data());
        EXPECT_EQ(res.size(), 1);
        EXPECT_EQ(collapse(best), collapse(res[0]));
    }

    // Checking against subsampling with 50% retention.
    {
        oiff::OptimizeFilter runner;
        runner.subsample_proportion = 0.5;
        auto res = runner.run_subsample(nobs, pvalues.data(), covariates.data());
        EXPECT_EQ(res.size(), runner.num_iterations);

        std::mt19937_64 rng(runner.random_seed);
        size_t keep = std::ceil(nobs * runner.subsample_proportion);
        std::vector<uint8_t> selected(nobs);
        std::vector<double> subp, subc;

        for (int it = 0; it < runner.num_iterations; ++it) {
            oiff::sample(nobs, keep, selected, rng);

            subp.clear();
            subc.clear();
            for (size_t i = 0; i < nobs; ++i) {
                if (selected[i]) {
                    subp.push_back(pvalues[i]);
                    subc.push_back(covariates[i]);
                }
            }

            auto manual = runner.run(subp.size(), subp.data(), subc.data());
            EXPECT_EQ(collapse(manual), collapse(res[it]));
        }

        runner.num_threads = 3;
        auto par = runner.run_subsample(nobs, pvalues.data(), covariates.data());
        for (int it = 0; it < runner.num_iterations; ++it) {
            EXPECT_EQ(collapse(par[it]), collapse(res[it]));
        }
    }
}

INSTANTIATE_TEST_SUITE_P(
    Scenarios,
    ScenarioTest,
    ::testing::Combine(
        ::testing::Values(100, 500, 1000), // number of observations.
        ::testing::Values(123456, 789, 0) // seed
    )
);
