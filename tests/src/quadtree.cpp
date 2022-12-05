#include <gtest/gtest.h>
#include "oiff/quadtree.hpp"

#include <vector>
#include <random>

class QuadTreeCountTest : public ::testing::TestWithParam<std::tuple<int, int, bool> > {};

TEST_P(QuadTreeCountTest, CountUnder) {
    auto param = GetParam();
    size_t nobs = std::get<0>(param);
    size_t seed = std::get<1>(param);
    bool has_ties = std::get<2>(param);

    std::vector<int> prank;
    std::vector<int> crank;

    if (has_ties) {
        std::mt19937_64 rng(seed);
        std::uniform_real_distribution dist;

        int counter = 0;
        for (size_t i = 0; i < nobs; ++i) {
            if (dist(rng) <= 0.9) { // occasional ties.
                ++counter;
            }
            prank.push_back(counter);
            crank.push_back(counter);
        }
    } else {
        for (size_t i = 0; i < nobs; ++i) {
            prank.push_back(i);
            crank.push_back(i);
        }
    }

    std::mt19937_64 rng(seed);
    std::shuffle(prank.begin(), prank.end(), rng);
    std::shuffle(crank.begin(), crank.end(), rng);

    auto tree = oiff::build_quad_tree(prank, crank);
    EXPECT_TRUE(tree.size() > 100);

    for (size_t i = 0; i < nobs; ++i) {
        auto observed = oiff::count_under(prank[i], crank[i], tree);

        // Manual counting.
        int expected = 0;
        for (size_t j = 0; j < nobs; ++j) {
            expected += (prank[j] <= prank[i] && crank[j] <= crank[i]);
        }

        ASSERT_EQ(observed, expected) << "failed for observation " << i;
    }
}

INSTANTIATE_TEST_SUITE_P(
    QuadTreeCount,
    QuadTreeCountTest,
    ::testing::Combine(
        ::testing::Values(100, 200, 500), // number of observations
        ::testing::Values(0, 12345, 987654321), // seed
        ::testing::Values(false, true) // ties
    )
);

