#include <gtest/gtest.h>
#include "oiff/rank.hpp"

#include <vector>
#include <random>

TEST(RankTest, Basic) {
    std::vector<double> input { 0.2, -1, 3.2, 3, -0.1, 1.2 };
    auto ranked = oiff::rankify(input.begin(), input.size());

    std::vector<int> expected{ 2, 0, 5, 4, 1, 3 };
    EXPECT_EQ(ranked, expected);
}

TEST(RankTest, Ties) {
    std::vector<double> input { 0.2, -1, 1.2, 0.2, -1, 0.2, 1.2 };
    auto ranked = oiff::rankify(input.begin(), input.size());

    std::vector<int> expected{ 1, 0, 2, 1, 0, 1, 2 };
    EXPECT_EQ(ranked, expected);
}

TEST(RankTest, Empty) {
    auto ranked = oiff::rankify((double*)NULL, 0);
    EXPECT_EQ(ranked.size(), 0);
}
