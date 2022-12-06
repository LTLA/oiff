#ifndef OIFF_H
#define OIFF_H

#include <vector>
#include <algorithm>

namespace oiff {

struct Node {
    int p_left, p_mid, p_right;
    int jump_left = -1, jump_right = -1;
};

void build(int p, int location, std::vector<Node>& tree) {
    auto& current = tree[location];
    if (current.p_right - current.p_left == 1) {
        return;
    }

    int jump;
    if (p < current.p_mid) {
        jump = current.jump_left;
        if (jump < 0) {
            jump = tree.size();
            current.jump_left = jump;
            tree.resize(jump + 1);

            const auto& current = tree[location]; // rebind due to possible allocation.
            auto& latest = tree.back();
            latest.p_left = current.p_left;
            latest.p_right = current.p_mid;
            latest.p_mid = latest.p_left + (latest.p_right - latest.p_left) / 2;
        }
    } else {
        jump = current.jump_right;
        if (jump < 0) {
            jump = tree.size();
            current.jump_right = jump;
            tree.resize(jump + 1);

            const auto& current = tree[location]; // rebind due to possible allocation.
            auto& latest = tree.back();
            latest.p_left = current.p_mid;
            latest.p_right = current.p_right;
            latest.p_mid = latest.p_left + (latest.p_right - latest.p_left) / 2;
        }
    }

    build(p, jump, tree);
    return;
}

template<typename P>
int find_boundary(int location, const std::vector<Node>& tree, const P* pvalues, const int* accumulated, P threshold) {
    const auto& current = tree[location];
    auto curp = pvalues[current.p_left];
    int cumulative = accumulated[current.p_left];

    // Quitting if it doesn't work at its most optimistic, i.e., assigning
    // all points in the interval to the lowest p-value.
    if (curp > threshold * cumulative) {
        return 0;
    }

    // Returning the cumulative number of discoveries if we're at a terminus.
    // This terminus is the boundary of significance for the given FDR threshold;
    // anything lower will end up getting this adjusted p-value or better due
    // to the cumulative minimum in the BH method.
    if (current.p_right - current.p_left == 1) {
        return cumulative;
    }

    // Checking if the boundary lies on the right side first; if it's
    // successful, we don't bother searching the left, because the right's
    // discoveries will be a superset of those on the left.
    if (current.jump_right >= 0) {
        auto right_hits = find_boundary(current.jump_right, tree, pvalues, accumulated, threshold);
        if (right_hits) {
            return right_hits;
        }
    }

    // Otherwise, searching the left side for the boundary.
    if (current.jump_left >= 0) {
        return find_boundary(current.jump_left, tree, pvalues, accumulated, threshold);
    }

    return 0;
}

template<class Iterator>
std::vector<int> order(Iterator start, size_t n) {
    std::vector<int> temp;
    temp.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        temp.push_back(i);
    }

    std::sort(temp.begin(), temp.end(), [&](int left, int right) -> bool {
        return *(start + left) < *(start + right);
    });

    return temp;
}

template<typename P, typename C>
C find_optimal_filter(size_t n, const P* pvalues, const C* covariates, P fdr_threshold) {
    if (!n) {
        return 0;
    }

    // Finding all unique p-values.
    std::vector<int> prank(n);
    std::vector<P> uniq_p;
    std::vector<int> uniq_count;
    {
        auto porder = order(pvalues, n);
        uniq_p.reserve(n);
        uniq_p.push_back(*(pvalues + porder[0]));

        uniq_count.reserve(n); 
        uniq_count.push_back(1);

        int rank = 0;
        for (size_t i = 1; i < n; ++i) {
            const auto& current = *(pvalues + porder[i]);
            if (current != uniq_p.back()) {
                uniq_p.push_back(current);
                uniq_count.push_back(1);
                ++rank;
            } else {
                ++uniq_count.back();
            }
            prank[porder[i]] = rank;
        }
    }
    
    // Iteratively build tree and query for the FDR threshold at increasing
    // covariate thresholds.
    std::vector<Node> tree(1);
    tree.reserve(n * 2);
    tree.back().p_left = 0;
    tree.back().p_right = uniq_p.size();
    tree.back().p_mid = uniq_p.size() / 2;

    auto corder = order(covariates, n);
    size_t i = 0;
    int max_hits = -1;
    C cov_threshold = 0;

    while (i < n) {
        auto current = covariates[corder[i]];
        do {
            build(pvalues[corder[i]], 0, tree);
            ++i;
        } while (i < n && covariates[corder[i]] == current);

        // Dividing by 'i' to get the Bonferroni adjusted threshold, which is
        // then internally multiplied by the rank to get the BH threshold.
        auto hits = find_boundary(0, tree, uniq_p.data(), uniq_count.data(), fdr_threshold / i);
        if (hits > max_hits) {
            max_hits = hits;
            cov_threshold = current;
        }
    }

    return cov_threshold;
}

}

#endif
