#ifndef OIFF_H
#define OIFF_H

#include <vector>
#include <algorithm>

namespace oiff {

struct Node {
    int p_left, p_mid, p_right;
    int jump_left = -1, jump_right = -1;
    int count = 0;
};

void build(int p, int location, std::vector<Node>& tree) {
    auto& current = tree[location];
    ++current.count;
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
std::pair<int, int> find_boundary(int location, const std::vector<Node>& tree, const P* pvalues, int accumulated, P threshold) {
    const auto& current = tree[location];
    auto curp = pvalues[current.p_left];
    int accumulated2 = accumulated + current.count;

    // Quitting if it doesn't work at its most optimistic, i.e., assigning
    // all points in the interval to the lowest p-value.
    if (curp > threshold * accumulated2) {
        return std::make_pair(0, accumulated2);
    }

    // Returning the cumulative number of discoveries if we're at a terminus.
    // This terminus is the boundary of significance for the given FDR threshold;
    // anything lower will end up getting this adjusted p-value or better due
    // to the cumulative minimum in the BH method.
    if (current.p_right - current.p_left == 1) {
        return std::make_pair(accumulated2, accumulated2);
    }

    // Searching the left side for the boundary. We need to do this first
    // in order to get the cumulative sum before searching the right.
    int best = 0;
    if (current.jump_left >= 0) {
        auto left = find_boundary(current.jump_left, tree, pvalues, accumulated, threshold);
        best = left.first;
        accumulated = left.second;
    }

    // Now checking the right side.
    if (current.jump_right >= 0) {
        auto right = find_boundary(current.jump_right, tree, pvalues, accumulated, threshold);
        if (right.first) {
            best = right.first;
        }
    }

    return std::make_pair(best, accumulated2);
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
std::pair<C, int> find_optimal_filter(size_t n, const P* pvalues, const C* covariates, P fdr_threshold) {
    if (!n) {
        return std::pair<C, int>(0, 0);
    }

    // Finding all unique p-values.
    std::vector<int> prank(n, -1);
    std::vector<P> uniq_p;
    {
        auto porder = order(pvalues, n);
        uniq_p.reserve(n);

        int rank = 0;
        auto first = *(pvalues + porder[0]);
        if (first <= fdr_threshold) {
            uniq_p.push_back(first);
            prank[porder[0]] = rank;
        }

        for (size_t i = 1; i < n; ++i) {
            const auto& current = *(pvalues + porder[i]);
            if (current > fdr_threshold) {
                break;
            }
            if (current != uniq_p.back()) {
                uniq_p.push_back(current);
                ++rank;
            }
            prank[porder[i]] = rank;
        }
    }
    
    // Iteratively build tree and query for the FDR threshold at increasing
    // covariate thresholds.
    std::vector<Node> tree(1);
    tree.reserve(uniq_p.size() * 2);
    tree.back().p_left = 0;
    tree.back().p_right = uniq_p.size();
    tree.back().p_mid = uniq_p.size() / 2;

    auto corder = order(covariates, n);
    size_t i = 0;
    int max_hits = -1;
    C cov_threshold = 0;

    while (i < n) {
        auto curcov = covariates[corder[i]];
        do {
            auto curprank = prank[corder[i]];
            if (curprank >= 0) {
                build(curprank, 0, tree);
            }
            ++i;
        } while (i < n && covariates[corder[i]] == curcov);

        // Dividing by 'i' to get the Bonferroni adjusted threshold, which is
        // then internally multiplied by the rank to get the BH threshold.
        auto hits = find_boundary(0, tree, uniq_p.data(), 0, fdr_threshold / i);

        // Using >= so that we favor a more relaxed filter, everything else being equal.
        if (hits.first >= max_hits) {
            max_hits = hits.first;
            cov_threshold = curcov;
        }
    }

    return std::pair<C, int>(cov_threshold, max_hits);
}

}

#endif
