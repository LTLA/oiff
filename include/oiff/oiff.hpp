#ifndef OIFF_H
#define OIFF_H

#include <vector>
#include <algorithm>
#include <random>
#include <thread>
#include <cstdint>
#include <cmath>

/**
 * @file oiff.hpp
 * @brief Optimizing the independent filter for FDR control.
 */

/**
 * @namespace oiff
 * @brief Optimizing the independent filter for FDR control.
 */
namespace oiff {

/**
 * @cond
 */
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

/* In this recursion, the return value is a pair containing:
 * 
 * - The largest number of discoveries at the specified FDR threshold
 *   within the interval represented by the current node.
 * - The number of hypotheses with p-values at or below the largest
 *   p-value in the interval represented by the current node.
 *
 * The first value is our desired result, but the second value needs
 * to be returned to allow calculation of a possibly larger number of 
 * discoveries at the parent node of the current node.
 */
template<typename P>
std::pair<int, int> find_boundary(int location, const std::vector<Node>& tree, const P* pvalues, int accumulated, P threshold) {
    const auto& current = tree[location];
    auto lowestp = pvalues[current.p_left];
    int accumulated2 = accumulated + current.count;

    // Quitting if it doesn't work at its most optimistic, i.e., assigning
    // all points in the interval to the lowest p-value.
    if (lowestp > threshold * accumulated2) {
        return std::make_pair(0, accumulated2);
    }

    // Returning the cumulative number of discoveries if we're at a terminal
    // node. This is potentially the boundary of significance for the given FDR
    // threshold; anything with lower p-values will end up getting this
    // adjusted p-value (or better) due to the cumulative minimum in the BH
    // method. Of course, we still need to check with the parent to see if
    // there might be more discoveries at an even larger p-value threshold.
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

    // Now checking the right side. If it gives us more discoveries, we switch
    // 'best' to its number of discoveries. No need to update 'accumulated' as
    // this should end up being equal to 'accumulated2' anyway.
    if (current.jump_right >= 0) {
        auto right = find_boundary(current.jump_right, tree, pvalues, accumulated, threshold);
        if (right.first) {
            best = right.first;
        }
    }

    return std::make_pair(best, accumulated2);
}

template<class Iterator>
std::vector<int> order(Iterator start, size_t n, bool reverse = false) {
    std::vector<int> temp;
    temp.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        temp.push_back(i);
    }

    if (reverse) {
        std::sort(temp.begin(), temp.end(), [&](int left, int right) -> bool {
            return *(start + left) > *(start + right);
        });
    } else {
        std::sort(temp.begin(), temp.end(), [&](int left, int right) -> bool {
            return *(start + left) < *(start + right);
        });
    }

    return temp;
}

template<class Engine>
double quick_uniform01(Engine& eng) {
    // Don't need to worry about checks for the edge case of returning 1
    // due to numerical imprecision after division; this won't break the 
    // use of this function in sample().
    return static_cast<double>(eng() - Engine::min()) / (static_cast<double>(Engine::max() - Engine::min()) + 1.0);
}

template<class Engine> 
void sample(size_t n, size_t s, std::vector<uint8_t>& chosen, Engine& eng) {
    std::fill(chosen.begin(), chosen.end(), static_cast<uint8_t>(0));
    for (size_t i = 0; i < n && s; ++i) {
        const double threshold = static_cast<double>(s)/(n - i);
        if (threshold >= 1 || quick_uniform01(eng) <= threshold) {
            chosen[i] = 1;
            --s;
        }
    }
}

template<typename T, typename Chosen>
void slice_vector(const std::vector<T>& input, const std::vector<Chosen>& keep, std::vector<T>& output) {
    output.clear();
    for (auto i : input) {
        if (keep[i]) {
            output.push_back(i);
        }
    }
}
/**
 * @endcond
 */

/**
 * @brief Optimize the filter threshold to maximize discoveries.
 *
 * This class contains the core functionality in the **oiff** library.
 * Methods are provided to identify the filter threshold that maximizes the number of discoveries at a nominal FDR threshold.
 * This can also be done via rounds of subsampling to mitigate loss of FDR control from aggressive inclusion of discoveries.
 */
struct OptimizeFilter {
    /**
     * The FDR threshold to be used.
     */
    double fdr_threshold = 0.05;

    /**
     * The number of threads to use during subsampling iterations.
     * Only used in `run_subsample()`.
     */
    int num_threads = 1;

    /**
     * The number of subsampling iterations to perform.
     * Only used in `run_subsample()`.
     */
    int num_iterations = 100;

    /**
     * The proportion of the hypotheses to subsample in each iteration.
     * Only used in `run_subsample()`.
     */
    double subsample_proportion = 0.1;

    /**
     * The random sed to use during subsampling.
     * Only used in `run_subsample()`.
     */
    uint64_t random_seed = 42;

    /**
     * Whether the filtering should retain larger values of the filter statistic.
     * By default, the filtering involves retaining values at or below the returned filter threshold.
     */
    bool retain_larger = false;

public:
    /**
     * @tparam F Arithmetic type for the filter statistic.
     * @brief Chosen filter thresholds that maximize the number of discoveries.
     */
    template<typename F>
    struct Thresholds {
        /**
         * First (i.e., most stringent) filter threshold that achieves the maximum number of discoveries.
         * If `OptimizeFilter::retain_larger` is false, this is the smallest filter threshold;
         * otherwise it is the largest.
         */
        F first = 0;

        /**
         * Last (i.e., most relaxed) filter threshold that achieves the maximum number of discoveries.
         * If `OptimizeFilter::retain_larger` is false, this is the largest filter threshold;
         * otherwise it is the smallest.
         */
        F last = 0;

        /**
         * An intermediate filter threshold that achieves the maximum number of discoveries.
         * This is usually somewhere between `first` and `last`,
         * and is close to the median of all filter thresholds that achieve maximum discoveries.
         * (It is exactly the median if the number of candidate thresholds is odd;
         * otherwise it is the closest threshold that is more relaxed than the median.) 
         */
        F middle = 0;

        /**
         * Maximum number of discoveries.
         */
        int number = 0;
    };

private:
    template<typename P, typename C>
    static Thresholds<C> find_optimal_filter(
        size_t n_full, 
        const P* pvalues, 
        const C* covariates, 
        const std::vector<int>& porder, 
        const std::vector<int>& corder, 
        P fdr_threshold, 
        std::vector<std::pair<C, int> >& candidates) 
    {
        // Without subsetting, n_working and n_full are the same. However, we want
        // to easily re-use the same code after subsetting, where 'porder' and
        // 'corder' contains a subset of the indices, while 'pvalues' and
        // 'covariates' still contain the full dataset to avoid a copy... Hence the
        // need to distinguish between n_full and n_working in the subsetted case.
        size_t n_working = porder.size();
        if (!n_working) {
            return Thresholds<C>();
        }

        // Finding all unique p-values and their relative ranks. Note that 
        // 'prank' is still full-length to enable easy indexing from 'corder'.
        std::vector<int> prank(n_full, -1);
        std::vector<P> uniq_p;
        {
            uniq_p.reserve(n_working);

            int rank = 0;
            auto first = *(pvalues + porder[0]);
            if (first <= fdr_threshold) {
                uniq_p.push_back(first);
                prank[porder[0]] = rank;
            }

            for (size_t i = 1; i < n_working; ++i) {
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
        
        // Iteratively build the tree and query for the FDR threshold. Each
        // iteration is done by increasing the covariate threshold and adding
        // more observations with that covariate value to the tree.
        std::vector<Node> tree(1);
        tree.reserve(uniq_p.size() * 2);
        tree.back().p_left = 0;
        tree.back().p_right = uniq_p.size();
        tree.back().p_mid = uniq_p.size() / 2;

        Thresholds<C> output;
        int& max_hits = output.number;
        candidates.clear();
        max_hits = -1;
        size_t i = 0;

        while (i < n_working) {
            size_t start = i;
            auto curcov = covariates[corder[i]];
            do {
                auto curprank = prank[corder[i]];
                if (curprank >= 0) {
                    build(curprank, 0, tree);
                }
                ++i;
            } while (i < n_working && covariates[corder[i]] == curcov);

            // Dividing by 'i' to get the Bonferroni adjusted threshold, which is
            // then internally multiplied by the rank to get the BH threshold.
            auto hits = find_boundary(0, tree, uniq_p.data(), 0, fdr_threshold / i);

            if (hits.first >= max_hits) {
                if (hits.first > max_hits) {
                    candidates.clear();
                    max_hits = hits.first;
                }
                candidates.emplace_back(curcov, i - start);
            }
        }

        if (candidates.size()) {
            output.first = candidates.front().first;
            output.last = candidates.back().first;

            // Computing the weighted median here.
            int total = 0; 
            for (const auto& x : candidates) {
                total += x.second;
            }
            int target = total / 2;

            total = 0;
            output.middle = output.first;
            for (const auto& x : candidates) {
                total += x.second;
                if (total > target) {
                    output.middle = x.first;
                    break;
                }
            }
        }

        return output;
    }

public:
    /**
     * @tparam P Floating-point type for the p-values.
     * @tparam F Arithmetic type for the filter statistic.
     *
     * @param n Number of hypotheses.
     * @param[in] pvalues Pointer to the start of an array of length `n` containing p-values for each hypothesis.
     * @param[in] filters Pointer to the start of an array of length `n` containing the filter statistic for each hypothesis.
     * 
     * @return A `Thresholds` object containing information about the chosen filter thresholds.
     */
    template<typename P, typename F>
    Thresholds<F> run(size_t n, const P* pvalues, const F* filters) const {
        auto porder = order(pvalues, n);
        auto corder = order(filters, n, retain_larger);
        std::vector<std::pair<F, int> > candidates;
        return find_optimal_filter(n, pvalues, filters, porder, corder, fdr_threshold, candidates);
    }

    /**
     * @tparam P Floating-point type for the p-values.
     * @tparam F Arithmetic type for the filter statistic.
     *
     * @param n Number of hypotheses.
     * @param[in] pvalues Pointer to the start of an array of length `n` containing p-values for each hypothesis.
     * @param[in] filters Pointer to the start of an array of length `n` containing the filter statistic for each hypothesis.
     * 
     * @return Vector of `Thresholds` objects, where each entry represents the results of one subsampling iteration.
     * Each object contains thresholds for a random subsample of hypotheses.
     */
    template<typename P, typename F>
    std::vector<Thresholds<F> > run_subsample(size_t n, const P* pvalues, const F* filters) const {
        auto porder = order(pvalues, n);
        auto corder = order(filters, n, retain_larger);

        std::mt19937_64 rng(random_seed);
        size_t keep = std::ceil(n * subsample_proportion);
        std::vector<Thresholds<F> > output(num_iterations);

        auto executor = [&](
            int iteration, 
            const std::vector<uint8_t>& s, 
            std::vector<int>& p, 
            std::vector<int>& c, 
            std::vector<std::pair<F, int> >& buffer) 
        -> void {
            slice_vector(porder, s, p);
            slice_vector(corder, s, c);
            output[iteration] = find_optimal_filter(n, pvalues, filters, p, c, fdr_threshold, buffer);
        };

        if (num_threads == 1) {
            std::vector<uint8_t> selected(n);
            std::vector<int> porder2, corder2;
            std::vector<std::pair<F, int> > candidates;
            for (int it = 0; it < num_iterations; ++it) {
                sample(n, keep, selected, rng);
                executor(it, selected, porder2, corder2, candidates);
            }

        } else {
            std::vector<std::vector<uint8_t> > selected(num_threads, std::vector<uint8_t>(n));
            std::vector<std::vector<int> > porder2(num_threads), corder2(num_threads);
            std::vector<std::vector<std::pair<F, int> > > candidates(num_threads);
            std::vector<std::thread> threads(num_threads);

            for (int it = 0; it < num_iterations; ++it) {
                size_t thread_id = it % num_threads;
                if (it >= num_threads) {
                    threads[thread_id].join();
                }
                sample(n, keep, selected[thread_id], rng); // RNG'ing is done in serial for simplicity.

                threads[thread_id] = std::thread(
                    [&](int iteration, size_t tid) -> void { 
                        executor(iteration, selected[tid], porder2[tid], corder2[tid], candidates[tid]); 
                    },
                    it, thread_id
                );
            }

            // Joining any threads that are still running.
            for (int t = 0; t < std::min(num_iterations, num_threads); ++t) {
                threads[t].join();
            }
        }

        return output;
    }
};

}

#endif
