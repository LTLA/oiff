#ifndef OIFF_QUADTREE_HPP
#define OIFF_QUADTREE_HPP

#include <vector>

namespace oiff {

template<typename T>
struct QuadTreeNode {
    T corner_p;
    T corner_c;
    T width_p;
    T width_c;

    // Sacrifice some memory for speed by pre-computing the midpoints.
    T mid_p;
    T mid_c;

    int count = 0;

    int jump_ = -1;
    int jump_p = -1;
    int jump_c = -1;
    int jump_pc = -1;
};

template<typename T>
void recursive_build_tree(T p, T c, int location, std::vector<QuadTreeNode<T> >& tree) {
    auto& current = tree[location];
    ++current.count;
    if (current.width_p == 1 && current.width_c == 1) {
        return;
    }

    int jump;
    if (p < current.mid_p) {
        if (c < current.mid_c) {
            jump = current.jump_;
            if (jump < 0) {
                jump = tree.size();
                current.jump_ = jump;
                tree.resize(jump + 1);
                const auto& current = tree[location]; // rbind due to allocation.
                auto& latest = tree.back();
                latest.corner_p = current.corner_p;
                latest.corner_c = current.corner_c;
                latest.width_p = current.mid_p - current.corner_p;
                latest.width_c = current.mid_c - current.corner_c;
                latest.mid_p = latest.corner_p + latest.width_p / 2;
                latest.mid_c = latest.corner_c + latest.width_c / 2;
            }
        } else {
            jump = current.jump_c;
            if (jump < 0) {
                jump = tree.size();
                current.jump_c = jump;
                tree.resize(jump + 1);
                const auto& current = tree[location]; // rbind due to allocation.
                auto& latest = tree.back();
                latest.corner_p = current.corner_p;
                latest.corner_c = current.mid_c;
                latest.width_p = current.mid_p - current.corner_p;
                latest.width_c = current.corner_c + current.width_c - current.mid_c;
                latest.mid_p = latest.corner_p + latest.width_p / 2;
                latest.mid_c = latest.corner_c + latest.width_c / 2;
            }
        }
    } else {
        if (c < current.mid_c) {
            jump = current.jump_p;
            if (jump < 0) {
                jump = tree.size();
                current.jump_p = jump;
                tree.resize(jump + 1);
                const auto& current = tree[location]; // rbind due to allocation.
                auto& latest = tree.back();
                latest.corner_p = current.mid_p;
                latest.corner_c = current.corner_c;
                latest.width_p = current.corner_p + current.width_p - current.mid_p;
                latest.width_c = current.mid_c - current.corner_c;
                latest.mid_p = latest.corner_p + latest.width_p / 2;
                latest.mid_c = latest.corner_c + latest.width_c / 2;
            }
        } else {
            jump = current.jump_pc;
            if (jump < 0) {
                jump = tree.size();
                current.jump_pc = jump;
                tree.resize(jump + 1);
                const auto& current = tree[location]; // rbind due to allocation.
                auto& latest = tree.back();
                latest.corner_p = current.mid_p;
                latest.corner_c = current.mid_c;
                latest.width_p = current.corner_p + current.width_p - current.mid_p;
                latest.width_c = current.corner_c + current.width_c - current.mid_c;
                latest.mid_p = latest.corner_p + latest.width_p / 2;
                latest.mid_c = latest.corner_c + latest.width_c / 2;
            }
        }
    }

    recursive_build_tree(p, c, jump, tree);
    return;
}

template<typename T>
using QuadTree = std::vector<QuadTreeNode<T> >;

template<typename T>
QuadTree<T> build_quad_tree(const std::vector<T>& prank, const std::vector<T>& crank) {
    size_t nobs = prank.size();
    QuadTree<T> tree(1);
    tree.reserve(nobs * 2);

    // Filling in the starting point.
    T min_p = 1, max_p = 1, min_c = 1, max_c = 1;
    if (!prank.empty()) {
        min_p = prank[0];
        max_p = prank[0];
        min_c = crank[0];
        min_c = crank[0];
    }

    for (size_t i = 1; i < nobs; ++i) {
        if (min_p > prank[i]) {
            min_p = prank[i];
        } else if (max_p < prank[i]) {
            max_p = prank[i];
        }
    }

    for (size_t i = 1; i < nobs; ++i) {
        if (min_c > crank[i]) {
            min_c = crank[i];
        } else if (max_c < crank[i]) {
            max_c = crank[i];
        }
    }

    tree[0].corner_c = min_c;
    tree[0].corner_p = min_p;
    tree[0].width_c = max_c - min_c + 1; // +1 as we need to include the max value.
    tree[0].width_p = max_p - min_p + 1;
    tree[0].mid_c = tree[0].corner_c + tree[0].width_c / 2;
    tree[0].mid_p = tree[0].corner_p + tree[0].width_p / 2;

    // Filling out the tree.
    for (size_t i = 0; i < nobs; ++i) {
        recursive_build_tree(prank[i], crank[i], 0, tree);
    }

    return tree;
}

template<typename T>
void recursive_count_under(T p, T c, int location, const QuadTree<T>& tree, int& count) {
    const auto& current = tree[location];
    if (current.width_p == 1 && current.width_c == 1) {
        count += current.count;
        return;
    }

    bool low_c = c < current.mid_c;
    bool low_p = p < current.mid_p;

    if (low_c || low_p) {
        if (current.jump_ >= 0) {
            recursive_count_under(p, c, current.jump_, tree, count);
        }

        if (!low_p) {
            if (current.jump_p >= 0) {
                recursive_count_under(p, c, current.jump_p, tree, count);
            }
        }
 
        if (!low_c) {
            if (current.jump_c >= 0) {
                recursive_count_under(p, c, current.jump_c, tree, count);
            }
        }

        return;
    }

    bool range_c = c < current.corner_c + current.width_c;
    bool range_p = p < current.corner_p + current.width_p;

    if (range_c || range_p) {
        if (current.jump_ >= 0) {
            count += tree[current.jump_].count;
        }

        if (current.jump_c >= 0) {
            if (range_c) {
                recursive_count_under(p, c, current.jump_c, tree, count);
            } else {
                count += tree[current.jump_c].count;
            }
        }

        if (current.jump_p >= 0) {
            if (range_p) {
                recursive_count_under(p, c, current.jump_p, tree, count);
            } else {
                count += tree[current.jump_p].count;
            }
        }

        if (current.jump_pc >= 0) {
            recursive_count_under(p, c, current.jump_pc, tree, count);
        }

        return;
    }

    count += current.count;
    return;
}

template<typename T>
int count_under(T p, T c, const QuadTree<T>& tree) {
    int count = 0;
    recursive_count_under(p, c, 0, tree, count);
    return count;
}

}

#endif
