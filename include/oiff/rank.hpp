#ifndef OIFF_RANK_HPP
#define OIFF_RANK_HPP

#include <vector>
#include <algorithm>

namespace oiff {

template<typename T = int, class Iterator>
std::vector<T> rankify(Iterator start, size_t n) {
    std::vector<T> temp;
    temp.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        temp.push_back(i);
    }

    std::sort(temp.begin(), temp.end(), [&](T left, T right) -> bool {
        return *(start + left) < *(start + right);
    });

    std::vector<T> output(n);
    if (n) {
        T rank = 0;
        auto last = *(start + temp[0]);
        for (size_t i = 0; i < n; ++i) {
            const auto& current = *(start + temp[i]);
            if (current != last) {
                last = current;
                ++rank;
            }
            output[temp[i]] = rank;
        }
    }

    return output;
}

}

#endif
