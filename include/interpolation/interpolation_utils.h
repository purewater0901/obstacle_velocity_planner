#ifndef OBSTACLE_VELOCITY_PLANNER_INTERPOLATION_UTILS_H
#define OBSTACLE_VELOCITY_PLANNER_INTERPOLATION_UTILS_H

#include <array>
#include <stdexcept>
#include <vector>
#include <string>

namespace interpolation_utils
{
    inline bool isIncreasing(const std::vector<double> & x)
    {
        if (x.empty()) {
            throw std::invalid_argument("Points is empty.");
        }

        for (size_t i = 0; i < x.size() - 1; ++i) {
            if (x.at(i) >= x.at(i + 1)) {
                return false;
            }
        }

        return true;
    }

    inline bool isNotDecreasing(const std::vector<double> & x)
    {
        if (x.empty()) {
            throw std::invalid_argument("Points is empty.");
        }

        for (size_t i = 0; i < x.size() - 1; ++i) {
            if (x.at(i) > x.at(i + 1)) {
                return false;
            }
        }

        return true;
    }

    inline void validateInput(
            const std::vector<double> & base_keys, const std::vector<double> & base_values,
            const std::vector<double> & query_keys)
    {
        // when vectors are empty
        if (base_keys.empty() || base_values.empty() || query_keys.empty()) {
            throw std::invalid_argument("Points is empty.");
        }

        // when size of vectors are less than 2
        if (base_keys.size() < 2 || base_values.size() < 2) {
            throw std::invalid_argument(
                    "The size of points is less than 2. base_keys.size() = " + std::to_string(base_keys.size()) +
                    ", base_values.size() = " + std::to_string(base_values.size()));
        }

        // when indices are not sorted
        if (!isIncreasing(base_keys) || !isNotDecreasing(query_keys)) {
            throw std::invalid_argument("Either base_keys or query_keys is not sorted.");
        }

        // when query_keys is out of base_keys (This function does not allow exterior division.)
        if (query_keys.front() < base_keys.front() || base_keys.back() < query_keys.back()) {
            throw std::invalid_argument("query_keys is out of base_keys");
        }

        // when sizes of indices and values are not same
        if (base_keys.size() != base_values.size()) {
            throw std::invalid_argument("The size of base_keys and base_values are not the same.");
        }
    }
}  // namespace interpolation_utils

#endif //OBSTACLE_VELOCITY_PLANNER_INTERPOLATION_UTILS_H
