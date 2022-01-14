#ifndef OBSTACLE_VELOCITY_PLANNER_LINEAR_INTERPOLATION_H
#define OBSTACLE_VELOCITY_PLANNER_LINEAR_INTERPOLATION_H

#include "interpolation/interpolation_utils.h"

#include <vector>

namespace interpolation
{
    double lerp(const double src_val, const double dst_val, const double ratio);

    std::vector<double> lerp(
            const std::vector<double> & base_keys, const std::vector<double> & base_values,
            const std::vector<double> & query_keys);
}  // namespace interpolation

#endif //OBSTACLE_VELOCITY_PLANNER_LINEAR_INTERPOLATION_H
