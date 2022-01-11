#ifndef OBSTACLE_VELOCITY_PLANNER_SBOUNDARY_H
#define OBSTACLE_VELOCITY_PLANNER_SBOUNDARY_H

#include <limits>
#include <vector>

class SBoundary
{
public:
    SBoundary(const double _max_s, const double _min_s) : max_s(_max_s), min_s(_min_s) {}
    SBoundary() : max_s(std::numeric_limits<double>::max()), min_s(0.0) {}

    double max_s = std::numeric_limits<double>::max();
    double min_s = 0.0;
};

using SBoundaries = std::vector<SBoundary>;

#endif //OBSTACLE_VELOCITY_PLANNER_SBOUNDARY_H
