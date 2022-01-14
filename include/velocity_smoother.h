#ifndef OBSTACLE_VELOCITY_PLANNER_VELOCITY_SMOOTHER_H
#define OBSTACLE_VELOCITY_PLANNER_VELOCITY_SMOOTHER_H

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Core>

#include "osqp_interface.h"

class VelocitySmoother
{
public:
    VelocitySmoother(const double over_v_weight, const double over_a_weight, const double over_j_weight, const double jerk_weight);

private:
    double over_v_weight_;
    double over_a_weight_;
    double over_j_weight_;
    double jerk_weight_;

    osqp::OSQPInterface qp_solver_;
};

#endif //OBSTACLE_VELOCITY_PLANNER_VELOCITY_SMOOTHER_H
