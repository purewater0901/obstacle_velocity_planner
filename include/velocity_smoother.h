#ifndef OBSTACLE_VELOCITY_PLANNER_VELOCITY_SMOOTHER_H
#define OBSTACLE_VELOCITY_PLANNER_VELOCITY_SMOOTHER_H

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Core>

#include "osqp_interface.h"

class VelocitySmoother
{
public:
    struct OptimizationData
    {
        std::vector<double> s;
        std::vector<double> v_max;
        double v0;
        double a0;
        double a_max;
        double a_min;
        double j_max;
        double j_min;
    };

    struct OptimizationResult
    {
        std::vector<double> t;
        std::vector<double> s;
        std::vector<double> v;
        std::vector<double> a;
        std::vector<double> j;
    };

    VelocitySmoother(const double over_v_weight, const double over_a_weight, const double over_j_weight, const double jerk_weight);

    OptimizationResult filterVelocity(const OptimizationData & data);

    std::vector<double> forwardJerkFilter(const OptimizationData & data, const double a_start);
    std::vector<double> backwardJerkFilter(const OptimizationData & data, const double a_stop);

private:
    double over_v_weight_;
    double over_a_weight_;
    double over_j_weight_;
    double jerk_weight_;

    osqp::OSQPInterface qp_solver_;

};

#endif //OBSTACLE_VELOCITY_PLANNER_VELOCITY_SMOOTHER_H
