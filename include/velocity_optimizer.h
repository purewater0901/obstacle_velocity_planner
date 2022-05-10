#ifndef OBSTACLE_VELOCITY_PLANNER_VELOCITY_OPTIMIZER_H
#define OBSTACLE_VELOCITY_PLANNER_VELOCITY_OPTIMIZER_H

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Core>

#include "SBoundary.h"
#include "osqp_interface.h"

class VelocityOptimizer
{
public:
    struct OptimizationData
    {
        size_t N;
        double dt;
        double s0;
        double v0;
        double a0;
        double v_max;
        double a_max;
        double a_min;
        double limit_a_min;
        double j_max;
        double j_min;
        double lim_j_max;
        double lim_j_min;
        double t_dangerous;
        double t_idling;
        SBoundaries s_boundary;
    };

    struct OptimizationResult
    {
        std::vector<double> t;
        std::vector<double> s;
        std::vector<double> v;
        std::vector<double> a;
        std::vector<double> j;
    };

    VelocityOptimizer(const double max_s_weight, const double max_v_weight,
                      const double over_s_safety_weight, const double over_s_ideal_weight,
                      const double over_v_weight,
                      const double over_a_weight, const double over_j_weight);

    OptimizationResult optimize(const OptimizationData & data);

private:
    double max_s_weight_;
    double max_v_weight_;
    double over_s_safety_weight_;
    double over_s_ideal_weight_;
    double over_v_weight_;
    double over_a_weight_;
    double over_j_weight_;

    osqp::OSQPInterface qp_solver_;
};

#endif //OBSTACLE_VELOCITY_PLANNER_VELOCITY_OPTIMIZER_H
