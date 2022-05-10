#include <iostream>
#include <chrono>
#include <vector>
#include <memory>
#include <eigen3/Eigen/Core>

#include "matplotlibcpp.h"
#include "osqp_interface.h"
#include "SBoundary.h"
#include "velocity_optimizer.h"
#include "velocity_smoother.h"
#include "interpolation/linear_interpolation.h"

int main() {

    const int N = 50;
    const double dt = 0.2;
    const double ds = 0.1;
    const double v_max = 10.0;
    const double a_max = 1.0;
    const double a_min = -1.0;
    const double limit_a_min = -3.0;
    const double j_max = 1.0;
    const double j_min = -1.0;
    const double s0 = 0.0;
    const double v0 = 3.0;
    const double a0 = 0.0;
    const double t_dangerous = 1.0;
    const double t_idling = 2.0;
    const double max_s_weight = 100.0;
    const double max_v_weight = 1.0;
    const double over_s_safety_weight = 1000000.0;
    const double over_s_ideal_weight = 50.0;
    const double over_v_weight = 500000.0;
    const double over_a_weight = 5000.0;
    const double over_j_weight = 10000.0;

    std::shared_ptr<VelocityOptimizer> velocity_optimizer_ptr_;
    velocity_optimizer_ptr_ =
            std::make_shared<VelocityOptimizer>(max_s_weight, max_v_weight, over_s_safety_weight, over_s_ideal_weight,
                                                over_v_weight, over_a_weight, over_j_weight);

    SBoundaries s_boundary(N);
    const double s_start = 6.0;
    const double v_obj = 0.0;
    const double s_max = 100; //s0 + v_max * (N-1) * dt;
    for(size_t i=0; i<N; ++i) {
        s_boundary.at(i).max_s = s_max;
    }
    for(size_t i=10; i<20; ++i) {
        s_boundary.at(i).max_s = s_start + v_obj * (i - 10) * dt;
        s_boundary.at(i).is_object = true;
    }
    for(size_t i=30; i<N; ++i) {
        s_boundary.at(i).max_s = s_boundary.at(29).max_s + v_max * (i - 29) * dt;
        s_boundary.at(i).is_object = false;
    }

    std::vector<double> s_safety_bound(N);
    for(size_t i=0; i<N; ++i) {
        s_safety_bound.at(i) = s_boundary.at(i).max_s;
    }

    VelocityOptimizer::OptimizationData data;
    data.N = N;
    data.dt = dt;
    data.s0 = 0.0;
    data.v0 = v0;
    data.a0 = a0;
    data.v_max = v_max;
    data.a_max = a_max;
    data.a_min = a_min;
    data.limit_a_min = limit_a_min;
    data.j_max = j_max;
    data.j_min = j_min;
    data.t_dangerous = t_dangerous;
    data.t_idling = t_idling;
    data.s_boundary = s_boundary;

    // Velocity Optimizer for Obstacle Avoidance
    const auto optimization_start_time = std::chrono::system_clock::now();
    const auto optimized_result = velocity_optimizer_ptr_->optimize(data);
    const auto optimization_end_time = std::chrono::system_clock::now();
    const double calculation_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
            optimization_end_time - optimization_start_time)
                                            .count() * 1.0e-6;

    std::cout << "Optimization Time: " << calculation_time << "[ms]" << std::endl;

    // Transformation from t to s
    std::vector<double> opt_positions = {0.0};
    std::vector<double> opt_velocity = {data.v0};
    for (size_t i = 1; i < optimized_result.s.size(); ++i) {
        const double prev_s = opt_positions.back();
        const double current_s = std::max(optimized_result.s.at(i), 0.0);
        const double current_v = std::max(optimized_result.v.at(i), 0.0);
        if (prev_s >= current_s) {
            continue;
        }
        opt_positions.push_back(current_s);
        opt_velocity.push_back(current_v);
    }

    if (opt_positions.size() == 1) {
        std::cout << "[Velocity Optimizer Result]: Optimized Trajectory is too small" << std::endl;
        return -1;
    }

    std::cout << "[Velocity Optimizer Result]: Trajectory Length: " << opt_positions.back() << std::endl;
    std::vector<double> query_positions;
    for (double s = 0; s <= opt_positions.back(); s+=ds) {
        query_positions.push_back(s);
    }
    auto resampled_opt_velocity =
            interpolation::lerp(opt_positions, opt_velocity, query_positions);
    resampled_opt_velocity.back() = 0.0;

    /*
     * Velocity Smoother
     */
    const double over_v_weight_smoother = 10000.0;
    const double over_a_weight_smoother = 500.0;
    const double over_j_weight_smoother = 200.0;
    const double jerk_weight = 0.0;
    const auto velocity_smoother_ptr_ = std::make_shared<VelocitySmoother>(over_v_weight_smoother, over_a_weight_smoother, over_j_weight_smoother, jerk_weight);

    std::vector<double> tmp_max_vels(resampled_opt_velocity.size(), 8.0);
    tmp_max_vels.back() = 0.0;

    VelocitySmoother::OptimizationData smoother_data_forward;
    smoother_data_forward.s = query_positions;
    smoother_data_forward.v_max = resampled_opt_velocity;
    smoother_data_forward.v0 = v0;
    smoother_data_forward.a0 = a0;
    smoother_data_forward.a_max = a_max;
    smoother_data_forward.a_min = a_min;
    smoother_data_forward.j_max = j_max;
    smoother_data_forward.j_min = j_min;
    const double a_stop_accel = 0.0;
    const double a_stop_decel = 0.0;

    const auto forward_filtered_vel = velocity_smoother_ptr_->forwardJerkFilter(smoother_data_forward, a_stop_accel);

    VelocitySmoother::OptimizationData smoother_data_backward = smoother_data_forward;
    smoother_data_backward.v0 = 0.0;
    smoother_data_backward.a0 = 0.0;
    const auto backward_filtered_vel = velocity_smoother_ptr_->backwardJerkFilter(smoother_data_forward, a_stop_decel);

    const auto merged_filtered_vel = velocity_smoother_ptr_->mergeFilteredTrajectory(smoother_data_forward, forward_filtered_vel, backward_filtered_vel);

    VelocitySmoother::OptimizationData smoother_data = smoother_data_forward;
    smoother_data.v_max = merged_filtered_vel;
    const auto smoothed_result = velocity_smoother_ptr_->optimize(smoother_data, a_stop_decel);
    const auto l2_smoothed_result = velocity_smoother_ptr_->optimizeL2(smoother_data, a_stop_decel);

    std::cout << "Finish Optimization" << std::endl;

    // Visualization
    std::vector<double> max_vels(N, v_max);
    matplotlibcpp::figure_size(1500, 900);
    matplotlibcpp::subplot(2, 1, 1);
    matplotlibcpp::named_plot("trajectory", optimized_result.t, optimized_result.s);
    matplotlibcpp::named_plot("Safety Bound", optimized_result.t, s_safety_bound);
    matplotlibcpp::legend();
    matplotlibcpp::title("Position");
    matplotlibcpp::subplot(2, 1, 2);
    matplotlibcpp::named_plot("maximum_velocity", optimized_result.s, max_vels);
    matplotlibcpp::named_plot("optimal_velocity", query_positions, resampled_opt_velocity);
    //matplotlibcpp::named_plot("forward_velocity", query_positions, forward_filtered_vel);
    //matplotlibcpp::named_plot("backward_velocity", query_positions, backward_filtered_vel);
    matplotlibcpp::named_plot("merged_velocity", query_positions, merged_filtered_vel);
    matplotlibcpp::named_plot("smoothed_velocity", query_positions, smoothed_result.v);
    //matplotlibcpp::named_plot("smoothed_acceleration", query_positions, smoothed_result.a);
    matplotlibcpp::named_plot("smoothed_jerk", query_positions, smoothed_result.j);
    // matplotlibcpp::named_plot("l2_smoothed_velocity", query_positions, l2_smoothed_result.v);
    //matplotlibcpp::named_plot("l2_smoothed_acceleration", query_positions, l2_smoothed_result.a);
    //matplotlibcpp::named_plot("l2_smoothed_jerk", query_positions, l2_smoothed_result.j);
    matplotlibcpp::title("Velocity");
    matplotlibcpp::legend();
    matplotlibcpp::show();

    return 0;
}
