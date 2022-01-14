#include <iostream>
#include <chrono>
#include <vector>
#include <memory>
#include <eigen3/Eigen/Core>

#include "matplotlibcpp.h"
#include "osqp_interface.h"
#include "SBoundary.h"
#include "velocity_optimizer.h"

int main() {

    const int N = 25;
    const double dt = 0.3;
    const double v_max = 10.0;
    const double a_max = 1.0;
    const double a_min = -1.0;
    const double j_max = 1.0;
    const double j_min = -1.0;
    const double s0 = 0.0;
    const double v0 = 5.0;
    const double a0 = 0.0;
    const double over_v_weight = 1000.0;
    const double over_a_weight = 1000.0;
    const double over_j_weight = 1000.0;

    std::shared_ptr<VelocityOptimizer> velocity_optimizer_ptr_;
    velocity_optimizer_ptr_ =
            std::make_shared<VelocityOptimizer>(over_v_weight, over_a_weight, over_j_weight);

    SBoundaries s_boundaries(N);
    for(size_t i=0; i<7; ++i) {
        s_boundaries.at(i).max_s = 8.5;
    }
    for(size_t i=7; i<16; ++i) {
        s_boundaries.at(i).max_s = 15.0;
    }
    for(size_t i=16; i<N; ++i) {
        s_boundaries.at(i).max_s = 25.0;
    }

    std::vector<double> s_lim(N);
    for(size_t i=0; i<N; ++i) {
        s_lim.at(i) = s_boundaries.at(i).max_s;
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
    data.j_max = j_max;
    data.j_min = j_min;
    data.s_boundaries = s_boundaries;

    // Velocity Optimizer for Obstacle Avoidance
    const auto optimization_start_time = std::chrono::system_clock::now();
    const auto optimized_result = velocity_optimizer_ptr_->optimize(data);
    const auto optimization_end_time = std::chrono::system_clock::now();
    const double calculation_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
            optimization_end_time - optimization_start_time)
                                            .count() *
                                    1.0e-6;

    std::cout << "Optimization Time: " << calculation_time << "[ms]" << std::endl;

    std::vector<double> max_vels(N, v_max);
    matplotlibcpp::figure_size(1200, 700);
    // matplotlibcpp::named_plot("optimal_velocity", opt_pos, opt_vel);
    // matplotlibcpp::named_plot("maximum_velocity", opt_pos, max_vels);
    matplotlibcpp::named_plot("trajectory", optimized_result.t, optimized_result.s);
    matplotlibcpp::named_plot("obstacle", optimized_result.t, s_lim);
    matplotlibcpp::title("Result");
    matplotlibcpp::legend();
    matplotlibcpp::show();

    return 0;
}
