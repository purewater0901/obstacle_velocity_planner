#include <iostream>
#include <chrono>
#include <vector>
#include <eigen3/Eigen/Core>

#include "matplotlibcpp.h"
#include "osqp_interface.h"
#include "SBoundary.h"

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

    osqp::OSQPInterface qp_solver;
    qp_solver.updateMaxIter(20000);
    qp_solver.updateRhoInterval(0);  // 0 means automatic
    //qp_solver.updateEpsRel(1.0e-4);  // def: 1.0e-4
    //qp_solver.updateEpsAbs(1.0e-8);  // def: 1.0e-4
    qp_solver.updateVerbose(false);

    // Variables: s_i, v_i, a_i, j_i, sigma_i, gamma_i, delta_i
    const int IDX_S0 = 0;
    const int IDX_V0 = N;
    const int IDX_A0 = 2*N;
    const int IDX_J0 = 3*N;
    const int IDX_SIGMA0 = 4*N;
    const int IDX_GAMMA0 = 5*N;
    const int IDX_DELTA0 = 6*N;
    const int l_variables = 7*N;
    const int l_constraints = 4*N + 3 * (N-1) + 3;

    // Start
    const auto ts = std::chrono::system_clock::now();

    // the matrix size depends on constraint numbers.
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(l_constraints, l_variables);

    std::vector<double> lower_bound(l_constraints, 0.0);
    std::vector<double> upper_bound(l_constraints, 0.0);

    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(l_variables, l_variables);
    std::vector<double> q(l_variables, 0.0);

    // Object Function
    // max v_i  -> min - v_i * dt
    for(size_t i=0; i<N; ++i) {
        q.at(IDX_V0 + i) = -dt;
    }

    // min sigma_i^2 * dt + gamma_i^2 * dt + delta_i^2 * dt
    for(size_t i=0; i<N; ++i) {
        P(IDX_SIGMA0 + i, IDX_SIGMA0 + i) += 1000;
        P(IDX_GAMMA0 + i, IDX_GAMMA0 + i) += 1000;
        P(IDX_DELTA0 + i, IDX_DELTA0 + i) += 1000;
    }

    // Constraint
    size_t constr_idx = 0;

    // Position Constraint: s_min < s_i < s_max
    for (size_t i=0; i<N; ++i, ++constr_idx) {
        A(constr_idx, IDX_S0 + i) = 1.0; // s_i
        upper_bound.at(constr_idx) = s_boundaries.at(i).max_s;
        lower_bound.at(constr_idx) = s_boundaries.at(i).min_s;
    }

    // Soft Velocity Constraint: 0 < v_i - sigma_i < v_max
    for (size_t i=0; i<N; ++i, ++constr_idx) {
        A(constr_idx, IDX_V0 + i) = 1.0; // v_i
        A(constr_idx, IDX_SIGMA0 + i) = -1.0; // sigma_i
        upper_bound.at(constr_idx) = v_max;
        lower_bound.at(constr_idx) = 0.0;
    }

    // Soft Acceleration Constraint: a_min < a_i - gamma_i < a_max
    for (size_t i=0; i<N; ++i, ++constr_idx) {
        A(constr_idx, IDX_A0 + i) = 1.0; // a_i
        A(constr_idx, IDX_GAMMA0 + i) = -1.0; // gamma_i
        upper_bound.at(constr_idx) = a_max;
        lower_bound.at(constr_idx) = a_min;
    }

    // Soft Jerk Constraint: j_min < j_i - delta_i < j_max
    for (size_t i=0; i<N; ++i, ++constr_idx) {
        A(constr_idx, IDX_J0 + i) = 1.0; // j_i
        A(constr_idx, IDX_DELTA0+ i) = -1.0; // delta_i
        upper_bound.at(constr_idx) = j_max;
        lower_bound.at(constr_idx) = j_min;
    }

    // Dynamic Constraint
    // s_i+1 = s_i + v_i * dt + 0.5 * a_i * dt^2 + 1/6 * j_i * dt^3
    for(size_t i=0; i<N-1; ++i, ++constr_idx) {
        A(constr_idx, IDX_S0+i+1) = 1.0; // s_i+1
        A(constr_idx, IDX_S0 + i) = -1.0; // -s_i
        A(constr_idx, IDX_V0 + i) = -dt; // -v_i*dt
        A(constr_idx, IDX_A0 + i) = -0.5 * dt * dt; // -0.5 * a_i * dt^2
        A(constr_idx, IDX_J0 + i) = -1.0/6.0 * dt * dt * dt; // -1.0/6.0 * j_i * dt^3
        upper_bound.at(constr_idx) = 0.0;
        lower_bound.at(constr_idx) = 0.0;
    }

    // v_i+1 = v_i + a_i * dt + 0.5 * j_i * dt^2
    for(size_t i=0; i<N-1; ++i, ++constr_idx) {
        A(constr_idx, IDX_V0 + i+1) = 1.0; // v_i+1
        A(constr_idx, IDX_V0 + i) = -1.0; // -v_i
        A(constr_idx, IDX_A0 + i) = -dt; // -a_i * dt
        A(constr_idx, IDX_J0 + i) = -0.5 * dt *dt; // -0.5 * j_i * dt^2
        upper_bound.at(constr_idx) = 0.0;
        lower_bound.at(constr_idx) = 0.0;
    }

    // a_i+1 = a_i + j_i * dt
    for(size_t i=0; i<N-1; ++i, ++constr_idx) {
        A(constr_idx, IDX_A0 + i+1) = 1.0; // a_i+1
        A(constr_idx, IDX_A0 + i) = -1.0; // -a_i
        A(constr_idx, IDX_J0 + i) = -dt; // -j_i * dt
        upper_bound.at(constr_idx) = 0.0;
        lower_bound.at(constr_idx) = 0.0;
    }

    // initial condition
    {
        A(constr_idx, IDX_S0) = 1.0;  // s0
        upper_bound[constr_idx] = s0;
        lower_bound[constr_idx] = s0;
        ++constr_idx;

        A(constr_idx, IDX_V0) = 1.0;  // v0
        upper_bound[constr_idx] = v0;
        lower_bound[constr_idx] = v0;
        ++constr_idx;

        A(constr_idx, IDX_A0) = 1.0;  // a0
        upper_bound[constr_idx] = a0;
        lower_bound[constr_idx] = a0;
    }

    // execute optimization
    const auto result = qp_solver.optimize(P, A, q, lower_bound, upper_bound);
    const std::vector<double> optval = std::get<0>(result);

    const auto tf1 = std::chrono::system_clock::now();
    const double dt_ms1 =
            std::chrono::duration_cast<std::chrono::nanoseconds>(tf1 - ts).count() * 1.0e-6;
    std::cout << "optimization time = " << dt_ms1 <<  "[ms]" << std::endl;

    // get velocity & acceleration
    std::vector<double> opt_pos(N);
    std::vector<double> opt_vel(N);
    std::vector<double> opt_acc(N);
    std::vector<double> opt_jerk(N);
    for (size_t i = 0; i < N; ++i) {
        opt_pos.at(i) = optval.at(IDX_S0 + i);
        opt_vel.at(i) = optval.at(IDX_V0 + i);
        opt_acc.at(i) = optval.at(IDX_A0 + i);
        opt_jerk.at(i) = optval.at(IDX_J0 + i);
        std::cout << "s[" << i << "]: " << opt_pos.at(i)
        << " v[" << i << "]: " << opt_vel.at(i)
        << " a[" << i << "]: " << opt_acc.at(i)
        << " j[" << i << "]: " << opt_jerk.at(i) << std::endl;
    }

    const int status_val = std::get<3>(result);
    if (status_val != 1) {
        std::cout <<"optimization failed : " << qp_solver.getStatusMessage().c_str() << std::endl;
        return -1;
    }

    std::vector<double> opt_time(N, 0.0);
    for (size_t i=1; i<N; ++i) {
        opt_time.at(i) = opt_time.at(i-1) + dt;
    }
    std::vector<double> max_vels(N, v_max);
    matplotlibcpp::figure_size(1200, 700);
    // matplotlibcpp::named_plot("optimal_velocity", opt_pos, opt_vel);
    // matplotlibcpp::named_plot("maximum_velocity", opt_pos, max_vels);
    matplotlibcpp::named_plot("trajectory", opt_time, opt_pos);
    matplotlibcpp::named_plot("obstacle", opt_time, s_lim);
    matplotlibcpp::title("Result");
    matplotlibcpp::legend();
    matplotlibcpp::show();

    return 0;
}
