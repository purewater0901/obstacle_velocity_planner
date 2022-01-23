#include "velocity_optimizer.h"

VelocityOptimizer::VelocityOptimizer(const double max_s_weight, const double max_v_weight,
                                     const double over_s_safety_weight, const double over_s_ideal_weight,
                                     const double over_v_weight,
                                     const double over_a_weight, const double over_j_weight)
                                     : max_s_weight_(max_s_weight),
                                       max_v_weight_(max_v_weight),
                                       over_s_safety_weight_(over_s_safety_weight),
                                       over_s_ideal_weight_(over_s_ideal_weight),
                                       over_v_weight_(over_v_weight),
                                       over_a_weight_(over_a_weight),
                                       over_j_weight_(over_j_weight)
{
    qp_solver_.updateMaxIter(200000);
    qp_solver_.updateRhoInterval(0);  // 0 means automatic
    qp_solver_.updateEpsRel(1.0e-4);  // def: 1.0e-4
    qp_solver_.updateEpsAbs(1.0e-8);  // def: 1.0e-4
    qp_solver_.updateVerbose(false);
}


VelocityOptimizer::OptimizationResult VelocityOptimizer::optimize(const OptimizationData & data)
{
    const size_t N = data.N;
    const double dt = data.dt;
    const double s0 = data.s0;
    const double v0 = data.v0;
    const double a0 = data.a0;
    const double v_max = std::max(data.v_max, 0.1);
    const double a_max = data.a_max;
    const double a_min = data.a_min;
    const double j_max = data.j_max;
    const double j_min = data.j_min;
    const auto s_safety = data.s_safety;
    const auto s_ideal = data.s_ideal;

    // Variables: s_i, v_i, a_i, j_i, over_s_safety_i, over_s_ideal_i, over_v_i, over_a_i, over_j_i
    const int IDX_S0 = 0;
    const int IDX_V0 = N;
    const int IDX_A0 = 2 * N;
    const int IDX_J0 = 3 * N;
    const int IDX_OVER_S_SAFETY0= 4 * N;
    const int IDX_OVER_S_IDEAL0= 5 * N;
    const int IDX_OVER_V0 = 6 * N;
    const int IDX_OVER_A0= 7 * N;
    const int IDX_OVER_J0= 8 * N;
    const int l_variables = 9 * N;
    const int l_constraints = 5 * N + 3 * (N - 1) + 3;

    // the matrix size depends on constraint numbers.
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(l_constraints, l_variables);
    std::vector<double> lower_bound(l_constraints, 0.0);
    std::vector<double> upper_bound(l_constraints, 0.0);

    // Object Variables
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(l_variables, l_variables);
    std::vector<double> q(l_variables, 0.0);

    // Object Function
    // min (s_i - s_ideal_i)^2 + (v_i - v_max)^2
    // + over_s_safety^2 + over_s_ideal^2 + over_v_ideal^2 + over_a_ideal^2
    // + over_j_ideal^2
    for (size_t i = 0; i < N; ++i) {
        q.at(IDX_S0 + i) += -s_ideal.at(i).max_s * max_s_weight_;
        q.at(IDX_V0 + i) += -v_max * max_v_weight_;
    }

    for (size_t i = 0; i < N; ++i) {
        P(IDX_OVER_S_SAFETY0 + i, IDX_OVER_S_SAFETY0 + i) += over_s_safety_weight_;
        P(IDX_OVER_S_IDEAL0 + i, IDX_OVER_S_IDEAL0 + i) += over_s_ideal_weight_;
        P(IDX_OVER_V0 + i, IDX_OVER_V0 + i) += over_v_weight_;
        P(IDX_OVER_A0 + i, IDX_OVER_A0 + i) += over_a_weight_;
        P(IDX_OVER_J0 + i, IDX_OVER_J0 + i) += over_j_weight_;

        P(IDX_S0+ i, IDX_S0+ i) += max_s_weight_;
        P(IDX_V0+ i, IDX_V0+ i) += max_v_weight_;
    }

    // Constraint
    size_t constr_idx = 0;

    // Safety Position Constraint: s_safety_min < s_i - over_s_safety_i < s_safety_max
    for (size_t i = 0; i < N; ++i, ++constr_idx) {
        A(constr_idx, IDX_S0 + i) = 1.0;  // s_i
        A(constr_idx, IDX_OVER_S_SAFETY0 + i) = -1.0;  // over_s_safety_i
        upper_bound.at(constr_idx) = s_safety.at(i).max_s;
        lower_bound.at(constr_idx) = s_safety.at(i).min_s;
    }

    // Ideal Position Constraint: s_ideal_min < s_i - over_s_ideal_i < s_ideal_max
    for (size_t i = 0; i < N; ++i, ++constr_idx) {
        A(constr_idx, IDX_S0 + i) = 1.0;  // s_i
        A(constr_idx, IDX_OVER_S_IDEAL0 + i) = -1.0;  // over_s_ideal_i
        upper_bound.at(constr_idx) = s_ideal.at(i).max_s;
        lower_bound.at(constr_idx) = s_ideal.at(i).min_s;
    }

    // Soft Velocity Constraint: 0 < v_i - over_v_i < v_max
    for (size_t i = 0; i < N; ++i, ++constr_idx) {
        A(constr_idx, IDX_V0 + i) = 1.0;       // v_i
        A(constr_idx, IDX_OVER_V0 + i) = -1.0;  // over_v_i
        upper_bound.at(constr_idx) = v_max;
        lower_bound.at(constr_idx) = 0.0;
    }

    // Soft Acceleration Constraint: a_min < a_i - over_a_i < a_max
    for (size_t i = 0; i < N; ++i, ++constr_idx) {
        A(constr_idx, IDX_A0 + i) = 1.0;       // a_i
        A(constr_idx, IDX_OVER_A0 + i) = -1.0;  // over_a_i
        upper_bound.at(constr_idx) = a_max;
        lower_bound.at(constr_idx) = a_min;
    }

    // Soft Jerk Constraint: j_min < j_i - over_j_i < j_max
    for (size_t i = 0; i < N; ++i, ++constr_idx) {
        A(constr_idx, IDX_J0 + i) = 1.0;       // j_i
        A(constr_idx, IDX_OVER_J0 + i) = -1.0;  // over_j_i
        upper_bound.at(constr_idx) = j_max;
        lower_bound.at(constr_idx) = j_min;
    }

    // Dynamic Constraint
    // s_i+1 = s_i + v_i * dt + 0.5 * a_i * dt^2 + 1/6 * j_i * dt^3
    for (size_t i = 0; i < N - 1; ++i, ++constr_idx) {
        A(constr_idx, IDX_S0 + i + 1) = 1.0;                    // s_i+1
        A(constr_idx, IDX_S0 + i) = -1.0;                       // -s_i
        A(constr_idx, IDX_V0 + i) = -dt;                        // -v_i*dt
        A(constr_idx, IDX_A0 + i) = -0.5 * dt * dt;             // -0.5 * a_i * dt^2
        A(constr_idx, IDX_J0 + i) = -1.0 / 6.0 * dt * dt * dt;  // -1.0/6.0 * j_i * dt^3
        upper_bound.at(constr_idx) = 0.0;
        lower_bound.at(constr_idx) = 0.0;
    }

    // v_i+1 = v_i + a_i * dt + 0.5 * j_i * dt^2
    for (size_t i = 0; i < N - 1; ++i, ++constr_idx) {
        A(constr_idx, IDX_V0 + i + 1) = 1.0;         // v_i+1
        A(constr_idx, IDX_V0 + i) = -1.0;            // -v_i
        A(constr_idx, IDX_A0 + i) = -dt;             // -a_i * dt
        A(constr_idx, IDX_J0 + i) = -0.5 * dt * dt;  // -0.5 * j_i * dt^2
        upper_bound.at(constr_idx) = 0.0;
        lower_bound.at(constr_idx) = 0.0;
    }

    // a_i+1 = a_i + j_i * dt
    for (size_t i = 0; i < N - 1; ++i, ++constr_idx) {
        A(constr_idx, IDX_A0 + i + 1) = 1.0;  // a_i+1
        A(constr_idx, IDX_A0 + i) = -1.0;     // -a_i
        A(constr_idx, IDX_J0 + i) = -dt;      // -j_i * dt
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
    const auto result = qp_solver_.optimize(P, A, q, lower_bound, upper_bound);
    const std::vector<double> optval = std::get<0>(result);

    const int status_val = std::get<3>(result);
    if (status_val != 1) {
        std::cerr << "optimization failed : " << qp_solver_.getStatusMessage().c_str() << std::endl;
    }

    std::vector<double> opt_time(N);
    std::vector<double> opt_pos(N);
    std::vector<double> opt_vel(N);
    std::vector<double> opt_acc(N);
    std::vector<double> opt_jerk(N);
    for (size_t i = 0; i < N; ++i) {
        opt_time.at(i) = i * dt;
        opt_pos.at(i) = optval.at(IDX_S0 + i);
        opt_vel.at(i) = optval.at(IDX_V0 + i);
        opt_acc.at(i) = optval.at(IDX_A0 + i);
        opt_jerk.at(i) = optval.at(IDX_J0 + i);
        std::cout << "opt_jerk[" << i << "]: " << opt_jerk.at(i) << std::endl;
    }

    OptimizationResult optimized_result;
    optimized_result.t = opt_time;
    optimized_result.s = opt_pos;
    optimized_result.v = opt_vel;
    optimized_result.a = opt_acc;
    optimized_result.j = opt_jerk;

    return optimized_result;
}
