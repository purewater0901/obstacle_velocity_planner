#include "velocity_smoother.h"

VelocitySmoother::VelocitySmoother(const double over_v_weight,
                                   const double over_a_weight,
                                   const double over_j_weight,
                                   const double jerk_weight)
                                   : over_v_weight_(over_v_weight),
                                     over_a_weight_(over_a_weight),
                                     over_j_weight_(over_j_weight),
                                     jerk_weight_(jerk_weight)
{
    qp_solver_.updateMaxIter(200000);
    qp_solver_.updateRhoInterval(0);  // 0 means automatic
    qp_solver_.updateEpsRel(1.0e-4);  // def: 1.0e-4
    qp_solver_.updateEpsAbs(1.0e-8);  // def: 1.0e-4
    qp_solver_.updateVerbose(false);
}

std::vector<double> VelocitySmoother::forwardJerkFilter(const OptimizationData & data, const double a_start)
{
    const double v0 = data.v0;
    const double a0 = data.a0;
    const double a_max = data.a_max;
    const double j_max = data.j_max;

    auto applyLimits = [&data, &a_start](double & v, double & a, size_t i) {
        double v_lim = data.v_max.at(i);
        static constexpr double ep = 1.0e-5;
        if (v > v_lim + ep) {
            v = v_lim;
            a = 0.0;

            if (v_lim < 1e-3 && i < data.v_max.size() - 1) {
                double next_v_lim = data.v_max.at(i + 1);
                if (next_v_lim >= 1e-3) {
                    a = a_start;  // start from stop velocity
                }
            }
        }

        if (v < 0.0) {
            v = a = 0.0;
        }
    };

    auto output_vel = data.v_max;

    double current_vel = v0;
    double current_acc = a0;
    applyLimits(current_vel, current_acc, 0);

    output_vel.front() = current_vel;
    for (size_t i = 1; i < data.v_max.size(); ++i) {
        const double ds = data.s.at(i) - data.s.at(i-1);
        const double max_dt = std::pow(6.0 * ds / j_max, 1.0 / 3.0);  // assuming v0 = a0 = 0.
        const double dt = std::min(ds / std::max(current_vel, 1.0e-6), max_dt);

        if (current_acc + j_max * dt >= a_max) {
            const double tmp_jerk = std::min((a_max - current_acc) / dt, j_max);
            current_vel = current_vel + current_acc * dt + 0.5 * tmp_jerk * dt * dt;
            current_acc = a_max;
        } else {
            current_vel = current_vel + current_acc * dt + 0.5 * j_max * dt * dt;
            current_acc = current_acc + j_max * dt;
        }
        applyLimits(current_vel, current_acc, i);
        output_vel.at(i) = current_vel;
    }
    return output_vel;
}

std::vector<double> VelocitySmoother::backwardJerkFilter(const OptimizationData & data, const double a_stop)
{
    auto input_data = data;
    std::reverse(input_data.v_max.begin(), input_data.v_max.end());
    auto filtered_vel = forwardJerkFilter(input_data, a_stop);
    std::reverse(filtered_vel.begin(), filtered_vel.end());
    return filtered_vel;
}

std::vector<double> VelocitySmoother::mergeFilteredTrajectory(const OptimizationData & data, const std::vector<double> & forward_filtered, const std::vector<double> & backward_filtered) const
{
    const double v0 = data.v0;
    const double a0 = data.a0;
    const double a_min = data.a_min;
    const double j_min = data.j_min;
    const auto s_vec = data.s;

    auto merged = forward_filtered;

    auto getVx = [](const std::vector<double> &v_lim, int i) {
        return v_lim.at(i);
    };

    size_t i = 0;

    if (getVx(backward_filtered, 0) < v0) {
        double current_vel = v0;
        double current_acc = a0;
        while (getVx(backward_filtered, i) < current_vel && current_vel <= getVx(forward_filtered, i) &&
               i < merged.size() - 1) {
            merged.at(i) = current_vel;

            const double ds = s_vec.at(i+1) - s_vec.at(i);
            const double max_dt =
                    std::pow(6.0 * ds / std::fabs(j_min), 1.0 / 3.0);  // assuming v0 = a0 = 0.
            const double dt = std::min(ds / std::max(current_vel, 1.0e-6), max_dt);

            if (current_acc + j_min * dt < a_min) {
                const double tmp_jerk = std::max((a_min - current_acc) / dt, j_min);
                current_vel = current_vel + current_acc * dt + 0.5 * tmp_jerk * dt * dt;
                current_acc = std::max(current_acc + tmp_jerk * dt, a_min);
            } else {
                current_vel = current_vel + current_acc * dt + 0.5 * j_min * dt * dt;
                current_acc = current_acc + j_min * dt;
            }
            ++i;
        }
    }

    // take smaller velocity point
    for (; i < merged.size(); ++i) {
        merged.at(i) = (getVx(forward_filtered, i) < getVx(backward_filtered, i))
                              ? forward_filtered.at(i)
                              : backward_filtered.at(i);
    }
    return merged;
}

VelocitySmoother::OptimizationResult VelocitySmoother::optimize(const OptimizationData & data, const double a_stop_decel)
{
    const double v0 = data.v0;
    const double a0 = data.a0;
    const double a_max = data.a_max;
    const double a_min = data.a_min;
    const double j_max = data.j_max;
    const double j_min = data.j_min;

    const size_t N = data.s.size();
    std::vector<double> interval_dist_arr(N-1);
    for(size_t i=0; i<N-1; ++i) {
        interval_dist_arr.at(i) = data.s.at(i+1) - data.s.at(i);
    }

    std::vector<double> v_max_arr = data.v_max;

    /*
     * x = [
     *      b[0], b[1], ..., b[N],               : 0~N
     *      a[0], a[1], .... a[N],               : N~2N
     *      delta[0], ..., delta[N],             : 2N~3N
     *      sigma[0], sigma[1], ...., sigma[N],  : 3N~4N
     *      gamma[0], gamma[1], ..., gamma[N]    : 4N~5N
     *     ]
     *
     * b[i]  : velocity^2
     * delta : 0 < b[i]-delta[i] < max_vel[i]*max_vel[i]
     * sigma : a_min < a[i] - sigma[i] < a_max
     * gamma : jerk_min < pseudo_jerk[i] * ref_vel[i] - gamma[i] < jerk_max
     */
    const uint32_t IDX_B0 = 0;
    const uint32_t IDX_A0 = N;
    const uint32_t IDX_DELTA0 = 2 * N;
    const uint32_t IDX_SIGMA0 = 3 * N;
    const uint32_t IDX_GAMMA0 = 4 * N;

    const uint32_t l_variables = 5 * N;
    const uint32_t l_constraints = 4 * N + 1;

    // the matrix size depends on constraint numbers.
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(l_constraints, l_variables);

    std::vector<double> lower_bound(l_constraints, 0.0);
    std::vector<double> upper_bound(l_constraints, 0.0);

    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(l_variables, l_variables);
    std::vector<double> q(l_variables, 0.0);

    /**************************************************************/
    /**************************************************************/
    /**************** design objective function *******************/
    /**************************************************************/
    /**************************************************************/

    // jerk: d(ai)/ds * v_ref -> minimize weight * ((a1 - a0) / ds * v_ref)^2 * ds
    constexpr double ZERO_VEL_THR_FOR_DT_CALC = 0.3;
    const double smooth_weight = jerk_weight_;
    for (size_t i = 0; i < N - 1; ++i) {
        const double ref_vel = v_max_arr.at(i);
        const double interval_dist = std::max(interval_dist_arr.at(i), 0.0001);
        const double w_x_ds_inv = (1.0 / interval_dist) * ref_vel;
        P(IDX_A0 + i, IDX_A0 + i) += smooth_weight * w_x_ds_inv * w_x_ds_inv * interval_dist;
        P(IDX_A0 + i, IDX_A0 + i + 1) -= smooth_weight * w_x_ds_inv * w_x_ds_inv * interval_dist;
        P(IDX_A0 + i + 1, IDX_A0 + i) -= smooth_weight * w_x_ds_inv * w_x_ds_inv * interval_dist;
        P(IDX_A0 + i + 1, IDX_A0 + i + 1) += smooth_weight * w_x_ds_inv * w_x_ds_inv * interval_dist;
    }

    for (size_t i = 0; i < N; ++i) {
        const double v_max = std::max(v_max_arr.at(i), 0.1);
        q.at(IDX_B0 + i) =
                -1.0 / (v_max * v_max);  // |v_max_i^2 - b_i|/v_max^2 -> minimize (-bi) * ds / v_max^2
        if (i < N - 1) {
            q.at(IDX_B0 + i) *= std::max(interval_dist_arr.at(i), 0.0001);
        }
        P(IDX_DELTA0 + i, IDX_DELTA0 + i) += over_v_weight_;  // over velocity cost
        P(IDX_SIGMA0 + i, IDX_SIGMA0 + i) += over_a_weight_;  // over acceleration cost
        P(IDX_GAMMA0 + i, IDX_GAMMA0 + i) += over_j_weight_;  // over jerk cost
    }

    /**************************************************************/
    /**************************************************************/
    /**************** design constraint matrix ********************/
    /**************************************************************/
    /**************************************************************/

    /*
    NOTE: The delta allows b to be negative. This is actually invalid because the definition is b=v^2.
    But mathematically, the strict b>0 constraint may make the problem infeasible, such as the case of
    v=0 & a<0. To avoid the infeasibility, we allow b<0. The negative b is dealt as b=0 when it is
    converted to v with sqrt. If the weight of delta^2 is large (the value of delta is very small),
    b is almost 0, and is not a big problem.
    */

    size_t constr_idx = 0;

    // Soft Constraint Velocity Limit: 0 < b - delta < v_max^2
    for (size_t i = 0; i < N; ++i, ++constr_idx) {
        A(constr_idx, IDX_B0 + i) = 1.0;       // b_i
        A(constr_idx, IDX_DELTA0 + i) = -1.0;  // -delta_i
        upper_bound[constr_idx] = v_max_arr.at(i) * v_max_arr.at(i);
        lower_bound[constr_idx] = 0.0;
    }

    // Soft Constraint Acceleration Limit: a_min < a - sigma < a_max
    for (size_t i = 0; i < N; ++i, ++constr_idx) {
        A(constr_idx, IDX_A0 + i) = 1.0;       // a_i
        A(constr_idx, IDX_SIGMA0 + i) = -1.0;  // -sigma_i

        constexpr double stop_vel = 1e-3;
        if (v_max_arr.at(i) < stop_vel) {
            // Stop Point
            upper_bound[constr_idx] = a_stop_decel;
            lower_bound[constr_idx] = a_stop_decel;
        } else {
            upper_bound[constr_idx] = a_max;
            lower_bound[constr_idx] = a_min;
        }
    }

    // Soft Constraint Jerk Limit: jerk_min < pseudo_jerk[i] * ref_vel[i] - gamma[i] < jerk_max
    // -> jerk_min * ds < (a[i+1] - a[i]) * ref_vel[i] - gamma[i] * ds < jerk_max * ds
    for (size_t i = 0; i < N - 1; ++i, ++constr_idx) {
        const double ref_vel = std::max(v_max_arr.at(i), ZERO_VEL_THR_FOR_DT_CALC);
        const double ds = interval_dist_arr.at(i);
        A(constr_idx, IDX_A0 + i) = -ref_vel;     // -a[i] * ref_vel
        A(constr_idx, IDX_A0 + i + 1) = ref_vel;  //  a[i+1] * ref_vel
        A(constr_idx, IDX_GAMMA0 + i) = -ds;      // -gamma[i] * ds
        upper_bound[constr_idx] = j_max * ds;     //  jerk_max * ds
        lower_bound[constr_idx] = j_min * ds;     //  jerk_min * ds
    }

    // b' = 2a ... (b(i+1) - b(i)) / ds = 2a(i)
    for (size_t i = 0; i < N - 1; ++i, ++constr_idx) {
        A(constr_idx, IDX_B0 + i) = -1.0;                            // b(i)
        A(constr_idx, IDX_B0 + i + 1) = 1.0;                         // b(i+1)
        A(constr_idx, IDX_A0 + i) = -2.0 * interval_dist_arr.at(i);  // a(i) * ds
        upper_bound[constr_idx] = 0.0;
        lower_bound[constr_idx] = 0.0;
    }

    // initial condition
    {
        A(constr_idx, IDX_B0) = 1.0;  // b0
        upper_bound[constr_idx] = v0 * v0;
        lower_bound[constr_idx] = v0 * v0;
        ++constr_idx;

        A(constr_idx, IDX_A0) = 1.0;  // a0
        upper_bound[constr_idx] = a0;
        lower_bound[constr_idx] = a0;
        ++constr_idx;
    }

    // execute optimization
    const auto result = qp_solver_.optimize(P, A, q, lower_bound, upper_bound);
    const std::vector<double> optval = std::get<0>(result);

    OptimizationResult opt_result;
    opt_result.t.resize(N);
    opt_result.s = data.s;
    opt_result.v.resize(N);
    opt_result.a.resize(N);
    opt_result.j.resize(N);

    size_t zero_vel_id = N-2;
    // get velocity & acceleration
    for (size_t i = 0; i < N-1; ++i) {
        const double b = optval.at(IDX_B0 + i);
        const double v = std::sqrt(std::max(b, 0.0));
        if(i > 10 && v < 1e-3) {
            zero_vel_id = i;
            break;
        }
        opt_result.v.at(i) = v;
        opt_result.a.at(i) = optval.at(IDX_A0 + i);
    }

    for (size_t i = zero_vel_id; i < N; ++i) {
        opt_result.v.at(i) = 0.0;
        opt_result.a.at(i) = a_stop_decel;
    }

    opt_result.t.front() = 0.0;
    opt_result.j.front() = 0.0;
    for(size_t i=0; i<zero_vel_id; ++i) {
        const double da = opt_result.a.at(i+1) - opt_result.a.at(i);
        const double ds = interval_dist_arr.at(i);
        const double dt = ds/opt_result.v.at(i);
        opt_result.t.at(i+1) = opt_result.t.at(i) + dt;
        opt_result.j.at(i+1) = da/std::max(dt, 1e-6);
    }
    for (size_t i = zero_vel_id; i < N-1; ++i) {
        opt_result.t.at(i+1) = opt_result.t.at(i);
        opt_result.j.at(i+1) = 0.0;
    }

    const int status_val = std::get<3>(result);
    if (status_val != 1) {
        std::cerr << "optimization failed : " << qp_solver_.getStatusMessage().c_str() << std::endl;
    }

    return opt_result;
}

VelocitySmoother::OptimizationResult VelocitySmoother::optimizeL2(const OptimizationData & data, const double a_stop_decel)
{
    const double v0 = data.v0;
    const double a0 = data.a0;
    const double a_max = data.a_max;
    const double a_min = data.a_min;

    const size_t N = data.s.size();
    std::vector<double> interval_dist_arr(N-1);
    for(size_t i=0; i<N-1; ++i) {
        interval_dist_arr.at(i) = data.s.at(i+1) - data.s.at(i);
    }

    std::vector<double> v_max_arr = data.v_max;

    /*
     * x = [b0, b1, ..., bN, |  a0, a1, ..., aN, | delta0, delta1, ..., deltaN, | sigma0, sigma1, ..., sigmaN | gamma0, gamma1, ...
     * gammaN] in R^{5N} b: velocity^2 a: acceleration delta: 0 < bi < v_max^2 + delta sigma: a_min <
     * ai - sigma < a_max
     */

    const uint32_t l_variables = 5 * N;
    const uint32_t l_constraints = 3 * N + 1;

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(
            l_constraints, l_variables);  // the matrix size depends on constraint numbers.

    std::vector<double> lower_bound(l_constraints, 0.0);
    std::vector<double> upper_bound(l_constraints, 0.0);

    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(l_variables, l_variables);
    std::vector<double> q(l_variables, 0.0);

    // design objective function
    for (unsigned int i = 0; i < N; ++i) {  // bi
        q[i] = -1.0;                          // |v_max^2 - b| -> minimize (-bi)
    }

    // pseudo jerk: d(ai)/ds -> minimize weight * (a1 - a0)^2
    for (unsigned int i = N; i < 2 * N - 1; ++i) {
        unsigned int j = i - N;
        const double w_x_ds_inv = jerk_weight_* (1.0 / std::max(interval_dist_arr.at(j), 0.0001));
        P(i, i) += w_x_ds_inv * w_x_ds_inv;
        P(i, i + 1) -= w_x_ds_inv * w_x_ds_inv;
        P(i + 1, i) -= w_x_ds_inv * w_x_ds_inv;
        P(i + 1, i + 1) += w_x_ds_inv * w_x_ds_inv;
    }

    for (unsigned int i = 2 * N; i < 3 * N; ++i) {  // over velocity cost
        P(i, i) += over_v_weight_;
    }

    for (unsigned int i = 3 * N; i < 4 * N; ++i) {  // over acceleration cost
        P(i, i) += over_a_weight_;
    }

    for (unsigned int i = 4 * N; i < 5 * N; ++i) {  // over position cost
        P(i, i) += 0.1;
    }

    /* design constraint matrix
    0 < b - delta < v_max^2
    NOTE: The delta allows b to be negative. This is actually invalid because the definition is b=v^2.
    But mathematically, the strict b>0 constraint may make the problem infeasible, such as the case of
    v=0 & a<0. To avoid the infeasibility, we allow b<0. The negative b is dealt as b=0 when it is
    converted to v with sqrt. If the weight of delta^2 is large (the value of delta is very small),
    b is almost 0, and is not a big problem.
    */
    for (unsigned int i = 0; i < N; ++i) {
        const int j = 2 * N + i;
        A(i, i) = 1.0;   // b_i
        A(i, j) = -1.0;  // -delta_i
        upper_bound[i] = v_max_arr[i] * v_max_arr[i];
        lower_bound[i] = 0.0;
    }

    // a_min < a - sigma < a_max
    for (unsigned int i = N; i < 2 * N; ++i) {
        const int j = 2 * N + i;
        A(i, i) = 1.0;   // a_i
        A(i, j) = -1.0;  // -sigma_i
        if (i != N && v_max_arr[i - N] < std::numeric_limits<double>::epsilon()) {
            upper_bound[i] = a_stop_decel;
            lower_bound[i] = a_stop_decel;
        } else {
            upper_bound[i] = a_max;
            lower_bound[i] = a_min;
        }
    }

    // b' = 2a ... (b(i+1) - b(i)) / ds = 2a(i) + gamma(i)
    for (unsigned int i = 2 * N; i < 3 * N - 1; ++i) {
        const unsigned int j = i - 2 * N;
        const double ds_inv = 1.0 / std::max(interval_dist_arr.at(j), 0.0001);
        A(i, j) = -ds_inv;     // b(i)
        A(i, j + 1) = ds_inv;  // b(i+1)
        A(i, j + N) = -2.0;    // a(i)
        A(i, j + 4*N) = 1.0;  // gamma(i)
        upper_bound[i] = 0.0;
        lower_bound[i] = 0.0;
    }

    // initial condition
    {
        const unsigned int i = 3 * N - 1;
        A(i, 0) = 1.0;  // b0
        upper_bound[i] = v0 * v0;
        lower_bound[i] = v0 * v0;

        A(i + 1, N) = 1.0;  // a0
        upper_bound[i + 1] = a0;
        lower_bound[i + 1] = a0;
    }

    const auto result = qp_solver_.optimize(P, A, q, lower_bound, upper_bound);

    // [b0, b1, ..., bN, |  a0, a1, ..., aN, |
    //  delta0, delta1, ..., deltaN, | sigma0, sigma1, ..., sigmaN]
    const std::vector<double> optval = std::get<0>(result);

    OptimizationResult opt_result;
    opt_result.t.resize(N);
    opt_result.s = data.s;
    opt_result.v.resize(N);
    opt_result.a.resize(N);
    opt_result.j.resize(N);

    size_t zero_vel_id = N-2;
    // get velocity & acceleration
    for (size_t i = 0; i < N-1; ++i) {
        const double b = optval.at(i);
        const double v = std::sqrt(std::max(b, 0.0));
        if(i > 10 && v < 1e-3) {
            zero_vel_id = i;
            break;
        }
        opt_result.v.at(i) = v;
        opt_result.a.at(i) = optval.at(N + i);
    }

    for (size_t i = zero_vel_id; i < N; ++i) {
        opt_result.v.at(i) = 0.0;
        opt_result.a.at(i) = a_stop_decel;
    }

    opt_result.t.front() = 0.0;
    opt_result.j.front() = 0.0;
    for(size_t i=0; i<zero_vel_id; ++i) {
        const double da = opt_result.a.at(i+1) - opt_result.a.at(i);
        const double ds = interval_dist_arr.at(i);
        const double dt = ds/opt_result.v.at(i);
        opt_result.t.at(i+1) = opt_result.t.at(i) + dt;
        opt_result.j.at(i+1) = da/std::max(dt, 1e-6);
    }
    for (size_t i = zero_vel_id; i < N-1; ++i) {
        opt_result.t.at(i+1) = opt_result.t.at(i);
        opt_result.j.at(i+1) = 0.0;
    }

    const int status_val = std::get<3>(result);
    if (status_val != 1) {
        std::cerr << "optimization failed : " << qp_solver_.getStatusMessage().c_str() << std::endl;
    }

    return opt_result;
}
