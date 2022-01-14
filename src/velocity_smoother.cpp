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

VelocitySmoother::OptimizationResult VelocitySmoother::filterVelocity(const OptimizationData & data)
{

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
