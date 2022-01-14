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