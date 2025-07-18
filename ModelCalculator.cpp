#include "ModelCalculator.h"
#include <cmath>

void ModelCalculator::SetParametrization(const std::vector<double>& p) {
    params = p;
}

CFFParams ModelCalculator::EvaluateCFF(double Q2, double xb, double t) const {
    CFFParams cff;
    // Toy polynomial model: linear in each variable
    cff.H  = params[0] + params[1]*Q2 + params[2]*xb + params[3]*t;
    cff.E  = params[4] + params[5]*Q2 + params[6]*xb + params[7]*t;
    cff.Ht = params[8] + params[9]*t;
    cff.Et = params[10] + params[11]*t;
    return cff;
}

void ModelCalculator::PredictObservables(const DataPoint& dp, double& bsa_model, double& xs_model) const {
    CFFParams cff = EvaluateCFF(dp.Q2, dp.xb, dp.t);
    bsa_model = cff.H * std::sin(dp.phi * M_PI / 180.0) + cff.E * std::cos(dp.phi * M_PI / 180.0);
    xs_model = std::abs(cff.H) + std::abs(cff.E) + cff.Ht * std::exp(-dp.t) + cff.Et * dp.t;
}
