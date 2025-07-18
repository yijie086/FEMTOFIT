#include "Chi2Calculator.h"
#include <cmath>

Chi2Calculator::Chi2Calculator(const std::vector<DataPoint>& data) : dataPoints(data) {}

double Chi2Calculator::operator()(const std::vector<double>& fit_params) {
    model.SetParametrization(fit_params);
    double chi2 = 0.0;
    for (const auto& dp : dataPoints) {
        double bsa_model, xs_model;
        model.PredictObservables(dp, bsa_model, xs_model);
        chi2 += std::pow((dp.BSA - bsa_model) / dp.BSA_err, 2);
        chi2 += std::pow((dp.XS - xs_model) / dp.XS_err, 2);
    }
    return chi2;
}
