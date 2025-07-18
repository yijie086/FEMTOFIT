#pragma once
#include "DataHandler.h"

struct CFFParams {
    double H;
    double E;
    double Ht;
    double Et;
};

class ModelCalculator {
public:
    void SetParametrization(const std::vector<double>& params);
    void PredictObservables(const DataPoint& dp, double& bsa_model, double& xs_model) const;

private:
    std::vector<double> params; // parameters for functional form
    CFFParams EvaluateCFF(double Q2, double xb, double t) const;
};
