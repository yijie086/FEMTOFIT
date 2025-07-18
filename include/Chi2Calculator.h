#pragma once
#include "ModelCalculator.h"

class Chi2Calculator {
public:
    Chi2Calculator(const std::vector<DataPoint>& data);
    double operator()(const std::vector<double>& fit_params);

private:
    std::vector<DataPoint> dataPoints;
    ModelCalculator model;
};
