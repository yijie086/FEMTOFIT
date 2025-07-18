#pragma once
#include "Chi2Calculator.h"

class Fitter {
public:
    explicit Fitter(const std::vector<DataPoint>& data);
    void Fit();

private:
    Chi2Calculator chi2;
};
