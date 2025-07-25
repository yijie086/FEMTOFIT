#pragma once

#include <map>
#include <vector>
#include <memory>
#include "DataHandler.h"
#include "TSpline.h"

class DispRelSolver {
public:
    DispRelSolver();

    // Store CFFs (real part vs xi) for a given t
    void addCFFs(double t, const std::vector<FitResult>& results);

    // Compute D-term using full dispersion relation (spline interpolated)
    double getDterm(double t) const;

    // Compute D-term at a specific xi
    double getDterm(double t, double xi_eval) const;

private:
    struct CFFData {
        std::vector<double> xis;
        std::vector<double> ReHs;
        std::shared_ptr<TSpline3> reH_spline;
        std::vector<FitResult> rawFits;  // original points
    };

    std::map<double, CFFData> cff_data_by_t;

    double dispersionIntegral(const TSpline3* spline, double xi_eval) const;
};
