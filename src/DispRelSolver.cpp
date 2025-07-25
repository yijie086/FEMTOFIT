#include "DispRelSolver.h"
#include <cmath>
#include <stdexcept>
#include <iostream>

DispRelSolver::DispRelSolver() {}

void DispRelSolver::addCFFs(double t, const std::vector<FitResult>& results) {
    CFFData data;
    data.rawFits = results;

    for (const auto& r : results) {
        data.xis.push_back(r.xi);
        data.ReHs.push_back(r.ReH);
    }

    if (data.xis.size() < 3) {
        std::cerr << "[WARN] Not enough points to fit spline at t = " << t << "\n";
        return;
    }

    data.reH_spline = std::make_shared<TSpline3>(
        ("spline_t_" + std::to_string(t)).c_str(),
        data.xis.data(), data.ReHs.data(), static_cast<int>(data.xis.size())
    );

    cff_data_by_t[t] = std::move(data);
}

double DispRelSolver::dispersionIntegral(const TSpline3* spline, double xi_eval) const {
    if (!spline) throw std::runtime_error("Invalid spline pointer");

    const int N = 1000;
    const double xi_min = 1e-4;
    const double xi_max = 1.0;
    const double dxi = (xi_max - xi_min) / N;

    double integral = 0.0;

    for (int i = 0; i < N; ++i) {
        double xi = xi_min + i * dxi + dxi / 2.0;

        double num = spline->Eval(xi);
        double denom = xi * xi - xi_eval * xi_eval;

        if (std::fabs(denom) < 1e-8) continue; // skip near-pole

        integral += xi * num / denom * dxi;
    }

    integral *= 2.0 * xi_eval / M_PI;

    return integral;
}

double DispRelSolver::getDterm(double t, double xi_eval) const {
    auto it = cff_data_by_t.find(t);
    if (it == cff_data_by_t.end()) throw std::runtime_error("No data for t = " + std::to_string(t));

    const auto& data = it->second;
    if (!data.reH_spline) throw std::runtime_error("Spline not initialized");

    double disp = dispersionIntegral(data.reH_spline.get(), xi_eval);

    double ReH_at_xi = data.reH_spline->Eval(xi_eval);

#ifdef DEBUG_DISPREL
    std::cerr << "[DEBUG] ReH(xi) = " << ReH_at_xi << ", Disp. integral = " << disp
              << ", D(t) = " << (ReH_at_xi - disp) << "\n";
#endif

    return ReH_at_xi - disp;
}

double DispRelSolver::getDterm(double t) const {
    auto it = cff_data_by_t.find(t);
    if (it == cff_data_by_t.end()) throw std::runtime_error("No data for t = " + std::to_string(t));

    const auto& data = it->second;
    if (data.rawFits.empty()) throw std::runtime_error("No fit data at t = " + std::to_string(t));

    double xi = data.rawFits.front().xi;
    return getDterm(t, xi);
}
