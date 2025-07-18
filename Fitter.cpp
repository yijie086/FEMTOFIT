#include "Fitter.h"
#include <iostream>
#include <vector>
#include <cmath>

Fitter::Fitter(const std::vector<DataPoint>& data) : chi2(data) {}

void Fitter::Fit() {
    std::vector<double> guess(12, 0.1); // 12 parameters for polynomial model
    double chi2val = chi2(guess);

    std::cout << "Initial chi2: " << chi2val << std::endl;
    std::cout << "Initial guess: ";
    for (auto v : guess) std::cout << v << " ";
    std::cout << std::endl;

    // Placeholder update — simulate slight change
    guess[0] += 0.01;
    chi2val = chi2(guess);

    std::cout << "Updated chi2: " << chi2val << std::endl;
    std::cout << "Updated guess: ";
    for (auto v : guess) std::cout << v << " ";
    std::cout << std::endl;
}
