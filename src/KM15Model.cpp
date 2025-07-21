#include "KM15Model.h"
#include "Gepard.h"

KM15Model::KM15Model() {}
KM15Model::~KM15Model() {}

void KM15Model::setKinematics(double Q2, double t, double xb, double E) {
    Gepard::call_KM15(Q2, t, xb, E, lastBSA, lastXS);
}

double KM15Model::getBSA() const {
    return lastBSA;
}

double KM15Model::getCrossSection() const {
    return lastXS;
}

std::string KM15Model::getName() const {
    return "KM15Model";
}
