#ifndef BMKModel_H
#define BMKModel_H

#include "TMath.h"
#include <cmath>

class BMKModel {
public:
    // Beam and target setup
    double q_beam, L_beam, L_target;
    double EB, xB, Q2, t, phi, theta_Tpol, phi_Tpol;

    // Derived kinematic quantities
    double xi, nu, y, eps, eps2, phi_BMK, t_min, K2, K, J, Ktild2, Ktilda, Jacob;
    double F1, F2, FF_comb1, FF_comb2, FF_comb3;
    double ImH, ImHt, ImE, ImEt;
    double ReH, ReHt, ReE, ReEt;

    bool VERB;

    // Constructor
    BMKModel(double rq_beam, double rL_beam, double rL_target,
             double rEB, double rxB, double rQ2, double rt, double rphi,
             double rtheta_Tpol = 0, double rphi_Tpol = 0);

    // Setters
    void setPrimaryVars(double rq_beam, double rL_beam, double rL_target,
                        double rEB, double rxB, double rQ2, double rt, double rphi,
                        double rtheta_Tpol = 0, double rphi_Tpol = 0);
    void setSecondaryVars(void);

    // Main observables
    double CrossSection(void);
    double BSA(void);

    // DVCS-related terms
    double T2(void);
    double BH2(void);
    double DVCS2(void);
    double BHDVCS(void);

    // Harmonic terms
    double c0_BH(void);
    double c1_BH(void);
    double c2_BH(void);
    double c0_I(void);
    double c1_I(void);
    double s1_I(void);
    double c0_DVCS(void);

    // BH denominators
    double BHP1(void);
    double BHP2(void);

    static constexpr double PI = TMath::Pi();
    static constexpr double alpha = 1. / 137.036;
    static constexpr double alpha3 = std::pow(alpha, 3.0);
    static constexpr double hbarc2 = 0.38938; // GeV^2 mbarn
    static constexpr double m = 0.000511;
    static constexpr double M = 0.93827;
    static constexpr double muP = 2.79285;
};

// Global controls and external CFF prototypes
extern bool hasH, hasHt, hasE, hasEt;
extern double renormImag, renormReal;

double GetF1(double T);
double GetF2(double T);
double GetGMP(double tau);
double GetGEP(double tau);
double GetImH(double xi, double t);
double GetImHt(double xi, double t);
double GetImE(double xi, double t);
double GetImEt(double xi, double t);
double GetReH(double xi, double t);
double GetReHt(double xi, double t);
double GetReE(double xi, double t);
double GetReEt(double xi, double t);

#endif // BMKModel_H
