#ifndef BMK_DVCS_H
#define BMK_DVCS_H

#include "TMath.h"
#include <cmath>

/*double PI = TMath::Pi();

double alpha = 1./137.036;
double alpha3 = TMath::Power(alpha,3.);
double hbarc2 = 0.38938;//GeV^{2} mbarn

double m = 0.000511;
double M = 0.93827;
double muP = 2.79285;*/

double GetF1_FX(double T);
double GetF2_FX(double T);
double GetGMP_FX(double tau);
double GetGEP_FX(double tau);

double GetImH_FX(double xi, double t);
double GetImHt_FX(double xi, double t);
double GetImE_FX(double xi, double t);
double GetImEt_FX(double xi, double t);

double GetReH_FX(double xi, double t);
double GetReHt_FX(double xi, double t);
double GetReE_FX(double xi, double t);
double GetReEt_FX(double xi, double t);

// extern bool hasH_FX, hasHt_FX, hasE_FX, hasEt_FX;

////////////////////////////////////////////////////////////////////////////////////////////////////////////

class BMK_DVCS {
public:
  double q_beam, L_beam, L_target;
  double EB, xB, Q2, t, phi, theta_Tpol,
      phi_Tpol; // primary variables ; EB on fixed target ; Trento convention

  double xi, nu, y, eps, eps2, phi_BMK, t_min, K2, K, J, Ktild2,
      Ktilda;   // secondary variables ;  phi_BMK = pi - phi_Trento
  double Jacob; // Jacobian from (xB,y) to (xB,Q2)
  double F1, F2, FF_comb1, FF_comb2, FF_comb3;
  double ImH, ImHt, ImE, ImEt;
  double ReH, ReHt, ReE, ReEt;

  void setUseWhichCFF(bool IsH = true, bool IsHt = true, bool IsE = true,
                      bool IsEt = true) {
    hasH_FX = IsH;
    hasHt_FX = IsHt;
    hasE_FX = IsE;
    hasEt_FX = IsEt;
  }

  void setImH(double val) {
    ImH = val;
    FitImH = true;
  }
  void setReH(double val) {
    ReH = val;
    FitReH = true;
  }
  void setImHt(double val) {
    ImHt = val;
    FitImHt = true;
  }
  void setReHt(double val) {
    ReHt = val;
    FitReHt = true;
  }
  void setImE(double val) {
    ImE = val;
    FitImE = true;
  }
  void setReE(double val) {
    ReE = val;
    FitReE = true;
  }
  void setImEt(double val) {
    ImEt = val;
    FitImEt = true;
  }
  void setReEt(double val) {
    ReEt = val;
    FitReEt = true;
  }

  void setUseDefaultCFF(bool useCFF = false) { IsDefaultCFF = useCFF; }
  bool VERB;

  BMK_DVCS(double rq_beam, double rL_beam, double rL_target, double rEB,
           double rxB, double rQ2, double rt, double rphi,
           double rtheta_Tpol = 0, double rphi_Tpol = 0);
  void setSecondaryVars(void);
  void setPrimaryVars(double rq_beam, double rL_beam, double rL_target,
                      double rEB, double rxB, double rQ2, double rt,
                      double rphi, double rtheta_Tpol = 0,
                      double rphi_Tpol = 0);

  double CrossSection(void);
  double TPolCrossSection(void);
  double BSA(void);
  double pBSA(void);
  double TLSA(void);
  double TLLSA(void);
  double TTSAx(void);
  double TTSAy(void);
  double TTSSAx(void);
  double TTSSAy(void);
  double BCA(void);
  double BCSA(void);
  double BC0SA(void);

  double BCLA(void);
  double BCLLA(void);
  double BCTxA(void);
  double BCTyA(void);

  double T2(void);
  double BH2(void);
  double DVCS2(void);
  double BHDVCS(void);

  double c0_BH(void);
  double c1_BH(void);
  double c2_BH(void);
  double c0_BH_LP(void);
  double c1_BH_LP(void);
  double c0_BH_TP(void);
  double c1_BH_TP(void);
  double s1_BH_TP(void);
  double BHP1(void);
  double BHP2(void);

  double c0_I(void);
  double c1_I(void);
  double s1_I(void);
  double c0_I_LP(void);
  double c1_I_LP(void);
  double s1_I_LP(void);
  double c0_I_TP(void);
  double c1_I_TP(void);
  double s1_I_TP(void);

  double c0_DVCS(void);
  double c0_DVCS_LP(void);
  double c0_DVCS_TP(void);

  double GetF1_FX(double T);

  double GetF2_FX(double T);
  double GetGMP_FX(double tau);
  double GetGEP_FX(double tau);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double GetImH_FX(double xi, double t);

  double GetImHt_FX(double xi, double t);

  double GetImE_FX(double xi, double t);

  double GetImEt_FX(double xi, double t);
  double GetReH_FX(double xi, double t);
  double GetReHt_FX(double xi, double t);

  double GetReE_FX(double xi, double t);

  double GetReEt_FX(double xi, double t);

  static constexpr double PI = TMath::Pi();
  static constexpr double alpha = 1. / 137.036;
  static constexpr double alpha3 = std::pow(alpha, 3.0);
  static constexpr double hbarc2 = 0.38938; // GeV^2 mbarn
  static constexpr double m = 0.000511;
  static constexpr double M = 0.93827;
  static constexpr double muP = 2.79285;

private:
  bool IsDefaultCFF = true;
  bool FitImH = false;
  bool FitImHt = false;
  bool FitImE = false;
  bool FitImEt = false;
  bool FitReH = false;
  bool FitReHt = false;
  bool FitReE = false;
  bool FitReEt = false;
  //
  bool hasH_FX = true;
  bool hasHt_FX = true;
  bool hasE_FX = true;
  bool hasEt_FX = true;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif // BMK_DVCS_H