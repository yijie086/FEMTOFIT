#ifndef FITTER_H
#define FITTER_H

#include <vector>
#include "BMK_DVCS.h"
#include "TMinuit.h"
#include "DataHandler.h"

#include "TGraph2D.h"
#include <TH1F.h>
#include "TCanvas.h"
#include <TLegend.h>
#include <TLatex.h>
#include "TGraph.h"
#include <TAxis.h> 
#include <TGraphErrors.h>

//std::vector<std::tuple<double, double, double>> fittedImH; // (xB, t, value)
//std::vector<std::tuple<double, double, double>> fittedReH;
// Repeat for ImE, ReE, etc.

class Fitter {
public:
    Fitter(const std::vector<InputDataPoint>& data);
    void setFitPrimaryVars(double rq_beam, double rL_beam, double rL_target,
                      double rEB, double rxB, double rQ2, double rt, double rtheta_Tpol = 0,
                      double rphi_Tpol = 0){
                        q_beam= rq_beam; 
                        L_beam = rL_beam;
                        L_target= rL_target;
                        EB= rEB;
                        xB= rxB;
                        Q2= rQ2;
                        t= rt;
                        theta_Tpol = rtheta_Tpol;
                        phi_Tpol = rphi_Tpol;
                      }
    void Fit();
    void GetCFFs();
    void PlotBSAandCrossSection(const std::string &tag);
    const std::vector<FitResult>& getFitResults() const { return fitResults; }
private:
    std::vector<InputDataPoint> _data;
    BMK_DVCS _dvcs;
   /*struct CFF{
    ImCFF
    ReCFF
   }*/
   static double q_beam ;
   static double L_beam ;
   static double L_target;
   static double EB;
   static double xB;
   static double Q2;
   static double t;
   static double theta_Tpol;
   static double phi_Tpol;
    std::vector<FitResult> fitResults;
    static void chi2Function(Int_t& npar, Double_t* grad, Double_t& fval, Double_t* par, Int_t flag);
    static Fitter* _instance;  // to access inside static chi2Function
};

#endif
