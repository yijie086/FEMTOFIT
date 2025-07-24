#include "Fitter.h"
#include <cmath>
#include <iostream>

// Static member definitions
double Fitter::q_beam = -1;
double Fitter::L_beam = 1;
double Fitter::L_target = 0;
double Fitter::EB = 5.75;
double Fitter::xB = 0.16;
double Fitter::Q2 = 1.33;
double Fitter::t = -0.21;
double Fitter::theta_Tpol = 0;
double Fitter::phi_Tpol = 0;

Fitter *Fitter::_instance = nullptr;

Fitter::Fitter(const std::vector<InputDataPoint> &data)
    : _data(data), _dvcs(-1, 1, 0, 7.546, 0.16, 1.33, -0.21, 0) {
  _instance = this;
}

void Fitter::Fit() {
  TMinuit minuit(1); // Fit ImH, ReH, ImE, ReE, ImHt, ReHt
  minuit.SetFCN(Fitter::chi2Function);

  double arglist[10];
  arglist[0] = 1;
  int ierflg = 0;
  minuit.mnexcm("SET ERR", arglist, 1, ierflg);

  minuit.DefineParameter(0, "ImH", 5.0, 0.001, 0, 10);
  //minuit.DefineParameter(1, "ReH", -7.0, 0.001, -10, 10);
  //minuit.DefineParameter(2, "ImE",  0.0, 0.01, -5,  5);
  // minuit.DefineParameter(2, "ImE",   0.0, 0.01, -20, 20);
  // minuit.DefineParameter(3, "ReE",   0.0, 0.01, -20, 20);
  // minuit.DefineParameter(4, "ImHt",  0.0, 0.01, -20, 20);
  // minuit.DefineParameter(5, "ReHt",  0.0, 0.01, -20, 20);

  minuit.Migrad();

  for (int i = 0; i < 2; ++i) {
    double val, err;

    minuit.GetParameter(i, val, err);
    std::cout << "Param " << i << ": " << val << " ± " << err << std::endl;
    // Calculate xi ≈ xB / (2 - xB) if needed
  }
  double val_ImH, err_ImH, val_ReH, err_ReH;
  minuit.GetParameter(0, val_ImH, err_ImH);
  minuit.GetParameter(1, val_ReH, err_ReH);
  double xi = xB / (2 - xB);

  FitResult result;
  result.xB = xB;
  result.t = t;
  result.xi = xi;
  result.ImH = val_ImH; // retrieved from minuit
  result.ReH = val_ReH;
  result.ImHerr = err_ImH; // error on ImH
  result.ReHerr = err_ReH; // error on ReH
  fitResults.push_back(result);
}

void Fitter::chi2Function(Int_t &npar, Double_t *, Double_t &fval,
                          Double_t *par, Int_t) {
  BMK_DVCS &dvcs = _instance->_dvcs;

  // Set all fitted CFFs
  //dvcs.setUseDefaultCFF(false); // ensure manual mode
  dvcs.setUseDefaultCFF(true);
  dvcs.setImH(par[0]);
  //dvcs.setReH(par[1]);
//dvcs.setImE(par[2]);  // <-- new
dvcs.setUseWhichCFF(true, true, true, true);
  /*dvcs.setImE(par[2]);
  dvcs.setReE(par[3]);
  dvcs.setImHt(par[4]);
  dvcs.setReHt(par[5]);*/

  double chi2 = 0.0;
//std::cout<<"adfddgs"<<std::endl;
  for (const auto &point : _instance->_data) {
    dvcs.setPrimaryVars(q_beam, L_beam, L_target, EB, xB, Q2, t, point.phi);
    double bsa_pred = dvcs.BSA();
    double cs_pred = dvcs.CrossSection();

    double delta_bsa = (point.bsa - bsa_pred) / point.bsa_err;
   // double delta_cs = (point.cs - cs_pred) / point.cs_err;

    chi2 += delta_bsa * delta_bsa ;//+ delta_cs * delta_cs;
  }
  std::cout << "par[0] (ImH): " << par[0] << ", par[1] (ReH): " << par[1]<< " chi2: " << chi2 << std::endl;


  fval = chi2;
}
void Fitter::PlotBSAandCrossSection(const std::string &tag) {
  const int n = _data.size();
  std::vector<double> phi(n), bsa(n), bsa_err(n), cs(n), cs_err(n);
  std::vector<double> bsa_fit(n), cs_fit(n);

  for (int i = 0; i < n; ++i) {
    phi[i] = _data[i].phi;
    bsa[i] = _data[i].bsa;
    bsa_err[i] = _data[i].bsa_err;
    cs[i] = _data[i].cs;
    cs_err[i] = _data[i].cs_err;

    _dvcs.setPrimaryVars(q_beam, L_beam, L_target, EB, xB, Q2, t, phi[i]);
    bsa_fit[i] = _dvcs.BSA();
    cs_fit[i] = _dvcs.CrossSection();
  }

  // ---- BSA Plot ----
  TCanvas *cBSA = new TCanvas("cBSA", "BSA", 800, 600);
  cBSA->SetGrid();

  TH1F *frameBSA = new TH1F("frameBSA", "", 100, 0, 360);
  frameBSA->SetMinimum(-0.3);
  frameBSA->SetMaximum(0.3);
  frameBSA->GetXaxis()->SetTitle("#phi [deg]");
  frameBSA->GetYaxis()->SetTitle("BSA (#sigma_{LU})");
  frameBSA->SetStats(0);
  frameBSA->Draw();

  TGraphErrors *gBSA = new TGraphErrors(n, phi.data(), bsa.data(), nullptr, bsa_err.data());
  gBSA->SetMarkerStyle(20);
  gBSA->SetMarkerColor(kBlue+1);
  gBSA->Draw("P SAME");

  TGraph *gBSA_fit = new TGraph(n, phi.data(), bsa_fit.data());
  gBSA_fit->SetLineColor(kRed);
  gBSA_fit->SetLineWidth(2);
  gBSA_fit->Draw("L SAME");

  TLatex labelBSA;
  labelBSA.SetNDC();
  labelBSA.SetTextSize(0.045);
  labelBSA.DrawLatex(0.15, 0.85, Form("-t = %.2f, x_{B} = %.2f, Q^{2} = %.2f", t, xB, Q2));

  cBSA->SaveAs(Form("BSA_fit_%s.png", tag.c_str()));

  // ---- Cross Section Plot ----
  TCanvas *cCS = new TCanvas("cCS", "Cross Section", 800, 600);
  cCS->SetGrid();

  double maxCS = *std::max_element(cs.begin(), cs.end()) * 1.3;
  TH1F *frameCS = new TH1F("frameCS", "", 100, 0, 360);
  frameCS->SetMinimum(0.1);
  frameCS->SetMaximum(maxCS);
  frameCS->GetXaxis()->SetTitle("#phi [deg]");
  frameCS->GetYaxis()->SetTitle("d^{4} #sigma GeV^{-4}");
  frameCS->SetStats(0);
  gPad->SetLogy();
  frameCS->Draw();
  
  TGraphErrors *gCS = new TGraphErrors(n, phi.data(), cs.data(), nullptr, cs_err.data());
  gCS->SetMarkerStyle(21);
  gCS->SetMarkerColor(kGreen+2);
  gCS->Draw("P SAME");

  TGraph *gCS_fit = new TGraph(n, phi.data(), cs_fit.data());
  gCS_fit->SetLineColor(kMagenta);
  gCS_fit->SetLineWidth(2);
  gCS_fit->Draw("L SAME");

  TLatex labelCS;
  labelCS.SetNDC();
  labelCS.SetTextSize(0.045);
  labelCS.DrawLatex(0.15, 0.85, Form("-t = %.2f, x_{B} = %.2f, Q^{2} = %.2f", t, xB, Q2));

  cCS->SaveAs(Form("CS_fit_%s.png", tag.c_str()));

  // Clean up
  delete gBSA;
  delete gBSA_fit;
  delete gCS;
  delete gCS_fit;
  delete cBSA;
  delete cCS;
  delete frameBSA;
  delete frameCS;
}
