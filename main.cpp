// STL
#include <iostream>

#include "DataHandler.h"
#include "Fitter.h"

/// models for fitting
#include "BMKModel.h"
#include "BMK_DVCS.h"
#include "KM15Model.h"
#include "Models.h"
#include "DispRelSolver.h"
#include <vector>

int main() {
  // DataHandler dh("data/sample_data.csv");
  // Fitter fitter(dh.GetData());
  // fitter.Fit();

  Models modelManager;
  // modelManager.registerModel(std::make_shared<KM15Model>());

  // auto model = modelManager.getModel("KM15Model");
  // model->setKinematics(2.0, -0.3, 0.1, 10.6);

  // std::cout << "BSA: " << model->getBSA() << std::endl;
  // std::cout << "Cross Section: " << model->getCrossSection() << std::endl;

  // Beam energy (GeV)
  double EB = 10.6;

  // Define kinematics: xB, Q2, t, phi
  double xB = 0.23;
  double Q2 = 2.2;
  double t = -0.322;
  double phi = 180.0; // degrees

  // Beam and target polarizations
  double q_beam = -1;  // electron
  double L_beam = 1;   // +1 or -1 for helicity
  double L_target = 0; // unpolarized

  // Create the DVCS object
  // BMKModel dvcs(q_beam, L_beam, L_target, EB, xB, Q2, t, phi);
  BMK_DVCS dvcs_FX(q_beam, L_beam, L_target, EB, xB, Q2, t, phi);

  // Optional: enable verbose output
  // dvcs.VERB = true;

  // Compute observables
  // double cross_sec = dvcs.CrossSection(); // nb/sr/GeV^4 or so
  // double bsa = dvcs.BSA();                // Beam Spin Asymmetry

  double cross_sec_FX = dvcs_FX.CrossSection(); // nb/sr/GeV^4 or so
  double bsa_FX = dvcs_FX.BSA();                // Beam Spin Asymmetry

  std::cout << "DVCS Kinematics:\n";
  std::cout << "xB = " << xB << ", Q2 = " << Q2 << " GeV^2, t = " << t
            << " GeV^2, phi = " << phi << " deg\n";
  std::cout << "Unpolarized Cross Section = "
            << "  Unpolarized Cross Section FX = " << cross_sec_FX << "nb\n";
  std::cout << "Beam Spin Asymmetry (BSA) = "
            << "  Beam Spin Asymmetry (BSA) FX = " << bsa_FX << "\n";


  /// chatgpt
  // Fitter should take the data BSA and Cross-section and fit simulatenously
  // return 0;

  std::vector<InputDataPoint> sampleData = {
      //{  0.0,   -7.447817e-17, 1.995165, 0.03 * std::abs(-7.447817e-17), 0.03
      //* 1.995165 },
      {20.0, 1.043311e-01, 1.531820, 0.03 * 1.043311e-01, 0.03 * 1.531820},
      {40.0, 1.937151e-01, 0.895475, 0.03 * 1.937151e-01, 0.03 * 0.895475},
      {60.0, 2.538035e-01, 0.531500, 0.03 * 2.538035e-01, 0.03 * 0.531500},
      {80.0, 2.749693e-01, 0.349975, 0.03 * 2.749693e-01, 0.03 * 0.349975},
      {100.0, 2.573233e-01, 0.257283, 0.03 * 2.573233e-01, 0.03 * 0.257283},
      {120.0, 2.102322e-01, 0.207848, 0.03 * 2.102322e-01, 0.03 * 0.207848},
      {140.0, 1.460321e-01, 0.181120, 0.03 * 1.460321e-01, 0.03 * 0.181120},
      {160.0, 7.423447e-02, 0.167758, 0.03 * 7.423447e-02, 0.03 * 0.167758},
      //{180.0,    0.000000e+00, 0.163688, 0.03 * 0.000000e+00, 0.03 * 0.163688
      //},
      {200.0, -7.423447e-02, 0.167758, 0.03 * 7.423447e-02, 0.03 * 0.167758},
      {220.0, -1.460321e-01, 0.181120, 0.03 * 1.460321e-01, 0.03 * 0.181120},
      {240.0, -2.102322e-01, 0.207848, 0.03 * 2.102322e-01, 0.03 * 0.207848},
      {260.0, -2.573233e-01, 0.257283, 0.03 * 2.573233e-01, 0.03 * 0.257283},
      {280.0, -2.749693e-01, 0.349975, 0.03 * 2.749693e-01, 0.03 * 0.349975},
      {300.0, -2.538035e-01, 0.531500, 0.03 * 2.538035e-01, 0.03 * 0.531500},
      {320.0, -1.937151e-01, 0.895475, 0.03 * 1.937151e-01, 0.03 * 0.895475},
      {340.0, -1.043311e-01, 1.531820, 0.03 * 1.043311e-01, 0.03 * 1.531820},
      //{360.0,   -7.447817e-17, 1.995165, 0.03 * std::abs(-7.447817e-17), 0.03
      //* 1.995165 }
  };

  /*std::vector<InputDataPoint> realsampleData = {
   {68.33  ,  8.9354 , 0.7223  , },
   {82.35  ,  5.7024 , 0.3494  , },
   {97.13  ,  3.5973 , 0.2525  , },
   {112.2  ,  2.9545 , 0.1979  , },
   {127.27 ,   2.604 , 0.2255  , },
   {142.33 ,   2.1674,  0.142  , },
   {157.4  ,  2.1833 , 0.1295  , },
   {172.47 ,   1.9659,  0.196  , },
   {187.53 ,   2.3694,  0.224  , },
   {202.6  ,  2.1443 , 0.1475  , },
   {217.67 ,   2.4837,  0.234  , },
   {232.73 ,   2.4709,  0.184  , },
   {247.8  ,  3.0839 , 0.2447  , },
   {262.87 ,   3.4295,  0.218  , },
   {277.65 ,   5.3907,  0.372  , },
   {291.67 ,   9.0739,  0.942  , },
  };*/

  Fitter fitter2(sampleData);
  fitter2.setFitPrimaryVars(q_beam, L_beam, L_target, EB, xB, Q2, t);
  fitter2.Fit();

  if (!fitter2.getFitResults().empty()) {
    const auto &res = fitter2.getFitResults().back(); // last result

    std::cout << "FitResult:"
              << " xB = " << res.xB << ", t = " << res.t << ", xi = " << res.xi
              << ", ImH = " << res.ImH << ", ReH = " << res.ReH << std::endl;
  }

  std::cout<< "running\n";
  ///
  EB = 5.75;
  DataHandler dh;
  auto kinematicBlocks = dh.loadData("/w/hallb-scshelf2102/clas12/singh/Softwares/FemtoFit/source/data/output.csv");
  //auto kinematicBlocks = dh.loadData("/work/clas12/singh/Softwares/newFemtoFit/FemtoFit/source/data/rgkdata1.csv");

  //Fitter fitter({});
  std::vector<FitResult> allFitResults;

  for (const auto &block : kinematicBlocks) {
    Fitter fitter(block.points);
    fitter.setFitPrimaryVars(q_beam, L_beam, L_target, EB, block.xB, block.Q2,
                             block.t);
    std::cout<<q_beam<<" "<<EB<<std::endl;               
    std::cout<<"block.t"<< block.t<<   "block.xB"<<block.xB<<" block.Q2"<< block.Q2<<std::endl;  
    for (int i=0;i<block.points.size();i++){
    std::cout<<"block.bsa.     "<< block.points[i].bsa<<   "block.xs.  "<<block.points[i].cs<<" block.phi.      "<< block.points[i].phi<<std::endl;
    }                  
    fitter.Fit();
    //fitter.PlotBSAandCrossSection(Form("CLAS6t%.2f_xB%.2f_Q2%.2f", block.t, block.xB, block.Q2));

    if (!fitter.getFitResults().empty()) allFitResults.push_back(fitter.getFitResults().back());
    
  }

  std::map<double, std::vector<FitResult>> resultsByT;

  for (const auto &res : allFitResults) {
    resultsByT[res.t].push_back(res);
        std::cout << "FitResult:"
              << " xB = " << res.xB << ", t = " << res.t << ", xi = " << res.xi
              << ", ImH = " << res.ImH << "+-" << res.ImHerr << ", ReH = " << res.ReH << "+-" << res.ReHerr << std::endl;
  }

 for (const auto &[tval, fitVec] : resultsByT) {
  TCanvas *cImH = new TCanvas("cImH", "ImH vs xi", 600, 600);
  cImH->SetGrid();

  // Create dummy histogram as frame
  TH1F *frame = new TH1F("frame", "", 100, 0.0, 0.3);  // dummy 100 bins
  //frame->SetMinimum(0.0);
  //frame->SetMaximum(10.0);
  frame->GetXaxis()->SetRangeUser(0.0, 0.3);
  frame->GetYaxis()->SetRangeUser(0.0, 10.0);
  frame->GetXaxis()->SetTitle("#xi");
  frame->GetYaxis()->SetTitle("ImH");
  frame->GetXaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleSize(0.045);
  frame->GetXaxis()->SetLabelSize(0.04);
  frame->GetYaxis()->SetLabelSize(0.04);
  frame->SetStats(0);
  frame->Draw();  // Draw only frame

  // Create TGraphErrors
  TGraphErrors *gImH = new TGraphErrors(fitVec.size());
  for (size_t i = 0; i < fitVec.size(); ++i) {
    gImH->SetPoint(i, fitVec[i].xi, fitVec[i].ImH);
    gImH->SetPointError(i, 0, fitVec[i].ImHerr);
  }

  gImH->SetMarkerStyle(20);
  gImH->SetMarkerColor(kBlue+1);
  gImH->SetLineColor(kBlue+1);
  gImH->SetLineWidth(2);
  gImH->Draw("P SAME");  // Plot over the frame

  // Legend
  TLegend *legend = new TLegend(0.6, 0.65, 0.78, 0.88);
  legend->AddEntry(gImH, "ImH", "lp");
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->Draw();
///
  TLatex latex;
  latex.SetNDC();  // Use normalized device coordinates (0 to 1)
  latex.SetTextSize(0.05);
  latex.SetTextAlign(13);  // Align left, bottom
  latex.DrawLatex(0.5, 0.85, Form("-t = %.2f GeV^{2}", tval));

  // Save
  cImH->SaveAs(Form("CLAS6_ImH_vs_xi_t%.2f.png", tval));

  delete gImH;
  delete legend;
  delete frame;
  delete cImH;
}
for (const auto &[tval, fitVec] : resultsByT) {
  TCanvas *cReH = new TCanvas("cReH", "ReH vs xi", 600, 600);
  cReH->SetGrid();

  // Create dummy histogram as frame
  TH1F *frame2 = new TH1F("frame2", "", 100, 0.0, 0.3);  // dummy 100 bins
  frame2->GetXaxis()->SetRangeUser(0.0, 0.3);
  frame2->GetYaxis()->SetRangeUser(-10, 0.0);
  frame2->GetXaxis()->SetTitle("#xi");
  frame2->GetYaxis()->SetTitle("ImH");
  frame2->GetXaxis()->SetTitleSize(0.045);
  frame2->GetYaxis()->SetTitleSize(0.045);
  frame2->GetXaxis()->SetLabelSize(0.04);
  frame2->GetYaxis()->SetLabelSize(0.04);
  frame2->SetStats(0);
  frame2->Draw();  // Draw only frame

  // Create TGraphErrors
  TGraphErrors *gReH = new TGraphErrors(fitVec.size());
  for (size_t i = 0; i < fitVec.size(); ++i) {
    gReH->SetPoint(i, fitVec[i].xi, fitVec[i].ReH);
    gReH->SetPointError(i, 0, fitVec[i].ReHerr);
  }

  gReH->SetMarkerStyle(20);
  gReH->SetMarkerColor(kBlue+1);
  gReH->SetLineColor(kBlue+1);
  gReH->SetLineWidth(2);
  gReH->Draw("P SAME");  // Plot over the frame

  // Legend
  TLegend *legend2 = new TLegend(0.6, 0.65, 0.78, 0.88);
  legend2->AddEntry(gReH, "ReH", "lp");
  legend2->SetBorderSize(0);
  legend2->SetFillStyle(0);
  legend2->Draw();
///
  TLatex latex2;
  latex2.SetNDC();  // Use normalized device coordinates (0 to 1)
  latex2.SetTextSize(0.05);
  latex2.SetTextAlign(13);  // Align left, bottom
  latex2.DrawLatex(0.5, 0.85, Form("-t = %.2f GeV^{2}", tval));

  // Save
  cReH->SaveAs(Form("CLAS6_ReH_vs_xi_t%.2f.png", tval));

  delete gReH;
  delete legend2;
  delete frame2;
  delete cReH;
}

// === After computing and adding D-terms ===
DispRelSolver drs;
std::vector<double> tvals, Dvals;

for (const auto& [tval, fits] : resultsByT) {
    drs.addCFFs(tval, fits);
    try {
        double D = drs.getDterm(tval);
        tvals.push_back(-tval);
        Dvals.push_back(D);
        std::cout << "t = " << -tval << ", D(t) = " << D << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Failed to compute D(t) for t = " << tval << ": " << e.what() << std::endl;
    }
}


//plotting D terms as function of xi and t
for (const auto& [tval, fitVec] : resultsByT) {
    TCanvas* cDxi = new TCanvas("cDxi", "D(t, #xi) vs #xi", 600, 600);
    cDxi->SetGrid();

    // Create frame histogram
    TH1F* frameXi = new TH1F("frameXi", "", 100, 0.0, 0.3);
    frameXi->GetXaxis()->SetTitle("#xi");
    frameXi->GetYaxis()->SetTitle("D(t, #xi)");
    frameXi->GetXaxis()->SetTitleSize(0.045);
    frameXi->GetYaxis()->SetTitleSize(0.045);
    frameXi->GetYaxis()->SetTitleOffset(0.6);
    frameXi->GetXaxis()->SetLabelSize(0.04);
    frameXi->GetYaxis()->SetLabelSize(0.04);
    frameXi->SetStats(0);
    frameXi->GetXaxis()->SetRangeUser(0.0, 0.3);
    frameXi->GetYaxis()->SetRangeUser(-15, 5);
    frameXi->Draw();

    TGraph* gD_vs_xi = new TGraph();
    int idx = 0;
    for (const auto& fit : fitVec) {
        try {
            double D = drs.getDterm(tval, fit.xi);
            gD_vs_xi->SetPoint(idx++, fit.xi, D);
        } catch (...) {}
    }

    gD_vs_xi->SetMarkerStyle(21);
    gD_vs_xi->SetMarkerColor(kGreen+2);
    gD_vs_xi->SetLineColor(kGreen+2);
    gD_vs_xi->SetLineWidth(2);
    gD_vs_xi->Draw("P SAME");

    // Add legend
    TLegend* legXi = new TLegend(0.6, 0.65, 0.78, 0.85);
    legXi->AddEntry(gD_vs_xi, "D(t, #xi)", "lp");
    legXi->SetBorderSize(0);
    legXi->SetFillStyle(0);
    legXi->Draw();

    // Annotate with t-value
    TLatex latexXi;
    latexXi.SetNDC();
    latexXi.SetTextSize(0.05);
    latexXi.SetTextAlign(13);
    latexXi.DrawLatex(0.5, 0.85, Form("-t = %.2f GeV^{2}", tval));

    // Save
    cDxi->SaveAs(Form("CLAS6_Dterm_vs_xi_t%.2f.png", tval));

    // Cleanup
    delete gD_vs_xi;
    delete legXi;
    delete frameXi;
    delete cDxi;
}
// === Plot D(t) ===
TCanvas* cD = new TCanvas("cD", "D-term vs t", 600, 600);
cD->SetGrid();

double tmin = *std::min_element(tvals.begin(), tvals.end());
double tmax = *std::max_element(tvals.begin(), tvals.end());

TH1F* frameD = new TH1F("frameD", "", 100, tmin - 0.05, tmax + 0.05);
frameD->GetXaxis()->SetTitle("-t [GeV^{2}]");
frameD->GetYaxis()->SetTitle("D(t)");
frameD->GetXaxis()->SetTitleSize(0.045);
frameD->GetYaxis()->SetTitleSize(0.045);
frameD->GetYaxis()->SetTitleOffset(0.6);
frameD->GetXaxis()->SetLabelSize(0.04);
frameD->GetYaxis()->SetLabelSize(0.04);
frameD->GetXaxis()->SetRangeUser(0.0, 0.6);
frameD->GetYaxis()->SetRangeUser(-15, 2.0);
frameD->SetStats(0);
frameD->Draw();

TGraph* gD = new TGraph(tvals.size());
for (size_t i = 0; i < tvals.size(); ++i) {
    gD->SetPoint(i, tvals[i], Dvals[i]);
}

gD->SetMarkerStyle(21);
gD->SetMarkerColor(kRed + 1);
gD->SetLineColor(kRed + 1);
gD->SetLineWidth(2);
gD->Draw("P SAME");

TLegend* legD = new TLegend(0.6, 0.65, 0.78, 0.85);
legD->AddEntry(gD, "D(t)", "lp");
legD->SetBorderSize(0);
legD->SetFillStyle(0);
legD->Draw();

cD->SaveAs("CLAS6_Dterm_vs_t.png");

// Cleanup
delete gD;
delete legD;
delete frameD;
delete cD;

}
