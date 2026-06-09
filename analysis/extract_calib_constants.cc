#include <array>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TLine.h"

namespace {

constexpr double kTargetDepositGeV = 2.396;

struct FitResult {
  bool ok = false;
  int status = -999;
  double amplitude = 0.0;
  double mean = 0.0;
  double meanErr = 0.0;
  double sigma = 0.0;
  double sigmaErr = 0.0;
  double constant = 0.0;
  double fitMin = 0.0;
  double fitMax = 0.0;
  double chi2 = 0.0;
  int ndf = 0;
  std::string histName;
  std::string reason;
};

FitResult FitHistogram(TH1 *h, const std::string &label) {
  FitResult out;
  out.histName = h ? h->GetName() : "null";

  if (!h) {
    out.reason = "histogram not found";
    return out;
  }
  if (h->GetEntries() < 50) {
    out.reason = "not enough entries";
    return out;
  }

  int peakBin = h->GetMaximumBin();
  // intADC_after 히스토그램에서 저ADC 스파이크(노이즈/엣지)를 피하기 위해
  // 5k ADC 이상 구간의 최대 bin을 우선 피크 후보로 사용한다.
  constexpr double kPeakSearchMinX = 5000.0;
  int peakBinAboveThreshold = -1;
  double peakContentAboveThreshold = -1.0;
  for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
    const double x = h->GetXaxis()->GetBinCenter(ib);
    if (x < kPeakSearchMinX) continue;
    const double c = h->GetBinContent(ib);
    if (c <= 0.0) continue;
    if (c > peakContentAboveThreshold) {
      peakContentAboveThreshold = c;
      peakBinAboveThreshold = ib;
    }
  }
  if (peakBinAboveThreshold > 0) peakBin = peakBinAboveThreshold;
  const double peakX = h->GetXaxis()->GetBinCenter(peakBin);
  double rms = h->GetRMS();
  const double xmin = h->GetXaxis()->GetXmin();
  const double xmax = h->GetXaxis()->GetXmax();
  const double width = xmax - xmin;
  if (rms <= 0.0) rms = width / 12.0;

  double fitMin = std::max(xmin, peakX - 1.5 * rms);
  double fitMax = std::min(xmax, peakX + 1.5 * rms);
  if (fitMax <= fitMin) {
    fitMin = xmin;
    fitMax = xmax;
  }

  TF1 fgaus("fgaus", "gaus", fitMin, fitMax);
  int fitStatus = h->Fit(&fgaus, "RQ0");
  if (fitStatus != 0) {
    fitMin = std::max(xmin, peakX - 2.5 * rms);
    fitMax = std::min(xmax, peakX + 2.5 * rms);
    fgaus.SetRange(fitMin, fitMax);
    fitStatus = h->Fit(&fgaus, "RQ0");
  }

  if (fitStatus != 0) {
    out.reason = "gaussian fit failed";
    return out;
  }

  out.status = fitStatus;
  out.fitMin = fitMin;
  out.fitMax = fitMax;
  out.amplitude = fgaus.GetParameter(0);
  out.mean = fgaus.GetParameter(1);
  out.meanErr = fgaus.GetParError(1);
  out.sigma = fgaus.GetParameter(2);
  out.sigmaErr = fgaus.GetParError(2);
  out.chi2 = fgaus.GetChisquare();
  out.ndf = fgaus.GetNDF();

  if (out.sigma <= 0.0) {
    out.reason = "non-positive sigma";
    return out;
  }

  if (out.mean <= 0.0) {
    out.reason = "non-positive mean";
    return out;
  }

  out.constant = kTargetDepositGeV / out.mean;
  out.ok = true;
  return out;
}

void DrawAndSave(TH1 *h, const FitResult &fit, const std::string &outPath, const std::string &title) {
  TCanvas c("c", "c", 1100, 900);
  c.Divide(1, 2);

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.045);

  // Top pad: full range view
  c.cd(1);
  auto *hFull = dynamic_cast<TH1 *>(h->Clone("h_full_view"));
  hFull->SetTitle((title + " (full range)").c_str());
  hFull->SetLineWidth(2);
  hFull->Draw("hist");

  if (fit.ok) {
    auto *fgaus = new TF1("fgaus_draw", "gaus", fit.fitMin, fit.fitMax);
    fgaus->SetNpx(1000);
    fgaus->SetParameters(fit.amplitude, fit.mean, fit.sigma);
    fgaus->SetLineColor(kRed + 1);
    fgaus->SetLineWidth(4);
    fgaus->DrawCopy("LSAME");
    delete fgaus;

    TLine lmin(fit.fitMin, 0.0, fit.fitMin, hFull->GetMaximum() * 0.95);
    TLine lmax(fit.fitMax, 0.0, fit.fitMax, hFull->GetMaximum() * 0.95);
    lmin.SetLineColor(kBlue + 1);
    lmax.SetLineColor(kBlue + 1);
    lmin.SetLineStyle(2);
    lmax.SetLineStyle(2);
    lmin.SetLineWidth(2);
    lmax.SetLineWidth(2);
    lmin.Draw("same");
    lmax.Draw("same");

    latex.DrawLatex(0.62, 0.86, Form("Fit: gaus [%.0f, %.0f]", fit.fitMin, fit.fitMax));
    latex.DrawLatex(0.62, 0.79, Form("Mean = %.3f #pm %.3f", fit.mean, fit.meanErr));
    latex.DrawLatex(0.62, 0.72, Form("Sigma = %.3f #pm %.3f", fit.sigma, fit.sigmaErr));
    latex.DrawLatex(0.62, 0.65, Form("#chi^{2}/NDF = %.2f / %d", fit.chi2, fit.ndf));
    latex.DrawLatex(0.62, 0.58, Form("Calib const = %.9f", fit.constant));
  } else {
    latex.SetTextColor(kRed + 1);
    latex.DrawLatex(0.62, 0.86, "Fit: FAILED");
    latex.DrawLatex(0.62, 0.79, ("Reason: " + fit.reason).c_str());
    latex.SetTextColor(kBlack);
  }

  // Bottom pad: zoom around fit range
  c.cd(2);
  auto *hZoom = dynamic_cast<TH1 *>(h->Clone("h_zoom_view"));
  hZoom->SetTitle((title + " (zoomed fit window)").c_str());
  hZoom->SetLineWidth(2);
  if (fit.ok) {
    const double zoomMin = std::max(hZoom->GetXaxis()->GetXmin(), fit.mean - 4.0 * fit.sigma);
    const double zoomMax = std::min(hZoom->GetXaxis()->GetXmax(), fit.mean + 4.0 * fit.sigma);
    hZoom->GetXaxis()->SetRangeUser(zoomMin, zoomMax);
  }
  hZoom->Draw("hist");

  if (fit.ok) {
    auto *fgausZoom = new TF1("fgaus_draw_zoom", "gaus", fit.fitMin, fit.fitMax);
    fgausZoom->SetNpx(1000);
    fgausZoom->SetParameters(fit.amplitude, fit.mean, fit.sigma);
    fgausZoom->SetLineColor(kRed + 1);
    fgausZoom->SetLineWidth(4);
    fgausZoom->DrawCopy("LSAME");
    delete fgausZoom;
  }

  c.Modified();
  c.Update();
  c.SaveAs(outPath.c_str());

  delete hFull;
  delete hZoom;
}

}  // namespace

int main() {
  gROOT->SetBatch(kTRUE);
  gSystem->mkdir("./calib_cont", kTRUE);

  const std::vector<std::pair<int, int>> runTowerMap = {
      {14127, 7}, {14128, 8}, {14129, 9}, {14130, 6}, {14131, 5},
      {14132, 4}, {14133, 1}, {14134, 2}, {14135, 3}};

  std::map<int, FitResult> towerC;
  std::map<int, FitResult> towerS;

  std::ofstream csv("./calib_cont/calib_constants_summary.csv");
  csv << "run,tower,component,hist_name,entries,mean,mean_err,sigma,sigma_err,constant,status\n";

  for (const auto &[run, tower] : runTowerMap) {
    const std::string inPath = "./Calib/Calib_Run_" + std::to_string(run) + ".root";
    TFile inFile(inPath.c_str(), "READ");
    if (inFile.IsZombie()) {
      std::cerr << "[WARN] cannot open " << inPath << '\n';
      csv << run << ",T" << tower << ",C,,0,0,0,0,0,0,FILE_OPEN_FAIL\n";
      csv << run << ",T" << tower << ",S,,0,0,0,0,0,0,FILE_OPEN_FAIL\n";
      continue;
    }

    const std::string histCName = "T" + std::to_string(tower) + "C_intADC_after";
    const std::string histSName = "T" + std::to_string(tower) + "S_intADC_after";

    auto *hC = dynamic_cast<TH1 *>(inFile.Get(histCName.c_str()));
    auto *hS = dynamic_cast<TH1 *>(inFile.Get(histSName.c_str()));

    if (hC) hC->SetDirectory(nullptr);
    if (hS) hS->SetDirectory(nullptr);
    inFile.Close();

    FitResult fitC = FitHistogram(hC, "Run " + std::to_string(run) + " T" + std::to_string(tower) + "C");
    FitResult fitS = FitHistogram(hS, "Run " + std::to_string(run) + " T" + std::to_string(tower) + "S");

    if (hC) {
      DrawAndSave(hC, fitC, "./calib_cont/Fit_Run_" + std::to_string(run) + "_T" + std::to_string(tower) + "C.png",
                  "Run " + std::to_string(run) + " / T" + std::to_string(tower) + "C_intADC_after");
    }
    if (hS) {
      DrawAndSave(hS, fitS, "./calib_cont/Fit_Run_" + std::to_string(run) + "_T" + std::to_string(tower) + "S.png",
                  "Run " + std::to_string(run) + " / T" + std::to_string(tower) + "S_intADC_after");
    }

    const double entriesC = hC ? hC->GetEntries() : 0.0;
    const double entriesS = hS ? hS->GetEntries() : 0.0;
    csv << run << ",T" << tower << ",C," << histCName << "," << entriesC << ","
        << fitC.mean << "," << fitC.meanErr << "," << fitC.sigma << "," << fitC.sigmaErr << "," << fitC.constant
        << "," << (fitC.ok ? "OK" : fitC.reason) << "\n";
    csv << run << ",T" << tower << ",S," << histSName << "," << entriesS << ","
        << fitS.mean << "," << fitS.meanErr << "," << fitS.sigma << "," << fitS.sigmaErr << "," << fitS.constant
        << "," << (fitS.ok ? "OK" : fitS.reason) << "\n";

    towerC[tower] = fitC;
    towerS[tower] = fitS;

    delete hC;
    delete hS;
  }

  csv.close();

  std::ofstream txt("./calib_cont/calib_constants_for_calib_DRC.txt");
  txt << std::fixed << std::setprecision(9);
  txt << "// Paste below lines into calib_DRC.cc\n";
  for (int t = 1; t <= 9; ++t) {
    const auto it = towerC.find(t);
    const double val = (it != towerC.end() && it->second.ok) ? it->second.constant : -1.0;
    txt << "float calib_T" << t << "C = " << val << ";\n";
  }
  txt << "\n";
  for (int t = 1; t <= 9; ++t) {
    const auto it = towerS.find(t);
    const double val = (it != towerS.end() && it->second.ok) ? it->second.constant : -1.0;
    txt << "float calib_T" << t << "S = " << val << ";\n";
  }
  txt.close();

  std::cout << "[DONE] Wrote outputs to ./calib_cont\n";
  std::cout << "       - calib_constants_summary.csv\n";
  std::cout << "       - calib_constants_for_calib_DRC.txt\n";
  std::cout << "       - Fit_Run_<run>_T<tower>{C,S}.png\n";

  return 0;
}
