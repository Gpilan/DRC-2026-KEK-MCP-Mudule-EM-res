#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

namespace {

struct RunEnergy {
  int run;
  double energy;
};

struct FitStat {
  bool ok = false;
  double mean = 0.0;
  double meanErr = 0.0;
  double sigma = 0.0;
  double sigmaErr = 0.0;
  double entries = 0.0;
  std::string reason;
};

struct RunComponentPoint {
  int run = -1;
  double energy = -1.0;
  std::string comp;  // C/S/Comb
  FitStat fit;
  double linearity = 0.0;
  double linearityErr = 0.0;
  double resolution = 0.0;
  double resolutionErr = 0.0;
};

struct EnergyComponentPoint {
  double energy = -1.0;
  std::string comp;
  int nRuns = 0;
  double mean = 0.0;
  double meanErr = 0.0;
  double sigma = 0.0;
  double sigmaErr = 0.0;
  double linearity = 0.0;
  double linearityErr = 0.0;
  double resolution = 0.0;
  double resolutionErr = 0.0;
};

FitStat FitGaussian(TFile &f, const std::string &histName) {
  FitStat out;
  auto *h = dynamic_cast<TH1 *>(f.Get(histName.c_str()));
  if (!h) {
    out.reason = "missing";
    return out;
  }
  if (h->GetEntries() < 50) {
    out.reason = "low_entries";
    return out;
  }

  out.entries = h->GetEntries();

  const int peakBin = h->GetMaximumBin();
  const double peakX = h->GetXaxis()->GetBinCenter(peakBin);
  const double peakY = h->GetBinContent(peakBin);
  const double xmin = h->GetXaxis()->GetXmin();
  const double xmax = h->GetXaxis()->GetXmax();

  double sigma0 = h->GetRMS();
  if (sigma0 <= 0.0) sigma0 = (xmax - xmin) / 12.0;

  double fitMin = std::max(xmin, peakX - 1.5 * sigma0);
  double fitMax = std::min(xmax, peakX + 1.5 * sigma0);
  if (fitMax <= fitMin) {
    fitMin = xmin;
    fitMax = xmax;
  }

  TF1 fgaus("fgaus_tmp", "gaus", fitMin, fitMax);
  fgaus.SetParameters(peakY, peakX, sigma0);
  int status = h->Fit(&fgaus, "RQ0");
  if (status != 0) {
    fitMin = std::max(xmin, peakX - 2.5 * sigma0);
    fitMax = std::min(xmax, peakX + 2.5 * sigma0);
    fgaus.SetRange(fitMin, fitMax);
    fgaus.SetParameters(peakY, peakX, sigma0);
    status = h->Fit(&fgaus, "RQ0");
  }

  if (status != 0) {
    out.reason = "gaus_fit_failed";
    return out;
  }

  out.ok = true;
  out.mean = fgaus.GetParameter(1);
  out.meanErr = fgaus.GetParError(1);
  out.sigma = fgaus.GetParameter(2);
  out.sigmaErr = fgaus.GetParError(2);
  return out;
}

void ScaleFitStat(FitStat &fit, double scale) {
  if (!fit.ok) return;
  fit.mean *= scale;
  fit.meanErr *= scale;
  fit.sigma *= scale;
  fit.sigmaErr *= scale;
}

void ComputeDerived(RunComponentPoint &p) {
  if (!p.fit.ok || p.energy <= 0.0 || p.fit.mean <= 0.0 || p.fit.sigma <= 0.0) return;
  p.linearity = p.fit.mean / p.energy;
  p.linearityErr = p.fit.meanErr / p.energy;

  p.resolution = p.fit.sigma / p.fit.mean;
  const double relSigmaErr = p.fit.sigmaErr / p.fit.sigma;
  const double relMeanErr = p.fit.meanErr / p.fit.mean;
  p.resolutionErr = p.resolution * std::sqrt(relSigmaErr * relSigmaErr + relMeanErr * relMeanErr);
}

EnergyComponentPoint WeightedMerge(const std::vector<RunComponentPoint> &pts, double energy, const std::string &comp) {
  EnergyComponentPoint out;
  out.energy = energy;
  out.comp = comp;

  double swMean = 0.0, swxMean = 0.0;
  double swSig = 0.0, swxSig = 0.0;
  for (const auto &p : pts) {
    if (!p.fit.ok || p.fit.meanErr <= 0.0 || p.fit.sigmaErr <= 0.0) continue;
    out.nRuns++;
    const double wMean = 1.0 / (p.fit.meanErr * p.fit.meanErr);
    const double wSig = 1.0 / (p.fit.sigmaErr * p.fit.sigmaErr);
    swMean += wMean;
    swxMean += wMean * p.fit.mean;
    swSig += wSig;
    swxSig += wSig * p.fit.sigma;
  }
  if (out.nRuns == 0 || swMean <= 0.0 || swSig <= 0.0) return out;

  out.mean = swxMean / swMean;
  out.meanErr = std::sqrt(1.0 / swMean);
  out.sigma = swxSig / swSig;
  out.sigmaErr = std::sqrt(1.0 / swSig);

  out.linearity = out.mean / out.energy;
  out.linearityErr = out.meanErr / out.energy;

  out.resolution = out.sigma / out.mean;
  const double relSigmaErr = out.sigmaErr / out.sigma;
  const double relMeanErr = out.meanErr / out.mean;
  out.resolutionErr = out.resolution * std::sqrt(relSigmaErr * relSigmaErr + relMeanErr * relMeanErr);
  return out;
}

constexpr int kColorC = kAzure + 2;
constexpr int kColorS = kRed;
constexpr int kColorSum = kBlack;
constexpr int kMarkerStyle = 20;  // circle
constexpr int kFitLineStyle = 2;    // dashed

void StyleGraph(TGraphErrors &g, int color) {
  g.SetLineColor(color);
  g.SetMarkerColor(color);
  g.SetMarkerStyle(kMarkerStyle);
  g.SetMarkerSize(1.0);
  g.SetLineWidth(2);
}

struct FitResult {
  bool ok = false;
  double p0 = 0.0;
  double p1 = 0.0;
  double p2 = 0.0;
  double chi2 = 0.0;
  int ndf = 0;
  int npar = 0;
};

FitResult FitAndDraw(TGraphErrors &g, const std::string &fitName, const char *formula, double xmin,
                     double xmax, int color) {
  FitResult out;
  if (g.GetN() < 2) return out;

  TF1 f(fitName.c_str(), formula, xmin, xmax);
  f.SetLineColor(color);
  f.SetLineStyle(kFitLineStyle);
  f.SetLineWidth(2);
  f.SetNpx(1000);
  const int status = g.Fit(&f, "Q0");
  if (status != 0) return out;

  out.ok = true;
  out.npar = f.GetNpar();
  out.p0 = f.GetParameter(0);
  if (out.npar > 1) out.p1 = f.GetParameter(1);
  if (out.npar > 2) out.p2 = f.GetParameter(2);
  out.chi2 = f.GetChisquare();
  out.ndf = f.GetNDF();

  TF1 *fDraw = dynamic_cast<TF1 *>(f.DrawCopy("LSAME"));
  if (fDraw) {
    fDraw->SetLineColor(color);
    fDraw->SetLineStyle(kFitLineStyle);
    fDraw->SetLineWidth(2);
    fDraw->SetNpx(1000);
  }
  return out;
}

// sigma/E = sqrt((p0/E)^2 + (p1/sqrt(E))^2 + p2^2)
FitResult FitResoEQuadratureAndDraw(TGraphErrors &g, const std::string &fitName, double xmin,
                                    double xmax, int color) {
  FitResult out;
  if (g.GetN() < 2) return out;

  TF1 f(fitName.c_str(), "TMath::Sqrt([0]*[0]/(x*x) + [1]*[1]/x + [2]*[2])", xmin, xmax);
  f.SetLineColor(color);
  f.SetLineStyle(kFitLineStyle);
  f.SetLineWidth(2);
  f.SetNpx(1000);
  f.SetParameter(0, 0.01);
  f.SetParameter(1, 0.15);
  f.SetParameter(2, 0.02);
  const int status = g.Fit(&f, "Q0");
  if (status != 0) return out;

  out.ok = true;
  out.npar = f.GetNpar();
  out.p0 = f.GetParameter(0);
  if (out.npar > 1) out.p1 = f.GetParameter(1);
  if (out.npar > 2) out.p2 = f.GetParameter(2);
  out.chi2 = f.GetChisquare();
  out.ndf = f.GetNDF();

  TF1 *fDraw = dynamic_cast<TF1 *>(f.DrawCopy("LSAME"));
  if (fDraw) {
    fDraw->SetLineColor(color);
    fDraw->SetLineStyle(kFitLineStyle);
    fDraw->SetLineWidth(2);
    fDraw->SetNpx(1000);
  }
  return out;
}

std::string LinearLegendText(const char *label, const FitResult &fit) {
  if (!fit.ok) return "";
  if (fit.npar == 2) {
    return Form("%s: %.3fx%+.3f (#chi^{2}/NDF = %.2f/%d)", label, fit.p0, fit.p1, fit.chi2, fit.ndf);
  }
  return Form("%s: %.3fx (#chi^{2}/NDF = %.2f/%d)", label, fit.p0, fit.chi2, fit.ndf);
}

std::string ResoLegendText(const char *label, const FitResult &fit) {
  if (!fit.ok) return "";
  return Form("%s: (%.3f/E) #oplus (%.3f/#sqrt{E}) #oplus %.3f (#chi^{2}/NDF = %.2f/%d)", label,
              fit.p0, fit.p1, fit.p2, fit.chi2, fit.ndf);
}

std::string InvSqrtLegendText(const char *label, const FitResult &fit) {
  if (!fit.ok) return "";
  return Form("%s: %+.4f/#sqrt{E} %+.4f (#chi^{2}/NDF = %.2f/%d)", label, fit.p0, fit.p1, fit.chi2,
              fit.ndf);
}

void AddColoredLegendEntry(TLegend &leg, TObject *obj, const char *text, const char *opt, int color) {
  leg.AddEntry(obj, text, opt);
  TList *entries = leg.GetListOfPrimitives();
  if (!entries || entries->IsEmpty()) return;
  if (auto *entry = dynamic_cast<TLegendEntry *>(entries->Last())) {
    entry->SetTextColor(color);
  }
}

void StyleCanvasFrame(TCanvas &c) {
  c.SetGrid();
  c.SetGridx();
  c.SetGridy();
}

void FinalizeCanvas(TCanvas &c) {
  c.Modified();
  c.Update();
}

TGraphErrors BuildGraph(const std::vector<EnergyComponentPoint> &pts, const std::string &comp,
                        const std::string &xKey, const std::string &yKey) {
  TGraphErrors g;
  int ip = 0;
  for (const auto &p : pts) {
    if (p.comp != comp || p.nRuns == 0) continue;
    double x = 0.0, y = 0.0, ex = 0.0, ey = 0.0;

    if (xKey == "E") x = p.energy;
    if (xKey == "invSqrtE") x = 1.0 / std::sqrt(p.energy);

    if (yKey == "Ereco") { y = p.mean; ey = p.meanErr; }
    if (yKey == "linearity") { y = p.linearity; ey = p.linearityErr; }
    if (yKey == "resolution") { y = p.resolution; ey = p.resolutionErr; }

    g.SetPoint(ip, x, y);
    g.SetPointError(ip, ex, ey);
    ip++;
  }
  return g;
}

std::vector<EnergyComponentPoint> MergeByEnergy(
    const std::vector<RunComponentPoint> &runPoints,
    const std::vector<std::pair<std::string, std::string>> &compHist) {
  std::vector<EnergyComponentPoint> ePoints;
  for (const auto &[comp, _] : compHist) {
    std::map<double, std::vector<RunComponentPoint>> grouped;
    for (const auto &p : runPoints) {
      if (p.comp == comp) grouped[p.energy].push_back(p);
    }
    for (const auto &[e, pts] : grouped) {
      ePoints.push_back(WeightedMerge(pts, e, comp));
    }
  }
  std::sort(ePoints.begin(), ePoints.end(),
            [](const EnergyComponentPoint &a, const EnergyComponentPoint &b) {
              if (a.comp == b.comp) return a.energy < b.energy;
              return a.comp < b.comp;
            });
  return ePoints;
}

void WritePerRunCsv(const std::vector<RunComponentPoint> &runPoints, const std::string &path) {
  std::ofstream out(path);
  out << "beamE,run,component,mean,meanErr,sigma,sigmaErr,linearity,linearityErr,resolution,"
         "resolutionErr,status\n";
  out << std::fixed << std::setprecision(6);
  for (const auto &p : runPoints) {
    out << p.energy << "," << p.run << "," << p.comp << ","
        << p.fit.mean << "," << p.fit.meanErr << ","
        << p.fit.sigma << "," << p.fit.sigmaErr << ","
        << p.linearity << "," << p.linearityErr << ","
        << p.resolution << "," << p.resolutionErr << ","
        << (p.fit.ok ? "OK" : p.fit.reason) << "\n";
  }
}

void WriteByEnergyCsv(const std::vector<EnergyComponentPoint> &ePoints, const std::string &path) {
  std::ofstream out(path);
  out << "beamE,component,nRuns,mean,meanErr,sigma,sigmaErr,linearity,linearityErr,resolution,"
         "resolutionErr\n";
  out << std::fixed << std::setprecision(6);
  for (const auto &p : ePoints) {
    out << p.energy << "," << p.comp << "," << p.nRuns << ","
        << p.mean << "," << p.meanErr << ","
        << p.sigma << "," << p.sigmaErr << ","
        << p.linearity << "," << p.linearityErr << ","
        << p.resolution << "," << p.resolutionErr << "\n";
  }
}

void DrawEnergyScanPlots(const std::vector<EnergyComponentPoint> &ePoints,
                         const std::string &tag) {
  const std::string outDir = "./energy_scan/";
  const std::string fitTag = tag.empty() ? "" : tag;
  constexpr double kXMinE = 0.0;
  constexpr double kXMaxE = 6.0;
  constexpr double kYMaxEreco = 6.0;
  constexpr double kYMaxReso = 0.5;

  TGraphErrors gE_C = BuildGraph(ePoints, "C", "E", "Ereco");
  TGraphErrors gE_S = BuildGraph(ePoints, "S", "E", "Ereco");
  TGraphErrors gE_Sum = BuildGraph(ePoints, "Sum", "E", "Ereco");
  TGraphErrors gLin_C = BuildGraph(ePoints, "C", "E", "linearity");
  TGraphErrors gLin_S = BuildGraph(ePoints, "S", "E", "linearity");
  TGraphErrors gLin_Sum = BuildGraph(ePoints, "Sum", "E", "linearity");
  TGraphErrors gResE_C = BuildGraph(ePoints, "C", "E", "resolution");
  TGraphErrors gResE_S = BuildGraph(ePoints, "S", "E", "resolution");
  TGraphErrors gResE_Sum = BuildGraph(ePoints, "Sum", "E", "resolution");
  TGraphErrors gResInv_C = BuildGraph(ePoints, "C", "invSqrtE", "resolution");
  TGraphErrors gResInv_S = BuildGraph(ePoints, "S", "invSqrtE", "resolution");
  TGraphErrors gResInv_Sum = BuildGraph(ePoints, "Sum", "invSqrtE", "resolution");

  StyleGraph(gE_C, kColorC);
  StyleGraph(gE_S, kColorS);
  StyleGraph(gE_Sum, kColorSum);
  StyleGraph(gLin_C, kColorC);
  StyleGraph(gLin_S, kColorS);
  StyleGraph(gLin_Sum, kColorSum);
  StyleGraph(gResE_C, kColorC);
  StyleGraph(gResE_S, kColorS);
  StyleGraph(gResE_Sum, kColorSum);
  StyleGraph(gResInv_C, kColorC);
  StyleGraph(gResInv_S, kColorS);
  StyleGraph(gResInv_Sum, kColorSum);

  constexpr double kEminFit = 0.4;
  constexpr double kEmaxFit = 6.2;
  constexpr double kInvSqrtMin = 0.35;
  constexpr double kInvSqrtMax = 1.55;

  TCanvas c1(("c1" + tag).c_str(), "Ereco vs Ebeam", 1000, 700);
  StyleCanvasFrame(c1);
  TMultiGraph mg1;
  mg1.Add(&gE_C, "P");
  mg1.Add(&gE_S, "P");
  mg1.Add(&gE_Sum, "P");
  mg1.Draw("A");
  mg1.GetXaxis()->SetTitle("E [GeV]");
  mg1.GetYaxis()->SetTitle("E_{reco} [GeV]");
  mg1.SetTitle("");
  const FitResult fitE_C =
      FitAndDraw(gE_C, ("fEreco_C" + fitTag).c_str(), "[0]*x+[1]", kEminFit, kEmaxFit, kColorC);
  const FitResult fitE_S =
      FitAndDraw(gE_S, ("fEreco_S" + fitTag).c_str(), "[0]*x+[1]", kEminFit, kEmaxFit, kColorS);
  const FitResult fitE_Sum =
      FitAndDraw(gE_Sum, ("fEreco_Sum" + fitTag).c_str(), "[0]*x+[1]", kEminFit, kEmaxFit, kColorSum);
  mg1.GetXaxis()->SetRangeUser(kXMinE, kXMaxE);
  mg1.GetYaxis()->SetRangeUser(0.0, kYMaxEreco);
  c1.Update();
  TLegend leg1(0.12, 0.60, 0.53, 0.89);
  leg1.SetTextSize(0.030);
  leg1.SetBorderSize(0);
  leg1.SetFillStyle(0);
  if (fitE_C.ok)
    AddColoredLegendEntry(leg1, &gE_C, LinearLegendText("C", fitE_C).c_str(), "lp", kColorC);
  if (fitE_S.ok)
    AddColoredLegendEntry(leg1, &gE_S, LinearLegendText("S", fitE_S).c_str(), "lp", kColorS);
  if (fitE_Sum.ok)
    AddColoredLegendEntry(leg1, &gE_Sum, LinearLegendText("Sum", fitE_Sum).c_str(), "lp", kColorSum);
  leg1.Draw();
  FinalizeCanvas(c1);
  c1.SaveAs((outDir + "linearity_Ereco_vs_Ebeam_overlay" + tag + ".png").c_str());

  TCanvas c2(("c2" + tag).c_str(), "Linearity ratio vs Ebeam", 1000, 700);
  StyleCanvasFrame(c2);
  TMultiGraph mg2;
  mg2.Add(&gLin_C, "P");
  mg2.Add(&gLin_S, "P");
  mg2.Add(&gLin_Sum, "P");
  mg2.Draw("A");
  mg2.GetXaxis()->SetTitle("E [GeV]");
  mg2.GetXaxis()->SetLimits(kXMinE, kXMaxE);
  mg2.GetYaxis()->SetTitle("E_{reco}/E_{beam}");
  mg2.GetYaxis()->SetRangeUser(0.8, 1.1);
  mg2.GetYaxis()->SetNdivisions(306, kFALSE);
  mg2.GetYaxis()->SetDecimals(kTRUE);
  gStyle->SetStripDecimals(kFALSE);
  mg2.SetTitle("");
  c2.Update();
  TLine lineUnity(gPad->GetUxmin(), 1.0, gPad->GetUxmax(), 1.0);
  lineUnity.SetLineColor(kGreen + 1);
  lineUnity.SetLineStyle(kFitLineStyle);
  lineUnity.SetLineWidth(2);
  lineUnity.Draw("SAME");
  TLegend leg2(0.12, 0.78, 0.30, 0.89);
  leg2.SetTextSize(0.035);
  leg2.SetBorderSize(0);
  leg2.SetFillStyle(0);
  AddColoredLegendEntry(leg2, &gLin_C, "C", "p", kColorC);
  AddColoredLegendEntry(leg2, &gLin_S, "S", "p", kColorS);
  AddColoredLegendEntry(leg2, &gLin_Sum, "Sum", "p", kColorSum);
  leg2.Draw();
  FinalizeCanvas(c2);
  c2.SaveAs((outDir + "linearity_ratio_vs_Ebeam_overlay" + tag + ".png").c_str());

  TCanvas c3(("c3" + tag).c_str(), "Resolution vs Ebeam", 1000, 700);
  StyleCanvasFrame(c3);
  TMultiGraph mg3;
  mg3.Add(&gResE_C, "P");
  mg3.Add(&gResE_S, "P");
  mg3.Add(&gResE_Sum, "P");
  mg3.Draw("A");
  mg3.GetXaxis()->SetTitle("E [GeV]");
  mg3.GetXaxis()->SetLimits(kXMinE, kXMaxE);
  mg3.GetYaxis()->SetTitle("#sigma/E");
  mg3.GetYaxis()->SetRangeUser(0.0, kYMaxReso);
  mg3.SetTitle("");
  const FitResult fitResE_C = FitResoEQuadratureAndDraw(gResE_C, ("fResE_C" + fitTag).c_str(),
                                                         kEminFit, kEmaxFit, kColorC);
  const FitResult fitResE_S = FitResoEQuadratureAndDraw(gResE_S, ("fResE_S" + fitTag).c_str(),
                                                         kEminFit, kEmaxFit, kColorS);
  const FitResult fitResE_Sum = FitResoEQuadratureAndDraw(gResE_Sum, ("fResE_Sum" + fitTag).c_str(),
                                                          kEminFit, kEmaxFit, kColorSum);
  TLegend leg3(0.25, 0.58, 0.75, 0.89);
  leg3.SetTextSize(0.030);
  leg3.SetBorderSize(0);
  leg3.SetFillStyle(0);
  if (fitResE_C.ok)
    AddColoredLegendEntry(leg3, &gResE_C, ResoLegendText("C", fitResE_C).c_str(), "lp", kColorC);
  if (fitResE_S.ok)
    AddColoredLegendEntry(leg3, &gResE_S, ResoLegendText("S", fitResE_S).c_str(), "lp", kColorS);
  if (fitResE_Sum.ok)
    AddColoredLegendEntry(leg3, &gResE_Sum, ResoLegendText("Sum", fitResE_Sum).c_str(), "lp", kColorSum);
  leg3.Draw();
  FinalizeCanvas(c3);
  c3.SaveAs((outDir + "resolution_vs_Ebeam_overlay" + tag + ".png").c_str());

  TCanvas c4(("c4" + tag).c_str(), "Resolution vs invSqrtE", 1000, 700);
  StyleCanvasFrame(c4);
  TMultiGraph mg4;
  mg4.Add(&gResInv_C, "P");
  mg4.Add(&gResInv_S, "P");
  mg4.Add(&gResInv_Sum, "P");
  mg4.Draw("A");
  mg4.GetXaxis()->SetTitle("1/#sqrt{E} [GeV^{-1/2}]");
  mg4.GetXaxis()->SetLimits(kInvSqrtMin, kInvSqrtMax);
  mg4.GetYaxis()->SetTitle("#sigma/E");
  mg4.GetYaxis()->SetRangeUser(0.0, kYMaxReso);
  mg4.SetTitle("");
  const FitResult fitResInv_C = FitAndDraw(gResInv_C, ("fResInv_C" + fitTag).c_str(), "[0]*x+[1]",
                                           kInvSqrtMin, kInvSqrtMax, kColorC);
  const FitResult fitResInv_S = FitAndDraw(gResInv_S, ("fResInv_S" + fitTag).c_str(), "[0]*x+[1]",
                                           kInvSqrtMin, kInvSqrtMax, kColorS);
  const FitResult fitResInv_Sum = FitAndDraw(gResInv_Sum, ("fResInv_Sum" + fitTag).c_str(),
                                             "[0]*x+[1]", kInvSqrtMin, kInvSqrtMax, kColorSum);
  TLegend leg4(0.02, 0.52, 0.43, 0.89);
  leg4.SetTextSize(0.030);
  leg4.SetBorderSize(0);
  leg4.SetFillStyle(0);
  if (fitResInv_C.ok)
    AddColoredLegendEntry(leg4, nullptr, InvSqrtLegendText("C", fitResInv_C).c_str(), "", kColorC);
  if (fitResInv_S.ok)
    AddColoredLegendEntry(leg4, nullptr, InvSqrtLegendText("S", fitResInv_S).c_str(), "", kColorS);
  if (fitResInv_Sum.ok)
    AddColoredLegendEntry(leg4, nullptr, InvSqrtLegendText("Sum", fitResInv_Sum).c_str(), "",
                          kColorSum);
  leg4.Draw();
  FinalizeCanvas(c4);
  c4.SaveAs((outDir + "resolution_vs_invSqrtE_overlay_fit" + tag + ".png").c_str());
}

std::vector<RunComponentPoint> FilterByMaxEnergy(const std::vector<RunComponentPoint> &runPoints,
                                                 double maxEnergy) {
  std::vector<RunComponentPoint> out;
  out.reserve(runPoints.size());
  for (const auto &p : runPoints) {
    if (p.energy <= maxEnergy) out.push_back(p);
  }
  return out;
}

}  // namespace

int main() {
  gROOT->SetBatch(kTRUE);
  gSystem->mkdir("./energy_scan", kTRUE);

  const std::vector<RunEnergy> mapping = {
      {14222, 3.0}, {14223, 3.5}, {14224, 4.0}, {14226, 4.0}, {14227, 1.5},
      {14228, 2.0}, {14229, 2.5}, {14230, 1.0}, {14231, 4.5}, {14232, 5.0},
      {14233, 5.5}, {14234, 5.8}, {14235, 0.5}, {14236, 0.5}, {14237, 2.0},
      {14238, 1.5}, {14239, 1.0}, {14240, 0.5}};

  const std::vector<std::pair<std::string, std::string>> compHist = {
      {"C", "totalEdepScaled_C_after"},
      {"S", "totalEdepScaled_S_after"},
      // Comb histogram stores C+S per event; Summation uses (C+S)/2.
      {"Sum", "totalEdepScaled_Comb_after"},
  };

  std::vector<RunComponentPoint> runPoints;
  for (const auto &re : mapping) {
    const std::string inPath = "./Calib/Calib_Run_" + std::to_string(re.run) + ".root";
    TFile inFile(inPath.c_str(), "READ");
    if (inFile.IsZombie()) {
      std::cerr << "[WARN] cannot open " << inPath << '\n';
      continue;
    }

    for (const auto &[comp, histName] : compHist) {
      RunComponentPoint p;
      p.run = re.run;
      p.energy = re.energy;
      p.comp = comp;
      p.fit = FitGaussian(inFile, histName);
      if (comp == "Sum") ScaleFitStat(p.fit, 0.5);
      ComputeDerived(p);
      runPoints.push_back(p);
    }
    inFile.Close();
  }

  const std::vector<EnergyComponentPoint> ePointsAll = MergeByEnergy(runPoints, compHist);
  WritePerRunCsv(runPoints, "./energy_scan/linearity_resolution_per_run.csv");
  WriteByEnergyCsv(ePointsAll, "./energy_scan/linearity_resolution_by_energy.csv");
  DrawEnergyScanPlots(ePointsAll, "");

  const std::vector<RunComponentPoint> runPointsLe5 = FilterByMaxEnergy(runPoints, 5.0);
  const std::vector<EnergyComponentPoint> ePointsLe5 = MergeByEnergy(runPointsLe5, compHist);
  WritePerRunCsv(runPointsLe5, "./energy_scan/linearity_resolution_per_run_le5GeV.csv");
  WriteByEnergyCsv(ePointsLe5, "./energy_scan/linearity_resolution_by_energy_le5GeV.csv");
  DrawEnergyScanPlots(ePointsLe5, "_le5GeV");

  std::cout << "[DONE] Prepared linearity/resolution outputs in ./energy_scan\n"
            << "       - linearity_resolution_per_run.csv\n"
            << "       - linearity_resolution_by_energy.csv\n"
            << "       - linearity_Ereco_vs_Ebeam_overlay.png\n"
            << "       - linearity_ratio_vs_Ebeam_overlay.png\n"
            << "       - resolution_vs_Ebeam_overlay.png\n"
            << "       - resolution_vs_invSqrtE_overlay_fit.png\n"
            << "       - linearity_resolution_per_run_le5GeV.csv\n"
            << "       - linearity_resolution_by_energy_le5GeV.csv\n"
            << "       - linearity_Ereco_vs_Ebeam_overlay_le5GeV.png\n"
            << "       - linearity_ratio_vs_Ebeam_overlay_le5GeV.png\n"
            << "       - resolution_vs_Ebeam_overlay_le5GeV.png\n"
            << "       - resolution_vs_invSqrtE_overlay_fit_le5GeV.png\n";
  return 0;
}
