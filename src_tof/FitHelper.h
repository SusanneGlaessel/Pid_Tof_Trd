#ifndef PID_SRC_FITHELPER_H_
#define PID_SRC_FITHELPER_H_

#include "Fitter.h"
#include "ParticleFit.h"
#include <TCutG.h>

TH2* cutTH2(std::shared_ptr<TH2D> hFull, TCutG* cut) {
  std::string name = hFull->GetName();
  name += "_cut";
  TH2* hCut = (TH2*) hFull->Clone(name.c_str());
  int nBinsX = hCut->GetNbinsX();
  int nBinsY = hCut->GetNbinsY();
  double x, y;

  for (int i = 1; i <= nBinsX; i++) {
    for (int j = 1; j <= nBinsY; j++) {
      x = hCut->GetXaxis()->GetBinCenter(i);
      y = hCut->GetYaxis()->GetBinCenter(j);
      if (!cut->IsInside(x, y)) {
        hCut->SetBinContent(i, j, 0.);
        hCut->SetBinError(i, j, 0.);
      }
    }
  }
  return hCut;
}

double asymmGauss(double x, std::array<double, 4> p) {
  double amp = p[0];
  double mean = p[1];
  double sig = p[2];
  double delta = p[3];
  if (x > mean) {
    sig *= 1 + delta;
  } else {
    sig *= 1 - delta;
  }
  return amp * exp(-0.5 * (x - mean) * (x - mean) / sig / sig);
}

struct FitParams {
  std::string name{""};
  float xmin{0}, xmax{0}, ymin{0}, ymax{0};
  float chi2_max{0};
  unsigned int pid;
  std::vector<std::string> func{"0", "pol9", "pol9"};
};

Pid::ParticleFit FitParticle(Pid::Fitter& tof, const FitParams& params, const std::shared_ptr<TH2>& hist) {

  std::string name = params.name;
  std::cout << "\n + " + name + "\n";

  tof.SetChi2Max(params.chi2_max);
  TF1 fit(("fit_" + name).c_str(), "gaus", params.ymin, params.ymax);
  fit.SetParNames("p0", "p1", "p2");
  fit.SetParLimits(0, 0., 4.e6);
  fit.SetParLimits(1, -.25, 0.03);
  fit.SetParLimits(2, 0., 0.7);

  TF1 fit_0((name + "_0").c_str(), params.func.at(0).c_str(), params.xmin, params.xmax);
  TF1 fit_1((name + "_1").c_str(), params.func.at(1).c_str(), params.xmin, params.xmax);
  TF1 fit_2((name + "_2").c_str(), params.func.at(2).c_str(), params.xmin, params.xmax);

  Pid::ParticleFit particle(params.pid);
  particle.SetParametrization({fit_0, fit_1, fit_2});
  particle.SetFitFunction(fit);
  particle.SetRange(params.xmin, params.xmax);
  particle.SetIsFitted();

  tof.AddParticle(particle, params.pid);
  tof.SetHisto2D(hist);
  tof.SetRangeX(params.xmin, params.xmax);
  tof.SetRangeY(params.ymin, params.ymax);
  tof.SetOutputFileName("pionpos.root");
  tof.Fit();
  particle = tof.GetParticle(0);
  tof.Clear();

  return particle;
}

Pid::ParticleFit FitPionsPos(Pid::Fitter& tof, const std::shared_ptr<TH2D>& hist) {

  FitParams pions;
  pions.xmin = 0.25;
  pions.xmax = 10.5;
  pions.ymin = -2.;
  pions.ymax = 2.;
  pions.pid = PidTofParticles::kPionPos;
  pions.name = "piplus";

  return FitParticle(tof, pions, hist);
}

#endif//PID_SRC_FITHELPER_H_
