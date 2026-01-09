#include "Fitter.h"
#include <iostream>

#include "TFile.h"
#include "TGraphErrors.h"
#include "TMath.h"

#include "Parameters.h"

// ClassImp(Pid::BaseFitterHelper);

namespace Pid {
/**
* Main function. Fitting TH2D bin-by-bin
*/
void Fitter::Fit() {
  const uint firstbin = histo2D_->GetXaxis()->FindBin(minx_);
  const uint lastbin = histo2D_->GetXaxis()->FindBin(maxx_);

  std::vector<std::vector<double>> params;
  std::vector<std::vector<double>> params_errors;
  std::vector<double> x;
  std::vector<double> chi2_x;
  std::vector<double> chi2_y;

  std::unique_ptr<TFile> f{TFile::Open(outfilename_, "recreate")};

  for (uint ibin = firstbin; ibin < lastbin; ++ibin) {
    std::shared_ptr<TH1> h1fit{histo2D_->ProjectionY(Form("py_%d", ibin), ibin, ibin)};

    std::vector<double> par;
    std::vector<double> par_err;

    const float mom = histo2D_->GetXaxis()->GetBinCenter(ibin);
    float chi2 = Fit1D(h1fit, par, par_err, mom);

    std::cout << mom << "  " << chi2 << std::endl;

    if (isnan(chi2) || isinf(chi2)) chi2 = -1.;

    chi2_x.push_back(mom);
    chi2_y.push_back(chi2);

    if (chi2 < 0. || chi2 > chi2_max_) continue;

    params.push_back(par);
    params_errors.push_back(par_err);
    x.push_back(mom);
  }

  Parameters p;
  p.SetParams(std::move(x), std::move(params), std::move(params_errors));
  //     p.SetParticles();
  TGraph chi2(chi2_x.size(), &(chi2_x[0]), &(chi2_y[0]));
  chi2.SetName("chi2");
  chi2.SetTitle("#chi^{2}/NDF;p (GeV/#it{c});#chi^{2}/NDF");
  chi2.Write();
  p.Parametrize(particles_);

  f->Close();
}

/**
* Constructs fit function as a sum of individual particle species. Parameters are also propagated
* @param h pointer to input histo
* @param par output: fit parameters 
* @param par_err output: fit parameters erorrs
* @param p track momentum
* @return chi2/NDF of the fit
*/
double Fitter::Fit1D(const std::shared_ptr<TH1>& h, std::vector<double>& par, std::vector<double>& par_err, double p) {
  auto f = ConstructFit1DFunction(p);//particles_.at(0).GetFunction(0.);

  h->Fit(f, "Q,M", "", miny_, maxy_);

  par = std::vector<double>(f->GetParameters(), f->GetParameters() + f->GetNpar());
  par_err = std::vector<double>(f->GetParErrors(), f->GetParErrors() + f->GetNpar());
  //
  //    int colors [6] = {kBlack,kGreen+2,kViolet,kOrange+2,kPink+2,kCyan+2};
  //    for (uint i = 0; i < particles_.size(); i++)
  //    {
  //      TF1 *func = (TF1*)particles_.at(i).GetFunction().Clone();
  //      func->SetLineColor(colors [i]);
  //      h->GetListOfFunctions()->Add(func);
  //    }
  //
  h->Write(Form("h_%f", p));

  return f->GetChisquare() / f->GetNDF();
}

/**
* Constructs fit function as a sum of individual particle species. Parameters are also propagated
* @param p track momentum
* @return pointer to TF1 function
*/
TF1* Fitter::ConstructFit1DFunction(double p) {
  TString sumname{""};
  std::vector<double> par{};

  uint iparticle{0};
  for (auto const& particle : particles_) {
    const TString name = particle.GetFunction().GetName();
    iparticle == 0 ? sumname = name : sumname += "+" + name;
    std::vector<double> par_i = particle.GetFunctionParams(p);

    par.insert(par.end(), par_i.begin(), par_i.end());
    iparticle++;
  }

  TF1* f;
  if (particles_.size() > 1)
    f = new TF1("fit1D", sumname, miny_, maxy_);
  else {
    f = new TF1;
    particles_.at(0).GetFunction().Copy(*f);
    //      f->SetName("fit1D");
    f->SetRange(miny_, maxy_);
  }
  f->SetParameters(&(par[0]));
  uint iparam_all{0};
  for (auto const& particle : particles_) {
    float pmin, pmax;
    bool notInRange = false;
    particle.GetRange(pmin, pmax);
    if (p < pmin || pmax < p) notInRange = true;
    for (uint iparam = 0; iparam < particle.GetNpar(); ++iparam, ++iparam_all) {
      if (notInRange) {
        f->FixParameter(iparam_all, 0.);
        continue;
      }

      if (particle.GetIsFixed(iparam)) {
        f->FixParameter(iparam_all, par.at(iparam_all));
        continue;
      }

      double parmin{-1.};
      double parmax{-1.};
      particle.GetFunction().GetParLimits(iparam, parmin, parmax);
      f->SetParLimits(iparam_all, parmin, parmax);
    }
  }
  //std::cout << sumname << std::endl;
  //std::cout << f->GetName() << " " << f->GetExpFormula() << std::endl;
  //for (auto ipar : par)
  //    std::cout << ipar << " ";
  //std::cout << std::endl;

  return f;
}
/**
* Clear everything
*/
void Fitter::Clear() {
  particles_.clear();
  particles_id_.clear();
  histo2D_.reset();
  minx_ = -1.;
  maxx_ = -1.;
  miny_ = -1.;
  maxy_ = -1.;
  chi2_max_ = 100.;
  outfilename_ = "out.root";
}

}// namespace Pid
