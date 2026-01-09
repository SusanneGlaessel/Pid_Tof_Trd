#include "Parameters.h"

#include "TFile.h"
#include "TGraphErrors.h"

// ClassImp(Pid::Parameters)

namespace Pid {

// -----   Default constructor   -------------------------------------------
Parameters::Parameters() = default;

/**
* Parametrize fit results for all particles
* @param particles vector of particles
*/
void Parameters::Parametrize(std::vector<ParticleFit>& particles) {
  std::vector<TGraphErrors> graphs;

  for (uint ipart = 0; ipart < particles.size(); ++ipart) {
    const uint nvar = particles.at(ipart).GetNpar();
    for (uint ivar = 0; ivar < nvar; ++ivar) {
      std::vector<double> y;
      std::vector<double> dy;

      for (uint ibin = 0; ibin < params_.size(); ++ibin) {
        const double par = params_.at(ibin).at(nvar * ipart + ivar);
        const double par_err = params_errors_.at(ibin).at(nvar * ipart + ivar);

        //              std::cout << par << " " <<  par_err << " " << std::endl;
        y.push_back(par);
        dy.push_back(par_err);
      }

      TF1& fit = particles.at(ipart).GetParametrizationFunction(ivar);
      auto* par = new TGraphErrors(params_.size(), &(x_[0]), &(y[0]), nullptr, &(dy[0]));
      par->SetName(fit.GetName());
      par->SetTitle(Form("%s;p (GeV/#it{c});%s", fit.GetName(), fit.GetName()));

      if (strcmp(fit.GetTitle(), "0") == 0) {
        std::function<double(const double*, const double*)> fitFunc = [=](const double* x, const double*) { return par->Eval(x[0]); };
        TF1 false_fit(fit.GetName(), fitFunc, fit.GetXmin(), fit.GetXmax(), 0);
        false_fit.SetNpx(par->GetN() * 5);
        false_fit.SetTitle("0");
        particles.at(ipart).SetParametrizationFunction(ivar, false_fit);
        par->GetListOfFunctions()->Add(false_fit.Clone());
      } else {
        if (particles.size() > 1) {
          float pmin, pmax;
          particles.at(ipart).GetRange(pmin, pmax);
          par->Fit(&fit, "Q,M", "", pmin, pmax);
        } else
          par->Fit(&fit, "Q,M,R");
      }
      par->Write();
      graphs.push_back(*par);
    }
  }
}

}// namespace Pid
