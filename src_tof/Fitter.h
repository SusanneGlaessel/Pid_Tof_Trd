/** @file   Fitter.h
    @class  Pid::Fitter
    @author Viktor Klochkov (klochkov44@gmail.com)
    @author Ilya Selyuzhenkov (ilya.selyuzhenkov@gmail.com)
    @brief  Class to fit 2D histograms
*/

#ifndef PidFitter_H
#define PidFitter_H 1

#include <utility>
#include <vector>

#include "TH1.h"
#include "TH2.h"
#include "TString.h"

#include "ParticleFit.h"

namespace Pid {

class Fitter {

 public:
  //     Fitter( std::vector <ParticleFit> &&particles ) : particles_(particles) {};
  Fitter() = default;
  ;

  void Fit();
  TF1* ConstructFit1DFunction(double p);
  double Fit1D(const std::shared_ptr<TH1>& h, std::vector<double>& par, std::vector<double>& par_err, double p);
  void Clear();

  void AddParticle(const ParticleFit& particle, uint id) {
    particles_.push_back(particle);
    particles_id_.push_back(id);
  }
  void SetHisto2D(std::shared_ptr<TH2> histo2D) { histo2D_ = std::move(histo2D); }
  void SetRangeX(double min, double max) { minx_ = min, maxx_ = max; }
  void SetRangeY(double min, double max) { miny_ = min, maxy_ = max; }
  void SetOutputFileName(TString name) { outfilename_ = std::move(name); }

  ParticleFit GetParticle(uint i) const { return particles_.at(i); };
  ParticleFit GetParticleSpecie(uint i) const {
    return particles_.at(std::find(particles_id_.begin(), particles_id_.end(), i) - particles_id_.begin());
  };

  void SetChi2Max(double chi2) { chi2_max_ = chi2; }

 private:
  std::vector<ParticleFit> particles_;
  std::vector<uint> particles_id_;
  std::shared_ptr<TH2> histo2D_{nullptr};

  TString outfilename_{"out.root"};

  double minx_{-1.};
  double maxx_{-1.};

  double miny_{-1.};
  double maxy_{-1.};

  double chi2_max_{100.};

  //     ClassDef(Fitter, 2);
};

}// namespace Pid
#endif// PidFitter_H
