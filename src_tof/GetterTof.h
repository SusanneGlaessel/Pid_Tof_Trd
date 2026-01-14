/** @file   GetterTof.h
    @class  Pid::GetterTof
    @author Viktor Klochkov (klochkov44@gmail.com)
    @author Ilya Selyuzhenkov (ilya.selyuzhenkov@gmail.com)
    @brief  Class to calculate Tof-PID probabilities
*/

#ifndef PidGetterTof_H
#define PidGetterTof_H 1

#include "ParticleFit.h"
#include "TObject.h"
#include <TCutG.h>
#include <TMultiGraph.h>
#include <TText.h>
#include <vector>

namespace Pid {

/**
 * @brief interface class for Pid getters
 */
class BaseGetter {

 public:
  virtual ~BaseGetter() = default;

  virtual double GetWeight(double var1, double var2, int pid) = 0;

  virtual std::map<int, double> GetWeights(double var1, double var2) = 0;

  virtual int GetPid(double var1, double var2, double purity) = 0;

  virtual void Streamer(TBuffer&){};
};

/**
 * @brief Bayesian Pid getter
 */
class GetterTof : public TObject, public BaseGetter {
 public:
  void AddParticle(const ParticleFit& particle, uint id) { species_[id] = particle; }
  void AddParticles(std::map<int, ParticleFit>&& species) { species_ = species; }

  std::map<int, double> GetBayesianProbability(double p, double m2);
  void SetRange(double min, double max) { minx_ = min, maxx_ = max; }

  std::map<uint, double> GetSigma(double p, double m2) {
    std::map<uint, double> sigma{};

    if (p > maxx_ || p < minx_)
      return sigma;

    for (auto& specie : species_) {
      sigma[specie.first] = abs(m2 - specie.second.GetMean(p)) / specie.second.GetSigma(p);
    }
    return sigma;
  }

  int GetPid(double p, double m2, double purity) override {
    auto prob = GetBayesianProbability(p, m2);
    for (auto& pid_prob : prob)
      if (pid_prob.second >= purity) return pid_prob.first;
    return 1;
  }

  double GetWeight(double var1, double var2, int pid) override {
    // not yet implemented
    return 1.0;
  }
  std::map<int, double> GetWeights(double var1, double var2) override {
    return GetBayesianProbability(var1, var2);
  }

  const ParticleFit& GetParticleFit(int pid) {
    auto it = species_.find(pid);
    if (it != species_.end()) {
      return it->second;
    }
    throw std::runtime_error("Particle " + std::to_string(pid) + " is not found!");
  }

 private:
  std::map<int, ParticleFit> species_{};
  double minx_{-100000.};
  double maxx_{100000.};

  ClassDefOverride(GetterTof, 1);
};

/**
 * @brief Pid getter based on graphical cut TCutG
 */
class CutGGetter : public TObject, public BaseGetter {

 public:
  void AddParticle(TCutG* cut, int pdgId) {
    if (cut) {
      auto insert_result = species_.insert({pdgId, {cut}});
      if (!insert_result.second) {
        (*insert_result.first).second.push_back(cut);
      }
      return;
    }

    throw std::logic_error("empty TCutG object");
  }

  int GetPid(double var1, double var2, double) override {

    for (const auto& specie : species_) {
      int pdgId = specie.first;
      auto specieCuts = specie.second;

      for (auto cut : specieCuts) {
        if (cut->IsInside(var1, var2)) {
          return pdgId;
        }
      }
    }

    return -1;
  }

  double GetWeight(double var1, double var2, int pid) override {
    return 1.0 * (GetPid(var1, var2, 1) == pid);
  }

  std::map<int, double> GetWeights(double var1, double var2) override {
    std::map<int, double> result;

    for (const auto& specie : species_) {
      result.insert({specie.first, GetWeight(var1, var2, specie.first)});
    }

    return {};
  }

  void Draw(Option_t* option = "") override {
    TObject::Draw(option);

    TMultiGraph mg("mg", "");
    TText pdgLabel;

    for (const auto& specie : species_) {
      auto specieCuts = specie.second;
      for (auto cut : specieCuts) mg.Add(cut);
    }

    mg.DrawClone(option);

    for (const auto& specie : species_) {
      int pdgId = specie.first;
      auto specieCuts = specie.second;

      double xc, yc;
      for (auto cut : specieCuts) {
        cut->Center(xc, yc);
        pdgLabel.DrawText(xc, yc, Form("%d", pdgId));
      }
    }
  }

 protected:
  std::map<int, std::vector<TCutG*>> species_{};

  ClassDefOverride(Pid::CutGGetter, 1)
};

}// namespace Pid
#endif// PidGetterTof_H
