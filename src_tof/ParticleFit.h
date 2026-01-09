/** @file   ParticleFit.h
    @class  Pid::ParticleFit
    @author Viktor Klochkov (klochkov44@gmail.com)
    @author Ilya Selyuzhenkov (ilya.selyuzhenkov@gmail.com)
    @brief  Class to store fit resuls for particle specie
*/

#ifndef PidParticleFit_H
#define PidParticleFit_H 1

#include <array>
#include <vector>

#include "ConstantsTof.h"
#include "TF1.h"

namespace Pid {

class ParticleFit {
 public:
  /**   Default constructor   **/
  ParticleFit();
  explicit ParticleFit(int type) : particle_type_(type){};
  virtual ~ParticleFit() = default;

  [[nodiscard]] std::vector<double> GetFunctionParams(float p) const;
  double Eval(double p, double m2);

  void SetParametrization(const std::vector<TF1>& parametrization) { parametrization_ = parametrization; }
  void SetParametrizationFunction(uint ivar, const TF1& func) { parametrization_.at(ivar) = func; }
  void SetFitFunction(const TF1& function) { function_ = function; }
  void SetRange(float min, float max) { minx_ = min, maxx_ = max; }
  void SetIsFitted(bool is = true) { isfitted_ = is; }
  void SetIsFixed(const std::vector<bool>& is) { isfixed_ = is; }

  [[nodiscard]] const TF1& GetFunction() const { return function_; }
  [[nodiscard]] uint GetNpar() const { return function_.GetNpar(); }
  TF1& GetParametrizationFunction(int ipar) { return parametrization_.at(ipar); }
  [[nodiscard]] bool GetIsFixed(uint ipar) const {
    if (ipar >= isfixed_.size())
      return false;
    return isfixed_.at(ipar);
  }

  [[nodiscard]] double GetSigma(float p) const { return parametrization_.at(PidFunction::kSigma).Eval(p); }
  [[nodiscard]] double GetMean(float p) const { return parametrization_.at(PidFunction::kMean).Eval(p); }
  [[nodiscard]] double GetIntegral(float p) const { return parametrization_.at(PidFunction::kA).Eval(p) / sqrt(2 * TMath::Pi() / GetSigma(p)); }

  void GetRange(float& min, float& max) const { min = minx_, max = maxx_; }

 private:
  TF1 function_;
  std::vector<TF1> parametrization_{};
  std::vector<bool> isfixed_{};

  float minx_{-1.};
  float maxx_{-1.};

  int npoints_{-1};
  int particle_type_{-1};

  bool isfitted_{false};

  ClassDef(ParticleFit, 2);
};

}// namespace Pid
#endif// PidParticleFit_H
