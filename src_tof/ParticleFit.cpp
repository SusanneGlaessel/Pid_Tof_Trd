#include "ParticleFit.h"

ClassImp(Pid::ParticleFit)

    namespace Pid {

  /// Default constructor   -------------------------------------------
  ParticleFit::ParticleFit() = default;

  /**
* Getter for vector of parameters for a ParticleFit
* @param p track momentum
* @return vector of parameters
*/
  std::vector<double> ParticleFit::GetFunctionParams(float p) const {
    std::vector<double> params;
    const uint npar = function_.GetNpar();

    if (parametrization_.size() != npar) {
      std::cout << "\n\nNumber of parameters is not equal to number of parameter functions!\nExiting...\n";
      exit(0);
    }

    if (!isfitted_)
      return params;

    for (uint i = 0; i < npar; ++i) {
      params.push_back(parametrization_.at(i).Eval(p));
    }

    return params;
  }

  /**
* Similar to TF1::Eval
* @param p track momentum
* @param m2 track mass square
* @return vector of parameters
*/
  double ParticleFit::Eval(double p, double m2) {
    if (p > maxx_ || p < minx_) return 0.;

    const uint npar = function_.GetNpar();
    if (parametrization_.size() != npar)
      exit(1);

    //  TF1* f2 = (TF1*) function_.Clone();
    //  f2->SetParameters(&(GetFunctionParams(p)[0]));
    //  double ret = f2->Eval(m2);
    //  delete f2;

    function_.SetParameters(&(GetFunctionParams(p)[0]));
    const double ret = function_.Eval(m2);

    return ret;
  }
}
