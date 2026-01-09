#include "GetterTof.h"

ClassImp(Pid::GetterTof)
    ClassImp(Pid::CutGGetter)

        namespace Pid {

  /**
* Constructs fit function as a sum of individual particle species. Parameters are also propagated
* @param p track momentum (or x-axis value, in general)
* @param m2 square mass for TOF (or y-axis value, in general)
* @return map with probabilities for all particle species
*/
  std::map<int, double> GetterTof::GetBayesianProbability(double p, double m2) {
    std::map<int, double> prob{};

    //  if (p > maxx_ || p < minx_)  //NOTE should be done in cuts maybe?
    //    return prob;

    double sum{0.};
    for (auto& specie : species_) {
      const double iprob = specie.second.Eval(p, m2);
      sum += iprob;
      prob[specie.first] = iprob;
    }
    //     std::cout << sum << std::endl;
    if (sum < 1.) return prob;// no particles in (p, m2) point

    for (auto& iprob : prob) {
      iprob.second /= sum;
    }

    return prob;
  }
}
