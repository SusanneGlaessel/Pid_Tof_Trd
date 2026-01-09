/** @file   ParticleProb.h
    @author Susanne Glaessel (glaessel@ikf.uni-frankfurt.de)
    @brief  Class to store probabilities for particle species, number of hits, truncation mode and probability mode
*/

#ifndef ParticleProb_H
#define ParticleProb_H

#include "ConstantsTrd.h"

#include "TString.h"
#include "TFile.h"
#include "TH2F.h"

using std::make_pair;
using std::array;

namespace PidTrd {

  class ParticleProb : public TObject {

  public:
    ParticleProb() = default;
    virtual ~ParticleProb() = default;

    ParticleProb(int type, int charge, int nhits, int truncmode, int probmode, TH2F* hprobabilities) { particle_type_ = type; charge_ = charge; nhits_ = nhits, truncmode_ = truncmode, probmode_ = probmode, hprobabilities_ = hprobabilities; };

    void Update(int type, int charge, int nhits, int truncmode, int probmode, TH2F* hprobabilities) { particle_type_ = type; charge_ = charge; nhits_ = nhits, truncmode_ = truncmode, probmode_ = probmode, hprobabilities_ = hprobabilities; };

    float Eval(float mom, float dEdx);

    int GetId()              {int type_pm = charge_ > 0 ? particle_type_ : particle_type_ + NumberOfPidsTrd-1;
                              return type_pm + nhits_ * 100 + truncmode_ * 1000 + probmode_ * 10000; };
    int GetType()            {return particle_type_; };
    int GetCharge()          {return charge_; };
    int GetNhits()           {return nhits_; };
    int GetTruncMode()       {return truncmode_; };
    int GetProbMode()        {return probmode_; };
    TH2F* GetProbabilities() {return hprobabilities_; };

  private:
    int particle_type_{-1};
    int charge_{0};
    int nhits_{0};
    int truncmode_{0};        
    int probmode_{0};
    TH2F* hprobabilities_;
    
    ClassDef(ParticleProb, 1);
  };
};
#endif
    
