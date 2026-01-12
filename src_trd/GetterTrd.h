/** @file   GetterTrd.h
    @author Susanne Glaessel (glaessel@ikf.uni-frankfurt.de)
    @brief  Class to calculate Trd-PID probabilities
*/

#ifndef PidGetterTrd_H
#define PidGetterTrd_H 1

#include "TObject.h"
#include "TString.h"
#include "ParticleProb.h"
#include "ContainerTrd.h"

namespace PidTrd {
  
  /**
   * @brief Pid Trd getter
   */  
  class GetterTrd : public TObject {
    
  public:
    GetterTrd() = default;
    virtual ~ GetterTrd() = default;

    void SetMinHits(int nhits_min) { nhits_min_ = nhits_min; }            // Min. number of required hits per track
    void SetTruncationMode(int trunc_mode) { trunc_mode_ = trunc_mode; }  // Modes for calculating the energy loss for up to 4 layers:
                                                                          // =0: <dEdx> average over all hits
                                                                          // =1-4: Select hits with lowest dEdx:
                                                                          //       =1: 1 hit, =2: 2 hits, =3: 3 hits, =4: 4 hits
    void SetProbabiltyMode(int prob_mode) { prob_mode_ = prob_mode; }     // =0: total probability - probability based on particle multiplicites
                                                                          // =1: likelihood - probability based on dEdx-distribution of particle species
    
    int GetMinHits() { return nhits_min_; }                               

    void AddParticlesProb(std::map<int, ParticleProb> particlesprob) {particles_prob_ = particlesprob;}
    void AddParticleProb(ParticleProb particleprob) {
      int type = particleprob.GetCharge() > 0 ? particleprob.GetType() : particleprob.GetType() + NumberOfPidsTrd-1;
      int nhits = particleprob.GetNhits();
      int truncmode = particleprob.GetTruncMode();
      int probmode = particleprob.GetProbMode();  
      particles_prob_[type + nhits * 100 + truncmode * 1000 + probmode * 10000] = particleprob;
    }
    
    std::map<int, ParticleProb> GetParticlesProb() { return particles_prob_; }
    ParticleProb GetParticleProb(int id) { return particles_prob_[id];}
    ParticleProb GetParticleProb(int type, int charge, int nhits, int truncmode, int probmode) {
      if (charge < 0) type += NumberOfPidsTrd-1;
      return particles_prob_[type + nhits * 100 + truncmode * 1000 + probmode * 10000];
    }
    
    std::map<int, float> GetTrdProbabilities(TrdContainer trdtrack, int ihit = -1);
    std::map<int, float> GetTrdProbabilitiesMulti(TrdContainer trdtrack);
    std::map<int, float> GetTrdProbabilitiesTotalAverage(float mom, float dEdx, int charge, int nhits, int probmode = 0);
    std::map<int, float> GetTrdProbabilitiesLikelihoodAverage(float mom, std::array<float,NumberOfTrdLayers> dEdx, int charge);
    int GetTrdPid(std::map<int, float> prob, float purity, int charge);
    
  private:

    std::map<int, ParticleProb> particles_prob_{};
    
    int nhits_min_{1};
    int trunc_mode_{0};
    int prob_mode_{0};

    array<TString, NumberOfProbMode> probnames_ = {"probT", "probL"};
    
    ClassDef(GetterTrd, 1);
  }; 
}// namespace PidTrd
#endif//PidTrd_Getter_H
