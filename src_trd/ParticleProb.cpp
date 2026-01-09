#include "ParticleProb.h"

ClassImp(PidTrd::ParticleProb)

namespace PidTrd {

/**
* Gets probability depending on momentum and dEdx
* @param momentum
* @param dEdx energy loss in the Trd
* @return probability for a particle specie
*/
    
  float ParticleProb::Eval(float mom, float dEdx)
  {
    Int_t binx = hprobabilities_->GetXaxis()->FindBin(mom);
    Int_t biny = hprobabilities_->GetYaxis()->FindBin(dEdx);
    return hprobabilities_->GetBinContent(binx,biny);
  }
}
    
