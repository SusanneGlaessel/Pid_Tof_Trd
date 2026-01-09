#include "GetterTrd.h"
#include <iostream>
using std::cout;
using std::endl;
using std::to_string;

ClassImp(PidTrd::GetterTrd)

namespace PidTrd {

  /**
* Gets probabilites for a hit or track for a selected truncation and probability mode. 
* @param trdtrack 
* @param ihit - for likelihood method (= -1 for total probability method)
* @return map with probabilities for a track or hit for all particle species
*/
  std::map<int, float> GetterTrd::GetTrdProbabilities(TrdContainer trdtrack, int ihit) {
    std::map<int, float> prob;

    float prob_bg = 1.0;
    float sum = 0.0;
    
    int nhits = trdtrack.GetNhitsTrd();
    float mom = trdtrack.GetP();
    float dEdx;

    if (prob_mode_ == 0) {
      if (ihit != -1)
	throw std::runtime_error("For probability mode 0 (total probability) function GetTrdProbabilities need to be used and ihit should not be set.");
      dEdx = trdtrack.GetdEdxTrack(trunc_mode_);
    }
    if (prob_mode_ == 1) {
      if (ihit == -1)
	throw std::runtime_error("For probability mode 1 (likelihood) function GetTrdProbabilitiesMulti need to be used.");
      dEdx = trdtrack.GetdEdxHits().at(ihit);
    }
    
    for(size_t ipid = 0; ipid < NumberOfPidsTrd - 1; ipid++) {
      Int_t ipid_pm = trdtrack.GetCharge() > 0 ? ipid : ipid + NumberOfPidsTrd-1;
      Int_t trunc_mode_getter;
      if (trunc_mode_ == 0) trunc_mode_getter = trdtrack.GetNhitsTrd() - 1;
      else if (trunc_mode_ > trdtrack.GetNhitsTrd()) trunc_mode_getter = trdtrack.GetNhitsTrd() - 1; // if (truncation > nhits) truncation = truncation of nhits
      else trunc_mode_getter = trunc_mode_ - 1;
      
      Int_t id = ipid_pm + (trdtrack.GetNhitsTrd() - 1) * 100 + trunc_mode_getter * 1000 + prob_mode_ * 10000;
      prob[pid_codes_trd_.at(ipid).first] = particles_prob_.find(id)->second.Eval(mom,dEdx);
      prob_bg -= prob[pid_codes_trd_.at(ipid).first];
    }
  
    if (prob_bg < 0) prob_bg = 0.0;
    prob[PidTrdParticles::kBgPos] = prob_bg;

    return prob;
  }

  /**
* Multiplies probabilites of selected hits of a track for a selected truncation and probability mode. 
* @param trdtrack 
* @return map with probabilities for all particle species for a track
*/
  
  std::map<int, float> GetterTrd::GetTrdProbabilitiesMulti(TrdContainer trdtrack) {
    std::map<int, float> prob;       
    for (const auto& pdg : pid_codes_trd_) 
      prob [pdg.first] = 1.0;
    
    std::map<int, float> prob_tmp;
    
    trdtrack.SelectHitIndices(trunc_mode_);

    for (int ihit = 0; ihit < NumberOfTrdLayers; ihit++) {
      if (trdtrack.GetHitsSelIndex().at(ihit) == false) continue;   
      prob_tmp = GetTrdProbabilities(trdtrack, ihit);

      for(size_t ipid = 0; ipid < NumberOfPidsTrd - 1; ipid++) {
	prob [pid_codes_trd_.at(ipid).first] *= prob_tmp [pid_codes_trd_.at(ipid).first];
      }
    }

    float prob_tot = 0.0;
    for(size_t ipid = 0; ipid < NumberOfPidsTrd - 1; ipid++)
      if (prob [pid_codes_trd_.at(ipid).first] >= 0. && prob[pid_codes_trd_.at(ipid).first] <= 1.) prob_tot += prob[pid_codes_trd_.at(ipid).first];
    
    for(size_t ipid = 0; ipid < NumberOfPidsTrd - 1; ipid++) {
      if (prob_tot > 0) prob[pid_codes_trd_.at(ipid).first] /= prob_tot;
      else prob[pid_codes_trd_.at(ipid).first]= 0.0;
    }

    if (prob_tot > 0)
      prob [PidTrdParticles::kBgPos] = 0.0;
    else
      prob [PidTrdParticles::kBgPos] = 1.0;
     
    return prob;
  }
  
   /**
* Gets probabilites for a track for the average dEdx with the total probability method. 
* @param momentum 
* @param average dEdx
* @param number of hits
* @return map with probabilities for a track for all particle species
*/
  
  std::map<int, float> GetterTrd::GetTrdProbabilitiesTotalAverage(float mom, float dEdx, int charge, int nhits, int prob_mode) {

    std::map<int, float> prob;
    float prob_bg = 1.0;
    
    for(size_t ipid = 0; ipid < NumberOfPidsTrd - 1; ipid++) {
      Int_t ipid_pm = charge > 0 ? ipid : ipid + NumberOfPidsTrd-1;
      Int_t trunc_mode_getter = nhits - 1;
      
      Int_t id = ipid_pm + (nhits - 1) * 100 + trunc_mode_getter * 1000 + prob_mode * 10000;
      prob[pid_codes_trd_.at(ipid).first] = particles_prob_.find(id)->second.Eval(mom,dEdx);
      prob_bg -= prob[pid_codes_trd_.at(ipid).first];
    }
    if (prob_bg < 0) prob_bg = 0.0;
    prob[PidTrdParticles::kBgPos] = prob_bg;
    return prob;
  }

    /**
* Gets probabilites for a track for the average dEdx with the likelihood method. 
* @param momentum 
* @param vector of dEdx
* @param number of hits
* @return map with probabilities for a track for all particle species

* Multiplies probabilites of selected hits of a track for a selected truncation and probability mode. 
* @param trdtrack 
* @return map with probabilities for all particle species for a track
*/
  
  std::map<int, float> GetterTrd::GetTrdProbabilitiesLikelihoodAverage(float mom, std::array<float, NumberOfTrdLayers> dEdx, int charge) {
 
    std::map<int, float> prob;       
    for (const auto& pdg : pid_codes_trd_) 
      prob [pdg.first] = 1.0;
    
    std::map<int, float> prob_tmp;
    int nhits = 0;

    for (int ihit = 0; ihit < NumberOfTrdLayers; ihit++) {
      if (dEdx.at(ihit) <= 0.0) continue;
      nhits ++;
    }
 
    for (int ihit = 0; ihit < NumberOfTrdLayers; ihit++) {
      if (dEdx.at(ihit) <= 0.0) continue;
      prob_tmp = GetTrdProbabilitiesTotalAverage(mom, dEdx.at(ihit), charge, nhits, 1);
      for(size_t ipid = 0; ipid < NumberOfPidsTrd - 1; ipid++) {
	prob [pid_codes_trd_.at(ipid).first] *= prob_tmp [pid_codes_trd_.at(ipid).first];
      }
    }

    float prob_tot = 0.0;
    for(size_t ipid = 0; ipid < NumberOfPidsTrd - 1; ipid++) 
      if (prob [pid_codes_trd_.at(ipid).first] >= 0. && prob [pid_codes_trd_.at(ipid).first] <= 1.) prob_tot += prob [pid_codes_trd_.at(ipid).first];
    for(size_t ipid = 0; ipid < NumberOfPidsTrd - 1; ipid++)  {
      if (prob_tot > 0) prob [pid_codes_trd_.at(ipid).first] /= prob_tot;
      else prob [pid_codes_trd_.at(ipid).first] = 0.0;
    }
    if (prob_tot > 0)
      prob [PidTrdParticles::kBgPos] = 0.0;
    else
      prob [PidTrdParticles::kBgPos] = 1.0;
    return prob;
  }
  
  /**
* Gets pid hypothesis by selecting the particle specie with the highest probability
* @param map with probabilities for all particle species
* @param minium selected purity
* @param charge
* @return pid
*/

  int GetterTrd::GetTrdPid(std::map<int, float> prob, float purity, int charge) {
    int pid;

    std::array<float, NumberOfPidsTrd> prob_vec;
    for (int i = 0; i < NumberOfPidsTrd; i++)
      prob_vec [i] = prob [pid_codes_trd_.at(i).first];
			       
    auto prob_max = std::max_element(std::begin(prob_vec), std::end(prob_vec));
    auto prob_max_index = std::distance(prob_vec.begin(), prob_max);
    if (*prob_max >= purity)
      pid = pid_codes_trd_.at(prob_max_index).first*charge;
    else
      pid = 1*charge;
    return pid;
  }
}
