#include "ContainerTrd.h"

/**
* Scales energy loss to the length of its trace through the Trd
*/

void TrdContainer::ScaleEnergyLossLength() {

  if (dEdx_is_scaled_ == true) {
     Warning("Exec", "Energy loss is already scaled.");
     return;
  }
  
  float pz = 0.0;
  if (TMath::Abs(mom_) > TMath::Abs(pT_))
    pz = TMath::Sqrt(mom_*mom_ - pT_*pT_);
  else 
    Warning("Exec", "Could not assign pz to the track, use pz=0. Energy loss will not be scaled.");

  if ( pz > 0 && mom_ > 0) {
    for (int ihit = 0 ; ihit < NumberOfTrdLayers; ihit++) 
      dEdx_hits_.at(ihit) *= pz / mom_;

    dEdx_is_scaled_ = true;
  }
}

/**
* Sorts hits of a track depending on its dEdx-value from lowest to highest
*/

std::array<float,NumberOfTrdLayers> TrdContainer::GetdEdxHitsSorted () {
      
  std::vector<float> dEdx;
  dEdx.clear();

  for (int ihit = 0; ihit < NumberOfTrdLayers; ihit++)
    if (dEdx_hits_.at(ihit) > 0)
      dEdx.push_back(dEdx_hits_.at(ihit));
  
  std::sort(dEdx.begin(), dEdx.end());

  for (int ihit = 0; ihit < dEdx.size(); ihit++)
    hits_sorted_.at(ihit) = dEdx.at(ihit);

  for (int ihit = dEdx.size(); ihit < NumberOfTrdLayers; ihit++)
    hits_sorted_.at(ihit) = 0.0;
  
  return hits_sorted_;
}

/**
* Gets number of hits for a track according to the truncation mode
* @param truncation mode
* @return number of selected hits
*/

int TrdContainer::GetNHitsSel(int trunc_mode) {

  if (trunc_mode == 0)
    nhits_sel_ = nhits_trd_;
  else
    nhits_sel_ = trunc_mode;
  if (nhits_sel_ > nhits_trd_)
    nhits_sel_ = nhits_trd_;
  
  return nhits_sel_;
}

/**
* Selects hits of track according to the truncation mode
* @param truncation mode
*/

void TrdContainer::SelectHitIndices(int trunc_mode) {

  if (nhits_sel_ == 0) nhits_sel_ = GetNHitsSel(trunc_mode);

  std::array<float, NumberOfTrdLayers> dEdx_hits_tmp = dEdx_hits_;
  
  for (int ihit = 0; ihit < NumberOfTrdLayers; ihit++)
    if (dEdx_hits_tmp.at(ihit) <= 0.0) dEdx_hits_tmp.at(ihit) = std::numeric_limits<float>::max();
  
  for (int ihit = 0; ihit < nhits_sel_; ihit++) {
    auto dEdx_min = std::min_element(std::begin(dEdx_hits_tmp), std::end(dEdx_hits_tmp));
    auto dEdx_min_index = std::distance(dEdx_hits_tmp.begin(), dEdx_min);
    hits_sel_index_.at(dEdx_min_index) = true;
    dEdx_hits_tmp.at(dEdx_min_index) = std::numeric_limits<float>::max();
  }
  
  indices_are_sel_  = true;
}

/**
* Calculates energy loss for a track for a selected truncation mode
* @param truncation mode
*/

void TrdContainer::CalculateEnergyLossTrack(int trunc_mode) {

  if (nhits_sel_ == 0) nhits_sel_ = GetNHitsSel(trunc_mode);
  
  std::vector<float> dEdx;
  dEdx.clear();

  for (int ihit = 0; ihit < NumberOfTrdLayers; ihit++)
    if (dEdx_hits_.at(ihit) > 0)
      dEdx.push_back(dEdx_hits_.at(ihit));
  
  std::sort(dEdx.begin(), dEdx.end());

  float dEdx_sum = 0;
  for (int ihit = 0 ; ihit < nhits_sel_; ihit++) {
    dEdx_sum += dEdx.at(ihit);
  }
  dEdx_track_.at(trunc_mode) = dEdx_sum / nhits_sel_;
}

/**
* Calculates energy loss for a track for a all truncation modes
*/

void TrdContainer::CalculateEnergyLossTrackAllModes() {

  std::vector<float> dEdx;
  dEdx.clear();

  for (int ihit = 0; ihit < NumberOfTrdLayers; ihit++)
    if (dEdx_hits_.at(ihit) > 0)
      dEdx.push_back(dEdx_hits_.at(ihit));

  std::sort(dEdx.begin(), dEdx.end());

  float dEdx_sum = 0;
  for (int ihit = 0 ; ihit < dEdx.size(); ihit++) {
    dEdx_sum += dEdx.at(ihit);
    dEdx_track_.at(ihit) = dEdx_sum / (ihit + 1);
  }
}
