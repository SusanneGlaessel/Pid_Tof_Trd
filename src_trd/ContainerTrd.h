/** @file   ContainerTrd.h
    @author Susanne Glaessel (glaessel@ikf.uni-frankfurt.de)
    @brief  Class to store Trd PID information
*/

#ifndef ContainerTrd_H
#define ContainerTrd_H

#include "ConstantsTrd.h"

#include "TMath.h"
#include <stdexcept>

class TrdContainer {
 public:
  TrdContainer() = default;
  
  TrdContainer(float mom, float pT, int charge, int nhits_trd, std::array<float,NumberOfTrdLayers> dEdx_hits, std::array<float, NumberOfTruncMode> dEdx_track, int nhits_sel, std::array<bool, NumberOfTrdLayers> hits_sel_index, int mc_pdg, bool dEdxIsScaled = false) : mom_(mom), pT_(pT), charge_(charge), nhits_trd_(nhits_trd), dEdx_hits_(dEdx_hits), dEdx_track_(dEdx_track), nhits_sel_(nhits_sel), hits_sel_index_ (hits_sel_index), mc_pdg_(mc_pdg), dEdx_is_scaled_(dEdxIsScaled) {}
  
  TrdContainer(const float mom, float pT, int charge, int nhits_trd, std::array<float,NumberOfTrdLayers> dEdx_hits, int mc_pdg, bool dEdxIsScaled = false) : mom_(mom), pT_(pT), charge_(charge), nhits_trd_(nhits_trd), dEdx_hits_(dEdx_hits), mc_pdg_(mc_pdg), dEdx_is_scaled_(dEdxIsScaled) {
    nhits_sel_ = 0;
    for (int ihit = 0; ihit < NumberOfTrdLayers; ihit++)
      hits_sel_index_.at(ihit) = false;
  }
  
  TrdContainer(float mom, float pT, int charge, int nhits_trd, std::array<float,NumberOfTrdLayers> dEdx_hits, bool dEdxIsScaled = false) : mom_(mom), pT_(pT), charge_(charge), nhits_trd_(nhits_trd), dEdx_hits_(dEdx_hits), dEdx_is_scaled_(dEdxIsScaled) {
    nhits_sel_ = 0;
    for (int ihit = 0; ihit < NumberOfTrdLayers; ihit++)
      hits_sel_index_.at(ihit) = false;
    mc_pdg_ = -2;
  }

  virtual ~TrdContainer() = default;

  void ScaleEnergyLossLength();
  void CalculateEnergyLossTrack(int trunc_mode); // Calculates energy loss of a track for a selected truncation mode
  void CalculateEnergyLossTrackAllModes();       // Calculates energy loss of a track for truncation modes 1-4
  void SelectHitIndices(int trunc_mode);         // Returns hit indices with dEdx > 0 selected in truncation mode
  int  GetNHitsSel(int trunc_mode);              // Returns number of hits with dEdx > 0 selected in truncation mode
  
  float GetP()       const { return mom_;       }
  int GetCharge()    const { return charge_;    }
  int GetNhitsTrd()  const { return nhits_trd_; }
  int GetMcPdg()     const {
    if (mc_pdg_ == -2) throw std::runtime_error("MC PDG not set.");
    else return mc_pdg_;}
  
  std::array<float,NumberOfTrdLayers> GetdEdxHits() const { return dEdx_hits_; }
  float GetdEdxTrack(int trunc_mode) {
    if (dEdx_track_.at(trunc_mode) == 0)
      CalculateEnergyLossTrack(trunc_mode);
    return dEdx_track_.at(trunc_mode);
  }
  std::array<float, NumberOfTrdLayers> GetdEdxHitsSorted (); // Orders hits from lowest to highest dEdx (hits with dEdx = 0 in last posisition)
  std::array<bool, NumberOfTrdLayers> GetHitsSelIndex ()   const {
    if ( indices_are_sel_  == false)
      throw std::runtime_error("Hits are not selected. Call function SelectHitIndices(trunc_mode) first.\n");
    else
      return hits_sel_index_;
  }
  
protected:
  float mom_{0.0};
  float pT_{0.0};;
  int charge_{0};
  int nhits_trd_{0};
  int mc_pdg_{-1};

  std::array<float,NumberOfTrdLayers> dEdx_hits_ = {0.0, 0.0, 0.0, 0.0};             // energy loss of 1-4 hits of a track
  std::array<float, NumberOfTruncMode> dEdx_track_ = {0.0, 0.0, 0.0, 0.0, 0.0};      // energy loss of a track for truncation modes 1-4

  std::array<float, NumberOfTrdLayers> hits_sorted_ = {0.0, 0.0, 0.0, 0.0};          // dEdx values of hits sorted from lowest to higest dEdx
  int nhits_sel_{0}; 
  std::array<bool, NumberOfTrdLayers> hits_sel_index_= {false, false, false, false}; // hits (in their original) order are selected according to truncation mode

  bool dEdx_is_scaled_{false};
  bool indices_are_sel_{false};
  
};

#endif//ContainerTrd_HPP
