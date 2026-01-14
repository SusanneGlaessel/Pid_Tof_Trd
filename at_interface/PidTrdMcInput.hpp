#ifndef PIDTRD_MCINPUT_HPP_
#define PIDTRD_MCINPUT_HPP_

#include "TH2F.h"

#include "AnalysisTree/Task.hpp"
#include "AnalysisTree/TaskManager.hpp"

#include "ContainerTrd.h"

using std::array;

class PidTrdMcInput : public AnalysisTree::Task {

 public:
  PidTrdMcInput(const std::string& outfilename) {
    outfilename_ = outfilename;
  };
  ~PidTrdMcInput() override = default;

  void Init() override;
  void InitMcHistograms();
  void OpenMcHistograms();
  void Exec() override;
  void Finish() override;

  void SetRecTracksName(const std::string& name) { rec_tracks_name_ = name; };
  void SetSimTracksName(const std::string& name) { sim_tracks_name_ = name; };
  void SetTrdTracksName(const std::string& name) { trd_tracks_name_ = name; };
  void SetUpdateMcHistos(const Bool_t update_mchisto) { update_mchisto_ = update_mchisto; };

 protected:
  float GetMomentum(const AnalysisTree::BranchChannel& trd_track);
  float GetPt(const AnalysisTree::BranchChannel& trd_particle);
  int GetMcPdg(const AnalysisTree::BranchChannel& rec_track);
  int GetCharge(const AnalysisTree::BranchChannel& rec_track);
  void GetEnergyLossHits(const AnalysisTree::BranchChannel& trd_track, int& nhits_trd, array<float, NumberOfTrdLayers>& dEdx);
  void CalculateEnergyLoss(TrdContainer track, float& dEdx_track, int trunc_mode);
  void CreateMcHistograms(int trunc_mode);
  void FillHistogram(TrdContainer track, float dEdx);
  void FillHistogram(TrdContainer track);
  bool IsRichElectron(const AnalysisTree::BranchChannel& rec_particle);

  AnalysisTree::Branch rec_tracks_;
  AnalysisTree::Branch sim_tracks_;
  AnalysisTree::Branch trd_tracks_;

  AnalysisTree::Matching* rec_to_sim_{nullptr};
  AnalysisTree::Matching* rec_to_trd_{nullptr};

  std::string rec_tracks_name_{"RecTracks"};   // Branch with input tracks
  std::string trd_tracks_name_{"TrdTracks"};   // Branch with Trd info (dEdx)
  std::string sim_tracks_name_{"SimParticles"};// Branch with Mc info

  std::vector<AnalysisTree::Field> trd_dEdx_field_{};

  TString outfilename_;
  TFile* outFile_;

  // MC histograms

  TString histnames_hits_all_pos_, histnames_hits_all_neg_;
  array<TString, NumberOfPidsTrd - 1> histnames_hits_pos_, histnames_hits_neg_;
  array<TString, NumberOfTrdLayers> histnames_all_pos_, histnames_all_neg_;
  array<array<TString, NumberOfPidsTrd - 1>, NumberOfTrdLayers> histnames_pos_, histnames_neg_;

  TString outfile_info_{""};

  TH2F* h2dEdx_p_pos_[NumberOfTruncMode - 1];
  TH2F* h2dEdx_p_neg_[NumberOfTruncMode - 1];
  array<TH2F*, (NumberOfPidsTrd - 1) * 2> h2dEdx_p_pdg_[NumberOfTrdLayers];
  TH2F* h2dEdx_hits_all_p_pos_;
  TH2F* h2dEdx_hits_all_p_neg_;
  array<TH2F*, (NumberOfPidsTrd - 1) * 2> h2dEdx_hits_all_p_pdg_;
  TH2F* h2dEdx_hits_p_pos_[NumberOfTruncMode - 1];
  TH2F* h2dEdx_hits_p_neg_[NumberOfTruncMode - 1];
  array<TH2F*, (NumberOfPidsTrd - 1) * 2> h2dEdx_hits_p_pdg_[NumberOfTrdLayers];

  Int_t nbins_mom_;
  Double_t bins_mom_[NbinsMax];
  Int_t nbins_dEdx_;
  Double_t bins_dEdx_[NbinsMax];

  Bool_t update_mchisto_{kFALSE};

  std::map<int, std::vector<TrdContainer>> trd_container_;//container with trd tracks
};

#endif//PIDTRD_AT_INTERFACE_PIDFILLER_HPP_
