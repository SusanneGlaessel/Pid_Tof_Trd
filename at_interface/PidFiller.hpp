#ifndef PID_INTERFACE_PIDFILLER_HPP_
#define PID_INTERFACE_PIDFILLER_HPP_

#include "ConstantsTof.h"
#include "ConstantsTrd.h"
#include "GetterTof.h"
#include "GetterTrd.h"
#include "ContainerTrd.h"

#include "AnalysisTree/Task.hpp"
#include "AnalysisTree/TaskManager.hpp"

#include "TH2F.h"

using std::make_pair;

class PidFiller : public AnalysisTree::Task {

 public:
  PidFiller(const std::string& pid_file_name_tof, const std::string& pid_file_name_trd, const std::string& getter_name_tof,  const std::string& getter_name_trd, int pid_mode);
  ~PidFiller() override { delete getter_tof_; delete getter_trd_; }

  void Init() override;
  void Exec() override;
  void Finish() override {
    auto* man = AnalysisTree::TaskManager::GetInstance();
    //    man->RemoveBranch(rec_tracks_name_);
  }

  void SetRecTracksName(const std::string& name) { rec_tracks_name_ = name; }
  void SetTofHitsName(const std::string& name)   { tof_hits_name_   = name; }
  void SetTrdTracksName(const std::string& name) { trd_tracks_name_ = name; };
  void SetRichRingsName(const std::string& name) { rich_rings_name_ = name; };

  // Settings for Trd Pid
  void SetMinHits(int nhits_min) { nhits_min_ = nhits_min; }            // Min. number of hits per track
  void SetTruncationMode(int trunc_mode) { trunc_mode_ = trunc_mode; }  // Calculation of energy loss for up to 4 layers:
                                                                        // =0: <dEdx> average over all hits
                                                                        // =1-4: Select hits with lowest dEdx:
                                                                        // =1: 1 hit, =2: 2 hits, =3: 3 hits, =4: 4 hits
  void SetProbabilityMode(int prob_mode) { prob_mode_ = prob_mode; }    // Probability for particle species i:
                                                                        // =0: total probability - probability based on particle multiplicites i
									// =1: likelihood - probability based on dEdx-distribution of particle
  
  void SetPurity(const float purity) { purity_ = purity; }
  
 protected:
  int signum(int x) const;


  float GetMomentumTrd(const AnalysisTree::BranchChannel& trd_particle);
  float GetPtTrd(const AnalysisTree::BranchChannel& trd_particle);
  int GetCharge(const AnalysisTree::BranchChannel& trd_particle);
  void GetEnergyLossHitsTrd(const AnalysisTree::BranchChannel& trd_track, int &nhits_trd, std::array<float, NumberOfTrdLayers> &dEdx);
  bool IsRichElectron(const AnalysisTree::BranchChannel& rec_particle);

  AnalysisTree::Branch rec_tracks_;
  AnalysisTree::Branch tof_hits_;
  AnalysisTree::Branch trd_tracks_;
  AnalysisTree::Branch rich_rings_;
  AnalysisTree::Branch ana_tracks_;
  AnalysisTree::Matching* rec_to_tof_{nullptr};
  AnalysisTree::Matching* rec_to_trd_{nullptr};
  AnalysisTree::Matching* rec_to_rich_{nullptr};

  // Tof fields
  AnalysisTree::Field qp_tof_field_;
  AnalysisTree::Field q_field_;
  AnalysisTree::Field mass2_field_;
  AnalysisTree::Field pid_tof_field_;
  std::vector<AnalysisTree::Field> prob_tof_field_{};

  // Trd fields
  std::vector<AnalysisTree::Field> trd_dEdx_field_{};  
  AnalysisTree::Field el_rich_field_;
  AnalysisTree::Field pid_trd_field_;
  AnalysisTree::Field nhits_trd_eloss_field_;
  std::vector<AnalysisTree::Field> prob_trd_field_{};

  std::vector<AnalysisTree::Matching*> in_matches_{};
  std::vector<AnalysisTree::Matching*> out_matches_{};

  std::string rec_tracks_name_{"VtxTracks"};    // Branch with input tracks
  std::string tof_hits_name_{"TofHits"};        // Branch with TOF info (m2)
  std::string trd_tracks_name_{"TrdTracks"};    // Branch with Trd info (dEdx)
  std::string rich_rings_name_{"RichRings"};    // Branch with Rich info
  std::string out_branch_name_{"RecParticles"}; // Output branch (based on VtxTracks) with pid info: probabilities and particle type hypothesis

  bool is_run_pid_tof_{false};
  bool is_run_pid_trd_{false};

  int trunc_mode_{0};
  int prob_mode_{0};
  float purity_{0.0};
  int nhits_min_{1};

  Pid::GetterTof* getter_tof_{nullptr};
  std::vector<std::pair<long long, std::string>> pid_codes_tof_{
      {PidTofParticles::kProton, "p"},
      {PidTofParticles::kPionPos, "pi"},
      {PidTofParticles::kKaonPos, "K"},
      {PidTofParticles::kDeutron, "d"},
      {PidTofParticles::kBgPos, "bg"}};

  PidTrd::GetterTrd* getter_trd_{nullptr};
};

#endif// PID_INTERFACE_PIDFILLER_HPP_
