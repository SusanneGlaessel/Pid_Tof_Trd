#ifndef PIDTRD_RUNGETTER_HPP_
#define PIDTRD_RUNGETTER_HPP_

#include "TH2F.h"

#include "ParticleProb.h"
#include "GetterTrd.h"

using std::make_pair;
using std::array;

class PidTrdRunGetter {

 public:
  PidTrdRunGetter(const std::string& mcfile_name, const std::string& getter_file, TString getter_name) {
    mcfile_name_ = mcfile_name;
    getter_file_ = getter_file;
    getter_name_ = getter_name;
  };
  virtual ~PidTrdRunGetter() = default;

  void Init();
  void OpenMcHistograms();
  void InitMcHistogramsOut();
  void Exec();
  void Finish();
  void SetWriteMcHistogramsOut(const Bool_t write_mchistograms_out) { write_mchistograms_out_ = write_mchistograms_out; }; // Probabilities are written and saved in histograms in addition to Getter

 protected:
  
  void CalculateProbabilities(int trunc_mode);
  void CalculateProbabilitiesTot(TH2F* h2all, TH2F* h2tmp, TH2F* &h2prob); // total probability - probability based on particle multiplicites
  void CalculateProbabilitiesLike(TH2F* h2tmp, TH2F* &h2prob);             // likelihood - probability based on dEdx-distribution of particle
 
  TString mcfile_name_{""};
  TFile *inFile_;
  TFile *outFile_;
  TString getter_file_;
  TString getter_name_;

  // MC histograms

  array<TString, NumberOfProbMode> probname_ = {"probT", "probL"};
  
  array<TString, NumberOfTrdLayers> histnames_all_pos_, histnames_all_neg_;
  array<array<TString, NumberOfPidsTrd - 1>, NumberOfTrdLayers> histnames_pos_, histnames_neg_;
  array<array<TString, NumberOfPidsTrd - 1>, NumberOfTrdLayers> histnames_prob_pos_, histnames_prob_neg_;

  TString histnames_hits_all_pos_, histnames_hits_all_neg_;
  array<TString, NumberOfPidsTrd - 1> histnames_hits_pos_, histnames_hits_neg_;
  
  Int_t nbins_mom_;
  Double_t bins_mom_ [NbinsMax];
  Int_t nbins_dEdx_;
  Double_t bins_dEdx_ [NbinsMax];

  Bool_t write_mchistograms_out_{kFALSE};

  PidTrd::GetterTrd getter_trd_{};
 
};

#endif//PIDTRD_RUNGETTER_HPP_
