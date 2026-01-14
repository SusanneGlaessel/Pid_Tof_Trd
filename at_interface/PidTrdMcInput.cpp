#include "PidTrdMcInput.hpp"
#include "ConstantsTrd.h"

using namespace AnalysisTree;

using std::to_string;

float PidTrdMcInput::GetMomentum(const AnalysisTree::BranchChannel& trd_track) {
  float momentum = 0;
  if (TMath::Abs(trd_track.Value(rec_tracks_.GetField("p")) > 0))
    momentum = trd_track.Value(rec_tracks_.GetField("p"));
  else if (TMath::Abs(trd_track.Value(rec_tracks_.GetField("p_out")) > 0))
    momentum = trd_track.Value(rec_tracks_.GetField("p_out"));
  else
    Warning("Exec", "Could not assign any momentum to the track, use p=0.");
  return momentum;
}

float PidTrdMcInput::GetPt(const AnalysisTree::BranchChannel& trd_track) {
  float pT = 0;
  if (TMath::Abs(trd_track.Value(trd_tracks_.GetField("pT")) > 0))
    pT = trd_track.Value(trd_tracks_.GetField("pT"));
  else if (TMath::Abs(trd_track.Value(trd_tracks_.GetField("pT_out")) > 0))
    pT = trd_track.Value(trd_tracks_.GetField("pT_out"));
  else
    Warning("Exec", "Could not assign any pT to the track, use pT=0.");
  return pT;
}

int PidTrdMcInput::GetMcPdg(const AnalysisTree::BranchChannel& rec_track) {
  const int sim_id = rec_to_sim_->GetMatch(rec_track.GetId());
  int pdg;
  if (sim_id < 0) pdg = -1;
  else
    pdg = sim_tracks_[sim_id].Value(sim_tracks_.GetField("pid"));
  return pdg;
}

int PidTrdMcInput::GetCharge(const AnalysisTree::BranchChannel& rec_track) {
  int charge = rec_track.Value(rec_tracks_.GetField("q"));
  return charge;
}

void PidTrdMcInput::GetEnergyLossHits(const AnalysisTree::BranchChannel& trd_track, int& nhits_trd, array<float, NumberOfTrdLayers>& dEdx) {

  nhits_trd = 0;
  for (int ihit = 0; ihit < NumberOfTrdLayers; ihit++) {
    if (trd_track.Value(trd_dEdx_field_.at(ihit)) > 0) {
      dEdx.at(ihit) = trd_track.Value(trd_dEdx_field_.at(ihit));
      nhits_trd++;
    } else
      dEdx.at(ihit) = 0;
  }
}

void PidTrdMcInput::FillHistogram(TrdContainer track) {

  float mom = track.GetP();
  int nhits_trd = track.GetNhitsTrd();
  int mc_pdg = track.GetMcPdg();
  int charge = track.GetCharge();

  // Histograms for dEdx of track
  for (int imode = 0; imode < nhits_trd; imode++) {

    float dEdx = track.GetdEdxTrack(imode);
    if (dEdx <= 0) continue;

    outFile_->cd(dirname_tracks_ + "/" + dirname_nhits_.at(nhits_trd - 1) + "/reco_info");
    if (charge > 0)
      h2dEdx_p_pos_[imode]->Fill(mom, dEdx);
    if (charge < 0)
      h2dEdx_p_neg_[imode]->Fill(mom, dEdx);

    for (int ipdg = 0; ipdg < NumberOfPidsTrd - 1; ipdg++) {
      if (TMath::Abs(mc_pdg) == pid_codes_trd_.at(ipdg).first) {
        outFile_->cd(dirname_tracks_ + "/" + dirname_nhits_.at(nhits_trd - 1) + "/reco_vs_sim_info/" + pid_codes_trd_.at(ipdg).second);
        if (charge > 0)
          h2dEdx_p_pdg_[imode].at(ipdg)->Fill(mom, dEdx);
        if (charge < 0)
          h2dEdx_p_pdg_[imode].at(NumberOfPidsTrd - 1 + ipdg)->Fill(mom, dEdx);
      }
    }
  }

  // Histograms for dEdx of hits
  std::array<float, NumberOfTrdLayers> dEdx_hits_sorted = track.GetdEdxHitsSorted();

  // All hits in one plot for every particle
  outFile_->cd(dirname_hits_ + "/reco_info");
  for (int ihit = 0; ihit < nhits_trd; ihit++) {
    if (charge > 0)
      h2dEdx_hits_all_p_pos_->Fill(mom, dEdx_hits_sorted.at(ihit));
    if (charge < 0)
      h2dEdx_hits_all_p_neg_->Fill(mom, dEdx_hits_sorted.at(ihit));
  }

  for (int ipdg = 0; ipdg < NumberOfPidsTrd - 1; ipdg++) {
    if (TMath::Abs(mc_pdg) == pid_codes_trd_.at(ipdg).first) {
      outFile_->cd(dirname_hits_ + "/reco_vs_sim_info/" + pid_codes_trd_.at(ipdg).second);
      for (int ihit = 0; ihit < nhits_trd; ihit++) {
        if (charge > 0)
          h2dEdx_hits_all_p_pdg_.at(ipdg)->Fill(mom, dEdx_hits_sorted.at(ihit));
        if (charge < 0)
          h2dEdx_hits_all_p_pdg_.at(NumberOfPidsTrd - 1 + ipdg)->Fill(mom, dEdx_hits_sorted.at(ihit));
      }
    }
  }
  for (int imode = 0; imode < nhits_trd; imode++) {
    outFile_->cd(dirname_hits_ + "/" + dirname_nhits_.at(nhits_trd - 1) + "/reco_info");
    for (int ihit = 0; ihit < imode + 1; ihit++) {
      if (charge > 0)
        h2dEdx_hits_p_pos_[imode]->Fill(mom, dEdx_hits_sorted.at(ihit));
      if (charge < 0)
        h2dEdx_hits_p_neg_[imode]->Fill(mom, dEdx_hits_sorted.at(ihit));
    }

    // Hits seperated by number of hits per track & truncation mode
    for (int ipdg = 0; ipdg < NumberOfPidsTrd - 1; ipdg++) {
      if (TMath::Abs(mc_pdg) == pid_codes_trd_.at(ipdg).first) {
        outFile_->cd(dirname_hits_ + "/" + dirname_nhits_.at(nhits_trd - 1) + "/reco_vs_sim_info/" + pid_codes_trd_.at(ipdg).second);
        for (int ihit = 0; ihit < imode + 1; ihit++) {
          if (charge > 0)
            h2dEdx_hits_p_pdg_[imode].at(ipdg)->Fill(mom, dEdx_hits_sorted.at(ihit));
          if (charge < 0)
            h2dEdx_hits_p_pdg_[imode].at(NumberOfPidsTrd - 1 + ipdg)->Fill(mom, dEdx_hits_sorted.at(ihit));
        }
      }
    }
  }
}

void PidTrdMcInput::CreateMcHistograms(int nhits) {

  TString histname;

  // Histgrams for dEdx of tracks
  for (int imode = 0; imode < NumberOfTruncMode - 1; imode++) {
    outFile_->cd(dirname_tracks_ + "/" + dirname_nhits_.at(nhits - 1) + "/reco_info");

    h2dEdx_p_pos_[imode] = (TH2F*) gDirectory->Get(histnames_all_pos_.at(imode));
    h2dEdx_p_neg_[imode] = (TH2F*) gDirectory->Get(histnames_all_neg_.at(imode));

    for (int ipdg = 0; ipdg < NumberOfPidsTrd - 1; ipdg++) {
      outFile_->cd(dirname_tracks_ + "/" + dirname_nhits_.at(nhits - 1) + "/reco_vs_sim_info/" + pid_codes_trd_.at(ipdg).second);

      histname = histnames_pos_.at(imode).at(ipdg);
      h2dEdx_p_pdg_[imode].at(ipdg) = (TH2F*) gDirectory->Get(histname);

      histname = histnames_neg_.at(imode).at(ipdg);
      h2dEdx_p_pdg_[imode].at(NumberOfPidsTrd - 1 + ipdg) = (TH2F*) gDirectory->Get(histname);
    }
  }

  // Histgrams for dEdx of hits
  outFile_->cd(dirname_hits_ + "/reco_info");
  h2dEdx_hits_all_p_pos_ = (TH2F*) gDirectory->Get(histnames_hits_all_pos_);
  h2dEdx_hits_all_p_neg_ = (TH2F*) gDirectory->Get(histnames_hits_all_neg_);

  for (int ipdg = 0; ipdg < NumberOfPidsTrd - 1; ipdg++) {
    outFile_->cd(dirname_hits_ + "/reco_vs_sim_info/" + pid_codes_trd_.at(ipdg).second);

    histname = histnames_hits_pos_.at(ipdg);
    h2dEdx_hits_all_p_pdg_.at(ipdg) = (TH2F*) gDirectory->Get(histname);

    histname = histnames_hits_neg_.at(ipdg);
    h2dEdx_hits_all_p_pdg_.at(NumberOfPidsTrd - 1 + ipdg) = (TH2F*) gDirectory->Get(histname);
  }

  for (int imode = 0; imode < NumberOfTruncMode - 1; imode++) {
    outFile_->cd(dirname_hits_ + "/" + dirname_nhits_.at(nhits - 1) + "/reco_info");

    h2dEdx_hits_p_pos_[imode] = (TH2F*) gDirectory->Get(histnames_all_pos_.at(imode));
    h2dEdx_hits_p_neg_[imode] = (TH2F*) gDirectory->Get(histnames_all_neg_.at(imode));

    for (int ipdg = 0; ipdg < NumberOfPidsTrd - 1; ipdg++) {
      outFile_->cd(dirname_hits_ + "/" + dirname_nhits_.at(nhits - 1) + "/reco_vs_sim_info/" + pid_codes_trd_.at(ipdg).second);

      histname = histnames_pos_.at(imode).at(ipdg);
      h2dEdx_hits_p_pdg_[imode].at(ipdg) = (TH2F*) gDirectory->Get(histname);

      histname = histnames_neg_.at(imode).at(ipdg);
      h2dEdx_hits_p_pdg_[imode].at(NumberOfPidsTrd - 1 + ipdg) = (TH2F*) gDirectory->Get(histname);
    }
  }

  auto it = trd_container_.find(nhits);
  if (it != trd_container_.end()) {
    for (auto track : it->second) {
      FillHistogram(track);
    }
  }
  outFile_->cd();
  outFile_->Write("", TObject::kOverwrite);
}

void PidTrdMcInput::InitMcHistograms() {

  outFile_ = new TFile(outfilename_, "RECREATE");

  // Create directories
  TDirectory *directory, *directory1, *directory2, *directory3;
  for (int idir = 0; idir < 2; idir++) {
    if (idir == 0) directory = outFile_->mkdir(dirname_tracks_);
    if (idir == 1) {
      directory = outFile_->mkdir(dirname_hits_);
      directory1 = directory->mkdir("reco_info");
      directory1 = directory->mkdir("reco_vs_sim_info");
      for (int ipdg = 0; ipdg < NumberOfPidsTrd - 1; ipdg++) {
        directory2 = directory1->mkdir(pid_codes_trd_.at(ipdg).second);
      }
    }
    for (int inhits = 0; inhits < NumberOfTrdLayers; inhits++) {
      directory1 = directory->mkdir(dirname_nhits_.at(inhits));
      directory2 = directory1->mkdir("reco_info");
      directory2 = directory1->mkdir("reco_vs_sim_info");
      for (int ipdg = 0; ipdg < NumberOfPidsTrd - 1; ipdg++) {
        directory3 = directory2->mkdir(pid_codes_trd_.at(ipdg).second);
      }
    }
  }

  // Initialize histograms
  TString dirname, histname, histtitle;
  TH2F* histogram;

  for (int idir = 0; idir < 2; idir++) {

    if (idir == 0) dirname = dirname_tracks_;

    if (idir == 1) {

      dirname = dirname_hits_;

      outFile_->cd(dirname + "/reco_info");
      histname = histnames_hits_all_pos_;
      histtitle = "dEdx hits vs. p_{rec} all+";
      histogram = new TH2F(histname, histtitle, nbins_mom_, bins_mom_, nbins_dEdx_, bins_dEdx_);
      histogram->Write();

      histname = histnames_hits_all_neg_;
      histtitle = "dEdx hits vs. p_{rec} all- ";
      histogram = new TH2F(histname, histtitle, nbins_mom_, bins_mom_, nbins_dEdx_, bins_dEdx_);
      histogram->Write();

      for (int ipdg = 0; ipdg < NumberOfPidsTrd - 1; ipdg++) {
        TString particlename = pid_codes_trd_.at(ipdg).second;

        outFile_->cd(dirname + "/reco_vs_sim_info/" + pid_codes_trd_.at(ipdg).second);

        histname = histnames_hits_pos_.at(ipdg);
        histtitle = Form("dEdx hits : p_{rec} (%s+)", particlename.Data());
        histogram = new TH2F(histname, histtitle, nbins_mom_, bins_mom_, nbins_dEdx_, bins_dEdx_);
        histogram->Write();

        histname = histnames_hits_neg_.at(ipdg);
        histtitle = Form("dEdx hits : p_{rec} (%s-)", particlename.Data());
        histogram = new TH2F(histname, histtitle, nbins_mom_, bins_mom_, nbins_dEdx_, bins_dEdx_);
        histogram->Write();
      }
    }

    for (int inhits = 0; inhits < NumberOfTrdLayers; inhits++) {
      for (int imode = 0; imode < NumberOfTruncMode - 1; imode++) {

        if (imode > inhits) continue;

        outFile_->cd(dirname + "/" + dirname_nhits_.at(inhits) + "/reco_info");

        histname = histnames_all_pos_.at(imode);
        histtitle = "dEdx vs. p_{rec} all+ " + to_string(inhits + 1) + " Trd Hits - " + histtitle_mode_.at(imode);
        histogram = new TH2F(histname, histtitle, nbins_mom_, bins_mom_, nbins_dEdx_, bins_dEdx_);
        histogram->Write();

        histname = histnames_all_neg_.at(imode);
        histtitle = "dEdx vs. p_{rec} all- " + to_string(inhits + 1) + " Trd Hits - " + histtitle_mode_.at(imode);
        histogram = new TH2F(histname, histtitle, nbins_mom_, bins_mom_, nbins_dEdx_, bins_dEdx_);
        histogram->Write();

        for (int ipdg = 0; ipdg < NumberOfPidsTrd - 1; ipdg++) {
          TString particlename = pid_codes_trd_.at(ipdg).second;

          outFile_->cd(dirname + "/" + dirname_nhits_.at(inhits) + "/reco_vs_sim_info/" + pid_codes_trd_.at(ipdg).second);

          histname = histnames_pos_.at(imode).at(ipdg);
          histtitle = Form("dEdx : p_{rec} (%s+) %d Trd Hits - %s", particlename.Data(), inhits + 1, histtitle_mode_.at(imode).Data());
          histogram = new TH2F(histname, histtitle, nbins_mom_, bins_mom_, nbins_dEdx_, bins_dEdx_);
          histogram->Write();

          histname = histnames_neg_.at(imode).at(ipdg);
          histtitle = Form("dEdx : p_{rec} (%s-) %d Trd Hits - %s", particlename.Data(), inhits + 1, histtitle_mode_.at(imode).Data());
          histogram = new TH2F(histname, histtitle, nbins_mom_, bins_mom_, nbins_dEdx_, bins_dEdx_);
          histogram->Write();
        }
      }
    }
  }
}

void PidTrdMcInput::OpenMcHistograms() {
  outFile_ = new TFile(outfilename_, "UPDATE");
}

void PidTrdMcInput::Init() {
  auto man = TaskManager::GetInstance();
  auto chain = man->GetChain();

  chain->InitPointersToBranches({});
  rec_tracks_ = chain->GetBranchObject(rec_tracks_name_);
  trd_tracks_ = chain->GetBranchObject(trd_tracks_name_);
  sim_tracks_ = chain->GetBranchObject(sim_tracks_name_);
  rec_to_trd_ = chain->GetMatching(rec_tracks_name_, trd_tracks_name_);
  rec_to_sim_ = chain->GetMatching(rec_tracks_name_, sim_tracks_name_);

  for (int i = 0; i < NumberOfTrdLayers; i++)
    trd_dEdx_field_.push_back(trd_tracks_.GetField(("energy_loss_" + to_string(i)).c_str()));

  Float_t bw_mom, bw_dEdx;

  if (update_mchisto_ == kFALSE) {
    bw_mom = BwMom;
    nbins_mom_ = NbinsMom;
    bins_mom_[0] = 0.;
    for (int i = 0; i < nbins_mom_; i++)
      bins_mom_[i + 1] = TMath::Power(10, -1.0 + i * bw_mom);

    bw_dEdx = BwdEdx;
    nbins_dEdx_ = NbinsdEdx;
    bins_dEdx_[0] = 0.;
    for (int i = 0; i < nbins_dEdx_; i++)
      bins_dEdx_[i + 1] = TMath::Power(10, -1.0 + i * bw_dEdx);
  }

  for (int imode = 0; imode < NumberOfTruncMode - 1; imode++) {
    histnames_all_pos_.at(imode) = "h2dEdx_p_pos_" + to_string(imode);
    histnames_all_neg_.at(imode) = "h2dEdx_p_neg_" + to_string(imode);
    for (int ipdg = 0; ipdg < NumberOfPidsTrd - 1; ipdg++) {
      TString particlename = pid_codes_trd_.at(ipdg).second;
      histnames_pos_.at(imode).at(ipdg) = "h2dEdx_p_" + particlename + "_pos_" + to_string(imode);
      histnames_neg_.at(imode).at(ipdg) = "h2dEdx_p_" + particlename + "_neg_" + to_string(imode);
    }

    histnames_hits_all_pos_ = "h2dEdx_p_pos";
    histnames_hits_all_neg_ = "h2dEdx_p_neg";
    for (int ipdg = 0; ipdg < NumberOfPidsTrd - 1; ipdg++) {
      TString particlename = pid_codes_trd_.at(ipdg).second;
      histnames_hits_pos_.at(ipdg) = "h2dEdx_p_" + particlename + "_pos";
      histnames_hits_neg_.at(ipdg) = "h2dEdx_p_" + particlename + "_neg";
    }
  }

  if (update_mchisto_ == kTRUE)
    OpenMcHistograms();
  else
    InitMcHistograms();
}

void PidTrdMcInput::Exec() {
  trd_container_.clear();

  for (int itrack = 0; itrack < rec_tracks_.size(); ++itrack) {
    const auto& rec_track = rec_tracks_[itrack];
    int itrd = rec_to_trd_->GetMatch(rec_track.GetId());
    if (itrd < 0) continue;
    const auto& trd_track = trd_tracks_[itrd];

    float mom = GetMomentum(trd_track);
    float pT = GetPt(trd_track);
    int mc_pdg = GetMcPdg(rec_track);
    int charge = GetCharge(rec_track);
    array<float, NumberOfTrdLayers> dEdx_hits = {0.0, 0.0, 0.0, 0.0};
    int nhits_trd = 0;

    GetEnergyLossHits(trd_track, nhits_trd, dEdx_hits);

    if (nhits_trd < 1) continue;

    TrdContainer trdtrack(mom, pT, charge, nhits_trd, dEdx_hits, mc_pdg);
    trdtrack.ScaleEnergyLossLength();
    trdtrack.CalculateEnergyLossTrackAllModes();

    auto it = trd_container_.find(nhits_trd);
    if (it != trd_container_.end()) {
      it->second.emplace_back(trdtrack);
    } else {
      trd_container_[nhits_trd] = {static_cast<TrdContainer>(trdtrack)};
    }
  }

  for (int inhits = 0; inhits < NumberOfTrdLayers; inhits++) {
    CreateMcHistograms(inhits + 1);
  }
}

void PidTrdMcInput::Finish() {
  outFile_->Close();
}
