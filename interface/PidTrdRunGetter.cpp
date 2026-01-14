#include "PidTrdRunGetter.hpp"
#include "ConstantsTrd.h"
using std::to_string;

void PidTrdRunGetter::CalculateProbabilities(int nhits) {

  TH2F* h2pos;
  TH2F* h2neg;
  TH2F* h2tmp;
  TH2F* h2prob;
  TString histname;
  PidTrd::ParticleProb particleprob;

  for (int truncmode = 0; truncmode < nhits + 1; truncmode++) {
    inFile_->cd(dirname_tracks_ + "/" + dirname_nhits_.at(nhits) + "/reco_info");
    h2pos = (TH2F*) gDirectory->Get(histnames_all_pos_.at(truncmode));
    h2neg = (TH2F*) gDirectory->Get(histnames_all_neg_.at(truncmode));

    for (int ipdg = 0; ipdg < NumberOfPidsTrd - 1; ipdg++) {
      TString particlename = pid_codes_trd_.at(ipdg).second.Data();
      int charge = 1;// charge = 1: positive particles

      int probmode = 0;// total probability - probability based on particle multiplicites
      inFile_->cd(dirname_tracks_ + "/" + dirname_nhits_.at(nhits) + "/reco_vs_sim_info/" + pid_codes_trd_.at(ipdg).second);
      h2tmp = (TH2F*) gDirectory->Get(histnames_pos_.at(truncmode).at(ipdg));
      h2prob = (TH2F*) h2tmp->Clone(histnames_prob_pos_.at(truncmode).at(ipdg));
      CalculateProbabilitiesTot(h2pos, h2tmp, h2prob);
      h2prob->SetTitle(Form("dEdx : p_{rec} probability (T) (%s+) %d Trd Hits - %s", pid_codes_trd_.at(ipdg).second.Data(), nhits + 1, histtitle_mode_.at(truncmode).Data()));
      particleprob.Update(ipdg, charge, nhits, truncmode, probmode, h2prob);
      getter_trd_.AddParticleProb(particleprob);
      if (write_mchistograms_out_ == kTRUE)
        outFile_->cd(dirname_nhits_.at(nhits) + "/" + particlename + "/" + probname_.at(0)) && h2prob->Write("", TObject::kOverwrite);

      probmode = 1;// likelihood - probability based on dEdx-distribution of particle
      inFile_->cd(dirname_hits_ + "/" + dirname_nhits_.at(nhits) + "/reco_vs_sim_info/" + pid_codes_trd_.at(ipdg).second);
      h2tmp = (TH2F*) gDirectory->Get(histnames_pos_.at(truncmode).at(ipdg));
      h2prob = (TH2F*) h2tmp->Clone(histnames_prob_pos_.at(truncmode).at(ipdg));
      CalculateProbabilitiesLike(h2tmp, h2prob);
      h2prob->SetTitle(Form("dEdx : p_{rec} probability (L) (%s+) %d Trd Hits - %s", pid_codes_trd_.at(ipdg).second.Data(), nhits + 1, histtitle_mode_.at(truncmode).Data()));
      particleprob.Update(ipdg, charge, nhits, truncmode, probmode, h2prob);
      getter_trd_.AddParticleProb(particleprob);
      if (write_mchistograms_out_ == kTRUE)
        outFile_->cd(dirname_nhits_.at(nhits) + "/" + particlename + "/" + probname_.at(1)) && h2prob->Write("", TObject::kOverwrite);

      charge = -1;// charge = -1: negative particles
      probmode = 0;
      inFile_->cd(dirname_tracks_ + "/" + dirname_nhits_.at(nhits) + "/reco_vs_sim_info/" + pid_codes_trd_.at(ipdg).second);
      h2tmp = (TH2F*) gDirectory->Get(histnames_neg_.at(truncmode).at(ipdg));
      h2prob = (TH2F*) h2tmp->Clone(histnames_prob_neg_.at(truncmode).at(ipdg));
      CalculateProbabilitiesTot(h2neg, h2tmp, h2prob);
      h2prob->SetTitle(Form("dEdx : p_{rec} probability (T) (%s-) %d Trd Hits - %s", pid_codes_trd_.at(ipdg).second.Data(), nhits + 1, histtitle_mode_.at(truncmode).Data()));
      particleprob.Update(ipdg, charge, nhits, truncmode, probmode, h2prob);
      getter_trd_.AddParticleProb(particleprob);
      if (write_mchistograms_out_ == kTRUE)
        outFile_->cd(dirname_nhits_.at(nhits) + "/" + particlename + "/" + probname_.at(0)) && h2prob->Write("", TObject::kOverwrite);

      probmode = 1;
      inFile_->cd(dirname_hits_ + "/" + dirname_nhits_.at(nhits) + "/reco_vs_sim_info/" + pid_codes_trd_.at(ipdg).second);
      h2tmp = (TH2F*) gDirectory->Get(histnames_neg_.at(truncmode).at(ipdg));
      h2prob = (TH2F*) h2tmp->Clone(histnames_prob_neg_.at(truncmode).at(ipdg));
      CalculateProbabilitiesLike(h2tmp, h2prob);
      h2prob->SetTitle(Form("dEdx : p_{rec} probability (L) (%s-) %d Trd Hits - %s", pid_codes_trd_.at(ipdg).second.Data(), nhits + 1, histtitle_mode_.at(truncmode).Data()));
      particleprob.Update(ipdg, charge, nhits, truncmode, probmode, h2prob);
      getter_trd_.AddParticleProb(particleprob);
      if (write_mchistograms_out_ == kTRUE)
        outFile_->cd(dirname_nhits_.at(nhits) + "/" + particlename + "/" + probname_.at(1)) && h2prob->Write("", TObject::kOverwrite);
    }
  }
  h2tmp->Delete();
}

void PidTrdRunGetter::CalculateProbabilitiesTot(TH2F* h2all, TH2F* h2tmp, TH2F*& h2prob) {
  h2prob->Divide(h2all);
}

void PidTrdRunGetter::CalculateProbabilitiesLike(TH2F* h2tmp, TH2F*& h2prob) {
  for (Int_t ibinx = 1; ibinx < h2tmp->GetNbinsX() + 1; ibinx++) {
    Double_t sumy = 0.;
    for (Int_t ibiny = 1; ibiny < h2tmp->GetNbinsY() + 1; ibiny++)
      sumy += h2tmp->GetBinContent(ibinx, ibiny);
    for (Int_t ibiny = 1; ibiny < h2tmp->GetNbinsY() + 1; ibiny++)
      if (sumy > 0) h2prob->SetBinContent(ibinx, ibiny, h2tmp->GetBinContent(ibinx, ibiny) / sumy);
  }
}

void PidTrdRunGetter::InitMcHistogramsOut() {

  TDirectory *directory, *directory1, *directory2;

  for (int imode = 0; imode < NumberOfTruncMode - 1; imode++) {
    for (int ipdg = 0; ipdg < NumberOfPidsTrd - 1; ipdg++) {
      TString particlename = pid_codes_trd_.at(ipdg).second;
      histnames_prob_pos_.at(imode).at(ipdg) = "h2dEdx_p_prob_" + particlename + "_pos_" + to_string(imode);
      histnames_prob_neg_.at(imode).at(ipdg) = "h2dEdx_p_prob_" + particlename + "_neg_" + to_string(imode);
    }
  }

  for (int inhits = 0; inhits < NumberOfTrdLayers; inhits++) {
    directory = outFile_->mkdir(dirname_nhits_.at(inhits));
    for (int ipdg = 0; ipdg < NumberOfPidsTrd - 1; ipdg++) {
      outFile_->cd();
      directory1 = directory->mkdir(pid_codes_trd_.at(ipdg).second);
      directory2 = directory1->mkdir(probname_.at(0));
      directory2 = directory1->mkdir(probname_.at(1));
    }
  }
}

void PidTrdRunGetter::OpenMcHistograms() {

  for (int imode = 0; imode < NumberOfTruncMode - 1; imode++) {
    histnames_all_pos_.at(imode) = "h2dEdx_p_pos_" + to_string(imode);
    histnames_all_neg_.at(imode) = "h2dEdx_p_neg_" + to_string(imode);
    for (int ipdg = 0; ipdg < NumberOfPidsTrd - 1; ipdg++) {
      TString particlename = pid_codes_trd_.at(ipdg).second;
      histnames_pos_.at(imode).at(ipdg) = "h2dEdx_p_" + particlename + "_pos_" + to_string(imode);
      histnames_neg_.at(imode).at(ipdg) = "h2dEdx_p_" + particlename + "_neg_" + to_string(imode);
    }
  }

  histnames_hits_all_pos_ = "h2dEdx_p_pos";
  histnames_hits_all_neg_ = "h2dEdx_p_neg";
  for (int ipdg = 0; ipdg < NumberOfPidsTrd - 1; ipdg++) {
    TString particlename = pid_codes_trd_.at(ipdg).second;
    histnames_hits_pos_.at(ipdg) = "h2dEdx_p_" + particlename + "_pos";
    histnames_hits_neg_.at(ipdg) = "h2dEdx_p_" + particlename + "_neg";
  }
}

void PidTrdRunGetter::Init() {

  inFile_ = new TFile(mcfile_name_, "READ");
  if (!inFile_ || !inFile_->IsOpen()) {
    throw std::runtime_error("Could not open input file: " + mcfile_name_);
  }

  OpenMcHistograms();

  if (write_mchistograms_out_ == kTRUE) {
    TString name = mcfile_name_;
    name.Replace(name.Index(".root"), 5, "");
    TString outfilename = Form("%s_probabilities.root", name.Data());
    outFile_ = new TFile(outfilename, "RECREATE");
    InitMcHistogramsOut();
  }
}

void PidTrdRunGetter::Exec() {

  for (int inhits = 0; inhits < NumberOfTrdLayers; inhits++)
    CalculateProbabilities(inhits);

  std::unique_ptr<TFile> outFile_getter{TFile::Open(getter_file_, "recreate")};
  outFile_getter->WriteObject(&getter_trd_, getter_name_);
  outFile_getter->Close();

  if (write_mchistograms_out_ == kTRUE) outFile_->Close();
}

void PidTrdRunGetter::Finish() {
  inFile_->Close();
  if (write_mchistograms_out_ == kTRUE) outFile_->Close();
}
