#include "PidFiller.hpp"

using std::string;
using namespace AnalysisTree;

PidFiller::PidFiller(const std::string& pid_file_name_tof, const std::string& pid_file_name_trd, const std::string& getter_name_tof, const std::string& getter_name_trd, int pid_mode) {
  std::cout << "PidFiller::PidFiller" << std::endl;

  std::unique_ptr<TFile> pid_file_tof; std::unique_ptr<TFile> pid_file_trd;

  if (pid_mode == 0) {
    is_run_pid_tof_ = true;
    is_run_pid_trd_ = true;
  }

  if (pid_mode == 1)
    is_run_pid_tof_ = true;

  if (pid_mode == 2)
    is_run_pid_trd_ = true;
  if (is_run_pid_tof_ == true) {
    pid_file_tof = std::unique_ptr<TFile>(TFile::Open(pid_file_name_tof.c_str(), "read"));
    if ((!pid_file_tof) || (pid_file_tof->IsZombie())) {
      throw std::runtime_error("No file or file is zombie: " + pid_file_name_tof);
    }
    getter_tof_ = pid_file_tof->Get<Pid::GetterTof>(getter_name_tof.c_str());
    if (getter_tof_ == nullptr) {
      throw std::runtime_error("PID Tof Getter is nullptr: " + getter_name_tof);
    }
  }

  if (is_run_pid_trd_ == true) {
    pid_file_trd = std::unique_ptr<TFile>(TFile::Open(pid_file_name_trd.c_str(), "read"));
    if ((!pid_file_trd) || (pid_file_trd->IsZombie())) {
      throw std::runtime_error("No file or file is zombie: " + pid_file_name_trd);
    }
    getter_trd_ = pid_file_trd->Get<PidTrd::GetterTrd>(getter_name_trd.c_str());
    if (getter_trd_ == nullptr) {
      throw std::runtime_error("PID Trd Getter is nullptr: " + getter_name_trd);
    }    
  }
}

int PidFiller::signum(int x) const {
  if (x > 0)
    return 1;
  else if (x == 0)
    return 0;
  else
    return -1;
}

 float PidFiller::GetMomentumTrd(const AnalysisTree::BranchChannel& trd_track) {
  float momentum = 0.0; 
  if (TMath::Abs(trd_track.Value(trd_tracks_.GetField("p")) > 0))
    momentum = trd_track.Value(trd_tracks_.GetField("p"));
  else if (TMath::Abs(trd_track.Value(trd_tracks_.GetField("p_out")) > 0))
    momentum = trd_track.Value(trd_tracks_.GetField("p_out"));
  else
    Warning("Exec", "Could not assign any momentum to the track, use p=0.");
    
  return momentum;
}

float PidFiller::GetPtTrd(const AnalysisTree::BranchChannel& trd_track) { 
  float pT = 0; 
  if (TMath::Abs(trd_track.Value(trd_tracks_.GetField("pT")) > 0))
    pT = trd_track.Value(trd_tracks_.GetField("pT"));
  else if (TMath::Abs(trd_track.Value(trd_tracks_.GetField("pT_out")) > 0))
    pT = trd_track.Value(trd_tracks_.GetField("pT_out"));
  else 
    Warning("Exec", "Could not assign any pT to the track, use pT=0.");
  return pT;
}

int PidFiller::GetCharge(const AnalysisTree::BranchChannel& rec_track) {
  return rec_track.Value(rec_tracks_.GetField("q"));
}

void PidFiller::GetEnergyLossHitsTrd(const AnalysisTree::BranchChannel& trd_track, int &nhits_trd, std::array<float, NumberOfTrdLayers> &dEdx) {
  
  nhits_trd = 0;
  
  for (int ihit = 0 ; ihit < NumberOfTrdLayers; ihit++) {
    if ( trd_track.Value(trd_dEdx_field_.at(ihit)) > 0) {
      dEdx.at(ihit) = trd_track.Value(trd_dEdx_field_.at(ihit));   
      nhits_trd ++;
    }
    else
      dEdx.at(ihit) = 0;  
  }
}

bool PidFiller::IsRichElectron(const AnalysisTree::BranchChannel& rec_track) {
  bool isElectron = false;
  int i_rich = rec_to_rich_->GetMatch(rec_track.GetId());
  if (i_rich > 0) {
    const auto& rich_ring = rich_rings_[i_rich];
    float Aaxis = rich_ring.Value(rich_rings_.GetField("axis_a"));
    float Baxis = rich_ring.Value(rich_rings_.GetField("axis_b"));
    Double_t dist  = 0;  // richRing->GetDistance();
	    
    Float_t mom = GetMomentumTrd(rec_track);
    Double_t MeanA    = 4.95;
    Double_t MeanB    = 4.54;
    Double_t RmsA     = 0.30;
    Double_t RmsB     = 0.22;
    Double_t RmsCoeff = 3.5;
    Double_t DistCut  = 1.;
	    
    if (mom < 5.) {
      if (fabs(Aaxis - MeanA) < RmsCoeff * RmsA && fabs(Baxis - MeanB) < RmsCoeff * RmsB && dist < DistCut)
	isElectron = true;
    }
    else {                     
      ///2 sigma
      Double_t polAaxis = 5.64791 - 4.24077 / (mom - 3.65494);
      Double_t polBaxis = 5.41106 - 4.49902 / (mom - 3.52450);
      if (Aaxis < (MeanA + RmsCoeff * RmsA) && Aaxis > polAaxis && Baxis < (MeanB + RmsCoeff * RmsB) && Baxis > polBaxis && dist < DistCut)              
	isElectron = true;
    }      
  }
  return isElectron;
}
 
void PidFiller::Init() {
  // Settings for Trd Pid
  if (is_run_pid_trd_ == true) {
    getter_trd_->SetMinHits(nhits_min_);
    getter_trd_->SetTruncationMode(trunc_mode_);
    getter_trd_->SetProbabiltyMode(prob_mode_);
  }

  auto man = TaskManager::GetInstance();
  auto chain = man->GetChain();
  
  rec_tracks_ = chain->GetBranch(rec_tracks_name_);
  in_branches_.emplace(rec_tracks_name_);
  auto conf = rec_tracks_.GetConfig().Clone(out_branch_name_, AnalysisTree::DetType::kParticle);

  // Set and configure branches for Tof Pid
  if (is_run_pid_tof_ == true) {
    tof_hits_ = chain->GetBranch(tof_hits_name_);
    rec_to_tof_ = chain->GetMatchPointers().find({rec_tracks_name_ + "2" + tof_hits_name_})->second;
  
    in_branches_.emplace(tof_hits_name_);

    std::vector<std::string> names{};
    for (const auto& pid : pid_codes_tof_) {
    names.push_back("prob_" + pid.second);
    }
    conf.AddFields<float>(names, "probability to be proton, pion, kaon etc");
  }

  // Set and configure branches for Trd Pid
  if (is_run_pid_trd_ == true) {
  
    trd_tracks_ = chain->GetBranchObject(trd_tracks_name_);
    rich_rings_ = chain->GetBranchObject(rich_rings_name_);
    rec_to_trd_ = chain->GetMatching(rec_tracks_name_, trd_tracks_name_);
    rec_to_rich_ = chain->GetMatching(rec_tracks_name_, rich_rings_name_);

    in_branches_.emplace(trd_tracks_name_);
    in_branches_.emplace(rich_rings_name_);

    conf.AddField<bool>("electron_rich", "rich electron hypothesis");
    conf.AddField<int>("pid_trd", "trd pid hypothesis");
    conf.AddField<int>("nhits_trd_eloss", "number of trd hits dEdx > 0");
    
    std::vector<std::string> names{};
    for (const auto& pid : pid_codes_trd_) {
      string name = Form("prob_trd_%s", pid.second.Data());
      names.push_back(Form("prob_trd_%s", pid.second.Data()));
    }
    conf.AddFields<float>(names, "trd probability to be proton, pion, kaon etc");
  }

  ana_tracks_ = Branch(conf);
  ana_tracks_.SetMutable();
  ana_tracks_.Freeze();

  rec_tracks_.Freeze();

  man->AddBranch(&ana_tracks_);

  // Add matching for ouput
  int imatch{0}; 
  std::vector<std::string> match_br;
  match_br.push_back("SimParticles");
  
  if (is_run_pid_tof_ == true)
    match_br.push_back("TofHits");

  if (is_run_pid_trd_ == true) {
    match_br.push_back("TrdTracks");
    match_br.push_back("RichRings");
  }
  
  out_matches_.assign(match_br.size(), nullptr);
  
  for (const auto& br : match_br) {
    in_matches_.emplace_back(chain->GetMatchPointers().find({rec_tracks_name_ + "2" + br})->second);
    man->AddMatching(out_branch_name_, br, out_matches_.at(imatch));
    imatch++;
  }

  // Set input variables for Tof Pid
  if (is_run_pid_tof_ == true) {
    qp_tof_field_ = tof_hits_.GetField("qp_tof");
    q_field_ = rec_tracks_.GetField("q");
    mass2_field_ = tof_hits_.GetField("mass2");
    pid_tof_field_ = ana_tracks_.GetField("pid");
    for (const auto& pid : pid_codes_tof_) {
      prob_tof_field_.push_back(ana_tracks_.GetField("prob_" + pid.second));
    }
  }

  // Set input variables for Trd Pid
  if (is_run_pid_trd_ == true) {
    nhits_trd_eloss_field_ = ana_tracks_.GetField("nhits_trd_eloss");
    pid_trd_field_ = ana_tracks_.GetField("pid_trd");
    el_rich_field_ = ana_tracks_.GetField("electron_rich");
    
    for (const auto& pid : pid_codes_trd_) {
      prob_trd_field_.push_back(ana_tracks_.GetField(Form("prob_trd_%s", pid.second.Data())));
    }  
    for (int i = 0; i < NumberOfTrdLayers; i++)
      trd_dEdx_field_.push_back(trd_tracks_.GetField(("energy_loss_" + std::to_string(i)).c_str()));
    
  }
}
 
void PidFiller::Exec() {
  ana_tracks_.ClearChannels();

  for (int i = 0; i < rec_tracks_.size(); ++i) {
    const auto& track = rec_tracks_[i];
    auto particle = ana_tracks_.NewChannel();
    particle.CopyContent(track);

    // Fill in Tof pid
    if (is_run_pid_tof_ == true) {
      auto q = track[q_field_];

      auto hit_id = rec_to_tof_->GetMatch(i);
      if (hit_id >= 0) {
	const auto& tof_hit = tof_hits_[hit_id];
	auto pq = tof_hit[qp_tof_field_];
	auto m2 = tof_hit[mass2_field_];

	auto pid = getter_tof_->GetPid(pq, m2, 0.5);
	if (pid == 1 && q < 0) {
	  pid = -1;
	}

	particle.SetValue(pid_tof_field_, pid);
	auto prob = getter_tof_->GetBayesianProbability(pq, m2);

	int specie{0};
	for (const auto& pdg : pid_codes_tof_) {
	  particle.SetValue(prob_tof_field_.at(specie++), prob[pdg.first * signum(q)]);// Think what to do in case of electrons and muons
	}
      } else {
	int specie{0};
	for (const auto& pdg : pid_codes_tof_) {
	  particle.SetValue(prob_tof_field_.at(specie++), -1.f);
	}
	particle.SetValue(pid_tof_field_, 2 * signum(q));
      }
    }

    // Fill in Trd pid
    if (is_run_pid_trd_ == true) {

      int itrd = rec_to_trd_->GetMatch(i);
      int nhits_trd = 0;
    
      if (itrd > -1) {
	const auto& trd_track = trd_tracks_[itrd];
	float mom = GetMomentumTrd(trd_track);
	float pT = GetPtTrd(trd_track);
	int charge = GetCharge(track);
      
	std::array<float, NumberOfTrdLayers> dEdx_hits = {0.0, 0.0, 0.0, 0.0};
	GetEnergyLossHitsTrd(trd_track, nhits_trd, dEdx_hits);

	TrdContainer trdtrack(mom, pT, charge, nhits_trd, dEdx_hits);
	trdtrack.ScaleEnergyLossLength();     
      
	int pidtrd_hypo;
	std::map<int, float> prob_trd;
	if (nhits_trd >= nhits_min_) {
	  if (prob_mode_ == 0) 
	    prob_trd = getter_trd_->GetTrdProbabilities(trdtrack);
	  if (prob_mode_ == 1) 
	    prob_trd = getter_trd_->GetTrdProbabilitiesMulti(trdtrack);

	  int specie{0};
	  for (const auto& pdg : pid_codes_trd_) 
	    particle.SetValue(prob_trd_field_.at(specie++), prob_trd[pdg.first]);
		
	  pidtrd_hypo = getter_trd_->GetTrdPid(prob_trd, purity_, charge);
	  particle.SetValue(pid_trd_field_, pidtrd_hypo);
	}
	else { 
	  particle.SetValue(pid_trd_field_, -2);
	  int specie{0};
	  for (const auto& pdg : pid_codes_trd_)
	    particle.SetValue(prob_trd_field_.at(specie++), -1.f);
	}
      }
      else { 
	particle.SetValue(pid_trd_field_, -2);
	int specie{0};
	for (const auto& pdg : pid_codes_trd_)
	  particle.SetValue(prob_trd_field_.at(specie++), -1.f);
      }

      particle.SetValue(nhits_trd_eloss_field_, nhits_trd);

      bool isElectron = IsRichElectron(track);
      particle.SetValue(el_rich_field_, isElectron);
    }
  }
  int i{0};
  for (auto& match : out_matches_) {
    auto m1 = in_matches_[i]->GetMatches(false);
    auto m2 = in_matches_[i]->GetMatches(true);
    match->SetMatches(m1, m2);
    i++;
  }
}
