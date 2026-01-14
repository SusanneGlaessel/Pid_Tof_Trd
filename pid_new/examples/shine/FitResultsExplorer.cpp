//
// Created by eugene on 29/05/2020.
//

#include <TFile.h>
#include <TKey.h>
#include <TTree.h>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include <ROOT/RDataFrame.hxx>

using namespace std;

struct BinKinematics_t {
  bool is_initialized{false};

  double p_lo{0.};
  double p_hi{0.};

  double pt_lo{0.};
  double pt_hi{0.};

  double phi_lo{0.};
  double phi_hi{0.};
};

struct Bin_t {
  string name{""};
  TDirectory* bin_dir{nullptr};

  TTree* track_tree{nullptr};
  TDirectory* qa_dir{nullptr};

  std::map<string, double> fit_parameters;

  BinKinematics_t bin_kin_info;
};

void readBinKinInfo(Bin_t& bin) {
  assert(bin.track_tree);

  ROOT::RDataFrame rdf(*bin.track_tree);

  auto p_lo = rdf.Min("p");
  auto p_hi = rdf.Max("p");
  auto pt_lo = rdf.Min("pt");
  auto pt_hi = rdf.Max("pt");
  auto phi_lo = rdf.Min("phi");
  auto phi_hi = rdf.Max("phi");

  bin.bin_kin_info.p_lo = *p_lo;
  bin.bin_kin_info.p_hi = *p_hi;
  bin.bin_kin_info.pt_lo = *pt_lo;
  bin.bin_kin_info.pt_hi = *pt_hi;
  bin.bin_kin_info.phi_lo = *phi_lo;
  bin.bin_kin_info.phi_hi = *phi_hi;
  bin.bin_kin_info.is_initialized = true;
}

vector<Bin_t> lookupBins(TDirectory* dir) {

  auto keys_list = dir->GetListOfKeys();
  vector<Bin_t> lookup_result;

  for (auto key_tobj : *keys_list) {
    auto key = static_cast<TKey*>(key_tobj);
    string key_name{key_tobj->GetName()};

    if (
        key_name.substr(0, strlen("bin_")) == "bin_" && key->ReadObject<TDirectory>() != nullptr) {
      std::cout << "Looking up " << key_name << std::endl;
      Bin_t bin;
      bin.name = key_name;
      bin.bin_dir = key->ReadObject<TDirectory>();

      bin.track_tree = bin.bin_dir->Get<TTree>("track_tree");
      bin.qa_dir = bin.bin_dir->Get<TDirectory>("qa");

      if (bin.qa_dir) {
        auto fit_parameters_ptr = bin.qa_dir->Get<map<string, double>>("fit_parameters");

        if (fit_parameters_ptr) {
          bin.fit_parameters = *fit_parameters_ptr;
        }
      }

      readBinKinInfo(bin);
      lookup_result.push_back(bin);
    }
  }

  return lookup_result;
}

vector<string> getMapKeys(const map<string, double>& map) {
  vector<string> result;
  transform(map.begin(), map.end(), back_inserter(result), [](const pair<string, double>& e) { return e.first; });
  return result;
}

int main() {

  string input_file_name{"/home/eugene/STORAGE/data/pid/vtx_tack_dumps/pbpb_13agev_16_026/tracks_binned.root"};
  string charge_str{"negative"};

  TFile input_file(input_file_name.c_str(), "read");
  auto tgt_dir = input_file.Get<TDirectory>(charge_str.c_str());

  auto bins_list = lookupBins(tgt_dir);

  BinKinematics_t bin_kin_info;

  TFile output_file("fit_results_explorer.root", "recreate");
  /* init output tree */
  map<string, pair<TTree*, double>> fit_param_trees;
  for (auto& bin : bins_list) {
    bin_kin_info = bin.bin_kin_info;

    for (auto& p : bin.fit_parameters) {
      auto par_name = p.first;

      auto par_tree_it = fit_param_trees.find(par_name);
      if (fit_param_trees.end() == par_tree_it) {
        auto par_tree = new TTree(par_name.c_str(), "");
        par_tree->Branch("p_lo", &bin_kin_info.p_lo);
        par_tree->Branch("p_hi", &bin_kin_info.p_hi);
        par_tree->Branch("pt_lo", &bin_kin_info.pt_lo);
        par_tree->Branch("pt_hi", &bin_kin_info.pt_hi);
        par_tree->Branch("phi_lo", &bin_kin_info.phi_lo);
        par_tree->Branch("phi_hi", &bin_kin_info.phi_hi);

        par_tree_it = fit_param_trees.emplace(par_name, make_pair(par_tree, 0.)).first;
        par_tree->Branch("par_value", &(par_tree_it->second.second));
      }

      auto par_value = p.second;
      par_tree_it->second.second = par_value;
      par_tree_it->second.first->Fill();
    }
  }

  for (auto& t : fit_param_trees) {
    t.second.first->Write();
  }

  return 0;
}
