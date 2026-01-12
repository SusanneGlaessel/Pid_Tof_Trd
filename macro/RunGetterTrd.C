#ifndef __CLING__

#include "GetterTrd.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TH2F.h>
#include <TRandom.h>
#include <map>
#include <vector>

void RunGetterTrd(TString InputFile);

int main(int argc, char** argv) {
  if (argc == 2) RunGetter(argv[1]);
  else
    RunGetter("pid_getter_trd.root");
  return 1;
}

#endif// __CLING__

void RunGetterTrd(TString InputFile = "pid_getter_trd.root") {

  TFile* f2{TFile::Open(InputFile)};
  std::unique_ptr<PidTrd::GetterTrd> getter{(PidTrd::GetterTrd*) f2->Get("pid_getter_trd")};

  std::unique_ptr<TRandom> r{new TRandom};

  std::vector<int> Pids = {PidTrdParticles::kBgPos, PidTrdParticles::kPionPos, PidTrdParticles::kKaonPos, PidTrdParticles::kProton, PidTrdParticles::kDeutron, PidTrdParticles::kTriton, PidTrdParticles::kHe3};
    
  std::map<int, TH2F*> h2Purity;
  for (auto pid : Pids) {
    h2Purity.emplace(pid, new TH2F(Form("h2Purity_%d", pid), Form("Purity for %d;p (GeV/#it{c});dE/dx (keV/cm)", pid), 1000, 0., 100, 1000, 0, 200));
  }
  TAxis* xAxis = h2Purity.at(PidTrdParticles::kBgPos)->GetXaxis();
  TAxis* yAxis = h2Purity.at(PidTrdParticles::kBgPos)->GetYaxis();
  uint nBinsX = xAxis->GetNbins();
  uint nBinsY = yAxis->GetNbins();
  int nhits = 2;
  int charge = 1;
  for (uint binX = 1; binX <= nBinsX; binX++) {
    float mom = xAxis->GetBinCenter(binX);
    for (uint binY = 1; binY <= nBinsY; binY++) {
      float dEdx = yAxis->GetBinCenter(binY);
      auto prob = getter->GetTrdProbabilitiesTotalAverage(mom, dEdx, charge, nhits);
      for (auto pid : Pids) 
        h2Purity.at(pid)->Fill(mom, dEdx, prob.at(pid));
    }
  }
  TCanvas* c = new TCanvas();
  c->Divide(3, 3);
  int i = 0;
  for (auto pid : Pids) {
    c->cd(++i);
    gPad->SetLogz();
    h2Purity.at(pid)->Draw("colz");
    h2Purity.at(pid)->GetYaxis()->SetTitleOffset(1.4);
    //h2Purity.at(pid)->SetMinimum(0.5);
  }
  
  // Calcluation of probabilities for test tracks
  
  for (int i=0; i<10; ++i) // Create test tracks
    {
      const float pT = 0.0;
      const float mom = r->Uniform(0.0, 30.0);
      const int charge = 1;
      
      std::array<float, NumberOfTrdLayers> dEdx_hits = {0.0, 0.0, 0.0, 0.0};
      float dEdx = 0.0;
      int nhits_trd = 0;
      for (int i = 0; i < NumberOfTrdLayers; i++) {
	dEdx_hits.at(i) = r->Uniform(1.5,40);
	if (dEdx_hits.at(i) > 0) {
	  dEdx += dEdx_hits.at(i);
	  nhits_trd ++;
	}
      }	

      dEdx /= nhits_trd;

      std::cout << std::endl;
      std::cout << "-- Track " << i << " -- " <<std::endl;
      std::cout << "momentum = " << mom << std::endl;
      std::cout << "dEdx (hits) = "<< dEdx_hits.at(0) << ", " <<dEdx_hits.at(1) << ", " <<dEdx_hits.at(2) << ", " <<dEdx_hits.at(3) << std::endl;
      std::cout << "dEdx (average) = " << dEdx << std::endl;
      std::cout << std::endl;

      // Simple functions for "average" mode (no truncation)
      
      float min_purity = 0.2; // minimum required purity for pid hypothesis
      
      std::cout << "- Total probability method -" << std::endl;

      auto probT = getter->GetTrdProbabilitiesTotalAverage(mom, dEdx, charge, nhits_trd);

      for (const auto& pdg : pid_codes_trd_) 
	printf("mom = %.4f dEdx (average) = %.4f  probability %-4s %-10d = %.4f\n", mom,  dEdx, pdg.second.Data(), pdg.first, probT[pdg.first]);
      std::cout << std::endl;

      int pid_hypothesis = getter->GetTrdPid(probT, min_purity, charge);

      std::cout << "pid hypothesis : " << pid_hypothesis << std::endl;
      std::cout << std::endl;
      
      std::cout << "- Likelihood method -" << std::endl;

      auto probL = getter->GetTrdProbabilitiesLikelihoodAverage(mom, dEdx_hits, charge);
      
      for (const auto& pdg : pid_codes_trd_) 
	printf("mom = %.4f dEdx (hits) = %.4f, %.4f, %.4f, %.4f  probability %-4s %-10d = %.4f\n", mom,  dEdx_hits.at(0),  dEdx_hits.at(1),  dEdx_hits.at(2),  dEdx_hits.at(3), pdg.second.Data(), pdg.first, probL[pdg.first]);
      std::cout << std::endl;
      
      pid_hypothesis = getter->GetTrdPid(probL, min_purity, charge);  

      std::cout << "pid hypothesis : " << pid_hypothesis << std::endl;
      std::cout << std::endl;
      
      // Alternative, more flexible way to run getter:
      // Functions with options for truncation mode
      /*
      float min_purity = 0.2; // minimum required purity for pid hypothesis

      //Options:
      //getter->SetTruncationMode(0); // Default truncation mode is "average", optional: set truncation mode (see README or GetterTrd.h)
      //getter->SetMinHits(1);        // Default minimum number of required hits per track is 1

      TrdContainer trdtrack(mom, pT, charge, nhits_trd, dEdx_hits);
	     
      std::cout << "- Total probability method -" << std::endl;

      getter->SetProbabiltyMode(0);
      auto probT = getter->GetTrdProbabilities(trdtrack); 

      for (const auto& pdg : pid_codes_trd_) 
	printf("mom = %.4f dEdx (average) = %.4f  probability %-4s %-10d = %.4f\n", mom,  dEdx, pdg.second.Data(), pdg.first, probT[pdg.first]);
      std::cout << std::endl;

      int pid_hypothesis = getter->GetTrdPid(probT, min_purity, charge);

      std::cout << "pid hypothesis : " << pid_hypothesis << std::endl;
      std::cout << std::endl;
      
      std::cout << "- Likelihood method -" << std::endl;

      getter->SetProbabiltyMode(1);
      auto probL = getter->GetTrdProbabilitiesMulti(trdtrack); 
      
      for (const auto& pdg : pid_codes_trd_) 
	printf("mom = %.4f dEdx (hits) = %.4f, %.4f, %.4f, %.4f  probability %-4s %-10d = %.4f\n", mom,  dEdx_hits.at(0),  dEdx_hits.at(1),  dEdx_hits.at(2),  dEdx_hits.at(3), pdg.second.Data(), pdg.first, probL[pdg.first]);
	std::cout << std::endl;

	pid_hypothesis = getter->GetTrdPid(probL, min_purity, charge);

	std::cout << "pid hypothesis : " << pid_hypothesis << std::endl;
	std::cout << std::endl;
      */

    }
}
