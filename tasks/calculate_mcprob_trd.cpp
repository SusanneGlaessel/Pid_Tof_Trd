#include <TROOT.h>
#include <iostream>
#include "PidTrdRunGetter.hpp"
#include "TStopwatch.h"

int calculate_mcprob_trd(const std::string& mcfile_name) {
  
  Bool_t write_mchistos_out = kTRUE;

  TString macroname = "create_getter";
 
  const std::string getter_file = "pid_getter_trd.root";
  const std::string getter_name = "pid_getter_trd";
    
  TString outFile;
  if (write_mchistos_out == kTRUE) {
      TString name = mcfile_name;
      name.Replace(name.Index(".root"),5,"");
      outFile = Form("%s_probabilities.root", name.Data());
  }
  
  std::cout << std::endl;
  std::cout << "-I- " << macroname << ": Using input file " << mcfile_name << std::endl;

  TStopwatch timer;
  timer.Start();
  
  auto start = std::chrono::system_clock::now(); 
  ROOT::EnableImplicitMT(4);

  gROOT -> SetBatch (true);
  gROOT -> SetEscape (true);

  auto* pid_trd_rungetter = new PidTrdRunGetter(mcfile_name, getter_file, getter_name);
  pid_trd_rungetter->SetWriteMcHistogramsOut(write_mchistos_out);
  pid_trd_rungetter->Init();
  pid_trd_rungetter->Exec();
  pid_trd_rungetter->Finish();

  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  std::cout << std::endl << std::endl;
  std::cout << "Macro finished successfully." << std::endl;
  std::cout << "Getter is " << getter_file << std::endl;
  if (write_mchistos_out == kTRUE) std::cout << "Output file is " << outFile << std::endl;
  std::cout << "Real time " << rtime << " s, CPU time " << ctime << "s"
            << std::endl
            << std::endl;

  return EXIT_SUCCESS;
}

int main(int argc, char **argv) {

  
  if (argc < 2) {
    std::cout << "Wrong number of arguments! Please use: " << std::endl;
    std::cout << " ./create_getter filename_mc\n";
    return EXIT_FAILURE;
  }

  const std::string& mcfile_name = argv[1];
  calculate_mcprob_trd(mcfile_name);
  
  return 0;
}
