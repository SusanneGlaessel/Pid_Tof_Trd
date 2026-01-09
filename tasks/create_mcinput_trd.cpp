#include <TROOT.h>

#include "PidTrdMcInput.hpp"
#include "AnalysisTree/TaskManager.hpp"

using namespace AnalysisTree;

int create_mcinput_trd(const std::string& filelist, const std::string& outfilename) {

  auto start = std::chrono::system_clock::now();
  ROOT::EnableImplicitMT(4);

  gROOT -> SetBatch (true);
  gROOT -> SetEscape (true);
  
  Bool_t update_mchistos = kFALSE;
  
  auto* pid_trd_mcinput = new PidTrdMcInput(outfilename);
  //pid_trd_mcinput->SetRecTracksName("RecTracks"); //analysistree with tof-pid
  pid_trd_mcinput->SetRecTracksName("VtxTracks");   // analysistree without tof-pid
  pid_trd_mcinput->SetSimTracksName("SimParticles");
  pid_trd_mcinput->SetTrdTracksName("TrdTracks");
  pid_trd_mcinput->SetUpdateMcHistos(update_mchistos);

  auto* man = TaskManager::GetInstance();
  man->AddTask(pid_trd_mcinput);  
  //man->Init({filelist}, {"pTree"}); //analysistree with tof-pid
  man->Init({filelist}, {"rTree"});  // analysistree without tof-pid
  man->Run(-1);// -1 = all events
  man->Finish();
  man->ClearTasks();
  
  return EXIT_SUCCESS;
}

int main(int argc, char **argv) {
  
  if (argc < 3) {
    std::cout << "Wrong number of arguments! Please use: " << std::endl;
    std::cout << " ./create_mcinput filelist.txt outfilename\n";
    return EXIT_FAILURE;
  }

  const std::string& filelist    = argv[1];
  const std::string& outfilename = argv[2];
  create_mcinput_trd(filelist, outfilename);
  
  return 0;
}
