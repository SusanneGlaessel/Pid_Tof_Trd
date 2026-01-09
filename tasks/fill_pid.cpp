#include "PidFiller.hpp"

#include "AnalysisTree/TaskManager.hpp"

using namespace AnalysisTree;

void fill_pid(const std::string& filelist, int pid_mode, const std::string& pid_file_tof, const std::string& pid_file_trd, const std::string& output, int truncation_mode, int probability_mode, int trdhits_min, float purity) {
    
  if (pid_mode == 0) {
    std::cout << "\n -I- fill_pid : Run with Tof and Trd pid.\n" << std::endl;
    std::cout << " -I- fill_pid : Using input files " << pid_file_tof << " " << pid_file_trd<<std::endl;
    std::cout << std::endl;
  }
  if (pid_mode == 1) {
    std::cout << "\n -I- fill_pid : Run with Tof pid only.\n" << std::endl;
    std::cout << " -I- fill_pid : Using input file " << pid_file_tof <<std::endl;
    std::cout << std::endl;
  }
  if (pid_mode == 2) {
    std::cout << "\n -I- fill_pid : Run with Trd pid only.\n" << std::endl;
    std::cout << " -I- fill_pid : Using input file " << pid_file_trd <<std::endl;
    std::cout << std::endl;
  }

  if (pid_mode == 0 || pid_mode == 2) {
    std::cout << "\n************************" << std::endl;
    std::cout << "Settings for Trd pid: "   << std::endl;
    std::cout << "  truncation mode  : " << truncation_mode  << std::endl;
    std::cout << "  probability mode : " << probability_mode << std::endl;
    std::cout << "  min hits trd     : " << trdhits_min      << std::endl;
    if (purity != 0.0) std::cout << "  min purtiy       : " << purity           << std::endl; 
    std::cout << "************************\n" << std::endl;
  }

  
  auto* man = TaskManager::GetInstance();
  man->SetOutputName(output, "aTree");

  // The output file will be based on the input one with additional branches (and possibly with removed other branches)
  man->SetWriteMode(eBranchWriteMode::kCopyTree);

  // The branch VtxTracks, which is an input branch, will be also removed, because a new one, based on it, will be created.
  man->SetBranchesExclude({"VtxTracks"});

  // Initialize the Pid::Getter prepared at the fitting step: 1-st argument - file, 2-nd argument - name of the getter.
  auto* pid_task = new PidFiller(pid_file_tof, pid_file_trd, "pid_getter_tof", "pid_getter_trd", pid_mode);

  pid_task -> SetRecTracksName("VtxTracks");
  pid_task -> SetTofHitsName("TofHits");
  pid_task -> SetTrdTracksName("TrdTracks");
  pid_task -> SetRichRingsName("RichRings");

  if ( pid_mode == 0 || pid_mode == 2) {
    pid_task -> SetMinHits(trdhits_min);
    pid_task -> SetTruncationMode(truncation_mode);
    pid_task -> SetProbabilityMode(probability_mode);
    pid_task -> SetPurity(purity);
  }
  
  man->AddTask(pid_task);

  man->Init({filelist}, {"rTree"});
  man->Run(-1);
  man->Finish();
}

int main(int argc, char** argv) {
 
  if (argc < 5) {
    std::cout << std::endl;
    std::cout << "Wrong number of arguments! Please use: " << std::endl;
    std::cout << std::endl;
    
    std::cout << "For Tof and Trd pid: " << std::endl;
    std::cout << " ./fill_pid 0 filelist.txt outputfile pid_file_tof pid_file_trd truncation_mode probability_mode min_hits\n";
    std::cout << "or: " << std::endl;
    std::cout << " ./fill_pid 0 filelist.txt outputfile pid_file_tof pid_file_trd truncation_mode probability_mode min_hits purity\n";
    std::cout << std::endl;
    
    std::cout << "For Tof pid only: " << std::endl;
    std::cout << " ./fill_pid 1 filelist.txt outputfile pid_file_tof\n";
    std::cout << std::endl;
    
    std::cout << "For Trd pid only: " << std::endl;
    std::cout << " ./fill_pid 2 filelist.txt outputfile pid_file_trd truncation_mode probability_mode min_hits\n";
    std::cout << "or: " << std::endl;
    std::cout << " ./fill_pid 2 filelist.txt outputfile pid_file_trd truncation_mode probability_mode min_hits purity\n";
    std::cout << std::endl;
    
    return EXIT_FAILURE;
  }
   
  int pid_mode = atoi(argv[1]);

  if (pid_mode == 0 && argc < 9) {
    std::cout << std::endl;
    std::cout << "Wrong number of arguments for running Tof and Trd pid simultanously! Please use: " << std::endl;
    std::cout << " ./fill_pid 0 filelist.txt outputfile pid_file_tof pid_file_trd truncation_mode probability_mode min_hits\n";
    std::cout << "or: " << std::endl;
    std::cout << " ./fill_pid 0 filelist.txt outputfile pid_file_tof pid_file_trd truncation_mode probability_mode min_hits purity\n";
    std::cout << std::endl;
    return EXIT_FAILURE;
  }

  if (pid_mode == 2 && argc < 8) {
    std::cout << std::endl;
    std::cout << "Wrong number of arguments for running Trd pid only! Please use: " << std::endl;
    std::cout << " ./fill_pid 2 filelist.txt outputfile pid_file_trd truncation_mode probability_mode min_hits\n";
    std::cout << "or: " << std::endl;
    std::cout << " ./fill_pid 2 filelist.txt outputfile pid_file_trd truncation_mode probability_mode min_hits purity\n";
    std::cout << std::endl;
    return EXIT_FAILURE;
  }

  const std::string filelist    = argv[2];
  const std::string output_file = argv[3];

  std::string pid_file_tof = ""; std::string pid_file_trd = "";
  if (pid_mode == 0) {
    pid_file_tof = argv[4];
    pid_file_trd = argv[5];
  }
  if (pid_mode == 1) {
    pid_file_tof = argv[4];
  }
  if (pid_mode == 2) {
    pid_file_trd = argv[4];
  }

  int truncation_mode = -1; int probability_mode = -1;  int trdhits_min = -1; double purity = 0.0;
  if (pid_mode == 0 || pid_mode == 2) {
    int i = 0;
    if (pid_mode == 0) i = 1;
    truncation_mode = atoi(argv[5+i]);
    probability_mode = atoi(argv[6+i]);
    trdhits_min = atoi(argv[7+i]);
    if (argv[8+i]) purity = atof(argv[8+i]);
    else purity = 0.0;
  }
  fill_pid(filelist, pid_mode, pid_file_tof, pid_file_trd, output_file, truncation_mode, probability_mode, trdhits_min, purity);
  
  return EXIT_SUCCESS;
}
