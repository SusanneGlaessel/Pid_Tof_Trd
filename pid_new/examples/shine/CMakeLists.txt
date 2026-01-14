add_executable(RunFitShine RunFitShine.cpp)
target_link_libraries(RunFitShine pid_new_core pid_new_bb pid_new_model)

add_executable(FitBinnedTracks FitBinnedTracks.cpp)
target_link_libraries(FitBinnedTracks pid_new_core pid_new_bb pid_new_model ${ROOT_LIBRARIES} ROOTVecOps)

add_executable(FitResultsExplorer FitResultsExplorer.cpp)
target_link_libraries(FitResultsExplorer pid_new_core pid_new_bb pid_new_model ${ROOT_LIBRARIES} ROOTVecOps)
