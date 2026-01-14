

add_library(pid_new_core SHARED BaseParticleFitModel.cpp FitterHelper.cpp PdgHelper.cpp PdgHelper.h FitParameter.cpp FitParameter.h ParticleFitModel.cpp ParticleFitModel.h)
target_link_libraries(pid_new_core PUBLIC ROOT::RooFit ROOT::EG)
target_include_directories(pid_new_core PUBLIC
        $<BUILD_INTERFACE:${Pid_SOURCE_DIR}>)
set_target_properties(pid_new_core PROPERTIES CXX_STANDARD 11)