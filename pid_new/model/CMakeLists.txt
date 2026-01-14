include_directories(${CMAKE_CURRENT_SOURCE_DIR})
ROOT_GENERATE_DICTIONARY(G__PidModel
        AsymmetricGaussianPDF.h
        LINKDEF PidModelLinkDef.h
        OPTIONS -inlineInputHeader
        )

add_library(pid_new_model SHARED AsymmetricGaussianPDF.cpp G__PidModel.cxx)
target_include_directories(pid_new_model PUBLIC
        $<BUILD_INTERFACE:${Pid_SOURCE_DIR}>)
target_link_libraries(pid_new_model PUBLIC ROOT::RooFit)

