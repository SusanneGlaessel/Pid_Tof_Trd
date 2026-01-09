set(SOURCES
        ParticleProb.cpp	
	ContainerTrd.cpp
        GetterTrd.cpp
        )

string(REPLACE ".cpp" ".h" HEADERS "${SOURCES}")
list(APPEND HEADERS ConstantsTrd.h)
include_directories(${PROJECT_INCLUDE_DIRECTORIES} ${CMAKE_CURRENT_SOURCE_DIR})
add_library(PidTrd SHARED ${SOURCES} G__PidTrd.cxx)

message("${PROJECT_INCLUDE_DIRECTORIES}")

ROOT_GENERATE_DICTIONARY(G__PidTrd ${HEADERS} LINKDEF PidTrdLinkDef.h)
target_link_libraries(PidTrd ${ROOT_LIBRARIES})

install(
       FILES
       ${HEADERS}
       DESTINATION
       include
       COMPONENT
       Devel
)

set(PCM_FILE_NAME libPidTrd)

install(
       FILES
       "${CMAKE_CURRENT_BINARY_DIR}/${PCM_FILE_NAME}_rdict.pcm"
       "${CMAKE_CURRENT_BINARY_DIR}/${PCM_FILE_NAME}.rootmap"
       DESTINATION
       lib
       OPTIONAL
)














