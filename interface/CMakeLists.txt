set(SOURCES
	PidTrdRunGetter.cpp
  )

string(REPLACE ".cpp" ".hpp" HEADERS "${SOURCES}")
include_directories(${PROJECT_INCLUDE_DIRECTORIES} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src_trd)
add_library(PidInt SHARED ${SOURCES} G__PidInt.cxx)

message("${PROJECT_INCLUDE_DIRECTORIES}")

ROOT_GENERATE_DICTIONARY(G__PidInt ${HEADERS} LINKDEF PidIntLinkDef.h)
target_link_libraries(PidInt ${ROOT_LIBRARIES} PidTrd)

install(
       FILES
       ${HEADERS}
       DESTINATION
       include
       COMPONENT
       Devel
)

set(PCM_FILE_NAME libPidInt)

install(
       FILES
       "${CMAKE_CURRENT_BINARY_DIR}/${PCM_FILE_NAME}_rdict.pcm"
       "${CMAKE_CURRENT_BINARY_DIR}/${PCM_FILE_NAME}.rootmap"
       DESTINATION
       lib
       OPTIONAL
)
