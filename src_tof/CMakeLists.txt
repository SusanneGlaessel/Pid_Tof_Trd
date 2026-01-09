set(SOURCES
        Fitter.cpp
        ParticleFit.cpp
        Parameters.cpp
        GetterTof.cpp
        )

string(REPLACE ".cpp" ".h" HEADERS "${SOURCES}")
list(APPEND HEADERS ConstantsTof.h)

set(DICT_FILE_NAME G__${PROJECT_NAME})
set(PCM_FILE_NAME lib${PROJECT_NAME})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
ROOT_GENERATE_DICTIONARY(${DICT_FILE_NAME} ${HEADERS} LINKDEF PidLinkDef.h)
add_library(Pid SHARED ${SOURCES} ${HEADERS} ${DICT_FILE_NAME})
foreach (_HEADER ${HEADERS})
    get_filename_component(_HEADER_FILENAME ${_HEADER} NAME)
    get_filename_component(_HEADER_PATH ${_HEADER} REALPATH)
    set(_INCLUDE_BINARY_DIR ${CMAKE_BINARY_DIR}/include/pid)
    set(_HEADER_SYMLINK_PATH ${_INCLUDE_BINARY_DIR}/${_HEADER_FILENAME})
    add_custom_command(TARGET Pid PRE_BUILD
            COMMAND ${CMAKE_COMMAND} ARGS -E make_directory ${_INCLUDE_BINARY_DIR}
            COMMAND ${CMAKE_COMMAND} ARGS -E create_symlink ${_HEADER_PATH} ${_HEADER_SYMLINK_PATH}
            )
endforeach ()

target_link_libraries(Pid ${ROOT_LIBRARIES})
target_include_directories(Pid
        PRIVATE ${CMAKE_SOURCE_DIR}/src_tof
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include/pid>
        $<INSTALL_INTERFACE:include/pid>
        )

install(
        FILES
        "${CMAKE_CURRENT_BINARY_DIR}/${PCM_FILE_NAME}_rdict.pcm"
        "${CMAKE_CURRENT_BINARY_DIR}/${PCM_FILE_NAME}.rootmap"
        DESTINATION
        lib
        OPTIONAL
)

install(FILES ${HEADERS} DESTINATION include)
