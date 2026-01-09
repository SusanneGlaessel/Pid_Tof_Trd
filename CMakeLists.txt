cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(Pid)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_cpp_STANDARD 11)
set(PROJECT_VERSION 1.0)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules)
#
# Options
#
set(PID_BUNDLED_AT OFF CACHE BOOL "Get and build AnalysisTree")
set(PID_BUNDLED_AT_VERSION "v2.2.7" CACHE STRING "Bundled AnalysisTree version")
set(PID_BUNDLED_AT_GIT_SHALLOW ON CACHE BOOL "Use CMake GIT_SHALLOW option")

# Telling CMake where to find the ROOT installation.
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/lib")

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS Core MathCore RIO Hist Tree)

# Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})
add_definitions(${ROOT_CXX_FLAGS})

if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    option(Pid_BUILD_NEW "Build experimental PID" OFF)
endif ()

set(Pid_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# by default build optimized code with debug symbols
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE RELEASE)
endif ()

# in DEBUG mode make verbose Makefile
if (CMAKE_BUILD_TYPE MATCHES DEBUG)
    set(CMAKE_VERBOSE_MAKEFILE ON)
endif ()

if (Pid_BUILD_NEW)
    add_subdirectory(pid_new)
endif ()

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g -DDEBUG -D__DEBUG -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ftree-vectorize -ffast-math -DNODEBUG")
message(STATUS "Using CXX flags for ${CMAKE_BUILD_TYPE}: ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")

add_subdirectory(src_tof)
add_subdirectory(src_trd)
add_subdirectory(interface)
if(${PID_BUNDLED_AT})
    add_subdirectory(at_interface)
endif()
add_subdirectory(tasks)

# Packaging routine
# For complete explanation
# @see https://cmake.org/cmake/help/git-master/manual/cmake-packages.7.html#creating-packages
include(GenerateExportHeader)

generate_export_header(Pid)
set_property(TARGET Pid PROPERTY VERSION ${PROJECT_VERSION})
set_property(TARGET Pid PROPERTY SOVERSION ${PROJECT_VERSION})
set_property(TARGET Pid PROPERTY INTERFACE_Pid_MAJOR_VERSION 1)
set_property(TARGET Pid APPEND PROPERTY COMPATIBLE_INTERFACE_STRING INTERFACE_Pid_MAJOR_VERSION)

install(TARGETS Pid EXPORT PidTargets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include
        )

install(
        FILES
        ${HEADERS}
        "${CMAKE_CURRENT_BINARY_DIR}/pid_export.h"
        DESTINATION
        include/pid
        COMPONENT
        Devel
)

install(
        FILES
        "${CMAKE_CURRENT_BINARY_DIR}/${PCM_FILE_NAME}_rdict.pcm"
        "${CMAKE_CURRENT_BINARY_DIR}/${PCM_FILE_NAME}.rootmap"
        DESTINATION
        lib
        OPTIONAL
)

configure_file(cmake_modules/PidConfig.sh.in ${CMAKE_BINARY_DIR}/PidConfig.sh)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/Pid/PidConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion
)

export(EXPORT PidTargets
        FILE "${CMAKE_CURRENT_BINARY_DIR}/Pid/PidTargets.cmake")

set(ConfigPackageLocation lib/cmake/Pid)

set(Pid_INCLUDE_DIR "include")
set(Pid_LIBRARY_DIR "lib")

configure_package_config_file(
        PidConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/Pid/PidConfig.cmake
        INSTALL_DESTINATION ${ConfigPackageLocation}
        PATH_VARS Pid_INCLUDE_DIR Pid_LIBRARY_DIR
        #  [NO_SET_AND_CHECK_MACRO]
        #  [NO_CHECK_REQUIRED_COMPONENTS_MACRO]
        #  [INSTALL_PREFIX <path>]
)

install(EXPORT PidTargets
        FILE
        PidTargets.cmake
        DESTINATION
        ${ConfigPackageLocation}
        )

install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/Pid/PidConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/Pid/PidConfigVersion.cmake"
        DESTINATION
        ${ConfigPackageLocation}
        COMPONENT
        Devel
        )

install(FILES
        ${CMAKE_BINARY_DIR}/PidConfig.sh
        DESTINATION
        bin
        COMPONENT
        Devel
        )
