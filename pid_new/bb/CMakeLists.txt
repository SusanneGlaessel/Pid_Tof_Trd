
set(SOURCES
        BetheBlochHelper.cpp
        BetheBlochShine.cpp
        )
set(HEADERS
        BetheBlochHelper.h
        BetheBlochShine.h
        )


include_directories(${CMAKE_CURRENT_SOURCE_DIR})
ROOT_GENERATE_DICTIONARY(G__BetheBlochShine BetheBlochShine.h
        LINKDEF BetheBlochShineLinkDef.h
        OPTIONS -inlineInputHeader)
add_library(pid_new_bb_shine SHARED BetheBlochShine.cpp G__BetheBlochShine.cxx)


add_library(pid_new_bb SHARED BetheBlochShine.cpp BetheBlochHelper.cpp)
target_include_directories(pid_new_bb PUBLIC
        $<BUILD_INTERFACE:${Pid_SOURCE_DIR}>)
target_link_libraries(pid_new_bb PUBLIC pid_new_core)
