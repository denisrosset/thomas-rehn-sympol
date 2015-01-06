cmake_minimum_required (VERSION 2.6)
project (SymPol)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/external/cmake/Modules")

include_directories(${CMAKE_CURRENT_BINARY_DIR})

FIND_PACKAGE(GMP REQUIRED)
include_directories(${GMP_INCLUDE_DIR})

IF( NOT CMAKE_BUILD_TYPE )
    SET( CMAKE_BUILD_TYPE Debug CACHE STRING
            "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
             FORCE )
ENDIF( NOT CMAKE_BUILD_TYPE )

add_definitions(-Wall)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O3 -g -DBOOST_DISABLE_ASSERTS")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DBOOST_DISABLE_ASSERTS")

add_subdirectory(sympol)
add_subdirectory(external/lrslib-042c)
add_subdirectory(external/cddlib-094f)

exec_program( ${CMAKE_COMMAND} ARGS -E copy_directory ${PROJECT_SOURCE_DIR}/data ${PROJECT_BINARY_DIR}/data )

