cmake_minimum_required (VERSION 2.6)
project (SymPol)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

#TODO: remove this condition, once Subversion_WC_INFO doesn't print an error if it can't find
#      a working copy
IF(SVNCHECK)
	find_package(Subversion)
	if(Subversion_FOUND)
		Subversion_WC_INFO(${CMAKE_CURRENT_SOURCE_DIR} ER)
		set(SUBVERSION_REVISION ${ER_WC_REVISION})
	endif(Subversion_FOUND)
ENDIF(SVNCHECK)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/revision.h.in ${CMAKE_CURRENT_BINARY_DIR}/revision.h)
#set(CPACK_PACKAGE_VERSION_PATCH "${SUBVERSION_REVISION}")

include_directories(${CMAKE_CURRENT_BINARY_DIR})

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

