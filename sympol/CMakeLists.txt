INCLUDE (CheckSymbolExists) 
CHECK_SYMBOL_EXISTS(sysconf "stdlib.h;unistd.h" HAVE_SYSCONF_PROTO)
CHECK_SYMBOL_EXISTS(getrusage "sys/resource.h" HAVE_GETRUSAGE_PROTO) 

INCLUDE (CheckIncludeFile)
IF(NAUTY_NTL_CHECK)
  SET(CMAKE_REQUIRED_INCLUDES ${PROJECT_SOURCE_DIR}/include)
  CHECK_INCLUDE_FILE("nauty.h" HAVE_NAUTY)
  CHECK_INCLUDE_FILE("NTL/config.h" HAVE_NTL)
ENDIF(NAUTY_NTL_CHECK)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)



SET (Boost_FIND_REQUIRED TRUE)
SET (Boost_FIND_QUIETLY TRUE)
SET (Boost_USE_MULTITHREADED TRUE)
find_package(Boost 1.34.1 COMPONENTS program_options REQUIRED)

include_directories(../external/boost ../external/permlib ../external/lrslib-042c ../external/cddlib-094f ${CMAKE_CURRENT_BINARY_DIR})
add_definitions(-Woverloaded-virtual -ansi -pedantic)
# lrs flags
add_definitions(-DGMP)
# cdd flags
add_definitions(-DGMPRATIONAL)

IF(HAVE_NAUTY)
  # nauty flags
  add_definitions(-DINT_SIZE=32 -DWORDSIZE=32)
ENDIF(HAVE_NAUTY)

IF(HAVE_NAUTY)
    add_library(libnauty STATIC IMPORTED)
    set_property(TARGET libnauty PROPERTY
                 IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/external/nauty/libnauty.a)
ENDIF(HAVE_NAUTY)

IF(HAVE_NTL)
    add_library(libntl STATIC IMPORTED)
    set_property(TARGET libntl PROPERTY
                 IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/external/NTL/libntl.a)
ENDIF(HAVE_NTL)

add_executable (sympol 
  configuration.cpp
  polyhedronio.cpp 
  sympol.cpp 
  polyhedron.cpp
  polyhedrondatastorage.cpp
  qarray.cpp
  automorphismcomputation.cpp
  raycomputationlrs.cpp
  raycomputationcdd.cpp
  recursionstrategy.cpp
  symmetrycomputation.cpp
  symmetrycomputationdirect.cpp
  symmetrycomputationadm.cpp
  symmetrycomputationidm.cpp
  facesuptosymmetrylist.cpp
  yal/reportlevel.cpp
  yal/usagestats.cpp
)

set(TARGET_LIBS lrsgmp cddgmp gmpxx ${Boost_LIBRARIES})
IF(HAVE_NAUTY)
    MESSAGE("building with nauty")
    set(TARGET_LIBS ${TARGET_LIBS} libnauty)
ENDIF(HAVE_NAUTY)
IF(HAVE_NTL)
    MESSAGE("building with NTL")
    set(TARGET_LIBS ${TARGET_LIBS} libntl)
ENDIF(HAVE_NTL)

target_link_libraries(sympol ${TARGET_LIBS})

install (TARGETS sympol DESTINATION bin)

