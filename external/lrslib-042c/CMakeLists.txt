add_definitions(-DTIMES -DGMP)
add_library(lrsgmp SHARED lrslib.c lrsgmp.c)

target_link_libraries(lrsgmp ${GMP_LIBRARIES})

install (TARGETS lrsgmp DESTINATION lib)
