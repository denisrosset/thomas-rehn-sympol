add_definitions(-DTIMES -DGMP)
add_library(lrsgmp SHARED lrslib.c lrsgmp.c)

install (TARGETS lrsgmp DESTINATION lib)
