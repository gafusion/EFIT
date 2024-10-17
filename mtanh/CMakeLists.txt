######################################################################
#
# CMakeLists.txt for mtanh routines
#
######################################################################

set(MTANH_SOURCES
 mds_rg.f90 mds_mtanh.f90 mds_mtanh_rdts.f90 mdslib.f90
)
add_library(mtanh ${MTANH_SOURCES})

# Need to link MDS+ (found in io)
target_link_libraries(mtanh ${MDSPLUS_LIBRARIES})

install(TARGETS mtanh DESTINATION lib
  PERMISSIONS OWNER_READ OWNER_EXECUTE OWNER_WRITE
              GROUP_READ GROUP_EXECUTE GROUP_WRITE
              WORLD_READ WORLD_EXECUTE
)

# Build test executable
add_executable(test_mtanh test_get_mtanh_ts.f90)
# TODO: linking could be done more cleanly
target_link_libraries(test_mtanh mtanh ${MDSPLUS_LIBRARIES})
