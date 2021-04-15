######################################################################          
#                                                                               
# CMakeLists.txt for cerlib
#                                                                               
######################################################################

# Need to find MDS+ ?
#include_directories(
#)

set(CERLIB_SOURCES
 mds_rg.f90 mds_mtanh.f90 mds_mtanh_rdts.f90 mds_write_rg.f90  
)

add_library(cerlib ${CERLIB_SOURCES})
install(TARGETS cerlib DESTINATION lib
  PERMISSIONS OWNER_READ OWNER_EXECUTE OWNER_WRITE
              GROUP_READ GROUP_EXECUTE GROUP_WRITE
              WORLD_READ WORLD_EXECUTE
)

# Need a find MDS+
#target_link_libraries(cerlib libMdsLib_client)

# It would be nice to have tests
#add_subdirectory(tests)
