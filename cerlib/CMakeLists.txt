######################################################################
#
# CMakeLists.txt for cerlib
#
######################################################################

set(CERLIB_SOURCES
 mds_rg.f90 mds_mtanh.f90 mds_mtanh_rdts.f90 mdslib.f90
)

add_library(cerlib ${CERLIB_SOURCES})
install(TARGETS cerlib DESTINATION lib
  PERMISSIONS OWNER_READ OWNER_EXECUTE OWNER_WRITE
              GROUP_READ GROUP_EXECUTE GROUP_WRITE
              WORLD_READ WORLD_EXECUTE
)

# Need to link MDS+ (found in io)
if(${HAVE_MDSPLUS})
  #include_directories(${MDSPLUS_INCLUDE_DIRS})
  target_link_libraries(cerlib ${MDSPLUS_LIBRARIES})
endif()

# It would be nice to have dedicated tests
#add_subdirectory(tests)
