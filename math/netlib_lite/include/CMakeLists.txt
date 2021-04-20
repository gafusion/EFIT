######################################################################
#
# CMakeLists.txt for vpstd
#
# $Id: CMakeLists.txt 53 2011-09-29 18:28:47Z cary $
#
######################################################################

# include_directories(
#   ${Hdf5_INCLUDE_DIRS}/fortran
# )

set(LSODE_SOURCES
  dcfode.f dewset.f dgbfa.f dgbsl.f dgefa.f dgesl.f dintdy.f dprepj.f
  dsolsy.f dsrcom.f dstode.f dumach.f dvnorm.f iumach.f ixsav.f lsode.f
  lsode_module.f xerrwd.f
)

add_library(lsode ${LSODE_SOURCES})

install(FILES lsode_dummy.h
  DESTINATION include
)

# Do for convenience lib only
if (0)
install(TARGETS lsode
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
endif ()

