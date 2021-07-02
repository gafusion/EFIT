######################################################################
#
# CMakeLists.txt for lsode
#
# $Id: CMakeLists.txt 58 2011-10-04 21:15:02Z cary $
#
######################################################################

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)
set(LSODE_SOURCES
  isize_get_lsode.c xerrwv.f idamax.f
  dcfode.f dewset.f dgbfa.f dgbsl.f dgefa.f dgesl.f dintdy.f dprepj.f
  dsolsy.f dsrcom.f dstode.f dumach.f dvnorm.f iumach.f ixsav.f lsode.f
  lsode_module.f xerrwd.f
  lsode_r8.F
)

add_library(lsode ${LSODE_SOURCES})

install(FILES lsode_dummy.h ${CMAKE_CURRENT_BINARY_DIR}/lsode_mod.mod
  DESTINATION include
)

install(TARGETS lsode
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

