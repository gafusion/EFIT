######################################################################
#
# CMakeLists.txt for other
#
# $Id: CMakeLists.txt 53 2011-09-29 18:28:47Z cary $
#
######################################################################

set(NLOTHER_SOURCES
  r8tomsqz.F  tomsqz.f
)

add_library(nlother ${NLOTHER_SOURCES})

install(TARGETS nlother
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

