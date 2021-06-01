######################################################################
#
# CMakeLists.txt for nimtests
#
# Complete harness is made of:
#    nimTests.cmake, nimtest_template.sh.in, and nimharness.sh
#
######################################################################

#---------------------------------------------------------------------
# Initialization
#---------------------------------------------------------------------
include_directories(
  ${CMAKE_BINARY_DIR}/efit
)

# This is a shell script used to drive running efit
set(EFIT_EXEC ${CMAKE_BINARY_DIR}/efit/efit)
set(EFIT_SIZE 65)
set(EFIT_MODE 2)
set(EFIT_TIME 1)
set(EFIT_TEST_CONFIGFILE ${CMAKE_CURRENT_SOURCE_DIR}/run_efit.sh.in)

add_subdirectory(DIIID)
