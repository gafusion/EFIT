######################################################################
#
# CMakeLists.txt for nimtests
#
######################################################################

#---------------------------------------------------------------------
# Initialization
#---------------------------------------------------------------------
include_directories(
  ${CMAKE_BINARY_DIR}/green
  ${CMAKE_BINARY_DIR}/efit
)

set(GRID_SIZE 65)

# This is a shell script used to drive running efund
set(EFUND_EXEC ${CMAKE_BINARY_DIR}/green/efund)
set(EFUND_TEST_CONFIGFILE ${CMAKE_CURRENT_SOURCE_DIR}/run_efund.sh.in)
set(EFUND_RESULTS_CONFIGFILE ${CMAKE_CURRENT_SOURCE_DIR}/check_efund_results.sh.in)

# This is a shell script used to drive running efit
set(EFIT_EXEC ${CMAKE_BINARY_DIR}/efit/efit)
set(EFIT_MODE 2)
set(EFIT_TIME 1)
set(EFIT_TEST_CONFIGFILE ${CMAKE_CURRENT_SOURCE_DIR}/run_efit.sh.in)
set(EFIT_RESULTS_CONFIGFILE ${CMAKE_CURRENT_SOURCE_DIR}/check_efit_results.sh.in)

add_subdirectory(DIIID)
