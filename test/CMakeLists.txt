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
set(EFUND_TEST_CONFIGFILE ${CMAKE_CURRENT_SOURCE_DIR}/run_green.sh.in)
set(EFUND_RESULTS_CONFIGFILE ${CMAKE_CURRENT_SOURCE_DIR}/check_green_results.sh.in)

# This is a shell script used to drive running efit
if (${ENABLE_PARALLEL})
   set(EFIT_EXEC "srun --mpi=pmi2 -n 4 ${CMAKE_BINARY_DIR}/efit/efit")
else ()
   set(EFIT_EXEC ${CMAKE_BINARY_DIR}/efit/efit)
endif ()
set(EFIT_MODE 2)
set(EFIT_TIME 1)
set(EFIT_TEST_CONFIGFILE ${CMAKE_CURRENT_SOURCE_DIR}/run_efit.sh.in)
set(EFIT_RESULTS_CONFIGFILE ${CMAKE_CURRENT_SOURCE_DIR}/check_efit_results.sh.in)

add_subdirectory(DIIID)
