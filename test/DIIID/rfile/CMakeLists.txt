######################################################################
#
# CMakeLists.txt for efit tests
#
######################################################################

set(EFIT_LABEL 186610.02300)
set(EFIT_IN rfile.txt)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${EFIT_IN}
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/gin
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

# This test can only be run serially
# Must be run on execute node at NERSC to use MPI
if (${ENABLE_PARALLEL})
  set(EFIT_EXEC "srun --mpi=pmi2 -n 1 ${CMAKE_BINARY_DIR}/efit/efit")
endif()

configure_file(${EFIT_TEST_CONFIGFILE}
               ${CMAKE_CURRENT_BINARY_DIR}/run_efit.sh @ONLY)
configure_file(${EFIT_RESULTS_CONFIGFILE}
               ${CMAKE_CURRENT_BINARY_DIR}/check_efit_results.sh @ONLY)

add_test(NAME rfile COMMAND run_efit.sh)
add_test(NAME rfile-results COMMAND check_efit_results.sh)
