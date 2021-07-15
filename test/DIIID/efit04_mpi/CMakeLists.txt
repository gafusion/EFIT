######################################################################
#
# CMakeLists.txt for efit tests
#
######################################################################

set(RUNID 186610)
set(TIMES 02300 02400 02500 02600)

foreach( t ${TIMES})
  set(EFIT_IN k${RUNID}.${t})
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${EFIT_IN}
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
endforeach()
string(REPLACE ";" " " TIMES "${TIMES}")

configure_file(${EFIT_MPI_TEST_CONFIGFILE}
               ${CMAKE_CURRENT_BINARY_DIR}/run_efit_mpi.sh @ONLY)
configure_file(${EFIT_RESULTS_CONFIGFILE}
               ${CMAKE_CURRENT_BINARY_DIR}/check_efit_results.sh @ONLY)

add_test(NAME efit04_mpi COMMAND run_efit_mpi.sh)
add_test(NAME efit04_mpi-results COMMAND check_efit_results.sh)
