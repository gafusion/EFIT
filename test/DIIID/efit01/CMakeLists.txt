######################################################################
#
# CMakeLists.txt for efit tests
#
######################################################################

set(EFIT_LABEL 186610.02300)
set(EFIT_IN k${EFIT_LABEL})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${EFIT_IN}
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

configure_file(${EFIT_TEST_CONFIGFILE}
               ${CMAKE_CURRENT_BINARY_DIR}/run_efit.sh @ONLY)
configure_file(${EFIT_RESULTS_CONFIGFILE}
               ${CMAKE_CURRENT_BINARY_DIR}/check_efit_results.sh @ONLY)

add_test(NAME efit01 COMMAND run_efit.sh)
add_test(NAME efit01-results COMMAND check_efit_results.sh)
