######################################################################
#
# CMakeLists.txt for efit test
#
######################################################################

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/mhdin.dat 
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

configure_file(${EFUND_TEST_CONFIGFILE} 
               ${CMAKE_CURRENT_BINARY_DIR}/run_efund.sh @ONLY)
configure_file(${EFUND_RESULTS_CONFIGFILE}
               ${CMAKE_CURRENT_BINARY_DIR}/check_efund_results.sh @ONLY)

add_test(NAME efund COMMAND run_efund.sh)
add_test(NAME efund-results COMMAND check_efund_results.sh)
