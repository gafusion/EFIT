######################################################################
#
# CMakeLists.txt for efit test
#
######################################################################

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/mhdin.dat 
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

configure_file(${EFUND_TEST_CONFIGFILE} 
               ${CMAKE_CURRENT_BINARY_DIR}/run_green.sh @ONLY)
configure_file(${EFUND_RESULTS_CONFIGFILE}
               ${CMAKE_CURRENT_BINARY_DIR}/check_green_results.sh @ONLY)

add_test(NAME ITER-green COMMAND run_green.sh)
add_test(NAME ITER-green-results COMMAND check_green_results.sh)
