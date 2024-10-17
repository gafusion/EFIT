######################################################################
#
# CMakeLists.txt for mtanh test
#
######################################################################

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run_mtanh.sh.in
               ${CMAKE_CURRENT_BINARY_DIR}/run_mtanh.sh @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/check_results.sh.in
               ${CMAKE_CURRENT_BINARY_DIR}/check_results.sh @ONLY)

add_test(NAME DIIID-mtanh COMMAND run_mtanh.sh)
add_test(NAME DIIID-mtanh-results COMMAND check_results.sh)
