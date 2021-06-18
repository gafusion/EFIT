######################################################################
#
# CMakeLists.txt for nimtests
#
# Complete harness is made of:
#    nimTests.cmake, nimtest_template.sh.in, and nimharness.sh
#
######################################################################


set(EFIT_LABEL 186610.02300)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/k${EFIT_LABEL} 
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

configure_file(${EFIT_TEST_CONFIGFILE} 
               ${CMAKE_CURRENT_BINARY_DIR}/run_efit.sh @ONLY)
add_test(NAME efit01 COMMAND run_efit.sh)

