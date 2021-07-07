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
add_test(NAME efit02er COMMAND run_efit.sh)
