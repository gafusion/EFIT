######################################################################
#
# CMakeLists.txt for efit tests
#
######################################################################

set(EFIT_LABEL 186610.00000)
set(EFIT_NEXT 186610.00001)
set(EFIT_IN rfile.txt)
set(RFILE_TEST True)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${EFIT_IN}
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/gin
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

# This test can only be run serially
# Must still be run on execute node at NERSC to use MPI
if (${ENABLE_PARALLEL})
  set(NPROC "1")
endif ()

configure_file(${EFIT_TEST_CONFIGFILE}
               ${CMAKE_CURRENT_BINARY_DIR}/run_efit.sh @ONLY)
configure_file(${EFIT_RESULTS_CONFIGFILE}
               ${CMAKE_CURRENT_BINARY_DIR}/check_efit_results.sh @ONLY)

add_test(NAME DIIID-rfile COMMAND run_efit.sh)
add_test(NAME DIIID-rfile-results COMMAND check_efit_results.sh)

if (${ENABLE_HDF5})
  set(RUNID 186610)
  set(EFIT_MODE 1)
  set(EFIT_IN ${RUNID}.h5)
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${EFIT_IN}
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
  configure_file(${EFIT_TEST_CONFIGFILE_HDF5}
                 ${CMAKE_CURRENT_BINARY_DIR}/run_efit_hdf5.sh @ONLY)
  add_test(NAME DIIID-rfile-hdf5 COMMAND run_efit_hdf5.sh)
  add_test(NAME DIIID-rfile-hdf5-results COMMAND check_efit_results.sh)
endif ()
