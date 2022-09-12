######################################################################
#
# CMakeLists.txt for efit tests
#
######################################################################

set(EFIT_LABEL 091101.00060)
set(EFIT_NEXT 091101.00061)
set(EFIT_IN rfile.txt)
set(RFILE_TEST True)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${EFIT_IN}
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/gin
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

# This test can only be run serially
# Must be run on execute node at NERSC to use MPI
if (${ENABLE_PARALLEL})
  set(NPROC "1")
endif ()

configure_file(${EFIT_TEST_CONFIGFILE}
               ${CMAKE_CURRENT_BINARY_DIR}/run_efit.sh @ONLY)
configure_file(${EFIT_RESULTS_CONFIGFILE}
               ${CMAKE_CURRENT_BINARY_DIR}/check_efit_results.sh @ONLY)

add_test(NAME ITER-rfile-modify COMMAND run_efit.sh)
add_test(NAME ITER-rfile-modify-results COMMAND check_efit_results.sh)

if (${ENABLE_HDF5})
  set(RUNID 091101)
  set(EFIT_MODE 1)
  set(EFIT_IN ${RUNID}.h5)
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${EFIT_IN}
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
  configure_file(${EFIT_TEST_CONFIGFILE_HDF5}
                 ${CMAKE_CURRENT_BINARY_DIR}/run_efit_hdf5.sh @ONLY)
  add_test(NAME ITER-rfile-modify-hdf5 COMMAND run_efit_hdf5.sh)
  add_test(NAME ITER-rfile-modify-hdf5-results COMMAND check_efit_results.sh)
  set(EFIT_MODE 2)
endif ()
set(RFILE_TEST False)
