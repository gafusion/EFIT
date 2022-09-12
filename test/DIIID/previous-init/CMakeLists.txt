######################################################################
#
# CMakeLists.txt for efit tests
#
######################################################################

set(RUNID 186610)
if (${ENABLE_PARALLEL})
  set(TIMES 02400 02500 02600 02700 02800 02900 03000 03100)
  set(EFIT_COUNT 8)
else ()
  set(TIMES 02400 02500)
  set(EFIT_COUNT 2)
endif ()

foreach (t ${TIMES})
  set(EFIT_LABEL ${EFIT_LABEL} ${RUNID}.${t})
  set(EFIT_IN ${EFIT_IN} k${RUNID}.${t})
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/k${RUNID}.${t}
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
endforeach ()
string(REPLACE ";" " " EFIT_LABEL "${EFIT_LABEL}")
string(REPLACE ";" "\n" EFIT_IN "${EFIT_IN}")

configure_file(${EFIT_TEST_CONFIGFILE}
               ${CMAKE_CURRENT_BINARY_DIR}/run_efit.sh @ONLY)
configure_file(${EFIT_RESULTS_CONFIGFILE}
               ${CMAKE_CURRENT_BINARY_DIR}/check_efit_results.sh @ONLY)

add_test(NAME DIIID-previous-init COMMAND run_efit.sh)
add_test(NAME DIIID-previous-init-results COMMAND check_efit_results.sh)

if (${ENABLE_HDF5})
  set(EFIT_MODE 1)
  if (${ENABLE_PARALLEL})
    set(EFIT_IN ${RUNID}_mpi.h5)
  else ()
    set(EFIT_IN ${RUNID}.h5)
  endif ()
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${EFIT_IN}
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
  configure_file(${EFIT_TEST_CONFIGFILE_HDF5}
                 ${CMAKE_CURRENT_BINARY_DIR}/run_efit_hdf5.sh @ONLY)
  add_test(NAME DIIID-previous-init-hdf5 COMMAND run_efit_hdf5.sh)
  add_test(NAME DIIID-previous-init-hdf5-results COMMAND check_efit_results.sh)
endif ()
