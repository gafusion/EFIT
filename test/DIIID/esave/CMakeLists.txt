######################################################################
#
# CMakeLists.txt for efit tests
#
######################################################################

set(RUNID 186610)
if (${ENABLE_PARALLEL})
  set(TIMES 02400 02500 02600 02700)
  set(EFIT_COUNT 4)
  set(STEPS 00000 00001 00002 00003)
else ()
  set(TIMES 02400)
  set(STEPS 00000)
endif ()

foreach (t ${TIMES})
  set(EFIT_LABEL ${EFIT_LABEL} ${RUNID}.${t})
  set(EFIT_IN ${EFIT_IN} k${RUNID}.${t})
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/k${RUNID}.${t}
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
endforeach ()
string(REPLACE ";" " " EFIT_LABEL "${EFIT_LABEL}")
string(REPLACE ";" "\n" EFIT_IN "${EFIT_IN}")
foreach (s ${STEPS})
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/esave${s}.dat
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
endforeach ()
if (NOT ${ENABLE_PARALLEL})
  file(RENAME ${CMAKE_CURRENT_BINARY_DIR}/esave${STEPS}.dat
       ${CMAKE_CURRENT_BINARY_DIR}/esave.dat)
endif ()

configure_file(${EFIT_TEST_CONFIGFILE}
               ${CMAKE_CURRENT_BINARY_DIR}/run_efit.sh @ONLY)
configure_file(${EFIT_RESULTS_CONFIGFILE}
               ${CMAKE_CURRENT_BINARY_DIR}/check_efit_results.sh @ONLY)

add_test(NAME DIIID-esave COMMAND run_efit.sh)
add_test(NAME DIIID-esave-results COMMAND check_efit_results.sh)

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
  add_test(NAME DIIID-esave-hdf5 COMMAND run_efit_hdf5.sh)
  add_test(NAME DIIID-esave-hdf5-results COMMAND check_efit_results.sh)
endif ()
