######################################################################
#
# CMakeLists.txt for efit tests
#
######################################################################

set(RUNID 186610)
set(FIRST_TIME 2400)
set(STEP_SIZE 100)
if (${ENABLE_PARALLEL})
  set(TIMES 02400 02500 02600 02700)
  set(EFIT_COUNT 4)
else ()
  set(TIMES 02400)
endif ()
if (${HAVE_NETCDF})
  set(outputs "g a k m")
else ()
  set(outputs "g a k")
endif ()

foreach (t ${TIMES})
  set(EFIT_LABEL ${EFIT_LABEL} ${RUNID}.${t})
endforeach ()
string(REPLACE ";" " " EFIT_LABEL "${EFIT_LABEL}")

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/efit_snap.dat
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

set(EFIT_MODE 5)
if (${ENABLE_PARALLEL})
  set(NPROC_SAVE ${NPROC})
  set(NPROC 1)
endif ()
configure_file(${EFIT_TEST_CONFIGFILE_SNAP}
               ${CMAKE_CURRENT_BINARY_DIR}/run_efit_kfile.sh @ONLY)
configure_file(${EFIT_RESULTS_CONFIGFILE_SNAP}
               ${CMAKE_CURRENT_BINARY_DIR}/check_efit_kfile_results.sh @ONLY)
if (${ENABLE_PARALLEL})
  set(NPROC ${NPROC_SAVE})
endif ()

add_test(NAME DIIID-efit01-kfile COMMAND run_efit_kfile.sh)
add_test(NAME DIIID-efit01-kfile-results COMMAND check_efit_kfile_results.sh)

set(EFIT_MODE 3)
configure_file(${EFIT_TEST_CONFIGFILE_SNAP}
               ${CMAKE_CURRENT_BINARY_DIR}/run_efit_snap.sh @ONLY)
configure_file(${EFIT_RESULTS_CONFIGFILE_SNAP}
               ${CMAKE_CURRENT_BINARY_DIR}/check_efit_snap_results.sh @ONLY)

add_test(NAME DIIID-efit01-snap COMMAND run_efit_snap.sh)
add_test(NAME DIIID-efit01-snap-results COMMAND check_efit_snap_results.sh)
