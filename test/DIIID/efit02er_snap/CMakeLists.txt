######################################################################
#
# CMakeLists.txt for efit tests
#
######################################################################

set(RUNID 186610)
set(TIME 2400)
set(STEP 100)
if (${ENABLE_PARALLEL})
  set(TIMES 02400 02500 02600 02700)
  set(EFIT_TIME 4)
else()
  set(TIMES 02400)
endif()

foreach( t ${TIMES})
  set(EFIT_LABEL ${EFIT_LABEL} ${RUNID}.${t})
endforeach()
string(REPLACE ";" " " EFIT_LABEL "${EFIT_LABEL}")
string(REPLACE ";" "\n" EFIT_IN "${EFIT_IN}")

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/efit_snap.dat
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

set(EFIT_MODE 5)
configure_file(${EFIT_TEST_CONFIGFILE_SNAP}
               ${CMAKE_CURRENT_BINARY_DIR}/run_efit_kfile.sh @ONLY)
configure_file(${EFIT_RESULTS_CONFIGFILE_SNAP}
               ${CMAKE_CURRENT_BINARY_DIR}/check_efit_kfile_results.sh @ONLY)

set(EFIT_MODE 3)
add_test(NAME DIIID-efit02er-kfile COMMAND run_efit_kfile.sh)
add_test(NAME DIIID-efit02er-kfile-results COMMAND check_efit_kfile_results.sh)

configure_file(${EFIT_TEST_CONFIGFILE_SNAP}
               ${CMAKE_CURRENT_BINARY_DIR}/run_efit_snap.sh @ONLY)
configure_file(${EFIT_RESULTS_CONFIGFILE_SNAP}
               ${CMAKE_CURRENT_BINARY_DIR}/check_efit_snap_results.sh @ONLY)
add_test(NAME DIIID-efit02er-snap COMMAND run_efit_snap.sh)
add_test(NAME DIIID-efit02er-snap-results COMMAND check_efit_snap_results.sh)
set(EFIT_MODE 2)
