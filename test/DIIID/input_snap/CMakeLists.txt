######################################################################
#
# CMakeLists.txt for efit tests
#
######################################################################

set(RUNID 186610)
# This could be changed to mode 5 for a more efficient test when outputs are compared in detail
set(EFIT_MODE 7) # Mode 3 should also work
set(SNAP_EXT "EFIT01_186610")
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
# not used
set(EFIT_IN "")

foreach (t ${TIMES})
  set(EFIT_LABEL ${EFIT_LABEL} ${RUNID}.${t})
endforeach ()
string(REPLACE ";" " " EFIT_LABEL "${EFIT_LABEL}")

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/efit_snap.dat_${SNAP_EXT}
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

configure_file(${EFIT_INPUT_CONFIGFILE}
               ${CMAKE_CURRENT_BINARY_DIR}/efit.input @ONLY)
configure_file(${EFIT_TEST_CONFIGFILE_INPUT}
               ${CMAKE_CURRENT_BINARY_DIR}/run_efit_input.sh @ONLY)
configure_file(${EFIT_RESULTS_CONFIGFILE_SNAP}
               ${CMAKE_CURRENT_BINARY_DIR}/check_efit_snap_results.sh @ONLY)

add_test(NAME DIIID-input-snap COMMAND run_efit_input.sh)
add_test(NAME DIIID-input-snap-results COMMAND check_efit_snap_results.sh)
