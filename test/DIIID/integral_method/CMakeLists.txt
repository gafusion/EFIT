######################################################################
#
# CMakeLists.txt for efit tests
#
######################################################################

set(RUNID 186610)
if (${ENABLE_PARALLEL})
  set(TIMES 02400 02500 02600 02700)
  set(EFIT_COUNT 4)
else ()
  set(TIMES 02400)
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

add_test(NAME DIIID-integral-method COMMAND run_efit.sh)
add_test(NAME DIIID-integral-method-results COMMAND check_efit_results.sh)
