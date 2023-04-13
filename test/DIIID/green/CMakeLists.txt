######################################################################
#
# CMakeLists.txt for efit test
#
######################################################################

if (${TEST_EFUND})
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/mhdin.dat 
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

  configure_file(${EFUND_TEST_CONFIGFILE} 
                 ${CMAKE_CURRENT_BINARY_DIR}/run_green.sh @ONLY)
  configure_file(${EFUND_RESULTS_CONFIGFILE}
                 ${CMAKE_CURRENT_BINARY_DIR}/check_green_results.sh @ONLY)

  add_test(NAME DIIID-green COMMAND run_green.sh)
  add_test(NAME DIIID-green-results COMMAND check_green_results.sh)
else ()
  if (NOT ${GRID_SIZE} EQUAL 65)
    message(FATAL_ERROR "Green Function Tables of this size must be made")
  endif ()

  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/mhdin.dat
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/output/ec6565.ddd
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/output/ep6565.ddd
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/output/fc6565.ddd
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/output/ra6565.ddd
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/output/re6565.ddd
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/output/rfcoil.ddd
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/output/rv6565.ddd
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
endif ()
