######################################################################
#
# CMakeLists.txt for efit tests
#
######################################################################

set(CURRENT_MACHINE_LINK_DIR ${CMAKE_CURRENT_BINARY_DIR}/)
if (${TEST_EFUND})
  set(CURRENT_MACHINE_BUILT_DIR ${CMAKE_CURRENT_BINARY_DIR}/green/)
endif ()

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/lim.dat 
     DESTINATION ${CURRENT_MACHINE_LINK_DIR})

add_subdirectory(green)
add_subdirectory(efit01)
add_subdirectory(efit02)
