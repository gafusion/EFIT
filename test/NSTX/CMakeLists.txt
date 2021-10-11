######################################################################
#
# CMakeLists.txt for efit tests
#
######################################################################

# TO DO: add test option that skips efund

set(CURRENT_MACHINE_LINK_DIR ${CMAKE_CURRENT_BINARY_DIR}/)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/lim.dat 
     DESTINATION ${CURRENT_MACHINE_LINK_DIR})

add_subdirectory(green)
#add_subdirectory(efit01)
add_subdirectory(efit02)
#add_subdirectory(rfile)
