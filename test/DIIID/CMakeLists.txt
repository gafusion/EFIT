######################################################################
#
# CMakeLists.txt for efit tests
#
######################################################################

# TODO: add test option that skips efund

set(CURRENT_MACHINE_LINK_DIR ${CMAKE_CURRENT_BINARY_DIR}/)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/btcomp.dat 
     DESTINATION ${CURRENT_MACHINE_LINK_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/ccoil.ddd 
     DESTINATION ${CURRENT_MACHINE_LINK_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/ccomp.dat 
     DESTINATION ${CURRENT_MACHINE_LINK_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/dcoef.dat 
     DESTINATION ${CURRENT_MACHINE_LINK_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/fitweight.dat 
     DESTINATION ${CURRENT_MACHINE_LINK_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/icomp.dat 
     DESTINATION ${CURRENT_MACHINE_LINK_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/lim.dat 
     DESTINATION ${CURRENT_MACHINE_LINK_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/n1coil.ddd 
     DESTINATION ${CURRENT_MACHINE_LINK_DIR})

add_subdirectory(green)
add_subdirectory(efit01)
add_subdirectory(efit02)
add_subdirectory(efit02er)
add_subdirectory(rfile)
add_subdirectory(kineticEFIT)
add_subdirectory(varyped)
