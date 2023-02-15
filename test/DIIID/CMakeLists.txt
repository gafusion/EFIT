######################################################################
#
# CMakeLists.txt for efit tests
#
######################################################################

set(CURRENT_MACHINE_LINK_DIR ${CMAKE_CURRENT_BINARY_DIR}/)
if (${TEST_EFUND})
  set(CURRENT_MACHINE_BUILT_DIR ${CMAKE_CURRENT_BINARY_DIR}/green/)
endif ()

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
add_subdirectory(vacuum)
add_subdirectory(esave)
add_subdirectory(restart)
add_subdirectory(previous-init)
add_subdirectory(input)
add_subdirectory(integral_method)
add_subdirectory(ece)
add_subdirectory(vsloop)
add_subdirectory(fitvs)
if (${HAVE_SNAP})
  add_subdirectory(efit01_snap)
  if (${HAVE_MSE})
    add_subdirectory(efit02_snap)
    add_subdirectory(efit02er_snap)
  endif ()
  add_subdirectory(input_snap)
endif ()
