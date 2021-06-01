
include_directories(
  ${CMAKE_BINARY_DIR}
)

set(GREEN_SOURCES
    exparm.f90 bfgrid.f90 cacoil.f90 cecoil.f90 filech.f90 
    coilsp.f90 consta.f90 cvesel.f90 fcoil.f90 vas_sub1.f90
    fshift.f90 input.f90 mprobe.f90 nio.f90
    pmodel.f90 rogowl.f90 siloop.f90 efund_read.f90 dprobe.f90
)

add_library(green ${GREEN_SOURCES})
install(TARGETS green DESTINATION lib
  PERMISSIONS OWNER_READ OWNER_EXECUTE OWNER_WRITE
              GROUP_READ GROUP_EXECUTE GROUP_WRITE
              WORLD_READ WORLD_EXECUTE
)

add_executable(efund  efundud.f90)
target_link_libraries(efund green)


# Add command and target to make green functions
# e.g., `make makeD3DGreensFunctions`  
# (and you can also `make help | grep -i green` to find it)
add_custom_target(makeD3DGreensFunctions)
add_custom_command(TARGET makeD3DGreensFunctions
    COMMAND ${CMAKE_COMMAND} -E copy 
        ${CMAKE_SOURCE_DIR}/share/data/mhdin_diii-d.dat 
        ${CMAKE_BINARY_DIR}/mhdin.dat
    COMMAND /bin/rm -f dprobe.dat
    COMMAND green/efund
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# This is to make sure the input file is in the right place for the above
# command
#file(COPY ${CMAKE_SOURCE_DIR}/share/data/mhdin_diii-d.dat 
#     DESTINATION ${CMAKE_BINARY_DIR}/)
#file(RENAME ${CMAKE_BINARY_DIR}/mhdin_diii-d.dat ${CMAKE_BINARY_DIR}/mhdin.dat)


