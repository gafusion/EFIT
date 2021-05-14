
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
