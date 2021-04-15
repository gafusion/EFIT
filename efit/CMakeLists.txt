
include_directories(
  ${CMAKE_BINARY_DIR}/nimloc
)

set(EFITLIB_SOURCES
  modules-efit.f90 get_eparmdud.f03 psicald.f90 lsolved.f90 
  # subs_d_pt_linux.f90  # Unique fortran technique
  subs_d_pt.f90
  bunemad.f90 getdiau.f90 getfnmdud.f90 ecom1-mods.f90 zplined.f90 dgglse.f90
  ecom2-mods.f90 vas_sub1.f90 get_opt_input.f03 data_input.f03
  matrix.f90 shapesurf.f90 getsets.f03 prtout.F90 writeout.f03 currnt.f90
  fcurrt.f90 boundud129.f90 chkeud129.f90 splineud129.f90 weqdud129.f03
  msels_data.f90 msels_hist.f90 pltdatud129.f90 getecdud129.f03 getneud129.f90
  ppbasisfuncud129.f90 ffbasisfuncud129.f90 wwbasisfuncud129.f90 tearingd129.f90
  expdataud129.f90 
  efitnotifier.c
)

add_library(efitlib ${EFITLIB_SOURCES})
target_link_libraries(efitlib green cerlib)
install(TARGETS efitlib DESTINATION lib
  PERMISSIONS OWNER_READ OWNER_EXECUTE OWNER_WRITE
              GROUP_READ GROUP_EXECUTE GROUP_WRITE
              WORLD_READ WORLD_EXECUTE
)

set(aux_libs)

if(Linux)
add_library(lin-dp lin-dp.f90)
set(aux_libs ${aux_libs} lin-dp)
endif()

#if (${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
#elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Lahey")
#endif ()

add_executable(efit  efitdud129.f90)
target_link_libraries(efit efitlib ${aux_libs})
