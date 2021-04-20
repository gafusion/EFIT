
include_directories(
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/io
)

# This calls out to routines needed in d3lib
set(EFITLIB_SOURCES
  modules-efit.f90 get_eparmdud.f90 psicald.f90 lsolved.f90 
  subs_d_pt.f90 efit_tables.f90 
  getdat.F90 getdiau.f90 
  bunemad.f90 getfnmdud.f90 ecom1-mods.f90 zplined.f90 
  ecom2-mods.f90 vas_sub1.f90 get_opt_input.F90 data_input.f90
  matrix.f90 shapesurf.f90 getsets.F90 prtout.F90 writeout.F90 currnt.f90
  fcurrt.f90 boundud129.f90 chkeud129.f90 splineud129.f90 weqdud129.f90
  msels_data.f90 msels_hist.f90 pltdatud129.f90 getecdud129.F90 getneud129.F90
  ppbasisfuncud129.f90 ffbasisfuncud129.f90 wwbasisfuncud129.f90 
  expdataud129.f90 
  efitnotifier.c
)
# TODO:  This has some strange code that I need to sort out. 
# tearingd129.f90

add_library(efitlib ${EFITLIB_SOURCES})
target_link_libraries(efitlib green)
# Diagnostic libs
target_link_libraries(efitlib cerlib)
#target_link_libraries(efitlib cerlib d3lib)
#target_link_libraries(efitlib cerlib mselib)
# Math libs
target_link_libraries(efitlib lsode r8slatec)
# IO libs
target_link_libraries(efitlib dissplalib ptdatalib)
if (NOT ${HAVE_NETCDF})
target_link_libraries(efitlib netcdflib)
endif()
if (NOT ${HAVE_MDSPLUS})
target_link_libraries(efitlib mdspluslib)
endif()
if (${HAVE_HDF5})
target_link_libraries(efitlib vshdf5lib)
endif()

install(TARGETS efitlib DESTINATION lib
  PERMISSIONS OWNER_READ OWNER_EXECUTE OWNER_WRITE
              GROUP_READ GROUP_EXECUTE GROUP_WRITE
              WORLD_READ WORLD_EXECUTE
)

#if (${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
#elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Lahey")
#endif ()

add_executable(efit  efitdud129.F90)
message(STATUS "math_libs: " ${math_libs})
message(STATUS "io_libs: " ${io_libs})
target_link_libraries(efit efitlib ${io_libs} ${math_libs})
