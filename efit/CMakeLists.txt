
include_directories(
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/io
)

# This calls out to routines needed in d3lib
set(EFITLIB_SOURCES
  modules-efit.F90 get_eparmdud.f90 psicald.f90 lsolved.f90 
  subs_d_pt.F90 efit_tables.f90 
  getdat.F90 getdiau.f90 read_namelist.F90
  bunemad.f90 getfnmdud.f90 ecom1-mods.f90 zplined.f90 
  ecom2-mods.f90 vas_sub1.f90 data_input.f90
  matrix.f90 shapesurf.f90 getsets.F90 prtout.F90 writeout.F90 currnt.f90
  fcurrt.f90 boundud129.f90 chkeud129.f90 splineud129.f90 weqdud129.f90
  msels_data.f90 msels_hist.f90 pltdatud129.f90 getecdud129.F90 getneud129.F90
  ppbasisfuncud129.f90 ffbasisfuncud129.f90 wwbasisfuncud129.f90 
  expdataud129.f90 
  cyclic.F90 efit_bdata.F90 set_routines.F90 pressure.F90 get_routines.F90
  efitnotifier.c
  ${CMAKE_CURRENT_BINARY_DIR}/set_expath.F90
)
# stubs are superior to ifdefs in this case.
if (${ENABLE_PARALLEL})
    set(EFITLIB_SOURCES ${EFITLIB_SOURCES} mpi_parallel.f90)
else()
    set(EFITLIB_SOURCES ${EFITLIB_SOURCES} mpi_serial.f90)
endif()


#set_expath.F90 is generated -- see top level CMakeLists.txt

# TODO:  This has some strange code that I need to sort out. 
# tearingd129.f90

add_library(efitlib ${EFITLIB_SOURCES})
add_dependencies(efitlib config_f)
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

add_executable(efit  efitdud129.F90 get_opt_input.F90)
target_link_libraries(efit efitlib ${io_libs} ${math_libs})

#add_executable(efitai  efitai.F90 get_opt_input_args.F90)
#add_executable(efitai  efitai.F90)
#target_link_libraries(efitai efitlib ${io_libs} ${math_libs})
