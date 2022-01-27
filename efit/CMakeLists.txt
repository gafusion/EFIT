
include_directories(
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/io
  ${CMAKE_BINARY_DIR}/io
)

# This calls out to routines needed in d3lib
set(EFITLIB_SOURCES
  modules-efit.F90 get_eparmdud.F90 psical.f90 lsolved.f90 
  subs_d_pt.F90 efit_tables.F90 
  getdat.F90 getdia.f90 read_namelist.F90
  bunemad.f90 getfnmdud.f90 ecom1-mods.f90 zplined.f90 
  ecom2-mods.f90 curve2d_mod.f90 vas_sub1.f90 data_input.F90
  matrix.F90 shapesurf.F90 getsets.F90 prtout.F90 writeout.F90 currnt.f90
  fcurrt.F90 bound.F90 chke.f90 spline.f90 weq.f90
  msels_data.f90 msels_hist.f90 pltdat.F90 getecd.F90 getne.F90
  ppbasisfunc.f90 ffbasisfunc.f90 wwbasisfunc.f90 
  expdata.f90  autoknot.F90 efit_utils.F90
  cyclic.F90 set_routines.F90 pressure.F90 get_routines.F90
  green.F90 fit.F90 beta_li.F90 flux_current.F90 weight.F90
  efitnotifier.c wmeasure.F90
  ${CMAKE_CURRENT_BINARY_DIR}/set_expath.F90
)

#set_expath.F90 is generated -- see top level CMakeLists.txt

# TODO:  This has some strange code that I need to sort out. 
# tearing.f90

add_library(efitlib ${EFITLIB_SOURCES})
add_dependencies(efitlib config_f)
target_link_libraries(efitlib green)
# Diagnostic libs
target_link_libraries(efitlib cerlib)
if(${HAVE_MDSPLUS})
  if(DEFINED USE_MDS)
    add_library(mselib STATIC IMPORTED)
    set_target_properties(mselib PROPERTIES IMPORTED_LOCATION ${MSE_LIB})
    add_library(d3lib STATIC IMPORTED)
    set_target_properties(d3lib PROPERTIES IMPORTED_LOCATION ${D3_LIB})
    target_link_libraries(efitlib mselib d3lib)
  else()
    target_link_libraries(efitlib ptdatalib)
  endif()
else()
  target_link_libraries(efitlib mdspluslib ptdatalib)
endif()
# Math libs
target_link_libraries(efitlib lsode r8slatec)
# IO libs
target_link_libraries(efitlib efitversion dissplalib)
if(NOT ${HAVE_NETCDF})
  target_link_libraries(efitlib netcdflib)
endif()
if (${HDF5_FOUND})
  target_link_libraries(efitlib hdf5lib)
endif()

install(TARGETS efitlib DESTINATION lib
  PERMISSIONS OWNER_READ OWNER_EXECUTE OWNER_WRITE
              GROUP_READ GROUP_EXECUTE GROUP_WRITE
              WORLD_READ WORLD_EXECUTE
)

#if(${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
#elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "Lahey")
#endif()

add_executable(efit efit.F90 get_opt_input.F90)
# TODO: linking could be done more cleanly
#target_link_libraries(efit efitlib)
target_link_libraries(efit efitlib ${io_libs} ${math_libs})

# Premature...
#add_executable(efitai efitai.F90 get_opt_input_args.F90)
#add_executable(efitai efitai.F90)
#target_link_libraries(efitai efitlib ${io_libs} ${math_libs})
