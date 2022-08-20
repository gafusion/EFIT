
include_directories(
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/io
  ${CMAKE_BINARY_DIR}/io
)

# This calls out to routines needed in d3lib
set(EFITLIB_SOURCES
  modules-efit.F90 set_eparm.F90 psical.f90 mat_solve.f90 
  tables.F90 read_namelist.F90 set_defaults.f90
  buneman.f90 set_filenames.f90 ecom1-mods.f90 zpline.f90 
  ecom2-mods.f90 curve2d_mod.f90 data_input.F90 set_init.f90
  response_matrix.F90 shapesurf.F90 setup_data_fetch.F90 current.f90
  external_current.F90 boundary.F90 spline.f90 pltdat.F90
  write_a.f90 write_g.f90 write_m.F90 write_ot.f90 print_terminal.F90
  ppbasisfunc.f90 ffbasisfunc.f90 wwbasisfunc.f90 eebasisfunc.f90
  get_exp_data.f90 get_kinetic_data.f90 auto_knot.F90
  ece.F90 msels.F90 stark.F90 chkerr.f90
  cyclic.F90 set_basis.f90 pressure.F90 utils.F90
  green.F90 fit.F90 beta_li.F90 pflux.f90 update_parameters.F90
  ${CMAKE_CURRENT_BINARY_DIR}/set_extvars.F90
)
#set_extvars.F90 is generated -- see top level CMakeLists.txt
if(${HAVE_SNAP})
  set(EFITLIB_SOURCES
    ${EFITLIB_SOURCES}
    getdat.F90 measurements.F90 diamagnetic.f90
    write_k.F90
  )
endif()

add_library(efitlib ${EFITLIB_SOURCES})
add_dependencies(efitlib config_f)
target_link_libraries(efitlib green)

# Diagnostic libs
# snap mode is only set up for DIIID currently and requires 2
# additional external libraries
if(${HAVE_SNAP})
  target_link_libraries(efitlib d3lib)
  if(${HAVE_MDSPLUS})
    if(${HAVE_MSE})
      target_link_libraries(efitlib cerlib mselib)
    else()
      target_link_libraries(efitlib cerlib)
    endif()
  endif()
#else()
#  target_link_libraries(efitlib ptdatalib)
endif()
#if(NOT ${HAVE_MDSPLUS})
#  target_link_libraries(efitlib mdspluslib)
#endif()

# Math libs
target_link_libraries(efitlib lsode r8slatec)
# IO libs
target_link_libraries(efitlib dissplalib)
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
