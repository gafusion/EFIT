######################################################################          
#                                                                               
# CMakeLists.txt for external io libraries
#                                                                               
######################################################################

# Compile HDF5 wrappers library
if (${HDF5_FOUND})
  include_directories(
    ${HDF5_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}/static
  )
  add_library(hdf5lib hdf5_fcapi.f90)
  install(TARGETS hdf5lib DESTINATION lib
    PERMISSIONS OWNER_READ OWNER_EXECUTE OWNER_WRITE
                GROUP_READ GROUP_EXECUTE GROUP_WRITE
                WORLD_READ WORLD_EXECUTE
  )
  target_link_libraries(hdf5lib ${io_libs})
  set(io_libs hdf5lib ${io_libs})
endif()

# Make dummy plotting routines library
# TODO: will this code be resurected or can it be removed?
add_library(dissplalib  disspla_dummy.F)
install(TARGETS dissplalib DESTINATION lib
  PERMISSIONS OWNER_READ OWNER_EXECUTE OWNER_WRITE
              GROUP_READ GROUP_EXECUTE GROUP_WRITE
              WORLD_READ WORLD_EXECUTE
)

# Make dummy netcdf routines if not linked
if (NOT ${HAVE_NETCDF})
  add_library(netcdflib  netcdf_dummy.F)
  install(TARGETS netcdflib DESTINATION lib
    PERMISSIONS OWNER_READ OWNER_EXECUTE OWNER_WRITE
                GROUP_READ GROUP_EXECUTE GROUP_WRITE
                WORLD_READ WORLD_EXECUTE
  )
  set(io_libs netcdflib ${io_libs})
endif ()

# Make dummy MDS+ routines if not linked
if (NOT ${HAVE_MDSPLUS})
  add_library(mdspluslib mdsfdummy.f)
  install(TARGETS mdspluslib DESTINATION lib
    PERMISSIONS OWNER_READ OWNER_EXECUTE OWNER_WRITE
                GROUP_READ GROUP_EXECUTE GROUP_WRITE
                WORLD_READ WORLD_EXECUTE
  )
  set(io_libs mdspluslib ${io_libs})
endif ()

# Make dummy ptdata routines if not linked
if (NOT DEFINED USE_MDS)
  add_library(ptdatalib  ptdata_dummy.F)
  install(TARGETS ptdatalib DESTINATION lib
    PERMISSIONS OWNER_READ OWNER_EXECUTE OWNER_WRITE
                GROUP_READ GROUP_EXECUTE GROUP_WRITE
                WORLD_READ WORLD_EXECUTE
  )
endif()
