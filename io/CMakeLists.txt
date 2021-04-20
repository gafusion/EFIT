######################################################################          
#                                                                               
# CMakeLists.txt for extlibs
#                                                                               
######################################################################


if (${HDF5_FOUND})
    include_directories(
      ${HDF5_INCLUDE_DIRS}
      ${HDF5_INCLUDE_DIRS}/static
    )
    add_library(vshdf5lib hdf5_fcapi.F)
    install(TARGETS iolib DESTINATION lib
      PERMISSIONS OWNER_READ OWNER_EXECUTE OWNER_WRITE
                  GROUP_READ GROUP_EXECUTE GROUP_WRITE
                  WORLD_READ WORLD_EXECUTE
    )
    target_link_libraries(vshdf5lib ${io_libs})
    set(io_libs vshdf5lib ${io_libs})
endif()

add_library(dissplalib  disspla_dummy.F)
install(TARGETS dissplalib DESTINATION lib
  PERMISSIONS OWNER_READ OWNER_EXECUTE OWNER_WRITE
              GROUP_READ GROUP_EXECUTE GROUP_WRITE
              WORLD_READ WORLD_EXECUTE
)
add_library(ptdatalib  ptdata_dummy.F)
install(TARGETS ptdatalib DESTINATION lib
  PERMISSIONS OWNER_READ OWNER_EXECUTE OWNER_WRITE
              GROUP_READ GROUP_EXECUTE GROUP_WRITE
              WORLD_READ WORLD_EXECUTE
)


if (NOT ${HAVE_NETCDF})
    add_library(netcdflib  netcdf_dummy.F)
    install(TARGETS netcdflib DESTINATION lib
      PERMISSIONS OWNER_READ OWNER_EXECUTE OWNER_WRITE
                  GROUP_READ GROUP_EXECUTE GROUP_WRITE
                  WORLD_READ WORLD_EXECUTE
    )
    set(io_libs netcdflib ${io_libs})
endif ()

if (NOT ${HAVE_MDSPLUS})
    add_library(mdspluslib mdsfdummy.f)
    install(TARGETS mdspluslib DESTINATION lib
      PERMISSIONS OWNER_READ OWNER_EXECUTE OWNER_WRITE
                  GROUP_READ GROUP_EXECUTE GROUP_WRITE
                  WORLD_READ WORLD_EXECUTE
    )
    set(io_libs mdspluslib ${io_libs})
endif ()
