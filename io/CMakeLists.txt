######################################################################          
#                                                                               
# CMakeLists.txt for external io libraries
#                                                                               
######################################################################

# Search and replace variables in the efit version file
set(REPLACE_VARS "CMAKE_Fortran_COMPILER_ID;CMAKE_Fortran_COMPILER_VERSION;CMAKE_Fortran_FLAGS;CMAKE_Fortran_BUILD_TYPE_FLAGS;GIT_BRANCH;GIT_COMMIT_HASH;GIT_URL;GIT_MODIFICATIONS")
file(READ efit_version.in.F90 FILE_CONTENT)
foreach(REPLACE_VAR ${REPLACE_VARS})
   string(LENGTH "${${REPLACE_VAR}}" SLEN)
   set(LINELEN 45)
   #message(STATUS "${REPLACE_VAR} ${${REPLACE_VAR}} ${SLEN}")
   if(${SLEN} GREATER ${LINELEN})
     set(STRSTART 0)
     set(STREND ${LINELEN})
     while(${STRSTART} LESS ${SLEN})
       string(SUBSTRING "${${REPLACE_VAR}}" ${STRSTART} ${LINELEN} SUBSTR)
       if(${STREND} GREATER ${SLEN})
         math(EXPR LINELEN "${SLEN} - ${STRSTART}")
         set(SUBSTR "${SUBSTR}")
       else()
         set(SUBSTR "${SUBSTR}\"// &\n     &  \"${REPLACE_VAR}")
       endif()
       string(REGEX REPLACE "${REPLACE_VAR}" "${SUBSTR}" FILE_CONTENT "${FILE_CONTENT}")
       math(EXPR STRSTART "${STRSTART} + ${LINELEN}")
       math(EXPR STREND "${STRSTART} + ${LINELEN}")
     endwhile()
   else()
     string(REGEX REPLACE "${REPLACE_VAR}" "${${REPLACE_VAR}}" FILE_CONTENT "${FILE_CONTENT}")
   endif()
endforeach()
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/efit_version.F90 "${FILE_CONTENT}")

# Compile efitversion library
add_library(efitversion efit_version.F90)
install(TARGETS efitversion DESTINATION lib
  PERMISSIONS OWNER_READ OWNER_EXECUTE OWNER_WRITE
              GROUP_READ GROUP_EXECUTE GROUP_WRITE
              WORLD_READ WORLD_EXECUTE
)
target_link_libraries(efitversion ${io_libs})
set(io_libs efitversion ${io_libs})

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
