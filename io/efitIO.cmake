######################################################################
#
# CMake to find IO dependencies
#
######################################################################

set(io_libs "")

#---------------------------------------------------------------------
#
# Find HDF5
# https://cmake.org/cmake/help/latest/module/FindHDF5.html
#
#---------------------------------------------------------------------
option(ENABLE_HDF5 "Enable HDF5" off)
set(HAVE_HDF5 False)   # Used in defines
if (${ENABLE_HDF5})
    set(HDF5_USE_STATIC_LIBRARIES TRUE)
    find_package(HDF5 COMPONENTS Fortran)
#    find_package(HDF5 COMPONENTS Fortran HL)
    if (${HDF5_FOUND})
      message(STATUS "Found HDF5")
      set(HAVE_HDF5 True)
      set(io_libs ${HDF5_LIBRARIES} ${Z_LIBRARIES} ${io_libs})
#      include_directories(${HDF5_INCLUDE_DIRS})
      if (UNIX)
        # Required for hdf5 version 1.8.11 and greater
        set(io_libs ${io_libs} ${CMAKE_DL_LIBS}) 
      endif ()
    endif()
endif()


#---------------------------------------------------------------------
#
# Find Netcdf
# https://cmake.org/cmake/help/latest/command/find_package.html
#
#---------------------------------------------------------------------
option(ENABLE_NETCDF "Enable NetCDF" off)
set(HAVE_NETCDF False)   # Used in defines
if (${ENABLE_NETCDF})
    set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_MODULE_PATH})
    find_package (NetCDF)
    if (${NetCDF_FOUND})
      message(STATUS "Found NetCDF")
      set(HAVE_NETCDF True)
      set(USE_NETCDF 1 CACHE BOOL "Whether NETCDF is linked")  # Used in directives
      include_directories(${NetCDF_INCLUDE_DIR})
      set(io_libs ${NetCDF_LIBRARIES} ${io_libs})
    endif()
endif()

#---------------------------------------------------------------------
#
# Find MDSplus
# https://cmake.org/cmake/help/latest/command/find_package.html
#
#---------------------------------------------------------------------
option(ENABLE_MDSPLUS "Enable MDS+" off)
set(HAVE_MDSPLUS False)   # Used in defines
if (${ENABLE_MDSPLUS})
  find_package(MDSPlus COMPONENTS 
               INSTALL_DIR "mdsplus"
               HEADERS "mdsdescrip.h" "mdslib.h"
               LIBRARIES "MdsLib"
               INCLUDE_SUBDIRS "include"
               LIBRARY_SUBDIRS "lib"
               )
  if (${MDSPLUS_FOUND})
    message(STATUS "Found MDS+")
    set(HAVE_MDSPLUS True)
    set(io_libs ${MDSPLUS_LIBRARIES} ${io_libs})
    # This may only be necessary for DIIID...
    if(NOT EXISTS ${MSE_LIB})
      message(STATUS "MSE_LIB not found, DIIID efit02 snap files will fail")
    else()
      if(NOT EXISTS ${D3_LIB})
        message(STATUS "D3_LIB not found, MSE is not being linked")
      else()
        set(USE_MSE TRUE)                 # ifdef
      endif()
    endif()
  endif()
endif()
