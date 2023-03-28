######################################################################
#
# CMake to find IO dependencies
#
######################################################################

set(io_libs "")

# Ideally, we would use the HDF5 installed cmake file, but in practice you
# can't do that with  most distros.  Rather than find_package all of those 
# dependencies, this does a simpler method based on find_library
# Based on a snippet found at the kitware gitlab site
macro(getAdditionalHdf5Libs)
  set(_additional_libs sz z)
  foreach(_additional_lib IN LISTS _additional_libs)
    if(HDF5_USE_STATIC_LIBRARIES)
      set(_libnames ${_additional_lib} lib${_additional_lib}.a)
    else()
      set(_libnames ${_additional_lib})
    endif(HDF5_USE_STATIC_LIBRARIES)
    set(_libvar "LIB_${_additional_lib}")
    find_library(${_libvar}
      NAMES ${_libnames}
      HINTS ${HDF5_ROOT}
      PATH_SUFFIXES lib Lib)
    if(NOT ${${_libvar}} STREQUAL "${_libvar}-NOTFOUND")
      list(APPEND HDF5_LIBRARIES ${${_libvar}})
    endif()
  endforeach()
list(REMOVE_DUPLICATES HDF5_LIBRARIES)
endmacro(getAdditionalHdf5Libs)

#---------------------------------------------------------------------
#
# Find HDF5
# https://cmake.org/cmake/help/latest/module/FindHDF5.html
#
#---------------------------------------------------------------------
option(ENABLE_HDF5 "Enable HDF5" off)
if(${ENABLE_HDF5})
  set(HDF5_USE_STATIC_LIBRARIES TRUE)
  find_package(HDF5 COMPONENTS Fortran HL)
  if(${HDF5_FOUND})
    getAdditionalHdf5Libs()
    set(USE_HDF5 1 CACHE BOOL "Whether HDF5 is linked")  # Used in directives
    include_directories(${HDF5_INCLUDE_DIRS})
    set(io_libs ${HDF5_LIBRARIES} ${io_libs})
    if(UNIX)
      # Required for hdf5 version 1.8.11 and greater
      set(io_libs ${io_libs} ${CMAKE_DL_LIBS}) 
    endif()
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
if(${ENABLE_NETCDF})
  set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_MODULE_PATH})
  find_package (NetCDF)
  if(${NetCDF_FOUND})
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
if(${ENABLE_MDSPLUS})
	find_package(MDSPLUS COMPONENTS 
               INSTALL_DIR "mdsplus"
               HEADERS "mdsdescrip.h" "mdslib.h"
	       LIBRARIES "MdsLib" "MdsShr" "MdsIpShr" "TdiShr" "TreeShr"
               INCLUDE_SUBDIRS "include"
               LIBRARY_SUBDIRS "lib"
               )
  if (${MDSPLUS_FOUND})
    message(STATUS "Found MDS+")
    set(HAVE_MDSPLUS True)
    set(USE_MDS True)                 # ifdef
    set(io_libs ${MDSPLUS_LIBRARIES} ${io_libs})
  endif()
endif()


#---------------------------------------------------------------------
#
# Setup Additional Libraries
# These are specialized and purpose built so they are taken directly
#   as input and not searched for as packages
#
#---------------------------------------------------------------------
if(NOT EXISTS ${D3_LIB})
  message(STATUS "D3_LIB not found, snap cannot be used")
  set(HAVE_SNAP FALSE)
else()
  add_library(d3lib STATIC IMPORTED)
  set_target_properties(d3lib PROPERTIES IMPORTED_LOCATION ${D3_LIB})
  set(HAVE_SNAP TRUE)
  set(USE_SNAP TRUE)  #ifdef
  if(NOT ${HAVE_MDSPLUS})
    message(STATUS "Without MDS+ no density or MSE signals available")
    set(HAVE_MSE FALSE)
  else()
    if(NOT EXISTS ${MSE_LIB})
      message(STATUS "MSE_LIB not found, some snap files will fail")
      set(HAVE_MSE FALSE)
    else()
      add_library(mselib STATIC IMPORTED)
      set_target_properties(mselib PROPERTIES IMPORTED_LOCATION ${MSE_LIB})
      set(HAVE_MSE TRUE)
      set(USE_MSE TRUE)  #ifdef
    endif()
  endif()
endif()
