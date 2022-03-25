# Try to find the MDS+ library and headers
# Usage of this module is as follows
#
# ===========================================================================
# Variables used by this module which can be used to change the default
# behaviour, and hence need to be set before calling find_package:
#
#  MDSPLUS_DIR
#       The preferred installation prefix for searching for MDSPLUS
#       Set this if the module has problems finding the proper MDSPLUS installation.
#
# If you don't supply MDSPLUS_DIR, the module will search on the standard
# system paths.
#
# ============================================================================
# Variables set by this module:
#
#  MDSPLUS_FOUND           System has MDSPLUS.
#
#  MDSPLUS_INCLUDE_DIRS    MDSPLUS include directories: not cached.
#
#  MDSPLUS_LIBRARIES       Link to these to use the MDSPLUS library: not cached.
#
# ===========================================================================
# If MDSPLUS is installed in a non-standard way, e.g. a non GNU-style install
# of <prefix>/{lib,include}, then this module may fail to locate the headers
# and libraries as needed. In this case, the following cached variables can
# be editted to point to the correct locations.
#
#  MDSPLUS_INCLUDE_DIR    The path to the MDSPLUS include directory: cached
#
#  MDSPLUS_LIBRARY        The path to the MDSPLUS library: cached
#
# You should not need to set these in the vast majority of cases
#

#-----------------------------------------------------------------------------

find_path(MDSPLUS_DIR
    NAMES include/mdsplus.h mdsplus/include/mdsplus.h include/mdsplus/mdsplus.h
    HINTS ENV MDSPLUS_DIR
    PATH_SUFFIXES lib lib64
    DOC "MDSPLUS root installation directory")

#-----------------------------------------------------------------------------

find_path(MDSPLUS_INCLUDE_DIR
    NAMES mdsplus.h mdsdescrip.h
    HINTS ${MDSPLUS_DIR} ENV MDSPLUS_DIR ENV CPATH
    PATH_SUFFIXES include include/mdsplus mdsplus/include
    DOC "Path to the MDSPLUS headers")

#-----------------------------------------------------------------------------

find_library(MDSPLUS_LIB
    NAMES MdsLib
    HINTS ${MDSPLUS_DIR} ENV MDSPLUS_DIR ENV LD_LIBRARY_PATH
        ENV LIBRARY_PATH ENV DYLD_LIBRARY_PATH
    PATH_SUFFIXES
        lib lib64 mdsplus lib/mdsplus lib64/mdsplus
        # system processor
        lib/mdsplus/${CMAKE_SYSTEM_PROCESSOR}
        lib64/mdsplus/${CMAKE_SYSTEM_PROCESSOR}
        lib/mdsplus/${CMAKE_SYSTEM_PROCESSOR}/lib
        lib64/mdsplus/${CMAKE_SYSTEM_PROCESSOR}/lib64
        lib64/mdsplus/${CMAKE_SYSTEM_PROCESSOR}/lib
        ${CMAKE_SYSTEM_PROCESSOR}/lib
        ${CMAKE_SYSTEM_PROCESSOR}/lib/mdsplus
        ${CMAKE_SYSTEM_PROCESSOR}/lib64
        ${CMAKE_SYSTEM_PROCESSOR}/lib64/mdsplus
    DOC "Path to the MDSPLUS library")

find_library(MdsShr_LIB
    NAMES MdsShr
    HINTS ${MDSPLUS_DIR} ENV MDSPLUS_DIR ENV LD_LIBRARY_PATH
        ENV LIBRARY_PATH ENV DYLD_LIBRARY_PATH
    PATH_SUFFIXES
        lib lib64 mdsplus lib/mdsplus lib64/mdsplus
        # system processor
        lib/mdsplus/${CMAKE_SYSTEM_PROCESSOR}
        lib64/mdsplus/${CMAKE_SYSTEM_PROCESSOR}
        lib/mdsplus/${CMAKE_SYSTEM_PROCESSOR}/lib
        lib64/mdsplus/${CMAKE_SYSTEM_PROCESSOR}/lib64
        lib64/mdsplus/${CMAKE_SYSTEM_PROCESSOR}/lib
        ${CMAKE_SYSTEM_PROCESSOR}/lib
        ${CMAKE_SYSTEM_PROCESSOR}/lib/mdsplus
        ${CMAKE_SYSTEM_PROCESSOR}/lib64
        ${CMAKE_SYSTEM_PROCESSOR}/lib64/mdsplus
    DOC "Path to the MDSPLUS library")

find_library(MdsIpShr_LIB
    NAMES MdsIpShr
    HINTS ${MDSPLUS_DIR} ENV MDSPLUS_DIR ENV LD_LIBRARY_PATH
        ENV LIBRARY_PATH ENV DYLD_LIBRARY_PATH
    PATH_SUFFIXES
        lib lib64 mdsplus lib/mdsplus lib64/mdsplus
        # system processor
        lib/mdsplus/${CMAKE_SYSTEM_PROCESSOR}
        lib64/mdsplus/${CMAKE_SYSTEM_PROCESSOR}
        lib/mdsplus/${CMAKE_SYSTEM_PROCESSOR}/lib
        lib64/mdsplus/${CMAKE_SYSTEM_PROCESSOR}/lib64
        lib64/mdsplus/${CMAKE_SYSTEM_PROCESSOR}/lib
        ${CMAKE_SYSTEM_PROCESSOR}/lib
        ${CMAKE_SYSTEM_PROCESSOR}/lib/mdsplus
        ${CMAKE_SYSTEM_PROCESSOR}/lib64
        ${CMAKE_SYSTEM_PROCESSOR}/lib64/mdsplus
    DOC "Path to the MDSPLUS library")

find_library(TdiShr_LIB
    NAMES TdiShr
    HINTS ${MDSPLUS_DIR} ENV MDSPLUS_DIR ENV LD_LIBRARY_PATH
        ENV LIBRARY_PATH ENV DYLD_LIBRARY_PATH
    PATH_SUFFIXES
        lib lib64 mdsplus lib/mdsplus lib64/mdsplus
        # system processor
        lib/mdsplus/${CMAKE_SYSTEM_PROCESSOR}
        lib64/mdsplus/${CMAKE_SYSTEM_PROCESSOR}
        lib/mdsplus/${CMAKE_SYSTEM_PROCESSOR}/lib
        lib64/mdsplus/${CMAKE_SYSTEM_PROCESSOR}/lib64
        lib64/mdsplus/${CMAKE_SYSTEM_PROCESSOR}/lib
        ${CMAKE_SYSTEM_PROCESSOR}/lib
        ${CMAKE_SYSTEM_PROCESSOR}/lib/mdsplus
        ${CMAKE_SYSTEM_PROCESSOR}/lib64
        ${CMAKE_SYSTEM_PROCESSOR}/lib64/mdsplus
    DOC "Path to the MDSPLUS library")

find_library(TreeShr_LIB
    NAMES TreeShr
    HINTS ${MDSPLUS_DIR} ENV MDSPLUS_DIR ENV LD_LIBRARY_PATH
        ENV LIBRARY_PATH ENV DYLD_LIBRARY_PATH
    PATH_SUFFIXES
        lib lib64 mdsplus lib/mdsplus lib64/mdsplus
        # system processor
        lib/mdsplus/${CMAKE_SYSTEM_PROCESSOR}
        lib64/mdsplus/${CMAKE_SYSTEM_PROCESSOR}
        lib/mdsplus/${CMAKE_SYSTEM_PROCESSOR}/lib
        lib64/mdsplus/${CMAKE_SYSTEM_PROCESSOR}/lib64
        lib64/mdsplus/${CMAKE_SYSTEM_PROCESSOR}/lib
        ${CMAKE_SYSTEM_PROCESSOR}/lib
        ${CMAKE_SYSTEM_PROCESSOR}/lib/mdsplus
        ${CMAKE_SYSTEM_PROCESSOR}/lib64
        ${CMAKE_SYSTEM_PROCESSOR}/lib64/mdsplus
    DOC "Path to the MDSPLUS library")

#-----------------------------------------------------------------------------

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set MDSPLUS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(MDSPLUS DEFAULT_MSG MDSPLUS_INCLUDE_DIR
	MDSPLUS_LIB MdsShr_LIB MdsIpShr_LIB TdiShr_LIB TreeShr_LIB)

#-----------------------------------------------------------------------------

if(MDSPLUS_FOUND)
    add_library(MDSPLUS INTERFACE)
    target_link_libraries(MDSPLUS INTERFACE ${MDSPLUS_LIBRARY})
    target_include_directories(MDSPLUS INTERFACE ${MDSPLUS_INCLUDE_DIR})
    get_filename_component(MDSPLUS_INCLUDE_DIRS ${MDSPLUS_INCLUDE_DIR} REALPATH)
    get_filename_component(MDSPLUS_LIBRARY ${MDSPLUS_LIB} REALPATH)
    get_filename_component(MdsShr_LIBRARY ${MdsShr_LIB} REALPATH)
    get_filename_component(MdsIpShr_LIBRARY ${MdsIpShr_LIB} REALPATH)
    get_filename_component(TdiShr_LIBRARY ${TdiShr_LIB} REALPATH)
    get_filename_component(TreeShr_LIBRARY ${TreeShr_LIB} REALPATH)
    set(MDSPLUS_LIBRARIES ${MDSPLUS_LIBRARY} ${MdsShr_LIBRARY} ${MdsIpShr_LIBRARY}
	     ${TdiShr_LIBRARY} ${TreeShr_LIBRARY})
endif()

mark_as_advanced(MDSPLUS_INCLUDE_DIR MDSPLUS_LIBRARY)
