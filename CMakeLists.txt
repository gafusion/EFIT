######################################################################
#
# CMakeLists.txt for EFIT infrastructure
#
######################################################################
# Project information
project(efit)
set(VERSION_MAJOR "0")
set(VERSION_MINOR "1")
set(VERSION_PATCH "0")
set(PROJECT_VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})
set(NIMROD_VERSION "${PROJECT_VERSION}-${PROJECT_REV}")
set(PROJECT_URL "https://gitlab.com/efit-ai/efit")

# Required version
cmake_minimum_required(VERSION 3.8.0)

#---------------------------------------------------------------------
# EFIT options
#---------------------------------------------------------------------
#TODO: Can we turn off?
#set(ENABLE_CERLIB FALSE CACHE BOOL "Whether to build cerlib ")

set(ENABLE_PRODUCTION FALSE CACHE BOOL "This is meant to be installed.")

set(TABLE_DIR ${CMAKE_BINARY_DIR} CACHE STRING "Location of tables")
set(INPUT_DIR ${CMAKE_BINARY_DIR} CACHE STRING "Location of input files")
set(STORE_DIR ${CMAKE_BINARY_DIR} CACHE STRING "Location of sources")

set(GEN_CONF_FILE ${PROJECT_BINARY_DIR}/efit/set_expath_conf.F90)
set(GEN_SRC_FILE ${PROJECT_BINARY_DIR}/efit/set_expath.F90)
configure_file(${PROJECT_SOURCE_DIR}/efit/set_expath_conf.F90.in
               ${GEN_CONF_FILE})

# Only update the src file if needed.  This reduces compilation
if(EXISTS ${GEN_SRC_FILE})
    message("File exists:  testing diff")
    execute_process( COMMAND ${CMAKE_COMMAND} -E compare_files 
                     ${GEN_CONF_FILE} ${GEN_SRC_FILE}
                     RESULT_VARIABLE compare_result)
    if( compare_result EQUAL 1)
        file(RENAME ${GEN_CONF_FILE} ${GEN_SRC_FILE})
    endif()
else()
    message("File does not exist")
    file(RENAME ${GEN_CONF_FILE} ${GEN_SRC_FILE})
endif()

#---------------------------------------------------------------------
# Set default values for performance issues if necessary
#---------------------------------------------------------------------
# -  See math/efitMath.cmake for math options
# -  See io/efitIO.cmake for IO options
set(ENABLE_OpenMP FALSE CACHE BOOL "Whether to enable OpenMP: Default is false")
set(ENABLE_SHARED_LIBS FALSE CACHE BOOL "Default is to build static")
set(ENABLE_MEM_PROF OFF CACHE BOOL "Whether to enable (limited) memory profiling")
set(ENABLE_NATIVE_ENDIAN FALSE CACHE BOOL "Whether to use native or big endian")
set(ENABLE_TIME_LEVEL1 OFF CACHE BOOL "Whether to add timings at level 1")
set(ENABLE_TIME_LEVEL2 OFF CACHE BOOL "Whether to add timings at level 2")
set(TRAP_FP_EXCEPTIONS FALSE CACHE BOOL "Whether to enable trapping of floating point exceptions")

#--------------------------------------------------------------------
# Initialization
#--------------------------------------------------------------------
enable_language(Fortran)

if (NOT DEFINED CMAKE_COLOR_MAKEFILE)
  set(CMAKE_COLOR_MAKEFILE TRUE)
endif ()
if (ENABLE_SHARED_LIBS)
  set(lib_type "LIBRARIES")
else ()
  set(lib_type "STLIBS")
endif ()

option(OBJ_MEM_PROF "enabling memory profiling" ${ENABLE_MEM_PROF})
option(TIME_LEVEL1 "do level 1 timing" ${ENABLE_TIME_LEVEL1})
option(TIME_LEVEL2 "do level 2 timing" ${ENABLE_TIME_LEVEL2})

#---------------------------------------------------------------------
# Fortran include file based on compiler.  Better than local.F
#---------------------------------------------------------------------
# always pick up config.h even for subdirs
include_directories( ${CMAKE_BINARY_DIR} )

# Format config.f
add_custom_command( 
  OUTPUT config.f
  COMMAND sed;'/^\\//;d';<;config.h;>;config.f
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR} 
)
add_custom_target(config_f DEPENDS config.f)
if (${ENABLE_PARALLEL})
  set (sharedir share-par)
  set (USEMPI TRUE)                 # ifdef
  message(STATUS "Building with MPI")
else()
  set (sharedir share-ser)
  set (USEMPI FALSE)
#  message(STATUS "Building without MPI")
endif()

install(FILES ${CMAKE_BINARY_DIR}/config.f 
              ${CMAKE_BINARY_DIR}/config.h DESTINATION ${sharedir}
        PERMISSIONS OWNER_READ OWNER_EXECUTE OWNER_WRITE
                    GROUP_READ GROUP_EXECUTE GROUP_WRITE
                    WORLD_READ WORLD_EXECUTE
)

#---------------------------------------------------------------------
# Set default fortran format to free and preprocess all files
#---------------------------------------------------------------------
# Worry about later
#set(CMAKE_Fortran_FORMAT FREE)
#set(CMAKE_Fortran_PREPROCESS ON)

#---------------------------------------------------------------------
# Set compiler specific definitions and flags. Flags are set to 
# good values for NIMROD. Valid build types are Debug, RelWithDebInfo 
# and Release. For the PGI compiler a custom build type FastBuild is 
# defined to aid develpment without the fastsse flag. Set endian and 
# underscore flags for all builds.
#---------------------------------------------------------------------
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "Cray")
  option(__cray "using Cray FC Compiler" on)
  set(CMAKE_Fortran_FLAGS "-h nosecond_underscore -s real64 ${CMAKE_Fortran_FLAGS}")
  if (NOT ${ENABLE_NATIVE_ENDIAN})
    set(CMAKE_Fortran_FLAGS "-h byteswapio ${CMAKE_Fortran_FLAGS}")
  endif ()
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  option(__gfortran "using gfortran FC Compiler" on)
  set(CMAKE_Fortran_FLAGS "-fno-second-underscore -fdefault-real-8 -fdefault-double-8 -fdec-math ${CMAKE_Fortran_FLAGS}")
  #set(CMAKE_Fortran_FLAGS "-Wimplicit-interface -Wimplicit-procedure ${CMAKE_Fortran_FLAGS}")
  # Don't allow automatic reallocation (F2003 standard) as it may be inefficient
  # inside loops, flag only needed for gfortran versions > 4.6  
  set(CMAKE_Fortran_FLAGS "-fno-realloc-lhs ${CMAKE_Fortran_FLAGS}")
  # Turn argument mismatch into a warning until fixed by updating MPI interface
  if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10.0)
    set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch ${CMAKE_Fortran_FLAGS}")
  endif ()
  if (NOT ${ENABLE_NATIVE_ENDIAN})
    set(CMAKE_Fortran_FLAGS "-fconvert=big-endian ${CMAKE_Fortran_FLAGS}")
  endif ()
  # Add trap flags and bounds checking if not debug and requested 
  # also enable backtraces upon error
  if (TRAP_FP_EXCEPTIONS)
    if (${CMAKE_BUILD_TYPE_UC} MATCHES "DEBUG") 
      set(CMAKE_Fortran_FLAGS "-ffpe-trap=invalid,zero,overflow ${CMAKE_Fortran_FLAGS}")
    else ()
      set(CMAKE_Fortran_FLAGS "-g -fbounds-check -finit-real=nan -fbacktrace -ffpe-trap=invalid,zero,overflow ${CMAKE_Fortran_FLAGS}")
    endif ()
  endif ()
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(CMAKE_C_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fbounds-check -finit-real=nan -fbacktrace")
  set(CMAKE_C_FLAGS_DEBUG "-O0 -g -fbounds-check")
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
  option(__ifort "using Intel FC Compiler" on)
  set(CMAKE_Fortran_FLAGS "-assume no2underscores,protect_parens -autodouble ${CMAKE_Fortran_FLAGS}")
  if (NOT ${ENABLE_NATIVE_ENDIAN})
    set(CMAKE_Fortran_FLAGS "-convert big_endian ${CMAKE_Fortran_FLAGS}")
  endif ()
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -check all -traceback")
  set(CMAKE_C_FLAGS_DEBUG "-O0 -g -traceback")
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
  option(__pgi "using PGI FC Compiler" on)
  set(CMAKE_Fortran_FLAGS "-Mextend -r8 ${CMAKE_Fortran_FLAGS}")
  if (NOT ${ENABLE_NATIVE_ENDIAN})
    set(CMAKE_Fortran_FLAGS "-byteswapio ${CMAKE_Fortran_FLAGS}")
  endif ()
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fastsse")
  set(CMAKE_C_FLAGS_RELEASE "-O3 -fastsse")
  set(CMAKE_Fortran_FLAGS_FASTBUILD "-O3")
  set(CMAKE_C_FLAGS_FASTBUILD "-O3")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -Mbounds")
  set(CMAKE_C_FLAGS_DEBUG "-g -O0 -Mbounds")
endif ()

# Remove fpic flags
string(REPLACE "-fpic" "" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
string(REPLACE "-fPIC" "" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")

# Write flags
message(STATUS "")
message(STATUS "--------- Final compiler flags ---------")
foreach (cmp C Fortran)
  foreach (bld FULL RELEASE RELWITHDEBINFO MINSIZEREL DEBUG)
    message(STATUS "CMAKE_${cmp}_FLAGS_${bld}= ${CMAKE_${cmp}_FLAGS_${bld}}")
  endforeach ()
  message(STATUS "CMAKE_${cmp}_FLAGS= ${CMAKE_${cmp}_FLAGS}")
endforeach ()
message(STATUS "")
set(CMAKE_Fortran_BUILD_TYPE_FLAGS "${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UC}}")
message(STATUS "CMAKE_Fortran_BUILD_TYPE_FLAGS= ${CMAKE_Fortran_BUILD_TYPE_FLAGS}")

#---------------------------------------------------------------------
# Always use rpath to greatest extent.
#---------------------------------------------------------------------
if (NOT DEFINED CMAKE_INSTALL_RPATH_USE_LINK_PATH)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif ()

#---------------------------------------------------------------------
# EFIT needs blas.  
# netlib_lite is overkill for efit, but we may use some of these 
# routines in the future.  Defer until later to decide how we 
# want to include it
#---------------------------------------------------------------------
include(math/efitMath.cmake)
add_subdirectory(math/netlib_lite)

#---------------------------------------------------------------------
# External IO libraries play a big role in what EFIT interfaces to
#---------------------------------------------------------------------
include(io/efitIO.cmake)
add_subdirectory(io)

#---------------------------------------------------------------------
# Add subdirectories 
#---------------------------------------------------------------------
# These are DIII-D libraries that we are trying to avoid
#add_subdirectory(d3lib)
#add_subdirectory(mselib)

add_subdirectory(cerlib)
add_subdirectory(green)
add_subdirectory(efit)

#---------------------------------------------------------------------
# Add testing
#---------------------------------------------------------------------
enable_testing()
add_subdirectory(test)


#---------------------------------------------------------------------
#
# Finalize (installation, packaging)
#
#---------------------------------------------------------------------
set(COMMON_INSTALL TRUE)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "NonIdeal Magnetohydrodynmics with Rotation Open Discussion")
set(CONFIG_FILES)  
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/config.h.cmake ${PROJECT_BINARY_DIR}/config.h)
