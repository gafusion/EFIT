######################################################################
#
# CMakeLists.txt for EFIT infrastructure
#
######################################################################
# Project information
project(efit-ai)
include(io/git.cmake)
set(PROJECT_VERSION ${GIT_BRANCH}-${GIT_COMMIT_HASH})
set(PROJECT_URL ${GIT_URL})

# Alternative if we would rather do things by hand...
#set(VERSION_MAJOR "0")
#set(VERSION_MINOR "1")
#set(VERSION_PATCH "0")
#set(PROJECT_VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})
#set(EFIT_VERSION "${PROJECT_VERSION}-${PROJECT_REV}")
#set(PROJECT_URL "https://gitlab.com/efit-ai/efit")

# Required version
cmake_minimum_required(VERSION 3.8.0)

# where to find the Find<PACKAGE>.cmake files
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/io)

#---------------------------------------------------------------------
# EFIT options
#---------------------------------------------------------------------
#TODO: Can we turn off?
#set(ENABLE_CERLIB FALSE CACHE BOOL "Whether to build cerlib ")

set(ENABLE_DOCS FALSE CACHE BOOL "Whether to build the docs or not.")
set(ENABLE_PRODUCTION FALSE CACHE BOOL "This is meant to be installed.")
set(TEST_EFUND True CACHE BOOL "Whether to run efund with tests.")

set(TABLE_DIR ${CMAKE_BINARY_DIR} CACHE STRING "Location of tables")
set(INPUT_DIR ${CMAKE_BINARY_DIR} CACHE STRING "Location of input files")
set(STORE_DIR ${CMAKE_BINARY_DIR} CACHE STRING "Location of sources")

set(IDEBUG "0" CACHE STRING "Debug output level (3 max)")
if (IDEBUG GREATER 0)
  set(DEBUG_LEVEL1 TRUE)
endif ()
if (IDEBUG GREATER 1)
  set(DEBUG_LEVEL2 TRUE)
endif ()
if (IDEBUG GREATER 2)
  set(DEBUG_LEVEL3 TRUE)
endif ()

set(JDEBUG "NONE" CACHE STRING "Further debug output from section")
if (JDEBUG MATCHES "MSEL")
  set(DEBUG_MSELS TRUE)
elseif (JDEBUG MATCHES "PLTS")
  set(DEBUG_PLTS TRUE) # write surface files for debugging/plotting. Serial only, not parallel
endif ()
#---------------------------------------------------------------------
# Set default values for performance issues if necessary
#---------------------------------------------------------------------
# -  See math/efitMath.cmake for math options
# -  See io/efitIO.cmake for IO options
set(ENABLE_PARALLEL FALSE CACHE BOOL "Whether to enable MPI: Default is false")
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
else ()
  set (sharedir share-ser)
  message(STATUS "Building without MPI")
endif ()

install(FILES ${CMAKE_BINARY_DIR}/config.f 
              ${CMAKE_BINARY_DIR}/config.h DESTINATION ${sharedir}
        PERMISSIONS OWNER_READ OWNER_EXECUTE OWNER_WRITE
                    GROUP_READ GROUP_EXECUTE GROUP_WRITE
                    WORLD_READ WORLD_EXECUTE
)

#---------------------------------------------------------------------
# OpenMP creates link errors
#---------------------------------------------------------------------
if (ENABLE_OpenMP)
    find_package(OpenMP)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif ()
endif ()
#---------------------------------------------------------------------
# Set default fortran format to free and preprocess all files
#---------------------------------------------------------------------
# Worry about later
#set(CMAKE_Fortran_FORMAT FREE)
#set(CMAKE_Fortran_PREPROCESS ON)

#---------------------------------------------------------------------
# Set compiler specific definitions and flags. Flags are set to 
# good values for EFIT. Valid build types are Default, Debug, MinSizeRel
# and Release. Set endian and underscore flags for all builds.
#---------------------------------------------------------------------
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "Cray")
  option(__cray "using Cray FC Compiler" on)
  set(CMAKE_Fortran_FLAGS "-s real64 ${CMAKE_Fortran_FLAGS}")
  if (NOT ${ENABLE_NATIVE_ENDIAN})
    set(CMAKE_Fortran_FLAGS "-h byteswapio ${CMAKE_Fortran_FLAGS}")
  endif ()
  set(CMAKE_Fortran_FLAGS_DEFAULT "-O2") # for reference
  set(CMAKE_C_FLAGS_DEFAULT "-O2") # for reference
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(CMAKE_C_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS_MINSIZEREL "-O3") # no small compile options available
  set(CMAKE_C_FLAGS_MINSIZEREL "-O3") # no small compile options available
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -Rb -ei0")
  set(CMAKE_C_FLAGS_DEBUG "-O0 -g")
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  option(__gfortran "using gfortran FC Compiler" on)
  set(CMAKE_Fortran_FLAGS "-fno-second-underscore -fdefault-real-8 -fdefault-double-8 ${CMAKE_Fortran_FLAGS}")
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
  set(CMAKE_Fortran_FLAGS_DEFAULT "-O0") # for reference
  set(CMAKE_C_FLAGS_DEFAULT "-O0") # for reference
  set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast")
  set(CMAKE_C_FLAGS_RELEASE "-Ofast -DNDEBUG")
  set(CMAKE_Fortran_FLAGS_MINSIZEREL "-Os")
  set(CMAKE_C_FLAGS_MINSIZEREL "-Os -DNDEBUG")
  set(CMAKE_Fortran_FLAGS_DEBUG "-Og -g -fbounds-check -finit-real=nan -fbacktrace")
  set(CMAKE_C_FLAGS_DEBUG "-Og -g -fbounds-check")
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
  option(__ifort "using Intel FC Compiler" on)
  set(CMAKE_Fortran_FLAGS "-assume no2underscores,protect_parens -autodouble ${CMAKE_Fortran_FLAGS}")
  if (NOT ${ENABLE_NATIVE_ENDIAN})
    set(CMAKE_Fortran_FLAGS "-convert big_endian ${CMAKE_Fortran_FLAGS}")
  endif ()
  set(CMAKE_Fortran_FLAGS_DEFAULT "-O2") # for reference
  set(CMAKE_C_FLAGS_DEFAULT "-O2") # for reference
  set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -ipo -fast")
  set(CMAKE_C_FLAGS_RELEASE "-Ofast -ipo -fast")
  set(CMAKE_Fortran_FLAGS_RELEASE "-Os")
  set(CMAKE_C_FLAGS_RELEASE "-Os")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -check all -traceback")
  set(CMAKE_C_FLAGS_DEBUG "-O0 -g -traceback")
elseif ((${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI") OR (${CMAKE_Fortran_COMPILER_ID} MATCHES "NVHPC"))
  option(__pgi "using PGI FC Compiler" on)
  set(CMAKE_Fortran_FLAGS "-Mextend -r8 ${CMAKE_Fortran_FLAGS}")
  if (NOT ${ENABLE_NATIVE_ENDIAN})
    set(CMAKE_Fortran_FLAGS "-byteswapio ${CMAKE_Fortran_FLAGS}")
  endif ()
  set(CMAKE_Fortran_FLAGS_DEFAULT "-O1") # from man page, contradicts online refs
  set(CMAKE_C_FLAGS_DEFAULT "-O1") # from man page, contradicts online refs
  set(CMAKE_Fortran_FLAGS_RELEASE "-O4 -fast")
  set(CMAKE_C_FLAGS_RELEASE "-O4 -fast")
  set(CMAKE_Fortran_FLAGS_MINSIZEREL "-O4 -fast") # no small compile options available
  set(CMAKE_C_FLAGS_MINSIZEREL "-O4 -fast") # no small compile options available
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -Mbounds")
  set(CMAKE_C_FLAGS_DEBUG "-g -O0 -Mbounds")
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Flang")
  option(__clang "using Clang FC Compiler" on)
  set(CMAKE_Fortran_FLAGS "-fdefault-real-8 ${CMAKE_Fortran_FLAGS}")
  #  set(CMAKE_Fortran_FLAGS "-fno-second-underscore -fdefault-real-8 -fdefault-double-8 ${CMAKE_Fortran_FLAGS}") # these flags are part of the standard but appear to be unused...
  if (NOT ${ENABLE_NATIVE_ENDIAN})
    set(CMAKE_Fortran_FLAGS "-byteswapio ${CMAKE_Fortran_FLAGS}")
  endif ()
  set(CMAKE_Fortran_FLAGS_DEFAULT "-O")
  set(CMAKE_C_FLAGS_DEFAULT "-O")
  set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast")
  set(CMAKE_C_FLAGS_RELEASE "-Ofast")
  set(CMAKE_Fortran_FLAGS_MINSIZEREL "-Oz")
  set(CMAKE_C_FLAGS_MINSIZEREL "-Oz")
  set(CMAKE_Fortran_FLAGS_DEBUG "-Og -g")
  #  set(CMAKE_Fortran_FLAGS_DEBUG "-Og -g -fbounds-check -finit-real=nan -fbacktrace") # these flags are part of the standard, but appear to be unused...
  set(CMAKE_C_FLAGS_DEBUG "-Og -g")
endif ()

# Remove fpic flags
string(REPLACE "-fpic" "" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
string(REPLACE "-fPIC" "" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")

message(STATUS "Build type is (blank is DEFAULT): ${CMAKE_BUILD_TYPE}")

# Write flags
message(STATUS "")
message(STATUS "--------- Final compiler flags ---------")
foreach (cmp C Fortran)
  foreach (bld DEFAULT RELEASE MINSIZEREL DEBUG)
    message(STATUS "CMAKE_${cmp}_FLAGS_${bld}= ${CMAKE_${cmp}_FLAGS_${bld}}")
  endforeach ()
  message(STATUS "CMAKE_${cmp}_FLAGS= ${CMAKE_${cmp}_FLAGS}")
endforeach ()
message(STATUS "")
set(CMAKE_Fortran_BUILD_TYPE_FLAGS "${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UC}}")
message(STATUS "CMAKE_Fortran_BUILD_TYPE_FLAGS= ${CMAKE_Fortran_BUILD_TYPE_FLAGS}")

#---------------------------------------------------------------------
# Need to configure set_extvars after all other flags are set.
#---------------------------------------------------------------------
set(GEN_CONF_FILE ${PROJECT_BINARY_DIR}/efit/set_extvars_conf.F90)
set(GEN_SRC_FILE ${PROJECT_BINARY_DIR}/efit/set_extvars.F90)
configure_file(${PROJECT_SOURCE_DIR}/efit/set_extvars_conf.F90.in
               ${GEN_CONF_FILE})

# Only update the src file if needed.  This reduces compilation
if (EXISTS ${GEN_SRC_FILE})
    message("File exists:  testing diff")
    execute_process( COMMAND ${CMAKE_COMMAND} -E compare_files 
                     ${GEN_CONF_FILE} ${GEN_SRC_FILE}
                     RESULT_VARIABLE compare_result)
    if (compare_result EQUAL 1)
        file(RENAME ${GEN_CONF_FILE} ${GEN_SRC_FILE})
    endif ()
else ()
    message("File does not exist")
    file(RENAME ${GEN_CONF_FILE} ${GEN_SRC_FILE})
endif ()

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
if (HAVE_MDSPLUS)
add_subdirectory(cerlib)
endif ()
add_subdirectory(green)
add_subdirectory(efit)

#---------------------------------------------------------------------
# Add testing
#---------------------------------------------------------------------
enable_testing()
add_subdirectory(test)

#---------------------------------------------------------------------
# Documentation
#---------------------------------------------------------------------
if (ENABLE_DOCS)
add_subdirectory(docs)
endif ()

#---------------------------------------------------------------------
#
# Finalize (installation, packaging)
#
#---------------------------------------------------------------------
set(COMMON_INSTALL TRUE)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Equilibrium Fitting Code")
set(CONFIG_FILES)  
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/config.h.cmake ${PROJECT_BINARY_DIR}/config.h)
