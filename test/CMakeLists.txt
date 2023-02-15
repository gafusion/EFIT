######################################################################
#
# CMakeLists.txt for nimtests
#
######################################################################

#---------------------------------------------------------------------
# Initialization
#---------------------------------------------------------------------
include_directories(
  ${CMAKE_BINARY_DIR}/green
  ${CMAKE_BINARY_DIR}/efit
)

set(GRID_SIZE 65)

# This is a shell script used to drive running efund
set(EFUND_EXEC ${CMAKE_BINARY_DIR}/green/efund)
set(EFUND_TEST_CONFIGFILE ${CMAKE_CURRENT_SOURCE_DIR}/run_green.sh.in)
set(EFUND_RESULTS_CONFIGFILE ${CMAKE_CURRENT_SOURCE_DIR}/check_green_results.sh.in)

# This is a shell script used to drive running efit
set(MPICMD "mpirun -n " CACHE STRING "Command to invoke ")
set(NPROC "4" CACHE STRING "Processors to use ")
if (NOT ${ENABLE_PARALLEL})
  set(MPICMD "")
  set(NPROC "") 
endif ()
set(EFIT_EXEC ${CMAKE_BINARY_DIR}/efit/efit)
set(EFIT_MODE 2)
set(EFIT_COUNT 1)
set(EFIT_SETUP_CONFIGFILE ${CMAKE_CURRENT_SOURCE_DIR}/efit.setup.in)
set(EFIT_INPUT_CONFIGFILE ${CMAKE_CURRENT_SOURCE_DIR}/efit.input.in)
set(EFIT_TEST_CONFIGFILE ${CMAKE_CURRENT_SOURCE_DIR}/run_efit.sh.in)
set(EFIT_TEST_CONFIGFILE_INPUT ${CMAKE_CURRENT_SOURCE_DIR}/run_efit_input.sh.in)
set(EFIT_TEST_CONFIGFILE_HDF5 ${CMAKE_CURRENT_SOURCE_DIR}/run_efit_hdf5.sh.in)
set(EFIT_TEST_CONFIGFILE_SNAP ${CMAKE_CURRENT_SOURCE_DIR}/run_efit_snap.sh.in)
set(EFIT_RESULTS_CONFIGFILE ${CMAKE_CURRENT_SOURCE_DIR}/check_efit_results.sh.in)
set(EFIT_RESULTS_CONFIGFILE_SNAP ${CMAKE_CURRENT_SOURCE_DIR}/check_efit_snap_results.sh.in)
set(RFILE_TEST False)
if (${HAVE_NETCDF})
  set(outputs "g a m")
else ()
  set(outputs "g a")
endif ()

set(IS_DIIID True)
set(IS_NSTX False)
add_subdirectory(DIIID)
set(IS_DIIID False)
set(IS_NSTX True)
add_subdirectory(NSTX)
set(IS_NSTX False)
add_subdirectory(ITER)
