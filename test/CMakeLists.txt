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
if (${ENABLE_PARALLEL})
    if (NOT DEFINED MPIEXEC)
       # mpiexec most common on laptops but we run
       # parallel at iris and nersc which are srun
       set(MPIEXEC "srun --mpi=pmi2")
    endif()
    set(NPROC "-n 4")
else ()
    set(MPIEXEC "")
    set(NPROC "") 
endif ()
set(EFIT_EXEC ${CMAKE_BINARY_DIR}/efit/efit)
set(EFIT_MODE 2)
set(EFIT_TIME 1)
set(EFIT_TEST_CONFIGFILE ${CMAKE_CURRENT_SOURCE_DIR}/run_efit.sh.in)
set(EFIT_TEST_CONFIGFILE_HDF5 ${CMAKE_CURRENT_SOURCE_DIR}/run_efit_hdf5.sh.in)
set(EFIT_RESULTS_CONFIGFILE ${CMAKE_CURRENT_SOURCE_DIR}/check_efit_results.sh.in)
set(RFILE_TEST False)

set(IS_DIIID True)
add_subdirectory(DIIID)
set(IS_DIIID False)
add_subdirectory(NSTX)
add_subdirectory(ITER)
