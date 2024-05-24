#!/bin/sh
# Note: before building on an ALCF machine, you must load the cray-hdf5 module:
#    module load cray-hdf5
# Note: if ENABLE_PARALLEL is on then you will have to submit a job to run the tests, e.g.
#    qsub -I -l select=1 -l filesystems=home:eagle -l walltime=00:05:00 -q debug
# Note: you must have ENABLE_HDF5 set to use the system NetCDF libraries

#rm -rf CMake*

if ! module list 2>&1 | grep PrgEnv-gnu 1>/dev/null 2>&1; then
  echo Wrong programming environment. Should be PrgEnv-gnu.
  exit
fi

arch=unknown
case ${HOST:0:7} in
  polaris*)

# Determine architecture from programming environment
    comp=gnu
    math_ver=12.3
    io_ver=12.3
    arch=x86_64
    crayroot=/opt/cray/pe

# System linear algebra
    SYSTEM_BLAS_SER_LIB=${crayroot}/libsci/default/${comp}/${math_ver}/${arch}/lib/libsci_gnu.so
    SYSTEM_LAPACK_SER_LIB=$SYSTEM_BLAS_SER_LIB

# System IO libs
#    # HDF5 location is set by modules so locations are not required
#    SYSTEM_HDF5_SER_DIR=${crayroot}/hdf5/default/${comp}/${io_ver}
#    SYSTEM_HDF5_PAR_DIR=${crayroot}/hdf5-parallel/default/${comp}/${io_ver}
    SYSTEM_NETCDF_SER_DIR=${crayroot}/netcdf/default/${comp}/${io_ver}
#    SYSTEM_NETCDF_PAR_DIR=${crayroot}/netcdf-hdf5parallel/default/${comp}/${io_ver}
    ;;
esac

cmake \
  -DMPICMD:STRING='mpiexec -np ' \
  -DBLAS_LIBRARIES:PATH=$SYSTEM_BLAS_SER_LIB \
  -DLAPACK_LIBRARIES:PATH=$SYSTEM_LAPACK_SER_LIB \
  -DENABLE_NETCDF:BOOL=ON \
  -DNetCDF_DIR=$SYSTEM_NETCDF_SER_DIR \
  -DNetCDF_INCLUDE_DIR:STRING=$SYSTEM_NETCDF_SER_DIR'/include' \
  -DENABLE_HDF5:BOOL=ON \
  -DENABLE_PARALLEL:BOOL=ON \
  -DCMAKE_BUILD_TYPE:STRING=Release \
  ..
