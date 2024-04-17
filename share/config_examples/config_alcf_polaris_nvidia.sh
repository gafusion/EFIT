#!/bin/sh

# should do module load PrgEnv-nvhpc cmake cray-hdf5
rm -rf CMake*

# System linear algebra (no path required, set by environment)
    SYSTEM_BLAS_SER_LIB=''
    SYSTEM_LAPACK_SER_LIB=$SYSTEM_BLAS_SER_LIB


cmake \
  -DMPICMD:STRING='mpiexec -np ' \
  -DBLAS_LIBRARIES:PATH=$SYSTEM_BLAS_SER_LIB \
  -DLAPACK_LIBRARIES:PATH=$SYSTEM_LAPACK_SER_LIB \
  -DENABLE_NETCDF:BOOL=ON \
  -DNetCDF_DIR:STRING='/opt/cray/pe/netcdf/default/nvidia/20.7' \
  -DNetCDF_INCLUDE_DIR:STRING='/opt/cray/pe/netcdf/default/nvidia/20.7/include' \
  -DENABLE_HDF5:BOOL=ON \
  -DENABLE_NETCDF:BOOL=ON \
  -DENABLE_PARALLEL:BOOL=ON \
  -DCMAKE_BUILD_TYPE:STRING=Release \
  ..
