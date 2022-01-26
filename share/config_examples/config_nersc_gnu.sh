#!/bin/sh

# Note: if ENABLE_PARALLEL is on then you will have to submit a job to run the tests, e.g.
#    salloc --nodes 1 --qos interactive --time 00:05:00 --constraint haswell
# Note: you must have ENABLE_HDF5 set to use the system NetCDF libraries

#rm -rf CMake*

if ! modulecmd bash list 2>&1 | grep PrgEnv-gnu 1>/dev/null 2>&1; then
  echo Wrong programming environment. Should be PrgEnv-gnu.
  exit
fi

arch=unknown
case $NERSC_HOST in
  cori*)

# Determine architecture from programming environment
    comp=gnu
    math_ver=8.1
    io_ver=8.2 # unclear why these are different
    arch=haswell
    if modulecmd bash list 2>&1 | grep craype-mic-knl 1>/dev/null 2>&1; then
      arch=mic_knl
    fi
    crayroot=/opt/cray/pe

# System linear algebra
    SYSTEM_BLAS_SER_LIB=${crayroot}/libsci/default/${comp}/${math_ver}/${arch}/lib/libsci_gnu.a
    SYSTEM_LAPACK_SER_LIB=$SYSTEM_BLAS_SER_LIB
#    MKL_ROOT_DIR=/opt/intel/compilers_and_libraries_2020/linux/mkl

# System IO libs
#    # System hdf5 ser library broken on head nodes 8/19
    SYSTEM_HDF5_SER_DIR=${crayroot}/hdf5/default/${comp}/${io_ver}
#    SYSTEM_HDF5_PAR_DIR=${crayroot}/hdf5-parallel/default/${comp}/${io_ver}
    SYSTEM_NETCDF_SER_DIR=${crayroot}/netcdf/default/${comp}/${io_ver}
#    SYSTEM_NETCDF_PAR_DIR=${crayroot}/netcdf-hdf5parallel/default/${comp}/${io_ver}
    ;;
esac

export PATH=$SYSTEM_HDF5_SER_DIR:${PATH}

cmake \
  -DMPICMD:STRING='srun' \
  -DBLAS_LIBRARIES:PATH=$SYSTEM_BLAS_SER_LIB \
  -DLAPACK_LIBRARIES:PATH=$SYSTEM_LAPACK_SER_LIB \
  -DENABLE_PARALLEL:BOOL=ON \
  -DENABLE_NETCDF:BOOL=ON \
  -DNetCDF_DIR=$SYSTEM_NETCDF_SER_DIR \
  -DENABLE_HDF5:BOOL=ON \
  -DCMAKE_BUILD_TYPE:STRING=Release \
  ..
