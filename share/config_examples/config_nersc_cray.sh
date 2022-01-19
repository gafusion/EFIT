#!/bin/sh

# Note: if ENABLE_PARALLEL is on then you will have to submit a job to run the tests, e.g.
#    salloc --nodes 1 --qos interactive --time 00:05:00 --constraint haswell
# Note: you must have ENABLE_HDF5 set to use the system NetCDF libraries

#rm -rf CMake*

if ! modulecmd bash list 2>&1 | grep PrgEnv-cray 1>/dev/null 2>&1; then
  echo Wrong programming environment. Should be PrgEnv-cray.
  exit
fi

arch=unknown
case $NERSC_HOST in
  cori*)

# Determine architecture from programming environment
    comp=crayclang
    ver=9.0
    arch=haswell
    if modulecmd bash list 2>&1 | grep craype-mic-knl 1>/dev/null 2>&1; then
      arch=mic_knl
    fi
    crayroot=/opt/cray/pe

# System linear algebra
    SYSTEM_BLAS_SER_LIB=${crayroot}/libsci/default/${comp}/${ver}/${arch}/lib/libsci_cray.a
    SYSTEM_LAPACK_SER_LIB=$SYSTEM_BLAS_SER_LIB
#    MKL_ROOT_DIR=/opt/intel/compilers_and_libraries_2020/linux/mkl

# System IO libs
#    # System hdf5 ser library broken on head nodes 8/19
    SYSTEM_HDF5_SER_DIR=${crayroot}/hdf5/default/${comp}/${ver}
#    SYSTEM_HDF5_PAR_DIR=${crayroot}/hdf5-parallel/default/${comp}/${ver}
    SYSTEM_NETCDF_SER_DIR=${crayroot}/netcdf/default/${comp}/${ver}
#    SYSTEM_NETCDF_PAR_DIR=${crayroot}/netcdf-hdf5parallel/default/${comp}/${ver}
    ;;
esac

export PATH=$SYSTEM_HDF5_SER_DIR:${PATH}

cmake \
-DBLAS_LIBRARIES:PATH=$SYSTEM_BLAS_SER_LIB \
-DLAPACK_LIBRARIES:PATH=$SYSTEM_LAPACK_SER_LIB \
-DENABLE_NETCDF:BOOL=ON \
-DNetCDF_DIR=$SYSTEM_NETCDF_SER_DIR \
-DENABLE_HDF5:BOOL=ON \
-DENABLE_PARALLEL:BOOL=ON \
-DCMAKE_BUILD_TYPE:STRING=Release \
..
