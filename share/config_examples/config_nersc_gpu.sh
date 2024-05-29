#!/bin/sh
# Note: before building on any NERSC machine, you must load the cray-hdf5 module:
#    module load cray-hdf5
# Note: you must have ENABLE_HDF5 set to use the system NetCDF libraries
# Note: if ENABLE_PARALLEL is on then you will have to submit a job to run the tests, e.g.
#    salloc --nodes 1 --qos interactive --time 00:05:00 --constraint haswell

#rm -rf CMake*

if ! [ ${PE_ENV,,} = 'nvidia' ]; then
  echo Wrong programming environment. Should be PrgEnv-nvidia.
  exit
fi

arch=unknown
case $NERSC_HOST in
  perlmutter*)

    # Determine architecture from programming environment
    comp=nvidia
    io_ver=23.3
    arch=x86_64
    crayroot=/opt/cray/pe

    # System linear algebra (no path required, set by environment)
    SYSTEM_BLAS_SER_LIB=''
    SYSTEM_LAPACK_SER_LIB=$SYSTEM_BLAS_SER_LIB

    # System IO libs
    # HDF5 location is set by modules so locations are not required
    SYSTEM_NETCDF_SER_DIR=${crayroot}/netcdf/default/${comp}/${io_ver}
    #SYSTEM_NETCDF_PAR_DIR=${crayroot}/netcdf-hdf5parallel/default/${comp}/${io_ver}
    ;;
esac

cmake \
  -DMPICMD:STRING='srun -n ' \
  -DBLAS_LIBRARIES:PATH=$SYSTEM_BLAS_SER_LIB \
  -DLAPACK_LIBRARIES:PATH=$SYSTEM_LAPACK_SER_LIB \
  -DENABLE_NETCDF:BOOL=ON \
  -DNetCDF_DIR=$SYSTEM_NETCDF_SER_DIR \
  -DENABLE_HDF5:BOOL=ON \
  -DENABLE_PARALLEL:BOOL=ON \
  -DCMAKE_BUILD_TYPE:STRING=Release \
  -DCMAKE_Fortran_FLAGS:STRING='-gpu=cc80,managed' \
  -DENABLE_OPENMP_NV:BOOL=ON \
  ..
