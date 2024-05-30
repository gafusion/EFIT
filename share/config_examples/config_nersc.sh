#!/bin/sh
# Note: before building on any NERSC machine, you must load the cray-hdf5 module:
#    module load cray-hdf5
# Note: you must have ENABLE_HDF5 set to use the system NetCDF libraries
# Note: if ENABLE_PARALLEL is on then you will have to submit a job to run the tests, e.g.
#    salloc --nodes 1 --qos interactive --time 00:05:00 --constraint haswell

#rm -rf CMake*

comp=${PE_ENV,,}
arch=x86_64
crayroot=/opt/cray/pe
case $NERSC_HOST in
  perlmutter*)
    case ${PE_ENV,,} in
      gnu)
    
        # Set architecture based on programming environment
        io_ver=12.3

        # System linear algebra (no path required, set by environment)
        SYSTEM_BLAS_SER_LIB=''
        SYSTEM_LAPACK_SER_LIB=$SYSTEM_BLAS_SER_LIB

        # System IO libs
        # HDF5 location is set by modules so locations are not required
        SYSTEM_NETCDF_SER_DIR=${crayroot}/netcdf/default/${comp}/${io_ver}
        #SYSTEM_NETCDF_PAR_DIR=${crayroot}/netcdf-hdf5parallel/default/${comp}/${io_ver}
        ;;

      nvidia)

        # Set architecture based on programming environment
        io_ver=23.3

        # System linear algebra (no path required, set by environment)
        SYSTEM_BLAS_SER_LIB=''
        SYSTEM_LAPACK_SER_LIB=$SYSTEM_BLAS_SER_LIB

        # System IO libs
        # HDF5 location is set by modules so locations are not required
        SYSTEM_NETCDF_SER_DIR=${crayroot}/netcdf/default/${comp}/${io_ver}
        #SYSTEM_NETCDF_PAR_DIR=${crayroot}/netcdf-hdf5parallel/default/${comp}/${io_ver}
        ;;

      aocc)

        # Set architecture based on programming environment
        io_ver=4.1

        # System linear algebra (no path required, set by environment)
        SYSTEM_BLAS_SER_LIB=''
        SYSTEM_LAPACK_SER_LIB=$SYSTEM_BLAS_SER_LIB

        # System IO libs
        # HDF5 location is set by modules so locations are not required
        SYSTEM_NETCDF_SER_DIR=${crayroot}/netcdf/default/${comp}/${io_ver}
        #SYSTEM_NETCDF_PAR_DIR=${crayroot}/netcdf-hdf5parallel/default/${comp}/${io_ver}
        ;;

      cray)

        # Set architecture based on programming environment
        comp=crayclang
        io_ver=17.0

        # System linear algebra (no path required, set by environment)
        SYSTEM_BLAS_SER_LIB=''
        SYSTEM_LAPACK_SER_LIB=$SYSTEM_BLAS_SER_LIB

        # System IO libs
        # HDF5 location is set by modules so locations are not required
        SYSTEM_NETCDF_SER_DIR=${crayroot}/netcdf/default/${comp}/${io_ver}
        #SYSTEM_NETCDF_PAR_DIR=${crayroot}/netcdf-hdf5parallel/default/${comp}/${io_ver}
        ;;

      *)
        echo "Unknown programming environment, trying to setup build anyway."
        ;;
    esac
    ;;

  *)
    echo "Unknown computer system, trying to setup build anyway."
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
  ..
#make
