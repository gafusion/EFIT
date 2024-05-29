#!/bin/sh
# Note: before building on an ALCF machine, you should load the cray-hdf5 module:
#     module load cray-hdf5
# Note: you must have ENABLE_HDF5 set to use the system NetCDF libraries
# Note: if ENABLE_PARALLEL is on then you will have to build from a compute node to run the tests, e.g.
#     qsub -I -l select=1 -l filesystems=home:eagle -l walltime=00:15:00 -q debug
#   with the following modules loaded before starting the build
#     module use /soft/modulefiles
#     module load spack-pe-base cmake

#rm -rf CMake*

comp=${PE_ENV,,}
arch=x86_64
crayroot=/opt/cray/pe
system=$(hostname -f 2>&1 | sed -n -e 's/^.*hsn.cm.//p')
env=$(module list 2>&1 | sed -n -e 's/^.*PrgEnv-//p')
case $system in
  polaris*)
    case ${PE_ENV,,} in
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
      
      gnu)
      
        # Set architecture based on programming environment
        math_ver=12.3
        io_ver=12.3
      
        # System linear algebra
        SYSTEM_BLAS_SER_LIB=${crayroot}/libsci/default/${comp}/${math_ver}/${arch}/lib/libsci_gnu.so
        SYSTEM_LAPACK_SER_LIB=$SYSTEM_BLAS_SER_LIB
      
        # System IO libs
        # HDF5 location is set by modules so locations are not required
        SYSTEM_NETCDF_SER_DIR=${crayroot}/netcdf/default/${comp}/${io_ver}
        #SYSTEM_NETCDF_PAR_DIR=${crayroot}/netcdf-hdf5parallel/default/${comp}/${io_ver}
        ;;
      
      cray)
      
        # Set architecture based on programming environment
        comp=crayclang
        math_ver=17.0
        io_ver=17.0
      
        # System linear algebra
        SYSTEM_BLAS_SER_LIB=${crayroot}/libsci/default/${comp}/${math_ver}/${arch}/lib/libsci_cray.so
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
  -DMPICMD:STRING='mpiexec -np ' \
  -DBLAS_LIBRARIES:PATH=$SYSTEM_BLAS_SER_LIB \
  -DLAPACK_LIBRARIES:PATH=$SYSTEM_LAPACK_SER_LIB \
  -DENABLE_NETCDF:BOOL=ON \
  -DNetCDF_DIR:STRING=$SYSTEM_NETCDF_SER_DIR \
  -DNetCDF_INCLUDE_DIR:STRING=$SYSTEM_NETCDF_SER_DIR'/include' \
  -DENABLE_HDF5:BOOL=ON \
  -DENABLE_PARALLEL:BOOL=ON \
  -DCMAKE_BUILD_TYPE:STRING=Release \
  ..
#make
