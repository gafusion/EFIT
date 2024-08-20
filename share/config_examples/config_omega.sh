# In order to build and run this version of EFIT on the system you must
#   execute the following module commands:
#
#    module load env/* mdsplus
#
# If you don't want MPI (slower in serial) simply remove the -DENABLE_PARALLEL... and set export FC=$MPIFC lines

#rm -rf CMake*

RUN_CMD='mpirun --bind-to none -n '
USE_OPENMP=OFF
env=$(module list 2>&1 | sed -n -e 's\^.*env/\\p')
case $env in
  gcc8*)

    export CC=/usr/bin/gcc
    export FC=/usr/bin/gfortran

    # System linear algebra (requires linking OpenMP libraries)
    SYSTEM_BLAS_SER_LIB=$BLAS_DIR'/lib/libopenblas.a'
    SYSTEM_LAPACK_SER_LIB=$SYSTEM_BLAS_SER_LIB
    USE_OPENMP=ON

    # System IO libs
    SYSTEM_NETCDF_DIR=$NETCDF_DIR
    SYSTEM_NETCDF_C_DIR=''
    SYSTEM_NETCDF_F_DIR=''

    # System physics libs
    SYSTEM_DD3_LIB='/fusion/usc/c8/src/d3lib/lib/libd3.a'
    SYSTEM_MSE_LIB='/fusion/projects/codes/mse/lib/V5_02-gcc8.x/libmse.a'
    ;;

  gcc11*)

    export CC=/opt/rh/gcc-toolset-11/root/usr/bin/gcc
    export FC=/opt/rh/gcc-toolset-11/root/usr/bin/gfortran

    # System linear algebra (requires linking OpenMP libraries)
    SYSTEM_BLAS_SER_LIB=$BLAS_DIR'/lib/libopenblas.a'
    SYSTEM_LAPACK_SER_LIB=$SYSTEM_BLAS_SER_LIB
    USE_OPENMP=ON

    # System IO libs
    SYSTEM_NETCDF_DIR=''
    SYSTEM_NETCDF_C_DIR=$NETCDF_DIR
    SYSTEM_NETCDF_F_DIR=$NETCDF_F_DIR

    # System physics libs (no compatible is 
    SYSTEM_DD3_LIB='/fusion/usc/c8/src/d3lib/lib/libd3.a'
    SYSTEM_MSE_LIB=''
    ;;

  nvhpc23*)

    # System linear algebra
    SYSTEM_BLAS_SER_LIB=$BLAS_DIR'/lib/libblas.a'
    SYSTEM_LAPACK_SER_LIB=$BLAS_DIR'/lib/liblapack.a'

    # System IO libs
    SYSTEM_NETCDF_DIR=''
    SYSTEM_NETCDF_C_DIR=$NETCDF_DIR
    SYSTEM_NETCDF_F_DIR=$NETCDF_F_DIR

    # System physics libs
    SYSTEM_DD3_LIB='/fusion/projects/codes/efit/c8/dev/d3lib_nvhpc23.11/libd3share.a'
    SYSTEM_MSE_LIB='/fusion/projects/codes/mse/lib/libmse.a'
    ;;

  intel2020*)

    export CC=/fusion/usc/opt/intel2020/bin/icc
    export FC=/fusion/usc/opt/intel2020/bin/ifort

    # System linear algebra
    SYSTEM_BLAS_SER_LIB='-L'$MKLROOT' -lmkl_core -lintlc -lmkl_intel_ilp64 -lmkl_sequential -lmkl_blas95_ilp64'
    SYSTEM_LAPACK_SER_LIB='-lmkl_lapack95_ilp64'

    # System IO libs
    SYSTEM_NETCDF_DIR=$NETCDF_DIR
    SYSTEM_NETCDF_C_DIR=''
    SYSTEM_NETCDF_F_DIR=''

    # System physics libs
    # ptdata calls hang when snap files are used (unclear why)
    #SYSTEM_DD3_LIB='/fusion/usc/c8/src/d3lib/lib/libd3.a'
    #SYSTEM_DD3_LIB='/fusion/projects/codes/efit/dev/d3lib_intel2020/libd3share.a'
    #SYSTEM_MSE_LIB='/fusion/projects/codes/mse/lib/V5_02-gcc8.x/libmse.a'
    SYSTEM_DD3_LIB=''
    SYSTEM_MSE_LIB=''
    ;;

  *)
    echo "Unknown programming environment, trying to setup build anyway."
    ;;
esac

export FC=$MPIFC

cmake \
  -DMPICMD:STRING="$RUN_CMD" \
  -DBLAS_LIBRARIES:PATH="$SYSTEM_BLAS_SER_LIB" \
  -DLAPACK_LIBRARIES:PATH="$SYSTEM_LAPACK_SER_LIB" \
  -DENABLE_OpenMP:BOOL=$USE_OPENMP \
  -DENABLE_NETCDF:BOOL=ON \
  -DNetCDF_DIR:PATH=$SYSTEM_NETCDF_DIR \
  -DNetCDF_C_DIR:PATH=$SYSTEM_NETCDF_C_DIR \
  -DNetCDF_FORTRAN_DIR:PATH=$SYSTEM_NETCDF_F_DIR\
  -DENABLE_HDF5:BOOL=ON \
  -DHDF5_ROOT:PATH=$HDF5_DIR\
  -DD3_LIB:PATH=$SYSTEM_DD3_LIB \
  -DMSE_LIB:PATH=$SYSTEM_MSE_LIB \
  -DENABLE_MDSPLUS:BOOL=ON \
  -DENABLE_PARALLEL:BOOL=ON \
  -DCMAKE_BUILD_TYPE:STRING=RELEASE \
  ..
#make
