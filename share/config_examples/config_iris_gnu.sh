# In order to run this version of EFIT on the system you must
#   execute the following module commands: (not required to build)
#
#    module switch gcc-4.7.2 gcc7/default
#
# If you don't want MPI (slower in serial) simply remove the FC=...
#   and -DENABLE_PARALLEL... lines

    module load cmake/3.8.2
    export CC=/fusion/projects/codes/gcc71-toolchain/bin/gcc
    export FC=/fusion/projects/codes/gcc71-toolchain/bin/mpifort

    cmake \
    -DBLAS_LIBRARIES:PATH='/usr/lib64/libblas.a' \
    -DLAPACK_LIBRARIES:PATH='/usr/lib64/liblapack.a' \
    -DENABLE_NETCDF:BOOL=ON \
    -DNetCDF_DIR:PATH='/fusion/usc/opt/netcdf/netcdf-4.4.1_gcc7' \
    -DENABLE_HDF5:BOOL=ON \
    -DENABLE_MDSPLUS:BOOL=ON \
    -DMDSPLUS_DIR:PATH='/fusion/usc/src/mdsplus/mdsplus-stable_release-7-96-9' \
    -DENABLE_PARALLEL:BOOL=ON \
    -DCMAKE_BUILD_TYPE:STRING=RELEASE \
    ..
