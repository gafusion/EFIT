# WARNING: there is currently no system HDF5 install that works
#   with this compiler version
# In order to run this version of EFIT on the system you must
#   execute the following module commands: (not required to build)
#
#    module switch gcc-4.7.2 gcc-9.2.0
#    module load mpich/3.2-gcc9.2.0
#
# If you don't want MPI (slower in serial) simply remove the FC=...
#   and -DENABLE_PARALLEL... lines

    module load cmake/3.8.2
    export CC=/fusion/usc/opt/gcc/gcc-9.2.0/bin/gcc
    export FC=/fusion/usc/opt/mpich/mpich-3.2/gcc-9.2.0/bin/mpifort

    cmake \
    -DBLAS_LIBRARIES:PATH='/usr/lib64/libblas.a' \
    -DLAPACK_LIBRARIES:PATH='/usr/lib64/liblapack.a' \
    -DENABLE_PARALLEL:BOOL=ON \
    -DENABLE_NETCDF:BOOL=ON \
    -DNetCDF_DIR:PATH='/fusion/usc/opt/env/gcc9.2' \
    -DCMAKE_BUILD_TYPE:STRING=RELEASE \
    ..
