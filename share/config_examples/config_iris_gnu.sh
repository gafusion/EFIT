# WARNING: there is currently no system HDF5 install that works
#   with this compiler version
# In order to run this version of EFIT on the system you must
#   execute the following module commands: (not required to build)
#
#    module switch gcc-4.7.2 gcc-9.2.0
#
# If you don't want MPI (slower in serial) simply remove the FC=...
#   and -DENABLE_PARALLEL... lines
#
# Iris does not appear to have MDS+ support for GNU at this time

    module load cmake/3.8.2
    export CC=/fusion/usc/opt/gcc/gcc-9.2.0/bin/gcc
    export FC=/fusion/usc/opt/mpich/mpich-3.2/gcc-9.2.0/bin/mpifort

    cmake \
    -DMPICMD:STRING='srun --mpi=pmi2' \
    -DBLAS_LIBRARIES:PATH='/usr/lib64/libblas.a' \
    -DLAPACK_LIBRARIES:PATH='/usr/lib64/liblapack.a' \
    -DENABLE_NETCDF:BOOL=ON \
    -DNetCDF_DIR:PATH='/fusion/usc/opt/env/gcc9.2' \
    -DENABLE_HDF5:BOOL=ON \
    -DHDF5_ROOT:PATH='/fusion/usc/opt/env/gcc9.2'\
    -DENABLE_PARALLEL:BOOL=ON \
    -DCMAKE_BUILD_TYPE:STRING=RELEASE \
    ..
