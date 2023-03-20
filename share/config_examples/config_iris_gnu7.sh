# In order to build and run this version of EFIT on the system you must
#   execute the following module commands:
#
#    module switch gcc-4.7.2 gcc7/default
#
# If you don't want MPI (slower in serial) simply remove the FC=...
#   and -DENABLE_PARALLEL... lines
#
# mselibs have not been build for GNU 7 yet

    module load cmake/3.8.2
    export CC=/fusion/projects/codes/gcc71-toolchain/bin/gcc
    export FC=/fusion/projects/codes/gcc71-toolchain/bin/mpifort

    cmake \
    -DMPICMD:STRING='srun --mpi=pmi2 -n ' \
    -DBLAS_LIBRARIES:PATH='/usr/lib64/libblas.a' \
    -DLAPACK_LIBRARIES:PATH='/usr/lib64/liblapack.a' \
    -DENABLE_NETCDF:BOOL=ON \
    -DNetCDF_DIR:PATH='/fusion/usc/opt/netcdf/netcdf-4.4.1_gcc7' \
    -DENABLE_HDF5:BOOL=ON \
    -DD3_LIB:PATH='/fusion/projects/codes/efit/dev/d3lib_gcc7/libd3share.a' \
    -DENABLE_MDSPLUS:BOOL=ON \
    -DMDSPLUS_DIR:PATH='/fusion/usc/src/mdsplus/mdsplus-stable_release-7-96-9' \
    -DENABLE_PARALLEL:BOOL=ON \
    -DCMAKE_BUILD_TYPE:STRING=RELEASE \
    ..
