# In order to build and run this version of EFIT on the system you must
#   execute the following module commands:
#
#    module purge
#    module load env/gcc9.2
#
# If you don't want MPI (slower in serial) simply remove the FC=...
#   and -DENABLE_PARALLEL... lines

    module load cmake/3.8.2
    export CC=/fusion/usc/opt/gcc/gcc-9.2.0/bin/gcc
    export FC=/fusion/usc/opt/mpich/mpich-3.2/gcc-9.2.0/bin/mpifort

    cmake \
    -DMPICMD:STRING='srun --mpi=pmi2 -n ' \
    -DBLAS_LIBRARIES:PATH='/usr/lib64/libblas.a' \
    -DLAPACK_LIBRARIES:PATH='/usr/lib64/liblapack.a' \
    -DENABLE_NETCDF:BOOL=ON \
    -DNetCDF_DIR:PATH='/fusion/usc/opt/env/gcc9.2' \
    -DENABLE_HDF5:BOOL=ON \
    -DHDF5_ROOT:PATH='/fusion/usc/opt/env/gcc9.2'\
    -DD3_LIB:PATH='/fusion/projects/codes/efit/dev/d3lib_gcc9.2.0/libd3share.a' \
    -DENABLE_MDSPLUS:BOOL=ON \
    -DMDSPLUS_DIR:PATH='/fusion/usc/src/mdsplus/mdsplus-stable_release-7-96-9' \
    -DMSE_LIB:PATH='/fusion/projects/codes/mse/lib/V4_05-gcc9.2/libmse.a' \
    -DENABLE_PARALLEL:BOOL=ON \
    -DCMAKE_BUILD_TYPE:STRING=RELEASE \
    ..
