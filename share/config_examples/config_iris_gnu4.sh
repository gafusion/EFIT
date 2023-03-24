# This version of EFIT is intended to be used with the default environment (old compilers)
#   and does not require any module changes.
#
# If you don't want MPI (slower in serial) simply remove the FC=...
#   and -DENABLE_PARALLEL... lines

    module load cmake/3.8.2
    export CC=/act/gcc-4.7.2/bin/gcc
    export FC=/fusion/usc/opt/mpich/mpich-3.2/gcc-4.7.2/bin/mpifort

    cmake \
    -DMPICMD:STRING='srun --mpi=pmi2 -n ' \
    -DBLAS_LIBRARIES:PATH='/usr/lib64/libblas.a' \
    -DLAPACK_LIBRARIES:PATH='/usr/lib64/liblapack.a' \
    -DENABLE_NETCDF:BOOL=ON \
    -DNetCDF_DIR:PATH='/fusion/usc/opt/netcdf/netcdf-4.4.1_mpich-3.2_gcc-4.7.2/' \
    -DENABLE_HDF5:BOOL=ON \
    -DHDF5_ROOT:PATH='/fusion/usc/opt/hdf5/hdf5-1.8.16-mpich-gcc-4.7.2/'\
    -DD3_LIB:PATH='/fusion/projects/codes/efit/dev/d3lib/libd3share.a' \
    -DMSE_LIB:PATH='/fusion/projects/codes/mse/lib/V4_05-gcc4.7.2/libmse.a' \
    -DENABLE_MDSPLUS:BOOL=ON \
    -DMDSPLUS_DIR:PATH='/fusion/usc/src/mdsplus/mdsplus-stable_release-7-96-9' \
    -DENABLE_PARALLEL:BOOL=ON \
    -DCMAKE_BUILD_TYPE:STRING=RELEASE \
    ..

