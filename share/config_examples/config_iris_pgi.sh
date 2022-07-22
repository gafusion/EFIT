# In order to build and run this version of EFIT on the system you must
#   execute the following module commands:
#
#    module purge
#    module load env/gcc9.2 pgf/18.7 mpich/3.2-pgf18.7 mse
#
# If you don't want MPI (slower in serial) simply remove the FC=...
#   and -DENABLE_PARALLEL... lines

    module load cmake/3.8.2
    module load hdf5/1.8.19-mpich3.2-pgf18.7
    export CC=/fusion/usc/opt/pgi/linux86-64/18.7/bin/pgcc
    export FC=/fusion/usc/opt/mpich/mpich-3.2/gcc7.1.0-pgf18.7/bin/mpifort

    cmake \
    -DMPICMD:STRING='srun --mpi=pmi2 -n ' \
    -DBLAS_LIBRARIES:PATH='/fusion/usc/opt/pgi/linux86-64/18.7/lib/libblas.a' \
    -DLAPACK_LIBRARIES:PATH='/fusion/usc/opt/pgi/linux86-64/18.7/lib/liblapack.a' \
    -DENABLE_NETCDF:BOOL=ON \
    -DNetCDF_DIR:PATH='/fusion/usc/opt/netcdf/netcdf-4.4.1_mpich-3.2_pgf-18.3/' \
    -DENABLE_HDF5:BOOL=TRUE \
    -DD3_LIB:PATH='/fusion/projects/codes/efit/dev/d3lib_gcc7/libd3share.a' \
    -DENABLE_MDSPLUS:BOOL=ON \
    -DMDSPLUS_DIR:PATH='/fusion/usc/src/mdsplus/mdsplus-stable_release-7-96-9' \
    -DMSE_LIB:PATH='/fusion/projects/codes/mse/lib/libmse.a' \
    -DENABLE_PARALLEL:BOOL=ON \
    -DCMAKE_BUILD_TYPE:STRING=RELEASE \
    ..

