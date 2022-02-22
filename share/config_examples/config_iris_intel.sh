# In order to run this version of EFIT on the system you must
#   execute the following module commands: (not required to build)
#
#    module switch gcc-4.7.2 gcc-9.2.0
#    module load intel/2018 mpich/3.2-intel2018
#
# If you don't want MPI (slower in serial) simply remove the FC=...
#   and -DENABLE_PARALLEL... lines
#
# Iris does not appear to have MDS+ support for Intel at this time

    module load cmake/3.8.2
    module load hdf5/1.8.19-mpich3.2-intel2018
    export CC=/fusion/usc/opt/intel2018/bin/icc
    export FC=/fusion/usc/opt/mpich/mpich-3.2/intel2018/bin/mpifort

    cmake \
    -DMPICMD:STRING='srun --mpi=pmi2' \
    -DBLAS_LIBRARIES:PATH='-L'$MKLROOT' -lmkl_core -lmkl_intel_ilp64 -lmkl_sequential -lmkl_blas95_ilp64' \
    -DLAPACK_LIBRARIES:PATH='-lmkl_lapack95_ilp64' \
    -DENABLE_NETCDF:BOOL=ON \
    -DNetCDF_DIR:PATH='/fusion/usc/opt/netcdf/netcdf-4.4.1_mpich-3.2_intel2018/' \
    -DENABLE_HDF5:BOOL=ON \
    -DENABLE_PARALLEL:BOOL=ON \
    -DCMAKE_BUILD_TYPE:STRING=RELEASE \
    ..
