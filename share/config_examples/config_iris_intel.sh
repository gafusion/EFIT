# In order to build and run this version of EFIT on the system you must
#   execute the following module commands:
#
#    module purge
#    module load env/gcc9.2 intel/2018 mpich/3.2-intel2018
#
# If you don't want MPI (slower in serial) simply remove the FC=...
#   and -DENABLE_PARALLEL... lines

    module load cmake/3.8.2
    module load hdf5/1.8.19-mpich3.2-intel2018
    export CC=/fusion/usc/opt/intel2018/bin/icc
    export FC=/fusion/usc/opt/mpich/mpich-3.2/intel2018/bin/mpifort

    cmake \
    -DMPICMD:STRING='srun --mpi=pmi2 -n ' \
    -DBLAS_LIBRARIES:PATH='-L'$MKLROOT' -lmkl_core -lmkl_intel_ilp64 -lmkl_sequential -lmkl_blas95_ilp64' \
    -DLAPACK_LIBRARIES:PATH='-lmkl_lapack95_ilp64' \
    -DENABLE_NETCDF:BOOL=ON \
    -DNetCDF_DIR:PATH='/fusion/usc/opt/netcdf/netcdf-4.4.1_mpich-3.2_intel2018/' \
    -DENABLE_HDF5:BOOL=ON \
    -DD3_LIB:PATH='/fusion/projects/codes/efit/dev/d3lib_gcc9.2.0/libd3share.a' \
    -DENABLE_MDSPLUS:BOOL=ON \
    -DMDSPLUS_DIR:PATH='/fusion/usc/src/mdsplus/mdsplus-stable_release-7-96-9' \
    -DMSE_LIB:PATH='/fusion/projects/codes/mse/lib/V4_05-gcc9.2.0/libmse.a' \
    -DENABLE_PARALLEL:BOOL=ON \
    -DCMAKE_BUILD_TYPE:STRING=RELEASE \
    ..
