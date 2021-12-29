# In order to run this version of EFIT on the system you must
#   execute the following module commands: (not required to build)
#
#    module switch gcc-4.7.2 gcc-9.2.0
#    module load intel/2018
#    module load mpich/3.2-intel2018
#
# If you don't want MPI (slower in serial) simply remove the FC=...
#   and -DENABLE_PARALLEL... lines

    module load cmake/3.8.2
    export CC=/fusion/usc/opt/intel2018/bin/icc
    export FC=/fusion/usc/opt/mpich/mpich-3.2/intel2018/bin/mpifort

    cmake \
    -DBLAS_LIBRARIES:PATH='-L'$MKLROOT' -lmkl_core -lmkl_intel_ilp64 -lmkl_sequential -lmkl_blas95_ilp64' \
    -DLAPACK_LIBRARIES:PATH='-lmkl_lapack95_ilp64' \
    -DENABLE_NETCDF:BOOL=ON \
    -DNetCDF_DIR:PATH='/fusion/usc/opt/netcdf/netcdf-4.4.1_mpich-3.2_intel2018/' \
    -DENABLE_PARALLEL:BOOL=ON \
    -DENABLE_MDSPLUS:BOOL=ON \
    -DMDSPLUS_DIR:PATH='/fusion/usc/src/mdsplus/mdsplus-stable_release-7-96-9' \
    -DMSE_LIB:PATH='/fusion/projects/codes/mse/lib/libmse.a' \
    -DCMAKE_BUILD_TYPE:STRING=RELEASE \
    ..
