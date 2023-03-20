# In order to build and run this version of EFIT on the system you must
#   execute the following module commands:
#
#    module load env/intel2020
#
# If you don't want MPI (slower in serial) simply remove the -DENABLE_PARALLEL... line
#   and set FC=/fusion/usc/opt/intel2020/bin/ifort
#
# ptdata calls hang when snap files are used (unclear why)

    export CC=/fusion/usc/opt/intel2020/bin/icc
    export FC=/fusion/usc/opt/mpich/mpich-3.2/intel2020/bin/mpifort

    cmake \
    -DMPICMD:STRING='srun --mpi=pmi2 -n ' \
    -DBLAS_LIBRARIES:PATH='-L'$MKLROOT' -lmkl_core -lintlc -lmkl_intel_ilp64 -lmkl_sequential -lmkl_blas95_ilp64' \
    -DLAPACK_LIBRARIES:PATH='-lmkl_lapack95_ilp64' \
    -DENABLE_NETCDF:BOOL=ON \
    -DNetCDF_DIR:PATH='/fusion/usc/opt/netcdf/netcdf-4.4.1_mpich-3.2_intel2020/' \
    -DENABLE_MDSPLUS:BOOL=ON \
    -DMDSPLUS_DIR:PATH='/fusion/usc/c8/opt/mdsplus/alpha/7.130.1' \
    -DENABLE_HDF5:BOOL=ON \
    -DENABLE_PARALLEL:BOOL=ON \
    -DCMAKE_BUILD_TYPE:STRING=RELEASE \
    ..

#    -DD3_LIB:PATH='/fusion/usc/c8/src/d3lib/lib/libd3.a' \
#    -DD3_LIB:PATH='/fusion/projects/codes/efit/dev/d3lib_intel2020/libd3share.a' \
