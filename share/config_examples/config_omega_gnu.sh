# In order to build and run this version of EFIT on the system you must
#   execute the following module commands:
#
#    module load env/gcc8.x mdsplus
#
# If you don't want MPI (slower in serial) simply remove the -DENABLE_PARALLEL... line
#   and set FC=/usr/bin/gfortran

    export CC=/usr/bin/gcc
    export FC=/fusion/usc/c8/opt/env/gcc-8.3.1/bin/mpifort

    cmake \
    -DMPICMD:STRING='srun --mpi=pmi2 -n ' \
    -DBLAS_LIBRARIES:PATH='/fusion/usc/c8/opt/env/gcc-8.3.1/lib/libopenblas.a' \
    -DLAPACK_LIBRARIES:PATH='/fusion/usc/c8/opt/env/gcc-8.3.1/lib/libopenblas.a' \
    -DENABLE_OpenMP:BOOL=ON \
    -DENABLE_NETCDF:BOOL=ON \
    -DNetCDF_DIR:PATH='/fusion/usc/c8/opt/env/gcc-8.3.1' \
    -DENABLE_HDF5:BOOL=ON \
    -DHDF5_ROOT:PATH='/fusion/usc/c8/opt/env/gcc-8.3.1'\
    -DD3_LIB:PATH='/fusion/usc/c8/src/d3lib/lib/libd3.a' \
    -DENABLE_MDSPLUS:BOOL=ON \
    -DMDSPLUS_DIR:PATH='/fusion/usc/c8/opt/mdsplus/alpha/7.130.1' \
    -DMSE_LIB:PATH='/fusion/projects/codes/mse/lib/V4_05-gcc8.x/libmse.a' \
    -DENABLE_PARALLEL:BOOL=ON \
    -DCMAKE_BUILD_TYPE:STRING=RELEASE \
    ..

