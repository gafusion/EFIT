# In order to build and run this version of EFIT on the system you must
#   execute the following module commands:
#
#    module load env/gcc11.x-mpi mdsplus
#
# If you don't want MPI (slower in serial) remove the -DENABLE_PARALLEL... line
#   and set:
# export FC=/opt/rh/gcc-toolset-11/root/usr/bin/gfortran
# -DBLAS_LIBRARIES:PATH='/fusion/usc/c8/opt/env/gcc-11.x/openblas/0.3.24/lib/libopenblas.a' \
# -DLAPACK_LIBRARIES:PATH='/fusion/usc/c8/opt/env/gcc-11.x/openblas/0.3.24/lib/libopenblas.a' \
# -DNetCDF_C_LIBRARY='/fusion/usc/c8/opt/env/gcc-11.x/netcdf-c/4.9.0' \
# -DNetCDF_DIR:PATH='/fusion/usc/c8/opt/env/gcc-11.x/netcdf-fortran/4.6.1' \
# -DHDF5_ROOT:PATH='/fusion/usc/c8/opt/env/gcc-11.x/hdf5/1.12.2'\

    export CC=/opt/rh/gcc-toolset-11/root/usr/bin/gcc
    export FC=/fusion/usc/c8/opt/env/gcc-11.x/mpich/4.1.2/bin/mpifort

    cmake \
    -DMPICMD:STRING='srun --mpi=pmi2 -n ' \
    -DBLAS_LIBRARIES:PATH='/fusion/usc/c8/opt/env/gcc-11.x/openblas-mpi/0.3.24/lib/libopenblas.a' \
    -DLAPACK_LIBRARIES:PATH='/fusion/usc/c8/opt/env/gcc-11.x/openblas-mpi/0.3.24/lib/libopenblas.a' \
    -DENABLE_OpenMP:BOOL=ON \
    -DENABLE_NETCDF:BOOL=ON \
    -DNetCDF_C_LIBRARY='/fusion/usc/c8/opt/env/gcc-11.x/netcdf-c-mpi/4.9.0' \
    -DNetCDF_DIR:PATH='/fusion/usc/c8/opt/env/gcc-11.x/netcdf-fortran-mpi/4.6.1' \
    -DENABLE_HDF5:BOOL=ON \
    -DHDF5_ROOT:PATH='/fusion/usc/c8/opt/env/gcc-11.x/hdf5-mpi/1.12.2'\
    -DD3_LIB:PATH='/fusion/usc/c8/src/d3lib/lib/libd3.a' \
    -DENABLE_MDSPLUS:BOOL=ON \
    -DMDSPLUS_DIR:PATH='/fusion/usc/c8/opt/mdsplus/alpha/7.130.1' \
    -DENABLE_PARALLEL:BOOL=ON \
    -DCMAKE_BUILD_TYPE:STRING=RELEASE \
    ..

#    -DMSE_LIB:PATH='/fusion/projects/codes/mse/lib/V4_05-gcc8.x/libmse.a' \
