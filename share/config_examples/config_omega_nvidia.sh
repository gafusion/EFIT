# In order to build and run this version of EFIT on the system you must
#   execute the following module commands: (not required to build)
#
#    module load env/nvhpc23.11-mpi mdsplus
#
# If you don't want MPI (slower in serial) simply remove the -DENABLE_PARALLEL... line
#   and set FC=$MPIF90
#
# Slurm also does not appear to be compatible with this compiler so you
#   need to run MPI jobs with mpirun at this time...

    export FC=$MPIF90

    cmake \
    -DMPICMD:STRING='mpirun -n ' \
    -DBLAS_LIBRARIES:PATH="${BLAS_DIR}/lib/libblas.a" \
    -DLAPACK_LIBRARIES:PATH="${BLAS_DIR}/lib/liblapack.a" \
    -DENABLE_NETCDF:BOOL=ON \
    -DNetCDF_C_DIR:PATH=$NETCDF_DIR \
    -DNetCDF_FORTRAN_DIR:PATH=$NETCDF_F_DIR\
    -DENABLE_HDF5:BOOL=TRUE \
    -DHDF5_ROOT:PATH=$HDF5_DIR \
    -DENABLE_MDSPLUS:BOOL=ON \
    -DMDSPLUS_DIR:PATH='/fusion/usc/c8/opt/mdsplus/stable/7.131.1' \
    -DD3_LIB:PATH='/fusion/projects/codes/efit/c8/dev/d3lib_nvhpc23.11/libd3share.a' \
    -DMSE_LIB:PATH='/fusion/projects/codes/mse/lib/libmse.a' \
    -DENABLE_PARALLEL:BOOL=ON \
    -DCMAKE_BUILD_TYPE:STRING=RELEASE \
    ..

