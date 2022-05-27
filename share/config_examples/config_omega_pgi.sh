# In order to build and run this version of EFIT on the system you must
#   execute the following module commands: (not required to build)
#
#    module load env/pgf20.11 mdsplus
#
# If you don't want MPI (slower in serial) simply remove the FC=...
#   and -DENABLE_PARALLEL... lines
#
# A compiler bug in NVIDIA's hpc_sdk 20.11 package prevents getdat from
#   being compiled.  hpc_sdk 21.1 or greater is required for snap modes.
#
# Slurm also does not appear to be compatible with this compiler so you
#   cannot run parallel builds at this time...

    export CC=/fusion/usc/opt/nvidia/hpc_sdk/Linux_x86_64/20.11/compilers/bin/pgcc

    cmake \
    -DMPICMD:STRING='srun --mpi=pmi2' \
    -DBLAS_LIBRARIES:PATH='/fusion/usc/c8/opt/env/nvf20.11/lib/libopenblas.a' \
    -DLAPACK_LIBRARIES:PATH='/fusion/usc/c8/opt/env/nvf20.11/lib/libopenblas.a' \
    -DENABLE_NETCDF:BOOL=ON \
    -DNetCDF_DIR:PATH='/fusion/usc/c8/opt/env/nvf20.11' \
    -DENABLE_HDF5:BOOL=TRUE \
    -DHDF5_ROOT:PATH='/fusion/usc/c8/opt/env/nvf20.11' \
    -DCMAKE_BUILD_TYPE:STRING=RELEASE \
    ..

#    export FC=/fusion/usc/c8/opt/env/nvf20.11/bin/mpifort
#    -DENABLE_PARALLEL:BOOL=ON \
#    -DD3_LIB:PATH='/fusion/projects/codes/efit/dev/d3lib_gcc9.2.0/libd3share.a' \
#    -DENABLE_MDSPLUS:BOOL=ON \
#    -DMDSPLUS_DIR:PATH='/fusion/usc/c8/opt/mdsplus/stable/7.131.1' \
#    -DMSE_LIB:PATH='/fusion/projects/codes/mse/lib/libmse.a' \
