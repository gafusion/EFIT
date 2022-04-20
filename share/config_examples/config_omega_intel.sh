# In order to build and run this version of EFIT on the system you must
#   execute the following module commands:
#
#    module load intel/2018
#
# *the module does not match the compiler used but that is the only way that a working compilation can be achieved on Omega right now...
#
# If you don't want MPI (slower in serial) simply remove the FC=...
#   and -DENABLE_PARALLEL... lines
#
# Omega does not appear to have MPI, NetCDF, HDF5, or MDS+ support for Intel at this time

    export CC=/fusion/usc/opt/intel2020/bin/icc

    cmake \
    -DMPICMD:STRING='srun --mpi=pmi2' \
    -DBLAS_LIBRARIES:PATH='-L'$MKLROOT' -lmkl_core -lintlc -lmkl_intel_ilp64 -lmkl_sequential -lmkl_blas95_ilp64' \
    -DLAPACK_LIBRARIES:PATH='-lmkl_lapack95_ilp64' \
    -DCMAKE_BUILD_TYPE:STRING=RELEASE \
    ..
