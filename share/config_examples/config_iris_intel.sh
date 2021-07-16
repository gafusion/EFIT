# In order to run this version of EFIT on the system you must
#   execute the following module commands: (not required to build)
#
#    module rm efit
#    module switch gcc-4.7.2 gcc-9.2.0
#    module load intel/2019
#    export LD_LIBRARY_PATH=/fusion/usc/opt/intel2019/compilers_and_libraries/linux/lib/intel64:$LD_LIBRARY_PATH

    module load cmake/3.8.2
    export CC=/fusion/usc/opt/intel2019/bin/icc

    cmake \
    -DBLAS_LIBRARIES:PATH='-L'$MKLROOT' -lmkl_core -lmkl_intel_ilp64 -lmkl_sequential -lmkl_blas95_ilp64' \
    -DLAPACK_LIBRARIES:PATH='-lmkl_lapack95_ilp64' \
    ..
