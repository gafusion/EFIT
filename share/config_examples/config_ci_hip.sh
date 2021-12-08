# CI machine at hip.txcorp.com

    cmake \
        -DCMAKE_BUILD_TYPE:STRING=Debug \
        -DCMAKE_COLOR_MAKEFILE:BOOL=TRUE \
        -DCMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
        -DBLAS_blas_LIBRARY:FILEPATH=/scratch/soft/lib/libfblas.a \
        -DLAPACK_lapack_LIBRARY:FILEPATH=/scratch/soft/lib/libflapack.a \
        ..
