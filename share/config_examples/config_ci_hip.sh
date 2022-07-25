#!/bin/bash
# CI machine at hip.txcorp.com

    export PATH=/scratch/soft/hdf5/bin:${PATH}
    export FC=/scratch/soft/mpich/bin/mpifort

    cmake \
        -DCMAKE_BUILD_TYPE:STRING=Debug \
        -DCMAKE_COLOR_MAKEFILE:BOOL=TRUE \
        -DCMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
        -DENABLE_PARALLEL:BOOL=ON \
        -DMPICMD:STRING='/scratch/soft/mpich/bin/mpiexec -n ' \
        -DBLAS_blas_LIBRARY:FILEPATH=/scratch/soft/lib/libfblas.a \
        -DBLAS_blas_LIBRARY:FILEPATH=/scratch/soft/lib/libfblas.a \
        -DLAPACK_lapack_LIBRARY:FILEPATH=/scratch/soft/lib/libflapack.a \
        -DSphinx_EXECUTABLE:FILEPATH=/scratch/soft/miniconda3/bin/sphinx-build \
        -DENABLE_HDF5:BOOL=TRUE \
        -DENABLE_DOCS:BOOL=TRUE \
        -DHDF5_DIR:PATH=/scratch/soft/hdf5 \
        ..

    # For future use:
