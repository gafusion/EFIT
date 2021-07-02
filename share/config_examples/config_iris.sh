    module load cmake/3.8.2
    module remove efit
    export CC=/fusion/usc/opt/gcc/gcc-9.2.0/bin/gcc

    cmake \
    -DBLAS_LIBRARIES:PATH='/usr/lib64/libblas.a' \
    -DLAPACK_LIBRARIES:PATH='/usr/lib64/liblapack.a' \
    -DENABLE_NETCDF:BOOL=ON \
    -DNetCDF_DIR:PATH='/fusion/usc/opt/env/gcc9.2' \
    ..
