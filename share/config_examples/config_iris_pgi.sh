# In order to run this version of EFIT on the system you must
#   execute the following module commands: (not required to build)
#
#    module rm efit
#    module switch gcc-4.7.2 gcc-9.2.0
#    module load pfg/18.7

    module load cmake/3.8.2
    export CC=/fusion/usc/opt/pgi/linux86-64/18.7/bin/pgcc

    cmake \
    -DBLAS_LIBRARIES:PATH='/fusion/usc/opt/pgi/linux86-64/18.7/lib/libblas.a' \
    -DLAPACK_LIBRARIES:PATH='/fusion/usc/opt/pgi/linux86-64/18.7/lib/liblapack.a' \
    -DENABLE_NETCDF:BOOL=ON \
    -DNetCDF_DIR:PATH='/fusion/usc/opt/env/gcc9.2-pgf18.3/' \
    ..
