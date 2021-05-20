    /bin/rm -rf CMake*

    cmake \
    -DENABLE_SHARED_LIBS:BOOL=TRUE \
    -DCMAKE_C_COMPILER:FILEPATH='gcc' \
    -DCMAKE_Fortran_COMPILER:FILEPATH='gfortran' \
    -DCMAKE_C_FLAGS:STRING='-fvisibility=default -fPIC -pipe' \
    -DCMAKE_Fortran_FLAGS:STRING='-fPIC -pipe' \
    -DHAVE_NETCDF:BOOL=TRUE \
    -DNETCDF_DIR:PATH='/fusion/usc/opt/env/gcc9.2' \
    -DBLAS_DIR:PATH='/fusion/usc/opt/env/gcc9.2' \
    -DBLAS_LIBS:PATH='/fusion/usc/opt/env/gcc9.2' \
    -DLAPACK_LIBRARIES:'PATH=/fusion/usc/opt/env/gcc9.2' \
    ..
