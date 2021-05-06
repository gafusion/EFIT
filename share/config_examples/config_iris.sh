

    cmake \
        -DCMAKE_INSTALL_PREFIX:PATH=$HOME/software \
        -DCMAKE_BUILD_TYPE:STRING=RelWithDebInfo \
        -DCMAKE_COLOR_MAKEFILE:BOOL=FALSE \
        -DCMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
        -DENABLE_SHARED_LIBS:BOOL=TRUE \
        -DCMAKE_C_COMPILER:FILEPATH='gcc' \
        -DCMAKE_Fortran_COMPILER:FILEPATH='gfortran' \
        -DCMAKE_C_FLAGS:STRING='-fvisibility=default -fPIC -pipe' \
        -DCMAKE_Fortran_FLAGS:STRING='-fPIC -pipe' \
        -DHDF5_DIR:PATH='$HOME/software/hdf5' \
        -DNETCDF_DIR:PATH='$HOME/software/netcdf' \
        -DMDSPLUS_DIR:PATH='/fusion/usc/src/mdsplus/mdsplus-stable_release-6-1-84' \
        -DBLAS_DIR:PATH=''\
        -DBLAS_LIBS:PATH=''\
        -DLAPACK_LIBRARIES:PATH=''\
        ..


#'/fusion/usc/opt/openblas/openblas-0.3.7'\
#'/fusion/usc/opt/openblas/openblas-0.3.7/lib'\
#'/fusion/usc/opt/openblas/openblas-0.3.7/lib'\