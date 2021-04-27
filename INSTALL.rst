


This gives installation instructions to EFIT using the CMake build system

Quickstart
==========

The easiest method to build is to ignore all of the third-party libraries::

    mkdir build
    cd build
    cmake ..

To build with all of the libraries and specify the compilers, one can use a
longer shell script to modify.  Here is an example `config.sh` with all of the 
bells and whistles turned on (TODO:  Need to test)::

    #!/bin/sh

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
        -DMDSPLUS_DIR:PATH='$HOME/software/mdsplus' \
        ..


For debugging, set the `CMAKE_BUILD_TYPE` to `Debug`.

