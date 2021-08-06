


This gives installation instructions to EFIT using the CMake build system

Quickstart
==========

Minimum requirements: cmake (>=3.8), fortran compiler, blas and lapack fortran
libraries

If these are already in your path or other standard locations, you can ignore
all other third-party libraries and simply build with::

    mkdir build
    cd build
    cmake ..
    make 
    make test

This will work on Ubuntu (18.04 and 20.04) as long as the following
packages have been installed (call apt-get to install)::
    git (optional)
    build-essential
    cmake
    gfortran
    libblas-dev or libopenblas-dev (or libopenblas-serial-dev)
    liblapack-dev or libopenblas-dev (or libopenblas-serial-dev)
(18.04 only builds with libopenblas-dev - not libblas-dev)
(16.04 and older requires non-standard CMake and gfortran versions)

It will also work on Mac with similar libraries (TO DO: list) (can be obtained with
brew install)

Configuring third-party libraries
=================================

To build with all of the libraries and specify the compilers, one can use a
longer shell script to store the required flags.  Here is an example 
`config.sh` with all of the bells and whistles turned on (assuming everything
is built in your $HOME/software directory)::

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
        -DBLAS_DIR:PATH='$HOME/software/blas' \
        -DLAPACK_DIR:PATH='$HOME/software/lapack' \
        -DENABLE_NETCDF:BOOL=ON \
        -DNetCDF_DIR:PATH='$HOME/software/netcdf' \
        ..

The following flags are still underdevelopment::

        -DHDF5_DIR:PATH='$HOME/software/hdf5' \
        -DMDSPLUS_DIR:PATH='$HOME/software/mdsplus' \

For debugging, set:: 

        -DCMAKE_BUILD_TYPE:STRIND=Debug

Config scripts for a number of systems have already been made and can be found
in the `share/config_examples/` directory.  The required dependencies for
these systems is described in those scripts as well (for best results read
before executing).  Please consider adding yours if it can help out other users.
