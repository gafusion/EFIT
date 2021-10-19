Installation
============

An editable tutorial document describing the installation process can be found at:
https://docs.google.com/document/d/1S8nYPoQVrPGiOJK7EuyNVBG61W51oXVofmNpmNoQo20

Public Installations
--------------------

Public installations of the EFIT-AI are available on the GA iris cluster, PPPL portal cluster
and nersc cori supercomputer

iris ::

    module switch gcc-4.7.2 gcc-9.2.0
    module load {intel/2018, pfg/18.7}
    module load {mpich/3.2-gcc9.2.0, mpich/3.2-intel2018, mpich/3.2-pgf18.7}
    export link_efit=/fusion/projects/codes/efit/efitai/efit_support_files/{device}/
    /fusion/projects/codes/efit/efitai/{gnu,intel,pgi}{_ser}/efit/efit {grid_size}

cori (intallation is pending ERCAP setup) (you will need to be added to the project repo in order to access these installations - email kruger@txcorp.com) ::

    module switch PrgEnv-gnu PrgEnv-{gnu,intel,cray}
    export link_efit=/global/common/software/efitai/efit_support_files/{device}/
    /global/common/software/efitai/efit/build_{gnu,intel,cray}/efit/efit {grid_size}

portal::

    module purge
    module load gcc/9.3.0
    export link_efit=/p/nstx/EFIT_GA/efit_support_files/{device}/
    /p/nstx/EFIT_GA/efit/build/efit/efit {grid_size}

Installing from source
----------------------

This gives installation instructions to EFIT using the CMake build system

Minimum requirements: cmake (>=3.8), fortran compiler, blas and lapack fortran
libraries

If these are already in your path or other standard locations, you can ignore
all other third-party libraries and simply build with::

    mkdir build
    cd build
    cmake ..
    make 

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
---------------------------------

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

        -DCMAKE_BUILD_TYPE:STRING=Debug

Config scripts for a number of supercomputers and compilers have already been made
and can be found in the `share/config_examples/` directory, including::

    config_iris_gnu.sh
    config_iris_intel.sh
    config_iris_pgi.sh
    config_nersc_gnu.sh
    config_nersc_intel.sh
    config_nersc_cray.sh
    config_portal.sh

The required environments for building on these systems are described in the scripts as well 
(for best results read before executing).

They can be used to install with the following commands::

    mkdir build
    cd build
    ../share/config_examples/config_{machine}_{compiler}.sh
    make

If you are trying to build for the first time on a different supercomputer or with a
different compiler, the best starting point is to change environment library paths from an 
existing configure script (e.g. try the most similar or iris_gnu.sh first) to match what
is available.  If you run into problems, contact a developer.

To ensure your build was successful, it is recommended that you run the included tests.
See RUN.rst for more info.

Once you have successfully built on a different system/compiler, please add your working
script to the collection to aid future users.

