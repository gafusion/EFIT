.. highlight:: bash

Installation
============

An editable tutorial document describing the installation process can be found at:
https://docs.google.com/document/d/1S8nYPoQVrPGiOJK7EuyNVBG61W51oXVofmNpmNoQo20

Public Installations
--------------------

Public installations of the EFIT-AI are available on the GA iris cluster, PPPL portal cluster
and nersc cori supercomputer

iris (loads DIII-D Green functions by default, to use others set the environment variable link_efit=/fusion/projects/codes/efit/efitai/efit_support_files/{machine}/ after loading the module)::

    module purge
    module load efitai/{pgi pgi_ser gnu gnu_ser intel intel_ser}
    efit {grid_size}

cori (intallation is pending ERCAP setup) (you will need to be added to the project repo in order to access these installations - email kruger@txcorp.com) ::

    module swap PrgEnv-${PE_ENV,,} PrgEnv-{gnu,intel,cray,nvidia,aocc}
    export link_efit=/global/common/software/efitai/efit_support_files/{device}/
    /global/common/software/efitai/{cori,perlmutter}/{gnu,intel,cray,nvidia,amd}{_ser}/efit/efit {grid_size}

portal::

    module purge
    module load gcc/9.3.0
    export link_efit=/p/nstx/EFIT_GA/efit_support_files/{device}/
    /p/nstx/EFIT_GA/efit/build/efit/efit {grid_size}

Installing from source
----------------------

This gives installation instructions to EFIT using the CMake build system.
`CMake <https://cmake.org>`__ provides many additional tools (GUI, TUI) that
allows for more advanced workflows.  Here, we provide the simplest workflow.

Minimum requirements: 
   + cmake (>=3.8) 
   + fortran compiler 
   + blas and lapack libraries

If these are already in your path or other standard locations, you can ignore
all other third-party libraries and simply build with::

    cd $EFIT_ROOT
    mkdir build
    cd build
    cmake ..
    make 

where ``EFIT_ROOT`` is an environment variable set to the EFIT source code
directory.
This will work on most systems and has been specifically tested on Mac 
and Ubuntu (18.04 and 20.04).   Note that in this example, the build directory
is a subdirectory of ``EFIT_ROOT``, but it can be located anywhere on the
filesystem.  The only thing that would change, is to replace the ``cmake ..`` with
``cmake $EFIT_ROOT``.  See the CMake documentation for more details.

For Ubuntu, the following packages need to be installed (call apt-get to
install)::

    git (optional)
    build-essential
    cmake
    gfortran
    libblas-dev or libopenblas-dev (or libopenblas-serial-dev)
    liblapack-dev or libopenblas-dev (or libopenblas-serial-dev)

(18.04 only builds with libopenblas-dev - not libblas-dev)
(16.04 and older requires non-standard CMake and gfortran versions)

For MacOS, XCode needs to be installed.  We recommend that one use brew to
install the gfortran library::

    brew install gcc

Once gfortran is in your path (follow instructions), then it can use the
Accelerate library from XCode for blas and lapack.   Here is an example CMake
configure script that is more complicated the just invoking the ``cmake ..`` as
shown above::

    #!/bin/sh

    rm -rf CMake*

    cmake \
        -DCMAKE_OSX_DEPLOYMENT_TARGET:STRING='10.9' \
        ..

additional configure options shown below can be added to any OS.


Configuring third-party libraries
---------------------------------

To build with all of the libraries and specify the compilers, one can use a
longer shell script to store the required flags.  Here is an example 
``config.sh`` with all of the bells and whistles turned on (assuming everything
is built in your $HOME/software directory)::

    #!/bin/sh

    export PATH='$HOME/software/hdf5':${PATH}
    
    cmake \
        -DENABLE_DOCS:BOOL=FALSE \
        -DCMAKE_INSTALL_PREFIX:PATH=$HOME/software \
        -DCMAKE_BUILD_TYPE:STRING=Release \
        -DCMAKE_COLOR_MAKEFILE:BOOL=FALSE \
        -DCMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
        -DENABLE_SHARED_LIBS:BOOL=TRUE \
        -DCMAKE_C_COMPILER:FILEPATH='gcc' \
        -DCMAKE_Fortran_COMPILER:FILEPATH='gfortran' \
        -DCMAKE_C_FLAGS:STRING='-fvisibility=default -fPIC -pipe' \
        -DCMAKE_Fortran_FLAGS:STRING='-fPIC -pipe' \
        -DBLAS_DIR:PATH='$HOME/software/blas' \
        -DLAPACK_DIR:PATH='$HOME/software/lapack' \
        -DENABLE_PARALLEL:BOOL=ON \
        -DMPICMD:STRING='mpirun -n ' \
        -DNPROC:STRING=2 \
        -DENABLE_NETCDF:BOOL=ON \
        -DNetCDF_DIR:PATH='$HOME/software/netcdf' \
        -DENABLE_HDF5:BOOL=ON \
        -DHDF5_ROOT:PATH='$HOME/software/hdf5' \
        -DENABLE_MDSPLUS:BOOL=ON \
        -DD3_LIB:PATH='/fusion/projects/codes/efit/dev/d3lib_gcc9.2.0/libd3share.a' \
        -DMSE_LIB:PATH='/fusion/projects/codes/mse/lib/libmse.a' \
        -DTEST_EFUND:BOOL=True \
        ..

For debugging, set:: 

        -DCMAKE_BUILD_TYPE:STRING=Debug

Config scripts for a number of supercomputers and compilers have already been made
and can be found in the ``share/config_examples/`` directory, including::

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

If you are trying to build for the first time on a different supercomputer
or with a different compiler, the best starting point is to change
environment library paths from an existing configure script (e.g. try the
most similar or ``iris_gnu.sh`` first) to match what is available.  If you run
into problems, contact a developer.

To ensure your build was successful, it is recommended that you run the included
tests.  See `quickstart <quickstart>`_ for more info.

Once you have successfully built on a different system/compiler, please add your
working script to the collection in ``$EFIT_ROOT/share/config_examples`` to aid
future users.

