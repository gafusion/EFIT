.. _quickstart:
.. highlight:: bash

Quick-start
===========

This gives instructions for running EFIT once the executables are built.

Running tests
-------------

If you are maintaining your own installation (e.g. not using a system install
or someone else's build) then it is recommended that you ensure all tests pass
before running any other cases.

The packaged tests are run by executing ``make test`` from the build directory.  
If you have a parallel build on a supercomputer, then you should launch an 
interactive session on a compute node before executing the tests.  By default,
parallel tests are executed with slrum, but if your system uses a different
program (e.g. mpirun, mpiexec) you need to configure the build with::

    -DMPCMD:STRING='mpirun -n ' 
    -DNPROC:STRING=2

to tell the test harness how to run the test.

For information about why a test failed, see 
``${EFIT_BUILDDIR}/Testing/Temporary/LastTest.log``.

You can also execute::

    make test VERBOSE=1

for more detailed running information.

To understand in detail how a test is run, one can change into the directory of
each test and examine the shell script that invokes that test.  For example, a
simple example of running EFIT can be seen by::

    cd ${EFIT_BUILDDIR}/test/DIIID/efit01
    ./run_efit.sh   # Execute script to run efit
    more ./run_efit.sh  # See how we execute in detail

where ``${EFIT_BUILDDIR}`` refers to the location of the EFIT build directory as
discussed in the `install <install>`_ section.

Executable location
-------------------

The key executables are::

    ${EFIT_BUILDDIR}/efit/efit
    ${EFIT_BUILDDIR}/green/efund

The efund executable generates the Green's function tables as we discuss next.

Forming Green's function tables
--------------------------------

An ``mhdin.dat`` file describing the machine geometry must be in the
same location that you run ``efund`` from.

Once this requirement is satisfied, the tables can be formed by calling
``efund``.  Example::

    green/efund

The tables are contained in files.  Example filenames of these files are 
``btcomp.dat, ccoil.ddd, ccomp.dat, n1coil.ddd, icomp.dat``.  

Running EFIT
------------

In order for EFIT to run it must know where the Green function tables and
machine data files can be found.  Default locations can be specified at compile
time by passing the flags ``-DINPUT_DIR=``, ``-DTABLE_DIR=``, and ``-DSTORE_DIR`` to
cmake.  Otherwise the default is set as the build directory.  This is superceded
at runtime by environment variables or input files.  

The precedence is decided by checking in that order:
     #.  efit.input file (setup namelist)
     #.  input files (efit_snap.dat, k-files, etc.)
     #.  environment variables
     #.  build defaults

The environment variables and setup namelist in efit.input are specified by ``link_efit`` and ``link_store``::

     TABLE_DIR=${link_efit}/green 
     INPUT_DIR=${link_efit}
     STORE_DIR=${links_store}


Input files can specify the ``TABLE_DIR``, ``INPUT_DIR``, and ``STORE_DIR``
variables in ``IN1``, ``INWANT``, and ``EFITIN`` namelists (kfiles, rfiles, xfiles,
and snap files).

Once the tables and experimnt data are properly located, EFIT can be run
by calling ``efit`` with the mesh size that matches the tables and an 
efit.input file or interactively.
Example::

    efit/efit 65
    2
    1
    kfile

For more information on the input variables see `namelists <namelists>`_.  

A description of the different input and output types can be found at
`EFIT IO Files <https://fusion.gat.com/theory/Efitiofiles>`__ (requires GA login).

Running EFIT with OMFIT
-----------------------

The recommedned tool for developing more compilicated workflows and analyzing results
is `OMFIT <https://omfit.io/>`__.  Interested users should refer to the documentation
on that page for setting up and using their software.  EFIT-AI (this version) can be
used in place of legacy EFIT within the 
`EFIT <https://omfit.io/modules/mod_EFIT.html/>`__, 
`EFITtime <https://omfit.io/modules/mod_EFITtime.html/>`__, and
`kineticEFITtime <https://omfit.io/modules/mod_kineticEFITtime.html/>`__, modules by
selecting the check-box in the GUI before executing.  Additional functionality and
integration is on going.  These modules have a variety of built in plotting options and 
more detailed analysis can be performed with the EFITviewer suite.
