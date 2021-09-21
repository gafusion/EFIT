Quick-start
===========

This gives instructions for running EFIT once the executables are built.

Running tests
-------------

If you are maintaining your own installation (e.g. not using a system install
or someone else's build) then it is recommended that you ensure all tests pass
before running any other cases.

The packaged tests are run by executing `make test` from the build directory.  
If you have a parallel build on a supercomputer, then you should launch an 
interactive session on a compute node before executing the tests.

For information about why a test failed, see 
`Testing/Temporary/LastTest.log`.

The complete ouput from each test can be found in the `test` subdirectrory of
the build (not the test source directory with the same name a level up).

Executable location
-------------------

By default the executables are placed in the `green/` and `efit/`
subdirectories of the build directory (not the source directories of the same
name).

Forming Green function tables
-----------------------------

An mhdin.dat file describing the machine geometry must be in the
same location that you run `efund` from.

Once this requirement is satisfied, the tables can be formed by calling
`efund` with the target mesh size, specified as a single integer describing
both R and Z dimensions.  Example::

    green/efund 129

This also creates a `dprobe.dat` file containing all the limiter and other
necessary information about the machine.

Running EFIT
------------

In order for EFIT to run it must know where the Green function tables and
machine data files can be found.  Default locations can be specified at compile
time by passing the flags `-DINPUT_DIR=`, `-DTABLE_DIR=`, and `-DSTORE_DIR` to
cmake.  Otherwise the default is set as the build directory.  This is superceded
at runtime by environment variables or input files.  Which has precedence is decided by checking in that order::
input files overwrite environment variables which overwrite the defaults
The environment variables are specified by `link_efit` (assigns 
`TABLE_DIR=${link_efit}/green` and `INPUT_DIR=${link_efit}`) and `link_store`
(assigns `STORE_DIR=${link_store}`).  Input files can specify the `TABLE_DIR`,
`INPUT_DIR`, and `STORE_DIR` variables in `IN1`, `INWANT`, and `EFITIN`
namelists (kfiles, rfiles, and snap files).

Once the tables and machine data are properly located, EFIT can be run
interactively by calling `efit` with the mesh size that matches the tables.
Example::

    efit/efit 65
    2
    1
    kfile

For more information on the input variables see https://fusion.gat.com/theory/Efitin1 (requires GA login).  A description of the different input and output
types can be found at https://fusion.gat.com/theory/Efitiofiles .

