



This gives instructions for running EFIT once the executables are built

Executable location
===================

By default the executables are placed in the `green/` and `efit/`
subdirectories of the build directory (not the source directories of the same
 name).

Forming Green function tables
=============================

An mhdin.dat file describing the machine geometry must be in the
same location that you run `efund` from.

Once this requirement is satisfied, the tables can be formed by calling
`efund` with the target mesh size, specified as a single integer describing
both R and Z dimensions.  Example::

    green/efund 129

This also creates a `dprobe.dat` file containing all the limiter and other
necessary information about the machine.

Running EFIT
===========

By default, EFIT checks the build directory for the Green function tables and
`dprobe.dat` file.  This can be set to a different location in the input
file, a `shot_table.txt` file, or with environment variables.  Examples::
    
    (To Do)

Once the tables are properly located, EFIT can be run interactively by calling
`efit` with the mesh size that matches the tables.  Example::

    efit/efit 65
    2
    1
    kfile

For more information on the input variables see https://fusion.gat.com/theory/Efitin1 (requires GA login).  A description of the different input and output
types can be found at https://fusion.gat.com/theory/Efitiofiles .

Running tests
=============

The packaged tests for the code can be run by executing `make test` from the
build directory.  

For more information about why a test failed, see 
`Testing/Temporary/LastTest.log`.

The complete ouput from each test can be found in the `test` subdirectrory of
the build (not the test source directory with the same name).
