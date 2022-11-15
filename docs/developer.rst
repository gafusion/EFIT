GitLab
=======

Day-to-day development happens at the `EFIT GitLab repository <https://gitlab.com/efit-ai/efit>`__.
There, you can find the history and development version of the source code,
`see or create issues <https://gitlab.com/efit-ai/efit/issues>`__,
`see or create merge requests <https://gitlab.com/efit-ai/efit/-/merge_requests>`__,
and more. 

For submitting the merge request, the developer should:
  
  + Set a reviewer to someone other than themself
  + Set the labels to the workflow stage
  + Before setting the workflow as ReadyToMerge, the developer should launch a
    pipeline (note that the play button has to be hit from the MR page) and
    ensure that the pipeline passes
  + The actual merge should always be done by someone other than the developer
    who submitted the MR


Slack and Mailing Lists
========================

Discussion occurs on the EFIT-AI slack channel.

Comments and questions regarding EFIT design and development should go to efit-physics@ice.txcorp.com.

Configuration issues should be reported to efit-physics@ice.txcorp.com.

For bugs and other issues, we encourage users to create a `GitLab <https://gitlab.com>`__ account and
`file issues directly <https://gitlab.com/efit-ai/efit/issues>`__;
this allows better tracking of long-term bug reports and feature requests.


Developer workflow
===================

Due to the small team, we do not use an advanced plan of gitlab which allows
formal workflow rules to be set up, but rather rely on conventions.  Our
conventions are:

  + All development is in a branch
  + Additional tests should be added to cover new functionality
  + Merge requests are submitted by developer
  + Merges to main must be done by another developer to ensure code review


Adding tests
=============

New tests will have many commonalities with the ones that already exist,
so a good starting point would be to copy an existing test directory and
replace the inputs and outputs.  All tests for a given experiment should use
the same Green tables that are in the test and if possible have the same
shot number.  Tests should include all of the following components when
applicable:

  + Input files (k, r, etc.) for 4 different time slices
  + Output files (g, m, a, etc.) for each time slice
  + OMAS files containing the inputs and outputs of both a single time
    slice and all 4 time slices

Note: OMAS files can be constructed with OMFIT using the `to_omas` script
in the EFITtime module and the `ods.save("/path/to/file.hdf5")` method.
The OMAS conversion does not know the expected dimensions and sizes for
arrays, however, so it cannot properly interpret the matrix syntax (e.g.
`CALPA(2,3) = ...`).  Therefore, all 2D arrays (`calpa`, `cgama`,
`ccoils`, etc.) should be flattened into 1D arrays with the same total
number of elements that EFIT is expecting (described in the `dprobe.dat`
file).

Once the correct input and output files are added, the following files
need to be edited to include the test:

  + ``$EFIT_ROOT/test/CMakeLists.txt`` (if a new experiment is being added)
  + ``$EFIT_ROOT/test/<experiment>/CMakeLists.txt`` add the subdirectory for
    the new test
  + ``$EFIT_ROOT/test/<experiment>/<new_test>/CMakeLists.txt`` can mostly be
    kept the same as other tests aside from changing the test name (shot
    number and time slice numbers should match the new test as well)


Fortran Guideines
=================

EFIT is a legacy fortran code originally written in the 1980's and before many
modern development tools and practices were in place.  Given this legacy,
modernizing EFIT while maintaining it's production status will take time.  As
developer's work on the code, the following rough guidelines are suggested (MR's
will not be rejected out-of-hand for violating these rules):


Use labels on blocks if possible
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fortran allows one to label `if` and `do` blocks to make them easier to read.
An example of an older style idiom is::

    if (a .eq. b .and. c .ne. d) then   
           <...>
    endif ! a .eq. b .and c . ne. d

The newer style idiom would be::

    in_bounds: if (a .eq. b .and. c .ne. d) then   
           <...>
    endif in_bounds

Here is an example for a do loop::

    i_loop:  do i = 1, nw
       j_loop: do j = 1, nh
          <...>
       enddo j_loop
    enddo i_loop

Using physics descriptions for loops and if blocks is ideal, but not always
possible or necessary.


Capitalization
~~~~~~~~~~~~~~

Some fortran codes capitalize intrinsic keywords, but EFIT does not.  So we use
`do` instead of `DO`, `subroutine` instead of `SUBROUTINE`, etc.


Variable typing
~~~~~~~~~~~~~~~~~

EFIT historically used fortran's implicit typing where reals and integers are
defined by the first character of the variable name; e.g., `error` would be a
real/double and `index` would be an integer.  We avoid this.  The 
rules of thumb is:

   + *Explicit typing are better than implicit typing*


Variable scoping
~~~~~~~~~~~~~~~~~

EFIT historically used common blocks which are no converted to global module
declarations.  Good programming practice is to use avoid this practice, and use
subroutine/function arguments to the extent possible.  The rules of thumb are:

   + *Explicit interfaces are better than implicit interfaces*
   + *Local declarations are better than global declarations*


