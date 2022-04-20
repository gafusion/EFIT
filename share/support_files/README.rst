This directory is intended to provide a central location for all the machine dependent files required to run efund and subsequently efit.  It has the same structure as the support file directories provided on GA clusters and NERSC, but does not contain any pre-computed Green function tables.  The logic used in efit to read such tables is described in detail in::
  efit/efit_tables.F90 -> set_table_dir

To summarize, the subdirectories of <machine>/green are labeled by a shot number corresponding to the first shot which should use the contained mhdin.dat (or pre-compiled tables) should be used.  Therefore, to run a choosen shot number, you should use the mhdin.dat (or tables) from the directory with the largest number that is less than the shot.
