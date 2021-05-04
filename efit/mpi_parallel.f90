!-------------------------------------------------------------------------------
!< interface to F77 MPI include file
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!> module for interface to F77 MPI include file
!  defines all MPI variables via parameter statements
!  use this module to define machine-specific MPI datatypes
!-------------------------------------------------------------------------------
MODULE mpi_efit
  IMPLICIT NONE

  INCLUDE "mpif.h"
END MODULE mpi_efit
