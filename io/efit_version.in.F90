!-----------------------------------------------------------------------
!     efit_version
!-----------------------------------------------------------------------
      SUBROUTINE efit_version(efitver)
      IMPLICIT NONE
      CHARACTER(8), INTENT(OUT) :: efitver
      CHARACTER(512) :: giturl,gitbranch,githash,lmods
      CHARACTER(512) :: fc,fcid,fcver
      CHARACTER(5096) :: fcflags
      INTEGER, DIMENSION(8) :: tval

      giturl= "GIT_URL"
      gitbranch= "GIT_BRANCH"
      githash= "GIT_COMMIT_HASH"
      lmods= "GIT_MODIFICATIONS"
      fcid= "CMAKE_Fortran_COMPILER_ID"
      fcver= "CMAKE_Fortran_COMPILER_VERSION"
      fcflags= "CMAKE_Fortran_FLAGS"
      fcflags= TRIM(fcflags)//" CMAKE_Fortran_BUILD_TYPE_FLAGS"
!-----------------------------------------------------------------------
!     write version #
!-----------------------------------------------------------------------
      WRITE(6,'(2a,/)') 'Git repository url ',TRIM(giturl)
      WRITE(6,'(2a,/)') 'Git branch ',TRIM(gitbranch)
      WRITE(6,'(4a,/)') 'Git hash ',TRIM(githash),' ',TRIM(lmods)
!-----------------------------------------------------------------------
!     simulation date and time
!-----------------------------------------------------------------------
      CALL DATE_AND_TIME(VALUES=tval)
      WRITE(6,'(a,i2,a,i2,a,i4,a,i2,a,i2,a,i2,/)')                      &
        'Simulation started on ',tval(2),'/',tval(3),'/',tval(1),       &
        ' at ',tval(5),':',tval(6),':',tval(7)
!-----------------------------------------------------------------------
!     compiler and flags
!-----------------------------------------------------------------------
      WRITE(6,'(4a,/)') 'Code compiled with ',TRIM(fcid),' ',TRIM(fcver)
      WRITE(6,'(2a,/)') ' and flags ',TRIM(fcflags)
!-----------------------------------------------------------------------
!     set output variables to efit match efit convention
!     use git hashes (first 8 or 10 chars) instead of compilation dates
!-----------------------------------------------------------------------
      efitver=githash(1:8)
!-----------------------------------------------------------------------
!     terminate
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE efit_version
