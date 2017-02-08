!
!  Module with routines to write into an MdPlus database.
! 
   module mds_write_rg
   
   use mds_rg

   contains


      subroutine mds_ld_error (ier,cfunc) 
      implicit none
      integer ier
      character(len=*) cfunc
      if (ier .ne. 0) write (*,*) ' Error on MDS load for ', cfunc
      end subroutine mds_ld_error

 
      ! Write ntimes pts of real data into MDSplus for node name.
      ! Written: 1/23/03 by rjg
      logical function Write_Mds_Real_Arr (name, ntimes, data, cunits_in)
      Implicit none
      Include 'mdslib.inc'


      Character(len=*), intent(in) :: name   ! Name of mds node, e.g. 'ALP_ERR'
      Integer, intent(in) ::  ntimes         ! Number of data points.
      Real, dimension(ntimes), intent(in) :: data ! Data
      Character(len=*), intent(in) :: cunits_in ! Units to be writtten with data

      Integer :: istat
      integer :: desc1, desc2
      character(len=512) :: expression
      character(len_trim(name)+1) :: node
      character(len_trim(cunits_in)+1) :: cunits
      integer :: ilen


      ! The basic call is:
      !   istat = MdsPut (node, expression, (descr1,descr2, ...), 0)
      !          where node is an Mdsplus node name
      !          expression is something to be done to the data
      !          there is one descriptor (descr1, etc) for every "$" sign
      !          in the expression and the final "0" is required for the
      !          fortran call.

      ! Remember to terminate every character string with a null, for C
      node = trim(name) // char(0)
      ! The TIME is expression means to use the node called TIME which is
      ! parallel to node "NODE".
      expression = "BUILD_SIGNAL(BUILD_WITH_UNITS($,$),,TIME)" // Char(0)
      desc1 = descr(IDTYPE_FLOAT,data,ntimes,0)
      cunits = trim(cunits_in) // char(0)
      ilen = len(cunits)-1
      desc2 = descr(IDTYPE_CSTRING, Cunits, 0, ilen, 0)
      
      ! Write the data
      istat = MdsPut (node, expression, desc1, desc2, 0 )

      Write_Mds_Real_Arr = MdsNoError(istat,'MdsPut',name)
      if (.not. Write_Mds_Real_Arr)                           &
         write (*,*) 'Write_Mds_Real_Arr: Could not write to ', node

      End function Write_Mds_Real_Arr


      ! This routine has been replaced and should no longer be needed.
      ! It does provide an example of how to obtain a multiplier and units
      ! stored in the tree.
      Subroutine Write_Mds_Variable_Ret (name, ntimes, data, cunits, ier)

      Implicit none

      Include 'mdslib.inc'


      Character(len=*) :: name        ! input:  name of mds variable, e.g. 'ALP_ERR'
      Integer      ntimes             ! input:  number of data point.
      Real         data(ntimes)        ! input:  data array
      Character*25 ctime              ! input:  time node

      Character(len=*) cunits
      Integer      ier                ! output: error code - even bad, odd good.
      Integer      istat
      Integer      i
      Integer      m
      integer      ilen

      istat = -999  


      m = len_trim(name)
      name = trim(name)

      ! Note that cunits is not actually used in this call.  The units
      ! used are those in the model tree.
      ! Just had a long talk with Justin, who explained all of this.
      ! The basic call is:
      !   istat = MdsPut (node, expression, (descr1,descr2, ...), 0)
      !          where node is an Mdsplus node name
      !          expression is something to be done to the data
      !          there is one descriptor (descr1, etc) for every "$" sign
      !          in the expression and the final "0" is required for the
      !          fortran call.
      !   In the example below, the expression uses the multiplier stored with
      !   the node to multiply the signal created by build_signal, which in
      !   turn uses a build with units.  The units are in the sub-node 
      !   "UNITS" of the desired node, which is why I have never been able
      !   to over-write them.  The timebase that is used is the data in the
      !   node "TIME".  Note that I really do not want to use the stored
      !   multiplier or units.
      istat = MdsPut(name(1:m)// Char(0), name(1:m)//":MULTIPLIER*"       &
     &               // "BUILD_SIGNAL(BUILD_WITH_UNITS($,"                &
     &               //  name(1:m)                                        &
     &               // ":UNITS),,time)" // Char(0),                      &
     &               descr(IDTYPE_FLOAT,data,ntimes,0), 0)

      ier = 0
      If (MdsError(istat,'MdsPut',name)) ier = 1

      End Subroutine Write_Mds_Variable_Ret



!---------------------------------------------------------------------
! Write ntimes pts of data into MDSplus for node name.
! Do not use multiplier
! Is this still used
      Subroutine Write_Mds_Rarray_No_Mult (name, ntimes, data, cunits, ier)
      Implicit none
      Include 'mdslib.inc'

      Character(len=*) :: name        ! input:  name of mds variable, e.g. 'ALP_ERR'
      Integer      ntimes             ! input:  number of data point.
      Real         data(ntimes)       ! input:  data array
      Character(len=*) cunits
      Integer      ier                ! output: error code - even bad, odd good.

      Integer      istat
      Integer      i
      Integer      m
      integer      ilen

!
      istat = -999  


      m = len_trim(name)
      name = name(1:m) // char(0)
      ilen = len_trim(cunits)
      cunits = cunits(1:ilen) // char(0) 
!         istat = mdsput (CFUNC, "BUILD_WITH_UNITS($,$)"//CHAR(0),          &
!     &                   descr(IDTYPE_FLOAT, array, isz, 0),               &
!     &                   descr(IDTYPE_CSTRING, Cunits, 0, 3),0)

      istat = mdsput (name, "BUILD_WITH_UNITS($,$)"//CHAR(0), &
     &              descr(IDTYPE_FLOAT, data, ntimes, 0),                 &
     &              descr(IDTYPE_CSTRING, Cunits, 0, ilen),0)
      if (MdsError (istat,'mdsput',name)) then
         ier = 1   
         call mds_ld_error (ier,name) 
      endif

      End Subroutine Write_Mds_Rarray_No_Mult 


      !---------------------------------------------------
      Subroutine Write_Mds_Real_Scalar (name, datum, cunits, ier)

      Implicit none

      Include 'mdslib.inc'

      Character(len=*) :: name        ! input:  name of mds variable, e.g. 'ALP_ERR'
      Real         datum              ! input:  datum (scalar)
      Character*20 cunits
      Integer      ier                ! output: error code - even bad, odd good.
      Integer      istat
      integer      ilen


      ilen = len_trim (cunits)
      cunits = trim(cunits) // char(0)
 
      istat = MdsPut(name//Char(0),                                        &
     &               "BUILD_WITH_UNITS($,$)"//Char(0),                    &
     &               descr(IDTYPE_FLOAT,  datum, 0),                       &
     &               descr(IDTYPE_CSTRING, cunits, 0, ilen),0)


!     How you would put data with no units.
!      istat = MdsPut(name//Char(0),                                        &
!     &               "$" // Char(0),                                       &
!     &               descr(IDTYPE_FLOAT, datum, 0), 0)



      ier = 0
      If (MdsError(istat,'MdsPut',name)) ier = 1
      
      End Subroutine Write_Mds_Real_Scalar 
  

!---------------------------------------------------------------------
! Write ntimes pts of data into MDSplus for node name.

      ! Disable until I write out units
 
      Subroutine Write_Mds_RScalar_Disabled (name, datum, ier)
      Implicit none

      Include 'mdslib.inc'


      Character(len=*) :: name        ! input:  name of mds variable, e.g. 'ALP_ERR'
      Real         datum              ! input:  datum (scalar)
      Integer      ier                ! output: error code - even bad, odd good.
      Integer      istat


      istat = MdsPut(name//Char(0),                                        &
     &               "$" // Char(0),                                       &
     &               descr(IDTYPE_FLOAT, datum, 0), 0)

      ier = 0
      If (MdsError(istat,'MdsPut',name)) ier = 1

      End Subroutine Write_Mds_RScalar_Disabled
  


      !---------------------------------------------------------------------
      ! Write an integer to MDSplus for node name.
      Subroutine Write_Mds_Int_Scalar (name, idatum, cunits, ier)
      Implicit none

      Include 'mdslib.inc'


      Character(len=*) :: name        ! input:  name of mds variable, e.g. 'ALP_ERR'
      Integer      idatum             ! input:  datum (scalar)
      Character(len=*) :: cunits      ! Units
      Integer      ier                ! output: error code - even bad, odd good.
      Integer      istat

      integer      ilen


      ilen = len_trim (cunits)
      cunits = trim(cunits) // char(0)
 
      istat = MdsPut(name//Char(0),                                        &
     &               "BUILD_WITH_UNITS($,$)"//Char(0),                    &
     &               descr(IDTYPE_LONG,  idatum, 0),                       &
     &               descr(IDTYPE_CSTRING, cunits, 0, ilen),0)


!      istat = MdsPut(name//Char(0),                                        &
!     &               "$" // Char(0),                                       &
!     &               descr(IDTYPE_LONG, idatum, 0), 0)

      ier = 0
      If (MdsError(istat,'MdsPut',name)) ier = 1

      End Subroutine Write_Mds_Int_Scalar
  


      !---------------------------------------------------------------------
      ! Write a character variable into MDSplus for node name.
      Subroutine Write_Mds_chars (name, cdatum, ier)
      Implicit none

      Include 'mdslib.inc'


      Character(len=*) :: name        ! input:  name of mds variable, e.g. 'ALP_ERR'
      Character*20 cdatum             ! input:  datum (scalar character string)
      Integer      ier                ! output: error code - 0 good,  1 bad.
      Integer      istat              ! error code - even bad, odd good.


      istat = MdsPut(name//Char(0),                                       &
     &               "$" // Char(0),                                      &
     &               descr(IDTYPE_CSTRING, cdatum, 0, 20), 0)

      ier = 0
      If (MdsError(istat,'MdsPut',name)) ier = 1

      End Subroutine Write_Mds_chars
  

!---------------------------------------------------------------------
! Write a real data value into MDSplus for node name.

      ! Disable until I write out units

      Subroutine Write_Mds_Real_Disabled (name, datum, ier)

      Implicit none

      Include 'mdslib.inc'

      Character(len=*) :: name        ! input:  name of mds variable, e.g. 'ALP_ERR'
      Real         datum              ! input:  datum (scalar)
      Integer      ier                ! output: error code - even bad, odd good.
      Integer      istat

      istat = MdsPut(name//Char(0),                                      &
     &               "$" // Char(0),                                     &
     &               descr(IDTYPE_FLOAT, datum, 0), 0)

      ier = 0
      If (MdsError(istat,'MdsPut',name)) ier = 1

!      istat = MdsPut(name // Char(0),                                    &
!     &                "$" // Char(0),                                    &
!     &               descr(IDTYPE_FLOAT, datum, 0), 0)



      Return
      End Subroutine Write_Mds_Real_Disabled
  

      ! Write an array of integers to mdsplus
      ! Written: circa 2000 by Terpstra
      ! Revised: 1/23/03 by rjg to write units
      logical function Write_Mds_Int_Arr (name, ntimes, idata, cunits)
      Implicit none
      Include 'mdslib.inc'

      Character(len=*), intent(in) :: name  ! name of mds node, e.g. 'ALP_ERR'
      Integer, intent(in) :: ntimes
      Character(len=*), intent(in) :: cunits
      Integer, dimension(ntimes), intent(in) :: idata   ! The data

      Integer      istat
      integer ::  desc1, desc2
      integer :: ilen

      desc1 =  descr(IDTYPE_LONG, idata,ntimes,0)
      ilen = len_trim(cunits)
      desc2 =  descr(IDTYPE_CSTRING, cunits(1:ilen)//char(0), 0, ilen ,0)


      istat = MdsPut (trim(name)//Char(0),                                 &
     &               "BUILD_SIGNAL(BUILD_WITH_UNITS($,$),,time)"//Char(0), &
     &               desc1, desc2, 0)

      Write_Mds_Int_Arr = (MdsNoError(istat,'MdsPut',name))
      if (.not. Write_Mds_Int_Arr)                                  &
         write (*,*) 'Write_Mds_Int_Arr: Could not write to ', name

      End Function Write_Mds_Int_Arr


!---------------------------------------------------------------------
! Write ntimes pts of data (unitless) into MDSplus for node name.
   
  ! I don't think this is being used
      Subroutine Write_Mds_Array_RETIRED (name, ntimes, data, ier)

      Implicit none

      Include 'mdslib.inc'


      Character(len=*) :: name        ! input:  name of mds variable, e.g. 'ALP_ERR'
      Real         data               ! input:  datum (scalar)
      Integer      ier                ! output: error code - even bad, odd good.
      Integer      ntimes
      Integer      istat


      istat = MdsPut(name // Char(0),                                     &
     &                "$" // Char(0),                                     &
     &               descr(IDTYPE_FLOAT, data, ntimes, 0), 0)


      ier = 0
      If (MdsError(istat,'MdsPut',name)) ier = 1


      End Subroutine Write_Mds_Array_RETIRED




      ! Write ntimes pts of character data into MDSplus for node name.
      ! New and not working yet
      ! Written: 1/24/03 by rjg
      ! Revised: 1/30/03 by rjg to use cdata_loc
!      logical function Write_Mds_Char_Arr (name, ntimes, cdata)
!      Implicit none
!      Include 'mdslib.inc'
!!
!      Character(len=*), intent(in) :: name        ! name of mds variable, e.g. 'ALP_ERR'
!      Integer, intent(in) :: ntimes
!      Character(len=*), dimension(ntimes), intent(in) :: cdata  ! Character data
!!
!      ! Want the elements in cdata_loc to be no longer than the longest 
!      ! element in cdata, but at least one character long (in order to
!      ! avoid code crashes.
!      character (max(1,maxval(len_trim(cdata(1:ntimes))))), dimension(ntimes) :: cdata_loc
!!
!      Integer :: istat
!      integer :: desc
!      integer :: length
!
!
!      length = len(cdata_loc(1))
!!
!! 
!!      print *,' max length in cdata = ', length
!
!      ! Load up Cdata_loc with the information in cdata
!      cdata_loc = cdata
!!
!      !  http://www.mdsplus.org/old/mdslib/mdslib.html#vararg
!      !  Creating a descriptor for a string array (length = space allocated 
!      !  for one member of the array, ntimes = number of elements in array): 
!!
!      desc = descr(IDTYPE_CSTRING, cdata_loc, ntimes, 0, length)
!!
!      istat = MdsPut(trim(name) // Char(0),                              &
!                    "$" // Char(0),                                     &
!                    desc, 0)
!
!      Write_Mds_Char_Arr = MdsNoError(istat,'MdsPut',name)
!      if (.not. Write_Mds_Char_Arr)                           &
!         write (*,*) 'Write_Mds_Char_Arr: Could not write to ', name
!
!      End function Write_Mds_Char_Arr
!
!
   end module mds_write_rg
