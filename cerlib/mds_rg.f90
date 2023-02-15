!
! Module with some common routines used to deal with mds database
! Take a look at https://fusion.gat.com/DIII-D/comp/analysis/mdsplus/f_api.html
!           or http://www.mdsplus.org/old/mdslib/mdslib.html
!  for info on the fortran api to mdsplus.
!
  Module mds_rg

!!!!!  use hpalias
  use mdslib
  implicit none

     logical, save :: connected_to_mds = .false.   ! True means we are
                                                   ! connected to mds server

  contains


     ! Subroutine to connect to mds server
     ! Written: 4/20/00 by rjg
     ! Revised: 1/12/01 by rjg to check if we are already connected
     ! Revised: 7/08/02 by rjg to add optional argument to return status
     ! Revised: 7/12/02 by rjg to solve logic problem introduced on 7/8/02
     subroutine wake_mds (lmds_awake)
                                                                       
     Implicit none

     logical, optional, intent(out) :: lmds_awake

     integer istat 

     ! If we are already connected, do not attempt to connect again to
     ! save wear and tear on the mds processes.
     if (.not. connected_to_mds) then

        ! mdsconnect returns the socket number of the connection.  This can be 
        ! odd or even.  It is the only exception to the rule for mds routines.
        ! Failure is indicated by a return value that is less than or equal to 
        ! zero.
        ! istat is the socket
        istat = Mdsconnect('atlas' // Char(0)) 
        connected_to_mds  = (istat .gt. 0)
        if (.not. connected_to_mds) print *, 'wake_mds: failed to connect'
     endif

     ! Return status of connection
     if (present(lmds_awake)) lmds_awake = connected_to_mds

     end subroutine wake_mds
     


     ! Subroutine to open a tree
     ! Return logical giving success or failure
     ! Written: 4/20/00 by rjg
     ! Revised: 4/27/01 by rjg to fix logic error on MdsError check
     subroutine wake_tree (tree, ishot, ltree_awake)
                                                                       
     Implicit none

     character(len=*), intent(in) :: tree  ! Name of tree to be opened;
                                           ! e.g., 'IONS' or 'ELECTRONS'
     integer, intent(in) :: ishot          ! Shot number of tree
     logical, intent(out) :: ltree_awake  ! TRUE means we successfully connected
                                      ! to mds server
                                      ! FALSE means we failed to connect
                                      ! to mds server.
     integer istat, ier, ilen
     
     ilen = len_trim(tree)
     istat = MdsOpen(tree(1:ilen)// Char(0), ishot) 
     ltree_awake = .not. MdsError (istat, 'MdsOpen', tree)

     end subroutine wake_tree
     


     ! Subroutine to close a tree
     ! Return logical giving success or failure
     ! Written: 4/20/00 by rjg
     subroutine close_tree (tree, ishot, ltree_awake)
                                                                       
     Implicit none

     character(len=*), intent(in) :: tree  ! Name of tree to be opened;
                                           ! e.g., 'IONS' or 'ELECTRONS'
     integer, intent(in) :: ishot          ! Shot number of tree
     logical, intent(out) :: ltree_awake  ! TRUE means we successfully connected
                                      ! to mds server
                                      ! FALSE means we failed to connect
                                      ! to mds server.
     integer istat, ier, ilen
     
     ilen = len_trim(tree)
     istat = MdsClose(tree(1:ilen)// Char(0), ishot) 
      Print *, 'After MdsClose: for shot ', ishot, ' tree  ', tree(1:ilen),  &
               '  istat = ', istat
     ltree_awake = MdsError (istat, 'MdsClose', tree)

     end subroutine close_tree
     


     ! Function to interpret status return from MDS function calls.
     ! Return TRUE if status is odd; else return FALSE.
     ! Written: 1/26/99 by rjg
     logical function good_status (status)
     implicit none
     integer*4 status

     integer num

     num = status/2
     num = status - 2*num

     if (num .eq. 1) then
        good_status = .true.
     else
        good_status = .false.
     endif

    END function good_status



     ! Written: 4/13/00 by rjg
     ! Return true if an error for istat
     logical function MdsError (istat,routine,name)
     implicit none
     integer, intent(in) :: istat   ! A return value from an mds function
     character(len=*) :: routine    ! Mds routine which returns the istat
     character(len=*) :: name       ! Mds node being loaded
     
     ! If istat is an even number, we have an error condition
     MdsError = (Mod(istat,2) .Eq. 0)
!     If (MdsError) then
!        print *,' Error from routine ', trim(routine), ' for ', trim(name)
!        print *,' istat was ', istat 
!     Endif
     End function MdsError 



     ! Written: 1/24/02 by rjg
     ! Return true if no error for istat
     logical function MdsNoError (istat,routine,name)
     implicit none
     integer, intent(in) :: istat   ! A return value from an mds function
     character(len=*) :: routine    ! Mds routine which returns the istat
     character(len=*) :: name       ! Mds node being loaded
     
     MdsNoError = .not. MdsError (istat, routine, name)

     End function MdsNoError 



     ! Routine to read string data from a node in mdsplus.
     ! We assume that we are connected to the server and that the proper
     ! tree has been opened.
     ! Written: 7/19/02 by rjg
     Subroutine Read_Mds_string (node, string, aok)

     Implicit none

     Character(len=*), intent(in) :: node ! Node name where desired
                                          ! data reside
                                          ! For example:
                                          ! '\top.nb30l:gas'
     Character(len=*), intent(out) :: string  ! The string read from mds+
     logical, intent(out) :: aok          ! TRUE if everything is okay

     Integer*4 :: istat              ! error code - even bad, odd good.
     Integer*4 :: idscr              ! Descriptor for the fortran api
     integer :: length  ! Length of string returned??

     ! Try to read back the string
     idscr = descr(IDTYPE_CSTRING, string, 0, len(string))
     istat = MDSVALUE(trim(node)//CHAR(0), idscr, 0, length)
     aok = .not. MdsError (istat,'Read_Mds_String',trim(node))

     end Subroutine Read_Mds_string 


     ! Routine to read real scalars from a node in mdsplus.
     ! We assume that we are connected to the server and that the proper
     ! tree has been opened.
     ! Written: 7/19/02 by rjg
     Subroutine Read_Mds_Real_Scalar (node, scalar, aok)

     Implicit none

     Character(len=*), intent(in) :: node ! Node name where desired
                                          ! data reside
                                          ! '\top.nb30l:nbvac_scalar'
     real, intent(out) :: scalar          ! The scalar read from mds+
     logical, intent(out) :: aok          ! TRUE if everything is okay

     Integer*4 :: istat              ! error code - even bad, odd good.
     Integer*4 :: idscr              ! Descriptor for the fortran api
     integer :: length  ! Length of something returned??


     ! Try to read back the string
     idscr = descr(IDTYPE_FLOAT, scalar, 0)
     istat = MDSVALUE(trim(node)//CHAR(0), idscr, 0, length)
     aok = .not. MdsError (istat,'Read_Mds_Real_Scalar',trim(node))

     end Subroutine Read_Mds_Real_Scalar


      subroutine mds_ld_error (ier,cfunc) 
      implicit none
      integer ier
      character(len=*) cfunc
      if (ier .ne. 0) write (*,*) ' Error on MDS load for ', cfunc
      end subroutine mds_ld_error

 
      ! Write ntimes pts of real data into MDSplus for node name.
      ! Written: 1/23/03 by rjg
      ! Revised: 03/11/08 by rjg to null out the node if ntimes is 0
      logical function Write_Mds_Real_Arr (name, ntimes, data, cunits_in)

      Implicit none

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

      
      if (ntimes .le. 0) then
          Write_Mds_Real_Arr = write_mds_null (name)
          return
      endif


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

      Character(len=*) :: name        ! input:  name of mds variable, e.g. 'ALP_ERR'
      Integer      ntimes             ! input:  number of data point.
      Real         data(ntimes)       ! input:  data array
      Character(len=*) cunits
      Integer      ier                ! output: error code - even bad, odd good.

      Integer      istat
      Integer      i
      Integer      m
      integer      ilen


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



!---------------------------------------------------------------------
! Write ntimes pts of data into MDSplus for node name.
! Do no use multiplier
! Revised: 03/11/08 by rjg to null out data if ntimes = 0

      Subroutine Dump_Mds_Rarray_No_Mult (name, ntimes, data, cunits, ier)
      use mdslib
      Implicit none

      Character(len=*) :: name   ! input:  name of mds variable, e.g. 'ALP_ERR'
      Integer      ntimes             ! input:  number of data point.
      Real         data(ntimes)       ! input:  data array
      Character(len=*) cunits
      Integer      ier             ! output: error code - even bad, odd good.

      Integer      istat
      Integer      i
      Integer      m
      integer      ilen

      logical :: lok


      ! If there are no valid times, null out the data in mdsplus
      if (ntimes .le. 0) then
         lok = write_mds_null (name)
         if (lok) then 
            ier = 0
         else
            ier = -20
         endif
         return
      endif


      m = len_trim(name)
      name = name(1:m) // char(0)
      ilen = len_trim(cunits)
      cunits = cunits(1:ilen) // char(0) 

      istat = mdsput (name, "BUILD_WITH_UNITS($,$)"//CHAR(0), &
     &              descr(IDTYPE_FLOAT, data, ntimes, 0),                 &
     &              descr(IDTYPE_CSTRING, Cunits, 0, ilen),0)
      if (MdsError (istat,'mdsput',name)) then
         ier = 1   
         call mds_ld_error (ier,name) 
      endif

      End Subroutine Dump_Mds_Rarray_No_Mult




! Write ntimes pts of real*8 data into MDSplus for node name.
! Do not use multiplier
!
! Revised: 12/01/06 by rjg to remove the use of r*16, which was needed
!                   for VMS code only (which had G-float and D-float)
! Revised: 03/11/08 by rjg to null out data if ntimes = 0
! Revised: 03/11/08 by rjg to change name and move into this file
      Subroutine Write_Mds_R8array_No_Mult (name, ntimes, data_r8, cunits, ier)
      use mdslib
      Implicit none

      Character(len=*), intent(inout) :: name   ! name of mds variable, e.g. 'ALP_ERR'
      Integer, intent(in) :: ntimes          ! number of data point.
      Real*8, intent(in) ::  data_r8(ntimes)    ! data array
      Character(len=*), intent(inout) :: cunits 
      Integer, intent(out) ::    ier         ! error code - even bad, odd good.


      Integer      istat
      Integer      i
      Integer      m
      integer      ilen
      integer :: if
      logical :: lok


      ! If there are no valid times, null out the data in mdsplus
      if (ntimes .le. 0) then
         lok = write_mds_null (name)
         if (lok) then 
            ier = 0
         else
            ier = -20
         endif
         return
      endif

      m = len_trim(name)
      name = name(1:m) // char(0)
      ilen = len_trim(cunits)
      cunits = cunits(1:ilen) // char(0) 

      istat = mdsput (name, "BUILD_WITH_UNITS($,$)"//CHAR(0), &
     &              descr(IDTYPE_DOUBLE, data_r8, ntimes, 0),                 &
     &              descr(IDTYPE_CSTRING, Cunits, 0, ilen),0)
      if (MdsError (istat,'mdsput',name)) then
         ier = 1   
         call mds_ld_error (ier,name) 
      endif

      End Subroutine Write_Mds_R8array_No_Mult 



      ! Revised: 03/03/05 by rjg for better declaration of cunits
      !---------------------------------------------------
      Subroutine Write_Mds_Real_Scalar (name, datum, cunits_in, ier)

      Implicit none

      Character(len=*) :: name        ! input:  name of mds variable, e.g. 'ALP_ERR'
      Real         datum              ! input:  datum (scalar)
      Character(len=*), intent(in) :: cunits_in      ! input, 
      Integer      ier                ! output: error code - even bad, odd good.
      Integer      istat
      integer      ilen
      character(len(cunits_in)) :: cunits

      cunits = cunits_in

      ilen = len_trim (cunits)
      cunits = trim(cunits) // char(0)
 
      istat = MdsPut(name//Char(0),                                        &
     &               "BUILD_WITH_UNITS($,$)"//Char(0),                    &
     &               descr(IDTYPE_FLOAT,  datum, 0),                       &
     &               descr(IDTYPE_CSTRING, cunits, 0, ilen),0)


      ier = 0
      If (MdsError(istat,'MdsPut',name)) ier = 1
      
      End Subroutine Write_Mds_Real_Scalar 
  

!---------------------------------------------------------------------
! Write ntimes pts of data into MDSplus for node name.

      ! Disable until I write out units
 
      Subroutine Write_Mds_RScalar_Disabled (name, datum, ier)

      Implicit none


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



      ier = 0
      If (MdsError(istat,'MdsPut',name)) ier = 1

      End Subroutine Write_Mds_Int_Scalar
  


      !---------------------------------------------------------------------
      ! Write a character variable into MDSplus for node name.
      ! Revised: 12/01/06 by rjg to make cdatum variable length
      Subroutine Write_Mds_chars (name, cdatum, ier)
      Implicit none

      Character(len=*) :: name        ! input:  name of mds variable, e.g. 'ALP_ERR'
      Character(len=*) :: cdatum      ! input:  datum (scalar character string)
      Integer      ier                ! output: error code - 0 good,  1 bad.
      Integer      istat              ! error code - even bad, odd good.


      istat = MdsPut(name//Char(0),                                       &
     &               "$" // Char(0),                                      &
     &               descr(IDTYPE_CSTRING, cdatum, 0, len(cdatum)), 0)


      ier = 0
      If (MdsError(istat,'MdsPut',name)) ier = 1

      End Subroutine Write_Mds_chars
  

      !---------------------------------------------------------------------
      ! Write a null into MDSplus for node name in order to disable the node.
      ! See http://www.mdsplus.org/old/mdslib/mdslib.html for info on mdsvalue2
      !
      ! Note that as of 12/20/12, code was crashing on venus with the use of
      ! tcl and mdsvalue2.  This was because name was being clobbered by 
      ! something, possibly the tcl command.  With Sean Flanagan's help,
      ! I redid the mdsplus calls by using an mdsput to write the null.
      ! This fixed the  problem.
      !
      ! Written: 03/10/08 by rjg
      ! Revised: 03/11/08 by rjg to use tcl and to use mdsvalue2
      ! Revised: 12/21/12 by rjg to use put rather than tcl and to use mdsvalue2;
      !                          
      logical function Write_Mds_Null (name)
      Implicit none

      Character(len=*) :: name    ! input:  name of mds variable, e.g. 'ROTC'

      Integer :: istat                 ! error code - even bad, odd good.


      ! Write a null to the node specified name
      istat = MdsPut (trim(name)//char(0), char(0), 0 )
 
      ! Want this to be true if istat is odd
      write_mds_null = (mod(istat,2) .ne. 0)

      if (.not. write_mds_null)   then 
          print *,' Write_Mds_Null: error from write_mds_null'
          print *,' istat = ', istat
      endif

      End function Write_Mds_Null
  

      !---------------------------------------------------------------------
      ! Write a null into MDSplus for node name in order to disable the node.
      ! Written: 03/11/08 by rjg
      logical function Write_Mds_Null_works (name)
      Implicit none

      Character(len=*) :: name    ! input:  name of mds variable, e.g. 'ROTC'

      character(len=50) :: tclstring   ! TCL command sent to mdsplus
      Integer :: istat                 ! error code - even bad, odd good.
      integer :: desc                  ! Descriptor needed for fortran API

      ! Want to set value stored in name to undefined.  In IDL, I would just
      ! do an mdsput with '' as the value.  In fortran API, the only way I 
      ! have found to do this is by using a TCL command.
      ! If name is ROTC, want tclstring to look like:
      !     tcl('put ROTC ""') 

      tclstring = "tcl('put " // trim(name) // ' ""' // "')"

      desc = descr(IDTYPE_CSTRING, trim(tclstring), 0, len_trim(tclstring) )

      istat = mdsvalue (trim(tclstring)//char(0), desc, 0, 1)

      ! Want this to be true if istat is odd
      write_mds_null_works = (mod(istat,2) .ne. 0)

      if (.not. write_mds_null_works)                                    &
          print *,' Write_Mds_Null: error from write_mds_null'

      End function Write_Mds_Null_works
  


!---------------------------------------------------------------------
! Write a real data value into MDSplus for node name.

      ! Disable until I write out units

      Subroutine Write_Mds_Real_Disabled (name, datum, ier)

      Implicit none

      Character(len=*) :: name        ! input:  name of mds variable, e.g. 'ALP_ERR'
      Real         datum              ! input:  datum (scalar)
      Integer      ier                ! output: error code - even bad, odd good.
      Integer      istat

      istat = MdsPut(name//Char(0),                                      &
     &               "$" // Char(0),                                     &
     &               descr(IDTYPE_FLOAT, datum, 0), 0)

      ier = 0
      If (MdsError(istat,'MdsPut',name)) ier = 1


      Return
      End Subroutine Write_Mds_Real_Disabled
  

      ! Write an array of integers to mdsplus
      ! Written: circa 2000 by Terpstra
      ! Revised: 1/23/03 by rjg to write units
      ! Revised: 03/11/08 by rjg to null out node if ntimes = 0
      logical function Write_Mds_Int_Arr (name, ntimes, idata, cunits)

      Implicit none

      Character(len=*), intent(in) :: name  ! name of mds node, e.g. 'ALP_ERR'
      Integer, intent(in) :: ntimes
      Character(len=*), intent(in) :: cunits
      Integer, dimension(ntimes), intent(in) :: idata   ! The data

      Integer      istat
      integer ::  desc1, desc2
      integer :: ilen


      ! If no data, null out the node      
      if (ntimes .le. 0) then
          Write_Mds_Int_Arr = write_mds_null (name)
          return
      endif

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



      ! Write ntimes pts of character data into MDSplus for node name.
      ! New and not working yet
      ! Written: 1/24/03 by rjg
      ! Revised: 1/30/03 by rjg to use cdata_loc
      ! Revised: 10/27/06 by rjg to work with Portland 
      ! Revised: 03/11/08 by rjg to null out data if ntimes = 0
      !                      I cannot get this to work so I am abandoning it
      ! Revised: 07/10/09 by novi to fix call to invalid array index "0".
      logical function Write_Mds_Char_Arr (name, ntimes, cdata)

      Implicit none

      Character(len=*), intent(in) :: name  ! mds nodename, e.g. 'FIT_MODEL'
      Integer, intent(in) :: ntimes
      Character(len=*), dimension(ntimes), intent(in) :: cdata  ! Character data

      ! Declare array of character strings so that we can properly size
      ! the goodies sent to mdsplus.  Make sure each element has at least
      ! one character to avoid code crashes.
      character(len=max(1,maxval(len_trim(cdata)))), dimension(ntimes) :: mydata


      ! Want the elements in cdata_loc to be no longer than the longest 
      ! element in cdata, but at least one character long (in order to
      ! avoid code crashes.

      Integer :: istat
      integer :: desc
      integer :: length

      mydata = cdata
       
      length = len(mydata(1))
 
      !  Http://www.mdsplus.org/old/mdslib/mdslib.html#vararg
      !  Creating a descriptor for a string array (length = space allocated 
      !  for one member of the array, ntimes = number of elements in array): 
      desc = descr(IDTYPE_CSTRING, mydata, ntimes, 0, length)

      istat = MdsPut(trim(name) // Char(0),                              &
     &               "$" // Char(0),                                     &
     &               desc, 0)

      Write_Mds_Char_Arr = MdsNoError(istat,'MdsPut',name)
      if (.not. Write_Mds_Char_Arr)                           &
         write (*,*) 'Write_Mds_Char_Arr: Could not write to ', name

      End function Write_Mds_Char_Arr


  end Module mds_rg
!
