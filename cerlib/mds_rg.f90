
!
!  Module with some common routines used to deal with mds database
!
  Module mds_rg

  implicit none
!  include 'hpalias.inc'

     logical, save :: connected_to_mds = .false.   ! True means we are
                                                   ! connected to mds server

  contains


     ! Subroutine to connect to mds server
     ! Written: 4/20/00 by rjg
     ! Revised: 1/12/01 by rjg to check if we are already connected
     ! Revised: 7/08/02 by rjg to add optional argument to return status
     ! Revised: 7/12/02 by rjg to solve logic problem introduced on 7/8/02
     subroutine wake_mds (lmds_awake)
                                                                       
     implicit none
     Include 'mdslib.inc' 

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
                                                                       
     implicit none
     Include 'mdslib.inc' 
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
                                                                       
     implicit none
     Include 'mdslib.inc' 
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
     


!      ! Written: 4/13/00 by rjg
!      logical function MdsError (istat,routine,name)
!      implicit none

!      integer, intent(in) :: istat   ! A return value from an mds function
!      character(len=*) :: routine    ! Mds routine which returns the istat
!      character(len=*) :: name       ! Mds node being loaded
     
!      ! If istat is an even number, we have an error condition
!      MdsError = .not. good_status (istat)
!      If (MdsError) then
!         print *,' Error from routine ', trim(routine), ' for ', trim(name)
!         print *,' istat was ', istat
!      Endif

!      End function MdsError 


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
     If (MdsError) then
        print *,' Error from routine ', trim(routine), ' for ', trim(name)
        print *,' istat was ', istat 
     Endif
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

     Include 'mdslib.inc'

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

     Include 'mdslib.inc'

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
 

  end Module mds_rg

