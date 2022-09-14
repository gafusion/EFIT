!**********************************************************************
!>
!!    This subroutine sets the file names for eqdsk files
!!    
!!
!!    @param itimeu : time in microseconds 
!!
!!    @param let : letter designation of file 
!!
!!    @param ishot : shot number
!!
!!    @param itime : time in milliseconds
!!
!!    @param fname : file name 
!!
!**********************************************************************
      subroutine setfnmeq(itimeu,let,ishot,itime,fname)
      implicit none
      integer*4, intent(in) :: itimeu,ishot,itime
      character*(*), intent(in) :: let
      character*(*), intent(inout) ::  fname
      integer*4 i,iitime,length
      character*5 sec_msec
!
      length=len(fname)
      do i=1,length
        fname(i:i)=' '
      enddo
      iitime=itime
      if (iitime.lt.0) iitime=90000-iitime
      write(sec_msec,1000) iitime
 1000 format(i5.5)
      if (itimeu.le.0) then
        write(fname,1020) let,ishot,sec_msec
      else
        write(fname,1030) let,ishot,sec_msec,itimeu
      endif
 1020 format(a1,i6.6,'.',a5)
 1030 format(a1,i6.6,'.',a5,'_',i3.3)
      return
      end subroutine setfnmeq

!**********************************************************************
!>
!!    This subroutine sets the file names for diagnostic files
!!
!!    @param let : is letter designation of file   
!!
!!    @param ishot : is shot number   
!!
!!    @param itime : is time in milliseconds 
!!
!!    @param fname : is file name 
!!
!**********************************************************************
      subroutine setfnmd(let,ishot,itime,fname)
      implicit none
      integer*4, intent(in) :: ishot,itime
      character*(*), intent(in) :: let
      character*(*), intent(inout) ::  fname
      integer*4 i,iitime,length
      character*5 sec_msec
!
      length=len(fname)
      do i=1,length
        fname(i:i)=' '
      enddo
      iitime=itime
      if (iitime.lt.0) iitime=90000-iitime
      write(sec_msec,1000) iitime
 1000 format(i5.5)
      write(fname,1020) let,ishot,sec_msec
 1020 format(a1,i6.6,'.',a5)
      return
      end subroutine setfnmd

!**********************************************************************
!>
!!    This subroutine sets the file names for time history files
!!
!!    @param let : is letter designation of file  (2 chracters)   
!!
!!    @param ishot : is shot number   
!!
!!    @param itime : is time in milliseconds 
!!
!!    @param fname : is file name 
!!
!**********************************************************************
      subroutine setfnmt(let,ishot,itime,fname)
      implicit none
      integer*4, intent(in) :: ishot,itime
      character*(*), intent(in) :: let
      character*(*), intent(inout) ::  fname
      integer*4 i,iitime,length
      character*5 sec_msec
!
      length=len(fname)
      do i=1,length
        fname(i:i)=' '
      enddo
      iitime=itime
      if (iitime.lt.0) iitime=90000-iitime
      write(sec_msec,1000) iitime
 1000 format(i5.5)
      write(fname,1020) let,ishot,sec_msec
 1020 format(a2,i6.6,'.',a5)
      return
      end subroutine setfnmt

!**********************************************************************
!>
!!    This subroutine sets the file names for pl plotting files
!!    
!!
!!    @param itimeu : time in microseconds 
!!
!!    @param let : letter designation of file 
!!
!!    @param ishot : shot number
!!
!!    @param itime : time in milliseconds
!!
!!    @param fname : file name 
!!
!**********************************************************************
      subroutine setfnmpl(itimeu,let,ishot,itime,fname)
      implicit none
      integer*4, intent(in) :: itimeu,ishot,itime
      character*(*), intent(in) :: let
      character*(*), intent(inout) ::  fname
      integer*4 i,iitime,length
      character*5 sec_msec
!
      length=len(fname)
      do i=1,length
        fname(i:i)=' '
      enddo
      iitime=itime
      if (iitime.lt.0) iitime=90000-iitime
      write(sec_msec,1000) iitime
 1000 format(i5.5)
      if (itimeu.le.0) then
        write(fname,1020) let,ishot,sec_msec
      else
        write(fname,1030) let,ishot,sec_msec,itimeu
      endif
 1020 format(a2,i6.6,'.',a5)
 1030 format(a2,i6.6,'.',a5,'_',i3.3)
      return
      end subroutine setfnmpl

!**********************************************************************
!>
!!    This subroutine sets the file names for q plotting files
!!    
!!
!!    @param let : letter designation of file 
!!
!!    @param ishot : shot number
!!
!!    @param itime : time in milliseconds
!!
!!    @param istore : directory flag\n
!!                    =0  Default directory
!!                    =1  Central directory 
!!
!!    @param fname : file name 
!!
!**********************************************************************      
      subroutine setfnmq(let,ishot,itime,istore,fname)
      use extvars, only : store_dir,lstdir
      implicit none
      integer*4, intent(in) :: ishot,itime,istore
      character*(*), intent(in) :: let
      character*(*), intent(inout) ::  fname
      integer*4 i,iitime,length
      character*5 sec_msec
!
      length=len(fname)
      do i=1,length
        fname(i:i)=' '
      enddo
      iitime=itime
      if (iitime.lt.0) iitime=90000-iitime
      write(sec_msec,1000) iitime
 1000 format(i5.5)
      write(fname,1020) let,ishot,sec_msec
 1020 format(a1,i6.6,'.',a5)
      if (istore .eq. 1) fname = store_dir(1:lstdir)//fname
      return
      end subroutine setfnmq
