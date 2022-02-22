!**********************************************************************
!>
!!    This subroutine does...
!!
!!    @param let : is letter designation of file   
!!
!!    @param ishot : is shot number   
!!
!!    @param itime : is time in miliseconds 
!!
!!    @param fname : is file name 
!!
!**********************************************************************
      subroutine getfnmd(let,ishot,itime,fname)
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      character*(*) let
      character*(*) fname
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
      end subroutine getfnmd

!**********************************************************************
!>
!!    This subroutine does...
!!
!!    @param let : is letter designation of file  (2 chracters)   
!!
!!    @param ishot : is shot number   
!!
!!    @param itime : is time in miliseconds 
!!
!!    @param fname : is file name 
!!
!**********************************************************************
      subroutine getfnm2(let,ishot,itime,fname)
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      character*(*) let
      character*(*) fname
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
      end subroutine getfnm2

!**********************************************************************
!>
!!    This subroutine does ...
!!    
!!
!!    @param itimeu : time in microseconds 
!!
!!    @param let : letter designation of file 
!!
!!    @param ishot : shot number
!!
!!    @param itime : time in miliseconds
!!
!!    @param fname : file name 
!!
!**********************************************************************
      subroutine getfnmu(itimeu,let,ishot,itime,fname)
      character*(*) let
      character*(*) fname
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
      end subroutine getfnmu

!**********************************************************************
!>
!!    This subroutine does ...
!!    
!!
!!    @param itimeu : time in microseconds 
!!
!!    @param let : letter designation of file 
!!
!!    @param ishot : shot number
!!
!!    @param itime : time in miliseconds
!!
!!    @param fname : file name 
!!
!**********************************************************************
      subroutine getfnmu2(itimeu,let,ishot,itime,fname)
      character*(*) let
      character*(*) fname
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
      end subroutine getfnmu2

!**********************************************************************
!>
!!    This subroutine does ...
!!    
!!
!!    @param let : letter designation of file 
!!
!!    @param ishot : shot number
!!
!!    @param itime : time in miliseconds
!!
!!    @param istore : directory flag\n
!!                    =0  Default directory
!!                    =1  Central directory 
!!
!!    @param fname : file name 
!!
!**********************************************************************      
      subroutine setfnme(let,ishot,itime,istore,fname)
      use exvars
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      character*(*) let
      character*(*) fname
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
      end subroutine setfnme
