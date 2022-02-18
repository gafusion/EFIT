      subroutine getfnmd(let,ishot,itime,fname)
!----------------------------------------------------------------------
!--   let is letter designation of file                              --
!--   ishot is shot number                                           --
!--   itime is time in miliseconds                                   --
!--   fname is file name                                             --
!----------------------------------------------------------------------
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      character*(*) let
      character*(*) fname
      character*5 sec_msec
!
      length=len(fname)
      do 100 i=1,length
       fname(i:i)=' '
  100 continue
      iitime=itime
      if (iitime.lt.0) iitime=90000-iitime
      write(sec_msec,1000) iitime
 1000 format(i5.5)
      write(fname,1020) let,ishot,sec_msec
 1020 format(a1,i6.6,'.',a5)
      return
      end
      subroutine getfnm2(let,ishot,itime,fname)
!----------------------------------------------------------------------
!--   let is letter designation of file (2 chracters)                --
!--   ishot is shot number                                           --
!--   itime is time in miliseconds                                   --
!--   fname is file name                                             --
!----------------------------------------------------------------------
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      character*(*) let
      character*(*) fname
      character*5 sec_msec
!
      length=len(fname)
      do 100 i=1,length
       fname(i:i)=' '
  100 continue
      iitime=itime
      if (iitime.lt.0) iitime=90000-iitime
      write(sec_msec,1000) iitime
 1000 format(i5.5)
      write(fname,1020) let,ishot,sec_msec
 1020 format(a2,i6.6,'.',a5)
      return
      end
      subroutine getfnmu(itimeu,let,ishot,itime,fname)
!----------------------------------------------------------------------
!--   itimeu is time in microseconds                                 --
!--   let is letter designation of file                              --
!--   ishot is shot number                                           --
!--   itime is time in miliseconds                                   --
!--   fname is file name                                             --
!----------------------------------------------------------------------
      character*(*) let
      character*(*) fname
      character*5 sec_msec
!
      length=len(fname)
      do 100 i=1,length
       fname(i:i)=' '
  100 continue
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
      end

      subroutine getfnmu_pefit(itimeu,let,ishot,itime,fname)
!----------------------------------------------------------------------
!--   itimeu is time in microseconds                                 --
!--   let is letter designation of file                              --
!--   ishot is shot number                                           --
!--   itime is time in miliseconds                                   --
!--   fname is file name                                             --
!----------------------------------------------------------------------
      character*(*) let
      character*(*) fname
      character*5 sec_msec
!
      length=len(fname)
      do 100 i=1,length
       fname(i:i)=' '
  100 continue
      iitime=itime
      if (iitime.lt.0) iitime=90000-iitime
      write(sec_msec,1000) iitime
 1000 format(i5.5)
      if (itimeu.le.0) then
        write(fname,1020) 'pefitkfile/',let,ishot,sec_msec
      else
        write(fname,1030) 'pefitkfile/',let,ishot,sec_msec,itimeu
      endif
 1020 format(a11,a1,i6.6,'.',a5)
 1030 format(a11,a1,i6.6,'.',a5,'_',i3.3)
      return
      end

      subroutine setfnme(let,ishot,itime,istore,fname)
!----------------------------------------------------------------------
!--   let is letter designation of file                              --
!--   ishot is shot number                                           --
!--   itime is time in miliseconds                                   --
!--   istore	directory flag					     --
!--		=0	Default directory			     --
!--		=1	Central directory			     --
!--   fname is file name                                             --
!----------------------------------------------------------------------
      use exvars
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      character*(*) let
      character*(*) fname
      character*5 sec_msec
!
      length=len(fname)
      do 100 i=1,length
       fname(i:i)=' '
  100 continue
      iitime=itime
      if (iitime.lt.0) iitime=90000-iitime
      write(sec_msec,1000) iitime
 1000 format(i5.5)
      write(fname,1020) let,ishot,sec_msec
 1020 format(a1,i6.6,'.',a5)
      if (istore .eq. 1) fname = store_dir(1:lstdir)//fname
      return
      end
!
!   This routine is required if the CVS revision numbers are to 
!   survive an optimization.
!
!
!   $Date: 2008/11/18 22:47:04 $ $Author: radhakri $
!
      subroutine getfnmx_rev(i)
      CHARACTER*100 opt
      character*10 s 
      if( i .eq. 0) s =   &
      '@(#)$RCSfile: getfnmdud.f90,v $ $Revision: 1.1.2.2 $\000'
      return
      end
