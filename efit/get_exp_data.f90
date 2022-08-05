!**********************************************************************
!>
!!    read in limiter data from limiter file
!!    
!!
!!    @param limmode : 0 - Calculate xltype(_180) only \n
!!           1 - Calculate xltype and read in limiter data
!!
!!    @param xltype : Inward distance of 0 degree outside limiter in cm.
!!            Initialization is passed in.
!!
!!    @param xltype_180 : Inward distance of 180 degree outside limiter in cm.
!!            Initialization is passed in.
!!
!!    @param shape_ext : Flag for whether profile_ext data is used. 
!!
!**********************************************************************
      subroutine getlim(limmode,xltype,xltype_180,shape_ext)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      integer*4 limmode,ilimshot,ioerr
      logical shape_ext
      character*100 filin
      data lfile/36/
!
      read_lim: if ((.not.shape_ext .and. limitr.le.0).or.(limmode.eq.0)) then
!
      if (abs(xltype*xltype_180).le.0.0001_dp.and.ishot.gt.52072) then
        call limpos(ishot,poslim,pos180,ier)
        if(ishot.ge.54992.and.ishot.le.55015) pos180=2.
        if (ier.eq.1) then
          if(abs(xltype).le.0.0001_dp) xltype=poslim
          if(abs(xltype_180).le.0.0001_dp) xltype_180=pos180
        else
          if(abs(xltype).le.0.0001_dp) xltype=0.0
          if(abs(xltype_180).le.0.0001_dp) xltype_180=0.0
        endif
      endif
      if(limmode.eq.0) return
!
! --- read in limiter data
!
      limid=iabs(limitr)

      filin=input_dir(1:lindir)//'lim.dat'
      open(unit=lfile,access='sequential',status='old',file=filin)
      ilimshot=999999
      do while (ishot.lt.ilimshot)
        read (lfile,*,iostat=ioerr) ilimshot, limitr
        if(ioerr.ne.0) exit
        do i=1,limitr
          read(lfile,5010) xlim(i),ylim(i)
        enddo
      enddo
      close(unit=lfile)
!----------------------------------------------------------------
!--   Fixed poloidal limiter for shot > 88400                  --
!----------------------------------------------------------------
      get_pol: if (ishot.le.88400) then
      dell=0.01_dp*xltype
      limup=0
      limbot=0
      do i=1,limitr
        if (ishot.le.52360) then
          if((ylim(i).eq.0.3260_dp).and.(xlim(i).gt.rcentr)) limup=i
          if((ylim(i).eq.-0.3260_dp).and.(xlim(i).gt.rcentr)) limbot=i
        endif
        if (ishot.gt.52360.and.ishot.le.57590) then
          if((ylim(i).eq.0.2000_dp).and.(xlim(i).gt.rcentr)) limup=i
          if((ylim(i).eq.-0.2000_dp).and.(xlim(i).gt.rcentr)) limbot=i
        endif
        if (ishot.gt.57590) then
          if((ylim(i).eq.0.3655_dp).and.(xlim(i).gt.rcentr)) limup=i
          if((ylim(i).eq.-0.3655_dp).and.(xlim(i).gt.rcentr)) limbot=i
        endif
      enddo
      limup0=limup
      limbot0=limbot
      if ((limup.gt.0).and.(limbot.gt.0)) then
        xlim(limup:limbot)=xlim(limup:limbot)-dell
        xlim(limitr+5:limitr+10)=xlim(limitr+5:limitr+10)-dell
      endif
!
! --- calculate the length of filename filimt
      lfilimt=0
      do i=1,len(filimt)
        if(filimt(i:i).ne.' ') lfilimt=lfilimt+1
      enddo
!
      if (ishot.ge.54021) then
        open(unit=lfile,status='old',file=filimt(1:lfilimt)//'_180')
        if (ishot.le.76756) then
          read (lfile,4980) limitr_180
          read (lfile,5000) (xlim_180(i),ylim_180(i),i=1,limitr_180) &
            ,(xlim_180(i),ylim_180(i),i=limitr_180+2,limitr_180+13)
        else
          read (lfile,*)    limitr_180
          read (lfile,*)    (xlim_180(i),ylim_180(i),i=1,limitr_180)
        endif
        close(unit=lfile)
!----------------------------------------------------------------------
!--     no 180 degree limiter for shot greater than 76756, only fixed--
!--     poloidal limiter                                             --
!----------------------------------------------------------------------
        dell=0.01_dp*xltype_180
        limup=0
        limbot=0
        do i=1,limitr_180
          if (ishot.le.57590) then
            if ((ylim_180(i).eq.0.3655_dp).and.(xlim_180(i).gt.rcentr)) &
              limup=i
            if ((ylim_180(i).eq.-0.3655_dp).and.(xlim_180(i).gt.rcentr)) &
              limbot=i
          endif
          if (ishot.gt.57590) then
            if ((ylim_180(i).eq.0.2000_dp).and.(xlim_180(i).gt.rcentr)) &
              limup=i
            if ((ylim_180(i).eq.-0.2000_dp).and.(xlim_180(i).gt.rcentr)) &
              limbot=i
          endif
        enddo
        if ((limup>0).and.(limbot>0)) then
          xlim_180(limup:limbot)=xlim_180(limup:limbot)-dell
          xlim_180(limitr_180+5:limitr_180+10) &
            =xlim_180(limitr_180+5:limitr_180+10)-dell
        endif
        if (ishot.ge.76757) then
          limup=limup0
          limbot=limbot0
          xltype_180=0.0
        endif
        if (limup0*limup*limbot*limbot0.gt.0) then
          rgraphite=2.3756_dp
          limup=min(limup,limup0)
          limbot=max(limbot,limbot0)
          do i=limup,limbot
            xlim(i)=min(xlim(i),xlim_180(i),rgraphite)
          enddo
        endif
      endif
      endif get_pol
      endif read_lim
      limitr=limitr+1
      xlim(limitr)=xlim(1)
      ylim(limitr)=ylim(1)
      xlmin=minval(xlim(1:limitr))
      xlmax=maxval(xlim(1:limitr))
      ylmin=minval(ylim(1:limitr))
      ylmax=maxval(ylim(1:limitr))
!---------------------------------------------------------------------------
!--   set up plotting limiter points if needed                            --
!---------------------------------------------------------------------------
      mimitr=limitr
      xmim(1:limitr)=xlim(1:limitr)
      ymim(1:limitr)=ylim(1:limitr)
      if(iplim.ne.-1.or.ishot.le.76756) return
      open(unit=lfile,status='old',file=filimt(1:lfilimt)//'_plot')
      read (lfile,*) mimitr
      read (lfile,*) (xmim(i),ymim(i),i=1,mimitr)
      close(unit=lfile)
      dell=0.01_dp*xltype
      limup=0
      limbot=0
      do i=1,mimitr
        if((ymim(i).eq.0.3655_dp).and.(xmim(i).gt.rcentr)) limup=i
        if((ymim(i).eq.-0.3655_dp).and.(xmim(i).gt.rcentr)) limbot=i
      enddo
      if ((limup.gt.0).and.(limbot.gt.0)) then
        xmim(limup:limbot)=xmim(limup:limbot)-dell
      endif
!
      open(unit=lfile,status='old', &
           file=filimt(1:lfilimt)//'_180'//'_plot' )
      read (lfile,*) mimitr_180
      read (lfile,*) (xmim_180(i),ymim_180(i),i=1,mimitr_180)
      close(unit=lfile)
      if (limup*limbot.gt.0) then
        rgraphite=2.3756_dp
        do i=limup,limbot
          xmim(i)=min(xmim(i),xmim_180(i),rgraphite)
        enddo
      endif
      mimitr=mimitr+1
      xmim(mimitr)=xmim(1)
      ymim(mimitr)=ymim(1)
!
 4980 format (i5)
 5000 format (2e12.6)
 5010 format (2E14.6)
 6530 format ('limd3d.0',i2)
 6540 format ('lim042386.0',i2)
 6542 format ('lim060686.0',i2)
 6544 format ('lim061686.0',i2)
 6546 format ('lim080786.0',i2)
 6548 format ('lim021887.0',i2)
 6650 format ('lim012588.0',i2)
 6652 format ('lim081390.0',i2)
 6654 format ('lim910417.0',i2)
 6656 format ('lim920203.0',i2)
 6658 format ('lim930411.0',i2)
 6660 format ('lim960321.0',i2)
 6662 format ('lim970407.0',i2)
 6664 format ('lim000113.0',i2)
 6666 format ('lim060530.0',i2)
 6668 format ('lim091109.0',i2)
 6670 format ('lim151110.0',i2)
 6672 format ('lim170120.0',i2)
 6674 format ('lim170120.',i3)
 6676 format ('lim200330.0',i2)
 6678 format ('lim200330.',i3)
!
      return
      end subroutine getlim

!**********************************************************************
!>
!!    read in SXR detectors geometry from a file if needed
!!    
!!
!!    @param ishot : shot number
!!
!!    @param ixray :
!!
!**********************************************************************
      subroutine getsxr(ishot,ixray)
      include 'eparm.inc'
      use error_control, only: errctrl_msg
      use var_cxray
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!------------------------------------------------------------------------
!--   read in new SXR detectors geometry for shot > 80744              --
!------------------------------------------------------------------------
      if ((ishot.ge.80744).and.(iand(ixray,1).ne.0)) then
        open(unit=80,status='old', &
             file=input_dir(1:lindir)//'sxr94.dat',iostat=ioerr)
        if (ioerr.ne.0) then
          call errctrl_msg('getsxr','could not open SXR file')
          stop
        endif
        do i=1,4
          ii=(i-1)*16
          read (80,*) rxrayi,zxrayi
          do j=ii+1,ii+16
            rxray(j)=rxrayi
            zxray(j)=zxrayi
            read (80,*) idumdum,xangle(j)
          enddo
        enddo
      endif
!------------------------------------------------------------------------
!--   read in new Toroidal x-ray geometry for shot > 91000 (since '97) --
!------------------------------------------------------------------------
      if((ishot.lt.91000).and.(iand(ixray,2).ne.0)) ixray = ixray-2
      if ((ishot.ge.91000).and.(iand(ixray,2).ne.0)) then
        open(unit=80,status='old', &
             file=input_dir(1:lindir)//'sxrt97.dat')
        if (ioerr.ne.0) then
          call errctrl_msg('getsxr','could not open SXR file')
          stop
        endif
        read (80,*) rxrayi,zxrayi
        do j=nangle+1,nangle+ntangle
          rxray(j)=rxrayi
          zxray(j)=zxrayi
          read (80,*) idumdum,xangle(j)
        enddo
      endif
!      
      return
      end subroutine getsxr
