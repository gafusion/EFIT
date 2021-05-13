!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          block data routine to hold experiment-dependent data    **
!**          statements for variables that are in common blocks.     **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          10/01/97..........created by Q.Peng, data were          **
!**                            in block data efit_btdata originally  ** 
!**                                                                  **
!**********************************************************************
!
      block data exp_bdata
      use commonblocks,only: byringr,byringz
      include 'eparmdud129.inc'
!vas      include 'modules2.inc' !to avoid the clash with allocat..
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!      include 'ecomdu1.f90'
!      include 'ecomdu2.f90'
!
!      common/adp/byringr(nh2),byringz(nh2),ringr(6),ringz(6),ringap
      common/adp/ringr(6),ringz(6),ringap
!
!---D3D-----------------------D3D----------------------------D3D-----
!-- the parameters set between the d3d comments are specifically   --
!-- for the D3D tokamak                                            --
!--------------------------------------------------------------------
!oct14      data fcid/1. ,2. ,3. ,4. ,5. ,6. ,7. ,8. ,9. , &
!oct14                10.,11.,12.,13.,14.,15.,16.,17.,18./
!vasorg      data relip/1.68/,zelip/0.00/,aelip/0.60/,eelip/1.2/
!vasorg      data rcentr/1.6955/,rzero/1.6955/,emaxis/1.3/
!oct14      data rcentr/1.6955/
!vasorg      data aaslop/0.6/,drslop/0.003/,nslref/1/
!oct14      data nslref/1/
!oct14      data chordv/1.486,1.945,2.098/
!oct14      data chordr/.00, 0.1524/,zcentr/0./
!oct14      data rmajts/1.94/,zlibim/-0.127/
!vasorg      data rubaf/1.372/,rudom/1.0420/,rlbaf/1.6810/
!vasorg      data zubaf/1.310/,zudom/1.1624/,zlbaf/-1.339/
!oct15      data ecurrt(3)/-1.e10/,ecurrt(4)/-1.e10/, &
!oct15           ecurrt(5)/-1.e10/,ecurrt(6)/-1.e10/, &
!oct15           ecurrt(1)/0.0/,ecurrt(2)/0.0/
      data ringr/1.766,1.680,1.674,2*1.671,1.681/
      data ringz/2*-1.255,-1.258,-1.264,-1.327,-1.335/
!vasorg      data do_spline_fit/.true./
!vasorg      data use_alternate_pointnames/0/
!vasorg      data alternate_pointname_file/'/link/efit/pcsnames.dat'/
!---------------------------------------------------------------------
!--  3.5 degree correction, per Ferron, Lazarus, Snider   08/02/90  --
!---------------------------------------------------------------------
!oct14      data xangle/    120.,124.,128.,132.,136.,140.,144.,148., &
!oct14                     152.,156.,160.,164.,168.,172.,176.,180.,180., &
!oct14                      184.,188.,192.,196.,200.,204.,208.,212.,216., &
!oct14                      220.,224.,228.,232.,236.,240., &
!oct14                   294.5,290.5,286.5,282.5,278.5,274.5,270.5,266.5, &
!oct14                   262.5,258.5,254.5,250.5,246.5,242.5,238.5,234.5, &
!oct14                   234.5,230.5,226.5,222.5,218.5,214.5,210.5,206.5, &
!oct14                   202.5,198.5,194.5,190.5,186.5,182.5,178.5,174.5, &
!oct14                   0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./
!oct14      data zxray/-10.7,-10.7,-10.7,-10.7,-10.7,-10.7, &
!oct14                    -10.7,-10.7,-10.7,-10.7,-10.7,-10.7,-10.7,-10.7, &
!oct14                    -10.7,-10.7,-14.7,-14.7,-14.7,-14.7,-14.7,-14.7, &
!oct14                    -14.7,-14.7,-14.7,-14.7,-14.7,-14.7,-14.7,-14.7, &
!oct14                    -14.7,-14.7, &
!oct14                    130.1,130.1,130.1,130.1,130.1,130.1,130.1,130.1, &
!oct14                    130.1,130.1,130.1,130.1,130.1,130.1,130.1,130.1, &
!oct14                    132.6,132.6,132.6,132.6,132.6,132.6,132.6,132.6, &
!oct14                    132.6,132.6,132.6,132.6,132.6,132.6,132.6,132.6, &
!oct14                   0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./
!oct14      data rxray/248.9,248.9,248.9,248.9,248.9,248.9,248.9, &
!oct14                248.9,248.9,248.9,248.9,248.9,248.9,248.9,248.9,248.9, &
!oct14                248.9,248.9,248.9,248.9,248.9,248.9,248.9,248.9, &
!oct14                248.9,248.9,248.9,248.9,248.9,248.9,248.9,248.9, &
!oct14                197.6,197.6,197.6,197.6,197.6,197.6,197.6,197.6, &
!oct14                197.6,197.6,197.6,197.6,197.6,197.6,197.6,197.6, &
!oct14                194.5,194.5,194.5,194.5,194.5,194.5,194.5,194.5, &
!oct14                194.5,194.5,194.5,194.5,194.5,194.5,194.5,194.5, &
!oct14                   0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./
!---D3D-----------------------D3D----------------------------D3D-----
!-- end of D3D parameters setting                                  --
!--------------------------------------------------------------------
!vasorg      data omega/1.0/,errorq/1.0e-03/,jjmax/1/
!oct14      data ksxr0/10*0/,ksxr2/10*0/,idosxr/1/
      end

!------------------------------------------------------------------------
!--
!--   read in limiter data from limiter file
!--
!--   type 0 - Calculate xltype(_180) only
!--        1 - Calculate xltype and read in limiter data
!--   xltype - Inward distance of 0 degree outside limiter in cm.
!--            Initialization is passed in.
!--   xltype_180 -
!--            Inward distance of 180 degree outside limiter in cm.
!--            Initialization is passed in.
!--
!--   10/02/97 created
!--
!------------------------------------------------------------------------
      subroutine getlim(type,xltype,xltype_180)
      include 'eparmdud129.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      integer type,ilimshot
      character*100 filin
      namelist/lim/xlim,ylim,limitr
      data lfile/36/
!
      if ((limitr.gt.0).and.(type.ne.0)) go to 50240
!
! --- read in limiter data
!
      limid=iabs(limitr)

      filin=input_dir(1:lindir)//'lim.dat'
      open(unit=lfile,access='sequential', &
      status='old',file=filin,err=120)

  100 read (lfile,*, end=120) ilimshot, limitr
 
      do i=1,limitr
        read(lfile,5010) xlim(i),ylim(i)
      enddo
      if (ishot.ge.ilimshot) then
        goto 120
      else
        goto 100
      endif
  120 close(unit=lfile)
  
!----------------------------------------------------------------
!--  Fixed poloidal limiter for shot > 88400                   --
!----------------------------------------------------------------
      if (ishot.gt.88400) go to 50240
      dell=0.01_dp*xltype
      limup=0
      limbot=0
      do i=1,limitr
        if (ishot.le.52360) then
        if ((ylim(i).eq.0.3260_dp).and.(xlim(i).gt.rcentr)) limup=i
        if ((ylim(i).eq.-0.3260_dp).and.(xlim(i).gt.rcentr)) limbot=i
        endif
        if (ishot.gt.52360.and.ishot.le.57590) then
        if ((ylim(i).eq.0.2000_dp).and.(xlim(i).gt.rcentr)) limup=i
        if ((ylim(i).eq.-0.2000_dp).and.(xlim(i).gt.rcentr)) limbot=i
        endif
        if (ishot.gt.57590) then
        if ((ylim(i).eq.0.3655_dp).and.(xlim(i).gt.rcentr)) limup=i
        if ((ylim(i).eq.-0.3655_dp).and.(xlim(i).gt.rcentr)) limbot=i
        endif
      enddo
      limup0=limup
      limbot0=limbot
      if ((limup.gt.0).and.(limbot.gt.0)) then
          do i=limup,limbot
            xlim(i)=xlim(i)-dell
          enddo   
          do i=limitr+5,limitr+10
            xlim(i)=xlim(i)-dell
          enddo    
      endif
!
! --- calculate the length of filename filimt
      lfilimt=0
      do i=1,len(filimt)
         if (filimt(i:i).ne.' ') lfilimt=lfilimt+1
      enddo
!
      if (ishot.ge.54021) then
      open(unit=lfile,                       status='old', &
        file=filimt(1:lfilimt)//'_180'          )
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
!--  no 180 degree limiter for shot greater than 76756, only fixed   --
!--  poloidal limiter                                                --
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
          do i=limup,limbot
            xlim_180(i)=xlim_180(i)-dell
          enddo
          do i=limitr_180+5,limitr_180+10
            xlim_180(i)=xlim_180(i)-dell
          enddo
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
50240 continue
      !SEK:  I have no idea this works if limitr < 0 so just doing this hack
      !TODO
      limitr=ABS(limitr)
      !SEK:  I have no idea this works
      limitr=limitr+1
      xlim(limitr)=xlim(1)
      ylim(limitr)=ylim(1)
      xlmin=xlim(1)
      xlmax=xlmin
      ylmin=ylim(1)
      ylmax=ylmin
      do i=2,limitr
        xlmin=min(xlmin,xlim(i))
        xlmax=max(xlmax,xlim(i))
        ylmin=min(ylmin,ylim(i))
        ylmax=max(ylmax,ylim(i))
      enddo
!---------------------------------------------------------------------------
!--  set up plotting limiter points if needed                             --
!---------------------------------------------------------------------------
      mimitr=limitr
      do i=1,limitr
       xmim(i)=xlim(i)
       ymim(i)=ylim(i)
      enddo
      if (iplim.ne.-1.or.ishot.le.76756) go to 83456
      open(unit=lfile,                       status='old', &
        file=filimt(1:lfilimt)//'_plot')
      read (lfile,*)    mimitr
      read (lfile,*)    (xmim(i),ymim(i),i=1,mimitr)
      close(unit=lfile)
      dell=0.01_dp*xltype
      limup=0
      limbot=0
      do 82206 i=1,mimitr
        if ((ymim(i).eq.0.3655_dp).and.(xmim(i).gt.rcentr)) limup=i
        if ((ymim(i).eq.-0.3655_dp).and.(xmim(i).gt.rcentr)) limbot=i
82206 continue
      if ((limup.le.0).or.(limbot.le.0)) go to 82211
      do 82208 i=limup,limbot
        xmim(i)=xmim(i)-dell
82208 continue
82211 continue
!
      open(unit=lfile,                       status='old', &
        file=filimt(1:lfilimt)//'_180'//'_plot' )
      read (lfile,*)    mimitr_180
      read (lfile,*)    (xmim_180(i),ymim_180(i),i=1,mimitr_180)
      close(unit=lfile)
      if (limup*limbot.gt.0) then
      rgraphite=2.3756_dp
      do 82220 i=limup,limbot
        xmim(i)=min(xmim(i),xmim_180(i),rgraphite)
82220 continue
      endif
      mimitr=mimitr+1
      xmim(mimitr)=xmim(1)
      ymim(mimitr)=ymim(1)
83456 continue
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
      end
!------------------------------------------------------------------------
!--
!--   read in SXR detectors geometry from a file if needed
!--
!--   10/01/97 created
!--   05/12/98 Q.P. added Toroidal X-ray. Data provided by R.Snider.
!--
!------------------------------------------------------------------------
      subroutine getsxr(ishot,ixray)
      include 'eparmdud129.inc'
      use var_cxray
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!vas      common/cxray/rxray(nangle+ntangle),zxray(nangle+ntangle), &
!vas                   xangle(nangle+ntangle) &
!vas                   ,ksxr0(10),ksxr2(10),idosxr
!------------------------------------------------------------------------
!--  read in new SXR detectors geometry for shot > 80744               --
!------------------------------------------------------------------------
      if ((ishot.ge.80744).and.(iand(ixray,1).ne.0)) then
        open(unit=80,status='old', &
           file=input_dir(1:lindir)//'sxr94.dat')
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
!--  read in new Toroidal x-ray geometry for shot > 91000 (since '97)  --
!------------------------------------------------------------------------
      if ((ishot.lt.91000).and.(iand(ixray,2).ne.0)) ixray = ixray-2
      if ((ishot.ge.91000).and.(iand(ixray,2).ne.0)) then
        open(unit=80,status='old', &
           file=input_dir(1:lindir)//'sxrt97.dat')
        read (80,*) rxrayi,zxrayi
        do j=nangle+1,nangle+ntangle
           rxray(j)=rxrayi
           zxray(j)=zxrayi
           read (80,*) idumdum,xangle(j)
        enddo
      endif
!      
      return
      end
