!**********************************************************************
!>
!!    set_init sets the initial condition for the plasma (either
!!      the normalized poloidal flux or plasma current) and other
!!      required variables to restart from an existing solution
!!
!!
!!    @param jtime : time index
!!
!**********************************************************************
      subroutine set_init(jtime)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      integer*4, intent(in) :: jtime
      integer*4 :: i,j,kk,ioerr,mw,mh
      real*8 :: delcur,sumi,erho
      character(14) :: sfile
      character(1000) :: line
      logical :: file_stat
      integer*8 :: ndout(1)
      real*8 :: alpa(nw)
!
      pcurrt=0.0
      if(ivacum.eq.1) return
!-----------------------------------------------------------------------
!--   icinit=1 uniform and elliptical flux surfaces
!--          2 parabolic and elliptical
!--         -2 use ESAVE.DAT at first time and then previous solution
!--         -3 use existing solution (geqdsk_ext or hdf5) for every time
!--         -4 use existing solution for first time and then previous
!--        -12 parabolic and elliptical for first time and then previous
!-----------------------------------------------------------------------
      if(jtime.eq.1) isicinit=icinit
      if (isicinit.lt.0 .and. isicinit.ne.-3) then
        if (jtime.gt.1) then
          icinit=isicinit
          return
        endif
        if (isicinit.eq.-12) then
          icinit=2
        elseif (isicinit.eq.-4) then
          icinit=-3
        endif
      endif
!
      select case (icinit)
      case (1)
        if (aelip.le.0.0) &
          aelip=min(0.50_dp, &
                minval(sqrt((xlim(1:limitr)-relip)**2 &
                          +((ylim(1:limitr)-zelip)/eelip)**2)))
        delcur=1.
        sumi=0.0
        do i=1,nw
          do j=1,nh
            kk=(i-1)*nh+j
            pcurrt(kk)=delcur*zero(kk)/rgrid(i)
            sumi=sumi+pcurrt(kk)
          enddo
        enddo
        cratio=ipmeas(jtime)/sumi
        pcurrt=pcurrt*cratio*zero
!
      case (2)
        zelip=zelipss
        if (zelip.gt.1.e5_dp) then
          ! DIII-D specific values
          zelip=1.447310_dp*(expmpi(jtime,37)-expmpi(jtime,43)) &
               +0.692055_dp*(expmpi(jtime,57)-expmpi(jtime,53)) &
               +0.728045_dp*(silopt(jtime,27)-silopt(jtime,37)) &
               +2.047150_dp*(silopt(jtime,2) -silopt(jtime,11))
          zelip=zelip*1.e6_dp/ipmeas(jtime)
          zbound=zelip
          eelip=1.5_dp
!----------------------------------------------------------------
!--       set zelip=0.0 if bad signals              96/06/24   --
!----------------------------------------------------------------
          if (size(fwtmp2).ge.37) then
            if(abs(fwtmp2(37)).le.1.0e-30_dp) zelip=0.0
          endif
          if (size(fwtmp2).ge.43) then
            if(abs(fwtmp2(43)).le.1.0e-30_dp) zelip=0.0
          endif
          if (size(fwtmp2).ge.57) then
            if(abs(fwtmp2(57)).le.1.0e-30_dp) zelip=0.0
          endif
          if (size(fwtmp2).ge.53) then
            if(abs(fwtmp2(53)).le.1.0e-30_dp) zelip=0.0
          endif
          if (size(fwtsi).ge.27) then
            if(abs(fwtsi(27)).le.1.0e-30_dp)  zelip=0.0
          endif
          if (size(fwtsi).ge.37) then
            if(abs(fwtsi(37)).le.1.0e-30_dp)  zelip=0.0
          endif
          if (size(fwtsi).ge.2) then
            if(abs(fwtsi(2)).le.1.0e-30_dp)   zelip=0.0
          endif
          if (size(fwtsi).ge.11) then
            if(abs(fwtsi(11)).le.1.0e-30_dp)  zelip=0.0
          endif
        endif
!
        do i=1,nw
          do j=1,nh
            kk=(i-1)*nh+j
            erho=sqrt((rgrid(i)-relip)**2+((zgrid(j)-zelip)/eelip)**2)
            xpsi(kk)=(erho/aelip)**2
          enddo
        enddo
!
      case (-2)
        if (nproc.gt.1) then
          WRITE(sfile,fmt='(i5.5)') rank
          sfile='esave'//TRIM(sfile)//'.dat'
        else
          sfile='esave.dat'
        endif
        open(unit=nsave,form='unformatted',file=sfile, &
             status='old',iostat=ioerr)
        if (ioerr.eq.0) then 
          read (nsave) mw,mh
          if ((mw.ne.nw.or.mh.ne.nh)) then
            call errctrl_msg('set_init','esave dimensions do not match')
            stop
          endif
          read (nsave) xpsi
          read (nsave) brsp
          read (nsave) www
          read (nsave) emaxis, rmaxis, fcentr
          close(unit=nsave)
        endif
!
      case (-3)
        if ((nw_ext.ne.nw.or.nh_ext.ne.nh)) then
          call errctrl_msg('set_init', &
                           'previous case had different dimensions')
          stop
        endif

        simag=simag_ext
        psibry=psibry_ext

        ! compute auxiliary quantities
        sidif=simag-psibry
        xmin=minval(rbdry_ext(1:nbdry_ext))
        xmax=maxval(rbdry_ext(1:nbdry_ext))
        ymin=minval(zbdry_ext(1:nbdry_ext))
        ymax=maxval(zbdry_ext(1:nbdry_ext))

        ! normalize psi
        do i=1,nw_ext
         do j=1,nh_ext
          kk=(i-1)*nh_ext+j
          if (icutfp.eq.0) then
            xpsi(kk)=1.1_dp
            if((rgrid(i).lt.xmin).or.(rgrid(i).gt.xmax)) cycle
            if((zgrid(j).lt.ymin).or.(zgrid(j).gt.ymax)) cycle
            xpsi(kk)=(simag-psirz_ext(kk))/sidif
          else
            call errctrl_msg('set_init','icutfp/=0 not set up')
            stop
          endif
         enddo
        enddo

        ! setup weights
        if (iweigh.le.0) then
          do kk=1,nwnh
            if ((xpsi(kk).lt.0.0).or.(xpsi(kk).gt.1.0)) then
              www(kk)=0.0
            else
              www(kk)=zero(kk)
            endif
          enddo
        else
          call errctrl_msg('set_init','iweight>0 not set up')
          stop
        endif

        ! need to use computed fcoil currents
        brsp(1:nfsum)=fcoil_ext

        ! fix signs
        ! note: plasma_ext does not appear to be useful for this...
        if (ipmeas(jtime).gt.0.0) then
          pprime_ext=-pprime_ext
          ffprim_ext=-ffprim_ext
        endif

        ! get p' and ff' coefficients
        ! note: this could be read in from g-files directly, but needs
        !       to be fit from imas files, so we treat both consistently
        if (icurrt.eq.2) then
          ! apparently fitting polynomials works equally well with
          ! splines...? (kppfnc=6 and kfffnc=6) - matches esave
          call fitpp(pprime_ext,nw_ext,alpa,kppcur)
          brsp(1+nfsum:nfsum+kppcur)=alpa(1:kppcur)*darea
          call fitfp(ffprim_ext,nw_ext,alpa,kffcur)
          brsp(1+nfsum+kppcur:nfsum+kppcur+kffcur)= &
            alpa(1:kffcur)*darea/twopi/tmu
        else
          call errctrl_msg('set_init','icurrt/=2 not set up')
          stop
        endif
      case default
        call errctrl_msg('set_init','icinit value not recognized')
        stop
      end select
      return
      end subroutine set_init
