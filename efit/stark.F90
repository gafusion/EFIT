#include "config.f"
#if defined(USE_MSE)
!**********************************************************************
!>
!!    getstark obtains the internal pitch angles
!!    from polarimetry measurement using Wroblewski's routine
!!    
!!    WARNING: this subroutine uses both REAL*4 (used by mselib) and
!!             REAL*8 variables conversions must be handled carefully
!!
!!    @param ktime : number of time slices
!!
!**********************************************************************
      subroutine getstark(ktime)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      real*4 avem,tanham(ktime,nmselp),sigham(ktime,nmselp), &
         rrham(nmselp),zzham(nmselp), &
         sarkar,sarkaz,a1ham(nmselp), &
         a2ham(nmselp),a3ham(nmselp),a4ham(nmselp), &
         a5ham(nmselp),a6ham(nmselp),a7ham(nmselp),atime(ktime), &
         spatial_avg_ham(nmselp,ngam_vars,ngam_u,ngam_w), &
         hgain(nmselp),hslope(nmselp),hscale(nmselp), &
         hoffset(nmselp),max_beamOff, &
         tanham_uncor(ktime,nmselp)
      real*4 fv30lt,fv30rt,fv210lt,fv210rt

      atime(1:ktime) = real(time(1:ktime),r4)
      if (dtmsefull .gt. 0.0) then
        avem = real(dtmsefull,r4)/1000.0
      else
        avem = 2.0*iavem/1000.0
      endif
      max_beamOff = real(t_max_beam_off,r4)/1000.0
      call set_mse_beam_logic(mse_strict,max_beamOff,ok_210lt,ok_30rt)
      tanham = 0.0
      tanham_uncor = 0.0
      sigham = 0.0
      fv30lt = real(v30lt,r4)
      fv30rt = real(v30rt,r4)
      fv210rt = real(v210rt,r4)
      fv210lt = real(v210lt,r4)
      call set_cer_correction(mse_usecer,mse_certree, &
                              mse_use_cer330,mse_use_cer210)
      call set_mse_beam_on_vlevel(fv30lt,fv210rt,fv210lt,fv30rt)
      call stark2cer(ishot,atime,ktime,avem,msefitfun,tanham,sigham, &
                     rrham,zzham,a1ham,a2ham,a3ham,a4ham,a5ham,a6ham, &
                     a7ham,iergam,msebkp,mse_quiet,tanham_uncor)
      call get_mse_spatial_data(spatial_avg_ham)
      call get_mse_calibration(msefitfun,hgain,hslope,hscale,hoffset)
      kfixstark = 0
      do n=1,nmselp
        rmse_gain(n) = real(hgain(n),dp)
        rmse_slope(n) = real(hslope(n),dp)
        rmse_scale(n) = real(hscale(n),dp)
        rmse_offset(n) = real(hoffset(n),dp)
        if(mse_spave_on(n) .ne. 0) kfixstark = 1
        do i=1,ngam_vars
          do j=1,ngam_u
            spatial_avg_gam(n,i,j,1:ngam_w) =  &
              real(spatial_avg_ham(n,i,j,1:ngam_w),dp)
          enddo
        enddo

        do i=1,ktime
          tangam(i,n) = real(tanham(i,n),dp)
          tangam_uncor(i,n) = real(tanham_uncor(i,n),dp)
          siggam(i,n) = real(sigham(i,n),dp)
          rrgam(i,n) = real(rrham(n),dp)
          zzgam(i,n) = real(zzham(n),dp)
          starkar(i,n) = real(sarkar,dp)
          starkaz(i,n) = real(sarkaz,dp)
          a1gam(i,n) = real(a1ham(n),dp)
          a2gam(i,n) = real(a2ham(n),dp)
          a3gam(i,n) = real(a3ham(n),dp)
          a4gam(i,n) = real(a4ham(n),dp)
          a5gam(i,n) = real(a5ham(n),dp)
          a6gam(i,n) = real(a6ham(n),dp)
          a7gam(i,n) = real(a7ham(n),dp)
          a8gam(i,n) = 0.0
          if (abs(tangam(i,n)).le.1.e-10_dp.and. &
            abs(siggam(i,n)).le.100.0) then
            siggam(i,n) = 0.0
          elseif (iergam(n).gt.0) then
            siggam(i,n) = 0.0
          endif
        enddo
      enddo
      !
      ! kfixstark is zero when none of the individual channels is turned on
      ! in this case set kwaitmse to zero which turns the mse spacial average
      ! off globally
      !
      if (kfixstark .eq. 0) then
        kwaitmse = 0
      elseif (kwaitmse .eq. 0) then
        kwaitmse = 5
      endif

      return
      end subroutine getstark

! =========================================================

#if defined(USEMPI)

!**********************************************************************
!>
!!    wrapper around getstark subroutine to handle MPI communication\n
!!    NOTE : NO error condition returned
!!
!!    @param ktime : number of time slices
!!
!**********************************************************************
      subroutine getstark_mpi(ktime)

        include 'eparm.inc'
        include 'modules1.inc'
        implicit integer*4 (i-n), real*8 (a-h,o-z)
        include 'mpif.h'

        ! NOTE : 12 scalar used because distributing 12 arrays amongst processes
        ! Dimensions ZWORK 2-D array
        !   TANGAM (ntime,nstark)
        !   SIGGAM (ntime,nstark)
        !   RRGAM  (ntime,nstark)
        !   ZZGAM  (ntime,nstark)
        !   A1GAM  (ntime,nstark)
        !   A2GAM  (ntime,nstark)
        !   A3GAM  (ntime,nstark)
        !   A4GAM  (ntime,nstark)
        !   A5GAM  (ntime,nstark)
        !   A6GAM  (ntime,nstark)
        !   A7GAM  (ntime,nstark)
        !   0.0
        ! WARNING : nsize < nstark (OKAY)
        ! WARNING : A7GAM explicitly set to zero by original GETSTARK_MPI code
        ! KWAITMSE
        ! FWTGAM (nstark)
        integer*4 :: i,j,ktime_all,offset,nsize
        integer*4, dimension(:), allocatable :: tmp1,tmp2
        double precision :: zwork(nmselp*12,ntime)

#ifdef DEBUG_LEVEL1
        ! timing variables
        integer*4 :: clock0,clock1,clockmax,clockrate,ticks
        double precision :: secs
        integer*4 :: total_bytes
#endif

        nsize=12*nmselp
        zwork(:,:) = 0.0
        allocate(tmp1(nproc),tmp2(nproc))

        ! Process with rank == 0 gets data from PTDATA/MDS+ database by calling GETSTARK
        timing_rank: if (rank == 0) then
          ! NOTE : Need to retrive data for ALL times
          ktime_all = sum(dist_data)
          ! TODO: it's not obvious why this info is always output for mpi
          ! runs, but not serial, so I am making it a debug option (for now)
#ifdef DEBUG_LEVEL1
          call system_clock(count_max=clockmax,count_rate=clockrate)
          call system_clock(count=clock0)
#endif
          call getstark(ktime_all)
#ifdef DEBUG_LEVEL1
          call system_clock(count=clock1)
          ticks = clock1-clock0
          secs = real(ticks,dp)/real(clockrate,dp)
          write(*,"(' GETSTARK call ',f6.2,' sec')") secs
          call system_clock(count_max=clockmax,count_rate=clockrate)
          call system_clock(count=clock0)
#endif
        endif timing_rank
        
        ! Process with rank == 0 gets distributes data
        if (rank == 0) then
          ! Pack ZWORK data array
          ! NOTE : Transposing data arrays (packing as consecutive chunks into each column of ZWORK array)
          do j=1,nmselp
            zwork(j,1:ktime_all)           = tangam(1:ktime_all,j)
            zwork(j+nmselp,1:ktime_all)    = siggam(1:ktime_all,j)
            zwork(j+nmselp*2,1:ktime_all)  = rrgam(1:ktime_all,j)
            zwork(j+nmselp*3,1:ktime_all)  = zzgam(1:ktime_all,j)
            zwork(j+nmselp*4,1:ktime_all)  = a1gam(1:ktime_all,j)
            zwork(j+nmselp*5,1:ktime_all)  = a2gam(1:ktime_all,j)
            zwork(j+nmselp*6,1:ktime_all)  = a3gam(1:ktime_all,j)
            zwork(j+nmselp*7,1:ktime_all)  = a4gam(1:ktime_all,j)
            zwork(j+nmselp*8,1:ktime_all)  = a5gam(1:ktime_all,j)
            zwork(j+nmselp*9,1:ktime_all)  = a6gam(1:ktime_all,j)
            ! WARNING : A7GAM explicitly set to zero by original
            ! GETSTARK_MPI code
            zwork(j+nmselp*10,1:ktime_all) = a7gam(1:ktime_all,j)
            !zwork(j+nmselp*10,1:ktime_all) = 0.0
            ! NOTE : Do NOT actually need to pack this data since
            !        array explicitly set to zero and we could just set A8GAM
            !        to zero
            zwork(j+nmselp*11,1:ktime_all) = 0.0
          enddo
        endif
        ! Distribute chunks of ZWORK array to processes
        ! NOTE : We need to recalculate distribution data
        tmp1(:) = dist_data(:)*nsize
        tmp2(:) = dist_data_displs(:)*nsize
        ! SIZE = SIZEOF(DOUBLE) * SUM(DIST_DATA(2:)) * NSIZE bytes
#ifdef DEBUG_LEVEL1
        total_bytes = 8*sum(dist_data(2:))*nsize
#endif
        if (rank == 0) then
          ! NOTE : DIST_DATA and DIST_DATA_DISPLS should be saved between calls since part of MPI_INFO module
          call MPI_SCATTERV(zwork,tmp1,tmp2,MPI_DOUBLE_PRECISION, &
                            MPI_IN_PLACE,tmp1(rank+1), &
                            MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        else
          call MPI_SCATTERV(zwork,tmp1,tmp2,MPI_DOUBLE_PRECISION,zwork, &
                            tmp1(rank+1),MPI_DOUBLE_PRECISION,0, &
                            MPI_COMM_WORLD,ierr)
        endif
        ! Unpack ZWORK array data
        ! NOTE : Only processes with rank > 0 need to unpack data
        ! Determine local KTIME value before proceeding
        ktime = dist_data(rank+1)
        if (rank > 0) then
          do j=1,nmselp
            tangam(1:ktime,j) = zwork(j,1:ktime)
            siggam(1:ktime,j) = zwork(j+nmselp,1:ktime)
            rrgam(1:ktime,j)  = zwork(j+nmselp*2,1:ktime)
            zzgam(1:ktime,j)  = zwork(j+nmselp*3,1:ktime)
            a1gam(1:ktime,j)  = zwork(j+nmselp*4,1:ktime)
            a2gam(1:ktime,j)  = zwork(j+nmselp*5,1:ktime)
            a3gam(1:ktime,j)  = zwork(j+nmselp*6,1:ktime)
            a4gam(1:ktime,j)  = zwork(j+nmselp*7,1:ktime)
            a5gam(1:ktime,j)  = zwork(j+nmselp*8,1:ktime)
            a6gam(1:ktime,j)  = zwork(j+nmselp*9,1:ktime)
            a7gam(1:ktime,j)  = zwork(j+nmselp*10,1:ktime)
            a8gam(1:ktime,j)  = zwork(j+nmselp*11,1:ktime)
          enddo
        endif

        ! KWAITMSE
        ! NOTE : Necessary to send KWAITMSE to ALL processes since controls main loop defined in EFIT
        ! SIZE = SIZEOF(INTEGER) * (NPROC - 1)
#ifdef DEBUG_LEVEL1
        total_bytes = total_bytes + 4*(nproc-1)
#endif
        call MPI_BCAST(kwaitmse,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        
        ! Distribute chunks of FWTGAM array
        ! NOTE : We need to recalculate distribution data
        ! WARNING : Uncertain if FWTGAM should be broadcast or if doing so could possible cause issues
        ! SIZE = SIZEOF(DOUBLE) * NMTARK * (NPROC - 1)
#ifdef DEBUG_LEVEL1
        total_bytes = total_bytes + 8*nmselp*(nproc-1)
#endif

        !!call MPI_BCAST(msefitfun,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        !!call MPI_BCAST(msebkp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(iergam,nstark,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(rmse_gain,nmselp,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(rmse_slope,nmselp,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(rmse_scale,nmselp,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(rmse_offset,nmselp,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(mse_spave_on,nmselp,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        ! SPATIAL_AVG_GAM
        call MPI_BCAST(spatial_avg_gam,nstark*ngam_vars*ngam_u*ngam_w, &
                       MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        deallocate(tmp1,tmp2)
        
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
#ifdef DEBUG_LEVEL1
        timing_rank0: if (rank == 0) then
          call system_clock(count=clock1)
          ticks = clock1-clock0
          secs = real(ticks,dp)/real(clockrate,dp)
          write(*,"(' GETSTARK transfer ',i10,' bytes in ',f6.2,'sec')") &
                total_bytes,secs
        endif timing_rank0
#endif

      end subroutine getstark_mpi
#endif
#endif

!**********************************************************************
!>
!!    setstark obtains the response matrix
!!    for polarimetry measurement
!!
!!    @param jtime : time index
!!
!**********************************************************************
      subroutine setstark(jtime)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      dimension rsplt(2500),zsplt(2500),csplt(2500), &
                rrgamin(nstark),zzgamin(nstark)
      character*40 filenmme
      character*2 ich
      data nset/20/,cdum/1.0/
!---------------------------------------------------------------------
!--   try read in, first locally, then efit area                    --
!---------------------------------------------------------------------
      do iname=1,nset
        write(ich,'(i2)') iname
        filenmme='rs'//trim(ch1)//trim(ch2)//'_'//trim(ich)//'.ddd'
        open(unit=nffile, &
             status='old',form='unformatted', &
             file=filenmme,iostat=ioerr)
        if(ioerr.ne.0) &
          open(unit=nffile, &
               status='old',form='unformatted', &
               file=table_dir(1:ltbdir)//filenmme,iostat=ioerr)
        if(ioerr.ne.0) exit
        read(nffile,iostat=ioerr) kkstark
        if(ioerr.ne.0) exit
        if (kkstark.ne.nstark) then
          close(unit=nffile)
          cycle
        endif
        read(nffile,iostat=ioerr)  rrgamin
        if(ioerr.ne.0) exit
        read(nffile,iostat=ioerr)  zzgamin
        if(ioerr.ne.0) exit
        do ii=1,nstark
          if(abs(rrgamin(ii)-rrgam(jtime,ii)).gt.1.e-4_dp) exit
          if(abs(zzgamin(ii)-zzgam(jtime,ii)).gt.1.e-4_dp) exit
        enddo
        if(abs(rrgamin(ii)-rrgam(jtime,ii)).gt.1.e-4_dp) cycle
        if(abs(zzgamin(ii)-zzgam(jtime,ii)).gt.1.e-4_dp) cycle
        read(nffile,iostat=ioerr) rbrfc
        if(ioerr.ne.0) exit
        read(nffile,iostat=ioerr) rbzfc
        if(ioerr.ne.0) exit
        read(nffile,iostat=ioerr) gbrpc
        if(ioerr.ne.0) exit
        read(nffile,iostat=ioerr) gbzpc
        if(ioerr.ne.0) exit
        read(nffile,iostat=ioerr) rbrec
        if(ioerr.ne.0) exit
        read(nffile,iostat=ioerr) rbzec
        if(ioerr.ne.0) exit
        close(unit=nffile)
        return
      enddo
!---------------------------------------------------------------------
!--   response due to F coils                                       --
!---------------------------------------------------------------------
      isplit=10
      itot=isplit*isplit
      fitot=itot
      rbrfc=0.0
      rbzfc=0.0
      do k=1,nfcoil
        call splitc(isplit,rsplt,zsplt,csplt, &
                    rf(k),zf(k),wf(k),hf(k),af(k),af2(k),cdum)
        do mmm=1,nstark
          if (rrgam(jtime,mmm).le.1.e-8_dp)  cycle
          brct=0.0
          bzct=0.0
          r1=rrgam(jtime,mmm)
          do l=1,itot
            a=rsplt(l)
            z1=zzgam(jtime,mmm)-zsplt(l)
            brct=brct+br(a,r1,z1)
            bzct=bzct+bz(a,r1,z1)
          enddo
          kkm=fcid(k)
          rbrfc(mmm,kkm)=rbrfc(mmm,kkm)+fcturn(k)*brct/fitot*tmu
          rbzfc(mmm,kkm)=rbzfc(mmm,kkm)+fcturn(k)*bzct/fitot*tmu
        enddo
      enddo
!---------------------------------------------------------------------
!--   plasma response                                               --
!---------------------------------------------------------------------
      zzx=0.0
      drzg=(drgrid**2+dzgrid**2)*25.
      isplit=18
      kk=1
      do mmm=1,nstark
        gbrpc(mmm,kk)=0.0
        gbzpc(mmm,kk)=0.0
        if(rrgam(jtime,mmm).le.1.e-8_dp) cycle
        r=rrgam(jtime,mmm)
        do ii=1,nw
          a=rgrid(ii)
          do jj=1,nh
            z=zzgam(jtime,mmm)-zgrid(jj)
            kk=(ii-1)*nh+jj
            dist=(a-r)**2+z**2
            if (dist.lt.drzg) then
   68        call splitc(isplit,rsplt,zsplt,csplt, &
                  rgrid(ii),zgrid(jj),drgrid,dzgrid,zzx,zzx,cdum)
             gbrt=0.0
             gbzt=0.0
             itot=isplit*isplit
             do k=1,itot
               a=rsplt(k)
               z=zzgam(jtime,mmm)-zsplt(k)
               distt=(a-r)**2+z**2
               if (distt.lt.1.e-8_dp.and.isplit.lt.49) then
                 isplit=isplit+2
                 go to 68
               endif
               gbrt=gbrt+br(a,r,z)
               gbzt=gbzt+bz(a,r,z)
             enddo
             gbrpc(mmm,kk)=gbrt*tmu/itot
             gbzpc(mmm,kk)=gbzt*tmu/itot
            else
             gbrpc(mmm,kk)=br(a,r,z)*tmu
             gbzpc(mmm,kk)=bz(a,r,z)*tmu
            endif
          enddo
        enddo
      enddo
!---------------------------------------------------------------------
!--   E coils response                                              --
!---------------------------------------------------------------------
      isplit=10
      itot=isplit*isplit
      fitot=real(itot,dp)
      do m=1,nstark
        rbrec(m,1:nesum)=0.0
        rbzec(m,1:nesum)=0.0
        if(rrgam(jtime,m).le.1.e-8_dp) cycle
        r1=rrgam(jtime,m)
        do k=1,necoil
          brct=0.0
          bzct=0.0
          call splitc(isplit,rsplt,zsplt,csplt, &
                      re(k),ze(k),we(k),he(k),zzx,zzx,cdum)
          do l=1,itot
            a=rsplt(l)
            z1=zzgam(jtime,m)-zsplt(l)
            brct=brct+br(a,r1,z1)
            bzct=bzct+bz(a,r1,z1)
          enddo
          brct=brct*tmu/fitot
          bzct=bzct*tmu/fitot
          kkm=ecid(k)
          rbrec(m,kkm)=rbrec(m,kkm)+brct
          rbzec(m,kkm)=rbzec(m,kkm)+bzct
        enddo
      enddo
!
! --- write out rstarkxx.dat if flag IOUT contains 8.
!
      if (iand(iout,8).ne.0) then
        open(unit=nffile,status='old',form='unformatted',iostat=ioerr, &
             file='rstarkxx.dat')
        if(ioerr.eq.0) close(unit=nffile,status='delete')
        open(unit=nffile,status='new',form='unformatted', &
             file='rstarkxx.dat')
        kkstark=nstark
        write(nffile) kkstark
        rrgamin(1:nstark)=rrgam(jtime,1:nstark)
        zzgamin(1:nstark)=zzgam(jtime,1:nstark)
        write(nffile) rrgamin
        write(nffile) zzgamin
        write(nffile) rbrfc
        write(nffile) rbzfc
        write(nffile) gbrpc
        write(nffile) gbzpc
        write(nffile) rbrec
        write(nffile) rbzec
        close(unit=nffile)
      endif
      return
      end subroutine setstark

!**********************************************************************
!>
!!    fixstark adjusts the internal pitch angles
!!    based on spatial averaging data
!!
!!    @param jtime :
!!
!!    @param kerror :
!!
!**********************************************************************
      subroutine fixstark(jtime,kerror)
      use commonblocks,only: ct,wkt,bkrt,bkzt
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      real*8 ppcurr,fpcurr,fpecrr,seval

      integer*4, intent(in) :: jtime
      integer*4, intent(out) :: kerror
      integer*4 i,ichan,ii,j,ier,ip_sign,lkrt,lkzt,ltime
      real*8 brl,btl,bzl,delsi,psi_norm,rl,siii,ssibry, &
             ssimag,sumf,tglocal,tl,ttl,zl
      real*8 pds(6)
      real*8,dimension(:),allocatable :: bwork,cwork,dwork

      kerror = 0
!
      allocate(bwork(nw),cwork(nw),dwork(nw))
!
      if (keecur .gt. 0 .and. ifirst .eq. 0) then
        write(6,*) "Spatial averaging correction of MSE data"
        write(6,*) "not supported with Er fit at this time."
        write(6,*) "Going ahead with correction anyway,at this"
        write(6,*) "time, on channels requested."
!        write(6,*) "Calculating but not applying correction."
!        where(mse_spave_on .ne. 0) mse_spave_on = -1
      endif

      ltime = jtime
      if (jtime .lt. 0) then
        ltime = -jtime
!       do ichan=1,nstark
!         tangam(ltime,ichan) = save_tangam(ltime,ichan)
!       enddo
      endif

      if (ifirst .eq. 0) then
        ifirst = 1
        write(6,*)"Calculate pitch angle corrections", &
          " using spatial averaging"
        do ichan=1,nstark
          save_gam(ltime,ichan) = atan(tangam(ltime,ichan))
          save_tangam(ltime,ichan) = tangam(ltime,ichan)
        enddo
        return
      endif

      ip_sign = -ipmhd(ltime)/abs(ipmhd(ltime))

      call sets2d(psi,ct,rgrid,nw,bkrt,lkrt,zgrid,nh,bkzt,lkzt,wkt,ier)
  
      if (ipmeas(ltime).gt.0.0) then
        ssimag=-simag
        ssibry=-psibry
      else
        ssimag=simag
        ssibry=psibry
      endif

      nonzerotime: if (jtime .gt. 0.0) then
!---------------------------------------------------------------------
!--   set up P' and FF', then integration                           --
!--   ffprim = (RBt) * d/dpsi(RBt)                                  --
!---------------------------------------------------------------------
      select case (icurrt)
      case (1)
        pprime(1)=cratio*sbeta/darea/srma
        ffprim(1)=cratio*srma*2.*salpha/darea*twopi*tmu
        pprime(nw)=pprime(1)
        ffprim(nw)=ffprim(1)
      case (2,5)
        pprime(nw)=ppcurr(x111,kppcur)/darea
        ffprim(nw)=fpcurr(x111,kffcur)/darea*twopi*tmu
        pprime(1)=ppcurr(x000,kppcur)/darea
        ffprim(1)=fpcurr(x000,kffcur)/darea*twopi*tmu
        if (kfffnc.eq.8) then
          ffprec(nw)=fpecrr(x111,kffcur)/darea*twopi*tmu
          ffprec(1)=fpecrr(x000,kffcur)/darea*twopi*tmu
        else
          ffprec(nw)=0.0
          ffprec(1)=0.0
        endif
      case (4)
        call currnt(n222,jtime,n222,kerror)
        if (kerror.gt.0) return
        pprime(1)=cratio/darea/rzero
        ffprim(1)=rbetap*cratio*rzero*twopi*tmu/darea
        ffprim(nw)=ffprim(1)*gammaf
        pprime(nw)=pprime(1)*gammap
      end select

      do i=2,nw-1
        ii=nw-i+1
        siii=1.0_dp-1.0_dp/(nw-1)*(i-1)
        sigrid(ii)=siii
        select case (icurrt)
        case(1)
          pprime(ii)=pprime(1)
          ffprim(ii)=ffprim(1)
        case (2,5)
          pprime(ii)=ppcurr(siii,kppcur)/darea
          ffprim(ii)=fpcurr(siii,kffcur)/darea*twopi*tmu
          if (kfffnc.eq.8) then
            ffprec(ii)=fpecrr(siii,kffcur)/darea*twopi*tmu
          else
            ffprec(ii)=0.0
          endif
        case(4)
          pprime(ii)=(1.-siii**enp)**emp*(1.-gammap)+gammap
          ffprim(ii)=ffprim(1)*pprime(ii)
          pprime(ii)=pprime(1)*pprime(ii)
        end select
      enddo
      endif nonzerotime

      fpol(nw)=fbrdy*tmu
      sumf=fpol(nw)**2/2.
      delsi=-(psibry+simag)/(nw-1)
      do i=1,nw-1
        sumf=sumf+0.5_dp*delsi*(ffprim(nw-i+1)+ffprim(nw-i))
        if (sumf .ge. 0.0) then
          fpol(nw-i)=sqrt(2.*sumf)*fpol(nw)/abs(fpol(nw))
        else
          fpol(nw-i)=fpol(nw)
        endif
      enddo

      call zpline(nw,sigrid,fpol,bwork,cwork,dwork)

      do ichan=1,nmselp
        if (mse_spave_on(ichan) .ne. 0) then

          ttl = 0.0
          rl = rrgam(ltime,ichan)
          zl = zzgam(ltime,ichan)
          call seva2d(bkrt,lkrt,bkzt,lkzt,ct,rl,zl,pds,ier,n333)
          brl = -pds(3)/rl
          bzl = pds(2)/rl
          psi_norm = (ssimag-pds(1)/ip_sign)/(ssimag-ssibry)
          btl = seval(nw,abs(psi_norm),sigrid,fpol,bwork, &
                      cwork,dwork)/rl
          tglocal = (bzl*a1gam(ltime,ichan))/  &
            (btl*a2gam(ltime,ichan)+brl*a3gam(ltime,ichan) &
            +bzl*a4gam(ltime,ichan))


          do i=1,ngam_u
            do j=1,ngam_w
              rl = spatial_avg_gam(ichan,1,i,j)
              zl = spatial_avg_gam(ichan,2,i,j)
              call seva2d(bkrt,lkrt,bkzt,lkzt,ct,rl,zl,pds,ier,n333)
              brl = -pds(3)/rl
              bzl = pds(2)/rl
              psi_norm = (ssimag-pds(1)/ip_sign)/(ssimag-ssibry)
              btl = seval(nw,abs(psi_norm),sigrid,fpol,bwork, &
                          cwork,dwork)/rl
              tl = 0.0
              tl = tl+spatial_avg_gam(ichan,4,i,j)*btl
              tl = tl+spatial_avg_gam(ichan,5,i,j)*brl
              tl = tl+spatial_avg_gam(ichan,6,i,j)*bzl
              tl = spatial_avg_gam(ichan,3,i,j)*bzl/tl
              ttl = ttl+tl
              !if(jtime .lt. 0) write(7,'(I2,6F13.8)') ichan,rl,zl,btl,brl,bzl,tl
            enddo
          enddo
          !spatial_fix(ichan,ltime) = atan(cmgam(ichan,ltime))-atan(ttl)
          spatial_fix(ichan,ltime) = atan(tglocal)-atan(ttl)

          if (jtime.gt.0.and.mse_spave_on(ichan) .eq. 1) then
            tangam(ltime,ichan) = tan(save_gam(ltime,ichan) &
              -spatial_fix(ichan,ltime))
          endif
        endif
      enddo

      if (jtime .lt. 0) then
        ifirst = 0
      endif
      
      deallocate(bwork,cwork,dwork)

      return
      end subroutine fixstark
