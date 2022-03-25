#include "config.f"
!**********************************************************************
!>
!!    getstark obtains the internal pitch angles
!!    from polarimetry measurement using Wroblewski's routine
!!    
!!    WARNING: this subroutine uses both REAL*4 (used by mselib) and
!!             REAL*8 variables conversions must be handled carefully
!!
!!
!!    @param ktime : number of time slices
!!
!**********************************************************************
      subroutine getstark(ktime)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      real*4 avem,tanham(ktime,nmtark),sigham(ktime,nmtark), &
         rrham(nmtark),zzham(nmtark), &
         sarkar,sarkaz,a1ham(nmtark), &
         a2ham(nmtark),a3ham(nmtark),a4ham(nmtark), &
         a5ham(nmtark),a6ham(nmtark),a7ham(nmtark),atime(ktime), &
         spatial_avg_ham(nmtark,ngam_vars,ngam_u,ngam_w), &
         hgain(nmtark),hslope(nmtark),hscale(nmtark), &
         hoffset(nmtark),max_beamOff, &
         tanham_uncor(ktime,nmtark)
      real*4 fv30lt,fv30rt,fv210lt,fv210rt

      atime(1:ktime) = real(time(1:ktime),r4)
      if (dtmsefull .gt. 0.0) then
        avem = real(dtmsefull,r4) / 1000.0
      else
        avem = 2.0*iavem / 1000.0
      endif
      max_beamOff = real(t_max_beam_off,r4) / 1000.0
      call  set_mse_beam_logic(mse_strict,max_beamOff,ok_210lt,ok_30rt)
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
      do n=1,nmtark
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
          a8gam(i,n)=0.0
          if (abs(tangam(i,n)).le.1.e-10_dp.and. &
            abs(siggam(i,n)).le.1.e-10_dp) then
            fwtgam(n)=0.0
            siggam(i,n)=0.0
          elseif (abs(tangam(i,n)).le.1.e-10_dp.and. &
            abs(siggam(i,n)).le.100.0) then
            fwtgam(n)=0.0
            siggam(i,n)=0.0
          elseif (iergam(n).gt.0) then
            fwtgam(n)=0.0
            siggam(i,n)=0.0
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
        integer*4,dimension(:),allocatable :: tmp1,tmp2
        double precision :: zwork(nmtark*12,ntime)

#ifdef DEBUG_LEVEL1
        ! timing variables
        integer*4 :: clock0,clock1,clockmax,clockrate,ticks
        double precision :: secs
        integer*4 :: total_bytes
#endif
        
        nsize = 12*nmtark
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
          write (*,"(' GETSTARK call ',f6.2,' sec')") secs
          call system_clock(count_max=clockmax,count_rate=clockrate)
          call system_clock(count=clock0)
#endif
        endif timing_rank
        
        ! Process with rank == 0 gets distributes data
        if (rank == 0) then
          ! Pack ZWORK data array
          ! NOTE : Transposing data arrays (packing as consecutive chunks into each column of ZWORK array)
          do j=1,nmtark
            zwork(j,1:ktime_all)           = tangam(1:ktime_all,j)
            zwork(j+nmtark,1:ktime_all)    = siggam(1:ktime_all,j)
            zwork(j+nmtark*2,1:ktime_all)  = rrgam(1:ktime_all,j)
            zwork(j+nmtark*3,1:ktime_all)  = zzgam(1:ktime_all,j)
            zwork(j+nmtark*4,1:ktime_all)  = a1gam(1:ktime_all,j)
            zwork(j+nmtark*5,1:ktime_all)  = a2gam(1:ktime_all,j)
            zwork(j+nmtark*6,1:ktime_all)  = a3gam(1:ktime_all,j)
            zwork(j+nmtark*7,1:ktime_all)  = a4gam(1:ktime_all,j)
            zwork(j+nmtark*8,1:ktime_all)  = a5gam(1:ktime_all,j)
            zwork(j+nmtark*9,1:ktime_all)  = a6gam(1:ktime_all,j)
            ! WARNING : A7GAM explicitly set to zero by original
            ! GETSTARK_MPI code
            zwork(j+nmtark*10,1:ktime_all) = a7gam(1:ktime_all,j)
            !zwork(j+nmtark*10,1:ktime_all) = 0.0
            ! NOTE : Do NOT actually need to pack this data since
            !        array explicitly set to zero and we could just set A8GAM
            !        to zero
            zwork(j+nmtark*11,1:ktime_all) = 0.0
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
          do j=1,nmtark
            tangam(1:ktime,j) = zwork(j,1:ktime)
            siggam(1:ktime,j) = zwork(j+nmtark,1:ktime)
            rrgam(1:ktime,j)  = zwork(j+nmtark*2,1:ktime)
            zzgam(1:ktime,j)  = zwork(j+nmtark*3,1:ktime)
            a1gam(1:ktime,j)  = zwork(j+nmtark*4,1:ktime)
            a2gam(1:ktime,j)  = zwork(j+nmtark*5,1:ktime)
            a3gam(1:ktime,j)  = zwork(j+nmtark*6,1:ktime)
            a4gam(1:ktime,j)  = zwork(j+nmtark*7,1:ktime)
            a5gam(1:ktime,j)  = zwork(j+nmtark*8,1:ktime)
            a6gam(1:ktime,j)  = zwork(j+nmtark*9,1:ktime)
            a7gam(1:ktime,j)  = zwork(j+nmtark*10,1:ktime)
            a8gam(1:ktime,j)  = zwork(j+nmtark*11,1:ktime)
          enddo
        endif

        ! KWAITMSE
        ! NOTE : Necessary to send KWAITMSE to ALL processes since controls main loop defined in EFITD
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
        total_bytes = total_bytes + 8*nmtark*(nproc-1)
#endif
        call MPI_BCAST(fwtgam,nmtark,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        !fwtgam(:) = fwtgam_mpi(:,rank+1)

        !!call MPI_BCAST(msefitfun,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        !!call MPI_BCAST(msebkp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(iergam,nstark,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(rmse_gain,nmtark,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(rmse_slope,nmtark,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(rmse_scale,nmtark,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(rmse_offset,nmtark,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(mse_spave_on,nmtark,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        ! SPATIAL_AVG_GAM
        call MPI_BCAST(spatial_avg_gam,nstark*ngam_vars*ngam_u*ngam_w, &
                       MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        !! SWTGAM
        !tmp1(:) = dist_data(:)
        !tmp2(:) = 0
        !! SIZE = SIZEOF(DOUBLE) * SUM(DIST_DATA(2:)) * NSIZE bytes
        !!total_bytes = total_bytes + 8*sum(dist_data(2:))*nsize
        !if (rank == 0) then
        !  ! NOTE : DIST_DATA and DIST_DATA_DISPLS should be saved between calls since part of MPI_INFO module
        !  call MPI_SCATTERV(swtgam,tmp1,tmp2,MPI_DOUBLE_PRECISION,MPI_IN_PLACE,tmp1(rank+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        !else
        !  call MPI_SCATTERV(swtgam,tmp1,tmp2,MPI_DOUBLE_PRECISION,swtgam,tmp1(rank+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        !endif
        !!call MPI_BCAST(swtgam,sum(dist_data),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        deallocate(tmp1,tmp2)
        
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
#ifdef DEBUG_LEVEL1
        timing_rank0: if (rank == 0) then
          call system_clock(count=clock1)
          ticks = clock1-clock0
          secs = real(ticks,dp)/real(clockrate,dp)
          write (*,"(' GETSTARK transfer ',i10,' bytes in ',f6.2,'sec')") &
                total_bytes,secs
        endif timing_rank0
#endif

      end subroutine getstark_mpi

#endif
