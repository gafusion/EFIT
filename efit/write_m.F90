!**********************************************************************
!>
!!    Subroutine  writes measurements (MSE, magnetic)
!!    and some of their corresponding calculated values in
!!    a netCDF file of name m0sssss.nc/m0sssss.0tttt.
!!    
!!    The input variable IOUT decides the mode,
!!    IOUT.and.2.ne.0 - one file for all slices m0sssss.nc
!!    IOUT.and.4.ne.0 - one file for each slice m0sssss.0tttt
!!    IOUT.and.6.ne.0 - one file for all m0sssss.0tttt_nnnn
!!    tttt as the first time,nnnn as the number of slices
!!    
!!    WARNING: this subroutine uses both REAL*4 (to write files) and
!!             REAL*8 variables, conversions must be handled carefully
!!
!!
!!    @param ktime : total number of time slices 
!!
!!    @param ifirsttime : index of starting time  
!!
!!    @param ilast :  index of ending time  
!!
!!    @param itype :  1 called from main routine with time loop.\n
!!                    2 called from main routine out of time loop.\n
!!                    write all slices at one time    
!!
!**********************************************************************
      subroutine write_m(ktime,ifirsttime,ilast,itype)
      use set_kinds, only: r4
      include 'eparm.inc'
      include 'modules1.inc'
      implicit none
      include 'netcdf.inc'   ! from the netCDF package..
!                            ..this must be symlinked to local directory
      integer*4, intent(in) :: ktime,ifirsttime,ilast,itype
      integer*4 i,ijump,iitime,m,n,nmse,npress1,npresw1,kzeroj1
      integer*4 nceq,idim_time,idim_nstark,idim_nsilop,idim_magpri, &
                idim_nfsum,idim_nesum,idim_npress,idim_npresw, &
                idim_npcurn,idim_nitera,idim_nacoil,idim_1,idim_device, &
                id_device,id_shot,id_fungam,id_tangam,id_tangam_uncor, & 
                id_time,id_fgam,id_siggam,id_id_fwtgam,id_rrgam,id_zzgam, &
                id_a1gam,id_a2gam,id_a3gam,id_a4gam,id_a5gam,id_a6gam, &
                id_a7gam,id_a8gam,id_cmgam,id_chigam,id_msebkp, &
                id_mseport,id_gaingam,id_fwtgam,id_slopegam,id_offsetgam, &
                id_silopt,id_sigsil,id_fwtsi,id_csilop,id_id_saisil, &
                id_expmpi,id_sigmpi,id_fwtmp2,id_cmpr2,id_saimpi, &
                id_psiref,id_sigref,id_fwtref,id_csiref,id_saisref, &
                id_diamag,id_sigdia,id_fwtdia,id_cdflux,id_chidflux, &
                id_fccurt,id_scalegam,id_saisil,id_sigfcc,id_fwtfc, &
                id_ccbrsp,id_chifcc,id_eccurt,id_sigecc,id_fwtec, &
                id_cecurr,id_chiecc,id_curc79,id_curc139,id_curc199, &
                id_curiu30,id_curil30,id_curiu90,id_curil90, &
                id_curiu150,id_curil150,id_accurt,id_caccurt,id_plasma, &
                id_sigpasma,id_fwtpasma,id_ipmhd,id_chipasma,id_pressr, &
                id_rpress,id_zpress,id_sigpre,id_fwtpre,id_cpress, &
                id_saipre,id_presw,id_rpresw,id_zpresw,id_sigprw, &
                id_fwtprw,id_cpresw,id_saiprw,id_czmaxi,id_cchisq, &
                id_cerror,id_chifin,id_chitot,id_darea,id_xrsp, &
                id_sizeroj,id_vzeroj,idim_kzeroj
      real*8 xdum,vm3,betatnx
      character let
      integer*4 dim2(2),c11(2),cnn(2),imap(2),stride(2)
      character(len=4) last
      character(len=80) eqdsk,title
      character*30 sfname
! --- temporay variables to convert double to single
      real*4 zwork(ntime+nsilop+nstark+nfsum+nesum+magpri+npress), &
             zcsilop(nsilop,ntime), &
             zcmpr2(magpri,ntime),zccbrsp(nfsum,ntime), &
             zsilopt(ntime,nsilop),zexpmpi(ntime,magpri), &
             zfccurt(ntime,nfsum),zeccurt(ntime,nesum), &
             zaccurt(ntime,nacoil),ziter(nitera)
      real*4,dimension(:,:),allocatable :: zcmgam,zstark
      character*85 presstext
      character*109 preswtext
!-----------------------------------------------------------------------
!--   write out t(shot).(time)_X files                                --
!-----------------------------------------------------------------------
      plot_t: if ((kwripre.eq.11).and.(itype.eq.2)) then
        xdum=0.0
        call setfnmd('t',ishot,itime,sfname)
        sfname=sfname(1:13)//'_chi2'
        open(unit=74,status='old',file=sfname,iostat=ierr)
        if(ierr.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=ifirsttime,ilast
          write (74,*) time(i),chisq(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_error'
        open(unit=74,status='old',file=sfname,iostat=ierr)
        if(ierr.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=ifirsttime,ilast
          write (74,*) time(i),terror(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_j1ave'
        open(unit=74,status='old',file=sfname,iostat=ierr)
        if(ierr.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=ifirsttime,ilast
          write (74,*) time(i),cj1ave(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_li'
        open(unit=74,status='old',file=sfname,iostat=ierr)
        if(ierr.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=ifirsttime,ilast
          write (74,*) time(i),li(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_betat'
        open(unit=74,status='old',file=sfname,iostat=ierr)
        if(ierr.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=ifirsttime,ilast
          write (74,*) time(i),betat(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_q95'
        open(unit=74,status='old',file=sfname,iostat=ierr)
        if(ierr.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=ifirsttime,ilast
          write (74,*) time(i),q95(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_q0'
        open(unit=74,status='old',file=sfname,iostat=ierr)
        if(ierr.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=ifirsttime,ilast
          write (74,*) time(i),qm(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_q0'
        open(unit=74,status='old',file=sfname,iostat=ierr)
        if(ierr.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=ifirsttime,ilast
          write (74,*) time(i),qm(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_eout'
        open(unit=74,status='old',file=sfname,iostat=ierr)
        if(ierr.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=ifirsttime,ilast
          write (74,*) time(i),elong(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_vout'
        open(unit=74,status='old',file=sfname,iostat=ierr)
        if(ierr.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=ifirsttime,ilast
          vm3=volume(i)/1.e6_dp
          write (74,*) time(i),vm3,xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_betan'
        open(unit=74,status='old',file=sfname,iostat=ierr)
        if(ierr.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=ifirsttime,ilast
          pasman=ipmhd(i)/1.e4_dp/aminor(i)/abs(bcentr(i))
          pasman=pasman*rout(i)/100./rcentr
          betatnx=betat(i)/pasman
          write (74,*) time(i),betatnx,xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_zts'
        open(unit=74,status='old',file=sfname,iostat=ierr)
        if(ierr.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=ifirsttime,ilast
          write (74,*) time(i),zuperts(i),xdum,xdum
        enddo
        close(unit=74)
      endif plot_t
!
      if((iand(iout,2).eq.0).and.(iand(iout,4).eq.0)) &
        return
!
      call NCPOPT(NCVERBOS)  ! get nc error message but not fatal
!
! --- time-dependent and the first slice and first time called, 
! --- create file and enter define mode
! --- if iout.and.2.ne.0 = create m0sssss.nc
! --- if iout.and.6.ne.0 = create m0sssss.0tttt_nnnn, 
! ---    while tttt is the first time and nnn is the number of slices.
      ijump=0
      if ((iand(iout,2).ne.0).and.(iand(iout,4).ne.0).and. &
          ((ifirsttime.eq.1).and.(itype.eq.1))) then
        let = 'm'
        iitime = time(ifirsttime)
        call setfnmd(let,ishot,iitime,eqdsk)
        write(last,'(i4.4)') ktime
        eqdsk = eqdsk(1:13)//'_'//last
        nceq = NCCRE(eqdsk,NCCLOB,ierr)
!
      elseif ((iand(iout,2).ne.0).and.(iand(iout,4).eq.0).and. &
              (ifirsttime.eq.1).and.(itype.eq.1)) then
        write(eqdsk,"('m',i6.6,'.nc')") ishot
        nceq = NCCRE(eqdsk,NCCLOB,ierr) ! create file, overwrite if exists
!
! --- creates one file for each slice
!
      elseif ((iand(iout,4).ne.0).and.(iand(iout,2).eq.0).and. &
              (itype.eq.1)) then
        let = 'm'
        iitime = time(ifirsttime)
        call setfnmeq(itimeu,let,ishot,iitime,eqdsk)
        nceq = NCCRE(eqdsk,NCCLOB,ierr)
!
! --- time-dependent but NOT the first slice NOR the first time,
! --- or single slice but NOT the first time, 
! --- skip define mode and go to write directly
!
      else
        ijump=1
      endif
      
      ijump0: if (ijump.eq.0) then
      if(ierr.ne.0) return
!
! --- nc file has been created successfully, now define variables, etc.
!
      title = 'EFIT measurement file Mssssss.ttttt/Mssssss.nc'
      call NCAPTC(nceq,NCGLOBAL,'title',NCCHAR,46,title,ierr)
! 
! --- define unlimited time dimension and scalar and array dimensions
!
      idim_time = NCDDEF(nceq,'dim_time',NCUNLIM,ierr)
      idim_1 = NCDDEF(nceq,'dim_scalar',1,ierr)
      idim_device = NCDDEF(nceq,'dim_device',10,ierr)
      if (klibim.gt.0) then
        nmse=nstark
      else
        nmse=nmselp
      endif
      idim_nstark = NCDDEF(nceq,'dim_nstark',nmse,ierr)
      allocate(zcmgam(nmse,ntime),zstark(ntime,nmse))
      idim_nsilop = NCDDEF(nceq,'dim_nsilop',nsilop,ierr)
      idim_magpri = NCDDEF(nceq,'dim_magpri',magpri,ierr)
      idim_nfsum = NCDDEF(nceq,'dim_nfcoil',nfsum,ierr)
      idim_nesum  = NCDDEF(nceq,'dim_nesum' ,nesum ,ierr)
      npress1 = npress
      if(npress.eq.0) npress1 = 1
      idim_npress = NCDDEF(nceq,'dim_npress',npress1,ierr)
      npresw1 = npresw
      if(npresw.eq.0) npresw1 = 1
      idim_npresw = NCDDEF(nceq,'dim_npresw',npresw1,ierr)
      idim_npcurn = NCDDEF(nceq,'dim_npcurn',npcurn,ierr)
      idim_nitera = NCDDEF(nceq,'dim_nitera',nitera,ierr)
      idim_nacoil = NCDDEF(nceq,'dim_nacoil',nacoil,ierr)
      kzeroj1 = kzeroj
      if(kzeroj.eq.0) kzeroj1 = 1
      idim_kzeroj = NCDDEF(nceq,'dim_kzeroj',kzeroj1,ierr)
      dim2(2) = idim_time
!-----------------------------------------------------------------------
!--   define variables
!-----------------------------------------------------------------------
! --- equilibrium specification
!
      id_device = NCVDEF(nceq,'device',NCCHAR,1,idim_device,ierr)
      ! more consistent with EFIT, but will break OMFIT omas conversion...
      !id_device = NCVDEF(nceq,'machine',NCCHAR,1,idim_device,ierr) !
      call NCAPTC(nceq,id_device,'long_name',NCCHAR,15, &
                  'machine name',ierr)
      id_shot = NCVDEF(nceq,'shot',NCLONG,1,idim_1,ierr)
      call NCAPTC(nceq,id_shot,'long_name',NCCHAR,11, &
                  'shot number',ierr)
!
      id_time = NCVDEF(nceq,'time',NCFLOAT,1,idim_time,ierr)
      call NCAPTC(nceq,id_time,'units',NCCHAR,4,'msec',ierr)
!
! --- MSE measurements
!
      id_fungam = NCVDEF(nceq,'msefitfun',NCLONG,1,idim_1,ierr)
      call NCAPTC(nceq,id_fungam,'long_name',NCCHAR,16, &
                  'MSE fit function',ierr)

      dim2(1) = idim_nstark
      id_tangam = NCVDEF(nceq,'tangam',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_tangam,'long_name',NCCHAR,35, &
                  'tangent of measured MSE pitch angle',ierr)
!
      id_tangam_uncor = NCVDEF(nceq,'tangam_uncor',NCFLOAT,&
                                2,dim2,ierr)
      call NCAPTC(nceq,id_tangam_uncor,'long_name',NCCHAR,54, &
            'tangent of measured MSE pitch angle w/o cer correction',ierr)
!
      id_fgam = NCVDEF(nceq,'fixgam',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_fgam,'long_name',NCCHAR,58, &
        'radians correction of tangam for spatial averaging effects',ierr)
!
      id_siggam = NCVDEF(nceq,'siggam',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_siggam,'long_name',NCCHAR,21, &
                  'uncertainty of tangam',ierr)
!
      id_fwtgam = NCVDEF(nceq,'fwtgam',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_fwtgam,'long_name',NCCHAR,31, &
                  'fitting weight for MSE channels',ierr)
!
      id_rrgam = NCVDEF(nceq,'rrgam',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_rrgam,'long_name',NCCHAR,22, &
                  'radius of MSE channels',ierr)
!
      id_zzgam = NCVDEF(nceq,'zzgam',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_zzgam,'long_name',NCCHAR,26, &
                  'Z position of MSE channels',ierr)
!
      id_a1gam = NCVDEF(nceq,'a1gam',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_a1gam,'long_name',NCCHAR,45, &
                  'viewing geometry coefficients of MSE channels',ierr)
!
      id_a2gam = NCVDEF(nceq,'a2gam',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_a2gam,'long_name',NCCHAR,45, &
                  'viewing geometry coefficients of MSE channels',ierr)
!
      id_a3gam = NCVDEF(nceq,'a3gam',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_a3gam,'long_name',NCCHAR,45, &
                  'viewing geometry coefficients of MSE channels',ierr)
!
      id_a4gam = NCVDEF(nceq,'a4gam',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_a4gam,'long_name',NCCHAR,45, &
                  'viewing geometry coefficients of MSE channels',ierr)
!
      id_a5gam = NCVDEF(nceq,'a5gam',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_a5gam,'long_name',NCCHAR,45, &
                  'viewing geometry coefficients of MSE channels',ierr)
!
      id_a6gam = NCVDEF(nceq,'a6gam',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_a6gam,'long_name',NCCHAR,45, &
                  'viewing geometry coefficients of MSE channels',ierr)
!
      id_a7gam = NCVDEF(nceq,'a7gam',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_a7gam,'long_name',NCCHAR,45, &
                  'viewing geometry coefficients of MSE channels',ierr)
!
      id_a8gam = NCVDEF(nceq,'a8gam',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_a8gam,'long_name',NCCHAR,45, &
                  'viewing geometry coefficients of MSE channels',ierr)
!
      id_cmgam = NCVDEF(nceq,'cmgam',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_cmgam,'long_name',NCCHAR,30, &
                  'calculated polarimetry signals',ierr)
!
      id_chigam = NCVDEF(nceq,'chigam',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_chigam,'long_name',NCCHAR,23, &
                  'chisq vs. polarimetries',ierr)
!
      id_msebkp = NCVDEF(nceq,'msebkp',NCFLOAT,1,idim_time,ierr) ! why is this a float?
      call NCAPTC(nceq,id_msebkp,'long_name',NCCHAR,30, &
                  'background substraction switch',ierr)
!     
      id_mseport = NCVDEF(nceq,'vport',NCLONG,1,idim_nstark,ierr)
      call NCAPTC(nceq,id_mseport,'long_name',NCCHAR,33, &
                  'view port locations of MSE system',ierr)
!
      if (msefitfun .eq. 1) then
        id_gaingam = NCVDEF(nceq,'mcal_gain',NCFLOAT,2,dim2,ierr)
        call NCAPTC(nceq,id_gaingam,'long_name',NCCHAR,38, &
                    'gain param for tangent offset function',ierr)
!
        id_slopegam = NCVDEF(nceq,'mcal_slope',NCFLOAT, &
                             2,dim2,ierr)
        call NCAPTC(nceq,id_slopegam,'long_name',NCCHAR,39, &
                     'slope param for tangent offset function',ierr)
!
        id_scalegam = NCVDEF(nceq,'mcal_scale',NCFLOAT, &
                             2,dim2,ierr)
        call NCAPTC(nceq,id_scalegam,'long_name',NCCHAR,39, &
                     'scale param for tangent offset function',ierr)
!
        id_offsetgam = NCVDEF(nceq,'mcal_offset',NCFLOAT, &
                              2,dim2,ierr)
        call NCAPTC(nceq,id_offsetgam,'long_name',NCCHAR,40, &
                    'offset param for tangent offset function',ierr)
!
      elseif (msefitfun .eq. 3) then
        id_gaingam = NCVDEF(nceq,'mcal3_gain',NCFLOAT,2,dim2,ierr)
        call NCAPTC(nceq,id_gaingam,'long_name',NCCHAR,37, &
                    'gain param for tangent slope function',ierr)
!
        id_slopegam = NCVDEF(nceq,'mcal3_phase',NCFLOAT,2,dim2,ierr)
        call NCAPTC(nceq,id_slopegam,'long_name',NCCHAR,38, &
                    'phase param for tangent slope function',ierr)
!
        id_scalegam = NCVDEF(nceq,'mcal3_btscale',NCFLOAT, &
                             2,dim2,ierr)
        call NCAPTC(nceq,id_scalegam,'long_name',NCCHAR,40, &
                    'btscale param for tangent slope function',ierr)
!
        id_offsetgam = NCVDEF(nceq,'mcal3_dc_offset',NCFLOAT, &
                              2,dim2,ierr)
        call NCAPTC(nceq,id_offsetgam,'long_name',NCCHAR,42, &
                    'dc_offset param for tangent slope function',ierr)
!
      elseif (msefitfun .eq. 4)then
        id_gaingam = NCVDEF(nceq,'mcal4_gain',NCFLOAT,2,dim2,ierr)
        call NCAPTC(nceq,id_gaingam,'long_name',NCCHAR,37, &
                    'gain param for tangent slope function',ierr)
!
        id_slopegam = NCVDEF(nceq,'mcal4_phase',NCFLOAT,2,dim2,ierr)
        call NCAPTC(nceq,id_slopegam,'long_name',NCCHAR,38, &
                    'phase param for tangent slope function',ierr)
!
        id_scalegam = NCVDEF(nceq,'mcal4_btscale',NCFLOAT, &
                             2,dim2,ierr)
        call NCAPTC(nceq,id_scalegam,'long_name',NCCHAR,40, &
                    'btscale param for tangent slope function',ierr)
!
        id_offsetgam = NCVDEF(nceq,'mcal4_dc_offset',NCFLOAT, &
                              2,dim2,ierr)
        call NCAPTC(nceq,id_offsetgam,'long_name',NCCHAR,42, &
                    'dc_offset param for tangent slope function',ierr)
      endif
!
! --- magnetic measurements
!
      dim2(1) = idim_nsilop
      id_silopt = NCVDEF(nceq,'silopt',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_silopt,'long_name',NCCHAR,19, &
                  'measured flux loops',ierr)
!
      id_sigsil = NCVDEF(nceq,'sigsil',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_sigsil,'long_name',NCCHAR,25, &
                  'uncertainty in flux loops',ierr)
!
      id_fwtsi = NCVDEF(nceq,'fwtsi',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_fwtsi,'long_name',NCCHAR,21, &
                  'weight for flux loops',ierr)
!
      id_csilop = NCVDEF(nceq,'csilop',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_csilop,'long_name',NCCHAR,21, &
                  'calculated flux loops',ierr)
!
      id_saisil = NCVDEF(nceq,'saisil',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_saisil,'long_name',NCCHAR,20, &
                  'chisq for flux loops',ierr)
!
      dim2(1) = idim_magpri
      id_expmpi = NCVDEF(nceq,'expmpi',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_expmpi,'long_name',NCCHAR,24, &
                  'measured magnetic probes',ierr)
!
      id_sigmpi = NCVDEF(nceq,'sigmpi',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_sigmpi,'long_name',NCCHAR,30, &
                  'uncertainty in magnetic probes',ierr)
!
      id_fwtmp2 = NCVDEF(nceq,'fwtmp2',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_fwtmp2,'long_name',NCCHAR,26, &
                  'weight for magnetic probes',ierr)
!
      id_cmpr2 = NCVDEF(nceq,'cmpr2',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_cmpr2,'long_name',NCCHAR,26, &
                  'calculated magnetic probes',ierr)
!
      id_saimpi = NCVDEF(nceq,'saimpi',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_saimpi,'long_name',NCCHAR,25, &
                  'chisq for magnetic probes',ierr)
!
      id_psiref = NCVDEF(nceq,'psiref',NCFLOAT,1,idim_time,ierr)
      call NCAPTC(nceq,id_psiref,'long_name',NCCHAR,28, &
                  'measured reference flux loop',ierr)
!
      id_sigref = NCVDEF(nceq,'sigref',NCFLOAT,1,idim_time,ierr)
      call NCAPTC(nceq,id_sigref,'long_name',NCCHAR,34, &
                  'uncertainty in reference flux loop',ierr)
!
      id_fwtref = NCVDEF(nceq,'fwtref',NCFLOAT,1,idim_time,ierr)
      call NCAPTC(nceq,id_fwtref,'long_name',NCCHAR,30, &
                  'weight for reference flux loop',ierr)
!
      id_csiref = NCVDEF(nceq,'csiref',NCFLOAT,1,idim_time,ierr)
      call NCAPTC(nceq,id_csiref,'long_name',NCCHAR,30, &
                  'calculated reference flux loop',ierr)
!
      id_saisref = NCVDEF(nceq,'saisref',NCFLOAT,1,idim_time,ierr)
      call NCAPTC(nceq,id_saisref,'long_name',NCCHAR,29, &
                  'chisq for reference flux loop',ierr)
!
      id_diamag = NCVDEF(nceq,'diamag',NCFLOAT,1,idim_time,ierr)
      call NCAPTC(nceq,id_diamag,'long_name',NCCHAR,25, &
                  'measured diamagnetic flux',ierr)
!
      id_sigdia = NCVDEF(nceq,'sigdia',NCFLOAT,1,idim_time,ierr)
      call NCAPTC(nceq,id_sigdia,'long_name',NCCHAR,31, &
                  'uncertainty of diamagnetic flux',ierr)
! renamed fwtdlc to match MDS+ and rt-EFIT
      id_fwtdia = NCVDEF(nceq,'fwtdia',NCFLOAT,1,idim_time,ierr)
      call NCAPTC(nceq,id_fwtdia,'long_name',NCCHAR,27, &
                  'weight for diamagnetic flux',ierr)
!
      id_cdflux = NCVDEF(nceq,'cdflux',NCFLOAT,1,idim_time,ierr)
      call NCAPTC(nceq,id_cdflux,'long_name',NCCHAR,27, &
                  'calculated diamagnetic flux',ierr)
!
      id_chidflux = NCVDEF(nceq,'chidflux',NCFLOAT,1,idim_time,ierr)
      call NCAPTC(nceq,id_chidflux,'long_name',NCCHAR,26, &
                  'chisq for diamagnetic flux',ierr)
!
! --- coil currents
!
      dim2(1) = idim_nfsum
      id_fccurt = NCVDEF(nceq,'fccurt',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_fccurt,'long_name',NCCHAR,36, &
                  'measured F-coil currents (Amp-turns)',ierr)
!
      id_sigfcc = NCVDEF(nceq,'sigfcc',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_sigfcc,'long_name',NCCHAR,30, &
                  'uncertainty in F-coil currents',ierr)
!
      id_fwtfc = NCVDEF(nceq,'fwtfc',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_fwtfc,'long_name',NCCHAR,26, &
                  'weight for F-coil currents',ierr)
!
      id_ccbrsp = NCVDEF(nceq,'ccbrsp',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_ccbrsp,'long_name',NCCHAR,38, &
                  'calculated F-coil currents (Amp-turns)',ierr)
!
      id_chifcc = NCVDEF(nceq,'chifcc',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_chifcc,'long_name',NCCHAR,25, &
                  'chisq for F-coil currents',ierr)
!
      dim2(1) = idim_nesum
      id_eccurt = NCVDEF(nceq,'eccurt',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_eccurt,'long_name',NCCHAR,30, &
                  'measured E-coil currents (Amp)',ierr)
!
      id_sigecc = NCVDEF(nceq,'sigecc',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_sigecc,'long_name',NCCHAR,30, &
                  'uncertainty in E-coil currents',ierr)
!
      id_fwtec = NCVDEF(nceq,'fwtec',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_fwtec,'long_name',NCCHAR,26, &
                  'weight for E-coil currents',ierr)
!
      id_cecurr = NCVDEF(nceq,'cecurr',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_cecurr,'long_name',NCCHAR,32, &
                  'calculated E-coil currents (Amp)',ierr)
!
      id_chiecc = NCVDEF(nceq,'chiecc',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_chiecc,'long_name',NCCHAR,25, &
                  'chisq for E-coil currents',ierr)
!
      id_curc79 = NCVDEF(nceq,'curc79',NCFLOAT,1,idim_time,ierr)
      call NCAPTC(nceq,id_curc79,'long_name',NCCHAR,17, &
                  'C coil 79 current',ierr)
!
      id_curc139 = NCVDEF(nceq,'curc139',NCFLOAT,1,idim_time,ierr)
      call NCAPTC(nceq,id_curc139,'long_name',NCCHAR,18, &
                  'C coil 139 current',ierr)
!
      id_curc199 = NCVDEF(nceq,'curc199',NCFLOAT,1,idim_time,ierr)
      call NCAPTC(nceq,id_curc199,'long_name',NCCHAR,18, &
                  'C coil 199 current',ierr)
!
      id_curiu30 = NCVDEF(nceq,'curiu30',NCFLOAT,1,idim_time,ierr)
      call NCAPTC(nceq,id_curiu30,'long_name',NCCHAR,23, &
                  'I coil 30 upper current',ierr)
!
      id_curil30 = NCVDEF(nceq,'curil30',NCFLOAT,1,idim_time,ierr)
      call NCAPTC(nceq,id_curil30,'long_name',NCCHAR,23, &
                  'I coil 30 lower current',ierr)
!
      id_curiu90 = NCVDEF(nceq,'curiu90',NCFLOAT,1,idim_time,ierr)
      call NCAPTC(nceq,id_curiu90,'long_name',NCCHAR,23, &
                  'I coil 90 upper current',ierr)
!
      id_curil90 = NCVDEF(nceq,'curil90',NCFLOAT,1,idim_time,ierr)
      call NCAPTC(nceq,id_curil90,'long_name',NCCHAR,23, &
                  'I coil 90 lower current',ierr)
!
      id_curiu150 = NCVDEF(nceq,'curiu150',NCFLOAT,1,idim_time,ierr)
      call NCAPTC(nceq,id_curiu150,'long_name',NCCHAR,24, &
                  'I coil 150 upper current',ierr)
!
      id_curil150 = NCVDEF(nceq,'curil150',NCFLOAT,1,idim_time,ierr)
      call NCAPTC(nceq,id_curil150,'long_name',NCCHAR,24, &
                  'I coil 150 lower current',ierr)
!
      dim2(1) = idim_nacoil
      id_accurt = NCVDEF(nceq,'accurt',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_accurt,'long_name',NCCHAR,24, &
                  'measured A coil currents',ierr)
!
      id_caccurt = NCVDEF(nceq,'caccurt',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_caccurt,'long_name',NCCHAR,26, &
                  'calculated A coil currents',ierr)
!
! --- plasma current
!
      id_plasma = NCVDEF(nceq,'plasma',NCFLOAT,1,idim_time,ierr)
      call NCAPTC(nceq,id_plasma,'long_name',NCCHAR,29, &
                  'measured plasma current (Amp)',ierr)
!
      id_sigpasma = NCVDEF(nceq,'sigpasma',NCFLOAT,1,idim_time,ierr)
      call NCAPTC(nceq,id_sigpasma,'long_name',NCCHAR,29, &
                  'uncertainty in plasma current',ierr)
! renamed fwtcur to match MDS+ and rt-EFIT
      id_fwtpasma = NCVDEF(nceq,'fwtpasma',NCFLOAT,1,idim_time,ierr)
      call NCAPTC(nceq,id_fwtpasma,'long_name',NCCHAR,25, &
                  'weight for plasma current',ierr)
!
      id_ipmhd = NCVDEF(nceq,'cpasma',NCFLOAT,1,idim_time,ierr)
      call NCAPTC(nceq,id_ipmhd,'long_name',NCCHAR,31, &
                  'calculated plasma current (Amp)',ierr)
!
      id_chipasma = NCVDEF(nceq,'chipasma',NCFLOAT,1,idim_time,ierr)
      call NCAPTC(nceq,id_chipasma,'long_name',NCCHAR,24, &
                  'chisq for plasma current',ierr)
!
! --- pressure
!
      dim2(1) = idim_npress
      id_pressr = NCVDEF(nceq,'pressr',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_pressr,'long_name',NCCHAR,57, &
         'measured pressure vs. normalized flux (kinetic fits only)',ierr)
!
      id_rpress = NCVDEF(nceq,'rpress',NCFLOAT,2,dim2,ierr)
      presstext(1:59) = &
         '<0 - input pressure profile vs. flux; >0 - R coordinates of'
      presstext(60:85) = 'input pressure profile (m)'
      call NCAPTC(nceq,id_rpress,'long_name',NCCHAR,85, &
                  presstext,ierr)
!
      id_zpress = NCVDEF(nceq,'zpress',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_zpress,'long_name',NCCHAR,43, &
                  'Z coordinates of input pressure profile (m)',ierr)
!
      id_sigpre = NCVDEF(nceq,'sigpre',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_sigpre,'long_name',NCCHAR,24, &
                  'uncertainty for pressure',ierr)
!
      id_fwtpre = NCVDEF(nceq,'fwtpre',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_fwtpre,'long_name',NCCHAR,19, &
                  'weight for pressure',ierr)
!
      id_cpress = NCVDEF(nceq,'cpress',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_cpress,'long_name',NCCHAR,59, &
         'calculated pressure vs. normalized flux (kinetic fits only)', &
         ierr)      
!
      id_saipre = NCVDEF(nceq,'saipre',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_saipre,'long_name',NCCHAR,17, &
                  'chisq of pressure',ierr)
!
      dim2(1) = idim_npresw
      id_presw = NCVDEF(nceq,'presw',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_presw,'long_name',NCCHAR,57, &
         'measured rotational pressure vs. normalized flux (option)',ierr)
!
      id_rpresw = NCVDEF(nceq,'rpresw',NCFLOAT,2,dim2,ierr)
      preswtext(1:49) = &
         '<0 - input rotational pressure profile vs. flux; '
      preswtext(50:109) = &
         '>0 - R coordinates of input rotational pressure profile (m)'
      call NCAPTC(nceq,id_rpresw,'long_name',NCCHAR,109, &
                  preswtext,ierr)
!
      id_zpresw = NCVDEF(nceq,'zpresw',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_zpresw,'long_name',NCCHAR,54, &
         'Z coordinates of input rotational pressure profile (m)',ierr)
!
      id_sigprw = NCVDEF(nceq,'sigprw',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_sigprw,'long_name',NCCHAR,35, &
                  'uncertainty for rotational pressure',ierr)
!
      id_fwtprw = NCVDEF(nceq,'fwtprw',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_fwtprw,'long_name',NCCHAR,30, &
                  'weight for rotational pressure',ierr)
!
      id_cpresw = NCVDEF(nceq,'cpresw',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_cpresw,'long_name',NCCHAR,59, &
         'calculated rotational pressure vs. normalized flux (option)', &
         ierr)      
!
      id_saiprw = NCVDEF(nceq,'saiprw',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_saiprw,'long_name',NCCHAR,28, &
                  'chisq of rotational pressure',ierr)
!
! --- Jt/R
!
      dim2(1) = idim_kzeroj
      id_vzeroj = NCVDEF(nceq,'vzeroj',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_vzeroj,'long_name',NCCHAR,53, &
         'measured Jt/R vs. normalized flux (kinetic fits only)',ierr)
!
      id_sizeroj = NCVDEF(nceq,'sizeroj',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_sizeroj,'long_name',NCCHAR,49, &
         'normalized flux locations corresponding to vzeroj',ierr)
!
! --- quality of fit parameters
!
      dim2(1) = idim_nitera
      id_czmaxi = NCVDEF(nceq,'czmaxi',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_czmaxi,'long_name',NCCHAR,21, &
                  'Zm (cm) vs. iteration',ierr)
!
      id_cchisq = NCVDEF(nceq,'cchisq',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_cchisq,'long_name',NCCHAR,19, &
                  'chisq vs. iteration',ierr)
!
      id_cerror = NCVDEF(nceq,'cerror',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_cerror,'long_name',NCCHAR,19, &
                  'error vs. iteration',ierr)
!
      id_chifin = NCVDEF(nceq,'chifin',NCFLOAT,1,idim_time,ierr)
      call NCAPTC(nceq,id_chifin,'long_name',NCCHAR,15, &
                  'chisq upon exit',ierr)
!
      id_chitot = NCVDEF(nceq,'chitot',NCFLOAT,1,idim_time,ierr)
      call NCAPTC(nceq,id_chitot,'long_name',NCCHAR,41, &
                  'more inclusive chisq (not typically used)',ierr)
!
! --- plasma coefficients
!
      id_darea = NCVDEF(nceq,'darea',NCFLOAT,1,idim_time,ierr)
      call NCAPTC(nceq,id_darea,'long_name',NCCHAR,33, &
                  'plasma coefficients normalization',ierr)
!
      dim2(1) = idim_npcurn
      id_xrsp = NCVDEF(nceq,'xrsp',NCFLOAT,2,dim2,ierr)
      call NCAPTC(nceq,id_xrsp,'long_name',NCCHAR,19, &
                  'plasma coefficients',ierr)
!
      call NCENDF(nceq,ierr)             ! leave define mode
!-----------------------------------------------------------------------
!--   write values of variables
!-----------------------------------------------------------------------
      call NCVPT1(nceq,id_shot,1,ishot,ierr)
      endif ijump0
      m = 1
      if (iand(iout,2).ne.0) &
        m = ifirsttime
      n = ilast-ifirsttime+1
      c11(1) = 1
      c11(2) = m
      cnn(2) = n
!
      itype1: if (itype.eq.1) then
!
! --- itype = 1, called by main routine within time loop,
! --- writes variables that are time-dependent but do not have time dimension
!
      cnn(1) = 10
      call NCVPTC(nceq,id_device,c11,cnn,device,10,ierr)
      cnn(1) = nmse
      zwork(1:nmse) = real(fwtgam(1:nmse),r4)
      call NCVPT(nceq,id_fwtgam,c11,cnn,zwork,ierr)
      zwork(1:nmse) = real(chigam(1:nmse),r4)
      call NCVPT(nceq,id_chigam,c11,cnn,zwork,ierr)
      call NCVPT(nceq,id_msebkp,m,n,real(msebkp,r4),ierr)
!
      cnn(1) = nsilop
      zwork(1:nsilop) = real(sigsil,r4)
      call NCVPT(nceq,id_sigsil,c11,cnn,zwork,ierr)
      zwork(1:nsilop) = real(fwtsi,r4)
      call NCVPT(nceq,id_fwtsi,c11,cnn,zwork,ierr)
      zwork(1:nsilop) = real(saisil,r4)
      call NCVPT(nceq,id_saisil,c11,cnn,zwork,ierr)
      cnn(1) = magpri
      zwork(1:magpri) = real(sigmpi,r4)
      call NCVPT(nceq,id_sigmpi,c11,cnn,zwork,ierr)
      zwork(1:magpri) = real(fwtmp2,r4)
      call NCVPT(nceq,id_fwtmp2,c11,cnn,zwork,ierr)
      zwork(1:magpri) = real(saimpi,r4)
      call NCVPT(nceq,id_saimpi,c11,cnn,zwork,ierr)
      cnn(1) = nfsum
      zwork(1:nfsum) = real(sigfcc,r4)
      call NCVPT(nceq,id_sigfcc,c11,cnn,zwork,ierr)
      zwork(1:nfsum) = real(fwtfc,r4)
      call NCVPT(nceq,id_fwtfc,c11,cnn,zwork,ierr)
      zwork(1:nfsum) = real(chifcc,r4)
      call NCVPT(nceq,id_chifcc,c11,cnn,zwork,ierr)
      cnn(1) = nesum
      zwork(1:nesum) = real(sigecc,r4)
      call NCVPT(nceq,id_sigecc,c11,cnn,zwork,ierr)
      zwork(1:nesum) = real(fwtec,r4)
      call NCVPT(nceq,id_fwtec,c11,cnn,zwork,ierr)
      zwork(1:nesum) = real(cecurr,r4)
      call NCVPT(nceq,id_cecurr,c11,cnn,zwork,ierr)
      zwork(1:nesum) = real(chiecc,r4)
      call NCVPT(nceq,id_chiecc,c11,cnn,zwork,ierr)
      call NCVPT(nceq,id_chiecc,c11,cnn,zwork,ierr)
      call NCVPT(nceq,id_sigref,m,n,real(sigref,r4),ierr)
      call NCVPT(nceq,id_fwtref,m,n,real(fwtref,r4),ierr)
      call NCVPT(nceq,id_csiref,m,n,real(csiref,r4),ierr)
      call NCVPT(nceq,id_saisref,m,n,real(saisref,r4),ierr)
      call NCVPT(nceq,id_fwtdia,m,n,real(fwtdlc,r4),ierr)
      call NCVPT(nceq,id_chidflux,m,n,real(chidflux,r4),ierr)
      call NCVPT(nceq,id_sigpasma,m,n,real(sigpasma,r4),ierr)
      call NCVPT(nceq,id_fwtpasma,m,n,real(fwtcur,r4),ierr)
      call NCVPT(nceq,id_chipasma,m,n,real(chipasma,r4),ierr)
!
!     Note that pressr and presw have been changed to calculated pressure 
!     from the measured ones. They are switched back in measurement file(s).
!
      cnn(1) = npress1
      zwork(1:npress1) = real(premea(1:npress1),r4)
      call NCVPT(nceq,id_pressr,c11,cnn,zwork,ierr)
      zwork(1:npress1) = real(pressr(1:npress1),r4)
      call NCVPT(nceq,id_cpress,c11,cnn,zwork,ierr)
      zwork(1:npress1) = real(rpress(1:npress1),r4)
      call NCVPT(nceq,id_rpress,c11,cnn,zwork,ierr)
      zwork(1:npress1) = real(zpress(1:npress1),r4)
      call NCVPT(nceq,id_zpress,c11,cnn,zwork,ierr)
      zwork(1:npress1) = real(fwtpre(1:npress1),r4)
      call NCVPT(nceq,id_fwtpre,c11,cnn,zwork,ierr)
      zwork(1:npress1) = real(sigpre(1:npress1),r4)
      call NCVPT(nceq,id_sigpre,c11,cnn,zwork,ierr)
      zwork(1:npress1) = real(saipre2(1:npress1),r4)
      call NCVPT(nceq,id_saipre,c11,cnn,zwork,ierr)
      cnn(1) = npresw1
      zwork(1:npresw1) = real(premew(1:npresw1),r4)
      call NCVPT(nceq,id_presw,c11,cnn,zwork,ierr)
      zwork(1:npresw1) = real(presw(1:npresw1),r4)
      call NCVPT(nceq,id_cpresw,c11,cnn,zwork,ierr)
      zwork(1:npresw1) = real(rpresw(1:npresw1),r4)
      call NCVPT(nceq,id_rpresw,c11,cnn,zwork,ierr)
      zwork(1:npresw1) = real(zpresw(1:npresw1),r4)
      call NCVPT(nceq,id_zpresw,c11,cnn,zwork,ierr)
      zwork(1:npress1) = real(fwtprw(1:npress1),r4)
      call NCVPT(nceq,id_fwtprw,c11,cnn,zwork,ierr)
      zwork(1:npresw1) = real(sigprw(1:npresw1),r4)
      call NCVPT(nceq,id_sigprw,c11,cnn,zwork,ierr)
      zwork(1:npresw1) = real(saiprw2(1:npresw1),r4)
      call NCVPT(nceq,id_saiprw,c11,cnn,zwork,ierr)
!
      cnn(1) = kzeroj1
      zwork(1:kzeroj1) = real(vzeroj(1:kzeroj1),r4)
      call NCVPT(nceq,id_vzeroj,c11,cnn,zwork,ierr)
      zwork(1:kzeroj1) = real(sizeroj(1:kzeroj1),r4)
      call NCVPT(nceq,id_sizeroj,c11,cnn,zwork,ierr)
!
      cnn(1) = nitera
      ziter = real(czmaxi(1:nitera),r4)
      call NCVPT(nceq,id_czmaxi,c11,cnn,ziter,ierr)
      ziter = real(cchisq(1:nitera),r4)
      call NCVPT(nceq,id_cchisq,c11,cnn,ziter,ierr)
      ziter = real(cerror(1:nitera),r4)
      call NCVPT(nceq,id_cerror,c11,cnn,ziter,ierr)
!
      cnn(1) = kwcurn
      zwork(1:kwcurn)=real(brsp(nfsum+1:nfsum+kwcurn)/darea,r4)
      call NCVPT(nceq,id_xrsp,c11,cnn,zwork,ierr)
      call NCVPT(nceq,id_darea,m,n,real(darea,r4),ierr)
!
! --- following variables do NOT have time dimension.
!
      call NCVPT(nceq,id_mseport,1,nmse,mseport(1:nmse),ierr)
      call NCVPT(nceq,id_chifin,1,1,real(chifin,r4),ierr)
      call NCVPT(nceq,id_chitot,1,1,real(chitot,r4),ierr)
!
      endif itype1
!
      itype_iout: if (((itype.eq.1).and.(iand(iout,2).eq.0)).or. &
                      ((itype.eq.2).and.(iand(iout,2).ne.0))) then
!
! --- following variables have time dimension.
! --- if called by main routine within time loop and in individual mode,
! --- writes a slice of the variables according the time specified;
! --- if called by main routine out of time loop (itype = 2) and in 
! --- accumulative mode, writes variables in block.
!
      zwork(1:ntime) = real(time,r4)
      call NCVPT(nceq,id_time,m,n,zwork(ifirsttime),ierr)
!
      cnn(1) = nmse
      zcmgam = real(cmgam(1:nmse,:),r4)
      call NCVPT(nceq,id_cmgam,c11,cnn,zcmgam(1,ifirsttime),ierr)
!
      cnn(1) = nsilop
      zcsilop = real(csilop,r4)
      call NCVPT(nceq,id_csilop,c11,cnn,zcsilop(1,ifirsttime),ierr)
      cnn(1) = magpri
      zcmpr2 = real(cmpr2,r4)
      call NCVPT(nceq,id_cmpr2,c11,cnn,zcmpr2(1,ifirsttime),ierr)
      cnn(1) = nfsum
      zccbrsp = real(ccbrsp,r4)
      call NCVPT(nceq,id_ccbrsp,c11,cnn,zccbrsp(1,ifirsttime),ierr)
      zwork(1:ntime) = real(psiref,r4)
      call NCVPT(nceq,id_psiref,m,n,zwork(ifirsttime),ierr)
      zwork(1:ntime) = real(diamag,r4)
      call NCVPT(nceq,id_diamag,m,n,zwork(ifirsttime),ierr)
      zwork(1:ntime) = real(cdflux,r4)
      call NCVPT(nceq,id_cdflux,m,n,zwork(ifirsttime),ierr)
      zwork(1:ntime) = real(sigdia,r4)
      call NCVPT(nceq,id_sigdia,m,n,zwork(ifirsttime),ierr)
!
      zwork(1:ntime) = real(ipmeas,r4)
      call NCVPT(nceq,id_plasma,m,n,zwork(ifirsttime),ierr)
      zwork(1:ntime) = real(ipmhd,r4)
      call NCVPT(nceq,id_ipmhd,m,n,zwork(ifirsttime),ierr)
!
      zwork(1:ntime) = real(curc79,r4)
      call NCVPT(nceq,id_curc79,m,n,zwork(ifirsttime),ierr)
      zwork(1:ntime) = real(curc139,r4)
      call NCVPT(nceq,id_curc139,m,n,zwork(ifirsttime),ierr)
      zwork(1:ntime) = real(curc199,r4)
      call NCVPT(nceq,id_curc199,m,n,zwork(ifirsttime),ierr)
      zwork(1:ntime) = real(curiu30,r4)
      call NCVPT(nceq,id_curiu30,m,n,zwork(ifirsttime),ierr)
      zwork(1:ntime) = real(curil30,r4)
      call NCVPT(nceq,id_curil30,m,n,zwork(ifirsttime),ierr)
      zwork(1:ntime) = real(curiu90,r4)
      call NCVPT(nceq,id_curiu90,m,n,zwork(ifirsttime),ierr)
      zwork(1:ntime) = real(curil90,r4)
      call NCVPT(nceq,id_curil90,m,n,zwork(ifirsttime),ierr)
      zwork(1:ntime) = real(curiu150,r4)
      call NCVPT(nceq,id_curiu150,m,n,zwork(ifirsttime),ierr)
      zwork(1:ntime) = real(curil150,r4)
      call NCVPT(nceq,id_curil150,m,n,zwork(ifirsttime),ierr)
!
! --- following data and their corresponding ncdf variables have 
! --- reversed dimensions
!
      call NCVPT1(nceq,id_fungam,1,msefitfun,ierr)
!
      stride(1) = 1
      stride(2) = 1
      imap(2) = 4                         ! number of bytes in float
!     imap(2) = imap(2)*2                 ! number of bytes in double
      imap(1) = imap(2)*ntime
!
      cnn(1) = nmse
      zstark = real(tangam(:,1:nmse),r4)
      call NCVPTG(nceq,id_tangam,c11,cnn,stride,imap, &
                  zstark(ifirsttime,1),ierr)
      zstark = real(tangam_uncor(:,1:nmse),r4)
      call NCVPTG(nceq,id_tangam_uncor,c11,cnn,stride,imap, &
                  zstark(ifirsttime,1),ierr)
      do i=1,nmse
        zstark(:,i) = real(spatial_fix(i,:),r4)
      enddo
      call NCVPTG(nceq,id_fgam,c11,cnn,stride,imap, &
                  zstark(ifirsttime,1),ierr)
      zstark = real(siggam(:,1:nmse),r4)
      call NCVPTG(nceq,id_siggam,c11,cnn,stride,imap, &
                  zstark(ifirsttime,1),ierr)
      zstark = real(rrgam(:,1:nmse),r4)
      call NCVPTG(nceq,id_rrgam,c11,cnn,stride,imap, &
                  zstark(ifirsttime,1),ierr)
      zstark = real(zzgam(:,1:nmse),r4)
      call NCVPTG(nceq,id_zzgam,c11,cnn,stride,imap, &
                  zstark(ifirsttime,1),ierr)
      zstark = real(a1gam(:,1:nmse),r4)
      call NCVPTG(nceq,id_a1gam,c11,cnn,stride,imap, &
                  zstark(ifirsttime,1),ierr)
      zstark = real(a2gam(:,1:nmse),r4)
      call NCVPTG(nceq,id_a2gam,c11,cnn,stride,imap, &
                  zstark(ifirsttime,1),ierr)
      zstark = real(a3gam(:,1:nmse),r4)
      call NCVPTG(nceq,id_a3gam,c11,cnn,stride,imap, &
                  zstark(ifirsttime,1),ierr)
      zstark = real(a4gam(:,1:nmse),r4)
      call NCVPTG(nceq,id_a4gam,c11,cnn,stride,imap, &
                  zstark(ifirsttime,1),ierr)
      zstark = real(a5gam(:,1:nmse),r4)
      call NCVPTG(nceq,id_a5gam,c11,cnn,stride,imap, &
                  zstark(ifirsttime,1),ierr)
      zstark = real(a6gam(:,1:nmse),r4)
      call NCVPTG(nceq,id_a6gam,c11,cnn,stride,imap, &
                  zstark(ifirsttime,1),ierr)
      zstark = real(a7gam(:,1:nmse),r4)
      call NCVPTG(nceq,id_a7gam,c11,cnn,stride,imap, &
                  zstark(ifirsttime,1),ierr)
      zstark = real(a8gam(:,1:nmse),r4)
      call NCVPTG(nceq,id_a8gam,c11,cnn,stride,imap, &
                  zstark(ifirsttime,1),ierr)
      do i=1,ntime
        zstark(i,:) = real(rmse_gain(1:nmse),r4)
      enddo
      call NCVPTG(nceq,id_gaingam,c11,cnn,stride,imap, &
                  zstark(ifirsttime,1),ierr)
      do i=1,ntime
        zstark(i,:) = real(rmse_slope(1:nmse),r4)
      enddo
      call NCVPTG(nceq,id_slopegam,c11,cnn,stride,imap, &
                  zstark(ifirsttime,1),ierr)
      do i=1,ntime
        zstark(i,:) = real(rmse_scale(1:nmse),r4)
      enddo
      call NCVPTG(nceq,id_scalegam,c11,cnn,stride,imap, &
                  zstark(ifirsttime,1),ierr)
      do i=1,ntime
        zstark(i,:) = real(rmse_offset(1:nmse),r4)
      enddo
      call NCVPTG(nceq,id_offsetgam,c11,cnn,stride,imap, &
                  zstark(ifirsttime,1),ierr)
!
      cnn(1) = nsilop
      zsilopt = real(silopt,r4)
      call NCVPTG(nceq,id_silopt,c11,cnn,stride,imap, &
                  zsilopt(ifirsttime,1),ierr)
      cnn(1) = magpri
      zexpmpi = real(expmpi,r4)
      call NCVPTG(nceq,id_expmpi,c11,cnn,stride,imap, &
                  zexpmpi(ifirsttime,1),ierr)
      cnn(1) = nfsum
      zfccurt = real(fccurt,r4)
      call NCVPTG(nceq,id_fccurt,c11,cnn,stride,imap, &
                  zfccurt(ifirsttime,1),ierr)
      cnn(1) = nesum
      zeccurt = real(eccurt,r4)
      call NCVPTG(nceq,id_eccurt,c11,cnn,stride,imap, &
                  zeccurt(ifirsttime,1),ierr)
!
      cnn(1) = nacoil
      zaccurt = real(accurt,r4)
      call NCVPTG(nceq,id_accurt,c11,cnn,stride,imap, &
                  zaccurt(ifirsttime,1),ierr)
      zaccurt = real(caccurt,r4)
      call NCVPTG(nceq,id_caccurt,c11,cnn,stride,imap, &
                  zaccurt(ifirsttime,1),ierr)
!
      call NCCLOS(nceq,ierr)             ! close the file
      endif itype_iout
!
      return
      end subroutine write_m
