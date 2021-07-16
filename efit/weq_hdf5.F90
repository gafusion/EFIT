!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          Subroutine wmeasure writes measurements (MSE, magnetic) **
!**          and some of their corresponding calculated values in    **
!**          a netCDF file of name m0sssss.nc/m0sssss.0tttt.         **
!**                                                                  **
!**          The input variable IOUT decides the mode,               **
!**          IOUT.and.2.ne.0 - one file for all slices m0sssss.nc    **
!**          IOUT.and.4.ne.0 - one file for each slice m0sssss.0tttt **
!**          IOUT.and.6.ne.0 - one file for all m0sssss.0tttt_nnnn   **
!**               tttt as the first time,nnnn as the number of slices**
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**          ktime - total number of time slices                     **
!**          ifirsttime- index of starting time                          **
!**          ilast - index of ending time                            **
!**          itype - 1 called from main routine with time loop       **
!**                  2 called from main routine out of time loop,    **
!**                    write all slices at one time                  **
!**                                                                  **
!**     REFERENCES:                                                  **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          12-10-97 created by Q.Peng                              **
!**          15-06-98 Q.Peng added measured pressure as premea.      **
!**          01-11-98 Q.Peng added mseport (vport in stark_multi)    **
!**          01-20-2000 Q.P. following the convention of A & G files **
!**                          use mssssss.0tttt_nnn for sub-millisec. **
!**                                                                  **
!**********************************************************************
      subroutine wmeasure(ktime,ifirsttime,ilast,itype)
      use set_kinds
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!      include 'ecomdu1.f90'
      include 'netcdf.inc'   ! from the netCDF package..
!                            ..this must be symlinked to local directory
! --- common block ids is local to wmeasure
!
      common/ids/nceq,idim_time,idim_1,idim_nstark,idim_nsilop &
      ,idim_magpri,idim_nfcoil,idim_nesum,idim_npress,idim_npcurn &
      ,idim_nitera &
      ,id_shot,id_time,id_tangam,id_cmgam,id_siggam,id_rrgam,id_zzgam &
      ,id_a1gam,id_a2gam,id_a3gam,id_a4gam,id_a5gam,id_a6gam,id_a7gam &
      ,id_a8gam,id_fgam, id_tangam_uncor &
      ,id_fwtgam,id_msebkp,id_mseport &
      ,id_silopt,id_csilop,id_fwtsi,id_expmpi,id_fwtmp2,id_cmpr2 &
      ,id_fccurt,id_ccbrsp,id_fwtfc,id_eccurt,id_cecurr,id_fwtec &
      ,id_diamag,id_cdflux,id_sigdia,id_plasma,id_cpasma &
      ,id_pressr,id_cpress,id_rpress,id_zpress,id_sigpre,id_saipre &
      ,id_chigam,id_saimpi,id_saisil,id_czmaxi,id_cchisq,id_cerror &
      ,id_fungam,id_gaingam,id_slopegam,id_scalegam,id_offsetgam
      dimension xrsp(npcurn)
!sri  to make it compatiable to mpi version
!      character eqdsk*20,let,title*80,last*3
      character let,title*80
      character(len=4) last
      character(len=80) eqdsk
      character*30 sfname
      integer dim2(2),c11(2),cnn(2),imap(2),stride(2)
! --- temporay variables to convert double to single
      real*4 zwork(ntime+nsilop+nstark+nfcoil+nesum+magpri+npress) &
      ,zcmgam(nstark,ntime),zcsilop(nsilop,ntime) &
      ,zcmpr2(magpri,ntime),zccbrsp(nfcoil,ntime),zstark(ntime,nstark) &
      ,zsilopt(ntime,nsilop),zexpmpi(ntime,magpri) &
      ,zfccurt(ntime,nfcoil),zeccurt(ntime,nesum) 
!vas f90 modifi
      character*85 vastext
!-----------------------------------------------------------------------
!--  write out t(shot).(time)_X files                                    --
!-----------------------------------------------------------------------
      if ((kwripre.eq.11).and.(itype.eq.2)) then
        call getfnmd('t',ishot,itime,sfname)
        sfname=sfname(1:13)//'_chi2'
          open(unit=74,status='old',file=sfname,err=19920)
          close(unit=74,status='delete')
19920     continue
          open(unit=74,status='new',file=sfname                       )
        do i=ifirsttime,ilast
          write (74,*) time(i),tsaisq(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_error'
          open(unit=74,status='old',file=sfname,err=19922)
          close(unit=74,status='delete')
19922     continue
          open(unit=74,status='new',file=sfname                       )
        do i=ifirsttime,ilast
          write (74,*) time(i),terror(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_j1ave'
          open(unit=74,status='old',file=sfname,err=19930)
          close(unit=74,status='delete')
19930     continue
          open(unit=74,status='new',file=sfname                       )
        do i=ifirsttime,ilast
          write (74,*) time(i),cj1ave(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_li'
          open(unit=74,status='old',file=sfname,err=19940)
          close(unit=74,status='delete')
19940     continue
          open(unit=74,status='new',file=sfname                       )
        do i=ifirsttime,ilast
          write (74,*) time(i),ali(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_betat'
          open(unit=74,status='old',file=sfname,err=19950)
          close(unit=74,status='delete')
19950     continue
          open(unit=74,status='new',file=sfname                       )
        do i=ifirsttime,ilast
          write (74,*) time(i),betat(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_q95'
          open(unit=74,status='old',file=sfname,err=19960)
          close(unit=74,status='delete')
19960     continue
          open(unit=74,status='new',file=sfname                       )
        do i=ifirsttime,ilast
          write (74,*) time(i),qpsib(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_q0'
          open(unit=74,status='old',file=sfname,err=19970)
          close(unit=74,status='delete')
19970     continue
          open(unit=74,status='new',file=sfname                       )
        do i=ifirsttime,ilast
          write (74,*) time(i),qqmagx(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_q0'
          open(unit=74,status='old',file=sfname,err=19980)
          close(unit=74,status='delete')
19980     continue
          open(unit=74,status='new',file=sfname                       )
        do i=ifirsttime,ilast
          write (74,*) time(i),qqmagx(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_eout'
          open(unit=74,status='old',file=sfname,err=19990)
          close(unit=74,status='delete')
19990     continue
          open(unit=74,status='new',file=sfname                       )
        do i=ifirsttime,ilast
          write (74,*) time(i),eout(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_vout'
          open(unit=74,status='old',file=sfname,err=19992)
          close(unit=74,status='delete')
19992     continue
          open(unit=74,status='new',file=sfname                       )
        do i=ifirsttime,ilast
          vm3=vout(i)/1.e6_dp
          write (74,*) time(i),vm3,xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_betan'
          open(unit=74,status='old',file=sfname,err=19995)
          close(unit=74,status='delete')
19995     continue
          open(unit=74,status='new',file=sfname                       )
        do i=ifirsttime,ilast
          pasman=cpasma(i)/1.e4_dp/aout(i)/abs(bcentr(i))
          pasman=pasman*rout(i)/100./rcentr
          betatnx=betat(i)/pasman
          write (74,*) time(i),betatnx,xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_zts'
          open(unit=74,status='old',file=sfname,err=19998)
          close(unit=74,status='delete')
19998     continue
          open(unit=74,status='new',file=sfname                       )
        do i=ifirsttime,ilast
          write (74,*) time(i),zuperts(i),xdum,xdum
        enddo
        close(unit=74)
      endif
!
      if ((iand(iout,2).eq.0).and.(iand(iout,4).eq.0)) &
           return
!
      call NCPOPT (NCVERBOS)              ! get nc error message but not fatal
!
! --- time-dependent and the first slice and first time called, 
! --- create file and enter define mode
! --- if iout.and.2.ne.0 = create m0sssss.nc
! --- if iout.and.6.ne.0 = create m0sssss.0tttt_nnnn, 
! ---    while tttt is the first time and nnn is the number of slices.

      if ((iand(iout,2).ne.0).and.(iand(iout,4).ne.0).and. &
           ((ifirsttime.eq.1).and.(itype.eq.1))) then
         let = 'm'
         iitime = time(ifirsttime)
         call getfnmd(let,ishot,iitime,eqdsk)
         if (ktime.le.9)                       write(unit=last,fmt=1030) ktime
         if ((ktime.gt.9).and.(ktime.le.99))   write(unit=last,fmt=1040) ktime
         if ((ktime.gt.99).and.(ktime.le.999)) write(unit=last,fmt=1050) ktime
         if (ktime.gt.999)                     write(unit=last,fmt=1060) ktime
 1030    format ('000',i1)
 1040    format ('00',i2)
 1050    format ('0',i3)
 1060    format (i4)
         eqdsk = eqdsk(1:13)//'_'//last
         nceq = NCCRE(eqdsk,NCCLOB,ierr)              
!
      elseif ((iand(iout,2).ne.0).and.(iand(iout,4).eq.0).and. &
              (ifirsttime.eq.1).and.(itype.eq.1)) then
         if (ishot.le.99999) write(unit=eqdsk,fmt=1010) ishot
         if (ishot.gt.99999) write(unit=eqdsk,fmt=1020) ishot
 1010    format ('m0',i5,'.nc')
 1020    format ('m',i6,'.nc')
         nceq = NCCRE(eqdsk,NCCLOB,ierr) ! create file, overwrite if exists
!
! --- creates one file for each slice
!
      elseif ((iand(iout,4).ne.0).and.(iand(iout,2).eq.0).and. &
              (itype.eq.1)) then
         let = 'm'
         iitime = time(ifirsttime)
         call getfnmu(itimeu,let,ishot,iitime,eqdsk)
         nceq = NCCRE(eqdsk,NCCLOB,ierr)  
!
! --- time-dependent but NOT the first slice NOR the first time,
! --- or single slice but NOT the first time, 
! --- skip define mode and go to write directly
!
      else
         go to 100
      endif
      
      if (ierr.ne.0) return
!
! --- nc file has been created successfully, now define variables, etc.
!
      title = 'EFIT measurement file Mssssss.ttttt/Mssssss.nc'
      call NCAPTC (nceq,NC_GLOBAL,'title',NCCHAR,46,title,ierr)
! 
! --- define unlimited time dimension and scalar and array dimensions
!
      idim_time = NCDDEF (nceq,'dim_time',NCUNLIM,ierr)
      idim_1 = NCDDEF (nceq,'dim_scalar',1,ierr)
      idim_nstark = NCDDEF (nceq,'dim_nstark',nstark,ierr)
      idim_nsilop = NCDDEF (nceq,'dim_nsilop',nsilop,ierr)
      idim_magpri = NCDDEF (nceq,'dim_magpri',magpri,ierr)
      idim_nfcoil = NCDDEF (nceq,'dim_nfcoil',nfcoil,ierr)
      idim_nesum  = NCDDEF (nceq,'dim_nesum' ,nesum ,ierr)
      npress1 = npress
      if (npress.eq.0) npress1 = 1
      idim_npress = NCDDEF (nceq,'dim_npress',npress1,ierr)
      idim_npcurn = NCDDEF (nceq,'dim_npcurn',npcurn,ierr)
      idim_nitera = NCDDEF (nceq,'dim_nitera',nitera,ierr)
      dim2(2) = idim_time
!
! --- define variables
!
      id_shot = NCVDEF (nceq,'shot',NCLONG,1,idim_1,ierr)
      call NCAPTC (nceq,id_shot,'long_name',NCCHAR,11, &
                  'shot number',ierr)
!
      id_time = NCVDEF (nceq,'time',NCFLOAT,1,idim_time,ierr)
      call NCAPTC (nceq,id_time,'units',NCCHAR,4,'msec',ierr)
!
! --- MSE measurements
!
      id_fungam = NCVDEF (nceq,'msefitfun',NCLONG,1,idim_1,ierr)
      call NCAPTC (nceq,id_fungam,'long_name',NCCHAR,16, &
                  'MSE fit function',ierr)

      dim2(1) = idim_nstark
      id_tangam = NCVDEF (nceq,'tangam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_tangam,'long_name',NCCHAR,35, &
      'tangent of measured MSE pitch angle',ierr)
!
      id_tangam_uncor = NCVDEF (nceq,'tangam_uncor',NCFLOAT,&
                             2,dim2,ierr)
      call NCAPTC (nceq,id_tangam_uncor,'long_name',NCCHAR,54, &
      'tangent of measured MSE pitch angle w/o cer correction',ierr)
!
      id_fgam = NCVDEF (nceq,'fixgam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_fgam,'long_name',NCCHAR,58, &
      'radians correction of tangam for spatial averaging effects',ierr)
!
      id_cmgam = NCVDEF (nceq,'cmgam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_cmgam,'long_name',NCCHAR,30, &
      'calculated polarimetry signals',ierr)
!
      id_siggam = NCVDEF (nceq,'siggam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_siggam,'long_name',NCCHAR,21, &
      'uncertainty of tangam',ierr)
!
      id_rrgam = NCVDEF (nceq,'rrgam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_rrgam,'long_name',NCCHAR,22, &
      'radius of MSE channels',ierr)
!
      id_zzgam = NCVDEF (nceq,'zzgam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_zzgam,'long_name',NCCHAR,26, &
      'Z position of MSE channels',ierr)
!
      id_a1gam = NCVDEF (nceq,'a1gam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_a1gam,'long_name',NCCHAR,45, &
      'viewing geometry coefficients of MSE channels',ierr)
!
      id_a2gam = NCVDEF (nceq,'a2gam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_a2gam,'long_name',NCCHAR,45, &
      'viewing geometry coefficients of MSE channels',ierr)
!
      id_a3gam = NCVDEF (nceq,'a3gam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_a3gam,'long_name',NCCHAR,45, &
      'viewing geometry coefficients of MSE channels',ierr)
!
      id_a4gam = NCVDEF (nceq,'a4gam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_a4gam,'long_name',NCCHAR,45, &
      'viewing geometry coefficients of MSE channels',ierr)
!
      id_a5gam = NCVDEF (nceq,'a5gam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_a5gam,'long_name',NCCHAR,45, &
      'viewing geometry coefficients of MSE channels',ierr)
!
      id_a6gam = NCVDEF (nceq,'a6gam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_a6gam,'long_name',NCCHAR,45, &
      'viewing geometry coefficients of MSE channels',ierr)
!
      id_a7gam = NCVDEF (nceq,'a7gam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_a7gam,'long_name',NCCHAR,45, &
      'viewing geometry coefficients of MSE channels',ierr)
!
      id_a8gam = NCVDEF (nceq,'a8gam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_a8gam,'long_name',NCCHAR,45, &
      'viewing geometry coefficients of MSE channels',ierr)
!
      id_fwtgam = NCVDEF (nceq,'fwtgam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_fwtgam,'long_name',NCCHAR,31, &
      'fitting weight for MSE channels',ierr)
!
      id_msebkp = NCVDEF (nceq,'msebkp',NCFLOAT,1,idim_time,ierr)
      call NCAPTC (nceq,id_msebkp,'long_name',NCCHAR,30, &
      'background substraction switch',ierr)
!     
      id_mseport = NCVDEF (nceq,'vport',NCLONG,1,idim_nstark,ierr)
      call NCAPTC (nceq,id_mseport,'long_name',NCCHAR,33, &
      'view port locations of MSE system',ierr)
!
      if( msefitfun .eq. 1)then

          id_gaingam = NCVDEF (nceq,'mcal_gain',NCFLOAT,2,dim2,ierr)
          call NCAPTC (nceq,id_gaingam,'long_name',NCCHAR,38, &
      'gain param for tangent offset function',ierr)

          id_slopegam = NCVDEF (nceq,'mcal_slope',NCFLOAT, &
                              2,dim2,ierr)
          call NCAPTC (nceq,id_slopegam,'long_name',NCCHAR,39, &
      'slope param for tangent offset function',ierr)

          id_scalegam = NCVDEF (nceq,'mcal_scale',NCFLOAT, &
                              2,dim2,ierr)
          call NCAPTC (nceq,id_scalegam,'long_name',NCCHAR,39, &
      'scale param for tangent offset function',ierr)

          id_offsetgam = NCVDEF (nceq,'mcal_offset',NCFLOAT, &
                              2,dim2,ierr)
          call NCAPTC (nceq,id_offsetgam,'long_name',NCCHAR,40, &
      'offset param for tangent offset function',ierr)


      else if( msefitfun .eq. 3)then

          id_gaingam = NCVDEF (nceq,'mcal3_gain',NCFLOAT,2,dim2,ierr)
          call NCAPTC (nceq,id_gaingam,'long_name',NCCHAR,37, &
      'gain param for tangent slope function',ierr)

          id_slopegam = NCVDEF (nceq,'mcal3_phase',NCFLOAT,2,dim2,ierr)
          call NCAPTC (nceq,id_slopegam,'long_name',NCCHAR,38, &
      'phase param for tangent slope function',ierr)

          id_scalegam = NCVDEF (nceq,'mcal3_btscale',NCFLOAT, &
                              2,dim2,ierr)
          call NCAPTC (nceq,id_scalegam,'long_name',NCCHAR,40, &
      'btscale param for tangent slope function',ierr)

          id_offsetgam = NCVDEF (nceq,'mcal3_dc_offset',NCFLOAT, &
                              2,dim2,ierr)
          call NCAPTC (nceq,id_offsetgam,'long_name',NCCHAR,42, &
      'dc_offset param for tangent slope function',ierr)

      else if( msefitfun .eq. 4)then

          id_gaingam = NCVDEF (nceq,'mcal4_gain',NCFLOAT,2,dim2,ierr)
          call NCAPTC (nceq,id_gaingam,'long_name',NCCHAR,37, &
      'gain param for tangent slope function',ierr)

          id_slopegam = NCVDEF (nceq,'mcal4_phase',NCFLOAT,2,dim2,ierr)
          call NCAPTC (nceq,id_slopegam,'long_name',NCCHAR,38, &
      'phase param for tangent slope function',ierr)

          id_scalegam = NCVDEF (nceq,'mcal4_btscale',NCFLOAT, &
                              2,dim2,ierr)
          call NCAPTC (nceq,id_scalegam,'long_name',NCCHAR,40, &
      'btscale param for tangent slope function',ierr)

          id_offsetgam = NCVDEF (nceq,'mcal4_dc_offset',NCFLOAT, &
                              2,dim2,ierr)
          call NCAPTC (nceq,id_offsetgam,'long_name',NCCHAR,42, &
      'dc_offset param for tangent slope function',ierr)

      endif
!
! --- Magnetic measurements
!
      dim2(1) = idim_nsilop
      id_silopt = NCVDEF (nceq,'silopt',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_silopt,'long_name',NCCHAR,19, &
      'measured flux loops',ierr)
!
      id_csilop = NCVDEF (nceq,'csilop',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_csilop,'long_name',NCCHAR,21, &
      'calculated flux loops',ierr)
!
      id_fwtsi = NCVDEF (nceq,'fwtsi',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_fwtsi,'long_name',NCCHAR,21, &
      'weight for flux loops',ierr)
!
      dim2(1) = idim_magpri
      id_expmpi = NCVDEF (nceq,'expmpi',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_expmpi,'long_name',NCCHAR,24, &
      'measured magnetic probes',ierr)
!
      id_fwtmp2 = NCVDEF (nceq,'fwtmp2',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_fwtmp2,'long_name',NCCHAR,26, &
      'weight for magnetic probes',ierr)
!
      id_cmpr2 = NCVDEF (nceq,'cmpr2',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_cmpr2,'long_name',NCCHAR,26, &
      'calculated magnetic probes',ierr)
!
      dim2(1) = idim_nfcoil
      id_fccurt = NCVDEF (nceq,'fccurt',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_fccurt,'long_name',NCCHAR,30, &
      'measured F-coil currents (Amp)',ierr)
!
      id_ccbrsp = NCVDEF (nceq,'ccbrsp',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_ccbrsp,'long_name',NCCHAR,32, &
      'calculated F-coil currents (Amp)',ierr)
!
      id_fwtfc = NCVDEF (nceq,'fwtfc',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_fwtfc,'long_name',NCCHAR,26, &
      'weight for F-coil currents',ierr)
!
      dim2(1) = idim_nesum
      id_eccurt = NCVDEF (nceq,'eccurt',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_eccurt,'long_name',NCCHAR,30, &
      'measured E-coil currents (Amp)',ierr)
!
      id_cecurr = NCVDEF (nceq,'cecurr',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_cecurr,'long_name',NCCHAR,32, &
      'calculated E-coil currents (Amp)',ierr)
!
      id_fwtec = NCVDEF (nceq,'fwtec',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_fwtec,'long_name',NCCHAR,26, &
      'weight for E-coil currents',ierr)
!
      id_diamag = NCVDEF (nceq,'diamag',NCFLOAT,1,idim_time,ierr)
      call NCAPTC (nceq,id_diamag,'long_name',NCCHAR,25, &
      'measured diamagnetic flux',ierr)
!
      id_cdflux = NCVDEF (nceq,'cdflux',NCFLOAT,1,idim_time,ierr)
      call NCAPTC (nceq,id_cdflux,'long_name',NCCHAR,27, &
      'calculated diamagnetic flux',ierr)
!
      id_sigdia = NCVDEF (nceq,'sigdia',NCFLOAT,1,idim_time,ierr)
      call NCAPTC (nceq,id_sigdia,'long_name',NCCHAR,31, &
      'uncertainty of diamagnetic flux',ierr)
!
! --- plasma current
!
      id_plasma = NCVDEF (nceq,'plasma',NCFLOAT,1,idim_time,ierr)
      call NCAPTC (nceq,id_plasma,'long_name',NCCHAR,29, &
      'measured plasma current (Amp)',ierr)
!
      id_cpasma = NCVDEF (nceq,'cpasma',NCFLOAT,1,idim_time,ierr)
      call NCAPTC (nceq,id_cpasma,'long_name',NCCHAR,31, &
      'calculated plasma current (Amp)',ierr)
!
      id_darea = NCVDEF (nceq,'darea',NCFLOAT,1,idim_time,ierr)
      call NCAPTC (nceq,id_darea,'long_name',NCCHAR,33, &
      'plasma coefficients normalization',ierr)
!
      dim2(1) = idim_npcurn
      id_xrsp = NCVDEF (nceq,'xrsp',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_xrsp,'long_name',NCCHAR,19, &
      'plasma coefficients',ierr)
!
! --- pressure
!
      dim2(1) = idim_npress
      id_pressr = NCVDEF (nceq,'pressr',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_pressr,'long_name',NCCHAR,57, &
      'measured pressure vs. normalized flux (kinetic fits only)',ierr)      
!
      id_cpress = NCVDEF (nceq,'cpress',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_cpress,'long_name',NCCHAR,59, &
      'calculated pressure vs. normalized flux (kinetic fits only)', &
      ierr)      
!
      id_rpress = NCVDEF (nceq,'rpress',NCFLOAT,2,dim2,ierr)
!vas f90 modifi.
vastext(1:59) = '<0 - input pressure profile vs. flux; >0 - R coordinates of'
vastext(60:85)= 'input pressure profile (m)'

      call NCAPTC (nceq,id_rpress,'long_name',NCCHAR,86, &
!vas      '<0 - input pressure profile vs. flux; >0 - R coordinates of input &
!vas       pressure profile (m)',ierr)
      vastext,ierr)
!
      id_zpress = NCVDEF (nceq,'zpress',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_zpress,'long_name',NCCHAR,43, &
      'Z coordinates of input pressure profile (m)',ierr)
!
      id_sigpre = NCVDEF (nceq,'sigpre',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_sigpre,'long_name',NCCHAR,24, &
      'uncertainty for pressure',ierr)
!
      id_saipre = NCVDEF (nceq,'saipre',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_saipre,'long_name',NCCHAR,17, &
      'chisq of pressure',ierr)

!
! --- fitting parameters
!
      dim2(1) = idim_nstark
      id_chigam = NCVDEF (nceq,'chigam',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_chigam,'long_name',NCCHAR,23, &
      'chisq vs. polarimetries',ierr)
!
      dim2(1) = idim_magpri
      id_saimpi = NCVDEF (nceq,'saimpi',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_saimpi,'long_name',NCCHAR,25, &
      'chisq vs. magnetic probes',ierr)
!
      dim2(1) = idim_nsilop
      id_saisil = NCVDEF (nceq,'saisil',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_saisil,'long_name',NCCHAR,19, &
      'chisq vs. PSI loops',ierr)
!
      dim2(1) = idim_nitera
      id_czmaxi = NCVDEF (nceq,'czmaxi',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_czmaxi,'long_name',NCCHAR,21, &
      'Zm (cm) vs. iteration',ierr)
!
      id_cchisq = NCVDEF (nceq,'cchisq',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_cchisq,'long_name',NCCHAR,19, &
      'chisq vs. iteration',ierr)
!
      id_cerror = NCVDEF (nceq,'cerror',NCFLOAT,2,dim2,ierr)
      call NCAPTC (nceq,id_cerror,'long_name',NCCHAR,19, &
      'error vs. iteration',ierr)
!
      call NCENDF (nceq,ierr)             ! leave define mode
!
! --- write values of variables
!
      call NCVPT1 (nceq,id_shot,1,ishot,ierr)
 100  continue
      m = 1
      if (iand(iout,2).ne.0) &
      m = ifirsttime
      n = ilast-ifirsttime+1
      c11(1) = 1
      c11(2) = m
      cnn(2) = n
!
      if (itype.eq.1) then
!
! --- itype = 1, called by main routine within time loop,
! --- writes variables that are time-dependent but do not have time dimension
!
      cnn(1) = nstark
      do j=1, nstark
         zwork(j) = real(fwtgam(j))
      enddo
      call NCVPT (nceq,id_fwtgam,c11,cnn,zwork,ierr)

      call NCVPT (nceq,id_msebkp,m,n,real(msebkp),ierr)
!
      cnn(1) = nsilop
      do j=1, nsilop
         zwork(j) = real(fwtsi(j))
      enddo
      call NCVPT (nceq,id_fwtsi,c11,cnn,zwork,ierr)
      cnn(1) = magpri
      do j=1, magpri
         zwork(j) = real(fwtmp2(j))
      enddo
      call NCVPT (nceq,id_fwtmp2,c11,cnn,zwork,ierr)
      cnn(1) = nfcoil
      do j=1, nfcoil
         zwork(j) = real(fwtfc(j))
      enddo
      call NCVPT (nceq,id_fwtfc,c11,cnn,zwork,ierr)
      cnn(1) = nesum
      do j=1, nesum
         zwork(j) = real(fwtec(j))
      enddo
      call NCVPT (nceq,id_fwtec,c11,cnn,zwork,ierr)
!
!     Note that pressr has been changed to calculated pressure from the 
!     measured one. It is switched back in measurement file.
!
      cnn(1) = npress1
      do j=1, npress1
         zwork(j) = real(premea(j))
      enddo
      call NCVPT (nceq,id_pressr,c11,cnn,zwork,ierr)
      do j=1, npress1
         zwork(j) = real(pressr(j))
      enddo
      call NCVPT (nceq,id_cpress,c11,cnn,zwork,ierr)
      do j=1, npress1
         zwork(j) = real(rpress(j))
      enddo
      call NCVPT (nceq,id_rpress,c11,cnn,zwork,ierr)
      do j=1, npress1
         zwork(j) = real(zpress(j))
      enddo
      call NCVPT (nceq,id_zpress,c11,cnn,zwork,ierr)
      do j=1, npress1
         zwork(j) = real(sigpre(j))
      enddo
      call NCVPT (nceq,id_sigpre,c11,cnn,zwork,ierr)
      do j=1, npress1
         zwork(j) = real(saipre2(j))
      enddo
      call NCVPT (nceq,id_saipre,c11,cnn,zwork,ierr)
!
      cnn(1) = nesum
      do j=1, nesum
         zwork(j) = real(cecurr(j))
      enddo
      call NCVPT (nceq,id_cecurr,c11,cnn,zwork,ierr)
      cnn(1) = kwcurn
      do i=1,kwcurn
         zwork(i)=brsp(nfcoil+i)/darea
      enddo
      call NCVPT (nceq,id_darea,m,n,real(darea),ierr)
      call NCVPT (nceq,id_xrsp,c11,cnn,zwork,ierr)
!
      cnn(1) = nstark
      do i=1, nstark
         zwork(i) = real(chigam(i))
      enddo
      call NCVPT (nceq,id_chigam,c11,cnn,zwork,ierr)
      cnn(1) = magpri
      do i=1, magpri
         zwork(i) = real(saimpi(i))
      enddo
      call NCVPT (nceq,id_saimpi,c11,cnn,zwork,ierr)
      cnn(1) = nsilop
      do i=1, nsilop
         zwork(i) = real(saisil(i))
      enddo
      call NCVPT (nceq,id_saisil,c11,cnn,zwork,ierr)
!
      cnn(1) = nitera
      do i=1, nitera
         zwork(i) = real(czmaxi(i))
      enddo
      call NCVPT (nceq,id_czmaxi,c11,cnn,zwork,ierr)
      do i=1, nitera
         zwork(i) = real(cchisq(i))
      enddo
      call NCVPT (nceq,id_cchisq,c11,cnn,zwork,ierr)
      do i=1, nitera
         zwork(i) = real(cerror(i))
      enddo
      call NCVPT (nceq,id_cerror,c11,cnn,zwork,ierr)
!
! --- following variables do NOT have time dimension.
!
      call NCVPT (nceq,id_mseport,1,nstark,mseport,ierr)
      endif
!
      if (((itype.eq.1).and.(iand(iout,2).eq.0)).or. &
          ((itype.eq.2).and.(iand(iout,2).ne.0))) then
!
! --- following variables have time dimension.
! --- if called by main routine within time loop and in individual mode,
! --- writes a slice of the variables according the time specified;
! --- if called by main routine out of time loop (itype = 2) and in 
! --- accumulative mode, writes variables in block.
!
      do i=1, ntime
         zwork(i) = real(time(i))
      enddo
      call NCVPT (nceq,id_time,m,n,zwork(ifirsttime),ierr)
      cnn(1) = nstark
      do i=1, nstark
      do j=1, ntime
         zcmgam(i,j) = real(cmgam(i,j))
      enddo
      enddo
      call NCVPT (nceq,id_cmgam,c11,cnn,zcmgam(1,ifirsttime),ierr)
      cnn(1) = nsilop
      do i=1, nsilop
      do j=1, ntime
         zcsilop(i,j) = real(csilop(i,j))
      enddo
      enddo
      call NCVPT (nceq,id_csilop,c11,cnn,zcsilop(1,ifirsttime),ierr)
      cnn(1) = magpri
      do i=1, magpri
      do j=1, ntime
         zcmpr2(i,j) = real(cmpr2(i,j))
      enddo
      enddo
      call NCVPT (nceq,id_cmpr2,c11,cnn,zcmpr2(1,ifirsttime),ierr)
      cnn(1) = nfcoil
      do i=1, nfcoil
      do j=1, ntime
         zccbrsp(i,j) = real(ccbrsp(i,j))
      enddo
      enddo
      call NCVPT (nceq,id_ccbrsp,c11,cnn,zccbrsp(1,ifirsttime),ierr)
      cnn(1) = nesum
      do i=1, ntime
         zwork(i) = real(diamag(i))
      enddo
      call NCVPT (nceq,id_diamag,m,n,zwork(ifirsttime),ierr)
      do i=1, ntime
         zwork(i) = real(sigdia(i))
      enddo
      call NCVPT (nceq,id_sigdia,m,n,zwork(ifirsttime),ierr)
      do i=1, ntime
         zwork(i) = real(cdflux(i))
      enddo
      call NCVPT (nceq,id_cdflux,m,n,zwork(ifirsttime),ierr)
!
      do i=1, ntime
         zwork(i) = real(pasmat(i))
      enddo
      call NCVPT (nceq,id_plasma,m,n,zwork(ifirsttime),ierr)
      do i=1, ntime
         zwork(i) = real(cpasma(i))
      enddo
      call NCVPT (nceq,id_cpasma,m,n,zwork(ifirsttime),ierr)
!
! --- following data and their corresponding ncdf variables have 
! --- reversed dimensions
!
      call NCVPT1 (nceq,id_fungam,1,msefitfun,ierr)

      stride(1) = 1
      stride(2) = 1
      cnn(1) = nstark
      imap(2) = 4                         ! number of bytes in float
!     imap(2) = imap(2)*2                 ! number of bytes in double
      imap(1) = imap(2)*ntime
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(tangam(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_tangam,c11,cnn,stride,imap, &
           zstark(ifirsttime,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(tangam_uncor(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_tangam_uncor,c11,cnn,stride,imap, &
           zstark(ifirsttime,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(spatial_fix(j,i))
      enddo
      enddo
      call NCVPTG (nceq,id_fgam,c11,cnn,stride,imap, &
           zstark(ifirsttime,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(siggam(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_siggam,c11,cnn,stride,imap, &
           zstark(ifirsttime,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(rrgam(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_rrgam,c11,cnn,stride,imap, &
           zstark(ifirsttime,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(zzgam(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_zzgam,c11,cnn,stride,imap, &
           zstark(ifirsttime,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(a1gam(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_a1gam,c11,cnn,stride,imap, &
           zstark(ifirsttime,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(a2gam(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_a2gam,c11,cnn,stride,imap, &
           zstark(ifirsttime,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(a3gam(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_a3gam,c11,cnn,stride,imap, &
           zstark(ifirsttime,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(a4gam(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_a4gam,c11,cnn,stride,imap, &
           zstark(ifirsttime,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(a5gam(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_a5gam,c11,cnn,stride,imap, &
           zstark(ifirsttime,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(a6gam(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_a6gam,c11,cnn,stride,imap, &
           zstark(ifirsttime,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(a7gam(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_a7gam,c11,cnn,stride,imap, &
           zstark(ifirsttime,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(a8gam(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_a8gam,c11,cnn,stride,imap, &
           zstark(ifirsttime,1),ierr)

      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(rmse_gain(j))
      enddo
      enddo
      call NCVPTG (nceq,id_gaingam,c11,cnn,stride,imap, &
           zstark(ifirsttime,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(rmse_slope(j))
      enddo
      enddo
      call NCVPTG (nceq,id_slopegam,c11,cnn,stride,imap, &
           zstark(ifirsttime,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(rmse_scale(j))
      enddo
      enddo
      call NCVPTG (nceq,id_scalegam,c11,cnn,stride,imap, &
           zstark(ifirsttime,1),ierr)
      do i=1, ntime
      do j=1, nstark
         zstark(i,j) = real(rmse_offset(j))
      enddo
      enddo
      call NCVPTG (nceq,id_offsetgam,c11,cnn,stride,imap, &
           zstark(ifirsttime,1),ierr)
!
      cnn(1) = nsilop
      do i=1, ntime
      do j=1, nsilop
         zsilopt(i,j) = real(silopt(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_silopt,c11,cnn,stride,imap, &
           zsilopt(ifirsttime,1),ierr)
      cnn(1) = magpri
      do i=1, ntime
      do j=1, magpri
         zexpmpi(i,j) = real(expmpi(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_expmpi,c11,cnn,stride,imap, &
           zexpmpi(ifirsttime,1),ierr)
      cnn(1) = nfcoil
      do i=1, ntime
      do j=1, nfcoil
         zfccurt(i,j) = real(fccurt(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_fccurt,c11,cnn,stride,imap, &
           zfccurt(ifirsttime,1),ierr)
      cnn(1) = nesum
      do i=1, ntime
      do j=1, nesum
         zeccurt(i,j) = real(eccurt(i,j))
      enddo
      enddo
      call NCVPTG (nceq,id_eccurt,c11,cnn,stride,imap, &
           zeccurt(ifirsttime,1),ierr)
!
      call NCCLOS (nceq,ierr)             ! close the file
      endif
!
      return
      end subroutine wmeasure

