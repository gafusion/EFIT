!**********************************************************************
!**                                                                  **
!**     getset performs inputing and initialization.                 **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          02/02/20..........first created                         **
!**                                                                  **
!**********************************************************************
      SUBROUTINE efund_getsizes
      
      USE exparm
      USE siloop
      USE cecoil
      USE fcoil
      USE pmodel
      USE consta
      USE input
      USE cacoil
      USE nio
      USE mprobe
      USE cvesel
      USE fshift

      integer*4:: istat,icycred_loopmax
      character(len=1000) :: line

      NAMELIST/machinein/nfcoil,nsilop,magpr2,nrogow,necoil,nesum, &
                         nfsum,nvsum,nvesel,nacoil,mgaus1,mgaus2,device, &
                         magpri322,magpri67,magprirdp,magudom,maglds, &
                         nnece,nnecein,neceo,mse315,mse45,mse15, &
                         mse1h,mse315_2,mse210,mpress,libim,nmsels, &
                         nnnte,ngam_vars,ngam_u,ngam_w,nlimit,nlimbd,nangle, &
                         ntangle,nfbcoil,mccoil,micoil,ndata,nwwcur, &
                         nffcur,nppcur,nercur,ntime,ndim,kxiter,mqwant, & 
                         mbdry,mbdry1,nxtram,nxtlim,nco2v,nco2r,modef, &
                         modep,modew,kubics,icycred_loopmax,nfourier, &
                         nsilol,nsilds

      device = 'DIII-D'
      nfcoil = 18
      nsilop = 44
      nsilol = -43
      nsilds = -1 ! default all psi loops to one side except 1
      magpr2 = 76
      nrogow = 1
      necoil = 122
      nesum = 6
      nfsum = 18
      nvsum = 24
      nvesel = 24
      nacoil = 1
      mgaus1 = 8
      mgaus2 = 10

      nnece=40
      nnecein=80
      neceo=1
      mse315=40
      mse45=0
      mse15=0
      mse1h=0
      mse315_2=0
      mse210=0
      mpress=201
      libim=32
      nmsels=16
      nnnte=801
      ngam_vars=9
      ngam_u=5
      ngam_w=3
      nlimit=160
      nlimbd=6
      nangle=64
      ntangle=12
      nfbcoil=12
      mccoil=6
      micoil=12     
      ndata=61
      nwwcur=32
      nffcur=32
      nppcur=32
      nercur=32
      npcurn=nffcur+nppcur
      nwcurn=nwwcur+npcurn
      nwcur2=nwcurn*2
      ntime=1001
      ndim=3200
      kxiter=515
!      mqwant=30 ! DIIID default
      mqwant=66 ! number used in NSTX
      mbdry=300
      mbdry1=110
      nxtram=10
      nxtlim=9
      nco2v=3
      nco2r=2
      modef=4
      modep=4
      modew=4
      kubics=4
      nfourier=5
      
      magpri67=-1
      magprirdp=-1
      magudom=-1
      maglds=-1
      
      OPEN(unit=nin,status='old',file='mhdin.dat',iostat=istat)

      READ (nin,machinein)

      if (istat>0) then
        backspace(nin)
        read(nin,fmt='(A)') line
        write(*,'(A)') 'Invalid line in namelist machinein: '//trim(line)
        stop
      endif

      CLOSE(nin)

      ! Handle if user probe array sizes aren't specified
      IF (trim(device)=='DIII-D') THEN
        mse315=11
        mse45=15
        mse15=10
        mse1h=4
        mse315_2=5
        mse210=24
      ENDIF

      IF (magpri67<0 .or. magprirdp<0 .or. magudom<0 .or. magudom<0) then 
        IF (trim(device)=='DIII-D') THEN
          magpri67=29
          magpri322=31
          magprirdp=8
          magudom=5
          maglds=3
        ELSE
          magpri67 = abs(magpri67)
          magprirdp = abs(magprirdp)
          magudom = abs(magudom)
          maglds = abs(maglds)
          magpri322 = magpr2 - magpri67 - magprirdp - magudom - maglds
        ENDIF
      ENDIF

      IF (nsilol<0 .or. nsilds<0) then 
        IF (trim(device)=='DIII-D') THEN
          nsilds = 3
          nsilol = 41
        ELSE
          nsilol = nsilop -1
          nsilds = 1
        ENDIF
      ENDIF

      ALLOCATE(rsi(nsilop),zsi(nsilop),wsi(nsilop),hsi(nsilop),&
               as(nsilop),as2(nsilop))
      rsi(:) = 0.0
      zsi(:) = 0.0
      wsi(:) = 0.0
      hsi(:) = 0.0
      as(:) = 0.0
      as2(:) = 0.0

      ALLOCATE(re(necoil),ze(necoil),he(necoil),we(necoil), &
               ecid(necoil),ecturn(necoil))
      re(:) = 0.0
      ze(:) = 0.0
      he(:) = 0.0
      we(:) = 0.0
      ecid(:) = 0.0
      ecturn(:) = 0.0

      ALLOCATE(rf(nfcoil),zf(nfcoil),wf(nfcoil),hf(nfcoil), &
               af(nfcoil),af2(nfcoil),turnfc(nfcoil), fcid(nfcoil), & 
               fcturn(nfcoil))
      rf(:) = 0.0
      zf(:) = 0.0
      wf(:) = 0.0
      hf(:) = 0.0
      af(:) = 0.0
      af2(:) = 0.0
      turnfc(:) = 0.0
      fcturn(:) = 0.0

      ALLOCATE(racoil(nacoil),zacoil(nacoil),wacoil(nacoil),hacoil(nacoil))
      racoil(:) = 0.0
      zacoil(:) = 0.0
      wacoil(:) = 0.0
      hacoil(:) = 0.0

      ALLOCATE(xmp2(magpr2),ymp2(magpr2),amp2(magpr2),smp2(magpr2))
      xmp2(:) = 0.0
      ymp2(:) = 0.0
      amp2(:) = 0.0
      smp2(:) = 0.0
      ALLOCATE(rvs(nvesel),zvs(nvesel),wvs(nvesel),hvs(nvesel),&
               avs(nvesel),avs2(nvesel),rsisvs(nvesel),vsid(nvesel))
      rvs(:) = 0.0
      zvs(:) = 0.0
      wvs(:) = 0.0
      hvs(:) = 0.0
      avs(:) = 0.0
      avs2(:) = 0.0
      rsisvs(:) = 0.0

      vsid(:) = 0.0

      ALLOCATE(nshiftrz(nfcoil))
      nshiftrz(:) = 0.

      ALLOCATE(rshift(nfcoil),zshift(nfcoil),pshift(nfcoil))
      rshift(:) = 0.0
      zshift(:) = 0.0
      pshift(:) = 0.0

      ALLOCATE(pmprobe(magpr2))
      pmprobe(:) = 0.
      END SUBROUTINE efund_getsizes

!**********************************************************************
!**                                                                  **
!**     getset performs inputing and initialization.                 **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**********************************************************************
      SUBROUTINE efund_getset

      USE exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,nesum,&
                      nfsum,nvsum,nvesel,nacoil,nw,nh,nwnh,device
      USE siloop
      USE cecoil
      USE fcoil
      USE pmodel
      USE consta
      USE input
      USE cacoil
      USE nio
      USE mprobe
      USE cvesel
      USE fshift
      use var_filech
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
      REAL*8, dimension (:), allocatable :: patmp2
      CHARACTER*10,dimension (:), allocatable:: mpnam2,lpname,vsname

      NAMELIST/in3/igrid,rleft,rright,zbotto,ztop,ifcoil, &
                   islpfc,iecoil,mpnam2,xmp2,ymp2,amp2,smp2,isize, &
                   rsi,zsi,wsi,hsi,as,as2,lpname,nsmp2,ivesel,rsisvs, &
                   vsname,turnfc,patmp2, &
                   iacoil,racoil,zacoil,wacoil,hacoil, &
                   rf,zf,fcid,wf,hf,wvs,hvs,avs,avs2,af,af2,fcturn, &
                   re,ze,ecid,ecturn,vsid,rvs,zvs,we,he, &
                   nshiftrz,rshift,zshift,pshift,pmprobe, &
                   nw,nh
!
      OPEN(unit=nin,status='old',file='mhdin.dat')
      OPEN(unit=nout,status='unknown',file='mhdout.dat')
!
      ALLOCATE(patmp2(magpr2),mpnam2(magpr2),lpname(nsilop),vsname(nvesel))
      vsname=''
      nw = 65
      nh = 65
!
      ifcoil=0
      igrid=0
      isize=0
      islpfc=0
      iecoil=0
      ivesel=0
      nsmp2=1
      rsi(1)=-1
      rf(1)=-1.
      re(1)=-1.
      rvs(1)=-1.
      wvs(1)=-1.
      nshiftrz(1:nfcoil)=0
      patmp2=0.
!---------------------------------------------------------------------
!--   isize=0      no finite size correction for flux loops         --
!--         1         finite size correction for flux loops         --
!--   islpfc=1     flux loops at F coils                            --
!---------------------------------------------------------------------
      READ (nin,in3)
!
      IF (.NOT. ALLOCATED(rgrid)) THEN
        ALLOCATE(rgrid(nw))
        rgrid(:) = 0.0
      ENDIF
      IF (.NOT. ALLOCATED(zgrid)) THEN
        ALLOCATE(zgrid(nh))
        zgrid(:) = 0.0
      ENDIF
!
      nwnh = nw * nh
!make the file names for green-table
      CALL inp_file_ch(nw,nh,ch1,ch2)
!----------------------------------------------------------------------
!--   READ f coil and psi loop dimensions                            --
!----------------------------------------------------------------------
      IF (rf(1).lt.0.0) THEN
        READ (nin,10000) (rf(i),zf(i),wf(i),hf(i),af(i),af2(i), &
                          i=1,nfcoil)
      ENDIF
      IF (rsi(1).lt.0.0) THEN
        READ (nin,10000) (rsi(i),zsi(i),wsi(i),hsi(i),as(i),as2(i), &
                          i=1,nsilop)
      ENDIF
      IF ((iecoil.gt.0).or.(ivesel.gt.0)) THEN
        IF (re(1).lt.0.0) THEN
          READ (nin,10020) (re(i),ze(i),we(i),he(i),ecid(i),i=1,necoil)
        ENDIF
        IF (ivesel.gt.0.and.rvs(1).lt.0.0) THEN
          IF (wvs(1).lt.0.0) THEN
            READ (nin,10000) (rvs(i),zvs(i),wvs(i),hvs(i),avs(i),avs2(i), &
                              i=1,nvesel)
          ELSE
            DO i=1,nvesel
              READ (nin,*) rvs(i),zvs(i)
            ENDDO
          ENDIF
        ENDIF
      ENDIF
!----------------------------------------------------------------------
!--   compute r and z arrays                                         --
!----------------------------------------------------------------------
      dr=(rright-rleft)/dble(nw-1)
      dz=(ztop-zbotto)/dble(nh-1)
      DO i=1,nw
        rgrid(i)=rleft+dr*(i-1)
      ENDDO 
      DO i=1,nh
        zgrid(i)=zbotto+dz*(i-1)
      ENDDO 

      WRITE (nout,in3)
      CLOSE(nin)
      CLOSE(nout)

      CALL dprobe_machinein(nfcoil,nsilop,magpr2,nrogow,necoil,nesum,&
                            nfsum,nvsum,nvesel,nacoil)
      
      CALL dprobe(mpnam2,lpname,patmp2)
      
      RETURN
10000 FORMAT (6e12.6)
10010 FORMAT (4e12.6)
10020 FORMAT (5e10.4)
      END SUBROUTINE efund_getset
