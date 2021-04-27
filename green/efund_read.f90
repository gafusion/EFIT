      SUBROUTINE efund_getsizes
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          getset performs inputing and initialization.            **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          02/02/20..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      
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

      NAMELIST/machinein/nfcoil,nsilop,magpr2,nrogow,necoil,nesum,&
                      nfsum,nvsum,nvesel,nacoil,mgaus1,mgaus2


      nfcoil = 18
      nsilop = 44
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

      OPEN(unit=nin,status='old',file='mhdin.dat' &
          )

      READ (nin,machinein,err=10)

   10 CLOSE(nin)

      allocate(rsi(nsilop),zsi(nsilop),wsi(nsilop),hsi(nsilop),&
               as(nsilop),as2(nsilop))
      rsi(:) = 0.0
      zsi(:) = 0.0
      wsi(:) = 0.0
      hsi(:) = 0.0
      as(:) = 0.0
      as2(:) = 0.0

      allocate(re(necoil),ze(necoil),he(necoil),we(necoil), &
               ecid(necoil),ecturn(necoil))
      re(:) = 0.0
      ze(:) = 0.0
      he(:) = 0.0
      we(:) = 0.0
      ecid(:) = 0.0
      ecturn(:) = 0.0

      allocate(rf(nfcoil),zf(nfcoil),wf(nfcoil),hf(nfcoil), &
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

      allocate(racoil(nacoil),zacoil(nacoil),wacoil(nacoil),hacoil(nacoil))
      racoil(:) = 0.0
      zacoil(:) = 0.0
      wacoil(:) = 0.0
      hacoil(:) = 0.0

      allocate(xmp2(magpr2),ymp2(magpr2),amp2(magpr2),smp2(magpr2))
      xmp2(:) = 0.0
      ymp2(:) = 0.0
      amp2(:) = 0.0
      smp2(:) = 0.0
      allocate(rvs(nvesel),zvs(nvesel),wvs(nvesel),hvs(nvesel),&
               avs(nvesel),avs2(nvesel),rsisvs(nvesel),vsid(nvesel))
      rvs(:) = 0.0
      zvs(:) = 0.0
      wvs(:) = 0.0
      hvs(:) = 0.0
      avs(:) = 0.0
      avs2(:) = 0.0
      rsisvs(:) = 0.0

      vsid(:) = 0.0

      allocate(nshiftrz(nfcoil))
      nshiftrz(:) = 0.

      allocate(rshift(nfcoil),zshift(nfcoil),pshift(nfcoil))
      rshift(:) = 0.0
      zshift(:) = 0.0
      pshift(:) = 0.0

      allocate(pmprobe(magpr2))
      pmprobe(:) = 0.
      END SUBROUTINE efund_getsizes

      SUBROUTINE efund_getset
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          getset performs inputing and initialization.            **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      USE exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,nesum,&
                      nfsum,nvsum,nvesel,nacoil,nw,nh,nwnh
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
!vas
      use var_filech
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
      REAL*8, dimension (:), allocatable :: patmp2
      CHARACTER*10,dimension (:), allocatable:: mpnam2,lpname,vsname

      NAMELIST/in3/igrid,rleft,rright,zbotto,ztop,ifcoil &
           ,islpfc,iecoil,mpnam2,xmp2,ymp2,amp2,smp2,isize,rsi,zsi,wsi &
           ,hsi,as,as2,lpname,nsmp2,ivesel,rsisvs,vsname,turnfc,patmp2 &
           ,iacoil,racoil,zacoil,wacoil,hacoil &
           ,rf,zf,fcid,wf,hf,wvs,hvs,avs,avs2,af,af2,fcturn &
           ,re,ze,ecid,ecturn,vsid,rvs,zvs,we,he &
           ,nshiftrz,rshift,zshift,pshift,pmprobe &
           ,nw,nh
!
      OPEN(unit=nin,status='old',file='mhdin.dat' &
          )
      OPEN(unit=nout,status='unknown',file='mhdout.dat' &
          )
!
      allocate(patmp2(magpr2),mpnam2(magpr2),lpname(nsilop),vsname(nvesel))
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
      DO i=1,nfcoil
        nshiftrz(i)=0
      ENDDO
!---------------------------------------------------------------------
!--  isize=0      no finite size correction for flux loops          --
!--        1         finite size correction for flux loops          --
!--  islpfc=1     flux loops at F coils                             --
!---------------------------------------------------------------------
      READ (nin,in3)
      !WRITE (nout,in3)
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
!vas
      call inp_file_ch(nw,nh,ch1,ch2)
!      print*,'file name : ', 'ep'//trim(ch1)// &
!                         trim(ch2)//'.ddd'
!
     print *, 'rf(1)=',rf
!----------------------------------------------------------------------
!-- READ f coil and psi loop dimensions                              --
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
  200 CONTINUE
!----------------------------------------------------------------------
!--  compute r and z arrays                                          --
!----------------------------------------------------------------------
      dr=(rright-rleft)/float(nw-1)
      dz=(ztop-zbotto)/float(nh-1)
      DO i=1,nw
         rgrid(i)=rleft+dr*(i-1)
      ENDDO 
      DO i=1,nh
         zgrid(i)=zbotto+dz*(i-1)
      ENDDO 
  300 CONTINUE
!
!vas
      WRITE (nout,in3)
      close(nin)
      close(nout)
      RETURN
10000 FORMAT (6e12.6)
10010 FORMAT (4e12.6)
10020 FORMAT (5e10.4)
      END SUBROUTINE efund_getset
