!**********************************************************************
!**
!**     efund_getsizes performs inputing and initialization from the
!**       mhdin.dat file
!**
!**       most default values set here are for DIII-D
!**
!**********************************************************************
      subroutine efund_getsizes
      
      use machine
      use grid
      use siloop
      use ecoil
      use fcoil
      use acoil
      use nio
      use mprobe
      use vessel

      implicit none
      integer*4:: istat,icycred_loopmax,kubics,npcurn,nwcur2,nwcurn
      character(len=1000) :: line

      namelist/machinein/nfcoil,nsilop,magpri,nrogow,necoil,nesum, &
                         nfsum,nvsum,nvesel,nacoil,device, &
                         nnece,nnecein,neceo,mpress,nmsels,nmselp,libim, &
                         nnnte,ngam_vars,ngam_u,ngam_w,nlimit,nlimbd,nangle, &
                         ntangle,nfbcoil,mccoil,micoil,ndata,nwwcur, &
                         nffcur,nppcur,nercur,ntime,ndim,kxiter,mqwant, & 
                         mbdry,mbdry1,nxtram,nxtlim,nco2v,nco2r,modef, &
                         modep,modew,kubics,icycred_loopmax,nfourier

      device = 'DIII-D'
      magpri = 76
      nacoil = 1
      necoil = 122
      nesum = 6
      nfcoil = 18
      nfsum = 18
      nrogow = 1
      nsilop = 44
      nvesel = 24
      nvsum = 24

      nnece=40
      nnecein=80
      neceo=1
      nmsels=16
      nmselp=69
      libim=32
      mpress=201
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
      nfourier=5
      
      open(unit=nin,status='old',file='mhdin.dat',iostat=istat)

      read (nin,machinein,iostat=istat)
      if (istat>0) then
        backspace(nin)
        read(nin,fmt='(A)') line
        write(*,'(A)') 'Invalid line in namelist machinein: '//trim(line)
        return
      endif

      close(nin)

      allocate(rsi(nsilop),zsi(nsilop),wsi(nsilop),hsi(nsilop),&
               as(nsilop),as2(nsilop))
      rsi = 0.0
      zsi = 0.0
      wsi = 0.0
      hsi = 0.0
      as = 0.0
      as2 = 0.0

      allocate(re(necoil),ze(necoil),he(necoil),we(necoil), &
               ecid(necoil),ecturn(necoil))
      re = 0.0
      ze = 0.0
      he = 0.0
      we = 0.0
      ecid = 0
      ecturn = 0.0

      allocate(rf(nfcoil),zf(nfcoil),wf(nfcoil),hf(nfcoil), &
               af(nfcoil),af2(nfcoil),turnfc(nfsum),fcid(nfcoil), & 
               fcturn(nfcoil))
      rf = 0.0
      zf = 0.0
      wf = 0.0
      hf = 0.0
      af = 0.0
      af2 = 0.0
      turnfc = 0.0
      fcturn = 0.0
      fcid = 0

      allocate(racoil(nacoil),zacoil(nacoil),wacoil(nacoil),hacoil(nacoil))
      racoil = 0.0
      zacoil = 0.0
      wacoil = 0.0
      hacoil = 0.0

      allocate(xmp2(magpri),ymp2(magpri),amp2(magpri),smp2(magpri))
      xmp2 = 0.0
      ymp2 = 0.0
      amp2 = 0.0
      smp2 = 0.0

      allocate(rvs(nvesel),zvs(nvesel),wvs(nvesel),hvs(nvesel),&
               avs(nvesel),avs2(nvesel),rsisvs(nvesel),vsid(nvesel))
      rvs = 0.0
      zvs = 0.0
      wvs = 0.0
      hvs = 0.0
      avs = 0.0
      avs2 = 0.0
      rsisvs = 0.0

      allocate(nshiftrz(nfcoil))
      nshiftrz = 0

      allocate(rshift(nfcoil),zshift(nfcoil),pshift(nfcoil))
      rshift = 0.0
      zshift = 0.0
      pshift = 0.0

      allocate(pmprobe(magpri))
      pmprobe = 0.

      open(unit=nout,status='unknown',file='mhdout.dat',delim='quote')
      write (nout,machinein)
      close(nout)
      end subroutine efund_getsizes

!**********************************************************************
!**                                                                  **
!**     getset performs inputing and initialization.                 **
!**                                                                  **
!**********************************************************************
      subroutine efund_getset

      use grid
      use siloop
      use ecoil
      use fcoil
      use acoil
      use nio
      use mprobe
      use vessel
      use utils, only: mgaus1,mgaus2
      implicit none
      integer*4 i,istat
      real*8 patmp2(magpri)
      character*10 mpnam2(magpri),lpname(nsilop),vsname(nvesel)
      character(1000) line

      namelist/in3/mpnam2,xmp2,ymp2,amp2,smp2,patmp2, &
                   rsi,zsi,wsi,hsi,as,as2,rsisvs,lpname, &
                   rvs,zvs,wvs,hvs,avs,avs2,vsid,vsname, &
                   racoil,zacoil,wacoil,hacoil, &
                   rf,zf,wf,hf,af,af2,fcid,fcturn,turnfc, &
                   re,ze,we,he,ecid,ecturn
      namelist/in5/rleft,rright,zbotto,ztop,mgaus1,mgaus2, &
                   nshiftrz,rshift,zshift,pshift,pmprobe,nsmp2, &
                   igrid,ifcoil,islpfc,iecoil,ivesel,iacoil,isize
!
      open(unit=nin,status='old',file='mhdin.dat')
!
      mpnam2=''
      lpname=''
      vsname=''
!
      ifcoil=0
      igrid=0
      islpfc=0
      iecoil=0
      ivesel=0
      isize=0
      mgaus1=8
      mgaus2=10
      nsmp2=1
      rleft=0.
      rright=0.
      ztop=0.
      zbotto=0.
      rsi(1)=-1
      rf(1)=-1.
      re(1)=-1.
      rvs(1)=-1.
      wvs(1)=-1.
      patmp2=0.
!---------------------------------------------------------------------
!--   isize=0      no finite size correction for flux loops         --
!--         1         finite size correction for flux loops         --
!--   islpfc=1     flux loops at F coils                            --
!---------------------------------------------------------------------
      read (nin,in3,iostat=istat)
      if (istat>0) then
        backspace(nin)
        read(nin,fmt='(A)') line
        write(*,'(A)') 'Invalid line in namelist in3: '//trim(line)
        return
      endif
      rewind(nin)
      read (nin,in5,iostat=istat)
      if (istat>0) then
        backspace(nin)
        read(nin,fmt='(A)') line
        write(*,'(A)') 'Invalid line in namelist in5: '//trim(line)
        return
      endif
!
      if (.not. allocated(rgrid)) then
        allocate(rgrid(nw))
        rgrid(:) = 0.0
      endif
      if (.not. allocated(zgrid)) then
        allocate(zgrid(nh))
        zgrid(:) = 0.0
      endif
!----------------------------------------------------------------------
!--   read f coil and psi loop dimensions                            --
!----------------------------------------------------------------------
      if (rf(1).lt.0.0) then
        read (nin,10000) (rf(i),zf(i),wf(i),hf(i),af(i),af2(i), &
                          i=1,nfcoil)
      endif
      if (rsi(1).lt.0.0) then
        read (nin,10000) (rsi(i),zsi(i),wsi(i),hsi(i),as(i),as2(i), &
                          i=1,nsilop)
      endif
      if ((iecoil.eq.1).or.(ivesel.eq.1)) then
        if (re(1).lt.0.0) then
          read (nin,10020) (re(i),ze(i),we(i),he(i),ecid(i),i=1,necoil)
        endif
        if (ivesel.eq.1.and.rvs(1).lt.0.0) then
          if (wvs(1).lt.0.0) then
            read (nin,10000) (rvs(i),zvs(i),wvs(i),hvs(i),avs(i),avs2(i), &
                              i=1,nvesel)
          else
            do i=1,nvesel
              read (nin,*) rvs(i),zvs(i)
            enddo
          endif
        endif
      endif
!----------------------------------------------------------------------
!--   compute r and z arrays                                         --
!----------------------------------------------------------------------
      dr=(rright-rleft)/dble(nw-1)
      dz=(ztop-zbotto)/dble(nh-1)
      do i=1,nw
        rgrid(i)=rleft+dr*(i-1)
      enddo 
      do i=1,nh
        zgrid(i)=zbotto+dz*(i-1)
      enddo 

      close(nin)
      open(unit=nout,status='unknown',file='mhdout.dat', &
           position='append',delim='quote')
      write (nout,in3)
      write (nout,in5)
      close(nout)
      return
10000 format (6e12.6)
10010 format (4e12.6)
10020 format (5e10.4)
      end subroutine efund_getset
