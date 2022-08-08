#include "config.f"
!**********************************************************************
!>
!!    get_constraints sets the pointnames and calls subroutines to
!!      query PTDATA or MDS+ to obtain the measurement constraints.
!!
!!
!!    @param nshot : shot number
!!
!!    @param times : first time requested (in seconds)
!!
!!    @param delt : length between time slices (in seconds)
!!
!!    @param np : number of time slices
!!
!!    @param iierr : error flag
!!
!**********************************************************************
      subroutine get_constraints(nshot,times,delt,np,iierr)
      use vtime_mod
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      character*10 nsingl(10),n1name,btcname, &
                   nc79name,nc139name,ncname(mccoil),niname(micoil)  !EJS(2014)

      integer*4 time_err,ioerr
      character*150 textline     !EJS(2014)
      character*10,dimension(:),allocatable :: ndenv,ndenr,fcname,ecname
      character*10 namedum
      character(len=1000) :: line
      real*8 dumbtc
      real*8,dimension(:),allocatable :: dumccc,dumcic
      logical read_btcshot

      NAMELIST/in3/mpnam2,xmp2,ymp2,amp2,smp2,rsi,zsi,wsi,hsi, &
                   as,as2,lpname,rsisvs,turnfc,patmp2,vsname, &
                   racoil,zacoil,wacoil,hacoil, &
                   rf,zf,fcid,wf,hf,wvs,hvs,avs,avs2,af,af2, &
                   re,ze,ecid,ecturn,vsid,rvs,zvs,we,he,fcturn
!
      namelist/in4/mpnam2,lpname,vsname,nsingl,n1name, & !JRF 
                   nc79name,nc139name,btcname,ndenv,ndenr, &
                   fcname,ecname
!
      ALLOCATE(ndenv(nco2v),ndenr(nco2r),fcname(nfcoil),ecname(nesum))

      nsingl(1) = 'IP        '
      nsingl(2) ='VLOOP     '
      nsingl(3) ='BCOIL     '
      nsingl(4) ='DIAMAG3   '   !diamagnetic flux 
      nsingl(5) ='RL01      '   !full rogowski 
      nsingl(6)='PINJ      '    !beam injection power
      n1name='N1COIL    '
      nc79name='C79       '
      nc139name='C139      '
      btcname ='BTI322    '
      ndenv(1)='DENV1     '
      ndenv(2)='DENV2     '
      ndenv(3)='DENV3     '
      ndenr(1)='DENR0     '
      ndenr(2)='DENMW     '
      fcname(1)='F1A       '
      fcname(2)='F2A       '
      fcname(3)='F3A       '
      fcname(4)='F4A       '
      fcname(5)='F5A       '
      fcname(6)='F6A       '
      fcname(7)='F7A       '
      fcname(8)='F8A       '
      fcname(9)='F9A       '
      fcname(10)='F1B       '
      fcname(11)='F2B       '
      fcname(12)='F3B       '
      fcname(13)='F4B       '
      fcname(14)='F5B       '
      fcname(15)='F6B       '
      fcname(16)='F7B       '
      fcname(17)='F8B       '
      fcname(18)='F9B       '
      ecname(1)='ECOILA    '
      ecname(2)='ECOILB    '
      ecname(3)='E567UP    '
      ecname(4)='E567DN    ' 
      ecname(5)='E89DN     '
      ecname(6)='E89UP     '
      data irdata/0/,baddat/0/
!----------------------------------------------------------------------
!--   read in pointnames
!----------------------------------------------------------------------
      open(unit=60,file=table_di2(1:ltbdi2)//'dprobe.dat', &
           status='old')
      read(60,in3,iostat=istat)
      if (istat>0) then
        backspace(nin)
        read(nin,fmt='(A)') line
        write(*,'(A)') 'Invalid line in namelist in3: '//trim(line)
        stop
      endif
      close(unit=60)
!
! !JRF The if statement here tests a flag from the snap file that
!     indicates whether the pointnames should be read from an alternate
!     file.
!
      if (use_alternate_pointnames .ne. 0) then
!
        open(unit=60,file=alternate_pointname_file,status='old')
        read(60,in4,iostat=istat)
        if (istat>0) then
          backspace(nin)
          read(nin,fmt='(A)') line
          write(*,'(A)') 'Invalid line in namelist in4: '//trim(line)
          stop
        endif
        close(unit=60)
      endif
!
!     times=times+0.000001
!     delt=delt+0.000001
      deltm=0.000001
      delfl=2.
      delno=.04
      delta=.001
      nptt=20
      do i=1,np
        time(i)=times+delt*(i-1)
      enddo
      krl01=0
      if(iierr.lt.0) krl01=1
      iierr=0
      i1 = 1
      i0 = 0
      r1 = 1.
!----------------------------------------------------------------------
!--   psi-loops
!----------------------------------------------------------------------
      do i=1,nsilop
        silopt(1:np,i)=0.
        ierpsi(i)=0
        call avdata(nshot,lpname(i),i1,ierpsi(i),silopt(1:np,i), &
                    np,times,delt,i0,r1,i1,psibit(i),iavem,time(1:np), &
                    ircfact,psi_rc(i),psircg(i), &
                    vrespsi(i),psi_k(i), &
                    t0psi(i),devpsi(1:np,i),navpsi(1:np,i),time_err)
        if (ierpsi(i).eq.3) then ! TODO: this doesn't appear to be possible
          iierr=1
          return
        endif
        if (i.eq.iabs(nslref)) then
          psiref(1:np)=silopt(1:np,i)
          silopt(1:np,i)=0.0
        endif
      enddo
      ierpsi(iabs(nslref))=0
      rnavpsi=navpsi
!----------------------------------------------------------------------
!--  plasma current                                                  --
!----------------------------------------------------------------------
      if (krl01.eq.0) then
        i=1
      else
        i=5
      endif
      if(use_alternate_pointnames .eq. 2) i = 5          !***JRF
      pasmat(1:np)=0.
      ierpla=0
      call avdata(nshot,nsingl(i),i1,ierpla,pasmat(1:np), &
                  np,times,delt,i0,r1,i1,bitip,iavem,time(1:np),ircfact, &
                  p_rc,prcg,vresp,p_k,t0p,devp(1:np), &
                  navp(1:np),time_err)

      rnavp=REAL(navp)
      if( (use_alternate_pointnames .eq. 1) .and. &      !JRF 
          (i .eq. 1) ) pasmat(1:np)=pasmat(1:np)*0.5e6
!----------------------------------------------------------------------
!--  loop voltage                                                    --
!----------------------------------------------------------------------
      vloopt(1:np)=0.
      ierlop=0
      call avdata(nshot,nsingl(2),i1,ierlop,vloopt(1:np), &
                  np,times,delt,i0,r1,i1,bitvl,iavev,time(1:np),ircfact, &
                  vl_rc,vlrcg,vresvl,vl_k,t0vl,devvl(1:np), &
                  navvl(1:np),time_err)
!
      if (use_alternate_pointnames .eq. 2) then    !JRF
        if( (ierpla .eq. 0) .and. (ierlop .eq. 0)) &
          pasmat(1:np) = pasmat(1:np) - vloopt(1:np) * 0.646/50.0 * 0.5e6
      endif
!---------------------------------------------------------------------
!--   Get density array from PTDATA or MDS+                         --
!---------------------------------------------------------------------
      denvt = 0.0
      denrt = 0.0
      s124411: if (nshot.lt.124411) then
      do i=1,3
        ierlop=0
        call avdata(nshot,ndenv(i),i1,ierlop,denvt(1:np,i), &
                    np,times,delt,i0,r1,i1,bitvl,iaved,time(1:np),ircfact, &
                    denv_rc(i),denvrcg(i),vresdenv(i), &
                    denv_k(i),t0denv(i),devdenv(1:np,i),navdenv(1:np,i), &
                    time_err)
        if(ierlop.eq.0) denvt(1:np,i)=denvt(1:np,i)*50.0
      enddo
      do i=1,2
        ierlop=0
        call avdata(nshot,ndenr(i),i1,ierlop,denrt(1:np,i), &
                    np,times,delt,i0,r1,i1,bitvl,iaved,time(1:np),ircfact, &
                    denr_rc(i),denrrcg(i),vresdenr(i), &
                    denr_k(i),t0denr(i),devdenr(1:np,i),navdenr(1:np,i), &
                    time_err)
        if(ierlop.eq.0) denrt(1:np,i)=denrt(1:np,i)*50.0
      enddo
      else s124411
!---------------------------------------------------------------------
!--   Get density array from MDS+                                   --
!---------------------------------------------------------------------
#if defined(USE_MDS)
      do i=1,3
        ierlop=0
        call amdata(nshot,ndenv(i),i1,ierlop,denvt(1:np,i), &
                    np,times,delt,i0,r1,i1,bitvl,iaved,time(1:np),ircfact)
        if(ierlop.eq.0) denvt(1:np,i)=denvt(1:np,i)*50.0
      enddo
      ierlop=0
      call amdata(nshot,ndenr(1),i1,ierlop,denrt(1:np,1), &
                  np,times,delt,i0,r1,i1,bitvl,iaved,time(1:np),ircfact)
      if(ierlop.eq.0) denrt(1:np,1)=denrt(1:np,1)*50.0
#else
      write(nttyo,*) "WARNING: density not set because MDS+ is missing"
#endif
      ierlop=0
      call avdata(nshot,ndenr(2),i1,ierlop,denrt(1:np,2), &
                  np,times,delt,i0,r1,i1,bitvl,iaved,time(1:np),ircfact, &
                  denr_rc(2),denrrcg(2),vresdenr(2), &
                  denr_k(2),t0denr(2),devdenr(1:np,2),navdenr(1:np,2), &
                  time_err)
      if(ierlop.eq.0) denrt(1:np,2)=denrt(1:np,2)*50.0
      endif s124411
!----------------------------------------------------------------------
!--   67-degree magnetic probes
!----------------------------------------------------------------------
      do i=1,magpri
        expmpi(1:np,i)=0.
        iermpi(i)=0
        sclmp=1.0
        call avdata(nshot,mpnam2(i),i1,iermpi(i),expmpi(1:np,i), &
                    np,times,delt,i0,sclmp,i1,bitmpi(i),iavem,time(1:np), &
                    ircfact,xmp_rc(i),xmprcg(i),vresxmp(i), &
                    xmp_k(i),t0xmp(i),devxmp(1:np,i),navxmp(1:np,i),time_err)
      enddo
      rnavxmp=navxmp
!--------------------------------------------------------------------
!--   New BT compensations for magnetic probes and flux loops      --
!--------------------------------------------------------------------
      have_btc: if (ibtcomp.gt.0) then
      open(unit=60,file=input_dir(1:lindir)//'btcomp.dat', &
           status='old',iostat=ioerr)
      if (ioerr.ne.0) then
        call errctrl_msg('get_consrtaints', &
                         input_dir(1:lindir)//'btcomp.dat not found')
        stop
      endif
      read_btcshot=.true.
      do while (read_btcshot)
        read (60,*,iostat=ioerr) ibtcshot, btcname
        if(is_iostat_end(ioerr)) exit
        if(ioerr.ne.0) cycle
        if (nshot.ge.ibtcshot) then
          bti322(1:np)=0.
          ierbtc=0
          call avdata(nshot,btcname,i1,ierbtc,bti322(1:np), &
                      np,times,delt,i0,r1,i1,bitbt,iavem,time(1:np), &
                      ircfact,bt_rc,btrcg,vresbt,bt_k, &
                      t0bt,devbt(1:np),navbt(1:np),time_err)
          if (ierbtc.ne.0) then
            bti322(1:np)=0.0
            exit
          endif
          do while (ioerr.eq.0)
            read (60,*,iostat=ioerr) namedum,dumbtc
            if (ioerr.ne.0) then
              read_btcshot=.false.
              exit
            endif
            do i=1,magpri
              if (mpnam2(i).eq.namedum) then
                expmpi(1:np,i)=expmpi(1:np,i)-dumbtc*bti322(1:np)
                cycle
              endif
            enddo
            do i=1,nsilop
              if (lpname(i).eq.namedum) then
                silopt(1:np,i)=silopt(1:np,i)-dumbtc*bti322(1:np)
                cycle
              endif
            enddo
          enddo
        endif
      enddo
      close(unit=60)
      endif have_btc
!---------------------------------------------------------------------
!--  correction to magnetic probes due to N1 and C Coil             --
!---------------------------------------------------------------------
      curccoi(1:np,1:mccoil)=0.
      have_coils: if (n1coil.gt.0.or.nccoil.gt.0) then
      open(unit=60,file=input_dir(1:lindir)//'ccomp.dat', &
           status='old')
!33000 read (60,*,err=34000,end=34000) ibtcshot,n1name,(ncname(i), &
!                                      i=1,3)
!---------------------------------------------  start of new section - EJS(2014)
!--  parse a line containing a shot number and unknown number of coil names
!--  assumes a format like the following:
!  156163      'N1COIL'  'RLC79'   'RLC139'   'RLC199'   'RLC259'  (etc.)
!---------------------------------------------------------------------
!33000 read(60,*,err=34000,end=34000) textline !-- read first line of file
33000 read(60,33002,err=34000,end=34000) textline !-- read first line of file
33002 format (a)
      read(textline,*,err=33000) ibtcshot !-- read shot number

      loc1 = index(textline,"'")  !-- find start of n=1 coil name
      ! if (loc1 .le. 0) go to 99  !-- check for end of string
      if (loc1 .le. 0) then
        iierr=1
        return
      endif
      textline = textline(loc1:len(textline)) !-- trim everything preceding

      read(textline,*) n1name   !-- read & drop n=1 coil name
      textline = textline(len(trim(n1name))+3:len(textline))
      nccomp = 0
      do while(index(textline,"'") .gt. 0) !-- loop to read C-coil names
        nccomp = nccomp+1   !-- count the number of coils

        loc1 = index(textline,"'")   !-- find start of coil name
        textline = textline(loc1:len(textline))!-- drop everything preceding
        ! write (6,*) 'CCOMP',textline
        read(textline,*) ncname(nccomp) !-- read & drop the coil name
        textline = textline(len(trim(ncname(nccomp)))+3:len(textline))

      enddo

      if (nccomp.gt.mccoil) then
        call errctrl_msg('get_constraints', &
                         'insufficient C-coil length (mccoil) for data')
        stop
      endif

      if(allocated(dumccc)) deallocate(dumccc)
      allocate(dumccc(nccomp))
!----------------------------------------------- end of new section - EJS(2014)
      if (nshot.ge.ibtcshot) then
!----------------------------------------------------------------------
!--     n1 Coil Current                                              --
!----------------------------------------------------------------------
        curtn1(1:np)=0.
        if (n1coil.gt.0) then
          iern1=0
          call avdata(nshot,n1name,i1,iern1,curtn1(1:np), &
                      np,times,delt,i0,r1,i1,bitn1,iavem,time(1:np), &
                      ircfact,xn1_rc,xn1rcg,vresxn1,xn1_k, &
                      t0xn1,devxn1(1:np),navxn1(1:np),time_err)
          if(iern1.ne.0) curtn1(1:np)=0.0
        endif
!----------------------------------------------------------------------
!--     C Coil Current                                               --
!----------------------------------------------------------------------
        if (nccoil.gt.0.and.nshot.ge.83350) then
          do k=1,nccomp      !EJS(2014)
            iercc=0
            call avdata(nshot,ncname(k),i1,iercc,curccoi(1:np,k), &
                        np,times,delt,i0,r1,i1,bitipc,iavem,time(1:np), &
                        ircfact,cc_rc(k),ccrcg(k),vrescc(k), &
                        cc_k(k),t0cc(k),devcc(1:np,k),navcc(1:np,k), &
                        time_err)
            if(iercc.ne.0) curccoi(1:np,k)=0.0
          enddo
          if(oldcomp) curccoi(1:np,3)=0.0
        endif
!
33200   read (60,*,err=34000,end=34000) namedum,dumbtc &
                                        ,(dumccc(k),k=1,nccomp) !EJS(2014)
        do i=1,magpri
          if (mpnam2(i).eq.namedum) then
            if (.not.oldcomp) then
              do j=1,np
                expmpi(j,i)=expmpi(j,i)-dumbtc*curtn1(j) &
                            -sum(dumccc(1:nccomp)*curccoi(j,1:nccomp)) !EJS(2014)
              enddo
            endif
            go to 33200
          endif
        enddo
!---------------------------------------------------- new section - EJS(2014)
! Compensate flux loops for coupling of individual C-coils.
! This was not needed previously, for antisymmetric (n=odd) C-coil pairs.
!-----------------------------------------------------------------------
        do i=1,nsilop
          if (lpname(i).eq.namedum) then
            if (.not.oldcomp) then
              do j=1,np
                silopt(j,i)=silopt(j,i)-dumbtc*curtn1(j) &
                            -sum(dumccc(1:nccomp)*curccoi(j,1:nccomp))
              enddo
            endif
            go to 33200
          endif
        enddo
!---------------------------------------------- end of new section - EJS(2014)
        go to 33200
      else
!------------------------------------------------------------------ EJS(2014)
! The following loop was intended to skip to the next shot number in the file,
! if the current shot number (nshot) is not in the range defined by ibtcshot.
! This loop may fail, if the number of intervening lines is not exactly  
! equal to magpol.
! ==> It is simpler and more reliable to go directly to 33000 and begin looking
! for a line that starts with a shot number.  
! The lines that need to be skipped begin with a probe name, not a shot number.
!------------------------------------------------------------------ 
!        do i=1,magpol
!          read(60,*,err=33000,end=34000) namedum,dumbtc &
!                                        ,(dumccc(k),k=1,3)
!        enddo
!------------------------------------------------------------------ 
        go to 33000
      endif
34000 continue
      close(unit=60)
      deallocate(dumccc)
      endif have_coils
!
      curc79(1:np)=curccoi(1:np,1)
      curc139(1:np)=curccoi(1:np,2)
      curc199(1:np)=curccoi(1:np,3)
      oldccomp=.false.
      if(oldcomp) oldccomp=.true.
!---------------------------------------------------------------------
!--  correction to magnetic probes due to I Coil                    --
!---------------------------------------------------------------------
      curicoi(1:np,1:micoil)=0.
      have_ic: if (nicoil.gt.0) then
      open(unit=60,file=input_dir(1:lindir)//'icomp.dat', &
           status='old')
!35000 read (60,*,err=36000,end=36000) ibtcshot,(niname(i),i=1,6)
!
!---------------------------------------------  start of new section - EJS(2014)
!--  parse a line containing a shot number and unknown number of coil names
!--  assumes a format like the following:
!  156163      'PCIU30'   'PCIU90'  'PCIU150'   'PCIL30'   'PCIL90'  (etc.)
!---------------------------------------------------------------------
!35000 read(60,*,err=36000,end=36000) textline !-- read first line of file
35000 read(60,35002,err=36000,end=36000) textline !-- read first line of file
35002 format (a)
      read(textline,*,err=35000) ibtcshot !-- read shot number

      nicomp = 0
      do while(index(textline,"'") .gt. 0) !-- loop to read I-coil names
        nicomp = nicomp+1   !-- count the number of coils

        loc1 = index(textline,"'")   !-- find start of coil name
        textline = textline(loc1:len(textline))!-- drop everything preceding
        ! write (6,*) 'ICOMP',textline
        ! read(textline,*) ncname(nicomp) !-- read & drop the coil name
        read(textline,*) niname(nicomp)    !-- read & drop the coil name
        textline = textline(len(trim(niname(nicomp)))+3:len(textline))

      enddo

      if (nicomp.gt.micoil) then
        call errctrl_msg('get_constraints', &
                         'insufficient I-coil length (micoil) for data')
        stop
      endif

      if(allocated(dumcic)) deallocate(dumcic)
      allocate(dumcic(nicomp))
!----------------------------------------------- end of new section - EJS(2014)
!      write (6,*) 'ICOMP', nshot,ibtcshot,nicomp,niname
      if (nshot.ge.ibtcshot) then
!----------------------------------------------------------------------
!--     I-Coil Currents                                                 --
!----------------------------------------------------------------------
        curicoi(1:np,1:nicomp)=0.
        if (nicoil.gt.0.and.nshot.ge.112962) then
          do k=1,nicomp      !EJS(2014)
            ieric(k)=0
            call avdata(nshot,niname(k),i1,ieric(k),curicoi(1:np,k), &
                        np,times,delt,i0,r1,i1,bitipc,iavem,time(1:np), &
                        ircfact,xic_rc(k),xicrcg(k), &
                        vresxic(k),xic_k(k),t0xic(k),devxic(1:np,k), &
                        navxic(1:np,k),time_err)
            if(ieric(k).ne.0.and.ieric(k).ne.-1) curicoi(1:np,k)=0.0
          enddo
        endif

35200   read (60,*,err=36000,end=36000) namedum,(dumcic(k),k=1,nicomp) !EJS(2014)
        do i=1,magpri
          if (mpnam2(i).eq.namedum) then
            oldexp=expmpi(1,i)
            do j=1,np
              expmpi(j,i)=expmpi(j,i)-sum(dumcic(1:nicomp)*curicoi(j,1:nicomp)) !EJS(2014)
            enddo
            go to 35200
          endif
        enddo
!---------------------------------------------------- new section - EJS(2014)
! Compensate flux loops for coupling of individual C-coils.
! This was not needed previously, for antisymmetric (n=odd) C-coil pairs.
!-----------------------------------------------------------------------
        do i=1,nsilop   
          if (lpname(i).eq.namedum) then
            do j=1,np
              silopt(j,i)=silopt(j,i)-sum(dumcic(1:nicomp)*curicoi(j,1:nicomp))
            enddo
            go to 35200
          endif
        enddo
!---------------------------------------------- end of new section - EJS(2014)
        go to 35200
      else
!------------------------------------------------------------------ EJS(2014)
! The following loop was intended to skip to the next shot number in the file,
! if the current shot number (nshot) is not in the range defined by ibtcshot.
! This loop may fail, if the number of intervening lines is not exactly  
! equal to magpol.
! ==> It is simpler and more reliable to go directly to 35000 and begin looking
! for a line that starts with a shot number.  
! The lines that need to be skipped begin with a probe name, not a shot number.
!------------------------------------------------------------------ 
!        do i=1,magpol
!          read(60,*,err=35000,end=36000) namedum,(dumcic(k),k=1,6)
!        enddo
!------------------------------------------------------------------ 
        go to 35000
      endif
36000 continue
      close(unit=60)
      deallocate(dumcic)
      endif have_ic
      curiu30(1:np)=curicoi(1:np,1)
      curiu90(1:np)=curicoi(1:np,2)
      curiu150(1:np)=curicoi(1:np,3)
      curil30(1:np)=curicoi(1:np,4)
      curil90(1:np)=curicoi(1:np,5)
      curil150(1:np)=curicoi(1:np,6)
!----------------------------------------------------------------------
!--   get toroidal B field                                           --
!----------------------------------------------------------------------
      call avdata(nshot,nsingl(3),i1,ierbto,bcentr(1:np), &
                  np,times,delt,i0,r1,i1,bitbto,iavem,time(1:np),ircfact, &
                  bc_rc,bcrcg,vresbc,bc_k,t0bc,devbc(1:np), &
                  navbc(1:np),time_err)
      if (time_err .eq. 1) then
        if (nvtime .eq. -1) then
          write(*,*) ''
          write(*,*) 'ERROR: BCOIL data unavailable'
          write(*,*) ''
        else
          write(*,*) ''
          write(*,*) 'ERROR: BCOIL data unavailable for following times'
          do j=1,nvtime-1
            write(*,'(F8.1)') vtime(j)
          enddo
        endif
        iierr = 1
        return
      endif
      rnavbc=navbc
!----------------------------------------------------------------------
!--   correct sign of toroidal magnetic field to be consistent with  --
!--   the actual sign consistent with a right-handed cylindrical     --
!--   coordinate system                01/09/86                      --
!----------------------------------------------------------------------
      bcentr(1:np)=bcentr(1:np)*tmu/rcentr*144.
      do i=1,nfcoil
        fccurt(1:np,i)=0.
        sclmp=1.0
        call avdata(nshot,fcname(i),i1,ierfc(i),fccurt(1:np,i), &
                    np,times,delt,i0,sclmp,i1,bitfc(i),iavem,time(1:np), &
                    ircfact,fc_rc(i),fcrcg(i),vresfc(i), &
                    fc_k(i),t0fc(i),devfc(1:np,i),navfc(1:np,i),time_err)
        fccurt(1:np,i)=fccurt(1:np,i)*turnfc(i)
        devfc(1:np,i) =devfc(1:np,i) *turnfc(i)
        bitfc(i)      =bitfc(i)      *turnfc(i)
        fc_k(i)       =fc_k(i)       *turnfc(i)
      enddo
      rnavfc=navfc
!----------------------------------------------------------------
!--   New E-coil connection after discharge 85700              --
!----------------------------------------------------------------
      do i=1,nesum
        if (nshot.le.85700.and.i.gt.2) cycle
        call avdata(nshot,ecname(i),i1,ierec(i),eccurt(1:np,i), &
                    np,times,delt,i0,r1,i1,bitec(i),iavem,time(1:np), &
                    ircfact,e_rc(i),ercg(i),vrese(i), &
                    e_k(i),t0e(i),deve(1:np,i),navec(1:np,i),time_err)
        if (time_err .eq. 1) then
          if (nvtime .eq. -1) then
            write(*,*) ''
            write(*,*) 'ERROR: ECOIL data unavailable'
            write(*,*) ''
          else
            write(*,*) ''
            write(*,*) &
              'ERROR: ECOIL data unavailable for following times'
            do j=1,nvtime-1
              write(*,'(F8.1)') vtime(j)
            enddo
          endif
          iierr = 1
          return
        endif
      enddo
      if (nshot.le.85700) then
        eccurt(1:np,3)=eccurt(1:np,1)
        eccurt(1:np,5)=eccurt(1:np,1)
        eccurt(1:np,4)=eccurt(1:np,2)
        eccurt(1:np,6)=eccurt(1:np,2)
      endif
      rnavec=navec
!----------------------------------------------------------------------
!--   uncompensated diamagnetic flux if compensated not available    --
!----------------------------------------------------------------------
      if (kcaldia.eq.0) then
        call avdiam(nshot,nsingl(4),i1,ierrdi,diamag(1:np), &
                    np,delt,i0,r1,i1,bitdia,iavem,time(1:np), &
                    sigdia,ierdia)
        if (ierdia(2).gt.0.and.ierdia(3).gt.0) then
          call avdata(nshot,nsingl(4),i1,ierrdi,diamag(1:np), &
                      np,times,delt,i0,r1,i1,bitdia,iavem,time(1:np), &
                      ircfact,diam_rc,diamrcg,vresdiam, &
                      diam_k,t0diam,devdiam(1:np),navdiam(1:np),time_err)
        endif
      elseif (kcaldia.eq.1) then
        call avdata(nshot,nsingl(4),i1,ierrdi,diamag(1:np), &
                    np,times,delt,i0,r1,i1,bitdia,iavem,time(1:np), &
                    ircfact,diam_rc,diamrcg,vresdiam, &
                    diam_k,t0diam,devdiam(1:np),navdiam(1:np),time_err)
      endif
      diamag(1:np)=1.0e-03*diamag(1:np)
      sigdia(1:np)=1.0e-03*abs(sigdia(1:np))
!------------------------------------------------------------------------
!--   get beam power                                                   --
!------------------------------------------------------------------------
      if (nshot.ge.53427) then
        call apdata(nshot,nsingl(6),i1,ierbim,pbinj(1:np), &
                    np,times,delt,i0,r1,i1,bitbim,iavem,time(1:np), &
                    beam_rc,beamrcg,vresbeam,beam_k, &
                    t0beam)
        if (ierbim.ne.0) then
          pbinj(1:np)=0.0
        else
          pbinj(1:np)=pbinj(1:np)*1.e+03
        endif
      endif
!
      return
      end subroutine get_constraints

!**********************************************************************
!>
!!    avdata gets data from PTDATA and optionally performs the
!!    average.
!!
!!
!!    @param nshot : shot number
!!
!!    @param name : the point name (10 ASCII characters)
!!
!!    @param mmm :
!!
!!    @param ierror : error flag
!!
!!    @param y : output data for each time
!!
!!    @param np : number of time slices
!!
!!    @param times : first time requested (in seconds)
!!
!!    @param delt : length between time slices (in seconds)
!!
!!    @param mm :
!!
!!    @param xxd :
!!
!!    @param nn :
!!
!!    @param bitvld :
!!
!!    @param kave : time window for averaging data (in milliseconds)
!!
!!    @param time : array of times requested
!!
!!    @param ircfact :
!!
!!    @param rcx :
!!
!!    @param rcgx :
!!
!!    @param vbitx :
!!
!!    @param zinhnox :
!!
!!    @param t0x :
!!
!!    @param stdevx :
!!
!!    @param navx :
!!
!!    @param ktime_err : error flag for time not found in database
!!
!**********************************************************************
      subroutine avdata(nshot,name,mmm,ierror,y, &
                        np,times,delt,mm,xxd,nn,bitvld,kave,time,ircfact, &
                        rcx,rcgx,vbitx,zinhnox,t0x,stdevx,navx, &
                        ktime_err)
      use vtime_mod
      use var_inaver
      use var_pcsys, only: do_spline_fit
      implicit none
      real*8 sevals
      integer*4, intent(in) :: nshot,mmm,np,mm,nn,kave,ircfact
      real*8, intent(in) :: time(np),delt,xxd,times
      character*10, intent(in) ::  name
      integer*4, intent(out) :: ierror,navx(np),ktime_err
      real*8, intent(out) :: y(np),bitvld,rcx,rcgx,vbitx,zinhnox, &
                             t0x,stdevx(np)
      integer*4 mave,kkk,npn,nnp,i,j,j_save
      real*8 xm5,xx,dtmin,tmin,tmax,bitvl,rcxx,rcgxx,vbitxx,zinhnoxx, &
             t0xx,dtave,delta_min,delta
      integer*4 navxx(ntims)
      real*8 w(npmax),xw(npmax),bw(ntims),cw(ntims),dw(ntims), &
             stdevxx(ntims)
      data xm5/0.00001/
!
      xx = xxd
      ierror = 1
      leave_units: if (iaveus.le.0) then
!----------------------------------------------------------------------
!--   milli-second averaging                                         --
!----------------------------------------------------------------------
      dtmin = 0.001001
      dtmin = min(dtmin,delt)
      dtmin = max(dtmin,xm5)
      mave = iabs(kave)
      do kkk=1,8
        tmin = times-(mave+10)*dtmin*kkk
        tmax = times+(np-1)*delt+(mave+10)*dtmin*kkk
        npn = (tmax-tmin)/dtmin + np*1.5
        npn = min0(npn,ntims)
        npn = max0(npn,10)
        !
        bitvl = 0.0
        if (name .ne. 'NONE      ') then !JRF
          call getdat_e &
            (nshot,name,mmm,ierror,xw,w,npn,tmin,tmax,mm,xx,bitvl,ircfact, &
            rcxx,rcgxx,vbitxx,zinhnoxx,t0xx)
          bitvld = bitvl
          rcx = rcxx
          rcgx = rcgxx
          vbitx = vbitxx
          zinhnox = zinhnoxx
          t0x = t0xx
          if((ierror .eq. -6).or.(ierror .eq. -7)) ierror = 0
        else
          ierror = 1
        endif
        !
        if(ierror .gt. 0) return
        if(npn.ge.2) go to 100
      enddo
      ierror = 1
      return
!
  100 continue
!------------------------------------------------------------------------
!--   Check valid range of time-slice data                             --
!------------------------------------------------------------------------
      ktime_err = 0
      nnp = np
      nvtime = -1
      if(time(1) .lt. xw(1)) ktime_err = 1
      if ((ktime_err .eq. 0) .and. (time(np) .gt. xw(npn))) then
        nvtime = 1
        nnp = 0
        do i=1,np
          if (time(i) .le. xw(npn)) then
            nnp = i
          else
            vtime(nvtime) = time(i)*1000.0
            nvtime = nvtime+1
            ktime_err = 1
          endif
        enddo
        if(nnp .eq. 0) nvtime = -1
      endif
      if(nnp .eq. 0) ktime_err = 1
      if (ktime_err .eq. 1) then
        ierror = 1
        return
      endif
      if (mave .ne. 0) then
        dtave = mave*dtmin*2.
        call smoothit2(xw(1:npn),w(1:npn),npn,dtave,stdevxx(1:npn), &
                       navxx(1:npn))
      endif
!
      if (do_spline_fit) then       !JRF
        call zplines(npn,xw(1:npn),w(1:npn),bw(1:npn),cw(1:npn),dw(1:npn))
        do i=1,nnp
          y(i)=sevals(npn,time(i),xw(1:npn),w(1:npn),bw(1:npn), &
                      cw(1:npn),dw(1:npn))
        enddo
      else
        do i=1,nnp
          delta_min = 1.0e30
          do j=1,npn
            delta = abs(xw(j) - time(i))
            if(delta .lt. delta_min) then
              j_save = j
              delta_min = delta
            endif
          enddo
          y(i) = w(j_save)
          stdevx(i) = stdevxx(j_save)
          navx(i) = navxx(j_save)
!          write(6,999) xw(j_save),name
!999       format(1x,'match at ',f15.8,'for ',a)
        enddo
      endif
!
      return
!
      endif leave_units
!--------------------------------------------------------------------------
!--   averaging in micro-seconds                                         --
!--------------------------------------------------------------------------
      dtmin = 0.000001
      mave = iabs(iaveus)
      do kkk=1,8
        tmin = times-(mave+100)*dtmin*kkk
        tmax = times+(np-1)*delt+(mave+100)*dtmin*kkk
        npn = (tmax-tmin)/dtmin + 1.5
        npn = (tmax-tmin)/dtmin + np*1.5
        npn = min0(npn,ntims)
        npn = max0(npn,10)
        bitvl=0.0
        if (name .ne. 'NONE      ') then
          call getdat_e &
          (nshot,name,mmm,ierror,xw,w,npn,tmin,tmax,mm,xx,bitvl,ircfact, &
          rcxx,rcgxx,vbitxx,zinhnoxx,t0xx)
          bitvld = bitvl
          rcx = rcxx
          rcgx = rcgxx
          vbitx = vbitxx
          zinhnox = zinhnoxx
          t0x = t0xx
          if((ierror .eq. -6).or.(ierror .eq. -7)) ierror = 0
        else
          ierror = 1
        endif
        if(ierror .gt. 0) return
        if(npn.ge.2) go to 1100
      enddo
      ierror = 1
      return
!
 1100 continue
!------------------------------------------------------------------------
!--   Check valid range of time-slice data                           --
!------------------------------------------------------------------------
      ktime_err = 0
      nnp = np
      nvtime = -1
      if(time(1) .lt. xw(1)) ktime_err = 1
      if ((ktime_err .eq. 0) .and. (time(np) .gt. xw(npn))) then
        nvtime = 1
        nnp = 0
        do i=1,np
          if (time(i) .le. xw(npn)) then
            nnp = i
          else
            vtime(nvtime) = time(i)*1000.0
            nvtime = nvtime+1
            ktime_err = 1
          endif
        enddo
        if(nnp .eq. 0) nvtime = -1
      endif
      if(nnp .eq. 0) ktime_err = 1
      if (ktime_err .eq. 1) then
        ierror = 1
        return
      endif
      if (mave .ne. 0) then
        dtave = mave*dtmin*2.
        call smoothit2(xw(1:npn),w(1:npn),npn,dtave,stdevxx(1:npn), &
                       navxx(1:npn))
      endif
!
      if (do_spline_fit) then !JRF
        call zplines(npn,xw(1:npn),w(1:npn),bw(1:npn),cw(1:npn),dw(1:npn))
        do i=1,nnp
          y(i)=sevals(npn,time(i),xw(1:npn),w(1:npn),bw(1:npn), &
                      cw(1:npn),dw(1:npn))
        enddo
      else
        do i=1,nnp
          delta_min = 1.0e30
          do j=1,npn
            delta = abs(xw(j) - time(i))
            if (delta .lt. delta_min) then
              j_save = j
              delta_min = delta
            endif
          enddo
          y(i) = w(j_save)
          stdevx(i) = stdevxx(j_save)
!          write(6,999) xw(j_save),name
        enddo
      endif
!
      return
      end subroutine avdata

!**********************************************************************
!>
!!    apdata gets the data and optionally performs the
!!    average
!!
!!
!!    @param nshot :
!!
!!    @param name : the point name (10 ASCII characters)
!!
!!    @param mmm :
!!
!!    @param ierror : error flag
!!
!!    @param y : output data for each time
!!
!!    @param np : number of time slices
!!
!!    @param times : first time requested (in seconds)
!!
!!    @param delt : length between time slices (in seconds)
!!
!!    @param mm :
!!
!!    @param xxd :
!!
!!    @param nn :
!!
!!    @param bitvld :
!!
!!    @param kave : time window for averaging data (in milliseconds)
!!
!!    @param time : array of times requested
!!
!!    @param rcx :
!!
!!    @param rcgx :
!!
!!    @param vbitx :
!!
!!    @param zinhnox :
!!
!!    @param t0x :
!!
!**********************************************************************
      subroutine apdata(nshot,name,mmm,ierror,y, &
                        np,times,delt,mm,xxd,nn,bitvld,kave,time, &
                        rcx,rcgx,vbitx,zinhnox,t0x)
      use vtime_mod, only: ntims,npmax
      use var_pcsys, only: do_spline_fit
      implicit none
      real*8 sevals
      integer*4, intent(in) :: nshot,mmm,np,mm,nn,kave
      real*8, intent(in) :: time(np),delt,xxd,times
      character*10, intent(in) ::  name
      integer*4, intent(out) :: ierror
      real*8, intent(out) :: y(np),bitvld,rcx,rcgx,vbitx,zinhnox,t0x
      integer*4 mave,kkk,npn,ircfac,ktime_err,nnp,i,j,j_save
      real*8 xm5,xx,dtmin,tmin,tmax,bitvl,rcxx,rcgxx,vbitxx,zinhnoxx, &
             t0xx,dtave,delta_min,delta
      real*8 w(npmax),xw(npmax),bw(ntims),cw(ntims),dw(ntims)
      data xm5/0.00001/
!
      xx = xxd
      ierror = 1
      dtmin = 0.001001
      dtmin= min(dtmin,delt)
      dtmin= max(dtmin,xm5)
      mave = iabs(kave)
      do kkk=1,8
        tmin = times-(mave+10)*dtmin*kkk
        tmax = times+(np-1)*delt+(mave+10)*dtmin*kkk
        npn = (tmax-tmin)/dtmin + 1.5
        npn = (tmax-tmin)/dtmin + np*1.5
        npn = min0(npn,ntims)
        npn = max0(npn,10)
        !
        bitvl = 0.0
        ircfac = 0
        if (name .ne. 'NONE      ') then  !JRF
          call getdat_e &
            (nshot,name,mmm,ierror,xw,w,npn,tmin,tmax,mm,xx,bitvl,ircfac, &
            rcxx,rcgxx,vbitxx,zinhnoxx,t0xx)
          bitvld=bitvl
          rcx=rcxx
          rcgx=rcgxx
          vbitx=vbitxx
          zinhnox=zinhnoxx
          t0x=t0xx
          if((ierror .eq. -6).or.(ierror .eq. -7)) ierror = 0
        else
          ierror = 1
        endif
        !
        if(ierror .gt. 0) return
        if(npn.ge.2) go to 100
      enddo
      ierror = 1
      return
!
  100 continue
!------------------------------------------------------------------------
!--   Check valid range of time-slice data                             --
!------------------------------------------------------------------------
      ktime_err = 0
      nnp = np
      if(time(1) .lt. xw(1)) ktime_err = 1
      if ((ktime_err .eq. 0) .and. (time(np) .gt. xw(npn))) then
        nnp = 0
        do i = 1, np
          if(time(i) .le. xw(npn)) nnp = i
        enddo
      endif
      if(nnp .eq. 0) ktime_err = 1
      if (ktime_err .eq. 1) then
        ierror=1
        return
      endif
      if (mave .ne. 0) then
        dtave = mave*dtmin*2.
        call smoothit(xw(1:npn),w(1:npn),npn,dtave)
      endif
!
      if (do_spline_fit) then !JRF
        call zplines(npn,xw(1:npn),w(1:npn),bw(1:npn),cw(1:npn),dw(1:npn))
        do i=1,nnp
          y(i)=sevals(npn,time(i),xw(1:npn),w(1:npn),bw(1:npn), &
                      cw(1:npn),dw(1:npn))
        enddo
      else
        do i=1,nnp
          delta_min = 1.0e30
          do j = 1,npn
            delta = abs(xw(j) - time(i))
            if (delta .lt. delta_min) then
              j_save = j
              delta_min = delta
            endif
          enddo
          y(i) = w(j_save)
        !            write(6,999) xw(j_save),name
        !999         format(1x,'match at ',f15.8,'for ',a)
        enddo
      endif
!
      return
      end subroutine apdata

! =========================================================

#if defined(USE_MDS)

!**********************************************************************
!>
!!    amdata gets the data and optionally performs the
!!    average from MDS+
!!
!!    WARNING: this subroutine uses both REAL*4 (used by MDS+) and
!!             REAL*8 variables conversions must be handled carefully
!!
!!
!!    @param nshot : shot number
!!
!!    @param name : the point name (10 ASCII characters)
!!
!!    @param mmm :
!!
!!    @param ierror : error flag
!!
!!    @param y : output data for each time
!!
!!    @param np : number of time slices
!!
!!    @param times : first time requested (in seconds)
!!
!!    @param delt : (unused)
!!
!!    @param mm :
!!
!!    @param xxd : (unused)
!!
!!    @param nn :
!!
!!    @param bitvld : (unused)
!!
!!    @param kave : time window for averaging data (in milliseconds)
!!
!!    @param time : array of times requested
!!
!!    @param ircfact :
!!
!**********************************************************************
      subroutine amdata(nshot,name,mmm,ierror,y, &
                        np,times,delt,mm,xxd,nn,bitvld,kave,time, &
                        ircfact)
      use set_kinds
      use var_pcsys, only: do_spline_fit
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      include 'mdslib.inc'
      character*10 name, MyTree
      real*8 y(np),time(np),delt,xxd,bitvld,times
      real*4, dimension(:), allocatable :: yw4,xw4
      real*8, dimension(:), allocatable :: yw,xw,bw,cw,dw,ew
      real*8 dtmin,dtave,delta_min,delta
      integer*4 :: stat,nshot,lenname,errallot,npn,mmm,ierror, &
                   np,mm,nn,kave,ircfact,ktime_err,nnp,mylen, &
                   i,j,j_save,dsc,f_dsc,t_dsc,ldum
      data dtmin/0.001001/
!
      ierror=0
      if(name .eq. 'NONE      ') return !JRF
!----------------------------------------------------------------------
!--   Get data from MDS+                                             --
!----------------------------------------------------------------------
      lenname = 0
      do i=1,len(name)
        if(name(i:i).ne.' ') lenname=lenname+1
      enddo
      stat = MdsConnect('atlas'//char(0))
      if (stat.eq.-1) then
        ierror = 1
        return
      endif
!
      MyTree = 'bci'
      stat = MdsOpen('bci'//char(0), nshot)
      if (mod(stat,2) .eq. 0) then
        ierror = 1
        return
      endif
!
      dsc = descr(IDTYPE_LONG, mylen, 0)  ! MDSPlus return a descriptor number
      stat = MdsValue('SIZE('//name(1:lenname)//')'//char(0), &
                      dsc, 0, ldum)
      if (mod(stat,2) .eq. 0) then
        ierror = 1
        return
      endif

      allocate(yw4(1:mylen),stat=errallot)
      allocate(xw4(1:mylen),stat=errallot)
      allocate(yw(1:mylen),stat=errallot)
      allocate(xw(1:mylen),stat=errallot)
      allocate(bw(1:mylen),stat=errallot)
      allocate(cw(1:mylen),stat=errallot)
      allocate(dw(1:mylen),stat=errallot)
      allocate(ew(1:mylen),stat=errallot)

! 20140905 tbt NOTE: This IF only checks the last allocate!
      if (errallot .ne. 0) then
        ierror = 1
        return
      endif
      npn=mylen
!
      f_dsc = descr(IDTYPE_FLOAT, yw4, mylen, 0) ! MDSPlus return a descriptor number
      t_dsc = descr(IDTYPE_FLOAT, xw4, mylen, 0) ! MDSPlus return a descriptor number
      stat=MdsValue(name(1:lenname)//char(0), f_dsc, 0, ldum)
      if (mod(stat,2) .eq. 0) then
        ierror = 1
        return
      endif
      stat=MdsValue('DIM_OF('//name(1:lenname)//',0)'//char(0), &
                    t_dsc, 0, ldum)
      yw(1:npn)=real(yw4(1:npn),dp)
      xw(1:npn)=real(xw4(1:npn),dp)/1000.
!
      mave = iabs(kave)
      if(ierror .gt. 0) return
      if (mylen.lt.2) then
        ierror = 1
        return
      endif
!------------------------------------------------------------------------
!--   Check valid range of time-slice data                             --
!------------------------------------------------------------------------
      ktime_err = 0
      nnp = np
      if(time(1) .lt. xw(1)) ktime_err = 1
      if ((ktime_err .eq. 0) .and. (time(np) .gt. xw(npn))) then
        nnp = 0
        do i=1,np
          if(time(i) .le. xw(npn)) nnp = i
        enddo
      endif
      if(nnp .eq. 0) ktime_err = 1
      if (ktime_err .eq. 1) then
        ierror=1
        return
      endif
      if (mave .ne. 0) then
        dtave = mave*dtmin*2.
        call smoothit(xw,yw,npn,dtave)
      endif
      !
      if (do_spline_fit) then       !JRF
        call zplines(npn,xw,yw,bw,cw,dw)
        do i=1,nnp
          y(i)=sevals(npn,time(i),xw,yw,bw,cw,dw)
        enddo
      else
        do i=1,nnp
          delta_min = 1.0e30
          do j=1,npn
            delta = abs(xw(j) - time(i))
            if (delta .lt. delta_min) then
              j_save = j
              delta_min = delta
            endif
          enddo
          y(i) = yw(j_save)
        !            write(6,999) xw(j_save),name
        !999         format(1x,'match at ',f15.8,'for ',a)
        enddo
      endif
!
      return
      end subroutine amdata

!**********************************************************************
!>
!!    gettanh gets the edge hyperbolic tangent fit parameters
!!    from MDS+
!!
!!
!!    @param ishot :
!!
!!    @param fitzts :
!!
!!    @param ktime :
!!
!!    @param time :
!!
!!    @param ztssym :
!!
!!    @param ztswid :
!!
!!    @param ptssym :
!!
!!    @param ztserr :
!!
!*********************************************************************
      subroutine gettanh(ishot,fitzts,ktime,time,ztssym,ztswid, &
                          ptssym,ztserr)
      implicit none
      integer*4, intent(in) :: ishot,ktime
      real*8, intent(in) ::  time(ktime)
      character*2, intent(in) :: fitzts
      real*8, intent(out) ::  ztssym(ktime),ztswid(ktime),ptssym(ktime)
      logical, intent(out) :: ztserr(ktime)
      integer*4 :: iishot,kktime
      character*2 ztsfit
      logical errzts(ktime)
      real*8 xw(ktime),bw(ktime),cw(ktime),dw(ktime)
!-----------------------------------------------------------------------
!--   Get edge pedestal tanh paramters                                --
!-----------------------------------------------------------------------
      if (fitzts.eq.'te') then
        ztsfit='TE'
        iishot=ishot
        kktime=ktime
        xw=time
        call Get_Mtanh_Ts(iishot,ztsfit,kktime,xw,bw,cw,dw,errzts)
        ztssym=bw
        ztswid=cw
        ptssym=dw
        ztserr=errzts
      endif
!
      return
      end subroutine gettanh

#endif

! =========================================================

!**********************************************************************
!>
!!    avdiam gets the compensated diamagnetic data and
!!    and optionally performs the average.
!!
!!
!!    @param nshot : shot number
!!
!!    @param name : the point name (10 ASCII characters)
!!
!!    @param mmm :
!!
!!    @param ierror : error flag
!!
!!    @param y : output data for each time
!!
!!    @param np : number of time slices
!!
!!    @param delt : length between time slices (in seconds)
!!
!!    @param mm :
!!
!!    @param xxd :
!!
!!    @param nn :
!!
!!    @param bitvl :
!!
!!    @param kave : time window for averaging data (in milliseconds)
!!
!!    @param time : array of times requested
!!
!!    @param sigmay :
!!
!!    @param ierdia :
!!
!**********************************************************************
      subroutine avdiam(nshot,name,mmm,ierror,y,np, &
                        delt,mm,xxd,nn,bitvl,kave,time,sigmay, &
                        ierdia)
      use vtime_mod, only: ntims
      implicit none
      real*8 sevals
      integer*4, intent(in) :: nshot,mmm,np,mm,nn,kave,ierdia(3) ! mm, mmm, and nn are unused
      real*8, intent(in) :: time(np),delt,xxd
      character*10, intent(in) :: name
      integer*4, intent(out) :: ierror
      real*8, intent(out) :: y(np),sigmay(np),bitvl
      integer*4 mave,npn,i
      real*8 xm5,dtmin,tavg,dtave
      real*8 w(ntims),xw(ntims),bw(ntims),cw(ntims),dw(ntims), &
             ew(ntims)
      data xm5/0.00001/
!
      ierror=0
      dtmin=0.001001
      dtmin=min(dtmin,delt)
      dtmin=max(dtmin,xm5)
      mave=iabs(kave)
      npn=ntims
      tavg=1.0 ! TODO: is this specific to DIII-D?
      call getdia(nshot,xw,npn,tavg,ierdia,w,ew)
      if (ierdia(2).gt.0.and.ierdia(3).gt.0) then
        ierror=1
        return
      endif
!------------------------------------------------------------------
!--   average data over mave ms                                  --
!------------------------------------------------------------------
      if (mave .ne. 0) then
        dtave=mave*dtmin*2.
        call smoothit(xw(1:npn),w(1:npn),npn,dtave)
      endif
!
      call zplines(npn,xw(1:npn),w(1:npn),bw(1:npn),cw(1:npn),dw(1:npn))
      do i=1,np
        y(i)=sevals(npn,time(i),xw(1:npn),w(1:npn),bw(1:npn), &
                    cw(1:npn),dw(1:npn))
      enddo
      call zplines(npn,xw(1:npn),ew(1:npn),bw(1:npn),cw(1:npn),dw(1:npn))
      do i=1,np
        sigmay(i)=sevals(npn,time(i),xw(1:npn),ew(1:npn),bw(1:npn), &
                    cw(1:npn),dw(1:npn))
      enddo
      return
      end subroutine avdiam

!**********************************************************************
!>
!!    This subroutine averages data over a sliding window of width
!!      timint in time twice.
!!
!!
!!    @param times : array of times that data has been sampled
!!
!!    @param datarr : array of data to be averaged
!!
!!    @param nts : number of elements in data array to apply average
!!
!!    @param timint : width of the time window to average data over
!!
!**********************************************************************
      subroutine smoothit(times,datarr,nts,timint)
      implicit none
      integer*4, intent(in) :: nts
      real*8, intent(in) :: times(nts),timint
      real*8, intent(inout) :: datarr(nts)
      integer*4 kount,mlow,mhigh,i
      real*8 dtt,val,dt,tlow,thigh
      real*8 work(nts)
      !
      if(timint .le. 0.) return
      dtt = timint*.5005
      do kount=1,2
        val = datarr(1)
        mlow = 1
        mhigh = 1
        do i=1,nts
          dt = min(dtt, (times(i)-times(1))*1.001)
          dt = min(dt, (times(nts)-times(i))*1.001)
          tlow = times(i) - dt
          thigh = times(i) + dt
          do while (times(mlow).lt.tlow .or. (mhigh.lt.nts &
                    .and. times(mhigh+1).le.thigh))
            if (times(mlow) .lt. tlow) then
              val = val - datarr(mlow)
              mlow = mlow + 1
            elseif (mhigh .lt. nts .and. times(mhigh+1) .le. thigh) then
              mhigh = mhigh + 1
              val = val + datarr(mhigh)
            endif
          enddo
          work(i)=val/(mhigh-mlow+1)
        enddo
        datarr(1:nts)=work(1:nts)
      enddo
      !
      return
      end subroutine smoothit

!**********************************************************************
!>
!!    This subroutine averages data over a sliding window of width
!!      timint in time twice and computes standard deviation of the 
!!      first average applied.
!!
!!
!!    @param times : array of times that data has been sampled
!!
!!    @param datarr : array of data to be averaged
!!
!!    @param nts : number of elements in data array to apply average
!!
!!    @param timint : width of the time window to average data over
!!
!!    @param stdev : standard deviation over the average window
!!
!!    @param nave : number of points in the average window
!!
!**********************************************************************
      subroutine smoothit2(times,datarr,nts,timint,stdev,nave)
      implicit none
      integer*4, intent(in) :: nts
      real*8, intent(in) :: times(nts),timint
      integer*4, intent(out) :: nave(nts)
      real*8, intent(out) :: stdev(nts)
      real*8, intent(inout) :: datarr(nts)
      integer*4 kount,mlow,mhigh,i
      real*8 dtt,val,val2,dt,tlow,thigh
      real*8 work(nts)
      !
      if(timint .le. 0.) return
      dtt = timint*.5005
      do kount=1,2
        val = datarr(1)
        val2 = datarr(1)**2
        mlow = 1
        mhigh = 1
        do i=1,nts
          dt = min(dtt, (times(i)-times(1))*1.001)
          dt = min(dt, (times(nts)-times(i))*1.001)
          tlow = times(i) - dt
          thigh = times(i) + dt
          do while (times(mlow).lt.tlow .or. (mhigh.lt.nts &
                    .and. times(mhigh+1).le.thigh))
            if (times(mlow) .lt. tlow) then
              val = val - datarr(mlow)
              val2 = val2 - datarr(mlow)**2
              mlow = mlow + 1
            elseif (mhigh .lt. nts .and. times(mhigh+1) .le. thigh) then
              mhigh = mhigh + 1
              val = val + datarr(mhigh)
              val2 = val2 + datarr(mhigh)**2
            endif
          enddo
          work(i)=val/(mhigh-mlow+1)
          if(kount .eq. 1) then   !-- calculate std dev based on raw data
            stdev(i) = val2/(mhigh-mlow+1)-(val/(mhigh-mlow+1))**2
            stdev(i) = sqrt(abs(stdev(i)))
            nave(i) = mhigh-mlow+1
          endif
        enddo
        datarr(1:nts)=work(1:nts)
      enddo
      !
      return
      end subroutine smoothit2

!**********************************************************************
!>
!!    the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
!!    for a cubic interpolating spline
!!    
!!    s(x) = y(i) + b(i)(x-x(i)) + c(i)(x-x(i))2 + d(i)(x-x(i))3
!!    
!!    for  x(i) .le. x .le. x(i+1)
!!
!!
!!    @param n : the number of data points or knots (n.ge.2)
!!
!!    @param x : the abscissas of the knots in strictly increasing order
!!
!!    @param y : the ordinates of the knots
!!
!!    @param b : s'(y)
!!
!!    @param c : s''(y)/2
!!
!!    @param d : s'''(y)/6 (derivative from the right)
!!
!**********************************************************************
      subroutine zplines(n, x, y, b, c, d)
      implicit none
      integer*4 n
      real*8 x(n), y(n), b(n), c(n), d(n)
      integer*4 nm1, ib, i
      real*8 t
!
      nm1 = n-1
      if( n .lt. 2 ) return
      cubic: if ( n .ge. 3 ) then
!
!     set up tridiagonal system
!
!     b = diagonal, d = offdiagonal, c = right hand side.
!
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do i = 2, nm1
        d(i) = x(i+1) - x(i)
        b(i) = 2.*(d(i-1) + d(i))
        c(i+1) = (y(i+1) - y(i))/d(i)
        c(i) = c(i+1) - c(i)
      enddo
!
!     end conditions.  third derivatives at  x(1)  and  x(n)
!     obtained from divided differences
!
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.
      c(n) = 0.
      if ( n .ne. 3 ) then
        c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
        c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
        c(1) = c(1)*d(1)**2/(x(4)-x(1))
        c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
      endif
!
!     forward elimination
!
      do i = 2, n
        t = d(i-1)/b(i-1)
        b(i) = b(i) - t*d(i-1)
        c(i) = c(i) - t*c(i-1)
      enddo
!
!     back substitution
!
      c(n) = c(n)/b(n)
      do ib = 1, nm1
        i = n-ib
        c(i) = (c(i) - d(i)*c(i+1))/b(i)
      enddo
!
!     c(i) is now the sigma(i) of the text
!
!     compute polynomial coefficients
!
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
      do i = 1, nm1
        b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))
        d(i) = (c(i+1) - c(i))/d(i)
        c(i) = 3.*c(i)
      enddo
      c(n) = 3.*c(n)
      d(n) = d(n-1)
      return
!
      endif cubic
      b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.
      d(1) = 0.
      b(2) = b(1)
      c(2) = 0.
      d(2) = 0.
      return
      end subroutine zplines

!**********************************************************************
!>
!!    this function evaluates the cubic spline function\n
!!    
!!    seval = y(i) + b(i)(u-x(i)) + c(i)(u-x(i))2 + d(i)(u-x(i))3\n
!!    
!!    where  x(i) .lt. u .lt. x(i+1), using horner's rule\n
!!    
!!    if  u .lt. x(1) then  i = 1  is used.\n
!!    if  u .ge. x(n) then  i = n  is used.\n
!!
!!    if  u  is not in the same interval as the previous call, then a
!!      binary search is performed to determine the proper interval.
!!
!!
!!    @param n : the number of data points
!!
!!    @param u : the abscissa at which the spline is to be evaluated
!!
!!    @param x : the arrays of data abscissas
!!
!!    @param y : the arrays of data ordinates
!!
!!    @param b : array of spline coefficients
!!
!!    @param c : array of spline coefficients
!!
!!    @param d : array of spline coefficients
!!
!**********************************************************************
      real*8 function sevals(n, u, x, y, b, c, d)
      implicit none
      integer*4 n
      real*8  u, x(n), y(n), b(n), c(n), d(n)
      integer*4 i, j, k
      real*8 dx
      data i/1/
      if ( i .ge. n ) i = 1
!
!  binary search
!
      if (( u .lt. x(i) ) .or. ( u .gt. x(i+1) )) then
        i = 1
        j = n+1
   20   k = (i+j)/2
        if ( u .lt. x(k) ) j = k
        if ( u .ge. x(k) ) i = k
        if ( j .gt. i+1 ) go to 20
      endif
!
!  evaluate spline
!
      dx = u - x(i)
      sevals= y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
      return
      end function sevals

! =========================================================

#if defined(USEMPI)

!**********************************************************************
!>
!!    wrapper around get_constraints subroutine to handle MPI comms
!!
!!
!!    @param nshot : shot number
!!
!!    @param times : first time requested (in seconds)
!!
!!    @param delt : length between time slices (in seconds)
!!
!!    @param ktime : number of time slices
!!
!!    @param istop : error flag
!!
!**********************************************************************
        subroutine get_constraints_mpi(nshot,times,delt,ktime,istop)
        use set_kinds
        include 'eparm.inc'
        include 'modules1.inc'
        implicit none
        include 'mpif.h'

        ! Dimension ZWORK 1-D array
        !   PSIBIT 1:nsilop
        !   BITMPI 1:magpri
        !   BITFC  1:nfcoil
        !   BITEC  1:nesum
        !   IERMPI to fix FWTMP2 1:magpri
        !   IERPSI to fix FWTSI 1:nsilop
        !   IERFC to fix FWTFC 1:nfcoil
        
        ! Number columns ZWORK2 2-D array
        !   EXPMPI 1:magpri
        !   SILOPT 1:nsilop
        !   FCCURT 1:nfcoil
        !   DENVT  1:nco2v
        !   DENRT  1:nco2r
        !   ECCURT 1:nesum
        integer*4, intent(in) :: nshot
        real*8, intent(in) :: times,delt
        integer*4, intent(inout) :: ktime
        integer*4, intent(out) :: istop
        integer*4 :: i,j,ktime_all,offset,stoff,endoff,nsize,nsize2
        integer*4, dimension(:), allocatable :: tmp1,tmp2
        double precision :: timeb_list(nproc), &
                zwork(4+nsilop+magpri+nfcoil+nesum+magpri+nsilop+nfcoil), &
                zwork2(18+magpri+nsilop+nfcoil+nco2v+nco2r+nesum,ntime)
                

#ifdef DEBUG_LEVEL1
        ! timing variables
        integer*4 :: clock0,clock1,clockmax,clockrate,ticks
        double precision :: secs
        integer*4 :: total_bytes
#endif

        zwork2 = 0.0
        nsize=4+nsilop+magpri+nfcoil+nesum+magpri+nsilop+nfcoil
        nsize2=18+magpri+nsilop+nfcoil+nco2v+nco2r+nesum
        allocate(tmp1(nproc),tmp2(nproc))

        ! NOTE : ktime contains the total number of time slices at this point
        ktime_all = ktime

        ! Process with rank == 0 gets data from PTDATA/MDS+ database by calling GETPTS
        timing: if (rank == 0) then
#ifdef DEBUG_LEVEL1
          call system_clock(count_max=clockmax,count_rate=clockrate)
          call system_clock(count=clock0)
#endif
          call get_constraints(nshot,times,delt,ktime,istop)
#ifdef DEBUG_LEVEL1
          call system_clock(count=clock1)
          ticks = clock1-clock0
          secs = real(ticks,dp)/real(clockrate,dp)
          write (*,"(' GETPTS call ',f6.2,' sec')") secs
#endif
        endif timing

        ! ALL processes get error information
        ! SIZE = SIZEOF(INT4) * (NPROC - 1) bytes
#ifdef DEBUG_LEVEL1
        total_bytes = 4*(nproc-1)
#endif
        call MPI_BCAST(istop,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        ! GETPTS error occurred if ISTOP == 1
        if (istop /= 0) return
        call MPI_BCAST(oldccomp,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

#ifdef DEBUG_LEVEL1
        timing_rank: if (rank == 0) then
          call system_clock(count_max=clockmax,count_rate=clockrate)
          call system_clock(count=clock0)
        endif timing_rank
#endif
        
        ! Each process computes local KTIME and TIMEB from distribution data
        ! NOTE : Avoids need to distribute KTIME and TIMEB information via individual MPI_SEND/MPI_RECV calls
        dist_data(:) = 0
        dist_data_displs(:) = 0
        ! (a) Distribute KTIME time steps
        i = 1
        do while (i <= ktime)
          do j=1,nproc
            if (i <= ktime) then
              dist_data(j) = dist_data(j)+1
              i = i+1
            endif
          enddo
        enddo
        ! (b) Compute displacements
        do i=2,nproc
          dist_data_displs(i) = dist_data_displs(i) + sum(dist_data(1:i-1))
        enddo
        ! Determine local KTIME and TIMEB values 
        ! NOTE : ktime has different definition from now on that it did
        ktime = dist_data(rank+1)
        timeb = times*1000.0+dist_data_displs(rank+1)*delt*1000.0
        
        ! ZWORK
        if (rank == 0) then
          ! Pack ZWORK array data
          zwork(1) = real(iavem,dp)   ! INT4  (1)
          zwork(2) = real(limitr,dp)  ! INT4  (1)
          zwork(3) = bitip         ! REAL8 (1)
          zwork(4) = rcentr        ! REAL8 (1)
          stoff = 4 + 1
          endoff = stoff + nsilop - 1
          zwork(stoff:endoff) = psibit(1:nsilop)  ! REAL8 (nsilop)
          stoff = endoff + 1
          endoff = endoff + magpri
          zwork(stoff:endoff) = bitmpi(1:magpri)  ! REAL8 (magpri)
          stoff = endoff + 1
          endoff = endoff + nfcoil
          zwork(stoff:endoff) = bitfc(1:nfcoil)   ! REAL8 (nfcoil)
          stoff = endoff + 1
          endoff = endoff + nesum
          zwork(stoff:endoff) = bitec(1:nesum)  ! REAL8 (nesum)
          stoff = endoff + 1
          endoff = endoff + magpri
          zwork(stoff:endoff) = iermpi(1:magpri) ! REAL8 (magpri)
          stoff = endoff + 1
          endoff = endoff + nsilop
          zwork(stoff:endoff) = ierpsi(1:nsilop) ! REAL8 (nsilop)
          stoff = endoff + 1
          endoff = endoff + nfcoil
          zwork(stoff:endoff) = ierfc(1:nfcoil) ! REAL8 (nfcoil)
! NOTE: all of the fwtmp2 are =1 at this point, for all ranks
! in order for the logic to be equivelant later, it is the error
! codes, iermpi that must be passed to other ranks
        endif
        ! Distribute ZWORK array to ALL processes
        ! SIZE = SIZEOF(DOUBLE) * NSIZE2 * (NPROC - 1) bytes
#ifdef DEBUG_LEVEL1
        total_bytes = total_bytes + 8*nsize*(nproc-1)
#endif
        call MPI_BCAST(zwork,nsize,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        ! Unpack ZWORK array data
        ! NOTE : Only processes with rank > 0 need to unpack data
        if (rank > 0) then
          iavem  = int(zwork(1))
          limitr = int(zwork(2))
          bitip  = zwork(3)
          rcentr = zwork(4)
          stoff = 4 + 1
          endoff = stoff + nsilop - 1
          psibit(1:nsilop) = zwork(stoff:endoff)
          stoff = endoff + 1
          endoff = endoff + magpri
          bitmpi(1:magpri) = zwork(stoff:endoff)
          stoff = endoff + 1
          endoff = endoff + nfcoil
          bitfc(1:nfcoil) = zwork(stoff:endoff)
          stoff = endoff + 1
          endoff = endoff + nesum
          bitec(1:nesum) = zwork(stoff:endoff)
          stoff = endoff + 1
          endoff = endoff + magpri
          iermpi(1:magpri) = zwork(stoff:endoff)
          stoff = endoff + 1
          endoff = endoff + nsilop
          ierpsi(1:nsilop) = zwork(stoff:endoff)
          stoff = endoff + 1
          endoff = endoff + nfcoil
          ierfc(1:nfcoil) = zwork(stoff:endoff)
        endif
        
        ! ZWORK2
        if (rank == 0) then
          ! Pack ZWORK2 array data
          zwork2(1,1:ktime_all)  = time(1:ktime_all)     ! REAL8 (ntime)
          zwork2(2,1:ktime_all)  = bcentr(1:ktime_all)   ! REAL8 (ntime)
          zwork2(3,1:ktime_all)  = pasmat(1:ktime_all)   ! REAL8 (ntime)
          zwork2(4,1:ktime_all)  = vloopt(1:ktime_all)   ! REAL8 (ntime)
          zwork2(5,1:ktime_all)  = psiref(1:ktime_all)   ! REAL8 (ntime)
          zwork2(6,1:ktime_all)  = diamag(1:ktime_all)   ! REAL8 (ntime)
          zwork2(7,1:ktime_all)  = sigdia(1:ktime_all)   ! REAL8 (ntime)
          zwork2(8,1:ktime_all)  = pbinj(1:ktime_all)    ! REAL8 (ntime)
          zwork2(9,1:ktime_all)  = curtn1(1:ktime_all)   ! REAL8 (ntime)
          zwork2(10,1:ktime_all) = curc79(1:ktime_all)   ! REAL8 (ntime)
          zwork2(11,1:ktime_all) = curc139(1:ktime_all)  ! REAL8 (ntime)
! NOTE: only adding those missing from the CURRENT (183350) k-files -MK
          zwork2(12,1:ktime_all) = curc199(1:ktime_all)  ! Addd by MK 2020.10.07
          zwork2(13,1:ktime_all) = curiu30(1:ktime_all)  ! Addd by MK 2020.10.07
          zwork2(14,1:ktime_all) = curiu90(1:ktime_all)  ! Addd by MK 2020.10.07
          zwork2(15,1:ktime_all) = curiu150(1:ktime_all) ! Addd by MK 2020.10.07
          zwork2(16,1:ktime_all) = curil30(1:ktime_all)  ! Addd by MK 2020.10.07
          zwork2(17,1:ktime_all) = curil90(1:ktime_all)  ! Addd by MK 2020.10.07
          zwork2(18,1:ktime_all) = curil150(1:ktime_all) ! Addd by MK 2020.10.07

          offset = 18
          do j=1,magpri
            zwork2(j+offset,1:ktime_all) = expmpi(1:ktime_all,j)  ! REAL8 (ntime,magpri)
          enddo
          offset = offset+magpri
          do j=1,nsilop
            zwork2(j+offset,1:ktime_all) = silopt(1:ktime_all,j)  ! REAL8 (ntime,nsilop)
          enddo
          offset = offset+nsilop
          do j=1,nfcoil
            zwork2(j+offset,1:ktime_all) = fccurt(1:ktime_all,j)  ! REAL8 (ntime,nfcoil)
          enddo
          offset = offset+nfcoil
          do j=1,nco2v
            zwork2(j+offset,1:ktime_all) = denvt(1:ktime_all,j)   ! REAL8 (ntime,nco2v)
          enddo
          offset = offset+nco2v
          do j=1,nco2r
            zwork2(j+offset,1:ktime_all) = denrt(1:ktime_all,j)   ! REAL8 (ntime,nco2r)
          enddo
          offset = offset+nco2r
          do j=1,nesum
            zwork2(j+offset,1:ktime_all) = eccurt(1:ktime_all,j)  ! REAL8 (ntime,nesum)
          enddo
        endif
        ! Distribute chunks of ZWORK2 array to processes
        tmp1(:) = dist_data(:)*nsize2
        tmp2(:) = dist_data_displs(:)*nsize2
        ! SIZE = SIZEOF(DOUBLE) * SUM(DIST_DATA(2:)) * NSIZE bytes
#ifdef DEBUG_LEVEL1
        total_bytes = total_bytes + 8*sum(dist_data(2:))*nsize2
#endif
        if (rank == 0) then
          call MPI_SCATTERV(zwork2,tmp1,tmp2,MPI_DOUBLE_PRECISION, &
                            MPI_IN_PLACE,tmp1(rank+1), &
                            MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        else
          call MPI_SCATTERV(zwork2,tmp1,tmp2,MPI_DOUBLE_PRECISION,zwork2, &
                            tmp1(rank+1),MPI_DOUBLE_PRECISION,0, &
                            MPI_COMM_WORLD,ierr)
        endif
        ! Unpack ZWORK2 array data
        ! NOTE : Only processes with rank > 0 need to unpack data
        if (rank > 0) then
          time(1:ktime)    = zwork2(1,1:ktime)
          bcentr(1:ktime)  = zwork2(2,1:ktime)
          pasmat(1:ktime)  = zwork2(3,1:ktime)
          vloopt(1:ktime)  = zwork2(4,1:ktime)
          psiref(1:ktime)  = zwork2(5,1:ktime)
          diamag(1:ktime)  = zwork2(6,1:ktime)
          sigdia(1:ktime)  = zwork2(7,1:ktime)
          pbinj(1:ktime)   = zwork2(8,1:ktime)
          curtn1(1:ktime)  = zwork2(9,1:ktime)
          curc79(1:ktime)  = zwork2(10,1:ktime)
          curc139(1:ktime) = zwork2(11,1:ktime)
! NOTE that I am only adding those missing from the CURRENT k-files -MK
          curc199(1:ktime) = zwork2(12,1:ktime) ! Addd by MK 2020.10.07
          curiu30(1:ktime) = zwork2(13,1:ktime) ! Addd by MK 2020.10.07
          curiu90(1:ktime) = zwork2(14,1:ktime) ! Addd by MK 2020.10.07
          curiu150(1:ktime)= zwork2(15,1:ktime) ! Addd by MK 2020.10.07
          curil30(1:ktime) = zwork2(16,1:ktime) ! Addd by MK 2020.10.07
          curil90(1:ktime) = zwork2(17,1:ktime) ! Addd by MK 2020.10.07
          curil150(1:ktime)= zwork2(18,1:ktime) ! Addd by MK 2020.10.07

          offset = 18
          do j=1,magpri
            expmpi(1:ktime,j) = zwork2(j+offset,1:ktime)
          enddo
          offset = offset+magpri
          do j=1,nsilop
            silopt(1:ktime,j) = zwork2(j+offset,1:ktime)
          enddo
          offset = offset+nsilop
          do j=1,nfcoil
            fccurt(1:ktime,j) = zwork2(j+offset,1:ktime)
          enddo
          offset = offset+nfcoil
          do j=1,nco2v
            denvt(1:ktime,j) = zwork2(j+offset,1:ktime)
          enddo
          offset = offset+nco2v
          do j=1,nco2r
            denrt(1:ktime,j) = zwork2(j+offset,1:ktime)
          enddo
          offset = offset+nco2r
          do j=1,nesum
            eccurt(1:ktime,j) = zwork2(j+offset,1:ktime)
          enddo
        endif
        
        deallocate(tmp1,tmp2)
        
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
#ifdef DEBUG_LEVEL1
        timing_rank0: if (rank == 0) then
          call system_clock(count=clock1)
          ticks = clock1-clock0
          secs = real(ticks,dp)/real(clockrate,dp)
          write (*,"(' GETPTS transfer ',i10,' bytes in ',f6.2,'sec')") &
                total_bytes,secs
        endif timing_rank0
#endif

      end subroutine get_constraints_mpi

#endif
