#include "config.f"
!**********************************************************************
!>
!!    getpts gets the magnetic data for use with EFIT
!!    and MFIT.
!!
!!    @param nshot : shot number
!!
!!    @param times : first time requested (in seconds)
!!
!!    @param delt : length between time slices (in seconds)
!!
!!    @param np : number of time slices
!!
!!    @param iierr :
!!
!**********************************************************************
      subroutine getpts(nshot,times,delt,np,iierr)
      use vtime_mod
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      character*10 nsingl(10),n1name,btcname &
                   ,nc79name,nc139name,ncname(mccoil),niname(micoil)  !EJS(2014)
!                         ---  or ---  ncname(mccoil),niname(micoil)  !EJS(2014)

      integer :: time_err
      character*150 textline     !EJS(2014)
      character*10,dimension(:),allocatable :: ndenv,ndenr,fcname,ecname
      character*10 namedum
      real*8 dumbtc
      real*8,dimension(:),allocatable :: dumccc,dumcic

      namelist/in3/mpnam2,xmp2,ymp2,amp2,smp2,rsi,zsi,wsi,hsi,as, &
        as2,lpname,rsisvs,vsname,turnfc,patmp2,racoil,zacoil, &
        wacoil,hacoil
!
      namelist/in4/mpnam2,lpname,vsname,nsingl,n1name,btcname, & !JRF 
                   nc79name,nc139name,btcname,ndenv,ndenr, &
                   fcname,ecname
!
      ALLOCATE(ndenv(nco2v),ndenr(nco2r),fcname(nfcoil),ecname(nesum))

      nsingl(1) = 'IP        '
      nsingl(2) ='VLOOP     '
      nsingl(3) ='BCOIL     '
      nsingl(4) ='DIAMAG3   '   !diamagnetic flux ... 
      nsingl(5) ='RL01      '     !full rogowski 
      nsingl(6)='PINJ      '     !beam injection power
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
!
      efitversion = 20201123
! NOTE this is only changed so serial/parallel k-files are identical
! no changes were made to getpts() only to getpts_mpi() - MK
!----------------------------------------------------------------------
!--   read in pointnames ...                                         --
!----------------------------------------------------------------------
      open(unit=60,file=table_di2(1:ltbdi2)//'dprobe.dat', &
           status='old')
      read(60,in3)
      close(unit=60)
!
! !JRF The if statement here tests a flag from the snap file that
!     indicates whether the pointnames should be read from an alternate
!     file.
!
      if(use_alternate_pointnames .ne. 0) then
!
         open(unit=60,file=alternate_pointname_file,status='old')
         read(60,in4)
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
      if (iierr.lt.0) then
        krl01=1
      endif
      iierr=0
      i1 = 1
      i0 = 0
      r1 = 1.
!----------------------------------------------------------------------
!--   psi-loops ...                                                  --
!----------------------------------------------------------------------
      do i=1,nsilop
        do j=1,np
          silopt(j,i)=0.
        enddo
        ierpsi(i)=0
        call avdata(nshot,lpname(i),i1,ierpsi(i),silopt(1,i), &
                    np,times,delt,i0,r1,i1,psibit(i),iavem,time, &
                    ircfact, do_spline_fit,psi_rc(i),psircg(i), &
                    vrespsi(i),psi_k(i), &
                    t0psi(i),devpsi(1,i),navpsi(1,i),time_err)
        if (ierpsi(i).eq.3) then
          iierr=1
          return
        endif
        if (i.eq.iabs(nslref)) then
          do j=1,np
            psiref(j)=silopt(j,iabs(nslref))
            silopt(j,iabs(nslref))=0.0
          enddo
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
      do j=1,np
        pasmat(j)=0.
      enddo
      ierpla=0
      call avdata(nshot,nsingl(i),i1,ierpla,pasmat(1), &
                  np,times,delt,i0,r1,i1,bitip,iavem,time,ircfact, &
                  do_spline_fit,p_rc,prcg,vresp,p_k,t0p,devp(1), &
                  navp(1),time_err)
      rnavp=REAL(navp)
      if( (use_alternate_pointnames .eq. 1) .and. &      !JRF 
          (i .eq. 1) ) then
         do j=1,np
            pasmat(j) = pasmat(j) * 0.5e6
         enddo
      endif
!
      i=2
      do j=1,np
        vloopt(j)=0.
      enddo
      ierlop=0
      call avdata(nshot,nsingl(i),i1,ierlop,vloopt(1), &
                  np,times,delt,i0,r1,i1,bitvl,iavev,time,ircfact, &
                  do_spline_fit,vl_rc,vlrcg,vresvl,vl_k,t0vl,devvl(1), &
                  navvl(1),time_err)
      do j=1,np
        vloopt(j)=vloopt(j)
      enddo
!
      if(use_alternate_pointnames .eq. 2) then    !JRF
         do j=1,np
            if( (ierpla .eq. 0) .and. (ierlop .eq. 0)) then
               pasmat(j) = pasmat(j) - vloopt(j) * 0.646/50.0 * 0.5e6
            endif
         enddo
      endif
!---------------------------------------------------------------------
!--   Get density array from PTDATA or MDS+                         --
!---------------------------------------------------------------------
      denvt = 0.0
      denrt = 0.0
      if (nshot.lt.124411) then
      do i=1,3
        ierlop=0
        call avdata(nshot,ndenv(i),i1,ierlop,denvt(1,i), &
                    np,times,delt,i0,r1,i1,bitvl,iaved,time,ircfact, &
                    do_spline_fit,denv_rc(i),denvrcg(i),vresdenv(i), &
                    denv_k(i),t0denv(i),devdenv(1,i),navdenv(1,i),time_err)
        if (ierlop.eq.0) then
          do j=1,np
            denvt(j,i)=denvt(j,i)*50.0
          enddo
        endif
      enddo
      do i=1,2
        ierlop=0
        call avdata(nshot,ndenr(i),i1,ierlop,denrt(1,i), &
                    np,times,delt,i0,r1,i1,bitvl,iaved,time,ircfact, &
                    do_spline_fit,denr_rc(i),denrrcg(i),vresdenr(i), &
                    denr_k(i),t0denr(i),devdenr(1,i),navdenr(1,i),time_err)
        if (ierlop.eq.0) then
          do j=1,np
            denrt(j,i)=denrt(j,i)*50.0
          enddo
        endif
      enddo
      else ! nshot.ge.124411
!---------------------------------------------------------------------
!--   Get density array from MDS+                                   --
!---------------------------------------------------------------------
      do i=1,3
        ierlop=0
        call amdata(nshot,ndenv(i),i1,ierlop,denvt(1,i), &
                    np,times,delt,i0,r1,i1,bitvl,iaved,time,ircfact, &
                    do_spline_fit)
        if (ierlop.eq.0) then
          do j=1,np
            denvt(j,i)=denvt(j,i)*50.0
          enddo
        endif
      enddo
      do i=1,1
        ierlop=0
        call amdata(nshot,ndenr(i),i1,ierlop,denrt(1,i), &
                    np,times,delt,i0,r1,i1,bitvl,iaved,time,ircfact, &
                    do_spline_fit)
        if (ierlop.eq.0) then
          do j=1,np
            denrt(j,i)=denrt(j,i)*50.0
          enddo
        endif
      enddo
      do i=2,2
        ierlop=0
        call avdata(nshot,ndenr(i),i1,ierlop,denrt(1,i), &
                    np,times,delt,i0,r1,i1,bitvl,iaved,time,ircfact, &
                    do_spline_fit,denr_rc(i),denrrcg(i),vresdenr(i), &
                    denr_k(i),t0denr(i),devdenr(1,i),navdenr(1,i),time_err)
        if (ierlop.eq.0) then
          do j=1,np
            denrt(j,i)=denrt(j,i)*50.0
          enddo
        endif
      enddo
      endif ! nshot.ge.124411
!----------------------------------------------------------------------
!--   67-degree magnetic probes ...                                  --
!----------------------------------------------------------------------
      do i=1,magpri
        do j=1,np
          expmpi(j,i)=0.
        enddo
        iermpi(i)=0
        sclmp=1.0
        call avdata(nshot,mpnam2(i),i1,iermpi(i),expmpi(1,i), &
                    np,times,delt,i0,sclmp,i1,bitmpi(i),iavem,time,ircfact, &
                    do_spline_fit,xmp_rc(i),xmprcg(i),vresxmp(i),xmp_k(i), &
                    t0xmp(i), devxmp(1,i),navxmp(1,i),time_err)
      enddo
      rnavxmp = navxmp
!--------------------------------------------------------------------
!--   New BT compensations for magnetic probes and flux loops      --
!--------------------------------------------------------------------
      if (ibtcomp.gt.0) then
      open(unit=60,file=input_dir(1:lindir)//'btcomp.dat', &
           status='old')
!     TODO: This code could be problematic... its intentions should be
!           checked
!31000 read (60,*,err=31000,end=32000) ibtcshot, btcname  !EJS(2014)
31000 read (60,*,end=32000) ibtcshot, btcname  !EJS(2014)

      if (nshot.ge.ibtcshot) then
        do j=1,np
          bti322(j)=0.
        enddo
        ierbtc=0
        call avdata(nshot,btcname,i1,ierbtc,bti322(1), &
                    np,times,delt,i0,r1,i1,bitbt,iavem,time,ircfact, &
                    do_spline_fit,bt_rc,btrcg,vresbt,bt_k,t0bt,devbt(1), &
                    navbt(1),time_err)
        if (ierbtc.ne.0) then
          do j=1,np
            bti322(j)=0.0
          enddo
          go to 32000
        endif
31200   read (60,*,err=32000,end=32000) namedum,dumbtc
        do i=1,magpri
         if (mpnam2(i).eq.namedum) then
           do j=1,np
             expmpi(j,i)=expmpi(j,i)-dumbtc*bti322(j)
           enddo
           go to 31200
         endif
        enddo
        do i=1,nsilop
         if (lpname(i).eq.namedum) then
           do j=1,np
             silopt(j,i)=silopt(j,i)-dumbtc*bti322(j)
           enddo
           go to 31200
         endif
        enddo
        go to 31200
      else
!------------------------------------------------------------------ EJS(2014)
! The following loop was intended to skip to the next shot number in the file,
! if the current shot number (nshot) is not in the range defined by ibtcshot.
! This loop may fail, if the number of intervening lines is not exactly  
! equal to magpol+nsilol.
! ==> It is simpler and more reliable to go directly to 31000 and begin looking
! for a line that starts with a shot number.  
! The lines that need to be skipped begin with a probe name, not a shot number.
!------------------------------------------------------------------ 
!        do i=1,magpol+nsilol
!          read(60,*,err=31000,end=32000) namedum,dumbtc
!        enddo
!------------------------------------------------------------------ 
        go to 31000
      endif
32000 continue
      close(unit=60)
      endif
!---------------------------------------------------------------------
!--  correction to magnetic probes due to N1 and C Coil             --
!---------------------------------------------------------------------
      do k=1,mccoil   ! or mccoil instead of 6 !EJS(2014)
       do j=1,np
         curccoi(j,k)=0.
       enddo
      enddo
      if (n1coil.gt.0.or.nccoil.gt.0) then
      open(unit=60,file=input_dir(1:lindir)//'ccomp.dat', &
           status='old'                                )
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
      !if (loc1 .le. 0) go to 99  !-- check for end of string
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
        !           write (6,*) 'CCOMP',textline
        read(textline,*) ncname(nccomp) !-- read & drop the coil name
        textline = textline(len(trim(ncname(nccomp)))+3:len(textline))

      enddo

      allocate(dumccc(nccomp))
!----------------------------------------------- end of new section - EJS(2014)
      if (nshot.ge.ibtcshot) then
!----------------------------------------------------------------------
!--     n1 Coil Current                                              --
!----------------------------------------------------------------------
        do j=1,np
          curtn1(j)=0.
        enddo
        if (n1coil.gt.0) then
          iern1=0
          call avdata(nshot,n1name,i1,iern1,curtn1(1), &
                      np,times,delt,i0,r1,i1,bitn1,iavem,time,ircfact, &
                      do_spline_fit,xn1_rc,xn1rcg,vresxn1,xn1_k, &
                      t0xn1,devxn1(1),navxn1(1),time_err)
          if (iern1.ne.0) then
           do j=1,np
            curtn1(j)=0.0
           enddo
          endif
        endif
!----------------------------------------------------------------------
!--     C Coil Current                                               --
!----------------------------------------------------------------------
        if (nccoil.gt.0.and.nshot.ge.83350) then
         do k=1,nccomp      !EJS(2014)
           iercc=0
           call avdata(nshot,ncname(k),i1,iercc,curccoi(1,k), &
                       np,times,delt,i0,r1,i1,bitipc,iavem,time,ircfact, &
                       do_spline_fit,cc_rc(k),ccrcg(k),vrescc(k), &
                       cc_k(k),t0cc(k), &
                       devcc(1,k),navcc(1,k),time_err)
           if (iercc.ne.0) then
            do j=1,np
              curccoi(j,k)=0.0
            enddo
           endif
         enddo
         if (oldcomp) then
           do j=1,np
             curccoi(j,3)=0.0
           enddo
         endif
        endif
!
33200   read (60,*,err=34000,end=34000) namedum,dumbtc &
                                        ,(dumccc(k),k=1,nccomp) !EJS(2014)
        do i=1,magpri
         if (mpnam2(i).eq.namedum) then
           do j=1,np
             if (.not.oldcomp) then
             expmpi(j,i)=expmpi(j,i)-dumbtc*curtn1(j)
             do k=1,nccomp     !EJS(2014)
              expmpi(j,i)=expmpi(j,i)-dumccc(k)*curccoi(j,k)
             enddo
             endif
           enddo
           go to 33200
         endif
        enddo
!---------------------------------------------------- new section - EJS(2014)
! Compensate flux loops for coupling of individual C-coils.
! This was not needed previously, for antisymmetric (n=odd) C-coil pairs.
!-----------------------------------------------------------------------
        do i=1,nsilop
         if (lpname(i).eq.namedum) then
           do j=1,np
             if (.not.oldcomp) then
             silopt(j,i)=silopt(j,i)-dumbtc*curtn1(j)
             do k=1,nccomp
              silopt(j,i)=silopt(j,i)-dumccc(k)*curccoi(j,k)
             enddo
             endif
           enddo
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
      endif
!
      do j=1,np
        curc79(j)=curccoi(j,1)
        curc139(j)=curccoi(j,2)
        curc199(j)=curccoi(j,3)
        curc259(j)=curccoi(j,4)
        curc319(j)=curccoi(j,5)
        curc19(j)=curccoi(j,6)
      enddo
      oldccomp=.false.
      if (oldcomp) oldccomp=.true.
!---------------------------------------------------------------------
!--  correction to magnetic probes due to I Coil                    --
!---------------------------------------------------------------------
      do k=1,micoil   ! or micoil instead of 12 !EJS(2014)
        do j=1,np
          curicoi(j,k)=0.
        enddo
      enddo
      if (nicoil.gt.0) then
      open(unit=60,file=input_dir(1:lindir)//'icomp.dat', &
           status='old'                                )
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
        !             write (6,*) 'ICOMP',textline
        !      read(textline,*) ncname(nicomp) !-- read & drop the coil name
        read(textline,*) niname(nicomp)    !-- read & drop the coil name
        textline = textline(len(trim(niname(nicomp)))+3:len(textline))

      enddo
      allocate(dumcic(nicomp))
!----------------------------------------------- end of new section - EJS(2014)
!      write (6,*) 'ICOMP', nshot,ibtcshot,nicomp,niname
      if (nshot.ge.ibtcshot) then
!----------------------------------------------------------------------
!--  I-Coil Currents                                                 --
!----------------------------------------------------------------------
        do k=1,nicomp      !EJS(2014)
         do j=1,np
          curicoi(j,k)=0.
         enddo
        enddo
        if (nicoil.gt.0.and.nshot.ge.112962) then
         do k=1,nicomp      !EJS(2014)
           ieric(k)=0
           call avdata(nshot,niname(k),i1,ieric(k),curicoi(1,k), &
                       np,times,delt,i0,r1,i1,bitipc,iavem,time,ircfact, &
                       do_spline_fit,xic_rc(k),xicrcg(k),vresxic(k), &
                       xic_k(k),t0xic(k),devxic(1,k),navxic(1,k),time_err)
           if (ieric(k).ne.0.and.ieric(k).ne.-1) then
            do j=1,np
              curicoi(j,k)=0.0
            enddo
           endif
         enddo
        endif

35200   read (60,*,err=36000,end=36000) namedum,(dumcic(k),k=1,nicomp)!EJS(2014)
        do i=1,magpri
         if (mpnam2(i).eq.namedum) then
           oldexp=expmpi(1,i)
           do j=1,np
             do k=1,nicomp     !EJS(2014)
              expmpi(j,i)=expmpi(j,i)-dumcic(k)*curicoi(j,k)
             enddo
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
             do k=1,nicomp
              silopt(j,i)=silopt(j,i)-dumcic(k)*curicoi(j,k)
             enddo
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
      endif
      do j=1,np
        curiu30(j)=curicoi(j,1)
        curiu90(j)=curicoi(j,2)
        curiu150(j)=curicoi(j,3)
        curil30(j)=curicoi(j,4)
        curil90(j)=curicoi(j,5)
        curil150(j)=curicoi(j,6)
        curiu210(j)=curicoi(j,7)
        curiu270(j)=curicoi(j,8)
        curiu330(j)=curicoi(j,9)
        curil210(j)=curicoi(j,10)
        curil270(j)=curicoi(j,11)
        curil330(j)=curicoi(j,12)
      enddo

!----------------------------------------------------------------------
!--        get toroidal B field                                      --
!----------------------------------------------------------------------
      call avdata(nshot,nsingl(3),i1,ierbto,bcentr, &
                  np,times,delt,i0,r1,i1,bitbto,iavem,time,ircfact, &
                  do_spline_fit,bc_rc,bcrcg,vresbc,bc_k,t0bc,devbc(1), &
                  navbc(1),time_err)
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
      do i=1,np
        bcentr(i)=bcentr(i)*tmu/rcentr*144.
      enddo
      do i=1,nfcoil
        do j=1,np
          fccurt(j,i)=0.
        enddo
        sclmp=1.0
        call avdata(nshot,fcname(i),i1,ierfc(i),fccurt(1,i), &
                    np,times,delt,i0,sclmp,i1,bitfc(i),iavem,time,ircfact, &
                    do_spline_fit,fc_rc(i),fcrcg(i),vresfc(i),fc_k(i),t0fc(i), &
                    devfc(1,i),navfc(1,i),time_err)
        do j=1,np
          fccurt(j,i)=fccurt(j,i)*turnfc(i)
          devfc(j,i) =devfc(j,i)*turnfc(i)
        enddo
        bitfc(i) =bitfc(i)*turnfc(i)
        fc_k(i)  =fc_k(i)*turnfc(i)
      enddo
      rnavfc = navfc
!----------------------------------------------------------------
!--   New E-coil connection after discharge 85700              --
!----------------------------------------------------------------
      do i=1,nesum
        if (nshot.le.85700.and.i.gt.2) cycle
        call avdata(nshot,ecname(i),i1,ierec(i),eccurt(1,i), &
                    np,times,delt,i0,r1,i1,bitec(i),iavem,time,ircfact, &
                    do_spline_fit,e_rc(i),ercg(i),vrese(i),e_k(i),t0e(i), &
                    deve(1,i),navec(1,i),time_err)
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
          do j=1,np
            eccurt(j,3)=eccurt(j,1)
            eccurt(j,5)=eccurt(j,1)
            eccurt(j,4)=eccurt(j,2)
            eccurt(j,6)=eccurt(j,2)
          enddo
      endif
      rnavec=navec
!----------------------------------------------------------------------
!--   uncompensated diamagnetic flux if compensated not available    --
!----------------------------------------------------------------------
      if (kcaldia.eq.0) then
        call avdiam(nshot,nsingl(4),i1,ierrdi,diamag(1), &
                  np,times,delt,i0,r1,i1,bitdia,iavem,time, &
                  sigdia,ierdia)
        if (ierdia(2).gt.0.and.ierdia(3).gt.0) then
          call avdata(nshot,nsingl(4),i1,ierrdi,diamag(1), &
                    np,times,delt,i0,r1,i1,bitdia,iavem,time,ircfact, &
                    do_spline_fit,diam_rc,diamrcg,vresdiam,diam_k, &
                    t0diam,devdiam(1),navdiam(1),time_err)
        endif
      endif
      if (kcaldia.eq.1) then
        call avdata(nshot,nsingl(4),i1,ierrdi,diamag(1), &
                    np,times,delt,i0,r1,i1,bitdia,iavem,time,ircfact, &
                    do_spline_fit,diam_rc,diamrcg,vresdiam,diam_k, &
                    t0diam,devdiam(1),navdiam(1),time_err)
      endif
      do i=1,np
        diamag(i)=1.0e-03*diamag(i)
        sigdia(i)=1.0e-03*abs(sigdia(i))
      enddo
!------------------------------------------------------------------------
!--   get beam power                                                   --
!------------------------------------------------------------------------
      if (nshot.ge.53427) then
        call apdata(nshot,nsingl(6),i1,ierbim,pbinj(1), &
                  np,times,delt,i0,r1,i1,bitbim,iavem,time, &
                  do_spline_fit,beam_rc,beamrcg,vresbeam,beam_k, &
                  t0beam,devbeam(1),navbeam(1))
        do i=1,np
          if (ierbim.ne.0) then
            pbinj(i)=0.0
          else
            pbinj(i)=pbinj(i)*1.e+03
          endif
        enddo
      endif
!
      return
      end


!**********************************************************************
!>
!!    avdata gets the data and optionally performs the
!!    average.
!!
!!    @param nshot : shot number
!!
!!    @param name : 
!!
!!    @param mmm :
!!
!!    @param ierror :
!!
!!    @param y :
!!
!!    @param np :
!!
!!    @param timesd :
!!
!!    @param deltd :
!!
!!    @param mm :
!!
!!    @param xxd :
!!
!!    @param nn :
!!
!!    @param bitvld :
!!
!!    @param kave :
!!
!!    @param time :
!!
!!    @param ircfact :
!!
!!    @param do_spline_fit :
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
                        np,timesd,deltd,mm,xxd,nn,bitvld,kave,time,ircfact, &
                        do_spline_fit,rcx,rcgx,vbitx,zinhnox,t0x,stdevx,navx, &
                        ktime_err)
      use eparm,only:ntime
      use vtime_mod
      parameter (ntims=8192)
      common/gggttt/w(ntims),xw(ntims),bw(ntims),cw(ntims),dw(ntims) &
                    ,ew(ntims),stdevxx(ntims),navxx(ntims)
      common/inaver/iavem,iaved,iavev,iaveus
      character*10 name
      integer, intent(out) :: ktime_err
      real*8 y(ntime),time(ntime),deltd,xxd,bitvld,timesd &
            , rcx,rcgx,vbitx,zinhnox,t0x &
            , stdevx(1)
      integer navx(1)
      data dtmin/0.001001/,xm5/0.00001/
      save dtmin,xm5
      logical*4 do_spline_fit
!
      delt=deltd
      times=timesd
      xx=xxd
      ierror=1
      if (iaveus.le.0) then
!----------------------------------------------------------------------
!--   milli-second averaging                                         --
!----------------------------------------------------------------------
      dtmin=0.001001
      dtmin= min (dtmin,delt)
      dtmin= max (dtmin,xm5)
      mave = iabs(kave)
      do kkk=1,8
        tmin = times-(mave+10)*dtmin*kkk
        tmax = times+(np-1)*delt+(mave+10)*dtmin*kkk
        npn = (tmax-tmin)/dtmin + 1.5
        npn = min0(npn,4000)
        npn = max0(npn,10)
        !
        bitvl=0.0
        if(name .ne. 'NONE      ') then !JRF
          call getdat_e &
            (nshot,name,mmm,ierror,xw,w,npn,tmin,tmax,mm,xx,bitvl,ircfact, &
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
        if (ierror .gt. 0) return
        if (npn.ge.2) go to 100
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
      if (time(1) .lt. xw(1)) then
        ktime_err = 1
      endif
      if ((ktime_err .eq. 0) .and. (time(np) .gt. xw(npn))) then
        nvtime = 1
        nnp = 0
        do i = 1, np
          if (time(i) .le. xw(npn)) then
            nnp = i
          else
            vtime(nvtime) = time(i)*1000.0
            nvtime = nvtime+1
            ktime_err = 1
          endif
        enddo
        if (nnp .eq. 0) nvtime = -1
      endif
      if (nnp .eq. 0) ktime_err = 1
      if (ktime_err .eq. 1) then
        ierror=1
        return
      endif
      if (mave .ne. 0) then
        dtave = mave*dtmin*2.
        call smoothit2(xw,w,npn,dtave,stdevxx,navxx)
      endif
!
      if (do_spline_fit) then       !JRF
         call zplines(npn,xw,w,bw,cw,dw)
         do i=1,nnp
             timenow=time(i)
             ynow=sevals(npn,timenow,xw,w,bw,cw,dw)
             y(i)=ynow
         enddo
      else
         do i=1,nnp
            timenow=time(i)
            delta_min = 1.0e30
            do j = 1,npn
               delta = abs(xw(j) - timenow)
               if(delta .lt. delta_min) then
                  j_save = j
                  delta_min = delta
               endif
            enddo
            y(i) = w(j_save)
            stdevx(i)=stdevxx(j_save)
            navx(i)=navxx(j_save)
!            write(6,999) xw(j_save),name
!999         format(1x,'match at ',f15.8,'for ',a)
         enddo
      endif
!
      return
!
      endif ! iaveus.le.0
!--------------------------------------------------------------------------
!--   averaging in micro-seconds                                         --
!--------------------------------------------------------------------------
      dtmin=0.000001
      mave = iabs(iaveus)
      do kkk=1,8
        tmin = times-(mave+100)*dtmin*kkk
        tmax = times+(np-1)*delt+(mave+100)*dtmin*kkk
        npn = (tmax-tmin)/dtmin + 1.5
        npn = min0(npn,4000)
        npn = max0(npn,10)
        bitvl=0.0
        if(name .ne. 'NONE      ') then
           call getdat_e &
           (nshot,name,mmm,ierror,xw,w,npn,tmin,tmax,mm,xx,bitvl,ircfact, &
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
        if (ierror .gt. 0) return
        if (npn.ge.2) go to 1100
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
      if (time(1) .lt. xw(1)) ktime_err = 1
      if ((ktime_err .eq. 0) .and. (time(np) .gt. xw(npn))) then
         nvtime = 1
         nnp = 0
         do i = 1, np
            if (time(i) .le. xw(npn)) then
               nnp = i
            else
               vtime(nvtime) = time(i)*1000.0
               nvtime = nvtime+1
               ktime_err = 1
            endif
         enddo
         if (nnp .eq. 0) nvtime = -1
      endif
      if (nnp .eq. 0) ktime_err = 1
      if (ktime_err .eq. 1) then
         ierror=1
         return
      endif
      if (mave .ne. 0) then
         dtave = mave*dtmin*2.
         call smoothit2(xw,w,npn,dtave,stdevxx,navxx)
      endif
!
      if (do_spline_fit) then !JRF
         call zplines(npn,xw,w,bw,cw,dw)
         do i=1,nnp
            timenow=time(i)
            ynow=sevals(npn,timenow,xw,w,bw,cw,dw)
            y(i)=ynow
         enddo
      else
         do i=1,nnp
            timenow=time(i)
            delta_min = 1.0e30
            do j = 1,npn
               delta = abs(xw(j) - timenow)
               if (delta .lt. delta_min) then
                  j_save = j
                  delta_min = delta
               endif
            enddo
            y(i) = w(j_save)
            stdevx(i) = stdevxx(j_save)
!            write(6,999) xw(j_save),name
         enddo
      endif
!
      return
      end


!**********************************************************************
!>
!!    amdata gets the data and optionally performs the
!!    average from MDS
!!
!!    @param nshot : shot number
!!
!!    @param name :
!!
!!    @param mmm :
!!
!!    @param ierror :
!!
!!    @param y :
!!
!!    @param np :
!!
!!    @param timesd :
!!
!!    @param deltd :
!!
!!    @param mm :
!!
!!    @param xxd :
!!
!!    @param nn :
!!
!!    @param bitvld :
!!
!!    @param kave :
!!
!!    @param time :
!!
!!    @param ircfact :
!!
!!    @param do_spline_fit :
!!
!**********************************************************************
      subroutine amdata(nshot,name,mmm,ierror,y, &
               np,timesd,deltd,mm,xxd,nn,bitvld,kave,time,ircfact, &
               do_spline_fit)
      include 'mdslib.inc'
      character*10 name, MyTree
      real*8 y(np),time(np),deltd,xxd,bitvld,timesd
      real*8, allocatable :: yw(:),xw(:),bw(:),cw(:),dw(:),ew(:)
      real*8 dtmin, xm5, dtave, xx, delt, times, delta_min, delta, &
             timenow, ynow
      integer :: status, nshot, lenname, errallot, npn, mmm, ierror, &
                 np, mm, nn, kave, ircfact, ktime_err, nnp, mylen, &
                 i, j, j_save, dsc, f_dsc, t_dsc
      logical*4 do_spline_fit
      data dtmin/0.001001/,xm5/0.00001/
      save dtmin,xm5
!
      delt=deltd
      times=timesd
      xx=xxd
      ierror=0
      if (name .eq. 'NONE      ') return !JRF
!----------------------------------------------------------------------
!--   Get data from MDS+                                             --
!----------------------------------------------------------------------
      lenname = 0
      do i=1,len(name)
        if (name(i:i).ne.' ') lenname=lenname+1
      enddo
      status= MdsConnect('atlas'//char(0))
      if (status.eq.-1) then
        ierror = 1
        return
      endif
!
      MyTree = 'bci'
      status = MdsOpen('bci'//char(0), nshot)
      if (mod(status,2) .eq. 0) then
        ierror = 1
        return
      endif
!
      dsc = descr(IDTYPE_LONG, mylen, 0)  ! MDSPlus return a descriptor number
      status = Mdsvalue('SIZE('//name(1:lenname)//')'//char(0), &
                        dsc, 0, 1)
      if (mod(status,2) .eq. 0) then
        ierror = 1
        return
      endif

!----------------------------------------------
!     20140905 tbt Getting with pgf90 14.6:
!       0: ALLOCATE: array already allocated
!       Putting in check and deallocation
      If (allocated(yw)) Deallocate(yw)
      If (allocated(xw)) Deallocate(xw)
      If (allocated(bw)) Deallocate(bw)
      If (allocated(cw)) Deallocate(cw)
      If (allocated(dw)) Deallocate(dw)
      If (allocated(ew)) Deallocate(ew)
!----------------------------------------------
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
      f_dsc = descr(IDTYPE_FLOAT, yw, mylen, 0) ! MDSPlus return a descriptor number
      t_dsc = descr(IDTYPE_FLOAT, xw, mylen, 0) ! MDSPlus return a descriptor number
      status=MdsValue(name(1:lenname)//char(0), f_dsc, 0, 1)
      if (mod(status,2) .eq. 0) then
        ierror = 1
        return
      endif
      status=MdsValue('DIM_OF('//name(1:lenname)//',0)'//char(0), &
                  t_dsc, 0, 1)
      xw(1:npn)=xw(1:npn)/1000.
!
      mave = iabs(kave)
      if (ierror .gt. 0) return
      if (mylen.lt.2) then
        ierror = 1
        return
      endif
!------------------------------------------------------------------------
!--   Check valid range of time-slice data                             --
!------------------------------------------------------------------------
      ktime_err = 0
      nnp = np
      if (time(1) .lt. xw(1)) ktime_err = 1
      if ((ktime_err .eq. 0) .and. (time(np) .gt. xw(npn))) then
        nnp = 0
        do i = 1, np
          if (time(i) .le. xw(npn)) nnp = i
        enddo
      endif
      if (nnp .eq. 0) ktime_err = 1
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
          timenow=time(i)
          ynow=sevals(npn,timenow,xw,yw,bw,cw,dw)
          y(i)=ynow
        enddo
      else
        do i=1,nnp
          timenow=time(i)
          delta_min = 1.0e30
          do j = 1,npn
            delta = abs(xw(j) - timenow)
            if(delta .lt. delta_min) then
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
      end


!**********************************************************************
!>
!!    apdata gets the data and optionally performs the
!!    average
!!
!!    @param nshot :
!!
!!    @param name :
!!
!!    @param mmm :
!!
!!    @param ierror :
!!
!!    @param y :
!!
!!    @param np :
!!
!!    @param timesd :
!!
!!    @param deltd :
!!
!!    @param mm :
!!
!!    @param xxd :
!!
!!    @param nn :
!!
!!    @param bitvld :
!!
!!    @param kave :
!!
!!    @param time :
!!
!!    @param do_spline_fit :
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
!**********************************************************************
      subroutine apdata(nshot,name,mmm,ierror,y, &
        np,timesd,deltd,mm,xxd,nn,bitvld,kave,time, &
        do_spline_fit,rcx,rcgx,vbitx,zinhnox,t0x,stdevx,navx)
      parameter (ntims=8192)
      common/gggttt/w(ntims),xw(ntims),bw(ntims),cw(ntims),dw(ntims) &
                    ,ew(ntims),stdevxx(ntims),navxx(ntims)
      real*8 y(np),time(np),deltd,xxd,bitvld,timesd &
             , rcx,rcgx,vbitx,zinhnox,t0x &
             , stdevx(1)
      integer navx(1)
      character*10 name
      data dtmin/0.001001/,xm5/0.00001/
      save dtmin,xm5
      logical*4 do_spline_fit
!
      delt=deltd
      times=timesd
      xx=xxd
      ierror=1
      dtmin=0.001001
      dtmin= min (dtmin,delt)
      dtmin= max (dtmin,xm5)
      mave = iabs(kave)
      do kkk=1,8
        tmin = times-(mave+10)*dtmin*kkk
        tmax = times+(np-1)*delt+(mave+10)*dtmin*kkk
        npn = (tmax-tmin)/dtmin + 1.5
        npn = min0(npn,4000)
        npn = max0(npn,10)
        !
        bitvl=0.0
        irfac = 0
        if(name .ne. 'NONE      ') then  !JRF
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
        if (ierror .gt. 0) return
        if (npn.ge.2) go to 100
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
      if (time(1) .lt. xw(1)) ktime_err = 1
      if ((ktime_err .eq. 0) .and. (time(np) .gt. xw(npn))) then
        nnp = 0
        do i = 1, np
          if (time(i) .le. xw(npn)) nnp = i
        enddo
      endif
      if (nnp .eq. 0) ktime_err = 1
      if (ktime_err .eq. 1) then
        ierror=1
        return
      endif
      if (mave .ne. 0) then
        dtave = mave*dtmin*2.
        call smoothit(xw,w,npn,dtave)
      endif
      !
      if (do_spline_fit) then       !JRF
        call zplines(npn,xw,w,bw,cw,dw)
        do i=1,nnp
          timenow=time(i)
          ynow=sevals(npn,timenow,xw,w,bw,cw,dw)
          y(i)=ynow
        enddo
      else
        do i=1,nnp
          timenow=time(i)
          delta_min = 1.0e30
          do j = 1,npn
            delta = abs(xw(j) - timenow)
            if(delta .lt. delta_min) then
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
      end


!**********************************************************************
!>
!!    gettanh gets the edge hyperbolic tangent fit parameters
!!    from MDS+
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
      parameter (ntims=8192)
      common/gggttt/w(ntims),xw(ntims),bw(ntims),cw(ntims),dw(ntims) &
                    ,ew(ntims),stdevxx(ntims),navxx(ntims)
      real*8 time(ktime),ztssym(ktime),ztswid(ktime),ptssym(ktime)
      character*2 fitzts,ztsfit
      logical ztserr(ktime)
      logical errzts(ntims)
      integer*4 iishot,kktime
!-----------------------------------------------------------------------
!--   Get edge pedestal tanh paramters                                --
!-----------------------------------------------------------------------
      if (fitzts.eq.'te') then
        ztsfit='TE'
        iishot=ishot
        kktime=ktime
        do i=1,ktime
         xw(i)=time(i)
        enddo
        call Get_Mtanh_Ts(iishot,ztsfit,kktime,xw,bw,cw,dw,errzts)
        do i=1,ktime
         ztssym(i)=bw(i)
         ztswid(i)=cw(i)
         ptssym(i)=dw(i)
         ztserr(i)=errzts(i)
        enddo
      endif
!
      return
      end


!**********************************************************************
!>
!!    avdata gets the compensated diamagnetic data and
!!    and optionally performs the average.
!!
!!    @param nshot : shot number
!!
!!    @param name :
!!
!!    @param mmm :
!!
!!    @param ierror :
!!
!!    @param y :
!!
!!    @param np :
!!
!!    @param timesd :
!!
!!    @param deltd :
!!
!!    @param mm :
!!
!!    @param xxd :
!!
!!    @param nn :
!!
!!    @param bitvl :
!!
!!    @param kave :
!!
!!    @param time :
!!
!!    @param sigmay :
!!
!!    @param ierdia :
!!
!**********************************************************************
      subroutine avdiam(nshot,name,mmm,ierror,y, &
                    np,timesd,deltd,mm,xxd,nn,bitvl,kave,time,sigmay, &
                        ierdia)
      parameter (ntims=8192)
      common/gggttt/w(ntims),xw(ntims),bw(ntims),cw(ntims),dw(ntims) &
                    ,ew(ntims),stdevxx(ntims),navxx(ntims)
      dimension ierdia(3)
      real*8 y(np),time(np),sigmay(np),deltd,xxd,timesd,bitvl,xlen(1)
      character*10 name
      data dtmin/0.001001/,xm5/0.00001/
      save dtmin,xm5
!
      delt=deltd
      ierror=0
      dtmin=0.001001
      dtmin=min(dtmin,delt)
      dtmin=max(dtmin,xm5)
      mave=iabs(kave)
      npn=ntims
      tavg=1.0
      call getdia(nshot,xw,npn,tavg,ierdia,w,ew)
      if (ierdia(2).gt.0.and.ierdia(3).gt.0) then
        ierror=1
        return
      endif
!------------------------------------------------------------------
!--   average data over mave ms                                  --
!------------------------------------------------------------------
      xlen=maxloc(xw)
      npn=xlen(1)
      if (mave.ne.0) then
        dtave=mave*dtmin*2.
        call smoothit(xw,w,npn,dtave)
      endif
!
      call zplines(npn,xw,w,bw,cw,dw)
      do i=1,np
          timenow=time(i)
          ynow=sevals(npn,timenow,xw,w,bw,cw,dw)
          y(i)=ynow
      enddo
      call zplines(npn,xw,ew,bw,cw,dw)
      do i=1,np
          timenow=time(i)
          ynow=sevals(npn,timenow,xw,ew,bw,cw,dw)
          sigmay(i)=ynow
      enddo
      return
      end


!**********************************************************************
!>
!!    zmooth smooths out the data.
!!
!!    @param y :
!!
!!    @param npts :
!!
!!    @param nave :
!!
!**********************************************************************
      subroutine zmooth(y,npts,nave)
      parameter (ntims=8192)
      common/gwork1/yave(ntims)
      dimension y(npts)
!
      if (nave.eq.0) return
      if (nave.ge.0) then
!---------------------------------------------------------------------
!--     uniform averaging                                           --
!---------------------------------------------------------------------
        imax=npts-nave
        ym=y(1)
        do m=1,npts
          jb=m-nave
          if (m.le.nave) jb=m
          je=m+nave
          if (m.gt.imax) je=m
          nt=0
          ynew=0.0
          do j=jb,je
            nt=nt+1
            ynew=ynew+y(j)
          enddo
          yave(m)=ynew/nt
        enddo
        do m=1,npts
          y(m)=yave(m)
        enddo
        return
      endif
!----------------------------------------------------------------------
!--   non-uniform weighted average                                   --
!----------------------------------------------------------------------
      nnave=iabs(nave)
      if (npts.lt.2*nnave+1) return
      do m=1,nnave
        yave(m)=0.0
        sum=0.0
        do i=1,nnave+1
          sumx=nnave-i+2
          sum=sum+sumx
          yave(m)=yave(m)+y(i+m-1)*sumx
        enddo
        yave(m)=yave(m)/sum
      enddo
      do m=nnave+1,npts-nnave
        sum=0.0
        yave(m)=0.0
        do i=1,nnave+1
          sumx=nnave-i+2
          sum=sum+sumx
          yave(m)=yave(m)+y(i+m-1)*sumx
        enddo
        do i=1,nnave
          sumx=nnave-i+1
          sum=sum+sumx
          yave(m)=yave(m)+y(m-i)*sumx
        enddo
        yave(m)=yave(m)/sum
      enddo
      do m=npts-nnave+1,npts
        sum=0.0
        yave(m)=0.0
        do i=1,nnave+1
          sumx=nnave-i+2
          sum=sum+sumx
          yave(m)=yave(m)+y(m-i+1)*sumx
        enddo
        yave(m)=yave(m)/sum
      enddo
      do m=1,npts
        y(m)=yave(m)
      enddo
      return
      end


!**********************************************************************
!>
!!    This subroutine smooths it??
!!
!!    @param times :
!!
!!    @param data :
!!
!!    @param nts :
!!
!!    @param timint :
!!
!**********************************************************************
      subroutine smoothit(times,data,nts,timint)
        use error_control, only: errctrl_msg
        parameter (ntims=8192)
        dimension work(ntims)
        dimension times(ntims),data(ntims)
        !
        if (times(nts) .ne. maxval(times)) then
          call errctrl_msg('smoothit', 'times do not match nts')
          stop
        endif
        if (timint .le. 0.) return
        dtt = timint*.5005
        do kount=1,2
          val = data(1)
          mlow = 1
          mhigh = 1
          do i=1,nts
            dt = amin1(dtt, (times(i)-times(1))*1.001)
            dt = amin1(dt, (times(nts)-times(i))*1.001)
            tlow = times(i) - dt
            thigh = times(i) + dt
10          if (times(mlow) .ge. tlow) go to 20
            val = val - data(mlow)
            mlow = mlow + 1
            go to 10
20          if (mhigh .ge. nts) go to 30
            if (times(mhigh+1) .gt. thigh) go to 30
            mhigh = mhigh + 1
            val = val + data(mhigh)
            go to 20
30          work(i)=val/(mhigh-mlow+1)
          enddo
          do i=1,nts
            data(i)=work(i)
          enddo
        enddo
        !
        return
      end


!**********************************************************************
!>
!!    This subroutine smooths it too??
!!
!!    @param times :
!!
!!    @param data :
!!
!!    @param nts :
!!
!!    @param timint :
!!
!**********************************************************************
      subroutine smoothit2(times,data,nts,timint,stdev,nave)
        use error_control, only: errctrl_msg
        parameter (ntims=8192)
        dimension work(ntims)
        dimension times(ntims),data(ntims),stdev(ntims),nave(ntims)
        !
        if (times(nts) .ne. maxval(times)) then
          call errctrl_msg('smoothit2', 'times do not match nts')
          stop
        endif
        if (timint .le. 0.) return
        dtt = timint*.5005
        do kount=1,2
          val = data(1)
          val2 = data(1)**2
          mlow = 1
          mhigh = 1
          do i=1,nts
            dt = amin1(dtt, (times(i)-times(1))*1.001)
            dt = amin1(dt, (times(nts)-times(i))*1.001)
            tlow = times(i) - dt
            thigh = times(i) + dt
10          if (times(mlow) .ge. tlow) go to 20
            val = val - data(mlow)
            val2 = val2 - data(mlow)**2
            mlow = mlow + 1
            go to 10
20          if (mhigh .ge. nts) go to 30
            if (times(mhigh+1) .gt. thigh) go to 30
            mhigh = mhigh + 1
            val = val + data(mhigh)
            val2 = val2 + data(mhigh)**2
            go to 20
30          work(i)=val/(mhigh-mlow+1)
            if(kount .eq. 1) then   !-- calculate std dev based on raw data
              stdev(i) = val2/(mhigh-mlow+1)-(val/(mhigh-mlow+1))**2
              stdev(i) = sqrt(abs(stdev(i)))
              nave(i) = mhigh-mlow+1
            endif
          enddo
          do i=1,nts
            data(i)=work(i)
          enddo
        enddo
        !
        return
      end


!**********************************************************************
!>
!!    the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
!!    for a cubic interpolating spline
!!    
!!    s(x) = y(i) + b(i)(x-x(i)) + c(i)(x-x(i))2 + d(i)(x-x(i))3
!!    
!!    for  x(i) .le. x .le. x(i+1)
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
      integer n
      real x(n), y(n), b(n), c(n), d(n)
      integer nm1, ib, i
      real t
!
      nm1 = n-1
      if ( n .lt. 2 ) return
      if ( n .ge. 3 ) then
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
      endif ! n .ge. 3
      b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.
      d(1) = 0.
      b(2) = b(1)
      c(2) = 0.
      d(2) = 0.
      return
      end


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
!!     if  u  is not in the same interval as the previous call, then a
!!    binary search is performed to determine the proper interval.
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
      real function sevals(n, u, x, y, b, c, d)
      integer n
      real  u, x(n), y(n), b(n), c(n), d(n)
      integer i, j, k
      real dx
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
      end


!*************************************************************************
!>
!!    This subroutine calculates the uncertainties for the magnetic
!!    diagnostics.  It is based on estimates described in
!!    DIII-D Physics Memo D3DPM 0202, "Estimating the Uncertainty of DIII-D
!!    Magnetic Data," by E.J. Strait (Aug. 30, 2002).
!!    
!!    The following sources of uncertainty are included. (Further explanation
!!    is given under the corresponding item numbers in the Physics Memo.)
!!    The individual terms are combined in quadrature to give
!!    the total uncertainty for each signal.
!!    
!!    1) Loop calibration   dS = a1 S
!!    2) Loop cal. - Long-term change  dS = a2 S
!!    3) Integrator calibration  dS = a3 S
!!    4) Int. cal. - Long-term change  dS = a4 S
!!    5) Integrator drift    dS = a5 K(RC/G) T
!!    6) Loop position    dS = a6 grad(S)
!!    7) Loop tilt angle   dS = a7 Bperp
!!    8) Bt pickup     dS = a8 dBt
!!    9) C-coil pickup   dS = a9 Cn
!!    10) Bp pickup in leads   dS = a10 K integral(Bperp^2)ds
!!    11) Bp pickup in ports   dS = a11 K Bport
!!    12) Detector nonlinearity  dS = a12 S/Bt
!!    13) Noise     dS = a13 (/<S^2> - <S>^2/)^0.5 / N^0.5
!!    14) Digitizer resolution  dS = a14 K(RC/G)(Vres/2) / N^0.5
!!    where
!!    S  = measured signal: flux, field, or current (in physics units)
!!    dS  = estimated uncertainty in the measured signal
!!    grad(S)  = gradient of measured flux or field (physics units/meter)
!!    N   = number of time samples averaged
!!    Vres  = one-bit resolution of the digitizer (volts)
!!    K   = inherent number (physics units/volt-second)
!!    RC  = integrator time constant (seconds)
!!    G   = integrator gain
!!    T   = time elapsed since the start of integration (seconds)
!!    dBt  = change in toroidal field since start of integration (seconds)
!!    Bperp  = poloidal field normal to axis of mag. probe or leads (Tesla)
!!    Cn  = current in a C-coil pair (Amps)
!!    integral(f)ds = integral along leads from probe to connector (meters)
!!    Bport  = poloidal field in port, perpendicular to leads (Tesla)
!!    an  = numerical coefficient: units vary with n,
!!    and values vary between types of signals
!!    
!!    Note that items 10 and 11 will not be implemented immediately due to
!!    the additional difficulty in calculating the path integral (10) and in
!!    estimating Bport outside the efit grid (11).
!!
!!    @param ishotx :  shot number
!!
!!    @param timexy : time slice
!!
!!    @param jtimex :
!!
!!    @param gradsmpx : Grad(S) of magnetic probe \n 
!!    s   = BR cost + BZ sint
!!
!!    @param gradsflx : Grad(S) of flux loop
!!
!!    @param bpermpx : B perpendicular to the magnetic probe
!!
!!    @param sigmafx : uncertainty of f coil
!!
!!    @param sigmabx : uncertainty of b coil
!!
!!    @param sigmaex : uncertainty of e coil
!!
!!    @param sigmaipx : 
!!
!!    @param sigmaflx : uncert. of flux loop
!!
!!    @param sigmampx : uncert. of magnetic probes
!!
!**********************************************************************
      subroutine magsigma(ishotx,timexy,jtimex,gradsmpx,gradsflx, &
                        bpermpx,sigmafx,sigmabx,sigmaex, &
                        sigmaipx,sigmaflx,sigmampx)

      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!
      dimension sigmafx(1),sigmabx(1)
      dimension sigmaex(1)
      dimension sigmaflx(1),gradsflx(1)
      dimension sigmampx(1),gradsmpx(1),bpermpx(1)

      !-----
      ! LOCAL VARIABLES
      !-----
      dimension dd(8)
      !vas equivalence  (dd(1), dpol) ,
      !vas 1  (dd(2), drdp) ,
      !vas 1  (dd(3), dfl) ,
      !vas 1  (dd(4), df67) ,
      !vas 1  (dd(5), dvv) ,
      !vas 1  (dd(6), dfc) ,
      !vas 1  (dd(7), dbe) ,
      !vas 1  (dd(8), dip)
      equivalence (dd(1), dpol), &
        (dd(2), drdp) , &
        (dd(3), dfl) , &
        (dd(4), df67) , &
        (dd(5), dvv) , &
        (dd(6), dfc) , &
        (dd(7), dbe) , &
        (dd(8), dip)

      ! rindex diagnostic  units

      ! 1 Pol. Mag. Probes   (Tesla)
      ! 2 RDP Mag. Probes   (Tesla)
      ! 3 F-coil Flux Loops (V-s/radian)
      ! 4 F6, F7 Flux Loops (V-s/radian)
      ! 5 Vessel Flux Loops (V-s/radian)
      ! 6 F-coil Rogowskis  (Amps)
      ! 7 E and B Rogowskis (Amps)
      ! 8 Ip Rogowskis   (Amps)

      !-----
      ! ARRAYS TO ACCUMULATE THE UNCERTAINTIES
      !-----
      real*8 :: sigmp(magpri),sigfl(nsilop),sigfc(nfcoil),sige(nesum)
      real*8 :: maskpol(magpri),maskrdp(magpri),maskff(nsilop), &
                maskinf(nsilop),maskoutf(nsilop),maskvv(nsilop)
      !-----
      ! MASKS FOR SUBSETS WITHIN ARRAYS
      !-----
!-- poloidal array probes
      maskpol = (/ 60*1., 16*0./)
!-- RDP probes
      maskrdp =  (/ 60*0., 16*1./)
!-- F-coil flux loops
      maskff = (/18*1.,7*0.,1., 0., 1.,7*0.,1.,0.,1.,6*0./)
!-- inner F-coil flux loops
      maskinf = (/ 5*1., 2*0., 7*1., 2*0., 2*1., 26*0. /)
!-- outer F-coil  loops
      maskoutf = (/5*0., 2*1., 7*0., 2*1., 2*0., &
                     7*0., 1., 0., 1., 7*0., 1., 0., 1., 6*0./)
!-- vacuum vessel flux loops
      maskvv=(/18*0.,7*1.,0., 1., 0., 7*1., 0., 1.,0.,6*1./)
      !-----
      ! INITIALIZE ARRAYS
      !-----
      sigmp = (/ (0., i=1, magpri) /)
      sigfl = (/ (0., i=1, nsilop) /)
      sigfc = (/ (0., i=1, nfcoil) /)
      sige = (/ (0., i=1, nesum)  /)
      sigb =   0.
      sigip =   0.

      !-----
      !  TESTING DIFFERENT CONTRIBUTIONS TO SIGMA
      !-----
      timex = timexy / 1000.
      if (ksigma .eq. 0 .or. ksigma .eq. 1) then
        !-------------------------------------------------------------------------
        !***** (1) LOOP CALIBRATION
        !-------------------------------------------------------------------------
        dd = (/.0019, .0019, 0., 0., 0., .0017, .0017, .003/)
        !-------------------------------------------------------------------------
        sigmp = sigmp + ( expmpi(jtimex,:) * dpol )**2
        sigfl = sigfl + 0
        sigfc = sigfc + ( fccurt(jtimex,:) * dfc  )**2
        sige  = sige  + ( eccurt(jtimex,:) * dbe  )**2
        sigb  = sigb  + ( bcentr(jtimex)   * dbe  )**2
        sigip = sigip + ( pasmat(jtimex)   * dip  )**2

      endif

      if (ksigma .eq. 0 .or. ksigma .eq. 2) then
        !-------------------------------------------------------------------------
        !***** (2) LOOP CALIBRATION: LONG-TERM
        !-------------------------------------------------------------------------
        dd = (/.002, .002, 0., 0., 0., .003, .003, .006/)
        !-------------------------------------------------------------------------
        sigmp = sigmp + ( expmpi(jtimex,:) * dpol )**2
        sigfl = sigfl + 0
        sigfc = sigfc + ( fccurt(jtimex,:) * dfc  )**2
        sige  = sige  + ( eccurt(jtimex,:) * dbe  )**2
        sigb  = sigb  + ( bcentr(jtimex)   * dbe  )**2
        sigip = sigip + ( pasmat(jtimex)   * dip  )**2

      endif

      if (ksigma .eq. 0 .or. ksigma .eq. 3) then
        !-------------------------------------------------------------------------
        !***** (3) INTEGRATOR CALIBRATION
        !-------------------------------------------------------------------------
        dd = (/.0013, .0013, .0013, .0013, .0013, .0013, .0013, .0013/)
        !-------------------------------------------------------------------------
        sigmp = sigmp + ( expmpi(jtimex,:) * dpol )**2
        sigfl = sigfl + ( silopt(jtimex,:) * dfl  )**2
        sigfc = sigfc + ( fccurt(jtimex,:) * dfc  )**2
        sige  = sige  + ( eccurt(jtimex,:) * dbe  )**2
        sigb  = sigb  + ( bcentr(jtimex)   * dbe  )**2
        sigip = sigip + ( pasmat(jtimex)   * dip  )**2

      endif


      if (ksigma .eq. 0 .or. ksigma .eq. 4) then
        !-------------------------------------------------------------------------
        !***** (4) INTEGRATOR CAL.: LONG-TERM
        !-------------------------------------------------------------------------
        dd = (/.0017, .0017, .0017, .0017, .0017, .0017, .0017, .0017/)
        !-------------------------------------------------------------------------
        sigmp = sigmp + ( expmpi(jtimex,:) * dpol )**2
        sigfl = sigfl + ( silopt(jtimex,:) * dfl  )**2
        sigfc = sigfc + ( fccurt(jtimex,:) * dfc  )**2
        sige  = sige  + ( eccurt(jtimex,:) * dbe  )**2
        sigb  = sigb  + ( bcentr(jtimex)   * dbe  )**2
        sigip = sigip + ( pasmat(jtimex)   * dip  )**2

      endif

      if (ksigma .eq. 0 .or. ksigma .eq. 5) then
        !-------------------------------------------------------------------------
        !***** (5) INTEGRATOR DRIFT
        !-------------------------------------------------------------------------
        dd = (/.0007, .0007, .0007, .0007, .0007, .0022, .0007, .0007/)
        !-------------------------------------------------------------------------
        sigmp = sigmp + (xmp_k  *xmprcg /vresxmp * (timex - t0xmp) &
          * dpol )**2
        sigfl = sigfl + ( psi_k * psircg/vrespsi * (timex - t0psi) &
          * dfl  )**2
        sigfc = sigfc + ( fc_k  * fcrcg /vresfc  * (timex - t0fc) &
          * dfc  )**2
        sige  = sige  + ( e_k   * ercg  /vrese   * (timex - t0e) &
          * dbe  )**2
        sigb  = sigb  + ( bc_k  * bcrcg /vresbc  * (timex - t0bc) &
          * dbe  )**2
        sigip = sigip + ( p_k   * prcg  /vresp   * (timex - t0p) &
          * dip  )**2

      endif

      if (ksigma .eq. 0 .or. ksigma .eq. 6) then
        !-------------------------------------------------------------------------
        !***** (6) LOOP POSITION
        !-------------------------------------------------------------------------
        dd = (/.0020, .0020, .0030, .0045, .0020, 0., 0., 0./)
        !-------------------------------------------------------------------------
        !vas f90 modifi
        !vas sigmp = sigmp + ( gradsmp(jtimex,:) *
        !vas 1  (maskpol * dpol + maskrdp * drdp ) )**2
        !vas sigfl = sigfl + ( gradsfl(jtimex,:) *
        !vas 1  (maskinf * dfl + maskoutf * df67 + maskvv * dvv ) )**2
        sigmp = sigmp + ( gradsmp(jtimex,:) * &
          (maskpol * dpol + maskrdp * drdp ) )**2
        sigfl = sigfl + ( gradsfl(jtimex,:) * &
          (maskinf * dfl + maskoutf * df67 + maskvv * dvv ) )**2
        sigfc = sigfc + 0
        sige  = sige  + 0
        sigb  = sigb  + 0
        sigip = sigip + 0

      endif

      if (ksigma .eq. 0 .or. ksigma .eq. 7) then
        !-------------------------------------------------------------------------
        !***** (7) LOOP TILT ANGLE - revised
        !-------------------------------------------------------------------------
        dd = (/.017, .017, 0., 0., 0., 0., 0., 0./)
        !-------------------------------------------------------------------------
        sigmp = sigmp + ( bpermp(jtimex,:) * dpol )**2
        sigfl = sigfl + 0
        sigfc = sigfc + 0
        sige  = sige  + 0
        sigb  = sigb  + 0
        sigip = sigip + 0

      endif

      if (ksigma .eq. 0 .or. ksigma .eq. 8) then
        !-------------------------------------------------------------------------
        !***** (8) BT PICKUP
        !-------------------------------------------------------------------------
        dd = (/.003, .003, .00044, .00044, .00044, 0., 0., 1.3e4/)
        !-------------------------------------------------------------------------
        sigmp = sigmp + ( bti322(jtimex) * dpol )**2
        sigfl = sigfl + ( bti322(jtimex) * dfl  )**2
        sigfc = sigfc + 0
        sige  = sige  + 0
        sigb  = sigb  + 0
        sigip = sigip + ( bti322(jtimex) * dip  )**2

      endif

      if (ksigma .eq. 0 .or. ksigma .eq. 9) then
        !-------------------------------------------------------------------------
        !**** (9a) C79 PICKUP
        !-------------------------------------------------------------------------
        dd  = (/2.3e-08,  1.4e-08,  5.1e-08,  5.1e-08, &
          3.6e-08, 0., 0., 0./)
        !-------------------------------------------------------------------------
        !vas f90 modifi
        sigmp = sigmp + ( curc79(jtimex) * &
          (maskpol * dpol + maskrdp * drdp ) )**2
        sigfl = sigfl + ( curc79(jtimex) * &
          (maskff * dfl   + maskvv * dvv   ) )**2

        !-------------------------------------------------------------------------
        !***** (9b) C139 PICKUP
        !-------------------------------------------------------------------------
        dd = (/8.1e-08,  2.3e-07,  4.1e-08,   4.1e-08, &
          3.2e-08, 0., 0., 0./)
        !-------------------------------------------------------------------------
        !vas f90 modifi
        sigmp = sigmp + ( curc139(jtimex) * &
          (maskpol * dpol + maskrdp * drdp ) )**2
        sigfl = sigfl + ( curc139(jtimex) * &
          (maskff * dfl   + maskvv * dvv   ) )**2

        !-------------------------------------------------------------------------
        !***** (9c) C199 PICKUP
        !-------------------------------------------------------------------------
        dd = (/3.1e-08,  1.1e-07,  5.2e-08,   5.2e-08, &
          3.9e-08, 0., 0., 0./)
        !-------------------------------------------------------------------------
        !vas f90 modifi
        sigmp = sigmp + ( curc199(jtimex) * &
          (maskpol * dpol + maskrdp * drdp ) )**2
        sigfl = sigfl + ( curc199(jtimex) * &
          (maskff * dfl   + maskvv * dvv   ) )**2
        sigfc = sigfc + 0
        sige  = sige  + 0
        sigb  = sigb  + 0
        sigip = sigip + 0

      endif

      if (ksigma .eq. 0 .or. ksigma .eq. 10) then
        !-------------------------------------------------------------------------
        !***** (10) BP PICKUP IN LEADS
        !-------------------------------------------------------------------------
        dd = (/.00016,  .00016,  0.,  0., .00016, 0., 0., 0./)
        !-------------------------------------------------------------------------
        Brleads = 0.
        sigmp = sigmp + ( Brleads *xmp_k  * dpol )**2
        sigfl = sigfl + ( Brleads * psi_k * &
          (maskff * dfl + maskvv * dvv ) )**2
        sigfc = sigfc + 0
        sige  = sige  + 0
        sigb  = sigb  + 0
        sigip = sigip + 0

      endif

      if (ksigma .eq. 0 .or. ksigma .eq. 11) then
        !-------------------------------------------------------------------------
        !***** (11) BP PICKUP IN PORT
        !-------------------------------------------------------------------------
        dd = (/.0002,  .0002,  0., 0., .0002, 0., 0., 0./)
        !-------------------------------------------------------------------------
        Bport = 0.
        sigmp = sigmp + ( Bport *xmp_k  * dpol )**2
        sigfl = sigfl + ( Bport * psi_k * &
          (maskff * dfl + maskvv * dvv ) )**2
        sigfc = sigfc + 0
        sige  = sige  + 0
        sigb  = sigb  + 0
        sigip = sigip + 0

      endif

      if (ksigma .eq. 0 .or. ksigma .eq. 12) then
        !-------------------------------------------------------------------------
        !***** (12) DETECTOR NONLINEARITY
        !-------------------------------------------------------------------------
        dd = (/.0025, .0025, 0., 0., 0., 0., 0., 0./)
        !-------------------------------------------------------------------------
        if (abs(bcentr(jtimex)) .gt. 0.1) then
          sigmp = sigmp + ( expmpi(jtimex,:) &
            / bcentr(jtimex) * dpol )**2
        endif
        sigfl = sigfl + 0
        sigfc = sigfc + 0
        sige  = sige  + 0
        sigb  = sigb  + 0
        sigip = sigip + 0

      endif

      if (ksigma .eq. 0 .or. ksigma .eq. 13) then
        !-------------------------------------------------------------------------
        !***** (13) NOISE
        !-------------------------------------------------------------------------
        dd = (/1., 1., 1., 1., 1., 1., 1., 1./)
        !-------------------------------------------------------------------------
        do i=1,magpri
          if (rnavxmp(jtimex,i).ne.0.0) &
            sigmp(i) = sigmp(i) + ( devxmp(jtimex,i) / &
            sqrt(rnavxmp(jtimex,i)) * dpol )**2
        enddo
        do i=1,nsilop
          if (rnavpsi(jtimex,i).ne.0.0) &
            sigfl(i) = sigfl(i) + ( devpsi(jtimex,i)/ &
            sqrt(rnavpsi(jtimex,i))* dfl  )**2
        enddo
        do i=1,nfcoil
          if (rnavfc(jtimex,i).ne.0.0) &
            sigfc(i) = sigfc(i) + ( devfc(jtimex,i) / &
            sqrt(rnavfc(jtimex,i)) * dfc  )**2
        enddo
        do i=1,nesum
          if (rnavec(jtimex,i).ne.0.0) &
            sige(i)  = sige(i)  + ( deve(jtimex,i)  / &
            sqrt(rnavec(jtimex,i))  * dbe  )**2
        enddo
        if (rnavbc(jtimex).ne.0.0) &
          sigb  = sigb  + ( devbc(jtimex)   / &
          sqrt(rnavbc(jtimex))   * dbe  )**2
        if (rnavp(jtimex).ne.0.0) &
          sigip = sigip + ( devp(jtimex)    / &
          sqrt(rnavp(jtimex))    * dip  )**2

      endif

      if (ksigma .eq. 0 .or. ksigma .eq. 14) then
        !-------------------------------------------------------------------------
        !***** (14) DIGITIZER RESOLUTION
        !-------------------------------------------------------------------------
        dd = (/1., 1., 1., 1., 1., 1., 1., 1./)
        !-------------------------------------------------------------------------
        do i=1,magpri
          if (rnavxmp(jtimex,i).ne.0.0) &
            sigmp(i) = sigmp(i) + (xmp_k(i)  *xmprcg(i)  /2/ &
            sqrt(rnavxmp(jtimex,i)) * dpol )**2
        enddo
        do i=1,nsilop
          if (rnavpsi(jtimex,i).ne.0.0) &
            sigfl(i) = sigfl(i) + ( psi_k(i) * psircg(i) /2/ &
            sqrt(rnavpsi(jtimex,i))* dfl  )**2
        enddo
        do i=1,nfcoil
          if (rnavfc(jtimex,i).ne.0.0) &
            sigfc(i) = sigfc(i) + ( fc_k(i)  * fcrcg(i)  /2/ &
            sqrt(rnavfc(jtimex,i)) * dfc  )**2
        enddo
        do i=1,nesum
          if (rnavec(jtimex,i).ne.0.0) &
            sige(i)  = sige(i)  + ( e_k(i)   * ercg(i)   /2/ &
            sqrt(rnavec(jtimex,i))  * dbe  )**2
        enddo
        if (rnavbc(jtimex).ne.0.0) &
          sigb  = sigb  + ( bc_k  * bcrcg  /2/ &
          sqrt(rnavbc(jtimex))   * dbe  )**2
        if (rnavp(jtimex).ne.0.0) &
          sigip = sigip + ( p_k   * prcg   /2/ &
          sqrt(rnavp(jtimex))    * dip  )**2

      !------------------------------------------------------------

      endif

      sigmamp(jtimex,:) = sqrt(sigmp)
      sigmafl(jtimex,:) = sqrt(sigfl)
      sigmaf(jtimex,:)  = sqrt(sigfc)
      sigmae(jtimex,:) = sqrt(sige)
      sigmab(jtimex)  = sqrt(sigb)
      sigmaip(jtimex) = sqrt(sigip)

      ! print *,'sigmamp'
      ! print 99,sigmamp
      ! print *,'sigmafl'
      ! print 99,sigmafl
      ! print *,'sigmaf'
      ! print 99,sigmaf
      ! print *,'sigmae'
      ! print 99,sigmae
      ! print *,'sigmab'
      ! print 99,sigmab
      ! print *,'sigmaip'
      ! print 99,sigmaip

!99    format(5e12.4)

      return
      end


! =========================================================

#if defined(USEMPI)

!**********************************************************************
!>
!!    MPI version of getpts subroutine...
!!
!!    @param nshot :
!!
!!    @param times :
!!
!!    @param delt :
!!
!!    @param ktime :
!!
!!    @param istop :
!!
!**********************************************************************
        subroutine getpts_mpi(nshot,times,delt,ktime,istop)
        use set_kinds
        include 'eparm.inc'
        include 'modules1.inc'
        implicit integer*4 (i-n), real*8 (a-h,o-z)
        include 'mpif.h'
        
        ! Number columns ZWORK 2-D array
        !   EXPMPI 1:magpri
        !   SILOPT 1:nsilop
        !   FCCURT 1:nfcoil
        !   DENVT  1:nco2v
        !   DENRT  1:nco2r
        !   ECCURT 1:nesum

        ! Dimension ZWORK2 1-D array
        !   PSIBIT 1:nsilop
        !   BITMPI 1:magpri
        !   BITFC  1:nfcoil
        !   BITEC  1:nesum ! added by MK
        !   IERMPI to fix FWTMP2 1:magpri! added by MK
        integer :: i,j,ktime_all,offset,nsize,nsize2
        integer,dimension(:),allocatable :: tmp1,tmp2
        double precision :: zwork(18+magpri+nsilop+nfcoil+nco2v+nco2r+nesum,ntime),&
                            zwork2(5+nsilop+magpri+nfcoil+nesum+magpri),timeb_list(nproc)

        ! TIMING >>>
        integer :: clock0,clock1,clockmax,clockrate,ticks
        double precision :: secs
        ! TIMING <<<
        integer :: total_bytes
        zwork(:,:) = 0.0
        nsize=18+magpri+nsilop+nfcoil+nco2v+nco2r+nesum
        nsize2=5+nsilop+magpri+nfcoil+nesum+magpri
        allocate(tmp1(nproc),tmp2(nproc))

        efitversion = 20201123

        ! Process with rank == 0 gets data from PTDATA/MDS+ database by calling GETPTS
        if (rank == 0) then
          ! NOTE : Need to retrive data for ALL times
          ktime_all = sum(dist_data)
          ! TIMING >>>
          call system_clock(count_max=clockmax,count_rate=clockrate)
          call system_clock(count=clock0)
          ! TIMING <<<
          call getpts(nshot,times,delt,ktime,istop)
          ! TIMING >>>
          call system_clock(count=clock1)
          ticks = clock1-clock0
          secs = real(ticks,dp)/real(clockrate,dp)
          write (*,"(' GETPTS call ',f6.2,' sec')") secs
          ! TIMING <<<
        endif
        ! ALL processes get error information
        ! SIZE = SIZEOF(INT4) * (NPROC - 1) bytes
        total_bytes = 4*(nproc-1)
        call MPI_BCAST(istop,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        ! GETPTS error occurred if ISTOP == 1
        if (istop /= 0) then
          return
        endif
        call MPI_BCAST(oldccomp,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

        ! TIMING >>>
        if (rank == 0) then
          call system_clock(count_max=clockmax,count_rate=clockrate)
          call system_clock(count=clock0)
        endif
        ! TIMING <<<
        
        ! Each process computes KTIME and TIMEB by computing distribution data
        ! NOTE : Avoids need to distribute KTIME and TIMEB information via individual MPI_SEND/MPI_RECV calls
        ktime_all = ktime
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
          do j=1,i-1
            dist_data_displs(i) = dist_data_displs(i) + dist_data(j)
          enddo
        enddo
        ! Determine local KTIME and TIMEB values
        ktime = dist_data(rank+1)
        timeb = times*1000.0+dist_data_displs(rank+1)*delt*1000.0
        
        ! ZWORK2
        if (rank == 0) then
          ! Pack ZWORK2 array data
          zwork2(1) = real(iavem,dp)   ! INT4  (1)
          zwork2(2) = real(limitr,dp)  ! INT4  (1)
          zwork2(3) = bitip         ! REAL8 (1)
          zwork2(4) = rcentr        ! REAL8 (1)
          if (oldccomp) then
             zwork2(5) = 1._dp
          else
             zwork2(5) = 0._dp
          endif
          offset = 5
          do i=1,nsilop
            zwork2(i+offset) = psibit(i)  ! REAL8 (nsilop)
          enddo
          offset = offset+nsilop
          do i=1,magpri
            zwork2(i+offset) = bitmpi(i)  ! REAL8 (magpri)
          enddo
          offset = offset+magpri
          do i=1,nfcoil
            zwork2(i+offset) = bitfc(i)   ! REAL8 (nfcoil)
          enddo
          offset = offset+nfcoil
          do i=1,nesum
            zwork2(i+offset) = bitec(i)  ! Added by MK 2020.10.07
          enddo
          offset = offset+nesum
          do i=1,magpri
            zwork2(i+offset) = iermpi(i) ! Added by MK 2020.10.07
! NOTE: all of the fwtmp2 are =1 at this point, for all ranks
! in order for the logic to be equivelant later, it is the error
! codes, iermpi that must be passed to other ranks
          enddo
        endif
        ! Distribute ZWORK2 array to ALL processes
        ! SIZE = SIZEOF(DOUBLE) * NSIZE2 * (NPROC - 1) bytes
        total_bytes = total_bytes + 8*nsize2*(nproc-1)
        call MPI_BCAST(zwork2,nsize2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        ! Unpack ZWORK2 array data
        ! NOTE : Only processes with rank > 0 need to unpack data
        if (rank > 0) then
          iavem  = int(zwork2(1))
          limitr = int(zwork2(2))
          bitip  = zwork2(3)
          rcentr = zwork2(4)
          if (abs(zwork2(5))<1.e-8_dp) then  
            oldccomp = .false.
          else
            oldccomp = .true.
          endif
          offset = 5
          do i=1,nsilop
            psibit(i) = zwork2(i+offset)
          enddo
          offset = offset+nsilop
          do i=1,magpri
            bitmpi(i) = zwork2(i+offset)
          enddo
          offset = offset+magpri
          do i=1,nfcoil
            bitfc(i) = zwork2(i+offset)
          enddo
          offset = offset+nfcoil
          do i=1,nesum
            bitec(i) = zwork2(i+offset) ! Added by MK 2020.10.07
          enddo
          offset = offset+nesum
          do i=1,magpri
            iermpi(i) = zwork2(i+offset) ! Added by MK 2020.10.07
          enddo
        endif
        
        ! ZWORK
        if (rank == 0) then
          ! Pack ZWORK array data
          do i=1,ktime_all
            zwork(1,i)  = time(i)     ! REAL8 (ntime)
            zwork(2,i)  = bcentr(i)   ! REAL8 (ntime)
            zwork(3,i)  = pasmat(i)   ! REAL8 (ntime)
            zwork(4,i)  = vloopt(i)   ! REAL8 (ntime)
            zwork(5,i)  = psiref(i)   ! REAL8 (ntime)
            zwork(6,i)  = diamag(i)   ! REAL8 (ntime)
            zwork(7,i)  = sigdia(i)   ! REAL8 (ntime)
            zwork(8,i)  = pbinj(i)    ! REAL8 (ntime)
            zwork(9,i)  = curtn1(i)   ! REAL8 (ntime)
            zwork(10,i) = curc79(i)   ! REAL8 (ntime)
            zwork(11,i) = curc139(i)  ! REAL8 (ntime)
! NOTE: only adding those missing from the CURRENT (183350) k-files -MK
            zwork(12,i) = curc199(i)  ! Addd by MK 2020.10.07
            zwork(13,i) = curiu30(i)  ! Addd by MK 2020.10.07
            zwork(14,i) = curiu90(i)  ! Addd by MK 2020.10.07
            zwork(15,i) = curiu150(i) ! Addd by MK 2020.10.07
            zwork(16,i) = curil30(i)  ! Addd by MK 2020.10.07
            zwork(17,i) = curil90(i)  ! Addd by MK 2020.10.07
            zwork(18,i) = curil150(i) ! Addd by MK 2020.10.07

            offset = 18
            do j=1,magpri
              zwork(j+offset,i) = expmpi(i,j)  ! REAL8 (ntime,magpri)
            enddo
            offset = offset+magpri
            do j=1,nsilop
              zwork(j+offset,i) = silopt(i,j)  ! REAL8 (ntime,nsilop)
            enddo
            offset = offset+nsilop
            do j=1,nfcoil
              zwork(j+offset,i) = fccurt(i,j)  ! REAL8 (ntime,nfcoil)
            enddo
            offset = offset+nfcoil
            do j=1,nco2v
              zwork(j+offset,i) = denvt(i,j)   ! REAL8 (ntime,nco2v)
            enddo
            offset = offset+nco2v
            do j=1,nco2r
              zwork(j+offset,i) = denrt(i,j)   ! REAL8 (ntime,nco2r)
            enddo
            offset = offset+nco2r
            do j=1,nesum
              zwork(j+offset,i) = eccurt(i,j)  ! REAL8 (ntime,nesum)
            enddo
          enddo
        endif
        ! Distribute chunks of ZWORK array to processes
        tmp1(:) = dist_data(:)*nsize
        tmp2(:) = dist_data_displs(:)*nsize
        ! SIZE = SIZEOF(DOUBLE) * SUM(DIST_DATA(2:)) * NSIZE bytes
        total_bytes = total_bytes + 8*sum(dist_data(2:))*nsize

        if (rank == 0) then
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
        if (rank > 0) then
          do i=1,ktime
            time(i)    = zwork(1,i)
            bcentr(i)  = zwork(2,i)
            pasmat(i)  = zwork(3,i)
            vloopt(i)  = zwork(4,i)
            psiref(i)  = zwork(5,i)
            diamag(i)  = zwork(6,i)
            sigdia(i)  = zwork(7,i)
            pbinj(i)   = zwork(8,i)
            curtn1(i)  = zwork(9,i)
            curc79(i)  = zwork(10,i)
            curc139(i) = zwork(11,i)
! NOTE that I am only adding those missing from the CURRENT k-files -MK
            curc199(i) = zwork(12,i) ! Addd by MK 2020.10.07
            curiu30(i) = zwork(13,i) ! Addd by MK 2020.10.07
            curiu90(i) = zwork(14,i) ! Addd by MK 2020.10.07
            curiu150(i)= zwork(15,i) ! Addd by MK 2020.10.07
            curil30(i) = zwork(16,i) ! Addd by MK 2020.10.07
            curil90(i) = zwork(17,i) ! Addd by MK 2020.10.07
            curil150(i)= zwork(18,i) ! Addd by MK 2020.10.07
        
            offset = 18
            do j=1,magpri
              expmpi(i,j) = zwork(j+offset,i)
            enddo
            offset = offset+magpri
            do j=1,nsilop
              silopt(i,j) = zwork(j+offset,i)
            enddo
            offset = offset+nsilop
            do j=1,nfcoil
              fccurt(i,j) = zwork(j+offset,i)
            enddo
            offset = offset+nfcoil
            do j=1,nco2v
              denvt(i,j) = zwork(j+offset,i)
            enddo
            offset = offset+nco2v
            do j=1,nco2r
              denrt(i,j) = zwork(j+offset,i)
            enddo
            offset = offset+nco2r
            do j=1,nesum
              eccurt(i,j) = zwork(j+offset,i)
            enddo
          enddo
        endif
        
        deallocate(tmp1,tmp2)
        
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
        ! TIMING >>>
        if (rank == 0) then
          call system_clock(count=clock1)
          ticks = clock1-clock0
          secs = real(ticks,dp)/real(clockrate,dp)
          write (*,"(' GETPTS transfer ',i10,' bytes in ',f6.2,'sec')") &
                total_bytes,secs
        endif
        ! TIMING <<<

      end subroutine getpts_mpi

!**********************************************************************
!>
!!    mpi version of getstark\n
!! NOTE : NO error condition returned
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
        integer :: i,j,ktime_all,nsize
        integer,dimension(:),allocatable :: tmp1,tmp2
        double precision :: zwork(nmtark*12,ntime)
        ! TIMING >>>
        integer :: clock0,clock1,clockmax,clockrate,ticks
        double precision :: secs
        ! TIMING <<<
        integer :: total_bytes
        nsize = 12*nmtark
        zwork(:,:) = 0.0
        allocate(tmp1(nproc),tmp2(nproc))

        ! Process with rank == 0 gets data from PTDATA/MDS+ database by calling GETSTARK
        if (rank == 0) then
          ! NOTE : Need to retrive data for ALL times
          ktime_all = sum(dist_data)
          ! TIMING >>>
          call system_clock(count_max=clockmax,count_rate=clockrate)
          call system_clock(count=clock0)
          ! TIMING <<<
          call getstark(ktime_all)
          ! TIMING >>>
          call system_clock(count=clock1)
          ticks = clock1-clock0
          secs = real(ticks,dp)/real(clockrate,dp)
          write (*,"(' GETSTARK call ',f6.2,' sec')") secs
          ! TIMING <<<
        endif

        ! TIMING >>>
        if (rank == 0) then
          call system_clock(count_max=clockmax,count_rate=clockrate)
          call system_clock(count=clock0)
        endif
        ! TIMING <<<
        
        ! Process with rank == 0 gets distributes data
        if (rank == 0) then
          ! Pack ZWORK data array
          ! NOTE : Transposing data arrays (packing as consecutive chunks into each column of ZWORK array)
          do i=1,ktime_all
            do j=1,nmtark
              zwork(j,i)           = tangam(i,j)
              zwork(j+nmtark,i)    = siggam(i,j)
              zwork(j+nmtark*2,i)  = rrgam(i,j)
              zwork(j+nmtark*3,i)  = zzgam(i,j)
              zwork(j+nmtark*4,i)  = a1gam(i,j)
              zwork(j+nmtark*5,i)  = a2gam(i,j)
              zwork(j+nmtark*6,i)  = a3gam(i,j)
              zwork(j+nmtark*7,i)  = a4gam(i,j)
              zwork(j+nmtark*8,i)  = a5gam(i,j)
              zwork(j+nmtark*9,i)  = a6gam(i,j)
              ! WARNING : A7GAM explicitly set to zero by original GETSTARK_MPI code
              zwork(j+nmtark*10,i) = a7gam(i,j)
              !zwork(j+nmtark*10,i) = 0.0
              ! NOTE : Do NOT actually need to pack this data since array explicitly set to zero and we could just set A8GAM to zero
              zwork(j+nmtark*11,i) = 0.0
            enddo
          enddo
        endif
        ! Distribute chunks of ZWORK array to processes
        ! NOTE : We need to recalculate distribution data
        tmp1(:) = dist_data(:)*nsize
        tmp2(:) = dist_data_displs(:)*nsize
        ! SIZE = SIZEOF(DOUBLE) * SUM(DIST_DATA(2:)) * NSIZE bytes
        total_bytes = 8*sum(dist_data(2:))*nsize
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
          do i=1,ktime
            do j=1,nmtark
              tangam(i,j) = zwork(j,i)
              siggam(i,j) = zwork(j+nmtark,i)
              rrgam(i,j)  = zwork(j+nmtark*2,i)
              zzgam(i,j)  = zwork(j+nmtark*3,i)
              a1gam(i,j)  = zwork(j+nmtark*4,i)
              a2gam(i,j)  = zwork(j+nmtark*5,i)
              a3gam(i,j)  = zwork(j+nmtark*6,i)
              a4gam(i,j)  = zwork(j+nmtark*7,i)
              a5gam(i,j)  = zwork(j+nmtark*8,i)
              a6gam(i,j)  = zwork(j+nmtark*9,i)
              a7gam(i,j)  = zwork(j+nmtark*10,i)
              a8gam(i,j)  = zwork(j+nmtark*11,i)
            enddo
          enddo
        endif

        ! KWAITMSE
        ! NOTE : Necessary to send KWAITMSE to ALL processes since controls main loop defined in EFITD
        ! SIZE = SIZEOF(INTEGER) * (NPROC - 1)
        total_bytes = total_bytes + 4*(nproc-1)
        call MPI_BCAST(kwaitmse,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        
        ! Distribute chunks of FWTGAM array
        ! NOTE : We need to recalculate distribution data
        ! WARNING : Uncertain if FWTGAM should be broadcast or if doing so could possible cause issues
        ! SIZE = SIZEOF(DOUBLE) * NMTARK * (NPROC - 1)
        total_bytes = total_bytes + 8*nmtark*(nproc-1)
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

        ! TEMP >>>
        ! SPATIAL_AVG_GAM
        call MPI_BCAST(spatial_avg_gam,nstark*ngam_vars*ngam_u*ngam_w, &
                       MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        ! TEMP <<<

        ! TEMP >>>
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
        ! TEMP <<<

        deallocate(tmp1,tmp2)
        
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
        ! TIMING >>>
        if (rank == 0) then
          call system_clock(count=clock1)
          ticks = clock1-clock0
          secs = real(ticks,dp)/real(clockrate,dp)
          write (*,"(' GETSTARK transfer ',i10,' bytes in ',f6.2,'sec')") &
                total_bytes,secs
        endif
        ! TIMING <<<

      end subroutine getstark_mpi

#endif
