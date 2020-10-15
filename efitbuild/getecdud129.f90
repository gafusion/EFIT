      subroutine getpts(nshot,times,delt,np,iierr)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          getpts gets the magnetic data for use with EFIT         **
!**          and MFIT.                                               **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/02/86..........revised for DIII-D                    **
!**          93/04/21..........revised for double precision option   **
!**      July 29, 2014.........revised GETDAT. No dimension(1) tbt   **
!**                   .........revised AMDATA to deallocate yw,xw,...**
!**                                                                  **
!**********************************************************************
      use vtime_mod
      include 'eparmdud129.f90'
      include 'modules1.f90'
!      include 'ecomdu1.f90'
      parameter (ntimem=ntime+10)
      common/cwork3/scrat(ntime),bscra(ntime),work(ntimem) &
                   ,cscra(ntime),dscra(ntime)
      character*10 nsingl(10),n1name,btcname &
                   ,nc79name,nc139name,ncname(mccoil),niname(micoil)  !EJS(2014)
!                         ---  or ---  ncname(mccoil),niname(micoil)  !EJS(2014)
      namelist/in3/mpnam2,xmp2,ymp2,amp2,smp2,rsi,zsi,wsi,hsi,as, &
        as2,lpname,rsisvs,vsname,turnfc,patmp2,racoil,zacoil, &
        wacoil,hacoil
!
      namelist/in4/mpnam2,lpname,vsname,nsingl,n1name,btcname, & !JRF 
                   nc79name,nc139name,btcname,ndenv,ndenr, &
                   fcname,ecname
!
      integer :: time_err
      character*150 textline     !EJS(2014)
      character*10 ndenv(nco2v),ndenr(nco2r),fcname(nfcoil), &
                   ecname(nesum),namedum
      real*8 dumccc(3),dumcic(6),dumbtc
      data nsingl(1) /'IP        '/, &
           nsingl(2) /'VLOOP     '/, &
           nsingl(3) /'BCOIL     '/, &
           nsingl(4) /'DIAMAG3   '/, &     !diamagnetic flux ... 
           nsingl(5) /'RL01      '/, &     !full rogowski 
           nsingl(6) /'PINJ      '/       !beam injection power
      data n1name   /'N1COIL    '/,nc79name/'C79       '/, &
           nc139name/'C139      '/
      data btcname  /'BTI322    '/
      data ndenv(1) /'DENV1     '/, &
           ndenv(2) /'DENV2     '/, &
           ndenv(3) /'DENV3     '/
      data ndenr(1) /'DENR0     '/,ndenr(2) /'DENMW     '/
      data fcname(1)/'F1A       '/,fcname(2)/'F2A       '/, &
           fcname(3)/'F3A       '/,fcname(4)/'F4A       '/, &
           fcname(5)/'F5A       '/,fcname(6)/'F6A       '/, &
           fcname(7)/'F7A       '/,fcname(8)/'F8A       '/, &
           fcname(9)/'F9A       '/,fcname(10)/'F1B       '/, &
           fcname(11)/'F2B       '/,fcname(12)/'F3B       '/, &
           fcname(13)/'F4B       '/,fcname(14)/'F5B       '/, &
           fcname(15)/'F6B       '/,fcname(16)/'F7B       '/, &
           fcname(17)/'F8B       '/,fcname(18)/'F9B       '/
      data ecname(1)/'ECOILA    '/,ecname(2)/'ECOILB    '/, &
           ecname(3)/'E567UP    '/,ecname(4)/'E567DN    '/, &
           ecname(5)/'E89DN     '/,ecname(6)/'E89UP     '/
      data irdata/0/,baddat/0/
!
      efitversion = 20201013 
! NOTE this is only changed so serial/parallel k-files are identical
! no changes were made to getpts() only to getpts_mpi() - MK
!----------------------------------------------------------------------
!-- read in pointnames ...                                           --
!----------------------------------------------------------------------
      open(unit=60,file=table_di2(1:ltbdi2)//'dprobe.dat', &
           status='old'                                )
      read(60,in3)
      close(unit=60)
!
! !JRF The if statement here tests a flag from the snap file that
!     indicates whether the pointnames should be read from an alternate
!     file.
!
      if(use_alternate_pointnames .ne. 0) then
!
         open(unit=60,file=alternate_pointname_file,status='old' &
                                             )
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
      do 2 i=1,np
    2   time(i)=times+delt*(i-1)
      krl01=0
      if (iierr.lt.0) then
        krl01=1
      endif
      iierr=0
      i1 = 1
      i0 = 0
      r1 = 1.
!----------------------------------------------------------------------
!--  psi-loops ...                                                   --
!----------------------------------------------------------------------
      do 10 i=1,nsilop
        do 5 j=1,np
   5    silopt(j,i)=0.
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
        do 7 j=1,np
        psiref(j)=silopt(j,iabs(nslref))
        silopt(j,iabs(nslref))=0.0
    7   continue
        endif
   10 continue
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
      do 25 j=1,np
   25   pasmat(j)=0.
      ierpla=0
      call avdata(nshot,nsingl(i),i1,ierpla,pasmat(1), &
                  np,times,delt,i0,r1,i1,bitip,iavem,time,ircfact, &
                  do_spline_fit,p_rc,prcg,vresp,p_k,t0p,devp(1), &
                  navp(1),time_err)
      rnavp=navp
      if( (use_alternate_pointnames .eq. 1) .and. &      !JRF 
          (i .eq. 1) ) then
         do j=1,np
            pasmat(j) = pasmat(j) * 0.5e6
         enddo
      endif
!
      i=2
      do 30 j=1,np
   30   vloopt(j)=0.
      ierlop=0
      call avdata(nshot,nsingl(i),i1,ierlop,vloopt(1), &
                  np,times,delt,i0,r1,i1,bitvl,iavev,time,ircfact, &
                  do_spline_fit,vl_rc,vlrcg,vresvl,vl_k,t0vl,devvl(1), &
                  navvl(1),time_err)
      do 32 j=1,np
   32   vloopt(j)=vloopt(j)
!
      if(use_alternate_pointnames .eq. 2) then    !JRF
         do j=1,np
            if( (ierpla .eq. 0) .and. (ierlop .eq. 0)) then
               pasmat(j) = pasmat(j) - vloopt(j) * 0.646/50.0 * 0.5e6
            endif
         enddo
      endif
!---------------------------------------------------------------------
!-- Get density array from PTDATA or MDS+                           --
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
        do 38 j=1,np
          denvt(j,i)=denvt(j,i)*50.0
   38   continue
        endif
      enddo
      do i=1,2
        ierlop=0
        call avdata(nshot,ndenr(i),i1,ierlop,denrt(1,i), &
                    np,times,delt,i0,r1,i1,bitvl,iaved,time,ircfact, &
                    do_spline_fit,denr_rc(i),denrrcg(i),vresdenr(i), &
                    denr_k(i),t0denr(i),devdenr(1,i),navdenr(1,i),time_err)
        if (ierlop.eq.0) then
        do 43 j=1,np
          denrt(j,i)=denrt(j,i)*50.0
   43   continue
        endif
      enddo
      else
!---------------------------------------------------------------------
!-- Get density array from MDS+                                     --
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
      endif
!----------------------------------------------------------------------
!-- 67-degree magnetic probes ...                                    --
!----------------------------------------------------------------------
      do 65 i=1,magpri
        do 62 j=1,np
   62   expmpi(j,i)=0.
        iermpi(i)=0
        sclmp=1.0
        call avdata(nshot,mpnam2(i),i1,iermpi(i),expmpi(1,i), &
                    np,times,delt,i0,sclmp,i1,bitmpi(i),iavem,time,ircfact, &
                    do_spline_fit,xmp_rc(i),xmprcg(i),vresxmp(i),xmp_k(i), &
                    t0xmp(i), devxmp(1,i),navxmp(1,i),time_err)
   65 continue
      rnavxmp = navxmp
!--------------------------------------------------------------------
!--      New BT compensations for magnetic probes and flux loops   --
!--------------------------------------------------------------------
      if (ibtcomp.gt.0) then
      open(unit=60,file=input_dir(1:lindir)//'btcomp.dat_156456', &
           status='old'                                )
31000 read (60,*,err=31000,end=32000) ibtcshot, btcname  !EJS(2014)

      if (nshot.ge.ibtcshot) then
        do 29 j=1,np
29        bti322(j)=0.
        ierbtc=0
        call avdata(nshot,btcname,i1,ierbtc,bti322(1), &
                    np,times,delt,i0,r1,i1,bitbt,iavem,time,ircfact, &
                    do_spline_fit,bt_rc,btrcg,vresbt,bt_k,t0bt,devbt(1), &
                    navbt(1),time_err)
        if (ierbtc.ne.0) then
          do 291 j=1,np
            bti322(j)=0.0
291       continue
          goto 32000
        endif
31200   read (60,*,err=32000,end=32000) namedum,dumbtc
        do i=1,magpri
         if (mpnam2(i).eq.namedum) then
           do j=1,np
             expmpi(j,i)=expmpi(j,i)-dumbtc*bti322(j)
           enddo
           goto 31200
         endif
        enddo
        do i=1,nsilop
         if (lpname(i).eq.namedum) then
           do j=1,np
             silopt(j,i)=silopt(j,i)-dumbtc*bti322(j)
           enddo
           goto 31200
         endif
        enddo
        goto 31200
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
        goto 31000
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
      open(unit=60,file=input_dir(1:lindir)//'ccomp.dat_168823', &
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
!----------------------------------------------- end of new section - EJS(2014)

!     write (6,*) 'CCOMP', nshot,ibtcshot,n1name,nccomp,ncname

      if (nshot.ge.ibtcshot) then
!----------------------------------------------------------------------
!--  n1 Coil Current                                                 --
!----------------------------------------------------------------------
        do 27 j=1,np
   27     curtn1(j)=0.
        if (n1coil.gt.0) then
          iern1=0
          call avdata(nshot,n1name,i1,iern1,curtn1(1), &
                      np,times,delt,i0,r1,i1,bitn1,iavem,time,ircfact, &
                      do_spline_fit,xn1_rc,xn1rcg,vresxn1,xn1_k, &
                      t0xn1,devxn1(1),navxn1(1),time_err)
          if (iern1.ne.0) then
           do 28 j=1,np
            curtn1(j)=0.0
   28      continue
          endif
        endif
!----------------------------------------------------------------------
!--  C Coil Current                                                  --
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
!           write (6,*) nshot,ibtcshot,mpnam2(i),namedum
           do j=1,np
             if (.not.oldcomp) then
             expmpi(j,i)=expmpi(j,i)-dumbtc*curtn1(j)
             do k=1,nccomp     !EJS(2014)
              expmpi(j,i)=expmpi(j,i)-dumccc(k)*curccoi(j,k)
             enddo
             endif
           enddo
           goto 33200
         endif
        enddo
!---------------------------------------------------- new section - EJS(2014)
! Compensate flux loops for coupling of individual C-coils.
! This was not needed previously, for antisymmetric (n=odd) C-coil pairs.
!-----------------------------------------------------------------------
        do i=1,nsilop
         if (lpname(i).eq.namedum) then
!           write (6,*) nshot,ibtcshot,lpname(i),namedum
           do j=1,np
             if (.not.oldcomp) then
             silopt(j,i)=silopt(j,i)-dumbtc*curtn1(j)
             do k=1,nccomp
              silopt(j,i)=silopt(j,i)-dumccc(k)*curccoi(j,k)
             enddo
             endif
           enddo
           goto 33200
         endif
        enddo
!---------------------------------------------- end of new section - EJS(2014)
        goto 33200
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
        goto 33000
      endif
34000 continue
      close(unit=60)
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
      open(unit=60,file=input_dir(1:lindir)//'icomp.dat_156456', &
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
!           write (6,*) nshot,ibtcshot,mpnam2(i),namedum
           oldexp=expmpi(1,i)
           do j=1,np
             do k=1,nicomp     !EJS(2014)
              expmpi(j,i)=expmpi(j,i)-dumcic(k)*curicoi(j,k)
             enddo
           enddo
           goto 35200
         endif
        enddo
!---------------------------------------------------- new section - EJS(2014)
! Compensate flux loops for coupling of individual C-coils.
! This was not needed previously, for antisymmetric (n=odd) C-coil pairs.
!-----------------------------------------------------------------------
        do i=1,nsilop   
         if (lpname(i).eq.namedum) then
!           write (6,*) nshot,ibtcshot,lpname(i),namedum
           do j=1,np
             do k=1,nicomp
              silopt(j,i)=silopt(j,i)-dumcic(k)*curicoi(j,k)
             enddo
           enddo
           goto 35200
         endif
        enddo
!---------------------------------------------- end of new section - EJS(2014)
        goto 35200
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
        goto 35000
      endif
36000 continue
      close(unit=60)
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
!--  correct sign of toroidal magnetic field to be consistent with   --
!--  the actual sign consistent with a right-handed cylindrical      --
!--  coordinate system                01/09/86                       --
!----------------------------------------------------------------------
      do 70 i=1,np
        bcentr(i)=bcentr(i)*tmu/rcentr*144.
   70 continue
      do 80 i=1,nfcoil
        do 75 j=1,np
   75   fccurt(j,i)=0.
        sclmp=1.0
        call avdata(nshot,fcname(i),i1,ierfc(i),fccurt(1,i), &
                    np,times,delt,i0,sclmp,i1,bitfc(i),iavem,time,ircfact, &
                    do_spline_fit,fc_rc(i),fcrcg(i),vresfc(i),fc_k(i),t0fc(i), &
                    devfc(1,i),navfc(1,i),time_err)
        do 77 j=1,np
        fccurt(j,i)=fccurt(j,i)*turnfc(i)
        devfc(j,i) =devfc(j,i)*turnfc(i)
   77 continue
      bitfc(i) =bitfc(i)*turnfc(i)
      fc_k(i)  =fc_k(i)*turnfc(i)
   80 continue
      rnavfc = navfc
!----------------------------------------------------------------
!--  New E-coil connection after discharge 85700, LLao 95/07/11--
!----------------------------------------------------------------
      do 90 i=1,nesum
        if (nshot.le.85700.and.i.gt.2) goto 83
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
            write(*,*) 'ERROR: ECOIL data unavailable for following times'
            do j=1,nvtime-1
              write(*,'(F8.1)') vtime(j)
            enddo
          endif
          iierr = 1
          return
        endif
   83   continue
   90 continue
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
!--  uncompensated diamagnetic flux if compensated not available     --
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
!--  get beam power                                                    --
!------------------------------------------------------------------------
      if (nshot.ge.53427) then
      call apdata(nshot,nsingl(6),i1,ierbim,pbinj(1), &
                  np,times,delt,i0,r1,i1,bitbim,iavem,time, &
                  do_spline_fit,beam_rc,beamrcg,vresbeam,beam_k, &
                  t0beam,devbeam(1),navbeam(1))
      do 150 i=1,np
        if (ierbim.ne.0) then
          pbinj(i)=0.0
        else
          pbinj(i)=pbinj(i)*1.e+03
        endif
  150 continue
      endif
!
      return
      end
      subroutine avdata(nshot,name,mmm,ierror,y, &
                        np,timesd,deltd,mm,xxd,nn,bitvld,kave,time,ircfact, &
                        do_spline_fit,rcx,rcgx,vbitx,zinhnox,t0x,stdevx,navx, &
                        ktime_err)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          avdata gets the data and optionally performs the        **
!**          average.                                                **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          09/06/86..........first created                         **
!**          93/04/21..........revised for double precision option   **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use eparmdud129,only:ntime
      use vtime_mod
      parameter (ntims=4096)
      common/gggttt/w(ntims),xw(ntims),bw(ntims),cw(ntims),dw(ntims) &
                    ,ew(ntims),stdevxx(ntims),navxx(ntims)
      common/inaver/iavem,iaved,iavev,iaveus
      character*10 name
      integer, intent(out) :: ktime_err
      real*8 y(1),time(ntime),deltd,xxd,bitvld,timesd &
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
      if (iaveus.gt.0) goto 1000
!----------------------------------------------------------------------
!-- milli-second averaging                                           --
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
!-- Check valid range of time-slice data                           --
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
         do 200 i=1,nnp
             timenow=time(i)
             ynow=sevals(npn,timenow,xw,w,bw,cw,dw)
             y(i)=ynow
  200    continue

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
 1000 continue
!--------------------------------------------------------------------------
!--  averaging in micro-seconds                                          --
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
!--     Check valid range of time-slice data                           --
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
                call smoothit2(xw,w,npn,dtave,stdevxx,navx)
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
               if(delta .lt. delta_min) then
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
      subroutine amdata(nshot,name,mmm,ierror,y, &
               np,timesd,deltd,mm,xxd,nn,bitvld,kave,time,ircfact, &
               do_spline_fit)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          amdata gets the data and optionally performs the        **
!**          average from MDS+                                       **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**        2006/06/14..........first created                         **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      include 'mdslib.inc'
      character*10 name, MyTree
      real*8 y(1),time(1),deltd,xxd,bitvld,timesd
      real*4, allocatable :: yw(:),xw(:),bw(:),cw(:),dw(:),ew(:)
      real*4 dtmin, xm5, dtave, xx, delt, times, delta_min, delta, &
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
!-- Get data from MDS+                                               --
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
      dsc = descr(IDTYPE_LONG, mylen, 0)
      status = Mdsvalue('SIZE('//name(1:lenname)//')'//char(0), &
                        dsc, 0, 1)
      if (mod(status,2) .eq. 0) then
        ierror = 1
        return
      endif


!----------------------------------------------
! 20140905 tbt Getting with pgf90 14.6:
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
      f_dsc = descr(IDTYPE_FLOAT, yw, mylen, 0)
      t_dsc = descr(IDTYPE_FLOAT, xw, mylen, 0)
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
      if (mylen.ge.2) go to 100
      ierror = 1
      return
!
  100 continue
!------------------------------------------------------------------------
!-- Check valid range of time-slice data                           --
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
        do 200 i=1,nnp
          timenow=time(i)
          ynow=sevals(npn,timenow,xw,yw,bw,cw,dw)
          y(i)=ynow
200     continue
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
      subroutine apdata(nshot,name,mmm,ierror,y, &
        np,timesd,deltd,mm,xxd,nn,bitvld,kave,time, &
        do_spline_fit,rcx,rcgx,vbitx,zinhnox,t0x,stdevx,navx)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          apdata gets the data and optionally performs the        **
!**          average.                                                **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          09/06/86..........first created                         **
!**          93/04/21..........revised for double precision option   **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      parameter (ntims=4096)
      common/gggttt/w(ntims),xw(ntims),bw(ntims),cw(ntims),dw(ntims) &
                    ,ew(ntims),stdevxx(ntims),navxx(ntims)
      real*8 y(1),time(1),deltd,xxd,bitvld,timesd &
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
  100   continue
  !------------------------------------------------------------------------
  !--     Check valid range of time-slice data                           --
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
        do 200 i=1,nnp
          timenow=time(i)
          ynow=sevals(npn,timenow,xw,w,bw,cw,dw)
          y(i)=ynow
200     continue
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
      subroutine gettanh(ishot,fitzts,ktime,time,ztssym,ztswid, &
                          ptssym,ztserr)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          gettanh gets the edge hyperbolic tangent fit parameters **
!**          from MDS+                                               **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          2001/03/09........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      parameter (ntims=4096)
      common/gggttt/w(ntims),xw(ntims),bw(ntims),cw(ntims),dw(ntims) &
                    ,ew(ntims),stdevxx(ntims),navxx(ntims)
      real*8 y(1),time(1),ztssym(1),ztswid(1),ptssym(1)
      character*2 fitzts,ztsfit
      logical ztserr(1)
      logical errzts(ntims)
      integer*4 iishot,kktime
!-----------------------------------------------------------------------
!-- Get edge pedestal tanh paramters                                  --
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
      subroutine avdiam(nshot,name,mmm,ierror,y, &
                    np,timesd,deltd,mm,xxd,nn,bitvl,kave,time,sigmay, &
                        ierdia)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          avdata gets the compensated diamagnetic data and        **
!**          and optionally performs the average.                    **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          09/06/86..........first created                         **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      parameter (ntims=4096)
      common/gggttt/w(ntims),xw(ntims),bw(ntims),cw(ntims),dw(ntims) &
                    ,ew(ntims),stdevxx(ntims),navxx(ntims)
      dimension ierdia(1)
      real*8 y(1),time(1),sigmay(1),deltd,xxd,timesd,bitvl
      character*10 name
      data dtmin/0.001001/,xm5/0.00001/
      save dtmin,xm5
!
      delt=deltd
      ierror=0
      dtmin=0.001001
      dtmin= min (dtmin,delt)
      dtmin= max (dtmin,xm5)
      mave=iabs(kave)
      npn = ntims
      tavg = 1.0
      call getdia(nshot,xw,npn,tavg,ierdia,w,ew)
      if (ierdia(2).gt.0.and.ierdia(3).gt.0) then
        ierror=1
        return
      endif
!------------------------------------------------------------------
!-- average data over mave ms                                    --
!------------------------------------------------------------------
      if (mave.ne.0) then
        dtave=mave*dtmin*2.
        call smoothit(xw,w,npn,dtave)
      endif
!
      call zplines(npn,xw,w,bw,cw,dw)
      do 200 i=1,np
          timenow=time(i)
          ynow=sevals(npn,timenow,xw,w,bw,cw,dw)
          y(i)=ynow
  200 continue
      call zplines(npn,xw,ew,bw,cw,dw)
      do 220 i=1,np
          timenow=time(i)
          ynow   =sevals(npn,timenow,xw,ew,bw,cw,dw)
          sigmay(i)=ynow
  220 continue
      return
      end
      subroutine zmooth(y,npts,nave)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          zmooth smooths out the data.                            **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          10/06/86..........first created                         **
!**                                                                  **
!**********************************************************************
      parameter (ntims=4096)
      common/gwork1/yave(ntims)
      dimension y(1)
!
      if (nave.eq.0) return
      if (nave.lt.0) go to 200
!---------------------------------------------------------------------
!-- uniform averaging                                               --
!---------------------------------------------------------------------
      imax=npts-nave
      ym=y(1)
      do 10 m=1,npts
        jb=m-nave
        if (m.le.nave) jb=m
        je=m+nave
        if (m.gt.imax) je=m
        nt=0
        ynew=0.0
        do 5 j=jb,je
        nt=nt+1
        ynew=ynew+y(j)
    5   continue
        yave(m)=ynew/nt
   10 continue
      do 20 m=1,npts
        y(m)=yave(m)
   20 continue
      return
  200 continue
!----------------------------------------------------------------------
!--  non-uniform weighted average                                    --
!----------------------------------------------------------------------
      nnave=iabs(nave)
      if (npts.lt.2*nnave+1) return
      do 320 m=1,nnave
      yave(m)=0.0
      sum=0.0
      do 300 i=1,nnave+1
        sumx=nnave-i+2
        sum=sum+sumx
        yave(m)=yave(m)+y(i+m-1)*sumx
  300 continue
      yave(m)=yave(m)/sum
  320 continue
      do 400 m=nnave+1,npts-nnave
        sum=0.0
        yave(m)=0.0
        do 330 i=1,nnave+1
          sumx=nnave-i+2
          sum=sum+sumx
          yave(m)=yave(m)+y(i+m-1)*sumx
  330   continue
        do 340 i=1,nnave
          sumx=nnave-i+1
          sum=sum+sumx
          yave(m)=yave(m)+y(m-i)*sumx
  340   continue
        yave(m)=yave(m)/sum
  400 continue
      do 500 m=npts-nnave+1,npts
        sum=0.0
        yave(m)=0.0
        do 440 i=1,nnave+1
          sumx=nnave-i+2
          sum=sum+sumx
          yave(m)=yave(m)+y(m-i+1)*sumx
  440   continue
        yave(m)=yave(m)/sum
  500 continue
      do 600 m=1,npts
        y(m)=yave(m)
  600 continue
      return
      end
      !
      subroutine smoothit(times,data,nts,timint)
        parameter (npmax=16384)
        dimension work(npmax)
        dimension times(2),data(2)
        !
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

      subroutine smoothit2(times,data,nts,timint,stdev,nave)
        parameter (npmax=16384)
        dimension work(npmax)
        dimension times(2),data(2),stdev(1),nave(1)
        !
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

      subroutine zplines(n, x, y, b, c, d)
      integer n
      real x(n), y(n), b(n), c(n), d(n)
!
!  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
!  for a cubic interpolating spline
!
!    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!
!    for  x(i) .le. x .le. x(i+1)
!
!  input..
!
!    n = the number of data points or knots (n.ge.2)
!    x = the abscissas of the knots in strictly increasing order
!    y = the ordinates of the knots
!
!  output..
!
!    b, c, d  = arrays of spline coefficients as defined above.
!
!  using  p  to denote differentiation,
!
!    y(i) = s(x(i))
!    b(i) = sp(x(i))
!    c(i) = spp(x(i))/2
!    d(i) = sppp(x(i))/6  (derivative from the right)
!
!cccccccccccccc
!  the accompanying function subprogram  seval  can be used
!  to evaluate the spline.
!
!
      integer nm1, ib, i
      real t
!
      nm1 = n-1
      if ( n .lt. 2 ) return
      if ( n .lt. 3 ) go to 50
!
!  set up tridiagonal system
!
!  b = diagonal, d = offdiagonal, c = right hand side.
!
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do 10 i = 2, nm1
         d(i) = x(i+1) - x(i)
         b(i) = 2.*(d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i) = c(i+1) - c(i)
   10 continue
!
!  end conditions.  third derivatives at  x(1)  and  x(n)
!  obtained from divided differences
!
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.
      c(n) = 0.
      if ( n .eq. 3 ) go to 15
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
!
!  forward elimination
!
   15 do 20 i = 2, n
         t = d(i-1)/b(i-1)
         b(i) = b(i) - t*d(i-1)
         c(i) = c(i) - t*c(i-1)
   20 continue
!
!  back substitution
!
      c(n) = c(n)/b(n)
      do 30 ib = 1, nm1
         i = n-ib
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
   30 continue
!
!  c(i) is now the sigma(i) of the text
!
!  compute polynomial coefficients
!
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
      do 40 i = 1, nm1
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.*c(i)
   40 continue
      c(n) = 3.*c(n)
      d(n) = d(n-1)
      return
!
   50 b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.
      d(1) = 0.
      b(2) = b(1)
      c(2) = 0.
      d(2) = 0.
      return
      end
      real function sevals(n, u, x, y, b, c, d)
      integer n
      real  u, x(n), y(n), b(n), c(n), d(n)
!
!cccccccccccccc
!cccccccccccccc
!  this subroutine evaluates the cubic spline function
!
!    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
!
!    where  x(i) .lt. u .lt. x(i+1), using horner's rule
!
!  if  u .lt. x(1) then  i = 1  is used.
!  if  u .ge. x(n) then  i = n  is used.
!
!  input..
!
!    n = the number of data points
!    u = the abscissa at which the spline is to be evaluated
!    x,y = the arrays of data abscissas and ordinates
!    b,c,d = arrays of spline coefficients computed by spline
!
!  if  u  is not in the same interval as the previous call, then a
!  binary search is performed to determine the proper interval.
!
      integer i, j, k
      real dx
      data i/1/
      if ( i .ge. n ) i = 1
      if ( u .lt. x(i) ) go to 10
      if ( u .le. x(i+1) ) go to 30
!
!  binary search
!
   10 i = 1
      j = n+1
   20 k = (i+j)/2
      if ( u .lt. x(k) ) j = k
      if ( u .ge. x(k) ) i = k
      if ( j .gt. i+1 ) go to 20
!
!  evaluate spline
!
   30 dx = u - x(i)
      sevals= y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
      return
      end
!
!**********************************************************************
! 
!  PRIMARY CALLING SEQUENCE USED BY REVIEW AND OTHER CODES.
!
      SUBROUTINE GETDAT(NSHOT,NAME,ICAL,IER,T,DATA,NP, &
      TMIN,TMAX,MOP,SCALE,jWAIT)

!----- 
! 1% tolerance for signal clipping added 8/14/89
! at present the clipping errors (ier = -1, -2) are used only by the
! yoka routine in NEWSPEC.
!-----
!
!     GETDAT CALLS PTDATA AND RETURNS SHOT DATA IN A MORE USEFUL WAY.
!
!     NSHOT IS THE SHOT NUMBER
!
!     NAME IS THE POINT NAME(10 ASCII CHARACTERS)
!
!     ICAL=0 RETURNS DATA IN DIGITIZER COUNTS
!          1 RETURNS CALIBRATED DATA IN PHYSICAL UNITS
!          2 RETURNS THE VOLTAGE INPUT TO THE DIGITIZER
!          3 RETURNS THE VOLTAGE AT THE INTEGRATOR OUTPUT
!          4 RETURNS THE INTEGRATED SIGNAL IN V-SEC
!          10-19 ARE THE SAME AS 0-9 EXCEPT FOR BASELINE ALGORITHM
!
!     ICAL=0 HAS NO BASELINE SUBTRACTION
!           1-9 USES BASELINE FROM PTDATA (OBTAINED FROM EARLY SAMPLES)
!           10-19 USE DIGITIZER MIDPOINT AS ZERO
!
!     IER IS AN ERROR RETURN FLAG.
!     IER = POSITIVE NUMBER IS THE ERROR RETURN FROM PTDATA.
!     IER = -1 IS FOR OVERFLOW OF THE DATA IN THE DIGITIZER
!     IER = -2 IS FOR UNDERFLOW OF THE DATA
!     IER = -3 IS FOR BASELINE(ZERO) OUTSIDE THE RANGE OF DIGITIZER
!     IER = -6 DATA ARE RETURNED WITH SOME LOSS OF TIME RESOLUTION 
!  (TOO MANY SAMPLES FOR BUFFER SIZE)
!     IER = -7 DATA MAY HAVE SOME LOSS OF TIME RESOLUTION 
!  (TOO MANY SAMPLES FOR BUFFER SIZE)
!  STATUS UNKNOWN DUE TO NECESSITY OF CALLINBG PTDATA BY TIME
!
!     T IS THE ARRAY IN WHICH THE TIMES OF THE DATA SAMPLES WILL BE RETURNED
!
!     DATA IS THE ARRAY IN WHICH THE DATA WILL BE RETURNED
!
!     NP IS THE MAXIMUM NUMBER OF DATA POINTS YOU WANT BACK
!       ACTUAL NUMBER IS RETURNED
!
!     TMIN, TMAX (IN SECONDS) DEFINE THE TIME INTERVAL ON WHICH YOU WANT
!       DATA.  ACTUAL VALUES ARE RETURNED.
!
!     THE POINTNAME DATA IS ALWAYS MULTIPLIED BY SCALE BEFORE USE OF
!
!     MOP--WHICH IS THE OP CODE.
!     IF MOP=0 THE POINT IS WRITTEN INTO DATA, OVERWRITING THE
!        PREVIOUS CONTENTS.
!     IF MOP=1,THE POINT IS ADDED INTO THE ARRAY DATA.
!     IF MOP=2,THE POINT IS SUBTRACTED FROM THE CURRENT CONTENTS OF DATA.
!     IF MOP=3 MULTIPLY CURRENT CONTENTS OF DATA BY POINTNAME
!     IF MOP=4 DIVIDE CURRENT CONTENTS OF DATA BY POINTNAME
!
!     IBFLAG = 0 PTDATA'S "ZERO OFFSET" VALUE IS USED FOR THE BASELINE.
!              1 SAMPLES 6 THROUGH 35 ARE AVERAGED TO PROVIDE THE BASELINE.
!                    IN THIS CASE, THE FIRST 40 SAMPLES ARE NOT RETURNED.
!
!

      parameter (npmax=262144)

      character*4 char4

! 20140905 tbt Dimension of 1 bombs when array bounds checked!
!     DIMENSION  DATA(1)
!     DIMENSION  FPDATA(1)

      DIMENSION  DATA(*)
      DIMENSION  FPDATA(npmax)

      dimension  int16(2),int32(2),iarr(50)
      dimension  iascii(12)
      CHARACTER*3  MONTH(12)
      CHARACTER  NAME*(*)
      character*12  pcst
      CHARACTER*10 PCSTIME
      CHARACTER  PHASE*4
      dimension  rarr(20+npmax)
      dimension  real32(100)
      REAL*8  REAL64(100)
      CHARACTER*10  SDATE, STIME

! 20140905 tbt Dimension of 1 bombs when array bounds checked!
!     DIMENSION  T(1)
      DIMENSION  T(*)

      REAL*8  TDUM(2+NPMAX)
      INTEGER*4  TEMP(npmax)
      REAL*8  TIME64(NPMAX)

! 20140905 tbt Dimension of 1 bombs when array bounds checked!
!     DIMENSION  TT(1)
      DIMENSION  TT(npmax)

      COMMON /TRAPPER/ IT_CHAN, ITRAPA

      equivalence (ichar,char4)
      equivalence (pcst,pcstime)
      EQUIVALENCE (TT,RARR(21)), (FPDATA(1),TEMP(1))

      DATA MONTH/'JAN','FEB','MAR','APR','MAY','JUN', &
         'JUL','AUG','SEP','OCT','NOV','DEC' /
      DATA PHASE /'.PLA'/
      data tolclip /.01/ !-- maximum fraction of clipped signals
      data kmax/10/
      data wait_time/30./

      LOGICAL INTCAL,efit
      character*10 cname, cname1, filnam
!----------------------------------------------------------------------

      INTCAL = .FALSE.
      efit = .false.
      rcfact = 1.0
      iwait = jwait

      GO TO 11111

!---------------------------------------------------------------
!
!  ALTERNATE ENTRY POINT FOR USE BY INTEGRATOR CALIBRATION CODE INTCAL.
!  THE ONLY DIFFERENCE IS THE RETURN OF ADDITIONAL PARAMETERS.
!
!---------------------------------------------------------------
      ENTRY GETDAT_I(NSHOT,NAME,ICAL,IER,T,DATA,NP, &
        TMIN,TMAX,MOP,SCALE,jWAIT, &
        RC_I,RCG_I,VPBIT_I,MZERO_I,NBITS_I,SDATE,STIME)

      INTCAL = .TRUE.
      efit = .false.
      rcfact = 1.0
      iwait = jwait

      GO TO 11111

!----------------------------------------------------------------------
!
! Alternate entry point for use by the EFIT code.
!
!----------------------------------------------------------------------

      ENTRY GETDAT_E(NSHOT,NAME,ICAL,IER,T,DATA,NP, &
      TMIN,TMAX,MOP,SCALE,bitone,ircfact,RC_E,RCG_E,VPBIT_E, &
       ZINHNO_E,T0_E)

      INTCAL = .FALSE.
      efit = .true.
      iwait = 0

        rcfact = 1.0

        IF (ircfact .eq. 1) THEN
                filnam='rcfact.dat'
                cname=name
                call str$upcase(cname,cname)
                open(unit=21,file=filnam,status='old',err=1000)
1001            read (21,9,end=1000,err=1000) cname1, fact
9               format(x,a10,f15.5)
                if(cname .ne. cname1)   go to 1001
                if(fact .gt. 0.) rcfact = fact
1000            close(unit=21)
        ENDIF

!----------------------------------------------------------------------
11111 kount=0

1      IASCII(1)=5
      INT16(1)=0
      INT32(1)=0
      REAL32(1)=50

! ---------------------- SET UP FOR THE PTDATA CALL ----------------------

      itype = 12  ! point type call, variable timing
      iarr(1) = npmax  ! number of points requested
      iarr(3) = 1  ! start point
      nint = 1
      iarr(4) = nint  ! point interval
      pzero = 0.0
!
!  PTDATA CALL ... GET ALL THE DATA!!
!  Request to get data by point sampling with PTDATA returning time array
!
!vas f90 modifi
!10 call ptdata(itype,nshot,%ref(phase),%ref(name),temp,ier,iarr,rarr,
!vas 1 iascii,int16,int32,real32)
10    call ptdata(itype,nshot,%ref(phase),%ref(name),temp,ier,iarr,rarr,&
        iascii,int16,int32,real32)

! write(13,91300)(rarr(klg),klg=540,640)
91300 format(x,1pe17.10)

      if (ier .eq. 4 .or. ier.eq.2) ier = 0
      VPBIT  = real32(51)
      RC     = real32(52)

!
!  CHECK DATA FORMAT, TRY AGAIN IF NECESSARY
!  Attempt to get data based on time sampling, using a computed
!  start and delta time.
!

      IF (ier .ne. 35) go to 30 ! ier=35 ... wrong call type for dfi

      itype = 2
!vas f90 modifi
!vas20 call ptdata(itype,nshot,%ref(phase),%ref(name),temp,ier,iarr,rarr,
!vas 1 iascii,int16,int32,real32)
20    call ptdata(itype,nshot,%ref(phase),%ref(name),temp,ier,iarr,rarr,&
        iascii,int16,int32,real32)
      if (ier .eq. 4 .or. ier.eq.2) ier = 0

      VPBIT  = real32(13)
      RC     = 1.0

30    continue

!
!  WAIT FOR THE DATA TO COME IN
!

!vas IF (iwait .eq. 1 .and. (ier .eq. 3 .or. ier .eq. 1)
!vas 1   .and. kount .lt. kmax .and. itrapa .eq. 0) THEN
      IF (iwait .eq. 1 .and. (ier .eq. 3 .or. ier .eq. 1) &
         .and. kount .lt. kmax .and. itrapa .eq. 0) THEN
! MPI >>>
        !call rev_wait(iwait, nshot)
! MPI <<<
        if (itrapa .eq. 0) go to 1
      ENDIF

      !-- AT THIS POINT ANY REMAINING POSITIVE ERROR VALUE IS FATAL.
      !-- NON-FATAL POSITIVE ERRORS HAVE BEEN HANDLED AS FOLLOWS:
      !  ier = 4  reset ier to 0
      !  ier = 2  non-fatal warning:  pointname has been revised
      !  ier = 1 or 3 looped back to start if iwait=1
      !  ier = 35  re-called ptdata w/o time array for this dfi

      IF (IER .gt. 0 .AND. IER.NE.2) THEN
70000   np = 2
        t(1)=tmin
        t(2)=tmax
        data(1)=0.
        data(2)=0.
        RETURN
      ENDIF

      !
      !  At this point, data is available.  Save some required information.
      !

      idfi = iarr(31)
      nret = iarr(2)
      nsampl = iarr(32)   ! number of 16-bit words
      nbits = iarr(33)   ! number of bits per word
      ZINHNO = rarr(4)                        ! inherent number
      RCG    = rarr(5)
      if(nbits .le.  8) nsampl = nsampl*2   ! correct number of
      if(nbits .gt. 16) nsampl = nsampl/2 !       available samples
      mzero = iarr(30)
      yzero = mzero
      IF (idfi.eq.158 .OR. IDFI.EQ.2158) THEN !-- NEW DFI WITH MORE ACCURATE OFFSET
        yzero = real32(14)
        mzero = anint(yzero)
      ENDIF

!--- WAVEFORMS HAVE AN ADDITIONAL OFFSET FOR CONVERSION TO PHYSICAL UNITS
! --- this offset makes sense only with ical=1 or 11,
! --- by convention we will add offset in ical=1 but not 11
!vas f90 modifi
!vas if ((idfi .eq. 125 .or. idfi .eq. 158 .OR. IDFI.EQ.2158)
!vas 1 .and. ical .eq. 1)  pzero = real32(8)
      if ((idfi .eq. 125 .or. idfi .eq. 158 .OR. IDFI.EQ.2158) &
        .and. ical .eq. 1)  pzero = real32(8)

        !--- Rely on PTDATA's time-type call if the
        !    pointname is too big for the buffer.
        !--- This approach may return less than optimum resolution if the requested
        !    interval crosses a clock domain boundary.
        !
        !  Must handle data from new (4-95) Plasma Control System separately.
        !  Assume that PCS data will be less than max qty asked for already.
        !  (Apr 5, 1994:  John Ferron said some waveforms could have as much
        !  as 90K samples.  Max data requested in 1st PTDATA call 128K.)
        !
        !  If data is NOT PCS, then compute appropriate start and delta times
        !  in order to do 2nd PTDATA call.

        if (idfi.eq.2201 .or. idfi.eq.2202 .or. idfi.eq.2203) then
          if (idfi.eq.2201) then  ! digitizer data, int
            vpbit = real32(5)
            yzero = real32(6)
            if (real32(2).ge.5) rc = real32(7)
          else if (idfi.eq.2202) then  ! processed data, int
            vpbit = -1.0
            pzero = real32(3)
            yzero = 0.0
          else if (idfi.eq.2203) then  ! processed data, real
            vpbit = -1.0
            pzero = real32(3)
            yzero = 0.0
        end if
        ichar = iascii(3)
        pcst(1:4) = char4
        ichar = iascii(4)
        pcst(5:8) = char4
        ichar = iascii(5)
        pcst(9:12) = char4
        iarr(1) = npmax
        iarr(3) = 1
        iarr(4) = nint
        iascii(1) = 0
        int16(1) = 0
        int32(1) = 0
        real32(1) = 0
        real64(1) = 0
        call ptdata64(64,nshot,%ref(phase),%ref(pcstime),time64,ker, &
          iarr,rarr,iascii,int16,int32,real32,real64,tdum)
        if (ker.ne.0 .and. ker.ne.2 .and. ker.ne.4) then
          ier = 1
          go to 70000
        end if
        nrettim = iarr(2)
        do k = 1,nrettim
          tt(k) = time64(k)
        end do
        go to 85000
      else IF (itype .eq. 12 .or. itype .eq. 2) THEN
      !
      !  at this point, is there more digitizer data available?
      !
      IF (nret .lt. nsampl) THEN
        itype = itype - 1  ! call ptdata by time
        rarr(1) = tmin   ! start time
        rarr(2) = (tmax-tmin)/(npmax-1) ! minimum delta-time
        iarr(1) = npmax   ! guaranteed to cover tmax-tmin
        IF (itype .eq. 11) THEN
          go to 10
        ELSE
          go to 20
        ENDIF    ! itype = 11
      endif     ! nret < nsample
      endif     ! itype = 12

      !--- ier = -7 for time type call ... possible loss of time resolution

      if ((itype .eq. 1 .or. itype .eq. 11)  .and. ier .le. 0) ier = -7

      !  FILL IN TIME ARRAY IF NECESSARY

      IF (itype .eq. 2) THEN
        dt = rarr(9)
        t0 = rarr(8) - dt
        DO i=1, nret
          tt(i) = t0 + i * dt
        ENDDO
      ENDIF

85000 continue

      !  FIND START AND STOP TIME OF DATA TO BE RETURNED
      !  --- THE OVERLAP OF THE REQUESTED AND DIGITIZED TIME INTERVALS

      tmind = tt(1)
      tmaxd = tt(nret)

      IF (tmax .le. tmin) THEN
        tmin = tmind
        tmax = tmaxd
      ENDIF

      tmin = amax1(tmin,tmind)
      tmax = amin1(tmax,tmaxd)

      IF (tmax .lt. tmin) THEN
        np = 0
        return
      ENDIF

      !  FIND THE FIRST AND LAST DATA POINTS TO BE RETURNED
      ! Traditionally, the efit code has used a different limit test here.  So that
      ! efit will continue to get the same answer, one of two different limit tests
      ! is used here.
      !
      nfirst = 1
      nlast  = 0
      ITOTALPTS = NRET
      IF (NSAMPL.LT.ITOTALPTS) ITOTALPTS = NSAMPL

      if (efit) then
        DO i=1,ITOTALPTS
          ttime = tt(i)
          if (ttime .lt. tmin) nfirst=i+1
          if (ttime .lt. tmax) nlast =i
        ENDDO
      else
        DO i=1,ITOTALPTS
          ttime = tt(i)
          if (ttime .lt. tmin) nfirst=i+1
          if (ttime .le. tmax) nlast =i
        ENDDO
      endif

      ndata = nlast - nfirst + 1

      !  DECIDE ON DATA POINT INTERVAL

      if (np .eq. 0) then
        np = 16384
      else
        np = amin0(np,npmax)
      endif

      nint = (ndata-1) / np + 1
      ! if (nint .gt. 1 .and. ier .le. 0) ier = -6
      if (nint.gt.1) ier = -6

      ! SOME OF THE ICH DATA  (ICHPWR) IS
      ! JUST LIKE DIGITIZER DATA EXCEPT IN
      ! REAL FORMAT

      IF(IDFI .EQ. 119 .or. idfi.eq.2119) THEN      !ICH STUFF EXISTS IN FLOATING FORMAT
        II=0
        DO I = NFIRST, NLAST, NINT
          II=II+1
          T(II) = TT(I)
          DATA(II) = FPDATA(I)
        ENDDO
        NP = II
        RETURN
      ENDIF

      !  CHECK FOR ZERO OFFSET OUT OF RANGE
      !  2013/04/03 Revised for new digitzers signed binary numbers  LL/ES

      ! min = 0
      !
      !  THIS NEXT LINE HAS BEEN MODIFIED TO (NBITS - 2) FOR TWO REASONS:
      !  (1)  THE FIRST BIT IS ACTUALLY BIT #0 (THEREFORE FOR A 32-BIT VARIABLE
      !       THE LARGEST UNSIGNED VALUE THAT THE VARIABLE CAN CONTAIN IS
      !       2**(NBITS) - 1
      !  (2)  FOR A "SIGNED" DATA VALUE, THE LAST BIT IS USED TO INDICATE THE
      !       SIGN.  THEREFORE, THE LARGEST VALUE A SIGNED VARIABLE CAN CONTAIN
      !       IS     2**(NBITS - 1) - 1
      !  THE CALCULATION IN #2 ABOVE WILL RESULT IN AN INTEGER OVERFLOW BECAUSE
      !  THE MULTIPLICATION IS DONE BEFORE THE SUBTRACTION.  AND IT IS THE
      !  RESULT OF THE MULTIPLICATION THAT WILL EXCEED THE ABILITY OF THE
      !  INTEGER TO HOLD IT.
      !
      !  THE NEW LINE IS A MODIFIED FORM TO MAKE SURE THAT INTEGER OVERFLOW
      !  DOES NOT OCCUR
      !
      IF (NBITS.GE.16) THEN  !  NECESSARY ONLY FOR 16-BIT DATA
        min =-(2 ** (nbits - 2) - 1) * 2 + 1
        max = (2 ** (nbits - 2) - 1) * 2 + 1
      ELSE    !  12-BIT DIG.S NEED THIS ONE
        min = 0
        MAX = 2**NBITS - 1
      END IF


      !--- WAVEFORMS ARE ALLOWED TO HAVE ZERO OFFSET OUTSIDE THE GENERATOR RANGE
      !vas f90 modifi
      !vas IF ((mzero .lt. min .or. mzero .gt. max) .AND.
      !vas 1 (idfi .ne. 125 .and. idfi .ne. 158 .AND. &
      IF ((mzero .lt. min .or. mzero .gt. max) .AND. &
        (idfi .ne. 125 .and. idfi .ne. 158 .AND. &
        IDFI.NE.2158)) THEN
        IF (ier .eq. 0 ) ier = -3
      ENDIF

      !  CALIBRATION CODES TO RETURN ABSOLUTE SIGNAL LEVEL, NO BASELINE SUBTRACTION
      !  MUST ASSUME DIGITIZER ZERO IS IN THE MIDDLE OF ITS RANGE, SINCE THIS
      !  INFORMATION IS NOT STORED IN THE DATA BASE!
      !  WAVEFORMS (idfi = 125, 158, 2158) DO NOT HAVE A BASELINE, SINCE WE ARCHIVE
      !  ONLY THE PROGRAMMED VALUE AND NOT THE ACTUAL VALUE.
      ical0 = ical
      IF (ical .ge. 10) THEN
        if (idfi .ne. 125 .and. idfi .ne. 158 &
          .AND. IDFI.NE.2158) mzero = 2**(nbits-1)
        yzero = mzero
        ical0 = ical-10
      ENDIF

      !   FILL UP PLOT ARRAYS

      nclip = 0
      ii = 0

      DO i=nfirst, nlast, nint

        ii = ii + 1
        !
        ! --- time array
        !
        t(ii) = tt(i)
        ! write(14,91300)t(ii)
        !
        ! --- data calibration
        !     some waveforms for PCS have data that is REAL*4
        !
        if (idfi.eq.2203) then
          y = fpdata(i)
        else
          y = temp(i)
        endif

        ! --- check for data out of digitizer range
        !      IF (TEMP(I) .GE. MAX .AND. IER.EQ.0)  IER = -1
        !      IF (TEMP(I) .LE. MIN .AND. IER.EQ.0)  IER = -2
        if(y .ge. max .or. y .le. min) nclip = nclip + 1

        ! --- offset for all but calibration code 0
        IF (ical .ne. 0) THEN
          y = y - yzero
        !  y = y - yzero - deltab   ! ICAL > 0
        ENDIF

        ! --  scale factors for the different calibration codes

        IF (ical0 .eq. 0) THEN
          y = y    ! ICAL = 0
        ELSE IF (ical0 .eq. 2) THEN
          y = -y * vpbit   ! ICAL = 2
        ELSE IF (ical0 .eq. 3) THEN
          y = -y * rcg / rc  ! ICAL = 3
        ELSE IF (ical0 .eq. 4) THEN
          y = -y * rcg   ! ICAL = 4
        ELSE
          y = -y * rcg * zinhno + pzero ! ICAL = 1 or other
        ENDIF

        y = scale * y * rcfact    ! user's scale factor

        ! --- arithmetic combinations as determined by mop code
        IF (mop .eq. 0) THEN
          data(ii) = y
        ELSE IF (mop .eq. 1) THEN
          data(ii) = data(ii) + y
        ELSE IF (mop .eq. 2) THEN
          data(ii) = data(ii) - y
        ELSE IF (mop .eq. 3) THEN
          data(ii) = data(ii) * y
        ELSE IF (mop .eq. 4) THEN
          if (y.ne.0) then
            data(ii) = data(ii) / y
          else
            data(ii) = 0.
          end if
        ELSE
          data(ii) = y
        ENDIF

      ENDDO



      ! --- actual number of points returned
      np = ii

      ! --- clipping tolerance:
      !     set error flag only if clipped fraction is outside tolerance
      if(np .gt. 0) then
        fclip = float(nclip)/float(np)
        if (fclip .gt. tolclip .and.  &
          ier .le. 0) ier = -1
      endif
      IF (INTCAL) then
        !
        !   ADDITIONAL ARGUMENTS FOR INTEGRATOR CALIBRATION CALL.
        ! THESE DUMMY ARGUMENTS MUST ONLY BE REFERENCED WHEN USING THE
        ! ENTRY POINT GETDAT_I.  OTHERWISE ACCESS VIOLATION ERRORS OCCUR.
        !
        RC_I = RC
        RCG_I = RCG
        VPBIT_I = VPBIT
        MZERO_I = MZERO
        NBITS_I = NBITS
        WRITE(SDATE,880) IARR(19),MONTH(IARR(18)), &
          IARR(20)-IARR(20)/100*100
        WRITE(STIME,881) IARR(15),IARR(16),IARR(17)
880     FORMAT(I2,'-',A3,'-',I2)
881     FORMAT(I2,':',I2,':',I2)
      !
      endif

      if(efit) then
        if (ical0 .eq. 0) then
          bitone = 1                              ! ical = 0
        else if (ical0 .eq. 2) then
          bitone = vpbit                          ! ical = 2
        else if (ical0 .eq. 3) then
          bitone = rcg / rc                       ! ical = 3
        else
          bitone = rcg * zinhno                   ! ical = 1 or other
        endif

        bitone = scale * abs(bitone)
        !--------------------------------------------------------------------
        !-- New returned values for new error matrix                       --
        !--------------------------------------------------------------------
        RC_E = RC
        RCG_E = RCG
        VPBIT_E = VPBIT
        ZINHNO_E = ZINHNO
        T0_E = TT(1)
      endif

      RETURN
      END

      subroutine magsigma(ishotx,timexy,jtimex,gradsmpx,gradsflx, &
                        bpermpx,sigmafx,sigmabx,sigmaex, &
                        sigmaipx,sigmaflx,sigmampx)
!*************************************************************************
!**
!**     MAIN PROGRAM: MHD FITTING CODE
!**
!**
!**     SUBPROGRAM DESCRIPTION:
!**     Calculates the experimental magnetic uncertainties
!**
!**     CALLING ARGUMENTS:
!**     ishot  = shot number
!**     time  = time slice
!**     gradsmp = Grad(S) of magnetic probe
!**     s   = BR cost + BZ sint
!**     grads_R = dBRdR cost + dBZdR sint (d/dR partial wrt R)
!**     grads_Z = dBRdZ cost + dBZdZ sint
!**     gradsmp  = sqrt (grads_R**2 + grads_Z**2)
!**     gradsfl  = Grad(S) of flux loop
!**     bpermp  = B perpendicular to the magnetic probe
!**     sigmaf  = uncertainty of f coil
!**     sigmab  = uncertainty of b coil
!**     sigmae  = uncertainty of e coil
!**     sigmafl  = uncert. of flux loop
!**     sigmamp  = uncert. of magnetic probes
!**     REFERENCES:
!**     (1)
!**     (2)
!**
!**     RECORD OF MODIFICATION:
!**     2003. ....first created  E.J. Strait
!**
!*************************************************************************
!**     This subroutine calculates the uncertainties for the magnetic 
!**     diagnostics.  It is based on estimates described in
!**     DIII-D Physics Memo D3DPM 0202, "Estimating the Uncertainty of DIII-D 
!**     Magnetic Data," by E.J. Strait (Aug. 30, 2002).
!**     
!**     The following sources of uncertainty are included. (Further explanation 
!**     is given under the corresponding item numbers in the Physics Memo.)
!**     The individual terms are combined in quadrature to give
!**     the total uncertainty for each signal.
!**     
!**     1) Loop calibration   dS = a1 S
!**     2) Loop cal. - Long-term change  dS = a2 S
!**     3) Integrator calibration  dS = a3 S
!**     4) Int. cal. - Long-term change  dS = a4 S
!**     5) Integrator drift    dS = a5 K(RC/G) T
!**     6) Loop position    dS = a6 grad(S)
!**     7) Loop tilt angle   dS = a7 Bperp
!**     8) Bt pickup     dS = a8 dBt
!**     9) C-coil pickup   dS = a9 Cn
!**     10) Bp pickup in leads   dS = a10 K integral(Bperp^2)ds
!**     11) Bp pickup in ports   dS = a11 K Bport
!**     12) Detector nonlinearity  dS = a12 S/Bt
!**     13) Noise     dS = a13 (/<S^2> - <S>^2/)^0.5 / N^0.5
!**     14) Digitizer resolution  dS = a14 K(RC/G)(Vres/2) / N^0.5
!**             where
!**     S  = measured signal: flux, field, or current (in physics units)
!**     dS  = estimated uncertainty in the measured signal
!**     grad(S)  = gradient of measured flux or field (physics units/meter)
!**     N   = number of time samples averaged
!**     Vres  = one-bit resolution of the digitizer (volts)
!**     K   = inherent number (physics units/volt-second)
!**     RC  = integrator time constant (seconds)
!**     G   = integrator gain
!**     T   = time elapsed since the start of integration (seconds)
!**     dBt  = change in toroidal field since start of integration (seconds)
!**     Bperp  = poloidal field normal to axis of mag. probe or leads (Tesla)
!**     Cn  = current in a C-coil pair (Amps)
!**     integral(f)ds = integral along leads from probe to connector (meters)
!**     Bport  = poloidal field in port, perpendicular to leads (Tesla)
!**     an  = numerical coefficient: units vary with n,
!**                       and values vary between types of signals
!**     
!**     Note that items 10 and 11 will not be implemented immediately due to
!**     the additional difficulty in calculating the path integral (10) and in 
!**     estimating Bport outside the efit grid (11).
!********************************************************************** 
!**     There are several classes of magnetic diagnostics.      
!**     They are represented in variable names by the following letters:
!**             mp  = magnetic probes
!**             fl  = flux loops
!**             f  = F-coil Rogowski loops
!**             b  = B-coil Rogowski loops
!**             e  = E-coil Rogowski loops
!**             ip  = plasma current Rogowski loops
!**
!********************************************************************** 


      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
!  include 'ecomdu1.f90'
!  include 'ecomdu2.f90'
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
      dimension sigmp(magpri)
      dimension sigfl(nsilop)
      dimension sigfc(nfcoil)
      dimension sige(nesum)
!-----
! MASKS FOR SUBSETS WITHIN ARRAYS
!-----
!-- poloidal array probes
      dimension maskpol(magpri)
        data maskpol/ 60*1., 16*0./
!-- RDP probes
      dimension maskrdp(magpri)
        data maskrdp/ 60*0., 16*1./
!-- F-coil flux loops
      dimension maskff(nsilop)
      data maskff/18*1.,7*0.,1., 0., 1.,7*0.,1.,0.,1.,6*0./
!-- inner F-coil flux loops
      dimension maskinf(nsilop)
      data maskinf/ 5*1., 2*0., 7*1., 2*0., 2*1., 26*0. /
!-- outer F-coil  loops
      dimension maskoutf(nsilop)
      data maskoutf/5*0., 2*1., 7*0., 2*1., 2*0., &
                     7*0., 1., 0., 1., 7*0., 1., 0., 1., 6*0./
!-- vacuum vessel flux loops
      dimension maskvv(nsilop)
        data maskvv/18*0.,7*1.,0., 1., 0., 7*1., 0., 1.,0.,6*1./
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
1       dd = (/.0019, .0019, 0., 0., 0., .0017, .0017, .003/)
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
2       dd = (/.002, .002, 0., 0., 0., .003, .003, .006/)
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
3       dd = (/.0013, .0013, .0013, .0013, .0013, .0013, .0013, .0013/)
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
4       dd = (/.0017, .0017, .0017, .0017, .0017, .0017, .0017, .0017/)
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
5       dd = (/.0007, .0007, .0007, .0007, .0007, .0022, .0007, .0007/)
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
6       dd = (/.0020, .0020, .0030, .0045, .0020, 0., 0., 0./)
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
7       dd = (/.017, .017, 0., 0., 0., 0., 0., 0./)
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
8       dd = (/.003, .003, .00044, .00044, .00044, 0., 0., 1.3e4/)
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
9       dd  = (/2.3e-08,  1.4e-08,  5.1e-08,  5.1e-08, &
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
10      dd = (/.00016,  .00016,  0.,  0., .00016, 0., 0., 0./)
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
11      dd = (/.0002,  .0002,  0., 0., .0002, 0., 0., 0./)
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
12      dd = (/.0025, .0025, 0., 0., 0., 0., 0., 0./)
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
13      dd = (/1., 1., 1., 1., 1., 1., 1., 1./)
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
14      dd = (/1., 1., 1., 1., 1., 1., 1., 1./)
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


100   sigmamp(jtimex,:) = sqrt(sigmp)
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

99    format(5e12.4)

      return
    end
    !****************************************************
    ! subroutine rev_wait() added by Qilong Ren
    !************************************************
    subroutine rev_wait(iwait, nshot)
      character*80 line
      data kmax/20/, wait_time/30./
      common/waitkount/k
      call shotno(line, lastsh, ier)
      !vas f90 modifi
      if((nshot .ne. lastsh .and. nshot .ne. lastsh+1) &
        .or. ier .ne. 0) then
        iwait = 0
        return
      endif
      call lib$wait(wait_time)
      k = k+1
      if(k .gt. kmax) iwait = 0
      return

      entry init_wait
      k = 1
      return
    end

! =========================================================

! MPI >>>
#if defined(USEMPI)

      subroutine getpts_mpi(nshot,times,delt,ktime,istop)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          getpts_mpi...                                           **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**       2020/09/18 ....... R.S. Bug fix, bcast oldccomp.           **
!**                          This removed small differences between  **
!**                          serial and parallel runs.               **
!**                                                                  **
!**********************************************************************
        include 'eparmdud129.f90'
        include 'modules1.f90'
        
        include 'mpif.h'
        
        ! Number columns ZWORK 2-D array
        !   EXPMPI 1:magpri
        !   SILOPT 1:nsilop
        !   FCCURT 1:nfcoil
        !   DENVT  1:nco2v
        !   DENRT  1:nco2r
        !   ECCURT 1:nesum
        parameter (nsize=18+magpri+nsilop+nfcoil+nco2v+nco2r+nesum)
        ! Dimension ZWORK2 1-D array
        !   PSIBIT 1:nsilop
        !   BITMPI 1:magpri
        !   BITFC  1:nfcoil
        !   BITEC  1:nesum ! added by MK
        !   IERMPI to fix FWTMP2 1:magpri! added by MK
        parameter (nsize2=5+nsilop+magpri+nfcoil+nesum+magpri)
        integer :: i,j,ktime_all,offset
        integer,dimension(:),allocatable :: tmp1,tmp2
        double precision :: zwork(nsize,ntime),zwork2(nsize2),timeb_list(nproc)
        ! TIMING >>>
        integer :: clock0,clock1,clockmax,clockrate,ticks
        double precision :: secs
        ! TIMING <<<
        integer :: total_bytes
        zwork(:,:) = 0.0
        allocate(tmp1(nproc),tmp2(nproc))

        efitversion = 20201013

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
          secs = dble(ticks)/dble(clockrate)
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
          zwork2(1) = dble(iavem)   ! INT4  (1)
          zwork2(2) = dble(limitr)  ! INT4  (1)
          zwork2(3) = bitip         ! REAL8 (1)
          zwork2(4) = rcentr        ! REAL8 (1)
          zwork2(5) = dble(oldccomp) ! added by MK 2020.10.07
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
          oldccomp = int(zwork2(5)) ! Added by MK 2020.10.07
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
          call MPI_SCATTERV(zwork,tmp1,tmp2,MPI_DOUBLE_PRECISION,MPI_IN_PLACE,tmp1(rank+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        else
          call MPI_SCATTERV(zwork,tmp1,tmp2,MPI_DOUBLE_PRECISION,zwork,       tmp1(rank+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
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
          secs = dble(ticks)/dble(clockrate)
          write (*,"(' GETPTS transfer ',i10,' bytes in ',f6.2,' sec')") total_bytes,secs
        endif
        ! TIMING <<<

      end subroutine getpts_mpi

! =========================

      ! NOTE : NO error condition returned
      subroutine getstark_mpi(ktime)

        include 'eparmdud129.f90'
        include 'modules1.f90'

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
        parameter (nsize=nmtark*12)
        integer :: i,j,ktime_all
        integer,dimension(:),allocatable :: tmp1,tmp2
        double precision :: zwork(nsize,ntime)
        ! TIMING >>>
        integer :: clock0,clock1,clockmax,clockrate,ticks
        double precision :: secs
        ! TIMING <<<
        integer :: total_bytes
        
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
          secs = dble(ticks)/dble(clockrate)
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
        total_bytes = total_bytes + 8*sum(dist_data(2:))*nsize
        if (rank == 0) then
          ! NOTE : DIST_DATA and DIST_DATA_DISPLS should be saved between calls since part of MPI_INFO module
          call MPI_SCATTERV(zwork,tmp1,tmp2,MPI_DOUBLE_PRECISION,MPI_IN_PLACE,tmp1(rank+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        else
          call MPI_SCATTERV(zwork,tmp1,tmp2,MPI_DOUBLE_PRECISION,zwork,tmp1(rank+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
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
        call MPI_BCAST(spatial_avg_gam,nstark*ngam_vars*ngam_u*ngam_w,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
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
          secs = dble(ticks)/dble(clockrate)
          write (*,"(' GETSTARK transfer ',i10,' bytes in ',f6.2,' sec')") total_bytes,secs
        endif
        ! TIMING <<<

      end subroutine getstark_mpi

#endif
! MPI <<<
