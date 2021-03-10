      subroutine fcurrt(jtime,iter,itertt,kerror)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          fcurrt computes the currents in the f coils.            **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**          24/07/85..........revised                               **
!**          16/08/90..........revised                               **
!**                                                                  **
!**********************************************************************
      use commonblocks,only: rfcpc
      include 'eparmdud129.f90'
      include 'modules2.f90'
      include 'modules1.f90'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      dimension afma(nfcoil,nfcoil),ifmatr(nfcoil),wfmatr(nfcoil) &
           ,wbry(msbdry),work(msbdr2),ut(msbdry,msbdry)
      dimension abry(msbdry,nfcoil+nvesel),bbry(msbdry),ainbry(nfcoil+nvesel,msbdry)
      dimension fcref(nfcoil)
      dimension pbry(msbdry),psbry(msbdry)
      integer, intent(inout) :: kerror
      data iskip/0/,idoit/0/

      kerror = 0

      if (ifcurr.gt.0) return
      if (itertt.le.1.and.icinit.lt.0) return
      if (islpfc.ne.1) go to 2000
!-----------------------------------------------------------------------
!--  flux loop on F coils, seldom used for DIII-D                     --
!-----------------------------------------------------------------------
      if (iskip.gt.0) go to 1200
      iskip=1
!vas
!vas      print*,'file name : ','fc'//trim(ch1)// &
!vas                         trim(ch2)//'.ddd'
      open(unit=nffile,status='old',form='unformatted', &
       file=table_dir(1:ltbdir)//'fc'//trim(ch1)//trim(ch2)//'.ddd')
        read (nffile) rfcfc
        read (nffile) rfcpc
        close(unit=nffile)
  !vas
  !vas      print*,'file name : ','fm'//trim(ch1)// &
  !vas                         trim(ch2)//'.ddd'
        open(unit=nffile,status='old',form='unformatted', &
       file=table_dir(1:ltbdir)//'fm'//trim(ch1)//trim(ch2)//'.ddd')
      read (nffile,err=1150) afma
      read (nffile) ifmatr
      close(unit=nffile)
      go to 1200
 1150 continue
      do 1170 i=1,nfcoil
      do 1170 j=1,nfcoil
        afma(i,j)=rfcfc(i,j)
 1170 continue
      m111=-1.0
      call decomp(nfcoil,nfcoil,afma,m111,ifmatr,wfmatr)
!vas
!vas      print*,'file name : ','fm'//trim(ch1)// &
!vas                         trim(ch2)//'.ddd'
      open(unit=nffile,status='old',form='unformatted', &
           file='fm'//trim(ch1)//trim(ch2)//'.ddd' &
              ,err=12927)
      close(unit=nffile,status='delete')
12927 continue
!vas
!vas      print*,'file name : ','fm'//trim(ch1)// &
!vas                          trim(ch2)//'.ddd'
      open(unit=nffile,status='new',form='unformatted', &
           file='fm'//trim(ch1)//trim(ch2)//'.ddd')
      write (nffile) afma
      write (nffile) ifmatr
      close(unit=nffile)
 1200 continue
      do 1320 i=1,nfcoil
        brsp(i)=0.0
        if (ivacum.gt.0) go to 1310
        do 1308 j=1,nwnh
        brsp(i)=brsp(i)+rfcpc(i,j)*pcurrt(j)
 1308   continue
 1310   continue
        if (iecurr.le.0) go to 1314
        do 1312 j=1,nesum
          brsp(i)=brsp(i)+rfcec(i,j)*ecurrt(j)
 1312   continue
 1314   continue
        if (ivesel.le.0.or.ifitvs.gt.0) go to 1318
        do 1316 j=1,nvesel
          brsp(i)=brsp(i)+rfcvs(i,j)*vcurrt(j)
 1316   continue
 1318   continue
        brsp(i)=csilop(i,jtime)-brsp(i)
        if (fitsiref) brsp(i)=brsp(i)+psiref(jtime)
 1320 continue
      call solve(nfcoil,nfcoil,afma,brsp,ifmatr)
      return
!-----------------------------------------------------------------------
!--  standard option, flux loops away from F coils                    --
!-----------------------------------------------------------------------
 2000 continue
      if (idoit.gt.0.and.itertt.gt.1) go to 3000
!-----------------------------------------------------------------------
!--  set up response matrix once for all                              --
!-----------------------------------------------------------------------
      idoit=0 ! Set idoit =0, this does save after setting to allocatable
      wsibry=psibry
!-----------------------------------------------------------------------
!--   get fixed boundary response from plasma                         --
!-----------------------------------------------------------------------
      if (nbdry.le.0.or.iconvr.ne.3) go to 2098
!-----------------------------------------------------------------------
!-- set up boundary fitting weights                                   --
!-----------------------------------------------------------------------
      z04=1.0e-04_dp
      fwtbdr=abs(errbry)*max(abs(sidif),z04)
      fwtbdr=1.0/fwtbdr
      do 2080 i=1,nbdry
        fwtbry(i)=fwtbdr*fwtbdry(i)
 2080 continue
 2098 continue
      if (idebug >=2) then
        write (106,*) 'FCURRT NBDRY,PSIBRY,WSIRY = ', nbdry,psibry,wsibry
        write (106,*) '       SIDIF,FWTBDR,ERRBRY= ', sidif,fwtbdr,errbry
        write (106,*) '       PSIBRY0= ', psibry0
      endif
!-----------------------------------------------------------------------
!--  set up response matrix, first fixed boundary, then flux loops    --
!--  and F coil currents                                              --
!-----------------------------------------------------------------------
      do 2500 nk=1,nfcoil
      nj=0
      if (nbdry.le.0.or.iconvr.ne.3) go to 2200
      do 2100 m=1,nbdry
        abry(m,nk)=fwtbry(m)*rbdrfc(m,nk)
 2100 continue
      nj=nbdry
 2200 continue
!-----------------------------------------------------------------------
!--  flux loop option                                                 --
!-----------------------------------------------------------------------
      do 2220 m=1,nsilop
        if (fwtsi(m).le.0.0) go to 2220
        nj=nj+1
        abry(nj,nk)=fwtsi(m)*rsilfc(m,nk)
 2220 continue
!-----------------------------------------------------------------------
!--  magnetic probes used only for vacuum error field analysis        --
!-----------------------------------------------------------------------
      if (ivacum.gt.0) then
      do 2222 m=1,magpri
        if (fwtmp2(m).le.0.0) go to 2222
        nj=nj+1
        abry(nj,nk)=fwtmp2(m)*rmp2fc(m,nk)
 2222 continue
      endif
!-----------------------------------------------------------------------
!--  F coil currents option                                           --
!-----------------------------------------------------------------------
      do 2225 m=1,nfcoil
        if (fwtfc(m).le.0.0) go to 2225
        nj=nj+1
        abry(nj,nk)=0.
        if (m.eq.nk) abry(nj,nk)=fwtfc(m)
 2225 continue
!-----------------------------------------------------------------------
!--  minimize coil current oscillations for fixed boundary option     --
!-----------------------------------------------------------------------
      if (nbdry.le.0.or.iconvr.ne.3) go to 2250
      do 2240 m=1,nfcoil
        nj=nj+1
        abry(nj,nk)=0.0
        if (nk.eq.m) abry(nj,nk)=cfcoil*fczero(m)
 2240 continue
!-----------------------------------------------------------------------
!--  constraint F coil currents to be symmetric if requested          --
!-----------------------------------------------------------------------
      if ((symmetrize).and.(cupdown.ne.0.0)) then
      do m=1,nfcoil/2
        nj=nj+1
        abry(nj,nk)=0.0
        if (nk.eq.m) abry(nj,nk)=cupdown
        if (nk.eq.(m+nfcoil/2)) abry(nj,nk)=-cupdown
      enddo
      endif
 2250 continue
!-----------------------------------------------------------------------
!--  F coil constraints Cf=x                                          --
!-----------------------------------------------------------------------
      if (kccoils.gt.0) then
      do m=1,kccoils
        nj=nj+1
        abry(nj,nk)=ccoils(nk,m)
      enddo
      endif
!-----------------------------------------------------------------------
!--  SOL constraints 2014/05/20 LL                                    --
!-----------------------------------------------------------------------
      if (nsol.gt.0) then
!-----------------------------------------------------------------------
!-- set up SOL fitting weights and response matrix                    --
!-----------------------------------------------------------------------
        z04s=1.0e-04_dp
        fwtsols=abs(errbry)*max(abs(sidif),z04s)
        fwtsols=1.0/fwtsols
        do i=1,nsol
          fwtsolw(i)=fwtsols*fwtsol(i)
          nj=nj+1
          abry(nj,nk)=fwtsolw(i)*rsolfc(i,nk)
        enddo
        if (idebug >= 2) write (6,*) 'FCURRT nsol,fwtsolw = ', nsol,(fwtsolw(i),i=1,nsol)
      endif
 2500 continue
!-----------------------------------------------------------------------
!--  optional vessel current model                                    --
!-----------------------------------------------------------------------
      neqn=nfcoil
      if (ifitvs.le.0) go to 2610
!
      if (nfourier.gt.1) then
        nuuu=nfourier*2+1
      else
        nuuu=nvesel
      endif
      do 2600 nkk=1,nuuu
      nk=nkk+nfcoil
      nj=0
      if (nbdry.le.0.or.iconvr.ne.3) go to 2510
      do 2505 m=1,nbdry
          if(nfourier.gt.1)then
          temp=0.
          do i=1,nvesel
          temp=temp+rbdrvs(m,i)*vecta(nkk,i)
          enddo
            abry(m,nk)=fwtbry(m)*temp
          else
            abry(m,nk)=fwtbry(m)*rbdrvs(m,nkk)
          endif
 2505 continue
      nj=nbdry
 2510 continue
      do 2520 m=1,nsilop
        if (fwtsi(m).le.0.0) go to 2520
        nj=nj+1
        if(nfourier.gt.1)then
          temp=0.
          do i=1,nvesel
          temp=temp+rsilvs(m,i)*vecta(nkk,i)
          enddo
          abry(nj,nk)=fwtsi(m)*temp
        else
           abry(nj,nk)=fwtsi(m)*rsilvs(m,nkk)
        endif
 2520 continue
!-----------------------------------------------------------------------
!-- magnetic probes used only for vacuum error analysis               --
!-----------------------------------------------------------------------
      if (ivacum.gt.0) then
      do 2522 m=1,magpri
        if (fwtmp2(m).le.0.0) go to 2522
        nj=nj+1
        if(nfourier.gt.1)then
          temp=0.
          do i=1,nvesel
          temp=temp+rmp2vs(m,i)*vecta(nkk,i)
          enddo
           abry(nj,nk)=fwtmp2(m)*temp
        else
           abry(nj,nk)=fwtmp2(m)*rmp2vs(m,nkk)
        endif
 2522 continue
      endif
      do 2525 m=1,nfcoil
        if (fwtfc(m).le.0.0) go to 2525
        nj=nj+1
        abry(nj,nk)=0.
 2525 continue
      if (nbdry.le.0.or.iconvr.ne.3) go to 2550
      do 2540 m=1,nfcoil
        nj=nj+1
        abry(nj,nk)=0.0
 2540 continue
 2550 continue
 2600 continue
      neqn=neqn+nvesel
 2610 continue
!-----------------------------------------------------------------------
!--  set up working matrix                                            --
!-----------------------------------------------------------------------
      do 2800 i=1,msbdry
      do 2800 j=1,msbdry
        ut(i,j)=0.0
 2800 continue
      do 2870 i=1,msbdry
 2870 ut(i,i) = 1.0
!-----------------------------------------------------------------------
!--  compute inverse of abry in least squares sence with singular     --
!--  value decomposition                                              --
!-----------------------------------------------------------------------
      if (idebug >= 2) write (106,*) 'FCURRT nj,neqn = ', nj,neqn
      call sdecm(abry,msbdry,nj,neqn,ut,msbdry,nj,wbry,work,ier)
      if (ier.eq.129) then
        kerror = 1
        call errctrl_msg('fcurrt','sdecm failed to converge')
        return
      end if
      cond=ier
      toler=1.0e-06_dp*wbry(1)
      do 2890 i=1,neqn
      do 2890 j=1,nj
        ainbry(i,j)=0.0
        do 2890 k=1,neqn
          t=0.0
          if (wbry(k).gt.toler) t=1.0/wbry(k)
          ainbry(i,j) = ainbry(i,j) + abry(i,k)*ut(k,j)*t
 2890   continue
!-----------------------------------------------------------------------
!--  compute and store coil reference currents, not need for          --
!--  non-fixed boundary when IFREF=-1                                 --
!-----------------------------------------------------------------------
      if (ifref.le.0.or.iconvr.ne.3) go to 3000
      do 2900 i=1,nj
        if (i.le.nbdry) then
          bbry(i)=fwtbry(i)
        else
          bbry(i)=0.0
        endif
 2900 continue
      do 2952 i=1,nfcoil
        fcref(i)=0.0
        do 2952 j=1,nj
          fcref(i)=fcref(i)+ainbry(i,j)*bbry(j)
 2952   continue
!-----------------------------------------------------------------------
!--  RHS of response matrix, start here if got inverse matrix already --
!-----------------------------------------------------------------------
 3000 continue
      nj=0
      if (nbdry.le.0.or.iconvr.ne.3) go to 3600
!-----------------------------------------------------------------------
!-- fixed boundary portion                                            --
!-----------------------------------------------------------------------
      do 3596 m=1,nbdry
        bbry(m)=0.0
        if (ivesel.le.0.or.ifitvs.gt.0) go to 3570
        do 3565 k=1,nvesel
          bbry(m)=bbry(m)+rbdrvs(m,k)*vcurrt(k)
 3565   continue
 3570   continue
        if (iecurr.le.0) go to 3590
        do 3580 k=1,nesum
          bbry(m)=bbry(m)+rbdrec(m,k)*ecurrt(k)
 3580   continue
 3590   continue
        if (iacoil.gt.0) then
          do 3592 k=1,nacoil
            bbry(m)=bbry(m)+rbdrac(m,k)*caccurt(jtime,k)
 3592     continue
        endif
        do 3594 k=1,nwnh
         bbry(m)=bbry(m)+rbdrpc(m,k)*pcurrt(k)
 3594   continue
        pbry(m)=bbry(m)
        bbry(m)=fwtbry(m)*(wsibry-bbry(m))
 3596 continue
      nj=nbdry
      if (idebug >= 2) then
        write (106,*) 'ivesel,iecurr,iacoil= ',ivesel,iecurr,iacoil
        write (106,*) 'pbry= ',(pbry(m),m=1,nbdry)
        write (106,*) 'bbry= ',(bbry(m),m=1,nbdry)
        write (106,*) 'pbry(24),bbry(24)= ',pbry(24),bbry(24)
        write (106,*) 'rgrid,zgrid= ',rgrid(75),zgrid(104)
        write (106,*) 'rbdry,zbdry= ',rbdry(24),zbdry(24)
        bdrav=abs(rbdrpc(24,1))
        bdrmax=abs(rbdrpc(24,1))
        mbbrmax=1
        do m=2,nwnh
          bdrav=bdrav+abs(rbdrpc(24,m))
          if (bdrmax.lt.abs(rbdrpc(24,m))) then
             bdrmax=abs(rbdrpc(24,m))
             mbdrmax=m
          endif
        enddo
        bdrav=bdrav/nwnh
        write (106,*) 'mbdrmax,bdrav,bdrmax= ',mbdrmax,bdrav,bdrmax
        write (106,*) 'rbdrpc(24,*)= ',(rbdrpc(24,m),m=1,nwnh)
        write (106,*) 'rbdrec(24,*)= ',(rbdrec(24,m),m=1,nesum)
        zdif=zbdry(24)-zgrid(104)
        bdrmax=psical(rbdry(24),rgrid(75),zdif)*tmu
        bdrmaxs=gridpc(9547,75)
        write (106,*) 'mbdrmax,bdrmaxs,bdrmax= ',mbdrmax,bdrmaxs,bdrmax
      endif
!-----------------------------------------------------------------------
!--  flux loops portion                                               --
!-----------------------------------------------------------------------
 3600 continue
      do 3630 m=1,nsilop
        if (fwtsi(m).le.0.0) go to 3630
        nj=nj+1
        bbry(nj)=0.0
        do 3610 j=1,nwnh
        bbry(nj)=bbry(nj)+gsilpc(m,j)*pcurrt(j)
 3610   continue
        if (iecurr.le.0) go to 3614
        do 3612 j=1,nesum
          bbry(nj)=bbry(nj)+rsilec(m,j)*ecurrt(j)
 3612   continue
 3614   continue
!-----------------------------------------------------------------------
!-- specify vessel currents ?                                         --
!-----------------------------------------------------------------------
        if (ivesel.le.0.or.ifitvs.gt.0) go to 3618
        do 3616 j=1,nvesel
          bbry(nj)=bbry(nj)+rsilvs(m,j)*vcurrt(j)
 3616   continue
 3618   continue
!-----------------------------------------------------------------------
!-- use advance divertor coil ?                                       --
!-----------------------------------------------------------------------
        if (iacoil.gt.0) then
          do 3621 j=1,nacoil
            bbry(nj)=bbry(nj)+rsilac(m,j)*caccurt(jtime,j)
 3621     continue
        endif
        if (fitsiref) then
        bbry(nj)=fwtsi(m)*(silopt(jtime,m)+psiref(jtime)-bbry(nj))
        else
        bbry(nj)=fwtsi(m)*(silopt(jtime,m)-bbry(nj))
        endif
 3630 continue
!-----------------------------------------------------------------------
!--  magnetic probes used only for vacuum error field analysis        --
!-----------------------------------------------------------------------
      if (ivacum.gt.0) then
      do 3650 m=1,magpri
        if (fwtmp2(m).le.0.0) go to 3650
        nj=nj+1
        bbry(nj)=0.0
        do 3640 j=1,nwnh
        bbry(nj)=bbry(nj)+gmp2pc(m,j)*pcurrt(j)
 3640   continue
        if (iecurr.le.0) go to 3644
        do 3642 j=1,nesum
          bbry(nj)=bbry(nj)+rmp2ec(m,j)*ecurrt(j)
 3642   continue
 3644   continue
!-----------------------------------------------------------------------
!-- specify vessel currents ?                                         --
!-----------------------------------------------------------------------
        if (ivesel.le.0.or.ifitvs.gt.0) go to 3648
        do 3646 j=1,nvesel
          bbry(nj)=bbry(nj)+rmp2vs(m,j)*vcurrt(j)
 3646   continue
 3648   continue
        if (iacoil.gt.0) then
        do 3649 j=1,nacoil
          bbry(nj)=bbry(nj)+rmp2ac(m,j)*caccurt(jtime,j)
 3649   continue
        endif
        bbry(nj)=fwtmp2(m)*(expmpi(jtime,m)-bbry(nj))
 3650 continue
      endif
!-----------------------------------------------------------------------
!-- F-coil currents specification                                     --
!-----------------------------------------------------------------------
      do 3660 m=1,nfcoil
        if (fwtfc(m).le.0.0) go to 3660
        nj=nj+1
        bbry(nj)=fccurt(jtime,m)*fwtfc(m)
 3660 continue
!-----------------------------------------------------------------------
!--  F coil current minimization                                      --
!-----------------------------------------------------------------------
      if (nbdry.le.0.or.iconvr.ne.3) go to 3665
      do 3662 m=1,nfcoil
        nj=nj+1
        bbry(nj)=0.0
 3662 continue
!-----------------------------------------------------------------------
!--  symmetrize F coil currents ?                                     --
!-----------------------------------------------------------------------
      if ((symmetrize).and.(cupdown.ne.0.0)) then
      do m=1,nfcoil/2
        nj=nj+1
        bbry(nj)=0.0
      enddo
      endif
 3665 continue
!-----------------------------------------------------------------------
!--  F coil constraints ?                                             --
!-----------------------------------------------------------------------
      if (kccoils.gt.0) then
        do m=1,kccoils
          nj=nj+1
          bbry(nj)=xcoils(m)
        enddo
      endif
!-----------------------------------------------------------------------
!-- SOL Constraints RHS                                               --
!-----------------------------------------------------------------------
      if (nsol.gt.0) then
        if (idebug >= 2) write (6,*) 'FCURRT wsisol = ', wsisol
        do m=1,nsol
          nj=nj+1
          bbry(nj)=0.0
          if (ivesel.gt.0.and.ifitvs.lt.0) then
            do k=1,nvesel
              bbry(nj)=bbry(nj)+rsolvs(m,k)*vcurrt(k)
            enddo
          endif
          if (iecurr.gt.0) then
             do k=1,nesum
                bbry(nj)=bbry(nj)+rsolec(m,k)*ecurrt(k)
             enddo
           endif
           do k=1,nwnh
             bbry(nj)=bbry(nj)+rsolpc(m,k)*pcurrt(k)
           enddo
           psbry(m)=bbry(nj)
           bbry(nj)=fwtsolw(m)*(wsisol-bbry(nj))
        enddo
      endif
!-----------------------------------------------------------------------
!--  now get F coil currents from precomputed inverse matrix          --
!-----------------------------------------------------------------------
      do 3670 i=1,nfcoil
        brsp(i)=0.0
        do 3670 j=1,nj
          brsp(i)=brsp(i)+ainbry(i,j)*bbry(j)
 3670 continue
      if (ifitvs.le.0) go to 3685
      do 3680 ii=1,nvesel
        i=ii+nfcoil
        vcurrt(ii)=0.0
        do 3680 j=1,nj
          vcurrt(ii)=vcurrt(ii)+ainbry(i,j)*bbry(j)
 3680 continue
!
       if(nfourier.gt.1) then
        do j=1,nvesel
        temp=0.
        do 2671 i=1,(nfourier*2+1)
          temp=temp+vcurrt(i)*vecta(i,j)
 2671   continue
        vcurrt(j)=temp
        enddo
       endif
!
 3685 continue
!-----------------------------------------------------------------------
!--  adjustments for various boundary flux options                    --
!-----------------------------------------------------------------------
      if (ifref.le.0.or.iconvr.ne.3) go to 4000
      sumif = 0.0
      sumifr = 0.0
!-----------------------------------------------------------------------
!-- sum of inner F-coils 1-5 A and B zero IFREF=1                     --
!-----------------------------------------------------------------------
      if (ifref.ne.1) go to 3700
      do 3690 i=1,5
 3690 sumif = sumif + brsp(i) + brsp(i+9)
      sumif = sumif + brsp(8) + brsp(17)
      do 3695 i=1,5
 3695 sumifr = sumifr + fcref(i) + fcref(i+9)
      sumifr = sumifr + fcref(8) + fcref(17)
      go to 3750
 3700 continue
!-----------------------------------------------------------------------
!--  sum of F coils selected through FCSUM vanish IFREF=2             --
!-----------------------------------------------------------------------
      if (ifref.ne.2) go to 3720
      do 3705 i=1,nfcoil
 3705 sumif=sumif+fcsum(i)*brsp(i)/turnfc(i)
      do 3710 i=1,nfcoil
 3710 sumifr = sumifr + fcref(i) * fcsum(i)/turnfc(i)
      go to 3750
 3720 continue
!----------------------------------------------------------------------
!--  choose boundary flux by minimize coil currents IFREF=3          --
!----------------------------------------------------------------------
      if (ifref.eq.3) then
      do 3725 i=1,nfcoil
 3725 sumif=sumif+fcref(i)*brsp(i)*fczero(i)
      do 3730 i=1,nfcoil
 3730 sumifr = sumifr + fcref(i)**2*fczero(i)
      go to 3750
      endif
 3750 continue
!-----------------------------------------------------------------------
!--  update boundary flux for IFREF=1-3                               --
!-----------------------------------------------------------------------
      if (ifref.le.3) then
      ssiref = sumif/sumifr
      do 3759 m=1,nfcoil
      silopt(jtime,m)=silopt(jtime,m)-ssiref
 3759 brsp(m) = brsp(m) - ssiref*fcref(m)
      wsibry=wsibry-ssiref
      wsisol=wsisol-ssiref
      endif
!------------------------------------------------------------------------
!--  fixed boundary flux specified through PSIBRY, IFREF=4             --
!------------------------------------------------------------------------
      if (ifref.eq.4) then
      ssiref=psibry0-psibry
      do 3770 m=1,nfcoil
        brsp(m)=brsp(m)+ssiref*fcref(m)
        silopt(jtime,m)=silopt(jtime,m)+ssiref
 3770 continue
      wsibry=wsibry+ssiref
      wsisol=wsisol+ssiref
      endif
      if (idebug >= 2) then
        write (106,*) '      PSIBRY0,PSIBRY,WSIBRY= ',psibry0,psibry,wsibry
      endif
!-----------------------------------------------------------------------
!--  done, estimated errors for fixed boundary calculations           --
!-----------------------------------------------------------------------
 4000 continue
      if (nbdry.gt.0.and.iconvr.eq.3) then
      erbmax=0.0
      erbave=0.0
      do 4100 i=1,nbdry
        xsibry=pbry(i)
        do 4050 m=1,nfcoil
          xsibry=xsibry+rbdrfc(i,m)*brsp(m)
 4050   continue
        erbloc(i)=abs((wsibry-xsibry)/sidif)
        erbmax=max(erbloc(i),erbmax)
        erbave=erbave+erbloc(i)
 4100 continue
      erbave=erbave/nbdry
      if (idebug /= 0) write (106,*) 'FCURRT erbmax,erbave,si = ', erbmax,erbave,wsibry
      if (idebug >= 2) then
         write (106,*) 'XSIBRY,PBRY(1),BRSP= ',xsibry,pbry(1),(brsp(m),m=1,nfcoil)
      endif
      endif
!
      if (nsol.gt.0.and.iconvr.eq.3) then
      erbsmax=0.0
      erbsave=0.0
      do i=1,nsol
        xsisol=psbry(i)
        do m=1,nfcoil
          xsisol=xsisol+rsolfc(i,m)*brsp(m)
        enddo
        erbsloc(i)=abs((wsisol-xsisol)/sidif)
        erbsmax=max(erbsloc(i),erbsmax)
        erbsave=erbsave+erbsloc(i)
      enddo
      erbsave=erbsave/nsol
      if (idebug /= 0) write (6,*) 'FCURRT erbsmax,erbsave,si = ', erbsmax,erbsave,wsisol
      endif
!
      return
      end
