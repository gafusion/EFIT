#include "config.f"
!**********************************************************************
!>
!!    fcurrt computes the currents in the f coils.
!!
!!
!!    @param jtime : time index
!!
!!    @param iter : current profile (outer) loop iteration index
!!
!!    @param itertt : total iteration index (current+equilibirum loops)
!!
!!    @param kerror : error flag
!!
!**********************************************************************
      subroutine fcurrt(jtime,iter,itertt,kerror)
      use commonblocks,only: rfcpc
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      integer*4, intent(inout) :: kerror
      dimension afma(nfcoil,nfcoil),ifmatr(nfcoil),wfmatr(nfcoil), &
                wbry(msbdry),work(msbdr2),ut(msbdry,msbdry)
      dimension abry(msbdry,nfcoil+nvesel),bbry(msbdry), &
                ainbry(nfcoil+nvesel,msbdry)
      dimension fcref(nfcoil)
      dimension pbry(msbdry),psbry(msbdry)
      integer*4, parameter :: islpfc=0 ! hardcoded option
      data iskip/0/,idoit/0/

      kerror = 0
      if(ifcurr.gt.0) return
      if(itertt.le.1.and.icinit.lt.0) return
      if(islpfc.eq.1) then
!-----------------------------------------------------------------------
!--     flux loop on F coils, no longer used...
!-----------------------------------------------------------------------
        if (iskip.le.0) then
          iskip=1

          open(unit=nffile,status='old',form='unformatted', &
          file=table_di2(1:ltbdi2)//'fc'//trim(ch1)//trim(ch2)//'.ddd')
          read (nffile) rfcfc
          read (nffile) rfcpc
          close(unit=nffile)
!vas
!vas      print*,'file name : ','fm'//trim(ch1)// &
!vas                         trim(ch2)//'.ddd'
          open(unit=nffile,status='old',form='unformatted', &
          file=table_di2(1:ltbdi2)//'fm'//trim(ch1)//trim(ch2)//'.ddd')
          read (nffile,iostat=ioerr) afma
          if (ioerr.eq.0) then
            read (nffile) ifmatr
            close(unit=nffile)
          else
            afma(1:nfcoil,1:nfcoil)=rfcfc(1:nfcoil,1:nfcoil)
            m111=-1.0
            call decomp(nfcoil,nfcoil,afma,m111,ifmatr,wfmatr)
!vas
!vas      print*,'file name : ','fm'//trim(ch1)// &
!vas                         trim(ch2)//'.ddd'
            open(unit=nffile,status='old',form='unformatted', &
                 file='fm'//trim(ch1)//trim(ch2)//'.ddd',iostat=ioerr)
            if(ioerr.eq.0) close(unit=nffile,status='delete')
!vas
!vas      print*,'file name : ','fm'//trim(ch1)// &
!vas                          trim(ch2)//'.ddd'
            open(unit=nffile,status='new',form='unformatted', &
                 file='fm'//trim(ch1)//trim(ch2)//'.ddd')
            write (nffile) afma
            write (nffile) ifmatr
            close(unit=nffile)
          endif
        endif

        do i=1,nfcoil
          brsp(i)=0.0
          if(ivacum.le.0) &
            brsp(i)=brsp(i)+sum(rfcpc(i,1:nwnh)*pcurrt(1:nwnh))
          if(iecurr.gt.0) &
            brsp(i)=brsp(i)+sum(rfcec(i,1:nesum)*ecurrt(1:nesum))
          if(ivesel.gt.0.and.ifitvs.le.0) &
            brsp(i)=brsp(i)+sum(rfcvs(i,1:nvesel)*vcurrt(1:nvesel))
          brsp(i)=csilop(i,jtime)-brsp(i)
          if(fitsiref) brsp(i)=brsp(i)+psiref(jtime)
        enddo
        call solve(nfcoil,nfcoil,afma,brsp,ifmatr)
        return
      endif
!-----------------------------------------------------------------------
!--    standard option, flux loops away from F coils                 --
!-----------------------------------------------------------------------
      first_iter: if (idoit.le.0.or.itertt.le.1) then
!-----------------------------------------------------------------------
!--    set up response matrix once for all                            --
!-----------------------------------------------------------------------
       idoit=0 !Set idoit=0, this does save after setting to allocatable
       wsibry=psibry
!-----------------------------------------------------------------------
!--    get fixed boundary response from plasma                        --
!-----------------------------------------------------------------------
       if (nbdry.gt.0.and.iconvr.eq.3) then
!-----------------------------------------------------------------------
!--      set up boundary fitting weights                              --
!-----------------------------------------------------------------------
         z04=1.0e-04_dp
         fwtbdr=abs(errbry)*max(abs(sidif),z04)
         fwtbdr=1.0/fwtbdr
         fwtbry(1:nbdry)=fwtbdr*fwtbdry(1:nbdry)
       endif
#ifdef DEBUG_LEVEL2
      write (106,*) 'FCURRT NBDRY,PSIBRY,WSIRY = ', nbdry,psibry,wsibry
      write (106,*) '       SIDIF,FWTBDR,ERRBRY= ', sidif,fwtbdr,errbry
      write (106,*) '       PSIBRY0= ', psibry0
#endif
!-----------------------------------------------------------------------
!--    set up response matrix, first fixed boundary, then flux loops  --
!--    and F coil currents                                            --
!-----------------------------------------------------------------------
       do nk=1,nfcoil
         nj=0
         if (nbdry.gt.0.and.iconvr.eq.3) then
           abry(1:nbdry,nk)=fwtbry(1:nbdry)*rbdrfc(1:nbdry,nk)
           nj=nbdry
         endif
!-----------------------------------------------------------------------
!--      flux loop option                                             --
!-----------------------------------------------------------------------
         do m=1,nsilop
           if(fwtsi(m).le.0.0) cycle
           nj=nj+1
           abry(nj,nk)=fwtsi(m)*rsilfc(m,nk)
         enddo
!-----------------------------------------------------------------------
!--      magnetic probes used only for vacuum error field analysis    --
!-----------------------------------------------------------------------
         if (ivacum.gt.0) then
           do m=1,magpri
             if(fwtmp2(m).le.0.0) cycle
             nj=nj+1
             abry(nj,nk)=fwtmp2(m)*rmp2fc(m,nk)
           enddo
         endif
!-----------------------------------------------------------------------
!--      F coil currents option                                       --
!-----------------------------------------------------------------------
         do m=1,nfcoil
           if(fwtfc(m).le.0.0) cycle
           nj=nj+1
           abry(nj,nk)=0.
           if (m.eq.nk) abry(nj,nk)=fwtfc(m)
         enddo
!-----------------------------------------------------------------------
!--      minimize coil current oscillations for fixed boundary option --
!-----------------------------------------------------------------------
         if (nbdry.gt.0.and.iconvr.eq.3) then
           do m=1,nfcoil
             nj=nj+1
             abry(nj,nk)=0.0
             if(nk.eq.m) abry(nj,nk)=cfcoil*fczero(m)
           enddo
!-----------------------------------------------------------------------
!--        constraint F coil currents to be symmetric if requested    --
!-----------------------------------------------------------------------
           if ((symmetrize).and.(cupdown.ne.0.0)) then
             do m=1,nfcoil/2
               nj=nj+1
               abry(nj,nk)=0.0
               if(nk.eq.m) abry(nj,nk)=cupdown
               if(nk.eq.(m+nfcoil/2)) abry(nj,nk)=-cupdown
             enddo
           endif
         endif
!-----------------------------------------------------------------------
!--      F coil constraints Cf=x                                      --
!-----------------------------------------------------------------------
         if (kccoils.gt.0) then
           do m=1,kccoils
             nj=nj+1
             abry(nj,nk)=ccoils(nk,m)
           enddo
         endif
!-----------------------------------------------------------------------
!--      SOL constraints                                              --
!-----------------------------------------------------------------------
         if (nsol.gt.0) then
!-----------------------------------------------------------------------
!--        set up SOL fitting weights and response matrix             --
!-----------------------------------------------------------------------
           z04s=1.0e-04_dp
           fwtsols=abs(errbry)*max(abs(sidif),z04s)
           fwtsols=1.0/fwtsols
           do i=1,nsol
             fwtsolw(i)=fwtsols*fwtsol(i)
             nj=nj+1
             abry(nj,nk)=fwtsolw(i)*rsolfc(i,nk)
           enddo
#ifdef DEBUG_LEVEL2
           write (6,*) 'FCURRT nsol,fwtsolw = ', &
                       nsol,(fwtsolw(i),i=1,nsol)
#endif
         endif
       enddo
!-----------------------------------------------------------------------
!--    optional vessel current model                                  --
!-----------------------------------------------------------------------
       neqn=nfcoil
       fit_vessel: if (ifitvs.gt.0) then
!
         if (nfourier.gt.1) then
           nuuu=nfourier*2+1
         else
           nuuu=nvesel
         endif
         do nkk=1,nuuu
           nk=nkk+nfcoil
           nj=0
           if (nbdry.gt.0.and.iconvr.eq.3) then
             do m=1,nbdry
               if (nfourier.gt.1) then
                 abry(m,nk)=fwtbry(m)* &
                   sum(rbdrvs(m,1:nvesel)*vecta(nkk,1:nvesel))
               else
                 abry(m,nk)=fwtbry(m)*rbdrvs(m,nkk)
               endif
             enddo
             nj=nbdry
           endif
           do m=1,nsilop
             if(fwtsi(m).le.0.0) cycle
             nj=nj+1
             if (nfourier.gt.1) then
               abry(nj,nk)=fwtsi(m)* &
                 sum(rsilvs(m,1:nvesel)*vecta(nkk,1:nvesel))
             else
               abry(nj,nk)=fwtsi(m)*rsilvs(m,nkk)
             endif
           enddo
!-----------------------------------------------------------------------
!--        magnetic probes used only for vacuum error analysis        --
!-----------------------------------------------------------------------
           if (ivacum.gt.0) then
             do m=1,magpri
               if(fwtmp2(m).le.0.0) cycle
               nj=nj+1
               if (nfourier.gt.1) then
                 abry(nj,nk)=fwtmp2(m)* &
                   sum(rmp2vs(m,1:nvesel)*vecta(nkk,1:nvesel))
               else
                 abry(nj,nk)=fwtmp2(m)*rmp2vs(m,nkk)
               endif
             enddo
           endif
           do m=1,nfcoil
             if(fwtfc(m).le.0.0) cycle
             nj=nj+1
             abry(nj,nk)=0.
           enddo
           if (nbdry.gt.0.and.iconvr.eq.3) then
             do m=1,nfcoil
               nj=nj+1
               abry(nj,nk)=0.0
             enddo
           endif
         enddo
         neqn=neqn+nvesel
       endif fit_vessel
!-----------------------------------------------------------------------
!--    set up working matrix                                          --
!-----------------------------------------------------------------------
       ut(1:msbdry,1:msbdry) = 1.0
!-----------------------------------------------------------------------
!--    compute inverse of abry in least squares sence with singular   --
!--    value decomposition                                            --
!-----------------------------------------------------------------------
#ifdef DEBUG_LEVEL2
       write (106,*) 'FCURRT nj,neqn = ', nj,neqn
#endif
       call sdecm(abry,msbdry,nj,neqn,ut,msbdry,nj,wbry,work,ier)
       if (ier.eq.129) then
         kerror = 1
         call errctrl_msg('fcurrt','sdecm failed to converge')
         return
       endif
       cond=ier
       toler=1.0e-06_dp*wbry(1)
       do i=1,neqn
         do j=1,nj
           ainbry(i,j)=0.0
           do k=1,neqn
             t=0.0
             if(wbry(k).gt.toler) t=1.0/wbry(k)
             ainbry(i,j) = ainbry(i,j) + abry(i,k)*ut(k,j)*t
           enddo
         enddo
       enddo
!-----------------------------------------------------------------------
!--    compute and store coil reference currents, not need for        --
!--    non-fixed boundary when IFREF=-1                               --
!-----------------------------------------------------------------------
       if (ifref.gt.0.and.iconvr.eq.3) then
         do i=1,nj
           if (i.le.nbdry) then
             bbry(i)=fwtbry(i)
           else
             bbry(i)=0.0
           endif
         enddo
         do i=1,nfcoil
           fcref(i)=sum(ainbry(i,1:nj)*bbry(1:nj))
         enddo
       endif
      endif first_iter
!-----------------------------------------------------------------------
!--   RHS of response matrix, start here if got inverse matrix already--
!-----------------------------------------------------------------------
      nj=0
      fixed_bdry: if (nbdry.gt.0.and.iconvr.eq.3) then
!-----------------------------------------------------------------------
!--     fixed boundary portion                                        --
!-----------------------------------------------------------------------
        do m=1,nbdry
          bbry(m)=0.0
          if(ivesel.gt.0.and.ifitvs.le.0) &
            bbry(m)=bbry(m)+sum(rbdrvs(m,1:nvesel)*vcurrt(1:nvesel))
          if(iecurr.gt.0) &
            bbry(m)=bbry(m)+sum(rbdrec(m,1:nesum)*ecurrt(1:nesum))
          if(iacoil.gt.0) &
            bbry(m)=bbry(m)+sum(rbdrac(m,1:nacoil)*caccurt(jtime,1:nacoil))
          bbry(m)=bbry(m)+sum(rbdrpc(m,1:nwnh)*pcurrt(1:nwnh))
          pbry(m)=bbry(m)
          bbry(m)=fwtbry(m)*(wsibry-bbry(m))
        enddo
        nj=nbdry
#ifdef DEBUG_LEVEL2
        write (106,*) 'ivesel,iecurr,iacoil= ',ivesel,iecurr,iacoil
        write (106,*) 'pbry= ',(pbry(m),m=1,nbdry)
        write (106,*) 'bbry= ',(bbry(m),m=1,nbdry)
        write (106,*) 'pbry(24),bbry(24)= ',pbry(24),bbry(24)
        write (106,*) 'rgrid,zgrid= ',rgrid(75*nw/129),zgrid(104*nh/129)
        write (106,*) 'rbdry,zbdry= ',rbdry(24),zbdry(24)
        bdrav=abs(rbdrpc(24,1))
        bdrmax=abs(rbdrpc(24,1))
        mbdrmax=1
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
        zdif=zbdry(24)-zgrid(104*nh/129)
        bdrmax=psical(rbdry(24),rgrid(75*nw/129),zdif)*tmu
        ! TODO: this goes out of bounds
!        bdrmaxs=gridpc(9547,75*nw/129)
        write (106,*) 'mbdrmax,bdrmax,bdrmaxs= ',mbdrmax,bdrmax!,bdrmaxs
#endif
      endif fixed_bdry
!-----------------------------------------------------------------------
!--   flux loops portion                                              --
!-----------------------------------------------------------------------
      do m=1,nsilop
        if (fwtsi(m).le.0.0) cycle
        nj=nj+1
        bbry(nj)=sum(gsilpc(m,1:nwnh)*pcurrt(1:nwnh))
        if(iecurr.gt.0) &
          bbry(nj)=bbry(nj)+sum(rsilec(m,1:nesum)*ecurrt(1:nesum))
!-----------------------------------------------------------------------
!--     specify vessel currents ?                                     --
!-----------------------------------------------------------------------
        if(ivesel.gt.0.and.ifitvs.le.0) &
          bbry(nj)=bbry(nj)+sum(rsilvs(m,1:nvesel)*vcurrt(1:nvesel))
!-----------------------------------------------------------------------
!--     use advance divertor coil ?                                   --
!-----------------------------------------------------------------------
        if(iacoil.gt.0) &
          bbry(nj)=bbry(nj)+sum(rsilac(m,1:nacoil)*caccurt(jtime,1:nacoil))
        if (fitsiref) then
          bbry(nj)=fwtsi(m)*(silopt(jtime,m)+psiref(jtime)-bbry(nj))
        else
          bbry(nj)=fwtsi(m)*(silopt(jtime,m)-bbry(nj))
        endif
      enddo
!-----------------------------------------------------------------------
!--   magnetic probes used only for vacuum error field analysis       --
!-----------------------------------------------------------------------
      vacuum: if (ivacum.gt.0) then
        do m=1,magpri
          if(fwtmp2(m).le.0.0) cycle
          nj=nj+1
          bbry(nj)=0.0
          bbry(nj)=bbry(nj)+sum(gmp2pc(m,1:nwnh)*pcurrt(1:nwnh))
          if(iecurr.gt.0) &
            bbry(nj)=bbry(nj)+sum(rmp2ec(m,1:nesum)*ecurrt(1:nesum))
!-----------------------------------------------------------------------
!--       specify vessel currents ?                                   --
!-----------------------------------------------------------------------
          if(ivesel.gt.0.and.ifitvs.le.0) &
            bbry(nj)=bbry(nj)+sum(rmp2vs(m,1:nvesel)*vcurrt(1:nvesel))
          if(iacoil.gt.0) &
            bbry(nj)=bbry(nj)+sum(rmp2ac(m,1:nacoil)*caccurt(jtime,1:nacoil))
          bbry(nj)=fwtmp2(m)*(expmpi(jtime,m)-bbry(nj))
        enddo
      endif vacuum
!-----------------------------------------------------------------------
!--   F-coil currents specification                                   --
!-----------------------------------------------------------------------
      do m=1,nfcoil
        if(fwtfc(m).le.0.0) cycle
        nj=nj+1
        bbry(nj)=fccurt(jtime,m)*fwtfc(m)
      enddo
!-----------------------------------------------------------------------
!--   F coil current minimization                                     --
!-----------------------------------------------------------------------
      if (nbdry.gt.0.and.iconvr.eq.3) then
        do m=1,nfcoil
          nj=nj+1
          bbry(nj)=0.0
        enddo
!-----------------------------------------------------------------------
!--     symmetrize F coil currents ?                                  --
!-----------------------------------------------------------------------
        if ((symmetrize).and.(cupdown.ne.0.0)) then
          do m=1,nfcoil/2
            nj=nj+1
            bbry(nj)=0.0
          enddo
        endif
      endif
!-----------------------------------------------------------------------
!--   F coil constraints ?                                            --
!-----------------------------------------------------------------------
      if (kccoils.gt.0) then
        do m=1,kccoils
          nj=nj+1
          bbry(nj)=xcoils(m)
        enddo
      endif
!-----------------------------------------------------------------------
!--   SOL Constraints RHS                                             --
!-----------------------------------------------------------------------
      if (nsol.gt.0) then
#ifdef DEBUG_LEVEL2
        write (6,*) 'FCURRT wsisol = ', wsisol
#endif
        do m=1,nsol
          nj=nj+1
          bbry(nj)=0.0
          if(ivesel.gt.0.and.ifitvs.lt.0) &
            bbry(nj)=bbry(nj)+sum(rsolvs(m,1:nvesel)*vcurrt(1:nvesel))
          if(iecurr.gt.0) &
            bbry(nj)=bbry(nj)+sum(rsolec(m,1:nesum)*ecurrt(1:nesum))
          bbry(nj)=bbry(nj)+sum(rsolpc(m,1:nwnh)*pcurrt(1:nwnh))
          psbry(m)=bbry(nj)
          bbry(nj)=fwtsolw(m)*(wsisol-bbry(nj))
        enddo
      endif
!-----------------------------------------------------------------------
!--   now get F coil currents from precomputed inverse matrix         --
!-----------------------------------------------------------------------
      do i=1,nfcoil
        brsp(i)=sum(ainbry(i,1:nj)*bbry(1:nj))
      enddo
      if (ifitvs.gt.0) then
        do ii=1,nvesel
          i=ii+nfcoil
          vcurrt(ii)=sum(ainbry(i,1:nj)*bbry(1:nj))
        enddo
!
        if(nfourier.gt.1) then
          do j=1,nvesel
            temp=sum(vcurrt(1:(nfourier*2+1))*vecta(1:(nfourier*2+1),j))
            vcurrt(j)=temp
          enddo
        endif
      endif
!-----------------------------------------------------------------------
!--   adjustments for various boundary flux options                   --
!-----------------------------------------------------------------------
      bdry_flux: if (ifref.gt.0.and.iconvr.eq.3) then
        select case (ifref)
        case(1)
!-----------------------------------------------------------------------
!--       sum of inner F-coils 1-5 A and B zero IFREF=1               --
!-----------------------------------------------------------------------
          sumif = sum(brsp(1:5) + brsp(10:14))
          sumifr = sum(fcref(1:5) + fcref(10:14))
          sumif = sumif + brsp(8) + brsp(17)
          sumifr = sumifr + fcref(8) + fcref(17)
!-----------------------------------------------------------------------
!--       sum of F coils selected through FCSUM vanish IFREF=2        --
!-----------------------------------------------------------------------
        case(2)
          sumif = sum(fcsum(1:nfcoil)*brsp(1:nfcoil)/turnfc(1:nfcoil))
          sumifr = sum(fcref(1:nfcoil)*fcsum(1:nfcoil)/turnfc(1:nfcoil))
        case(3)
!----------------------------------------------------------------------
!--     choose boundary flux by minimize coil currents IFREF=3       --
!----------------------------------------------------------------------
          sumif = sum(fcref(1:nfcoil)*brsp(1:nfcoil)*fczero(1:nfcoil))
          sumifr = sum(fcref(1:nfcoil)**2*fczero(1:nfcoil))
        end select
!-----------------------------------------------------------------------
!--     update boundary flux for IFREF=1-3                            --
!-----------------------------------------------------------------------
        if (ifref.le.3) then
          ssiref = sumif/sumifr
          silopt(jtime,1:nfcoil)=silopt(jtime,1:nfcoil)-ssiref
          brsp(1:nfcoil) = brsp(1:nfcoil) - ssiref*fcref(1:nfcoil)
          wsibry=wsibry-ssiref
          wsisol=wsisol-ssiref
        endif
!------------------------------------------------------------------------
!--     fixed boundary flux specified through PSIBRY, IFREF=4          --
!------------------------------------------------------------------------
        if (ifref.eq.4) then
          ssiref=psibry0-psibry
          brsp(1:nfcoil)=brsp(1:nfcoil)+ssiref*fcref(1:nfcoil)
          silopt(jtime,1:nfcoil)=silopt(jtime,1:nfcoil)+ssiref
          wsibry=wsibry+ssiref
          wsisol=wsisol+ssiref
        endif

#ifdef DEBUG_LEVEL2
        write (106,*) '      PSIBRY0,PSIBRY,WSIBRY= ',&
                      psibry0,psibry,wsibry
#endif
      endif bdry_flux
!-----------------------------------------------------------------------
!--   estimated errors for fixed boundary calculations                --
!-----------------------------------------------------------------------
      if (nbdry.gt.0.and.iconvr.eq.3) then
      erbmax=0.0
      erbave=0.0
      do i=1,nbdry
        xsibry=pbry(i)+sum(rbdrfc(i,1:nfcoil)*brsp(1:nfcoil))
        erbloc(i)=abs((wsibry-xsibry)/sidif)
        erbmax=max(erbloc(i),erbmax)
        erbave=erbave+erbloc(i)
      enddo
      erbave=erbave/nbdry
#ifdef DEBUG_LEVEL1
      write (6,*) 'FCURRT erbmax,erbave,si = ',erbmax,erbave,wsibry
#endif
#ifdef DEBUG_LEVEL2
      write (106,*) 'XSIBRY,PBRY(1),BRSP= ', &
                    xsibry,pbry(1),(brsp(m),m=1,nfcoil)
#endif
      endif
!
      if (nsol.gt.0.and.iconvr.eq.3) then
        erbsmax=0.0
        erbsave=0.0
        do i=1,nsol
          xsisol=psbry(i)+sum(rsolfc(i,1:nfcoil)*brsp(1:nfcoil))
          erbsloc(i)=abs((wsisol-xsisol)/sidif)
          erbsmax=max(erbsloc(i),erbsmax)
          erbsave=erbsave+erbsloc(i)
        enddo
        erbsave=erbsave/nsol
#ifdef DEBUG_LEVEL1
        write (6,*) 'FCURRT erbsmax,erbsave,si = ',erbsmax,erbsave,wsisol
#endif
      endif
!
      return
      end subroutine fcurrt

!**********************************************************************
!>
!!    vescur computes the currents induced in the vessel
!!    segments due to E coils and F coils
!!
!!
!!    @param jtime : time index
!!
!**********************************************************************
      subroutine vescur(jtime)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      integer*4, intent(in) :: jtime
!
      if(ivesel.eq.5) return
!
      if(ivesel.eq.1) return
!
      vcurrt(1:nvesel)=vloopt(jtime)/rsisvs(1:nvesel)
      return
      end subroutine vescur
