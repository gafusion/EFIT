#include "config.f"
!**********************************************************************
!>
!!    fcurrt computes the currents in the f coils.
!!
!!    @param jtime :
!!
!!    @param iter :
!!
!!    @param itertt :
!!
!!    @param kerror :
!!
!**********************************************************************
      subroutine fcurrt(jtime,iter,itertt,kerror)
      use commonblocks,only: rfcpc
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      dimension afma(nfcoil,nfcoil),ifmatr(nfcoil),wfmatr(nfcoil), &
                wbry(msbdry),work(msbdr2),ut(msbdry,msbdry)
      dimension abry(msbdry,nfcoil+nvesel),bbry(msbdry), &
                ainbry(nfcoil+nvesel,msbdry)
      dimension fcref(nfcoil)
      dimension pbry(msbdry),psbry(msbdry)
      integer*4, intent(inout) :: kerror
      data iskip/0/,idoit/0/

      kerror = 0
      if (ifcurr.gt.0) return
      if (itertt.le.1.and.icinit.lt.0) return
      if (islpfc.eq.1) then
!-----------------------------------------------------------------------
!--     flux loop on F coils, seldom used for DIII-D                  --
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
            do i=1,nfcoil
              do j=1,nfcoil
                afma(i,j)=rfcfc(i,j)
              enddo
            enddo
            m111=-1.0
            call decomp(nfcoil,nfcoil,afma,m111,ifmatr,wfmatr)
!vas
!vas      print*,'file name : ','fm'//trim(ch1)// &
!vas                         trim(ch2)//'.ddd'
            open(unit=nffile,status='old',form='unformatted', &
                 file='fm'//trim(ch1)//trim(ch2)//'.ddd',iostat=ioerr)
            if (ioerr.eq.0) close(unit=nffile,status='delete')
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
          if (ivacum.le.0) then
            do j=1,nwnh
              brsp(i)=brsp(i)+rfcpc(i,j)*pcurrt(j)
            enddo
          endif
          if (iecurr.gt.0) then
            do j=1,nesum
              brsp(i)=brsp(i)+rfcec(i,j)*ecurrt(j)
            enddo
          endif
          if (ivesel.gt.0.and.ifitvs.le.0) then
            do j=1,nvesel
              brsp(i)=brsp(i)+rfcvs(i,j)*vcurrt(j)
            enddo
          endif
          brsp(i)=csilop(i,jtime)-brsp(i)
          if (fitsiref) brsp(i)=brsp(i)+psiref(jtime)
        enddo
        call solve(nfcoil,nfcoil,afma,brsp,ifmatr)
        return
!-----------------------------------------------------------------------
!--     standard option, flux loops away from F coils                 --
!-----------------------------------------------------------------------
      endif
      if (idoit.le.0.or.itertt.le.1) then
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
         do i=1,nbdry
           fwtbry(i)=fwtbdr*fwtbdry(i)
         enddo
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
           do m=1,nbdry
             abry(m,nk)=fwtbry(m)*rbdrfc(m,nk)
           enddo
           nj=nbdry
         endif
!-----------------------------------------------------------------------
!--      flux loop option                                             --
!-----------------------------------------------------------------------
         do m=1,nsilop
           if (fwtsi(m).le.0.0) cycle
           nj=nj+1
           abry(nj,nk)=fwtsi(m)*rsilfc(m,nk)
         enddo
!-----------------------------------------------------------------------
!--      magnetic probes used only for vacuum error field analysis    --
!-----------------------------------------------------------------------
         if (ivacum.gt.0) then
           do m=1,magpri
             if (fwtmp2(m).le.0.0) cycle
             nj=nj+1
             abry(nj,nk)=fwtmp2(m)*rmp2fc(m,nk)
           enddo
         endif
!-----------------------------------------------------------------------
!--      F coil currents option                                       --
!-----------------------------------------------------------------------
         do m=1,nfcoil
           if (fwtfc(m).le.0.0) cycle
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
             if (nk.eq.m) abry(nj,nk)=cfcoil*fczero(m)
           enddo
!-----------------------------------------------------------------------
!--        constraint F coil currents to be symmetric if requested    --
!-----------------------------------------------------------------------
           if ((symmetrize).and.(cupdown.ne.0.0)) then
             do m=1,nfcoil/2
               nj=nj+1
               abry(nj,nk)=0.0
               if (nk.eq.m) abry(nj,nk)=cupdown
               if (nk.eq.(m+nfcoil/2)) abry(nj,nk)=-cupdown
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
!--      SOL constraints 2014/05/20 LL                                --
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
       if (ifitvs.gt.0) then
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
               if(nfourier.gt.1) then
                 temp=0.
                 do i=1,nvesel
                   temp=temp+rbdrvs(m,i)*vecta(nkk,i)
                 enddo
                 abry(m,nk)=fwtbry(m)*temp
               else
                 abry(m,nk)=fwtbry(m)*rbdrvs(m,nkk)
               endif
             enddo
             nj=nbdry
           endif
           do m=1,nsilop
             if (fwtsi(m).le.0.0) cycle
             nj=nj+1
             if(nfourier.gt.1) then
               temp=0.
               do i=1,nvesel
                 temp=temp+rsilvs(m,i)*vecta(nkk,i)
               enddo
               abry(nj,nk)=fwtsi(m)*temp
             else
               abry(nj,nk)=fwtsi(m)*rsilvs(m,nkk)
             endif
           enddo
!-----------------------------------------------------------------------
!--        magnetic probes used only for vacuum error analysis        --
!-----------------------------------------------------------------------
           if (ivacum.gt.0) then
             do m=1,magpri
               if (fwtmp2(m).le.0.0) cycle
               nj=nj+1
               if (nfourier.gt.1) then
                 temp=0.
                 do i=1,nvesel
                   temp=temp+rmp2vs(m,i)*vecta(nkk,i)
                 enddo
                 abry(nj,nk)=fwtmp2(m)*temp
               else
                 abry(nj,nk)=fwtmp2(m)*rmp2vs(m,nkk)
               endif
             enddo
           endif
           do m=1,nfcoil
             if (fwtfc(m).le.0.0) cycle
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
       endif
!-----------------------------------------------------------------------
!--    set up working matrix                                          --
!-----------------------------------------------------------------------
       do i=1,msbdry
         do j=1,msbdry
           ut(i,j)=0.0
         enddo
       enddo
       do i=1,msbdry
         ut(i,i) = 1.0
       enddo
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
       end if
       cond=ier
       toler=1.0e-06_dp*wbry(1)
       do i=1,neqn
         do j=1,nj
           ainbry(i,j)=0.0
           do k=1,neqn
             t=0.0
             if (wbry(k).gt.toler) t=1.0/wbry(k)
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
           fcref(i)=0.0
           do j=1,nj
             fcref(i)=fcref(i)+ainbry(i,j)*bbry(j)
           enddo
         enddo
       endif
      endif
!-----------------------------------------------------------------------
!--   RHS of response matrix, start here if got inverse matrix already--
!-----------------------------------------------------------------------
      nj=0
      if (nbdry.gt.0.and.iconvr.eq.3) then
!-----------------------------------------------------------------------
!--     fixed boundary portion                                        --
!-----------------------------------------------------------------------
        do m=1,nbdry
          bbry(m)=0.0
          if (ivesel.gt.0.and.ifitvs.le.0) then
            do k=1,nvesel
              bbry(m)=bbry(m)+rbdrvs(m,k)*vcurrt(k)
            enddo
          endif
          if (iecurr.gt.0) then
            do k=1,nesum
              bbry(m)=bbry(m)+rbdrec(m,k)*ecurrt(k)
            enddo
          endif
          if (iacoil.gt.0) then
            do k=1,nacoil
              bbry(m)=bbry(m)+rbdrac(m,k)*caccurt(jtime,k)
            enddo
          endif
          do k=1,nwnh
            bbry(m)=bbry(m)+rbdrpc(m,k)*pcurrt(k)
          enddo
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
      endif
!-----------------------------------------------------------------------
!--   flux loops portion                                              --
!-----------------------------------------------------------------------
      do m=1,nsilop
        if (fwtsi(m).le.0.0) cycle
        nj=nj+1
        bbry(nj)=0.0
        do j=1,nwnh
          bbry(nj)=bbry(nj)+gsilpc(m,j)*pcurrt(j)
        enddo
        if (iecurr.gt.0) then
          do j=1,nesum
            bbry(nj)=bbry(nj)+rsilec(m,j)*ecurrt(j)
          enddo
        endif
!-----------------------------------------------------------------------
!--     specify vessel currents ?                                     --
!-----------------------------------------------------------------------
        if (ivesel.gt.0.and.ifitvs.le.0) then
          do j=1,nvesel
            bbry(nj)=bbry(nj)+rsilvs(m,j)*vcurrt(j)
          enddo
        endif
!-----------------------------------------------------------------------
!--     use advance divertor coil ?                                   --
!-----------------------------------------------------------------------
        if (iacoil.gt.0) then
          do j=1,nacoil
            bbry(nj)=bbry(nj)+rsilac(m,j)*caccurt(jtime,j)
          enddo
        endif
        if (fitsiref) then
          bbry(nj)=fwtsi(m)*(silopt(jtime,m)+psiref(jtime)-bbry(nj))
        else
          bbry(nj)=fwtsi(m)*(silopt(jtime,m)-bbry(nj))
        endif
      enddo
!-----------------------------------------------------------------------
!--   magnetic probes used only for vacuum error field analysis       --
!-----------------------------------------------------------------------
      if (ivacum.gt.0) then
        do m=1,magpri
          if (fwtmp2(m).le.0.0) cycle
          nj=nj+1
          bbry(nj)=0.0
          do j=1,nwnh
            bbry(nj)=bbry(nj)+gmp2pc(m,j)*pcurrt(j)
          enddo
          if (iecurr.gt.0) then
            do j=1,nesum
              bbry(nj)=bbry(nj)+rmp2ec(m,j)*ecurrt(j)
            enddo
          endif
!-----------------------------------------------------------------------
!--       specify vessel currents ?                                   --
!-----------------------------------------------------------------------
          if (ivesel.gt.0.and.ifitvs.le.0) then
            do j=1,nvesel
              bbry(nj)=bbry(nj)+rmp2vs(m,j)*vcurrt(j)
            enddo
          endif
          if (iacoil.gt.0) then
            do j=1,nacoil
              bbry(nj)=bbry(nj)+rmp2ac(m,j)*caccurt(jtime,j)
            enddo
          endif
          bbry(nj)=fwtmp2(m)*(expmpi(jtime,m)-bbry(nj))
        enddo
      endif
!-----------------------------------------------------------------------
!--   F-coil currents specification                                   --
!-----------------------------------------------------------------------
      do m=1,nfcoil
        if (fwtfc(m).le.0.0) cycle
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
!--   now get F coil currents from precomputed inverse matrix         --
!-----------------------------------------------------------------------
      do i=1,nfcoil
        brsp(i)=0.0
        do j=1,nj
          brsp(i)=brsp(i)+ainbry(i,j)*bbry(j)
        enddo
      enddo
      if (ifitvs.gt.0) then
        do ii=1,nvesel
          i=ii+nfcoil
          vcurrt(ii)=0.0
          do j=1,nj
            vcurrt(ii)=vcurrt(ii)+ainbry(i,j)*bbry(j)
          enddo
        enddo
!
        if(nfourier.gt.1) then
          do j=1,nvesel
            temp=0.
            do i=1,(nfourier*2+1)
              temp=temp+vcurrt(i)*vecta(i,j)
            enddo
            vcurrt(j)=temp
          enddo
        endif
      endif
!-----------------------------------------------------------------------
!--   adjustments for various boundary flux options                   --
!-----------------------------------------------------------------------
      if (ifref.gt.0.and.iconvr.eq.3) then
        sumif = 0.0
        sumifr = 0.0
!-----------------------------------------------------------------------
!--     sum of inner F-coils 1-5 A and B zero IFREF=1                 --
!-----------------------------------------------------------------------
        select case (ifref)
        case(1)
          do i=1,5
            sumif = sumif + brsp(i) + brsp(i+9)
            sumifr = sumifr + fcref(i) + fcref(i+9)
          enddo
          sumif = sumif + brsp(8) + brsp(17)
          sumifr = sumifr + fcref(8) + fcref(17)
!-----------------------------------------------------------------------
!--     sum of F coils selected through FCSUM vanish IFREF=2          --
!-----------------------------------------------------------------------
        case(2)
          do i=1,nfcoil
            sumif = sumif + fcsum(i)*brsp(i)/turnfc(i)
            sumifr = sumifr + fcref(i)*fcsum(i)/turnfc(i)
          enddo
!----------------------------------------------------------------------
!--     choose boundary flux by minimize coil currents IFREF=3       --
!----------------------------------------------------------------------
        case(3)
          do i=1,nfcoil
            sumif = sumif + fcref(i)*brsp(i)*fczero(i)
            sumifr = sumifr + fcref(i)**2*fczero(i)
          enddo
        endselect
!-----------------------------------------------------------------------
!--     update boundary flux for IFREF=1-3                            --
!-----------------------------------------------------------------------
        if (ifref.le.3) then
          ssiref = sumif/sumifr
          do m=1,nfcoil
            silopt(jtime,m)=silopt(jtime,m)-ssiref
            brsp(m) = brsp(m) - ssiref*fcref(m)
          enddo
          wsibry=wsibry-ssiref
          wsisol=wsisol-ssiref
        endif
!------------------------------------------------------------------------
!--     fixed boundary flux specified through PSIBRY, IFREF=4          --
!------------------------------------------------------------------------
        if (ifref.eq.4) then
          ssiref=psibry0-psibry
          do m=1,nfcoil
            brsp(m)=brsp(m)+ssiref*fcref(m)
            silopt(jtime,m)=silopt(jtime,m)+ssiref
          enddo
          wsibry=wsibry+ssiref
          wsisol=wsisol+ssiref
        endif

#ifdef DEBUG_LEVEL2
        write (106,*) '      PSIBRY0,PSIBRY,WSIBRY= ',&
                      psibry0,psibry,wsibry
#endif
!-----------------------------------------------------------------------
!--     done, estimated errors for fixed boundary calculations        --
!-----------------------------------------------------------------------
      endif
      if (nbdry.gt.0.and.iconvr.eq.3) then
      erbmax=0.0
      erbave=0.0
      do i=1,nbdry
        xsibry=pbry(i)
        do m=1,nfcoil
          xsibry=xsibry+rbdrfc(i,m)*brsp(m)
        enddo
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
          xsisol=psbry(i)
          do m=1,nfcoil
            xsisol=xsisol+rsolfc(i,m)*brsp(m)
          enddo
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
