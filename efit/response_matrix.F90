#include "config.f"
!**********************************************************************
!>
!!    calculate the appropriate response matrix and
!!    invert it to get the plasma current strengths.\n
!!
!!    note that npcurn=nffcur+nppcur, 
!!    nrsmat=nfsum+npcurn+number of constraints.
!!
!!    @param jtime : time index
!!
!!    @param iter : current profile (outer) loop iteration index
!!
!!    @param ichisq : chisq flag
!!
!!    @param nniter : total iteration index (current+equilibirum loops)
!!
!!    @param kerror : error flag
!!
!**********************************************************************
      subroutine response_matrix(jtime,iter,ichisq,nniter,kerror)
      use commonblocks,only: c,wk,bkx,bky
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      real*8 pwcurr,prcurr

      integer*4, intent(in) :: jtime,iter,nniter
      integer*4, intent(out) :: ichisq,kerror
      integer*4 i,j,kk,m,mk,mp1,n,ne,nj,nk,nkk,nmw,ier,info,mcentral, &
                ncrsp,need,needs,nfedge,nfffff,nkb,nload,npedge,nsavdz,nsq
      real*8 brspmin,cdelznow,cm,cmbr,cmbz,cmv,drgam,erbot,erup,ework, &
             fwtbdr,fwtbpp,fwtxxa,fwtxxf,fwtxxo,fwtxxzj,pres0,prew0, &
             ptop0,pwop0,pwp0r2,saiold,saisq,scadelz,siedge,t,toler,xjj, &
             ysiwant
      character(128) tmpstr
      integer*4, dimension(mfnpcr)     :: ipvttmp
      real*8, dimension(2)             :: arspdet2
      real*8, dimension(mfnpcr)        :: wrsp,worktmp
      real*8, dimension(nrsmat,mfnpcr) :: arsp,arsptmp
      real*8, dimension(nrsmat)        :: b,brsold
      real*8, dimension(kzeroj)        :: rxxx,rxxxf,rxx2,rxxw
      real*8 work(nrsmat*2),vcurrto(nvesel)
      real*8 xpsfp(nffcur),xpspp(nppcur),xpspwp(nwwcur)
      real*8 crsp(4*(npcurn-2)+6+npcurn*npcurn,nrsmat)
      real*8 z(4*(npcurn-2)+6+npcurn*npcurn)
      real*8 pds(6)
      integer*4, parameter :: nnn=1,minite=8 ! TODO: should this be fixed?
      real*8, parameter :: ten24=1.e4_dp,z04=1.0e-04_dp

      kerror = 0
      csilopv = 0.0
      cmpr2v = 0.0
      arsp = 0.0
      crsp = 0.0
      z = 0.0
      if(iconvr.eq.3) return
!----------------------------------------------------------------------
!--   Variable fitdelz                                               --
!----------------------------------------------------------------------
      if(fitdelz) scadelz=scaledz
!----------------------------------------------------------------------
!--   set up fitting weight for boundary constraints                 --
!----------------------------------------------------------------------
      if (nbdry.gt.0) then
        fwtbdr=abs(errbry)*max(abs(sidif),z04)
        fwtbdr=1.0/fwtbdr
        do i=1,nbdry
          if (sigrbd(i).lt.1.e10_dp.and.sigzbd(i).lt.1.e10_dp) then
            call seva2d(bkx,lkx,bky,lky,c,rbdry(i),zbdry(i),pds,ier,n666)
            fwtbdr=sqrt((sigrbd(i)*pds(2))**2+(sigzbd(i)*pds(3))**2)
            fwtbdr=1.0/fwtbdr
          endif
          fwtbry(i)=fwtbdr*fwtbdry(i)
        enddo
      endif
!----------------------------------------------------------------------
!--   set up the response matrix arsp                                --
!----------------------------------------------------------------------
      ichisq=0
      nsq=2
      saiold=chisq(jtime)
      brsold=brsp
      if(ivesel.eq.3) vcurrto=vcurrt
!----------------------------------------------------------------------
!--   singular decomposition, first F-coil currents, set up arsp     --
!----------------------------------------------------------------------
      do nk=1,nfsum
        nj=0
        do m=1,nsilop
          if(fwtsi(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtsi(m)*rsilfc(m,nk)
        enddo
        do m=1,magpri
          if(fwtmp2(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtmp2(m)*rmp2fc(m,nk)
        enddo
        do m=1,nstark
          if(fwtgam(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtgam(m)*rgamfc(m,nk)
        enddo
        do m=1,nmsels
          if(fwtbmselt(jtime,m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
        do m=1,nmsels
          if(fwtemselt(jtime,m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
        do m=1,nece
          if(fwtece(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtece(m)*recefc(m,nk)
        enddo
        if (fwtecebz.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=fwtecebz*recebzfc(nk)
        endif
        if (fwtcur.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        if (fwtqa.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        if (fwtbp.gt.0.0) then
          do m=2,kffcur
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (fwtdlc.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        do m=1,nfsum
          if(fwtfc(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.
          if (m.eq.nk) arsp(nj,nk)=fwtfc(m)
        enddo
!--------------------------------------------------------------------------
!--     pressure data                                                    --
!--------------------------------------------------------------------------
        if (kprfit.gt.0.and.kdofit.ne.0.and.npress.gt.0) then
          do m=1,npress
            if(fwtpre(m).le.0.0) cycle
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
!--------------------------------------------------------------------
!--     P'(1)                                                      --
!--------------------------------------------------------------------
        if (kpressb.eq.2) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
!---------------------------------------------------------------------------
!--    rotational pressure data                                           --
!---------------------------------------------------------------------------
        if (kprfit.ge.3.and.npresw.gt.0) then
          do m=1,npresw
            if (fwtprw(m).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=0.0
            endif
          enddo
        endif
!---------------------------------------------------------------------------
!--     J(PSIWANT) constraint                                             --
!---------------------------------------------------------------------------
        if (kzeroj.gt.0) then
          do i=1,kzeroj
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
!----------------------------------------------------------------------------
!--     P', FF', and rotational constraints                                --
!----------------------------------------------------------------------------
        if (kcalpa.gt.0) then
          do j=1,kcalpa
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (kcgama.gt.0) then
          do j=1,kcgama
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (kcomega.gt.0) then
          do j=1,kcomega
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
!------------------------------------------------------------------------------
!--     Boundary constraints, data                                           --
!------------------------------------------------------------------------------
        if (nbdry.gt.0) then
          do j=1,nbdry
            if (fwtbdry(j).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=fwtbry(j)*rbdrfc(j,nk)
            endif
          enddo
        endif
!------------------------------------------------------------------------------
!--     E coil currents, data                                                --
!------------------------------------------------------------------------------
        if (iecurr.eq.2) then
          do m=1,nesum
            if (fwtec(m).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=0.
            endif
          enddo
        endif
!------------------------------------------------------------------------------
!--     Reference flux, data                                                 --
!------------------------------------------------------------------------------
        if (fitsiref) then
          if (fwtref.gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=0.0
          endif
        endif
!------------------------------------------------------------------------
!--     Summation of F-coils currents
!------------------------------------------------------------------------
        if (fitfcsum) then
          nj=nj+1
          arsp(nj,nk)=fwtfcsum(nk)
        endif
      enddo
!-----------------------------------------------------------------------
!--  plasma current P', FF', and Pw', set up response matrix arsp     --
!-----------------------------------------------------------------------
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        ysiwant=sizeroj(i)
        call setfp(ysiwant,xpsfp)
        call setpp(ysiwant,xpspp)
!----------------------------------------------------------------------
!--  local or flux surface averaged J constraint ?                   --
!----------------------------------------------------------------------
        if (rzeroj(i).gt.0.0) rxxx(i)=rzeroj(i)
        if (rzeroj(i).lt.0.0) rxxx(i)=rseps(1,jtime)/100.
        rxxxf(i)=rxxx(i)
        if (rzeroj(i).eq.0.0) then
          rxxx(i)=1./r1sdry(i)
          rxxxf(i)=r1sdry(i)/r2sdry(i)
        endif
        if (rxxx(i).le.0.0) then
          rxxx(i)=rcentr
          rxxxf(i)=rxxx(i)
        endif
        brspmin=max(ten24,abs(brsp(nfsum+1)))
        fwtxxzj=fwtxxj*1000./brspmin
        rotation: if (kvtor.gt.0) then
          call setpwp(ysiwant,xpspwp)
          if (nniter.gt.0) then
            if (kvtor.eq.2.or.kvtor.eq.3) then
              prew0=pwcurr(ysiwant,kwwcur)
              pres0=prcurr(ysiwant,kppcur)
              if (abs(pres0).gt.1.e-10_dp) then
                pwop0=prew0/pres0
              else
                pwop0=0.0
              endif
            endif
          endif
          if (rzeroj(i).gt.0.0) then
            rxx2(i)=(rzeroj(i)/rvtor)**2-1.
            rxxw(i)=rxx2(i)*rzeroj(i)
            if (nniter.gt.0) then
              if (kvtor.eq.2) then
                rxxw(i)=rxxw(i)*(1.+pwop0*rxx2(i))
                rxxx(i)=rxxx(i)*(1.-0.5_dp*(pwop0*rxx2(i))**2)
              endif
              if (kvtor.eq.3) then
                pwp0r2=pwop0*rxx2(i)
                ptop0=exp(pwp0r2)
                rxxw(i)=rxxw(i)*ptop0
                rxxx(i)=rxxx(i)*ptop0*(1.-pwp0r2)
              endif
            endif
          endif
          if (rzeroj(i).lt.0.0) then
            rxxw(i)=rseps(1,jtime)/100.
            rxx2(i)=(rxxw(i)/rvtor)**2-1.
            rxxw(i)=rxx2(i)*rxxw(i)
            if (nniter.gt.0) then
              if (kvtor.eq.2) then
                rxxw(i)=rxxw(i)*(1.+pwop0*rxx2(i))
                rxxx(i)=rxxx(i)*(1.-0.5_dp*(pwop0*rxx2(i))**2)
              endif
              if (kvtor.eq.3) then
                pwp0r2=pwop0*rxx2(i)
                ptop0=exp(pwp0r2)
                rxxw(i)=rxxw(i)*ptop0
                rxxx(i)=rxxx(i)*ptop0*(1.-pwp0r2)
              endif
            endif
          endif
          if (rzeroj(i).eq.0.0) then
            rxx2(i)=r2wdry/rvtor**2-1.
            rxxw(i)=rxx2(i)/r1sdry(i)
            if (nniter.gt.0) then
              if (kvtor.eq.2) then
                rxxw(i)=rxxw(i)+pwop0*r4wdry/r1sdry(i)
                rxxx(i)=rxxx(i)-0.5_dp*pwop0**2*r4wdry/r1sdry(i)
              endif
              if (kvtor.eq.3) then
                rxxx(i)=(rpwdry-pwop0*rp2wdry)/r1sdry(i)
                rxxw(i)=rp2wdry/r1sdry(i)
              endif
            endif
          endif
        endif rotation
       enddo
      endif
!----------------------------------------------------------------------
!--   start loop for plasma fitting parameters: P', FF', Pw'         --
!----------------------------------------------------------------------
      do nk=nfsum+1,nfnwcr
        n=nk-nfsum
        nj=0
        do m=1,nsilop
          if(fwtsi(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtsi(m)*rsilpc(m,n)
        enddo
        do m=1,magpri
          if(fwtmp2(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtmp2(m)*rmp2pc(m,n)
        enddo
        do m=1,nstark
          if(fwtgam(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtgam(m)*rgampc(m,n)
        enddo
        do m=1,nmsels
          if(fwtbmselt(jtime,m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtbmselt(jtime,m)*rmlspc(m,n)
        enddo
#ifdef DEBUG_MSELS
        m=3
        write (6,*) 'MATRIX jtime, m,n,fwtbmselt,rmlspc = ',jtime, &
             m,n,fwtbmselt(jtime,m),rmlspc(m,n)
        write (6,*) 'nj,nk,arsp= ',nj,nk,arsp(nj,nk)
#endif
        do m=1,nmsels
          if(fwtemselt(jtime,m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
        do m=1,nece
          if(fwtece(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtece(m)*recepc(m,n)
        enddo
        if (fwtecebz.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=fwtecebz*recebzpc(n)
        endif
        if (fwtcur.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=fwtcur*fgowpc(n)
        endif
        if (fwtqa.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        if (fwtbp.gt.0.0) then
          do m=2,kffcur
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (fwtdlc.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        do m=1,nfsum
          if(fwtfc(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.
        enddo
!--------------------------------------------------------------------
!--     pressure                                                   --
!--------------------------------------------------------------------
        if (kprfit.gt.0.and.kdofit.ne.0.and.npress.gt.0) then
          do m=1,npress
            if(fwtpre(m).le.0.0) cycle
            nj=nj+1
            arsp(nj,nk)=0.0
            if (n.le.kppcur) arsp(nj,nk)=rprepc(m,n)/sigpre(m)*fwtpre(m)
          enddo
        endif
!--------------------------------------------------------------------
!--     P'(1)                                                      --
!--------------------------------------------------------------------
        if (kpressb.eq.2) then
          nj=nj+1
          if (n.le.kppcur) arsp(nj,nk)=1./sigppb/darea
        endif
!---------------------------------------------------------------------
!--     rotational pressure                                         --
!---------------------------------------------------------------------
        if (kprfit.ge.3.and.npresw.gt.0) then
          do m=1,npresw
            if (fwtprw(m).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=0.0
              if (n.gt.kpcurn) then
                nmw=n-kpcurn
                arsp(nj,nk)=rprwpc(m,nmw)/sigprw(m)*fwtprw(m)
              endif
            endif
          enddo
        endif
!----------------------------------------------------------------------
!--     J(PSIWANT) constraint                                        --
!----------------------------------------------------------------------
        if (kzeroj.gt.0) then
          do i=1,kzeroj
            nj=nj+1
            ysiwant=sizeroj(i)
            if (n.le.kppcur) then
              call setpp(ysiwant,xpspp)
              xjj=xpspp(n)
              arsp(nj,nk)=rxxx(i)*fwtxxzj*xjj
            elseif (n.le.kpcurn) then
              call setfp(ysiwant,xpsfp)
              xjj=xpsfp(n-kppcur)
              arsp(nj,nk)=fwtxxzj/rxxxf(i)*xjj
            elseif (kvtor.gt.0) then
              xjj=xpspwp(n-kpcurn)
              arsp(nj,nk)=rxxw(i)*fwtxxzj*xjj
            endif
          enddo
        endif
!-------------------------------------------------------------------------
!--     p' and ff' constraints                                          --
!-------------------------------------------------------------------------
        if (kcalpa.gt.0) then
          brspmin=max(ten24,abs(brsp(nfsum+1)))
          fwtxxa=fwtxx*1000./brspmin
          do j=1,kcalpa
            nj=nj+1
            if (n.le.kppcur) then
              arsp(nj,nk)=calpa(n,j)*fwtxxa
            else
              arsp(nj,nk)=0.0
            endif
          enddo
        endif
        if (kcgama.gt.0) then
          brspmin=max(ten24,abs(brsp(nfsum+1)))
          fwtxxf=fwtxx*1000./brspmin
          do j=1,kcgama
            nj=nj+1
            if (n.le.kppcur) then
              arsp(nj,nk)=0.0
            elseif (n.le.kpcurn) then
              arsp(nj,nk)=cgama(n-kppcur,j)*fwtxxf
            else
              arsp(nj,nk)=0.0
            endif
          enddo
        endif
!---------------------------------------------------------------------
!--     rotational constraints                                      --
!---------------------------------------------------------------------
        if (kcomega.gt.0) then
          brspmin=max(ten24,abs(brsp(nfsum+1)))
          fwtxxo=fwtxx*1000./brspmin
          do j=1,kcomega
            nj=nj+1
            if (n.le.kpcurn) then
              arsp(nj,nk)=0.0
            else
              arsp(nj,nk)=comega(n,j)*fwtxxo
            endif
          enddo
        endif
!------------------------------------------------------------------------------
!--     Boundary constraints due to plasma contributions                     --
!------------------------------------------------------------------------------
        if (nbdry.gt.0) then
          do j=1,nbdry
            if (fwtbdry(j).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=fwtbry(j)*gbdrpc(j,n)
            endif
          enddo
        endif
!------------------------------------------------------------------------------
!--     E coil currents                                                      --
!------------------------------------------------------------------------------
        if (iecurr.eq.2) then
          do m=1,nesum
            if (fwtec(m).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=0.
            endif
          enddo
        endif
!------------------------------------------------------------------------------
!--     Reference flux                                                       --
!------------------------------------------------------------------------------
        if (fitsiref) then
          if (fwtref.gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=0.0
          endif
        endif
!-------------------------------------------------------------------------
!--     Summation of F-coils currents
!-------------------------------------------------------------------------
        if (fitfcsum) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
      enddo
      need=nfnwcr
!----------------------------------------------------------------------
!--   fit vessel currents                                            --
!----------------------------------------------------------------------
      fit_vessel: if (ivesel.eq.3) then
      if (nfourier.gt.1) then
        need=need+nfourier*2+1
      else
        need=need+nvesel
      endif
      do nk=nfnwcr+1,need
        mk=nk-nfnwcr
        nj=0
        do m=1,nsilop
          if(fwtsi(m).le.0.0) cycle
          nj=nj+1
          if (nfourier.gt.1) then
            arsp(nj,nk)=fwtsi(m)*sum(rsilvs(m,:)*vecta(mk,:))
          else
            arsp(nj,nk)=fwtsi(m)*rsilvs(m,mk)
          endif
        enddo
        do m=1,magpri
          if(fwtmp2(m).le.0.0) cycle
          nj=nj+1
          if (nfourier.gt.1) then
            arsp(nj,nk)=fwtmp2(m)*sum(rmp2vs(m,:)*vecta(mk,:))
          else
            arsp(nj,nk)=fwtmp2(m)*rmp2vs(m,mk)
          endif
        enddo
        do m=1,nstark
          if(fwtgam(m).le.0.0) cycle
          nj=nj+1
          if (nfourier.gt.1) then
            arsp(nj,nk)=fwtgam(m)*sum(rgamvs(m,:)*vecta(mk,:))
          else
            arsp(nj,nk)=fwtgam(m)*rgamvs(m,mk)
          endif
        enddo
        do m=1,nmsels
          if(fwtbmselt(jtime,m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
        do m=1,nmsels
          if(fwtemselt(jtime,m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
        do m=1,nece
          if(fwtece(m).le.0.0) cycle  ! add by qian for ece
          nj=nj+1
          if (nfourier.gt.1) then
            arsp(nj,nk)=fwtece(m)*sum(recevs(m,:)*vecta(mk,:))
          else
            arsp(nj,nk)=fwtece(m)*recevs(m,mk)
          endif
        enddo
        if (fwtecebz.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=fwtecebz*recebzvs(mk)
        endif
        if (fwtcur.gt.0.0) then
          nj=nj+1
          if (nfourier.gt.1) then
            arsp(nj,nk)=fwtcur*sum(vecta(mk,:))
          else
            arsp(nj,nk)=fwtcur
          endif
        endif
        if (fwtqa.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        if (fwtbp.gt.0.0) then
          do m=2,kffcur
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (fwtdlc.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        do m=1,nfsum
          if(fwtfc(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.
        enddo
!--------------------------------------------------------------------------
!--     pressure                                                         --
!--------------------------------------------------------------------------
        if (kprfit.gt.0.and.kdofit.ne.0.and.npress.gt.0) then
          do m=1,npress
            if(fwtpre(m).le.0.0) cycle
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
!--------------------------------------------------------------------
!--     P'(1)                                                      --
!--------------------------------------------------------------------
        if (kpressb.eq.2) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
!--------------------------------------------------------------------------
!--     rotational pressure                                              --
!--------------------------------------------------------------------------
        if (kprfit.ge.3.and.npresw.gt.0) then
          do m=1,npresw
            if (fwtprw(m).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=0.0
            endif
          enddo
        endif
        if (kzeroj.gt.0) then
          do i=1,kzeroj
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (kcalpa.gt.0) then
          do j=1,kcalpa
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (kcgama.gt.0) then
          do j=1,kcgama
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (kcomega.gt.0) then
          do j=1,kcomega
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
!------------------------------------------------------------------------------
!--     Boundary constraints                                                 --
!------------------------------------------------------------------------------
        if (nbdry.gt.0) then
          do j=1,nbdry
            if (fwtbdry(j).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=fwtbry(j)*rbdrvs(j,mk)
            endif
          enddo
        endif
!------------------------------------------------------------------------------
!--     E coil currents                                                      --
!------------------------------------------------------------------------------
        if (iecurr.eq.2) then
          do m=1,nesum
            if (fwtec(m).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=0.
            endif
          enddo
        endif
!------------------------------------------------------------------------------
!--     Reference flux                                                       --
!------------------------------------------------------------------------------
        if (fitsiref) then
          if (fwtref.gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=0.0
          endif
        endif
!------------------------------------------------------------------------------
!--     add by qian sum of vessel sement currents=IPV1-IP1                   --
!------------------------------------------------------------------------------
!       TODO: isumvesel appears to never be defined in efit...
!        if (isumvesel.gt.0) then
!            nj=nj+1
!            if (nfourier.gt.0) then
!              arsp(nj,nk)=fwtcur*sum(vecta(mk,:))/100.
!            else
!              arsp(nj,nk)=1.0
!            endif
!        endif
      enddo
      endif fit_vessel
!-----------------------------------------------------------------------
!--   boundary pressure term for kinetic fitting P(1)                 --
!--   P(1) is an additional fitting parameter in kinetic fitting      --
!-----------------------------------------------------------------------
      kinetic: if (kprfit.gt.0.and.kdofit.ne.0.and.npress.gt.0) then
      need=need+1
      nk=need
      nj=0
      do m=1,nsilop
        if(fwtsi(m).le.0.0) cycle
        nj=nj+1
        arsp(nj,nk)=0.0
      enddo
      do m=1,magpri
        if(fwtmp2(m).le.0.0) cycle
        nj=nj+1
        arsp(nj,nk)=0.0
      enddo
      do m=1,nstark
        if(fwtgam(m).le.0.0) cycle
        nj=nj+1
        arsp(nj,nk)=0.0
      enddo
      do m=1,nmsels
        if(fwtbmselt(jtime,m).le.0.0) cycle
        nj=nj+1
        arsp(nj,nk)=0.0
      enddo
      do m=1,nmsels
        if(fwtemselt(jtime,m).le.0.0) cycle
        nj=nj+1
        arsp(nj,nk)=0.0
      enddo
      do m=1,nece
        if(fwtece(m).le.0.0) cycle
        nj=nj+1
        arsp(nj,nk)=0.0
      enddo
      if (fwtecebz.gt.0.0) then
        nj=nj+1
        arsp(nj,nk)=0.0
      endif
      if (fwtcur.gt.0.0) then
        nj=nj+1
        arsp(nj,nk)=0.0
      endif
      if (fwtqa.gt.0.0) then
        nj=nj+1
        arsp(nj,nk)=0.0
      endif
      if (fwtbp.gt.0.0) then
        do m=2,kffcur
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
      endif
      if (fwtdlc.gt.0.0) then
        nj=nj+1
        arsp(nj,nk)=0.0
      endif
      do m=1,nfsum
        if(fwtfc(m).le.0.0) cycle
        nj=nj+1
        arsp(nj,nk)=0.
      enddo
      do m=1,npress
        if(fwtpre(m).le.0.0) cycle
        nj=nj+1
        arsp(nj,nk)=1./sigpre(m)*fwtpre(m)
      enddo
!-----------------------------------------------------------------------
!--   boundary pressure constraint on P'(1), no coupling to P(1)      --
!-----------------------------------------------------------------------
      if (kpressb.eq.2) then
        nj=nj+1
        arsp(nj,nk)=0.0
      endif
      if (kprfit.ge.3.and.npresw.gt.0) then
        do m=1,npresw
          if (fwtprw(m).gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=0.0
          endif
        enddo
      endif
      if (kzeroj.gt.0) then
        do i=1,kzeroj
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
      endif
      if (kcalpa.gt.0) then
        do j=1,kcalpa
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
      endif
      if (kcgama.gt.0) then
        do j=1,kcgama
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
      endif
      if (kcomega.gt.0) then
        do j=1,kcomega
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
      endif
      if (nbdry.gt.0) then
        do j=1,nbdry
          if (fwtbdry(j).gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=0.0
          endif
        enddo
      endif
!------------------------------------------------------------------------------
!--   E coil currents                                                        --
!------------------------------------------------------------------------------
      if (iecurr.eq.2) then
        do m=1,nesum
          if (fwtec(m).gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=0.
          endif
        enddo
      endif
!------------------------------------------------------------------------------
!--   Reference flux                                                         --
!------------------------------------------------------------------------------
      if (fitsiref) then
        if (fwtref.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
      endif
!-----------------------------------------------------------------------------
!--   Summation of Fcoils current
!-----------------------------------------------------------------------------
      if (fitfcsum) then
         nj = nj + 1
         arsp(nj,nk) = 0.0
      endif
      endif kinetic
!-----------------------------------------------------------------------
!--   DELZ rigid vertical shift   96/01                               --
!-----------------------------------------------------------------------
      delz_shift: if (fitdelz.and.nniter.ge.ndelzon) then
        need=need+1
        nsavdz=need
        nk=need
        nj=0
        do m=1,nsilop
          if(fwtsi(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtsi(m)*gsildz(m)*scadelz
        enddo
        do m=1,magpri
          if(fwtmp2(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtmp2(m)*gmp2dz(m)*scadelz
        enddo
        do m=1,nstark
          if(fwtgam(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtgam(m)*rgamdz(m)*scadelz
        enddo
        do m=1,nmsels
          if(fwtbmselt(jtime,m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtbmselt(jtime,m)*rmlsdz(m)*scadelz
        enddo
        do m=1,nmsels
          if(fwtemselt(jtime,m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtemselt(jtime,m)*relsdz(m)*scadelz
        enddo
        do m=1,nece
          if(fwtece(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtece(m)*recedz(m)*scaledz
        enddo
        if (fwtecebz.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=fwtecebz*recebzdz*scaledz
        endif
        if (fwtcur.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=fwtcur*fgowdz*scadelz
        endif
        if (fwtqa.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        if (fwtbp.gt.0.0) then
          do m=2,kffcur
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (fwtdlc.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        do m=1,nfsum
          if(fwtfc(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.
        enddo
!--------------------------------------------------------------------------
!--     pressure                                                         --
!--------------------------------------------------------------------------
        if (kprfit.gt.0.and.kdofit.ne.0.and.npress.gt.0) then
          do m=1,npress
            if(fwtpre(m).le.0.0) cycle
            nj=nj+1
            arsp(nj,nk)=rpredz(m)/sigpre(m)*fwtpre(m)*scadelz
          enddo
        endif
!--------------------------------------------------------------------------
!--     P'(1)                                                            --
!--------------------------------------------------------------------------
        if (kpressb.eq.2) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
!--------------------------------------------------------------------------
!--     rotational pressure                                              --
!--------------------------------------------------------------------------
        if (kprfit.ge.3.and.npresw.gt.0) then
          do m=1,npresw
            if (fwtprw(m).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=rprwdz(m)/sigprw(m)*fwtprw(m)*scadelz
            endif
          enddo
        endif
        if (kzeroj.gt.0) then
          do i=1,kzeroj
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (kcalpa.gt.0) then
          do j=1,kcalpa
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (kcgama.gt.0) then
          do j=1,kcgama
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (kcomega.gt.0) then
          do j=1,kcomega
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (nbdry.gt.0) then
          do j=1,nbdry
            if (fwtbdry(j).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=fwtbry(j)*gbdrdz(j)*scadelz
            endif
          enddo
        endif
        if (iecurr.eq.2) then
         do m=1,nesum
          if (fwtec(m).gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=0.
          endif
         enddo
        endif
        if (fitsiref) then
          if (fwtref.gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=0.0
          endif
        endif
!-----------------------------------------------------------------------------
!--     Summation of Fcoils current
!-----------------------------------------------------------------------------
        if (fitfcsum) then
           nj = nj + 1
           arsp(nj,nk) = 0.0
        endif
      endif delz_shift
!-----------------------------------------------------------------------
!--   set up response matrix for advanced divertor coil               --
!--   lacking boundary constraints                                    --
!-----------------------------------------------------------------------
      have_acoil: if (iacoil.gt.0) then
      nkb=need+1
      need=need+nacoil
      do nk=nkb,need
        mk=nk-nkb+1
        nj=0
        do m=1,nsilop
          if(fwtsi(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtsi(m)*rsilac(m,mk)
        enddo
        do m=1,magpri
          if(fwtmp2(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtmp2(m)*rmp2ac(m,mk)
        enddo
        do m=1,nstark
          if(fwtgam(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtgam(m)*rgamac(m,mk)
        enddo
        do m=1,nece
          if(fwtece(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtece(m)*receac(m,mk)
        enddo
        if (fwtecebz.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=fwtecebz*recebzac(mk)
        endif
        if (fwtcur.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        if (fwtqa.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        if (fwtbp.gt.0.0) then
          do m=2,kffcur
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (fwtdlc.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        do m=1,nfsum
          if(fwtfc(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.
        enddo
!--------------------------------------------------------------------------
!--     pressure                                                         --
!--------------------------------------------------------------------------
        if (kprfit.gt.0.and.kdofit.ne.0.and.npress.gt.0) then
          do m=1,npress
            if(fwtpre(m).le.0.0) cycle
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
!--------------------------------------------------------------------------
!--     P'(1)                                                            --
!--------------------------------------------------------------------------
        if (kpressb.eq.2) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
!---------------------------------------------------------------------------
!--     rotational pressure                                               --
!---------------------------------------------------------------------------
        if (kprfit.ge.3.and.npresw.gt.0) then
          do m=1,npresw
            if (fwtprw(m).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=0.0
            endif
          enddo
        endif
        if (kzeroj.gt.0) then
          do i=1,kzeroj
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (kcalpa.gt.0) then
          do j=1,kcalpa
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (kcgama.gt.0) then
          do j=1,kcgama
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (kcomega.gt.0) then
          do j=1,kcomega
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
!------------------------------------------------------------------------------
!--     E coil currents                                                      --
!------------------------------------------------------------------------------
        if (iecurr.eq.2) then
          do m=1,nesum
            if (fwtec(m).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=0.
            endif
          enddo
        endif
!------------------------------------------------------------------------------
!--     Reference flux                                                       --
!------------------------------------------------------------------------------
        if (fitsiref) then
          if (fwtref.gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=0.0
          endif
        endif
!-----------------------------------------------------------------------------
!--     Summation of Fcoils current
!-----------------------------------------------------------------------------
        if (fitfcsum) then
          nj = nj + 1
          arsp(nj,nk) = 0.0
        endif
      enddo
      endif have_acoil
!-----------------------------------------------------------------------
!--   set up response matrix for E coils                              --
!-----------------------------------------------------------------------
      ecoil: if (iecurr.eq.2) then
      nkb=need+1
      need=need+nesum
      do nk=nkb,need
        mk=nk-nkb+1
        nj=0
        do m=1,nsilop
          if(fwtsi(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtsi(m)*rsilec(m,mk)
        enddo
        do m=1,magpri
          if(fwtmp2(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtmp2(m)*rmp2ec(m,mk)
        enddo
        do m=1,nstark
          if(fwtgam(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtgam(m)*rgamec(m,mk)
        enddo
        do m=1,nmsels
          if(fwtbmselt(jtime,m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
        do m=1,nmsels
          if(fwtemselt(jtime,m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
        do m=1,nece
          if(fwtece(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtece(m)*receec(m,mk)
        enddo
        if (fwtecebz.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=fwtecebz*recebzec(mk)
        endif
        if (fwtcur.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        if (fwtqa.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        if (fwtbp.gt.0.0) then
          do m=2,kffcur
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (fwtdlc.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        do m=1,nfsum
          if(fwtfc(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.
        enddo
!--------------------------------------------------------------------------
!--     pressure                                                         --
!--------------------------------------------------------------------------
        if (kprfit.gt.0.and.kdofit.ne.0.and.npress.gt.0) then
          do m=1,npress
            if(fwtpre(m).le.0.0) cycle
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
!--------------------------------------------------------------------------
!--     P'(1)                                                            --
!--------------------------------------------------------------------------
        if (kpressb.eq.2) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
!---------------------------------------------------------------------------
!--     rotational pressure                                               --
!---------------------------------------------------------------------------
        if (kprfit.ge.3.and.npresw.gt.0) then
          do m=1,npresw
            if (fwtprw(m).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=0.0
            endif
          enddo
        endif
        if (kzeroj.gt.0) then
          do i=1,kzeroj
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (kcalpa.gt.0) then
          do j=1,kcalpa
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (kcgama.gt.0) then
          do j=1,kcgama
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (kcomega.gt.0) then
          do j=1,kcomega
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
!------------------------------------------------------------------------------
!--     Boundary constraints                                                 --
!------------------------------------------------------------------------------
        if (nbdry.gt.0) then
          do j=1,nbdry
            if (fwtbdry(j).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=fwtbry(j)*rbdrec(j,mk)
            endif
          enddo
        endif
!------------------------------------------------------------------------------
!--     E coil currents                                                      --
!------------------------------------------------------------------------------
        do m=1,nesum
          if(fwtec(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.
          if (m.eq.mk) arsp(nj,nk)=fwtec(m)
        enddo
!------------------------------------------------------------------------------
!--     Reference flux                                                       --
!------------------------------------------------------------------------------
        if (fitsiref) then
          if (fwtref.gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=0.0
          endif
        endif
!-----------------------------------------------------------------------------
!--     Summation of Fcoils current
!-----------------------------------------------------------------------------
        if (fitfcsum) then
          nj = nj + 1
          arsp(nj,nk) = 0.0
        endif
      enddo
      endif ecoil
!------------------------------------------------------------------------------
!--   fitting relative flux, set up response for fitted reference flux       --
!------------------------------------------------------------------------------
      fit_flux: if (fitsiref) then
        need=need+1
        nj=0
        nk=need
        do m=1,nsilop
          if(fwtsi(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=-fwtsi(m)*scalesir
        enddo
        do m=1,magpri
          if(fwtmp2(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
        do m=1,nstark
          if(fwtgam(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
        do m=1,nmsels
          if(fwtbmselt(jtime,m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
        do m=1,nmsels
          if(fwtemselt(jtime,m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
        do m=1,nece
          if(fwtece(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
        if (fwtecebz.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        if (fwtcur.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        if (fwtqa.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        if (fwtbp.gt.0.0) then
          do m=2,kffcur
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (fwtdlc.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        do m=1,nfsum
          if(fwtfc(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.
        enddo
!--------------------------------------------------------------------------
!--     pressure                                                         --
!--------------------------------------------------------------------------
        if (kprfit.gt.0.and.kdofit.ne.0.and.npress.gt.0) then
          do m=1,npress
            if(fwtpre(m).le.0.0) cycle
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
!--------------------------------------------------------------------------
!--     P'(1)                                                            --
!--------------------------------------------------------------------------
        if (kpressb.eq.2) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
!---------------------------------------------------------------------------
!--     rotational pressure                                               --
!---------------------------------------------------------------------------
        if (kprfit.ge.3.and.npresw.gt.0) then
          do m=1,npresw
            if (fwtprw(m).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=0.0
            endif
          enddo
        endif
!---------------------------------------------------------------------------
!--     J(PSIWANT)                                                        --
!---------------------------------------------------------------------------
        if (kzeroj.gt.0) then
          do i=1,kzeroj
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
!----------------------------------------------------------------------------
!--     P', FF', and rotational                                            --
!----------------------------------------------------------------------------
        if (kcalpa.gt.0) then
          do j=1,kcalpa
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (kcgama.gt.0) then
          do j=1,kcgama
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (kcomega.gt.0) then
          do j=1,kcomega
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
!------------------------------------------------------------------------------
!--     Boundary constraints                                                 --
!------------------------------------------------------------------------------
        if (nbdry.gt.0) then
          do j=1,nbdry
            if (fwtbdry(j).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=0.0
            endif
          enddo
        endif
!------------------------------------------------------------------------------
!--     E coil currents                                                      --
!------------------------------------------------------------------------------
        if (iecurr.eq.2) then
          do m=1,nesum
            if (fwtec(m).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=0.
            endif
          enddo
        endif
!------------------------------------------------------------------------------
!--     Reference flux                                                       --
!------------------------------------------------------------------------------
        if (fitsiref) then
          if (fwtref.gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=fwtref*scalesir
          endif
        endif
!-----------------------------------------------------------------------------
!--     Summation of Fcoils current                                         --
!-----------------------------------------------------------------------------
        if (fitfcsum) then
           nj = nj + 1
           arsp(nj,nk) = 0.0
        endif
      endif fit_flux
!------------------------------------------------------------------------------
!--   set up response matrix for ER, fitting ER                              --
!------------------------------------------------------------------------------
      er_fit: if (keecur.gt.0.and.kdomse.le.0) then
      needs=need
      need=need+keecur
      do nk=needs+1,need
        nkk=nk-needs
        nj=0
        do m=1,nsilop
          if(fwtsi(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
        do m=1,magpri
          if(fwtmp2(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
        do m=1,nstark
          if(fwtgam(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtgam(m)*rgamer(m,nkk)
        enddo
        do m=1,nmsels
          if(fwtbmselt(jtime,m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
        do m=1,nmsels
          if(fwtemselt(jtime,m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtemselt(jtime,m)*relser(m,nkk)
        enddo
        do m=1,nece
          if(fwtece(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
        if (fwtecebz.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        if (fwtcur.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        if (fwtqa.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        if (fwtbp.gt.0.0) then
          do m=2,kffcur
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (fwtdlc.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        do m=1,nfsum
          if(fwtfc(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.
        enddo
!--------------------------------------------------------------------------
!--     pressure                                                         --
!--------------------------------------------------------------------------
        if (kprfit.gt.0.and.kdofit.ne.0.and.npress.gt.0) then
          do m=1,npress
            if(fwtpre(m).le.0.0) cycle
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
!--------------------------------------------------------------------------
!--     P'(1)                                                            --
!--------------------------------------------------------------------------
        if (kpressb.eq.2) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
!--------------------------------------------------------------------------
!--     rotational pressure                                                         --
!--------------------------------------------------------------------------
        if (kprfit.ge.3.and.npresw.gt.0) then
          do m=1,npresw
            if (fwtprw(m).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=0.0
            endif
          enddo
        endif
        if (kzeroj.gt.0) then
          do i=1,kzeroj
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (kcalpa.gt.0) then
          do j=1,kcalpa
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (kcgama.gt.0) then
          do j=1,kcgama
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (kcomega.gt.0) then
          do j=1,kcomega
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (nbdry.gt.0) then
          do j=1,nbdry
            if (fwtbdry(j).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=0.0
            endif
          enddo
        endif
        if (iecurr.eq.2) then
          do m=1,nesum
            if (fwtec(m).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=0.
            endif
          enddo
        endif
        if (fitsiref) then
          if (fwtref.gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=0.0
          endif
        endif
!-----------------------------------------------------------------------------
!--     Summation of Fcoils current
!-----------------------------------------------------------------------------
        if (fitfcsum) then
          nj = nj + 1
          arsp(nj,nk) = 0.0
        endif
      enddo
      endif er_fit
!-----------------------------------------------------------------------
!--   Response for Pedge hyperbolic tangent                           --
!-----------------------------------------------------------------------
      pedge_htan: if (kedgep.ne.0) then
        need=need+1
        nk=need
        npedge=need
        nj=0
        do m=1,nsilop
          if(fwtsi(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtsi(m)*rsilpe(m)
        enddo
        do m=1,magpri
          if(fwtmp2(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtmp2(m)*rmp2pe(m)
        enddo
        do m=1,nstark
          if(fwtgam(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtgam(m)*rgampe(m)
        enddo
        do m=1,nmsels
          if(fwtbmselt(jtime,m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
        do m=1,nmsels
          if(fwtemselt(jtime,m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
        if (fwtcur.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=fwtcur*fgowpe
        endif
        if (fwtqa.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        if (fwtbp.gt.0.0) then
          do m=2,kffcur
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (fwtdlc.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        do m=1,nfsum
          if(fwtfc(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.
        enddo
!--------------------------------------------------------------------
!--     pressure                                                   --
!--------------------------------------------------------------------
        if (kprfit.gt.0.and.kdofit.ne.0.and.npress.gt.0) then
          do m=1,npress
            if(fwtpre(m).le.0.0) cycle
            nj=nj+1
            arsp(nj,nk)=rprepe(m)/sigpre(m)*fwtpre(m)
          enddo
        endif
!--------------------------------------------------------------------
!--     P'(1)                                                      --
!--------------------------------------------------------------------
        if (kpressb.eq.2) then
          nj=nj+1
          ! TODO: s1edge is not defined...
!          arsp(nj,nk)=1./sigppb/darea/cosh(s1edge)**2/pe_width/sidif
          arsp(nj,nk)=0.0
        endif
!---------------------------------------------------------------------
!--     rotational pressure                                         --
!---------------------------------------------------------------------
        if (kprfit.ge.3.and.npresw.gt.0) then
          do m=1,npresw
            if (fwtprw(m).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=0.0
            endif
          enddo
        endif
!----------------------------------------------------------------------
!--     J(PSIWANT) constraint                                        --
!----------------------------------------------------------------------
        if (kzeroj.gt.0) then
          do i=1,kzeroj
            nj=nj+1
            siedge=(sizeroj(i)-pe_psin)/pe_width
            arsp(nj,nk)=rxxx(i)*fwtxxzj/cosh(siedge)**2/pe_width/sidif
          enddo
        endif
!-------------------------------------------------------------------------
!--     p' and ff' constraints                                          --
!-------------------------------------------------------------------------
        if (kcalpa.gt.0) then
          do j=1,kcalpa
            nj=nj+1
            arsp(nj,nk)=calpa(kppcur+1,j)*fwtxxa
          enddo
        endif
        if (kcgama.gt.0) then
          do j=1,kcgama
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
!---------------------------------------------------------------------
!--     rotational constraints                                      --
!---------------------------------------------------------------------
        if (kcomega.gt.0) then
          do j=1,kcomega
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
!------------------------------------------------------------------------------
!--     Boundary constraints due to plasma contributions                     --
!------------------------------------------------------------------------------
        if (nbdry.gt.0) then
          do j=1,nbdry
            if (fwtbdry(j).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=fwtbry(j)*gbdrpe(j)
            endif
          enddo
        endif
!------------------------------------------------------------------------------
!--     E coil currents                                                      --
!------------------------------------------------------------------------------
        if (iecurr.eq.2) then
          do m=1,nesum
            if (fwtec(m).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=0.
            endif
          enddo
        endif
!------------------------------------------------------------------------------
!--     Reference flux                                                       --
!------------------------------------------------------------------------------
        if (fitsiref) then
          if (fwtref.gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=0.0
          endif
        endif
!-----------------------------------------------------------------------------
!--     Summation of Fcoils current
!-----------------------------------------------------------------------------
        if (fitfcsum) then
          nj = nj + 1
          arsp(nj,nk) = 0.0
        endif
      endif pedge_htan
!-----------------------------------------------------------------------
!--   Response for f2edge hyperbolic tangent                          --
!-----------------------------------------------------------------------
      f2edge_htan: if (kedgef.ne.0) then
        need=need+1
        nk=need
        nfedge=need
        nj=0
        do m=1,nsilop
          if(fwtsi(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtsi(m)*rsilfe(m)
        enddo
        do m=1,magpri
          if(fwtmp2(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtmp2(m)*rmp2fe(m)
        enddo
        do m=1,nstark
          if(fwtgam(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtgam(m)*rgamfe(m)
        enddo
        do m=1,nmsels
          if(fwtbmselt(jtime,m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=fwtbmselt(jtime,m)*rmlsfe(m)
        enddo
        do m=1,nmsels
          if(fwtemselt(jtime,m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.0
        enddo
        if (fwtcur.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=fwtcur*fgowfe
        endif
        if (fwtqa.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        if (fwtbp.gt.0.0) then
          do m=2,kffcur
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (fwtdlc.gt.0.0) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
        do m=1,nfsum
          if(fwtfc(m).le.0.0) cycle
          nj=nj+1
          arsp(nj,nk)=0.
        enddo
!--------------------------------------------------------------------
!--     pressure                                                   --
!--------------------------------------------------------------------
        if (kprfit.gt.0.and.kdofit.ne.0.and.npress.gt.0) then
          do m=1,npress
            if(fwtpre(m).le.0.0) cycle
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
!--------------------------------------------------------------------
!--     P'(1)                                                      --
!--------------------------------------------------------------------
        if (kpressb.eq.2) then
          nj=nj+1
          arsp(nj,nk)=0.0
        endif
!---------------------------------------------------------------------
!--     rotational pressure                                         --
!---------------------------------------------------------------------
        if (kprfit.ge.3.and.npresw.gt.0) then
          do m=1,npresw
            if (fwtprw(m).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=0.0
            endif
          enddo
        endif
!----------------------------------------------------------------------
!--     J(PSIWANT) constraint                                        --
!----------------------------------------------------------------------
        if (kzeroj.gt.0) then
          do i=1,kzeroj
            nj=nj+1
            siedge=(sizeroj(i)-fe_psin)/fe_width
            arsp(nj,nk)=fwtxxzj/rxxxf(i)/cosh(siedge)**2/fe_width/sidif
          enddo
        endif
!-------------------------------------------------------------------------
!--     p' and ff' constraints                                          --
!-------------------------------------------------------------------------
        if (kcalpa.gt.0) then
          do j=1,kcalpa
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
        if (kcgama.gt.0) then
          do j=1,kcgama
            nj=nj+1
            arsp(nj,nk)=cgama(kffcur+1,j)*fwtxxf
          enddo
        endif
!---------------------------------------------------------------------
!--     rotational constraints                                      --
!---------------------------------------------------------------------
        if (kcomega.gt.0) then
          do j=1,kcomega
            nj=nj+1
            arsp(nj,nk)=0.0
          enddo
        endif
!------------------------------------------------------------------------------
!--     Boundary constraints due to plasma contributions                     --
!------------------------------------------------------------------------------
        if (nbdry.gt.0) then
          do j=1,nbdry
            if (fwtbdry(j).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=fwtbry(j)*gbdrfe(j)
            endif
          enddo
        endif
!------------------------------------------------------------------------------
!--     E coil currents                                                      --
!------------------------------------------------------------------------------
        if (iecurr.eq.2) then
          do m=1,nesum
            if (fwtec(m).gt.0.0) then
              nj=nj+1
              arsp(nj,nk)=0.
            endif
          enddo
        endif
!------------------------------------------------------------------------------
!--     Reference flux                                                       --
!------------------------------------------------------------------------------
        if (fitsiref) then
          if (fwtref.gt.0.0) then
            nj=nj+1
            arsp(nj,nk)=0.0
          endif
        endif
!-----------------------------------------------------------------------------
!--     Summation of Fcoils current
!-----------------------------------------------------------------------------
        if (fitfcsum) then
          nj = nj + 1
          arsp(nj,nk) = 0.0
        endif
      endif f2edge_htan
!-----------------------------------------------------------------------
!--   stabilization constraint for dz, setup here for all fitting paras
!-----------------------------------------------------------------------
      if (fitdelz.and.nniter.ge.ndelzon) then
        if (stabdz.gt.0.0) then
          nj=nj+1
          do i=1,need
            arsp(nj,i )=0.0
          enddo
          arsp(nj,nsavdz)=stabdz
        endif
      endif
!-----------------------------------------------------------------------
!--   user defined constraint equations
!-----------------------------------------------------------------------
      if (kccoils.gt.0) then
         do m=1,kccoils
           nj=nj+1
           do i=1,nfsum
             arsp(nj,i)=ccoils(i,m)
           enddo
         enddo
      endif
!-----------------------------------------------------------------------
!--   set up the corresponding data, RHS                              --
!-----------------------------------------------------------------------
      nj=0
      do m=1,nsilop
        if(fwtsi(m).le.0.0) cycle
        nj=nj+1
        brsp(nj)=fwtsi(m)*silopt(jtime,m)
        if (iecurr.eq.1) then
          ework=0.0
          do ne=1,nesum
            ework=ework+rsilec(m,ne)*ecurrt(ne)
          enddo
          brsp(nj)=brsp(nj)-fwtsi(m)*ework
        endif
        if(ivesel.eq.0 .or. ivesel.eq.3) cycle
        ework=0.0
        do ne=1,nvesel
          ework=ework+rsilvs(m,ne)*vcurrt(ne)
        enddo
        brsp(nj)=brsp(nj)-fwtsi(m)*ework
      enddo
      do m=1,magpri
        if(fwtmp2(m).le.0.0) cycle
        nj=nj+1
        brsp(nj)=fwtmp2(m)*expmpi(jtime,m)
        if (iecurr.eq.1) then
          ework=0.0
          do ne=1,nesum
            ework=ework+rmp2ec(m,ne)*ecurrt(ne)
          enddo
          brsp(nj)=brsp(nj)-fwtmp2(m)*ework
        endif
        if(ivesel.eq.0 .or. ivesel.eq.3) cycle
        ework=0.0
        do ne=1,nvesel
          ework=ework+rmp2vs(m,ne)*vcurrt(ne)
        enddo
        brsp(nj)=brsp(nj)-fwtmp2(m)*ework
      enddo
      do m=1,nstark
        if(fwtgam(m).le.0.0) cycle
        nj=nj+1
        brsp(nj)=fwtgam(m)*rhsgam(jtime,m)
        if (iecurr.eq.1) then
          ework=0.0
          do ne=1,nesum
            ework=ework+rgamec(m,ne)*ecurrt(ne)
          enddo
          brsp(nj)=brsp(nj)-fwtgam(m)*ework
        endif
        if(ivesel.eq.0 .or. ivesel.eq.3) cycle
        ework=0.0
        do ne=1,nvesel
          ework=ework+rgamvs(m,ne)*vcurrt(ne)
        enddo
        brsp(nj)=brsp(nj)-fwtgam(m)*ework
      enddo
!-----------------------------------------------------------------------
!--   set up the corresponding MSE-LS data + MSE-LE                   --
!-----------------------------------------------------------------------
      do m=1,nmsels
        if(fwtbmselt(jtime,m).le.0.0) cycle
        nj=nj+1
        brsp(nj)=fwtbmselt(jtime,m)*(rhsmls(jtime,m)+bmselt(jtime,m)**2)
      enddo
      do m=1,nmsels
        if(fwtemselt(jtime,m).le.0.0) cycle
        nj=nj+1
        brsp(nj)=fwtemselt(jtime,m)*emselt(jtime,m)
      enddo
      do m=1,nece
        if(fwtece(m).le.0.0) cycle
        nj=nj+1
        brsp(nj)=fwtece(m)*brspece(jtime,m)
        if (iecurr.eq.1) then
          ework=0.0
          do ne=1,nesum
            ework=ework+receec(m,ne)*ecurrt(ne)
          enddo
          brsp(nj)=brsp(nj)-fwtece(m)*ework
        endif
        if(ivesel.eq.0 .or. ivesel.eq.3) cycle
        ework=0.0
        do ne=1,nvesel
          ework=ework+recevs(m,ne)*vcurrt(ne)
        enddo
        brsp(nj)=brsp(nj)-fwtece(m)*ework
      enddo
      if (fwtecebz.gt.0.0) then
        nj=nj+1
        brsp(nj)=fwtecebz*brspecebz(jtime)
        if (iecurr.eq.1) then
          ework=0.0
          do ne=1,nesum
            ework=ework+recebzec(ne)*ecurrt(ne)
          enddo
          brsp(nj)=brsp(nj)-fwtecebz*ework
        endif
        if (ivesel.eq.1 .or. ivesel.eq.2) then
          ework=0.0
          do ne=1,nvesel
            ework=ework+recebzvs(ne)*vcurrt(ne)
          enddo
          brsp(nj)=brsp(nj)-fwtecebz*ework
        endif
      endif
      if (fwtcur.gt.0.0) then
        nj=nj+1
        brsp(nj)=fwtcur*ipmeas(jtime)
      endif
      if (fwtqa.gt.0.0) then
        nj=nj+1
        do i=1,kppcur
          arsp(nj,nfsum+i)=fwtqa*rqajx(i)
        enddo
        do i=1,kffcur
          arsp(nj,nbase+i)=fwtqa*rqafx(i)
        enddo
        if (kedgep.gt.0) then
          arsp(nj,npedge)=fwtqa*rqape
        endif
        if (kedgef.gt.0) then
          arsp(nj,nfedge)=fwtqa*rqafe
        endif
        if (kvtor.gt.0) then
          do i=1,kwwcur
            arsp(nj,nfsum+kpcurn+i)=fwtqa*rqawx(i)
          enddo
        endif
        brsp(nj)=fwtqa/qvfit*ipmeas(jtime)/abs(ipmeas(jtime))
      endif
      if (fwtbp.gt.0.0) then
        fwtbpp=fwtbp/abs(brsold(nfsum+1)*brsold(nbase+2))
        do i=2,kffcur
          nj=nj+1
          arsp(nj,nfsum+1)=fwtbpp*brsold(nbase+i)
          arsp(nj,nfsum+i)=-fwtbpp*brsold(nbase+1)
          arsp(nj,nbase+1)=-fwtbpp*brsold(nfsum+i)
          arsp(nj,nbase+i)=fwtbpp*brsold(nfsum+1)
          brsp(nj)=fwtbpp*(brsold(nfsum+1)*brsold(nbase+i)- &
                           brsold(nfsum+i)*brsold(nbase+1))
        enddo
      endif
      if (fwtdlc.gt.0.0) then
        nj=nj+1
        arsp(nj,(nbase+1):(nbase+kffcur))=fwtdlc*rspdlc(1:kffcur)
        if (kedgef.gt.0) then
          arsp(nj,nfedge)=fwtdlc*rdlcfe
        endif
        brsp(nj)=-fwtdlc*diamag(jtime)
      endif
      do i=1,nfsum
        if(fwtfc(i).le.0.0) cycle
        nj=nj+1
        brsp(nj)=fccurt(jtime,i)*fwtfc(i)
      enddo
!------------------------------------------------------------------
!--   pressure                                                   --
!------------------------------------------------------------------
      if (kprfit.gt.0.and.kdofit.ne.0.and.npress.gt.0) then
        do i=1,npress
          if(fwtpre(i).le.0.0) cycle
          nj=nj+1
          brsp(nj)=pressr(i)/sigpre(i)*fwtpre(i)
        enddo
      endif
!--------------------------------------------------------------------
!--     P'(1)                                                      --
!--------------------------------------------------------------------
      if (kpressb.eq.2) then
        nj=nj+1
        brsp(nj)=prespb/sigppb
      endif
!-------------------------------------------------------------------
!--   rotational pressure                                         --
!-------------------------------------------------------------------
      if (kprfit.ge.3.and.npresw.gt.0) then
        do i=1,npresw
          if (fwtprw(i).gt.0.0) then
          nj=nj+1
          brsp(nj)=(presw(i)-preswb)/sigprw(i)*fwtprw(i)
          endif
        enddo
      endif
!--------------------------------------------------------------------
!--   J(PSIWANT)                                                   --
!--------------------------------------------------------------------
      if (kzeroj.gt.0) then
        do i=1,kzeroj
          nj=nj+1
          brsp(nj)=fwtxxzj*vzeroj(i)*darea*ipmeas(jtime)/carea
        enddo
      endif
!--------------------------------------------------------------------
!--   P' and FF'                                                   --
!--------------------------------------------------------------------
      if (kcalpa.gt.0) then
        do j=1,kcalpa
          nj=nj+1
          brsp(nj)=fwtxxa*xalpa(j)*darea
        enddo
      endif
      if (kcgama.gt.0) then
        do j=1,kcgama
          nj=nj+1
          brsp(nj)=fwtxxf*xgama(j)*darea
        enddo
      endif
!---------------------------------------------------------------------
!--   rotational                                                    --
!---------------------------------------------------------------------
      if (kcomega.gt.0) then
        do j=1,kcomega
          nj=nj+1
          brsp(nj)=xomega(j)*darea*fwtxxo
        enddo
      endif
!------------------------------------------------------------------------------
!--   Boundary constraints                                                   --
!------------------------------------------------------------------------------
      if (nbdry.gt.0) then
        do j=1,nbdry
          if (fwtbdry(j).gt.0.0) then
            nj=nj+1
            if (iecurr.eq.1) then
              brsp(nj)=sum(rbdrec(j,:)*ecurrt)
            else
              brsp(nj)=0.0
            endif
            brsp(nj)=fwtbry(j)*(psibry-brsp(nj))
          endif
        enddo
      endif
!------------------------------------------------------------------------------
!--   E coil currents                                                        --
!------------------------------------------------------------------------------
      if (iecurr.eq.2) then
        do i=1,nesum
          if (fwtec(i).gt.0.0) then
            nj=nj+1
            brsp(nj)=ecurrt(i)*fwtec(i)
          endif
        enddo
      endif
!------------------------------------------------------------------------------
!--   fitting relative flux                                                  --
!------------------------------------------------------------------------------
      if (fitsiref) then
        if (fwtref.gt.0.0) then
          nj=nj+1
          brsp(nj)=fwtref*psiref(jtime)
        endif
      endif
!----------------------------------------------------------------------
!--   Summation of F-coils current
!----------------------------------------------------------------------
      if (fitfcsum) then
        nj = nj + 1
        brsp(nj) = 0.0
      endif
!-----------------------------------------------------------------------
!--   stabilization constraint for dz                                 --
!-----------------------------------------------------------------------
      if (fitdelz.and.nniter.ge.ndelzon) then
        if (stabdz.gt.0.0) then
          nj=nj+1
          brsp(nj)=0.0
        endif
      endif
!-----------------------------------------------------------------------
!--   user defined constraint equations
!-----------------------------------------------------------------------
      if (kccoils.gt.0) then
        do j=1,kccoils
          nj=nj+1
          brsp(nj)=xcoils(j)
        enddo
      endif
!
      ncrsp=0
      nfffff=nfsum
      call ppcnst(ncrsp,crsp,z,nfffff)
      call ffcnst(ncrsp,crsp,z,nfffff)
      call wwcnst(ncrsp,crsp,z,nfffff)
      if (keecur.gt.0.and.kdomse.eq.0) then
        nfffff=needs-kppcur-kffcur
        needer=needs
        call eecnst(ncrsp,crsp,z,nfffff)
      endif
!---------------------------------------------------------------------
!--   preconditioning A matrix if need                              --
!---------------------------------------------------------------------
      if (scalea) then
        call dgeequ(nj,need,arsp,nrsmat,rowscale,colscale, &
                    rowcnd,colcnd,arspmax,infosc)
        do j=1,nj
          arsp(j,1:need)=arsp(j,1:need)*colscale(1:need)
        enddo
      endif

      if (ncrsp .le. 0) then
        call sdecm(arsp,nrsmat,nj,need,brsp,nrsmat,nnn,wrsp,work,ier)
        if (ier.eq.129) then
          kerror=1
          call errctrl_msg('response_matrix', &
                 'problem in decomposition, sdecm failed to converge')
          return
        end if

        !-----------------------------------------------------------------------
        !--  unfold fitting parameters                                        --
        !-----------------------------------------------------------------------
        if ( wrsp(need).eq.0 ) then
          kerror=1
          call errctrl_msg('response_matrix', &
                 'problem in decomposition, wrsp(need).eq.0')
          return
        end if
        condno=wrsp(1)/wrsp(need)
        toler=condin*wrsp(1)
        do i=1,need
          t=0.0
          if(wrsp(i).gt.toler) t=brsp(i)/wrsp(i)
          work(i)=t
        end do
        do i=1,need
          brsp(i)=sum(arsp(i,1:need)*work(1:need))
        end do

      else
        b=brsp
        info=0
        call dgglse(int(nj,8),int(need,8),int(ncrsp,8),arsp,int(nrsmat,8), &
                    crsp,int(4*(npcurn-2)+6+npcurn*npcurn,8),b,z,brsp,work,&
                    int(nrsmat*2,8),int(info,8),condno)
        if (info.gt.0) then ! special hack to info in dgglse
          kerror=1
          write(tmpstr,'(a,i4,a,i4,a)') 'A(',info,',',info, &
                                        ')=0 in dgglse, divide by zero.'
          call errctrl_msg('response_matrix',tmpstr)
          return
        else if (info.lt.0) then
          kerror=1
          call errctrl_msg('response_matrix', &
                           'calling argument in dgglse was bad')
          return
        endif
      endif
!----------------------------------------------------------------------
!--   rescale results if A is preconditioned                         --
!----------------------------------------------------------------------
      if(scalea) brsp(1:need)=brsp(1:need)*colscale(1:need)
      nload=nfnwcr
      if (ivesel.eq.3) then
        if (nfourier.gt.1) then
          do j=1,nvesel
            vcurrt(j)=sum(brsp((1+nfnwcr):(nfourier*2+1+nfnwcr)) &
                          *vecta(1:(nfourier*2+1),j))
          enddo
          nload=nload+(nfourier*2+1)
        else
          vcurrt=brsp((1+nfnwcr):(nvesel+nfnwcr))
          nload=nload+nvesel
        endif
      endif
!
      if (kprfit.gt.0.and.kdofit.gt.0) then
        nload=nload+1
        prbdry=brsp(nload)
      endif
      if (fitdelz.and.nniter.ge.ndelzon) then
        nload=nload+1
        cdelz(nniter)=brsp(nload)*scaledz*100.*relaxdz
      else
        cdelz(nniter)=0.0
      endif
      if (iacoil.eq.0) then
        caccurt(jtime,:)=0.0
      else
        caccurt(jtime,:)=brsp((1+nload):(nacoil+nload))
        nload=nload+nacoil
      endif
!------------------------------------------------------------------------------
!--   E coil currents                                                        --
!------------------------------------------------------------------------------
      if (iecurr.eq.2) then
        cecurr=brsp((1+nload):(nesum+nload))
        nload=nload+nesum
      endif
!------------------------------------------------------------------------------
!--   Reference flux                                                         --
!------------------------------------------------------------------------------
      if (fitsiref) then
        nload=nload+1
        csiref=brsp(nload)*scalesir
      endif
!------------------------------------------------------------------------------
!--   Hyperbolic tangents                                                    --
!------------------------------------------------------------------------------
      if(kedgep.gt.0) pedge=brsp(npedge)
      if(kedgef.gt.0) f2edge=brsp(nfedge)
!----------------------------------------------------------------------------
!--   Update poloidal flux due to vertical shift if necessary              --
!----------------------------------------------------------------------------
      if (ifitdelz.eq.2) then
       if (fitdelz.and.nniter.ge.ndelzon) then
        cdelznow=cdelz(nniter)/100.
        call sets2d(psi,c,rgrid,nw,bkx,lkx,zgrid,nh,bky,lky,wk,ier)
        do i=1,nw
         do j=1,nh
          kk=(i-1)*nh+j
          call seva2d(bkx,lkx,bky,lky,c,rgrid(i),zgrid(j),pds,ier,n333)
          psi(kk)=psi(kk)+cdelznow*pds(3)
         enddo
        enddo
       endif
      endif
!----------------------------------------------------------------------
!--   calculate the fitting figure of merit saisq                    --
!----------------------------------------------------------------------
      saisq=0.0
      do m=1,nsilop
        cm=sum(rsilfc(m,1:nfsum)*brsp(1:nfsum))
        if(ivesel.gt.0) cm=cm+sum(rsilvs(m,:)*vcurrt)
        if (iecurr.eq.1) then
          cm=cm+sum(rsilec(m,:)*ecurrt)
        elseif (iecurr.eq.2) then
          cm=cm+sum(rsilec(m,:)*cecurr)
        endif
        if(iacoil.gt.0) cm=cm+sum(rsilac(m,:)*caccurt(jtime,:))
        cmv=cm
        cm=cm+sum(rsilpc(m,1:kwcurn)*brsp((1+nfsum):(kwcurn+nfsum)))
        if(fitsiref) cm=cm-csiref
        if(kedgep.gt.0) cm=cm+rsilpe(m)*pedge
        if(kedgef.gt.0) cm=cm+rsilfe(m)*f2edge
        if (swtsi(m).ne.0.0) then
          saisil(m)=fwtsi(m)**nsq*(silopt(jtime,m)-cm)**2
          saisil(m)=saisil(m)/swtsi(m)**nsq
        else
          saisil(m)=0.0
        endif
        saisq=saisq+saisil(m)
        csilop(m,jtime)=cm
        csilopv(m,jtime)=cmv
      enddo
!
      do m=1,magpri
        cm=sum(rmp2fc(m,1:nfsum)*brsp(1:nfsum))
        if(ivesel.gt.0) cm=cm+sum(rmp2vs(m,:)*vcurrt)
        if (iecurr.eq.1) then
          cm=cm+sum(rmp2ec(m,:)*ecurrt)
        elseif (iecurr.eq.2) then
          cm=cm+sum(rmp2ec(m,:)*cecurr)
        endif
        if(iacoil.gt.0) cm=cm+sum(rmp2ac(m,:)*caccurt(jtime,:))
        cmv=cm
        cm=cm+sum(rmp2pc(m,1:kwcurn)*brsp((1+nfsum):(kwcurn+nfsum)))
        if(kedgep.gt.0) cm=cm+rmp2pe(m)*pedge
        if(kedgef.gt.0) cm=cm+rmp2fe(m)*f2edge
        if (swtmp2(m).ne.0.0) then
          saimpi(m)=fwtmp2(m)**nsq*(expmpi(jtime,m)-cm)**2
          saimpi(m)=saimpi(m)/swtmp2(m)**nsq
        else
          saimpi(m)=0.0
        endif
        saisq=saisq+saimpi(m)
        cmpr2(m,jtime)=cm
        cmpr2v(m,jtime)=cmv
      enddo
!
      MSE: if (kstark.gt.0) then
      chigamt=0.0
      do m=1,nstark
        chigam(m)=0.0
        cmgam(m,jtime)=0.0
        if(rrgam(jtime,m).le.0.0) cycle
        cmbr=sum(rbrfc(m,1:nfsum)*brsp(1:nfsum)) &
            +sum(rbrpc(m,1:kwcurn)*brsp((1+nfsum):(kwcurn+nfsum)))
        cmbz=sum(rbzfc(m,1:nfsum)*brsp(1:nfsum)) &
            +sum(rbzpc(m,1:kwcurn)*brsp((1+nfsum):(kwcurn+nfsum)))
        if (ivesel.gt.0) then
          cmbr=cmbr+sum(rbrvs(m,:)*vcurrt)
          cmbz=cmbz+sum(rbzvs(m,:)*vcurrt)
        endif
        if (iecurr.eq.1) then
          cmbr=cmbr+sum(rbrec(m,:)*ecurrt)
          cmbz=cmbz+sum(rbzec(m,:)*ecurrt)
        elseif (iecurr.eq.2) then
          cmbr=cmbr+sum(rbrec(m,:)*cecurr)
          cmbz=cmbz+sum(rbzec(m,:)*cecurr)
        endif
        if (iacoil.gt.0) then
          cmbr=cmbr+sum(rbrac(m,:)*caccurt(jtime,:))
          cmbz=cmbz+sum(rbzac(m,:)*caccurt(jtime,:))
        endif
        if (kedgep.gt.0) then
          cmbr=cmbr+rbrpe(m)*pedge
          cmbz=cmbz+rbzpe(m)*pedge
        endif
        if (kedgef.gt.0) then
          cmbr=cmbr+rbrfe(m)*f2edge
          cmbz=cmbz+rbzfe(m)*f2edge
        endif
        cm=a2gam(jtime,m)*btgam(m)+a3gam(jtime,m)*cmbr &
          +a4gam(jtime,m)*cmbz
        bzmsec(m)=cmbz
        if (keecur.le.0) then
          bzmse(m)=(tangam(jtime,m)*cm-a8gam(jtime,m)*cmbr) &
                   /a1gam(jtime,m)
          cm=(a1gam(jtime,m)*cmbz+a8gam(jtime,m)*cmbr)/cm
        else
          cerer(1:keecur)=brsp((1+needer):(keecur+needer))
          erup=sum(e1rbz(m,1:keecur)*cerer(1:keecur))
          erbot=sum((e2rbz(m,1:keecur)+e3rbr(m,1:keecur))*cerer(1:keecur))
          cm=cm-erbot
          bzmse(m)=(tangam(jtime,m)*cm+erup-a8gam(jtime,m)*cmbr) &
                   /a1gam(jtime,m)
          cm=(a1gam(jtime,m)*cmbz+a8gam(jtime,m)*cmbr-erup)/cm
          ermse(m)=-erup/a5gam(jtime,m)
        endif
        if (swtgam(m).ne.0.0) then
          chigam(m)=fwtgam(m)**nsq*(tangam(jtime,m)-cm)**2
          chigam(m)=chigam(m)/swtgam(m)**nsq
        else
          chigam(m)=0.0
        endif
        cmgam(m,jtime)=cm
      enddo
      chigamt=chigamt+sum(chigam(1:nmselp))
      chi2gamt(jtime)=chigamt
      chilibt=sum(chigam((nmselp+1):nstark))
      if (ishot.le.97400) then
        mcentral=15
      else
        mcentral=10
      endif
      do m=1,mcentral
        mp1=m+1
        drgam=rrgam(jtime,mp1)-rrgam(jtime,m)
        cjmse(m)=-(bzmse(mp1)-bzmse(m))/drgam/twopi/tmu
        cjmsec(m)=-(bzmsec(mp1)-bzmsec(m))/drgam/twopi/tmu
      enddo
      endif MSE
!
      tchimls=0.0
      do m=1,nmsels
        cm=sqrt(sum(rmlspc(m,1:kwcurn)*brsp((1+nfsum):(kwcurn+nfsum))) &
                -rhsmls(jtime,m))
        if (swtbmsels(m).ne.0.0) then
          chimls(m)=fwtbmselt(jtime,m)**nsq*(bmselt(jtime,m)-cm)**2
          chimls(m)=chimls(m)/swtbmsels(m)**nsq
        else
          chimls(m)=0.0
        endif
        cmmls(jtime,m)=cm
        tchimls=tchimls+chimls(m)
      enddo
      chi2mls(jtime)=tchimls
!
      tchiels=0.0
      do m=1,nmsels
        cm=sum(relser(m,1:keecur)*cecurr(1:keecur))
        if (swtemsels(m).ne.0.0) then
          chiels(m)=fwtemselt(jtime,m)**nsq*(emselt(jtime,m)-cm)**2
          chiels(m)=chiels(m)/swtemsels(m)**nsq
        else
          chiels(m)=0.0
        endif
        cmels(jtime,m)=cm
        tchiels=tchiels+chiels(m)
      enddo
!
      tchiece=0.0
      do m=1,nece
        cm=sum(recefc(m,1:nfsum)*brsp(1:nfsum)) &
          +sum(recepc(m,1:kwcurn)*brsp((1+nfsum):(kwcurn+nfsum)))
        if(ivesel.gt.0) cm=cm+sum(recevs(m,:)*vcurrt)
        if (iecurr.eq.1) then 
          cm=cm+sum(receec(m,:)*ecurrt)
        elseif (iecurr.eq.2) then
          cm=cm+sum(receec(m,:)*cecurr)
        endif
        if(iacoil.gt.0) cm=cm+sum(receac(m,:)*caccurt(jtime,:))
        if (swtece(m).ne.0.0) then
          chiece(m)=fwtece(m)**nsq*(brspece(jtime,m)-cm)**2
          chiece(m)=chiece(m)/swtece(m)**nsq
        else
          chiece(m)=0.0
        endif
        cmece(m,jtime)=cm
        tchiece=tchiece+chiece(m)
      enddo
!
      cm=sum(recebzfc(1:nfsum)*brsp(1:nfsum)) &
        +sum(recebzpc(1:kwcurn)*brsp((1+nfsum):(kwcurn+nfsum)))
      if(ivesel.gt.0) cm=cm+sum(recevs(m,:)*vcurrt) !TODO: should this be recebzvs?
      if (iecurr.eq.1) then
        cm=cm+sum(recebzec*ecurrt)
      elseif (iecurr.eq.2) then
        cm=cm+sum(recebzec*cecurr)
      endif
      if(iacoil.gt.0) cm=cm+sum(receac(m,:)*caccurt(jtime,:)) !TODO: should this be recebzac?
      if (swtecebz.ne.0.0) then
        chiecebz=fwtecebz**nsq*(brspecebz(jtime)-cm)**2
        chiecebz=chiecebz/swtecebz**nsq
      else
        chiecebz=0.0
      endif
      cmecebz(jtime)=cm
!
      cm=sum(brsp((1+nfsum):nfnwcr)*fgowpc(1:(nfnwcr-nfsum)))
      if(kedgep.gt.0) cm=cm+fgowpe*pedge
      if(kedgef.gt.0) cm=cm+fgowfe*f2edge
      ipmhd(jtime)=cm
      if (kfffnc.eq.8) then
        cjeccd=brsp(nfnwcr)*fgowpc(nfnwcr-nfsum)/1000.
      else
        cjeccd=0.0
      endif
      if(ivesel.eq.3) cm=cm+sum(vcurrt)
      if (swtcur.ne.0.0) then
        chipasma=(fwtcur/swtcur)**nsq*(ipmeas(jtime)-cm)**2
      else
        chipasma=0.0
      endif
      saisq=saisq+chipasma
!
      tchifcc=0.0
      do i=1,nfsum
        chifcc(i)=0.0
        if (fwtfc(i).gt.0.0) then
          chifcc(i)=fwtfc(i)**nsq*(brsp(i)-fccurt(jtime,i))**2
          chifcc(i)=chifcc(i)/swtfc(i)**nsq
        endif
        saisq=saisq+chifcc(i)
        tchifcc=tchifcc+chifcc(i)
      enddo
!
      if (iecurr.eq.2) then
        do i=1,nesum
          chiecc(i)=0.0
          if (fwtec(i).gt.0.0) then
            chiecc(i)=fwtec(i)**nsq*(cecurr(i)-ecurrt(i))**2
            chiecc(i)=chiecc(i)/swtec(i)**nsq
          endif
          saisq=saisq+chiecc(i)
        enddo
      endif
!
      if (fitsiref) then
        saisref=0.0
        if (fwtref.gt.0.0) then
          saisref=fwtref**nsq*(psiref(jtime)-csiref)**2
          saisref=saisref/swtsi(nslref)**nsq
        endif
        saisq=saisq+saisref
      endif
!--------------------------------------------------------------------------
!--   pressure                                                           --
!--------------------------------------------------------------------------
      chipre=0.0
      if (kprfit.gt.0.and.kdofit.ne.0.and.npress.gt.0) then
        do m=1,npress
          cm=sum(rprepc(m,1:kppcur)*brsp((1+nfsum):(kppcur+nfsum)))
          if(kedgep.gt.0) cm=cm+rprepe(m)*pedge
          cm=cm+prbdry
          if (fwtpre(m).gt.0.0) then
            saipre(m)=((cm-pressr(m))/sigpre(m))**2
          else
            saipre(m)=0.0
          endif
          saipre2(m)=saipre(m)  ! preserve saipre - changed later by pltout
          chipre=chipre+saipre(m)
          precal(m)=cm
        enddo
      endif
!--------------------------------------------------------------------------
!--   rotational pressure                                                --
!--------------------------------------------------------------------------
      chiprw=0.0
      rotational_pressure: if (kprfit.ge.3.and.npresw.gt.0) then
        do m=1,npresw
          cm=sum(rprwpc(m,1:kwwcur)*brsp((1+nfnpcr):(kwwcur+nfnpcr)))
          cm=cm+preswb
          if (fwtprw(m).gt.0.0) then
            saiprw(m)=((cm-presw(m))/sigprw(m))**2
          else
            saiprw(m)=0.0
          endif
          saiprw2(m)=saiprw(m)  ! preserve saiprw - changed later by pltout
          chiprw=chiprw+saiprw(m)
          prwcal(m)=cm
        enddo
      endif rotational_pressure
!--------------------------------------------------------------------------
!--   check iconvr=2 criteria to determine if fit should be stopped
!--------------------------------------------------------------------------
      if (iconvr.eq.2) then
       if ((nniter.ge.minite).or.((eouter.le.elomin).and.(fwtdlc.le.0.0))) then
        if ((nniter.ge.kcallece).or.(kfitece.le.0.0)) then
         if ((errorm.le.errmin).or.((eouter.le.elomin).and.(fwtdlc.le.0.0))) then
          if (saisq.le.saicon) then
           if (abs(saisq-saiold).le.0.10_dp .or. saisq.ge.saiold) then
            ! criteria satisfied, restore previous solution and stop fit 
            ichisq=1
            brsp=brsold
            if(ivesel.eq.3) vcurrt=vcurrto
            saisq=saiold
           endif
          endif
         endif
        endif
       endif
      endif
      chisq(jtime)=saisq
!--------------------------------------------------------------------------
!--   write status to fitout.dat
!--------------------------------------------------------------------------
      if (iand(iout,1).ne.0) then
        write (nout,7400) time(jtime),chipre,ipmhd(jtime), &
          nniter,condno,saisq,chigamt
        ! TODO: chimlst and chielst are never defined in efit...
!        write (nout,97400) time(jtime),chimlst,chielst
        write (nout,7445) need
        write (nout,7450) (brsp(i),i=1,need)
        write (nout,7450) (wrsp(i),i=1,need)
      endif

      return

 7400 format (/,2x,'time = ',e12.5,2x,'chipr = ',e12.5, &
              2x,'current = ',e12.5,/,2x,'it = ',i5, &
              1x,'condn = ',1pe11.4, &
              1x,'chisq = ',1pe11.4,1x,'chigam = ',1pe11.4, &
              /,2x,'chiecebz = ',1pe11.4,1x,'tchieceR+R- = ', &
              1pe11.4)
!97400 format (/,2x,'time = ',e12.5,2x,'chimls= ',e12.5, &
!              2x,'chiels = ',e12.5)
 7445 format (10x,'fitting parameters:   ',i5)
 7450 format (8(1x,e12.5,1x))
 7460 format (10x,'chi ip:',/,15x,e12.5)
      end subroutine response_matrix
