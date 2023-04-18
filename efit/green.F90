#include "config.f"
!**********************************************************************
!>
!!    green sets up the appropriate response functions for use
!!    with the routine response_matrix.
!!
!!    @param iflag :
!!
!!    @param jtime : time index
!!
!!    @param niter : iteration number
!!
!**********************************************************************
      subroutine green(iflag,jtime,niter)
      use commonblocks,only: c,bkx,bky
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      real*8 pwcurr,prcurr,erpote
      integer*4, intent(in) :: iflag,jtime,niter
      integer*4 i,j,jj,jjk,k,kk,m,ier
      real*8 factor,ecorrect,epotp,kffp1,pres0,prew0,ptop0,pwop0,pwp0r2, &
             qconst,rmggvt,rmsnow,rmvjj,rmvnow,rmvtor,siedge,upsi,upsi1, &
             wwwww,xmsinow,xn,xpsdb,xpsdd,xpsnow,xpzero,ypsi,zmsnow
      real*8 xpsii(nwcurn),xpsis(nwcurn),xpsisb(nwcurn),xsier(nercur), &
             pds(6),xnsi(nppcur)
!      
      kffp1=kffcur+1
      do jj=1,kwcurn
        if (iflag.ne.1) then
          rsilpc(:,jj)=0.0
          rmp2pc(:,jj)=0.0
          if (kstark.gt.0.or.kdomse.gt.0) then
            rbrpc(:,jj)=0.0
            rbzpc(:,jj)=0.0
          endif
          if(mmbmsels.gt.0.or.kdomsels.gt.0) rmlspc(1:nmsels,jj)=0.0
          if(kece.gt.0) recepc(:,jj)=0.0
          if(kecebz.gt.0) recebzpc(jj)=0.0
          if(nbdry.gt.0) gbdrpc(1:nbdry,jj)=0.0
        endif
        fgowpc(jj)=0.0
      enddo
      if (fitdelz.and.niter.ge.ndelzon) then
        gsildz=0.0
        gmp2dz=0.0
        fgowdz=0.0
        if (kstark.gt.0) then
          gbrdz=0.0
          gbzdz=0.0
        endif
        if(kece.gt.0) gecedz=0.0
        if(kecebz.gt.0) gecebzdz=0.0
        if(nbdry.gt.0) gbdrdz(1:nbdry)=0.0
      endif
      if(iflag.ne.1) rspdlc(1:kffcur)=0.0
!
      if (fwtdlc.gt.0.0) then
        upsi1=1.0
        call setff(upsi1,xpsisb)
      endif
!------------------------------------------------------------------
!--   Hyperbolic tangent term                                    --
!------------------------------------------------------------------
      if (kedgep.gt.0) then
        if (iflag.ne.1) then
          rsilpe=0.0
          rmp2pe=0.0
          if (kstark.gt.0.or.kdomse.gt.0) then
            rbrpe=0.0
            rbzpe=0.0
          endif
          if(nbdry.gt.0) gbdrpe(1:nbdry)=0.0
        endif
        fgowpe=0.0
      endif
!
      if (kedgef.gt.0) then
        if (iflag.ne.1) then
          rsilfe=0.0
          rmp2fe=0.0
          if (kstark.gt.0.or.kdomse.gt.0) then
            rbrfe=0.0
            rbzfe=0.0
          endif
          if(nbdry.gt.0) gbdrfe(1:nbdry)=0.0
          rdlcfe=0.0
        endif
        fgowfe=0.0
      endif
!
      do i=1,nw
        do j=1,nh
          kk=(i-1)*nh+j
          if ((xpsi(kk).ge.0.0).and.(xpsi(kk).le.1.0)) then
            call setpp(xpsi(kk),xpsii)
            do jj=1,kppcur
              factor=xpsii(jj)*rgrid(i)*www(kk)
!------------------------------------------------------------------------
!--           correction for rotation, p0' terms                       --
!------------------------------------------------------------------------
              if (niter.gt.1) then
                if (kvtor.eq.2) then
                  prew0=pwcurr(xpsi(kk),kwwcur)
                  pres0=prcurr(xpsi(kk),kppcur)
                  if (abs(pres0).gt.1.e-10_dp) then
                    pwop0=prew0/pres0
                    pwp0r2=pwop0*rgrvt(i)
                  else
                    pwop0=0.0
                    pwp0r2=0.0
                  endif
                  factor=factor*(1.-0.5_dp*pwp0r2**2)
                elseif (kvtor.eq.3) then
                  prew0=pwcurr(xpsi(kk),kwwcur)
                  pres0=prcurr(xpsi(kk),kppcur)
                  if (abs(pres0).gt.1.e-10_dp) then
                    pwop0=prew0/pres0
                    pwp0r2=pwop0*rgrvt(i)
                    ptop0=exp(pwp0r2)
                  else
                    pwop0=0.0
                    pwp0r2=0.0
                    ptop0=1.0
                  endif
                  factor=factor*ptop0*(1.-pwp0r2)
                endif
              endif
              if (iflag.ne.1) then
                rsilpc(:,jj)=rsilpc(:,jj)+gsilpc(:,kk)*factor
                rmp2pc(:,jj)=rmp2pc(:,jj)+gmp2pc(:,kk)*factor
                if (kstark.gt.0.or.kdomse.gt.0) then
                  rbrpc(:,jj)=rbrpc(:,jj)+gbrpc(:,kk)*factor
                  rbzpc(:,jj)=rbzpc(:,jj)+gbzpc(:,kk)*factor
                endif
                if(kece.gt.0) &
                  recepc(:,jj)=recepc(:,jj)+gecepc(:,kk)*factor
                if(kecebz.gt.0) &
                  recebzpc(jj)=recebzpc(jj)+gecebzpc(kk)*factor
!---------------------------------------------------------------------------
!--             Boundary constraints                                      --
!---------------------------------------------------------------------------
                if (nbdry.gt.0) then
                  do m=1,nbdry
                    gbdrpc(m,jj)=gbdrpc(m,jj)+rbdrpc(m,kk)*factor
                  enddo
                endif
              endif
              fgowpc(jj)=fgowpc(jj)+factor
            enddo
!-------------------------------------------------------------------------
!--         Hyperbolic tangent term for P'                              --
!-------------------------------------------------------------------------
            if (kedgep.gt.0) then
              siedge=(xpsi(kk)-pe_psin)/pe_width
              xpsnow=1./pe_width/sidif/cosh(siedge)**2
              factor=xpsnow*rgrid(i)*www(kk)
              if (iflag.ne.1) then
                rsilpe=rsilpe+gsilpc(:,kk)*factor
                rmp2pe=rmp2pe+gmp2pc(:,kk)*factor
                if (kstark.gt.0.or.kdomse.gt.0) then
                  rbrpe=rbrpe+gbrpc(:,kk)*factor
                  rbzpe=rbzpe+gbzpc(:,kk)*factor
                endif
                if (nbdry.gt.0) then
                  do m=1,nbdry
                    gbdrpe(m)=gbdrpe(m)+rbdrpc(m,kk)*factor
                  enddo
                endif
              endif
              fgowpe=fgowpe+factor
            endif
          endif
!-------------------------------------------------------------------------
!--       attached current for ff' ?                                    --
!-------------------------------------------------------------------------
          if (icutfp.eq.0) then
            upsi=xpsi(kk)
            wwwww=www(kk)
          else
            upsi=xpsi(kk)*xpsimin
            wwwww=zero(kk)
          endif
          if ((upsi.ge.0.0).and.(upsi.le.1.0)) then
            call setfp(upsi,xpsii)
            if(fwtdlc.gt.0.0) call setff(upsi,xpsis)
            do jj=kppcur+1,kpcurn
              jjk=jj-kppcur
              factor=xpsii(jjk)/rgrid(i)*wwwww
              if (iflag.ne.1) then
                rsilpc(:,jj)=rsilpc(:,jj)+gsilpc(:,kk)*factor
                rmp2pc(:,jj)=rmp2pc(:,jj)+gmp2pc(:,kk)*factor
                if (kstark.gt.0.or.kdomse.gt.0) then
                  rbrpc(:,jj)=rbrpc(:,jj)+gbrpc(:,kk)*factor
                  rbzpc(:,jj)=rbzpc(:,jj)+gbzpc(:,kk)*factor
                endif
                if(kece.gt.0) &
                  recepc(:,jj)=recepc(:,jj)+gecepc(:,kk)*factor
                if(kecebz.gt.0) &
                  recebzpc(jj)=recebzpc(jj)+gecebzpc(kk)*factor
!---------------------------------------------------------------------------
!--             Boundary constraints                                      --
!---------------------------------------------------------------------------
                if (nbdry.gt.0) then
                  do m=1,nbdry
                    gbdrpc(m,jj)=gbdrpc(m,jj)+rbdrpc(m,kk)*factor
                  enddo
                endif
              endif
              fgowpc(jj)=fgowpc(jj)+factor
              if (iflag.eq.1) cycle
              if (fwtdlc.gt.0.0) then
                xpsdd=xpsis(jjk)
                xpsdb=xpsisb(jjk)
                rspdlc(jjk)=(xpsdd-xpsdb) &
                           /rgrid(i)*www(kk)+rspdlc(jjk)
              endif
            enddo
!----------------------------------------------------------------------
!--         Hyperbolic tangent term                                  --
!----------------------------------------------------------------------
            if (kedgef.gt.0) then
              siedge=(upsi-fe_psin)/fe_width
              xpsnow=1./fe_width/sidif/cosh(siedge)**2
              if (icutfp.gt.0) xpsnow=xpsnow*xpsimin
              factor=xpsnow/rgrid(i)*wwwww
              if (iflag.ne.1) then
                rsilfe=rsilfe+gsilpc(:,kk)*factor
                rmp2fe=rmp2fe+gmp2pc(:,kk)*factor
                if (kstark.gt.0.or.kdomse.gt.0) then
                  rbrfe=rbrfe+gbrpc(:,kk)*factor
                  rbzfe=rbzfe+gbzpc(:,kk)*factor
                endif
                if (nbdry.gt.0) then
                  do m=1,nbdry
                    gbdrfe(m)=gbdrfe(m)+rbdrpc(m,kk)*factor
                  enddo
                endif
                if (fwtdlc.gt.0.0) then
                  xpsdd=tanh(siedge)
                  xpsdb=tfedge
                  rdlcfe=(xpsdd-xpsdb) &
                        /rgrid(i)*www(kk)+rdlcfe
                endif
              endif
              fgowfe=fgowfe+factor
            endif
          endif
!----------------------------------------------------------------------
!--       toroidal rotation terms                                    --
!----------------------------------------------------------------------
          rotation: if ((kvtor.gt.0).and.(xpsi(kk).ge.0.0).and.(xpsi(kk).le.1.0)) &
            then
            call setpwp(xpsi(kk),xpsii)
            do jj=kpcurn+1,kwcurn
              jjk=jj-kpcurn
              factor=xpsii(jjk)*rgsvt(i)*www(kk)
              if (niter.gt.1) then
                if (kvtor.eq.2) then
                  factor=factor*(1.+pwp0r2)
                elseif (kvtor.eq.3) then
                  factor=factor*ptop0
                endif
              endif
              if (iflag.ne.1) then
                rsilpc(:,jj)=rsilpc(:,jj)+gsilpc(:,kk)*factor
                rmp2pc(:,jj)=rmp2pc(:,jj)+gmp2pc(:,kk)*factor
                if (kstark.gt.0.or.kdomse.gt.0) then
                  rbrpc(:,jj)=rbrpc(:,jj)+gbrpc(:,kk)*factor
                  rbzpc(:,jj)=rbzpc(:,jj)+gbzpc(:,kk)*factor
                endif
                if(kece.gt.0) &
                  recepc(:,jj)=recepc(:,jj)+gecepc(:,kk)*factor
                if(kecebz.gt.0) &
                  recebzpc(jj)=recebzpc(jj)+gecebzpc(kk)*factor
!---------------------------------------------------------------------------
!--             Boundary constraints                                      --
!---------------------------------------------------------------------------
                if (nbdry.gt.0) then
                  do m=1,nbdry
                    gbdrpc(m,jj)=gbdrpc(m,jj)+rbdrpc(m,kk)*factor
                  enddo
                endif
              endif
              fgowpc(jj)=fgowpc(jj)+factor
            enddo
          endif rotation
        enddo
      enddo
!----------------------------------------------------------------------
!--   diamagnetic flux                                               --
!----------------------------------------------------------------------
      if (fwtdlc.gt.0.0) then
        do jj=1,kffcur
          rspdlc(jj)=rspdlc(jj)*sidif*twopi/fbrdy
          if (icutfp.gt.0) rspdlc(jj)=rspdlc(jj)/xpsimin
        enddo
        rdlcfe=rdlcfe*twopi/fbrdy
      endif
!-----------------------------------------------------------------------
!--   axial q constraint                                              --
!-----------------------------------------------------------------------
      qconst=twopi*emaxis*rmaxis**2/(emaxis**2+1.)/abs(fcentr)/darea
      xpzero=0.0
      call setpp(xpzero,xpsii)
      rmvtor=qconst*rmaxis
      rmvjj=rmaxis
      if(kvtor.gt.0) rmggvt=(rmaxis/rvtor)**2-1.
      if (niter.gt.1) then
        if (kvtor.eq.2) then
          prew0=pwcurr(xpzero,kwwcur)
          pres0=prcurr(xpzero,kppcur)
          if (abs(pres0).gt.1.e-10_dp) then
            pwop0=prew0/pres0
            pwp0r2=pwop0*rmggvt
          else
            pwop0=0.0
            pwp0r2=0.0
          endif
          rmvtor=rmvtor*(1.-0.5_dp*pwp0r2**2)
          rmvjj=rmvjj*(1.-0.5_dp*pwp0r2**2)
        elseif (kvtor.eq.3) then
          prew0=pwcurr(xpzero,kwwcur)
          pres0=prcurr(xpzero,kppcur)
          if (abs(pres0).gt.1.e-10_dp) then
            pwop0=prew0/pres0
            pwp0r2=pwop0*rmggvt
            ptop0=exp(pwp0r2)
          else
            pwop0=0.0
            pwp0r2=0.0
            ptop0=1.0
          endif
          rmvtor=rmvtor*ptop0*(1.-pwp0r2)
          rmvjj=rmvjj*ptop0*(1.-pwp0r2)
        endif
      endif
      rqajtor=rmvtor
      do jj=1,kppcur
        factor=xpsii(jj)
        rqajx(jj)=rmvtor*factor
        rjjjx(jj)=rmvjj*factor
      enddo
      if (kedgep.gt.0) then
        siedge=(xpzero-pe_psin)/pe_width
        rqapetor=rmvtor/pe_width/sidif
        rqape=rqapetor/cosh(siedge)**2
      endif
!
      call setfp(xpzero,xpsii)
      rmvtor=qconst/rmaxis
      rqaftor=rmvtor
      do jj=1,kffcur
        factor=xpsii(jj)
        rqafx(jj)=rmvtor*factor
        rjjfx(jj)=factor/rmaxis
      enddo
      if (kedgef.gt.0) then
        siedge=(xpzero-fe_psin)/fe_width
        rqafetor=rmvtor/fe_width/sidif
        rqafe=rqafetor/cosh(siedge)**2
      endif
!
      if (kvtor.gt.0) then
        call setpwp(xpzero,xpsii)
        rmvnow=rmaxis*rmggvt
        if (niter.gt.1) then
          if (kvtor.eq.2) then
            rmvnow=rmvnow*(1.+pwp0r2)
          endif
          if (kvtor.eq.3) then
            rmvnow=rmvnow*ptop0
          endif
        endif
        rmvtor=qconst*rmvnow
        do jj=1,kwwcur
          factor=xpsii(jj)
          rqawx(jj)=rmvtor*factor
          rjjwx(jj)=rmvnow*factor
        enddo
      endif
!----------------------------------------------------------------------
!--   response functions for MSE                                     --
!----------------------------------------------------------------------
      MSE: if (kstark.gt.0.or.kdomse.gt.0) then
        do m=1,nstark
          if (fwtgam(m).gt.0.0) then
            do jj=1,kwcurn
              rgampc(m,jj)=(a3gam(jtime,m)*tangam(jtime,m) &
                          - a8gam(jtime,m))*rbrpc(m,jj) &
                          +(a4gam(jtime,m)*tangam(jtime,m) &
                          - a1gam(jtime,m))*rbzpc(m,jj)
            enddo
            do jj=1,nfsum
              rgamfc(m,jj)=(a3gam(jtime,m)*tangam(jtime,m) &
                          - a8gam(jtime,m))*rbrfc(m,jj) &
                          +(a4gam(jtime,m)*tangam(jtime,m) &
                          - a1gam(jtime,m))*rbzfc(m,jj)
            enddo
            do jj=1,nesum
              rgamec(m,jj)=(a3gam(jtime,m)*tangam(jtime,m) &
                          - a8gam(jtime,m))*rbrec(m,jj) &
                          +(a4gam(jtime,m)*tangam(jtime,m) &
                          - a1gam(jtime,m))*rbzec(m,jj)
            enddo
            if (ivesel.eq.3) then
              do jj=1,nvesel
                rgamvs(m,jj)=(a3gam(jtime,m)*tangam(jtime,m) &
                            - a8gam(jtime,m))*rbrvs(m,jj) &
                            +(a4gam(jtime,m)*tangam(jtime,m) &
                            - a1gam(jtime,m))*rbzvs(m,jj)
              enddo
            endif
            rhsgam(jtime,m)=-tangam(jtime,m)*btgam(m)*a2gam(jtime,m)
            if (kedgep.gt.0) then
              rgampe(m)=(a3gam(jtime,m)*tangam(jtime,m) &
                       - a8gam(jtime,m))*rbrpe(m) &
                       +(a4gam(jtime,m)*tangam(jtime,m) &
                       - a1gam(jtime,m))*rbzpe(m)
            endif
            if (kedgef.gt.0) then
              rgamfe(m)=(a3gam(jtime,m)*tangam(jtime,m) &
                       - a8gam(jtime,m))*rbrfe(m) &
                       +(a4gam(jtime,m)*tangam(jtime,m) &
                       - a1gam(jtime,m))*rbzfe(m)
            endif
          endif
          if (keecur.gt.0) then
            rmsnow=rrgam(jtime,m)
            zmsnow=zzgam(jtime,m)
            call seva2d(bkx,lkx,bky,lky,c,rmsnow,zmsnow,pds,ier,n333)
            xmsinow=(simag-pds(1))/sidif ! TODO: simag is unset here...
            call seter(xmsinow,xsier)
            e1rbz(m,1:keecur)=a5gam(jtime,m)*pds(2)*xsier(1:keecur)
            e2rbz(m,1:keecur)=a7gam(jtime,m)*pds(2)*xsier(1:keecur)
            e3rbr(m,1:keecur)=a6gam(jtime,m)*pds(3)*xsier(1:keecur)
            rgamer(m,1:keecur)=-(e2rbz(m,1:keecur) &
                               + e3rbr(m,1:keecur))*tangam(jtime,m) &
                               + e1rbz(m,1:keecur)
          endif
        enddo
      endif MSE
!----------------------------------------------------------------------
!--   response functions for MSE-LS                                  --
!----------------------------------------------------------------------
#ifdef DEBUG_MSELS
      write (6,*) 'GREEN mmbmsels = ',mmbmsels
#endif
      if (mmbmsels.gt.0.or.kdomsels.gt.0) then
        do m=1,nmsels
          if ((fwtbmselt(jtime,m).gt.0.0).or.(kdomsels.ne.0)) then
            if (rrmselt(jtime,m).le.0.0) cycle
            call seva2d(bkx,lkx,bky,lky,c,rrmselt(jtime,m), &
                        zzmselt(jtime,m),pds,ier,n333)
            xn=(simag-pds(1))/sidif ! TODO: simag is unset here...
            if (xn.ge.1.0) then
              rmlspc(m,(kppcur+1):kpcurn)=fmlscut
            else
              call setff(xn,xnsi)
#ifdef DEBUG_LEVEL2
              write (6,*) 'GREEN MSELS m,xn,xnsi,tmu02,dsi,dA= ',m,xn, &
                          xnsi(1),xnsi(2),xnsi(3),tmu02,sidif,darea
#endif
              do jj=kppcur+1,kpcurn
                jjk=jj-kppcur
                rmlspc(m,jj)=-l1mselt(jtime,m)*sidif/darea*xnsi(jjk)*tmu02
                rmlspc(m,jj)=rmlspc(m,jj)/rrmselt(jtime,m)**2
              enddo
            endif
            rhsmls(jtime,m)=-l1mselt(jtime,m)* &
              (rcentr*bcentr(jtime)/rrmselt(jtime,m))**2
            ecorrect=1.0
            ! TODO: mmemsels never defined here, was mmbmsels intended?
!            if (mmemsels.gt.0) then
            if (mmbmsels.gt.0) then
              epotp=erpote(ypsi,keecur)
              ecorrect=ecorrect-l4mselt(jtime,m)*rrmselt(jtime,m)*epotp
            endif
            rhsmls(jtime,m)=rhsmls(jtime,m)-l2mselt(jtime,m)*ecorrect* &
              brmls(m)*btmls(m)
            rhsmls(jtime,m)=rhsmls(jtime,m)-l3mselt(jtime,m)*ecorrect**2* &
              brmls(m)**2
            rhsmls(jtime,m)=rhsmls(jtime,m)-(l3mselt(jtime,m)*ecorrect**2+ &
              l1mselt(jtime,m))*bzmls(m)**2
          endif
          ! TODO: see above
!          if (mmemsels.gt.0) then
          if (mmbmsels.gt.0) then
            call seter(xn,xsier)
            relser(m,1:keecur)=xsier(1:keecur)
          endif
        enddo
#ifdef DEBUG_LEVEL2
        m=3
        write (6,*) 'GREEN MSELS m,rmlspc,rhsmls= ',m,rmlspc(m,3), &
                    rhsmls(1,m)
#else
#ifdef DEBUG_MSELS
        m=3
        write (6,*) 'GREEN MSELS m,rmlspc,rhsmls= ',m,rmlspc(m,3), &
                    rhsmls(1,m)
#endif
#endif
      endif
!-----------------------------------------------------------------------
!--   response functions for q constraints at PSIWANT                 --
!-----------------------------------------------------------------------
      if (nqwant.gt.0) then
        do k=1,nqwant
          do jj=1,kwcurn
             fgowsw(jj,k)=0.0
          enddo
          do i=1,nw
            do j=1,nh
              kk=(i-1)*nh+j
              if ((xpsi(kk).ge.0.0).and.(xpsi(kk).le.siwantq(k))) &
                then
                call setpp(xpsi(kk),xpsii)
                do jj=1,kppcur
                  factor=xpsii(jj)*rgrid(i)*www(kk)
                  fgowsw(jj,k)=fgowsw(jj,k)+factor
                enddo
              endif
              if (icutfp.eq.0) then
                upsi=xpsi(kk)
                wwwww=www(kk)
              else
                upsi=xpsi(kk)*xpsimin
                wwwww=zero(kk)
              endif
              if ((upsi.lt.0.0).or.(upsi.gt.psiwant)) cycle
              call setfp(upsi,xpsii)
              do jj=kppcur+1,kpcurn
                jjk=jj-kppcur
                factor=xpsii(jjk)/rgrid(i)*wwwww
                fgowsw(jj,k)=fgowsw(jj,k)+factor
              enddo
            enddo
          enddo
        enddo
      endif
      if (iflag.eq.1) return
!--------------------------------------------------------------------------
!--   deltaz in fitting                                                  --
!--------------------------------------------------------------------------
      if (fitdelz.and.niter.ge.ndelzon) then
        do i=1,nw
          do j=1,nh
            kk=(i-1)*nh+j
            rdjdz(kk)=0.0
            if ((xpsi(kk).ge.0.0).and.(xpsi(kk).le.1.0)) then
              call setppp(xpsi(kk),xpsii)
              call seva2d(bkx,lkx,bky,lky,c,rgrid(i),zgrid(j),pds,ier,n333)
              do jj=1,kppcur
                factor=xpsii(jj)*rgrid(i)*www(kk)*pds(3)*brsp(nfsum+jj)
                factor=-factor/sidif
                rdjdz(kk)=rdjdz(kk)+factor
                gsildz=gsildz+gsilpc(:,kk)*factor
                gmp2dz=gmp2dz+gmp2pc(:,kk)*factor
                if (kstark.gt.0) then
                  gbrdz=gbrdz+gbrpc(:,kk)*factor
                  gbzdz=gbzdz+gbzpc(:,kk)*factor
                endif
                if(kece.gt.0) &
                  recedz=recedz+gecepc(:,kk)*factor
                if(kecebz.gt.0) &
                  recebzdz=recebzdz+gecebzpc(kk)*factor
!---------------------------------------------------------------------------
!--             Boundary constraints                                      --
!---------------------------------------------------------------------------
                if (nbdry.gt.0) then
                  do m=1,nbdry
                    gbdrdz(m)=gbdrdz(m)+rbdrpc(m,kk)*factor
                  enddo
                endif
                fgowdz=fgowdz+factor
              enddo
            endif
!-------------------------------------------------------------------------
!--         attached current for ff' ?                                  --
!-------------------------------------------------------------------------
            if (icutfp.eq.0) then
              upsi=xpsi(kk)
              wwwww=www(kk)
            else
              upsi=xpsi(kk)*xpsimin
              wwwww=zero(kk)
            endif
            if ((upsi.ge.0.0).and.(upsi.le.1.0)) then
              call setfpp(upsi,xpsii)
              do jj=kppcur+1,kpcurn
                jjk=jj-kppcur
                factor=xpsii(jjk)/rgrid(i)*wwwww*pds(3)*brsp(nfsum+jj)
                factor=-factor/sidif
                rdjdz(kk)=rdjdz(kk)+factor
                gsildz=gsildz+gsilpc(:,kk)*factor
                gmp2dz=gmp2dz+gmp2pc(:,kk)*factor
                if (kstark.gt.0) then
                  gbrdz=gbrdz+gbrpc(:,kk)*factor
                  gbzdz=gbzdz+gbzpc(:,kk)*factor
                endif
                if(kece.gt.0) &
                  recedz=recedz+gecepc(:,kk)*factor
                if(kecebz.gt.0) &
                  recebzdz=recebzdz+gecebzpc(kk)*factor
!---------------------------------------------------------------------------
!--             Boundary constraints                                      --
!---------------------------------------------------------------------------
                if (nbdry.gt.0) then
                  do m=1,nbdry
                    gbdrdz(m)=gbdrdz(m)+rbdrpc(m,kk)*factor
                  enddo
                endif
                fgowdz=fgowdz+factor
              enddo
            endif
!----------------------------------------------------------------------
!--         toroidal rotation contributions                      --
!----------------------------------------------------------------------
            tor_rot: if ((kvtor.gt.0).and.(xpsi(kk).ge.0.0).and.(xpsi(kk).le.1.0)) &
              then
              call setpwpp(xpsi(kk),xpsii)
              do jj=kpcurn+1,kwcurn
                jjk=jj-kpcurn
                factor=xpsii(jjk)*rgsvt(i)*www(kk)*pds(3)*brsp(nfsum+jj)
                factor=-factor/sidif
                rdjdz(kk)=rdjdz(kk)+factor
                gsildz=gsildz+gsilpc(:,kk)*factor
                gmp2dz=gmp2dz+gmp2pc(:,kk)*factor
                if (kstark.gt.0) then
                  gbrdz=gbrdz+gbrpc(:,kk)*factor
                  gbzdz=gbzdz+gbzpc(:,kk)*factor
                endif
                if(kece.gt.0) &
                  recedz=recedz+gecepc(:,kk)*factor
                if(kecebz.gt.0) &
                  recebzdz=recebzdz+gecebzpc(kk)*factor
!---------------------------------------------------------------------------
!--             Boundary constraints                                      --
!---------------------------------------------------------------------------
                if (nbdry.gt.0) then
                  do m=1,nbdry
                    gbdrdz(m)=gbdrdz(m)+rbdrpc(m,kk)*factor
                  enddo
                endif
                fgowdz=fgowdz+factor
              enddo
            endif tor_rot
          enddo
        enddo
!----------------------------------------------------------------------
!--     response functions for MSE                                   --
!----------------------------------------------------------------------
        if (kstark.gt.0) then
          do m=1,nstark
            if (fwtgam(m).le.0.0) cycle
              rgamdz(m)=(a3gam(jtime,m)*tangam(jtime,m) &
                       - a8gam(jtime,m))*gbrdz(m) &
                       +(a4gam(jtime,m)*tangam(jtime,m) &
                       - a1gam(jtime,m))*gbzdz(m)
          enddo
        endif
      endif
      return
      end subroutine green
