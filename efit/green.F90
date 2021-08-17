!**********************************************************************
!>
!!    green sets up the appropriate response functions for use
!!    with the routine matrix.
!!
!!    @param ifag :
!!
!!    @param jtime :
!!
!!    @param niter :
!!
!**********************************************************************
      subroutine green(ifag,jtime,niter)
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      dimension xpsii(nwcurn),xpsis(nwcurn),xpsisb(nwcurn)
      dimension xsier(nercur)
      dimension pds(6),xnsi(nppcur)
!
      kffp1=kffcur+1
      do jj=1,kwcurn
        if (ifag.ne.1) then
          do m=1,nsilop
            rsilpc(m,jj)=0.0
          enddo
          do m=1,magpri
            rmp2pc(m,jj)=0.0
          enddo
          if (kstark.gt.0.or.kdomse.gt.0) then
            do m=1,nstark
              rbrpc(m,jj)=0.0
              rbzpc(m,jj)=0.0
           enddo
          endif
          if (mmbmsels.gt.0.or.kdomsels.gt.0) then
            do m=1,nmsels
              rmlspc(m,jj)=0.0
            enddo
          endif
          if (kece.gt.0) then
            do m=1,nnece
              recepc(m,jj)=0.0
            enddo
          endif
          if (kecebz.gt.0) then
            recebzpc(jj)=0.0
          endif
          if (nbdry.gt.0) then
            do m=1,nbdry
              gbdrpc(m,jj)=0.0
            enddo
          endif
        endif
        fgowpc(jj)=0.0
      enddo
      if (fitdelz.and.niter.ge.ndelzon) then
        do m=1,nsilop
          gsildz(m)=0.0
        enddo
        do m=1,magpri
          gmp2dz(m)=0.0
        enddo
        fgowdz=0.0
        if (kstark.gt.0) then
          do m=1,nstark
            gbrdz(m)=0.0
            gbzdz(m)=0.0
          enddo
        endif
        if (kece.gt.0) then
          do m=1,nnece
            gecedz(m)=0.0
          enddo
        endif
        if (kecebz.gt.0) then
          gecebzdz=0.0
        endif
        if (nbdry.gt.0) then
          do m=1,nbdry
            gbdrdz(m)=0.0
          enddo
        endif
      endif
      if (ifag.ne.1) then
        do jj=1,kffcur
          rspdlc(jj)=0.0
        enddo
      endif
!
      if (fwtdlc.gt.0.0) then
        upsi1=1.0
        call setff(upsi1,xpsisb)
      endif
!------------------------------------------------------------------
!--   Hyperbolic tangent term                                    --
!------------------------------------------------------------------
      if (kedgep.gt.0) then
        if (ifag.ne.1) then
          do m=1,nsilop
            rsilpe(m)=0.0
          enddo
          do m=1,magpri
            rmp2pe(m)=0.0
          enddo
          if (kstark.gt.0.or.kdomse.gt.0) then
            do m=1,nstark
              rbrpe(m)=0.0
              rbzpe(m)=0.0
            enddo
          endif
          if (nbdry.gt.0) then
            do m=1,nbdry
              gbdrpe(m)=0.0
            enddo
          endif
        endif
        fgowpe=0.0
      endif
!
      if (kedgef.gt.0) then
        if (ifag.ne.1) then
          do m=1,nsilop
            rsilfe(m)=0.0
          enddo
          do m=1,magpri
            rmp2fe(m)=0.0
          enddo
          if (kstark.gt.0.or.kdomse.gt.0) then
            do m=1,nstark
              rbrfe(m)=0.0
              rbzfe(m)=0.0
            enddo
          endif
          if (nbdry.gt.0) then
            do m=1,nbdry
              gbdrfe(m)=0.0
            enddo
          endif
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
              if (ifag.ne.1) then
                do m=1,nsilop
                  rsilpc(m,jj)=rsilpc(m,jj) + gsilpc(m,kk)*factor
                enddo
                do m=1,magpri
                  rmp2pc(m,jj)=rmp2pc(m,jj) + gmp2pc(m,kk)*factor
                enddo
                if (kstark.gt.0.or.kdomse.gt.0) then
                  do m=1,nstark
                    rbrpc(m,jj)=rbrpc(m,jj) + gbrpc(m,kk)*factor
                    rbzpc(m,jj)=rbzpc(m,jj) + gbzpc(m,kk)*factor
                  enddo
                endif
                if (kece.gt.0) then
                  do m=1,nnece
                    recepc(m,jj)=recepc(m,jj) + gecepc(m,kk)*factor
                  enddo
                endif
                if (kecebz.gt.0) then
                  recebzpc(jj)=recebzpc(jj) + gecebzpc(kk)*factor
                endif
!---------------------------------------------------------------------------
!--             Boundary constraints                                      --
!---------------------------------------------------------------------------
                if (nbdry.gt.0) then
                  do m=1,nbdry
                    gbdrpc(m,jj)=gbdrpc(m,jj) + rbdrpc(m,kk)*factor
                  enddo
                endif
              endif
              fgowpc(jj)=fgowpc(jj) + factor
            enddo
!-------------------------------------------------------------------------
!--         Hyperbolic tangent term for P'                              --
!-------------------------------------------------------------------------
            if (kedgep.gt.0) then
              siedge=(xpsi(kk)-pe_psin)/pe_width
              xpsnow=1./pe_width/sidif/cosh(siedge)**2
              factor=xpsnow*rgrid(i)*www(kk)
              if (ifag.ne.1) then
                do m=1,nsilop
                  rsilpe(m)=rsilpe(m) + gsilpc(m,kk)*factor
                enddo
                do m=1,magpri
                  rmp2pe(m)=rmp2pe(m) + gmp2pc(m,kk)*factor
                enddo
                if (kstark.gt.0.or.kdomse.gt.0) then
                  do m=1,nstark
                    rbrpe(m)=rbrpe(m) + gbrpc(m,kk)*factor
                    rbzpe(m)=rbzpe(m) + gbzpc(m,kk)*factor
                  enddo
                endif
                if (nbdry.gt.0) then
                  do m=1,nbdry
                    gbdrpe(m)=gbdrpe(m) + rbdrpc(m,kk)*factor
                  enddo
                endif
              endif
              fgowpe=fgowpe + factor
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
            if (fwtdlc.gt.0.0) then
              call setff(upsi,xpsis)
            endif
            do jj=kppcur+1,kpcurn
              jjk=jj-kppcur
              factor=xpsii(jjk)/rgrid(i)*wwwww
              if (ifag.ne.1) then
                do m=1,nsilop
                  rsilpc(m,jj)=rsilpc(m,jj) + gsilpc(m,kk)*factor
                enddo
                do m=1,magpri
                  rmp2pc(m,jj)=rmp2pc(m,jj) + gmp2pc(m,kk)*factor
                enddo
                if (kstark.gt.0.or.kdomse.gt.0) then
                  do m=1,nstark
                    rbrpc(m,jj)=rbrpc(m,jj) + gbrpc(m,kk)*factor
                    rbzpc(m,jj)=rbzpc(m,jj) + gbzpc(m,kk)*factor
                  enddo
                endif
                if (kece.gt.0) then
                  do m=1,nece
                    recepc(m,jj)=recepc(m,jj) + gecepc(m,kk)*factor
                  enddo
                endif
                if (kecebz.gt.0) then
                  recebzpc(jj)=recebzpc(jj) + gecebzpc(kk)*factor
                endif
!---------------------------------------------------------------------------
!--             Boundary constraints                                      --
!---------------------------------------------------------------------------
                if (nbdry.gt.0) then
                  do m=1,nbdry
                    gbdrpc(m,jj)=gbdrpc(m,jj) + rbdrpc(m,kk)*factor
                  enddo
                endif
              endif
              fgowpc(jj)=fgowpc(jj)+factor
              if (ifag.eq.1) cycle
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
              if (ifag.ne.1) then
                do m=1,nsilop
                  rsilfe(m)=rsilfe(m) + gsilpc(m,kk)*factor
                enddo
                do m=1,magpri
                  rmp2fe(m)=rmp2fe(m) + gmp2pc(m,kk)*factor
                enddo
                if (kstark.gt.0.or.kdomse.gt.0) then
                  do m=1,nstark
                    rbrfe(m)=rbrfe(m) + gbrpc(m,kk)*factor
                    rbzfe(m)=rbzfe(m) + gbzpc(m,kk)*factor
                  enddo
                endif
                if (nbdry.gt.0) then
                  do m=1,nbdry
                    gbdrfe(m)=gbdrfe(m) + rbdrpc(m,kk)*factor
                  enddo
                endif
                if (fwtdlc.gt.0.0) then
                  xpsdd=tanh(siedge)
                  xpsdb=tfedge
                  rdlcfe=(xpsdd-xpsdb) &
                        /rgrid(i)*www(kk)+rdlcfe
                endif
              endif
              fgowfe=fgowfe + factor
            endif
          endif
!----------------------------------------------------------------------
!--       toroidal rotation terms                                    --
!----------------------------------------------------------------------
          if ((kvtor.gt.0).and.(xpsi(kk).ge.0.0).and.(xpsi(kk).le.1.0)) &
            then
            call setpwp(xpsi(kk),xpsii)
            do jj=kpcurn+1,kwcurn
              jjii=jj-kpcurn
              factor=xpsii(jjii)*rgsvt(i)*www(kk)
              if (niter.gt.1) then
                if (kvtor.eq.2) then
                  factor=factor*(1.+pwp0r2)
                elseif (kvtor.eq.3) then
                  factor=factor*ptop0
                endif
              endif
              if (ifag.ne.1) then
                do m=1,nsilop
                  rsilpc(m,jj)=rsilpc(m,jj) + gsilpc(m,kk)*factor
                enddo
                do m=1,magpri
                  rmp2pc(m,jj)=rmp2pc(m,jj) + gmp2pc(m,kk)*factor
                enddo
                if (kstark.gt.0.or.kdomse.gt.0) then
                  do m=1,nstark
                    rbrpc(m,jj)=rbrpc(m,jj) + gbrpc(m,kk)*factor
                    rbzpc(m,jj)=rbzpc(m,jj) + gbzpc(m,kk)*factor
                  enddo
                endif
                if (kece.gt.0) then
                  do m=1,nece
                    recepc(m,jj)=recepc(m,jj) + gecepc(m,kk)*factor
                  enddo
                endif
                if (kecebz.gt.0) then
                  recebzpc(jj)=recebzpc(jj) + gecebzpc(kk)*factor
                endif
!---------------------------------------------------------------------------
!--             Boundary constraints                                      --
!---------------------------------------------------------------------------
                if (nbdry.gt.0) then
                  do m=1,nbdry
                    gbdrpc(m,jj)=gbdrpc(m,jj) + rbdrpc(m,kk)*factor
                  enddo
                endif
              endif
              fgowpc(jj)=fgowpc(jj) + factor
            enddo
          endif
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
      if (kvtor.gt.0) rmggvt=(rmaxis/rvtor)**2-1.
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
        endif
        if (kvtor.eq.3) then
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
      rqajtor = rmvtor
      do jj=1,kppcur
        factor= xpsii(jj)
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
      rqaftor = rmvtor
      do jj=1,kffcur
        factor= xpsii(jj)
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
          factor= xpsii(jj)
          rqawx(jj)=rmvtor*factor
          rjjwx(jj)=rmvnow*factor
        enddo
      endif
!----------------------------------------------------------------------
!--   response functions for MSE                                     --
!----------------------------------------------------------------------
      if (kstark.gt.0.or.kdomse.gt.0) then
        do m=1,nstark
          if (fwtgam(m).gt.0.0) then
            do jj=1,kwcurn
              rgampc(m,jj)=(a3gam(jtime,m)*tangam(jtime,m) &
                          - a8gam(jtime,m))*rbrpc(m,jj) &
                          +(a4gam(jtime,m)*tangam(jtime,m) &
                          - a1gam(jtime,m))*rbzpc(m,jj)
            enddo
            do jj=1,nfcoil
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
            if (ifitvs.gt.0) then
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
            xmsinow=(simag-pds(1))/sidif
            call seter(xmsinow,xsier)
            do jj=1,keecur
              e1rbz(m,jj)=a5gam(jtime,m)*pds(2)*xsier(jj)
              e2rbz(m,jj)=a7gam(jtime,m)*pds(2)*xsier(jj)
              e3rbr(m,jj)=a6gam(jtime,m)*pds(3)*xsier(jj)
              rgamer(m,jj)=-(e2rbz(m,jj) &
                           + e3rbr(m,jj))*tangam(jtime,m) &
                           + e1rbz(m,jj)
            enddo
          endif
        enddo
      endif
!----------------------------------------------------------------------
!--   response functions for MSE-LS                                  --
!----------------------------------------------------------------------
      if (jdebug.eq.'MSEL') write (6,*) 'GREEN mmbmsels = ',mmbmsels
        if (mmbmsels.gt.0.or.kdomsels.gt.0) then
        do m=1,nmsels
          if ((fwtbmselt(jtime,m).gt.0.0).or.(kdomsels.ne.0)) then
            if (rrmselt(jtime,m).le.0.0) cycle
            call seva2d(bkx,lkx,bky,lky,c,rrmselt(jtime,m), &
                        zzmselt(jtime,m),pds,ier,n333)
            xn=(simag-pds(1))/sidif
            if (xn.ge.1.0) then
              do jj=kppcur+1,kpcurn
                rmlspc(m,jj)=fmlscut
              enddo
            else
              call setff(xn,xnsi)
              if (idebug>=2) then
                write (6,*) 'GREEN MSELS m,xn,xnsi,tmu02,dsi,dA= ',m,xn, &
                            xnsi(1),xnsi(2),xnsi(3),tmu02,sidif,darea
              endif
              do jj=kppcur+1,kpcurn
                mjj=jj-kppcur
                rmlspc(m,jj)=-l1mselt(jtime,m)*sidif/darea*xnsi(mjj)*tmu02
                rmlspc(m,jj)=rmlspc(m,jj)/rrmselt(jtime,m)**2
              enddo
            endif
            rhsmls(jtime,m)=-l1mselt(jtime,m)* &
              (rcentr*bcentr(jtime)/rrmselt(jtime,m))**2
            ecorrect=1.0
            if (mmemsels.gt.0) then
              epotp=erpote(ypsi,keecur)
              ecorrect=ecorrect-l4mselt(jtime,m)*rrmselt(jtime,m)*epotp
            endif
            rhsmls(jtime,m)=rhsmls(jtime,m)- l2mselt(jtime,m)*ecorrect* &
              brmls(m)*btmls(m)
            rhsmls(jtime,m)=rhsmls(jtime,m)- l3mselt(jtime,m)*ecorrect**2* &
              brmls(m)**2
            rhsmls(jtime,m)=rhsmls(jtime,m)-(l3mselt(jtime,m)*ecorrect**2+ &
              l1mselt(jtime,m))*bzmls(m)**2
          endif
          if (mmemsels.gt.0) then
            call seter(xn,xsier)
            do jj=1,keecur
              relser(m,jj)=xsier(jj)
            enddo
          endif
        enddo
        if (idebug>=2.or.jdebug.eq.'MSEL') then
          m=3
          write (6,*) 'GREEN MSELS m,rmlspc,rhsmls= ',m,rmlspc(m,3), &
                      rhsmls(1,m)
        endif
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
                  fgowsw(jj,k)=fgowsw(jj,k) + factor
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
                jji=jj-kppcur
                factor=xpsii(jji)/rgrid(i)*wwwww
                fgowsw(jj,k)=fgowsw(jj,k)+factor
              enddo
            enddo
          enddo
        enddo
      endif
      if (ifag.eq.1) return
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
                factor=xpsii(jj)*rgrid(i)*www(kk)*pds(3)*brsp(nfcoil+jj)
                factor=-factor/sidif
                rdjdz(kk)=rdjdz(kk)+factor
                do m=1,nsilop
                  gsildz(m)=gsildz(m) + gsilpc(m,kk)*factor
                enddo
                do m=1,magpri
                  gmp2dz(m)=gmp2dz(m) + gmp2pc(m,kk)*factor
                enddo
                if (kstark.gt.0) then
                  do m=1,nstark
                    gbrdz(m)=gbrdz(m) + gbrpc(m,kk)*factor
                    gbzdz(m)=gbzdz(m) + gbzpc(m,kk)*factor
                  enddo
                endif
                if (kece.gt.0) then
                  do m=1,nece
                    recedz(m)=recedz(m) + gecepc(m,kk)*factor
                  enddo
                endif
                if (kecebz.gt.0) then
                  recebzdz=recebzdz + gecebzpc(kk)*factor
                endif
!---------------------------------------------------------------------------
!--             Boundary constraints                                      --
!---------------------------------------------------------------------------
                if (nbdry.gt.0) then
                  do m=1,nbdry
                    gbdrdz(m)=gbdrdz(m) + rbdrpc(m,kk)*factor
                  enddo
                endif
                fgowdz=fgowdz + factor
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
                factor=xpsii(jjk)/rgrid(i)*wwwww*pds(3)*brsp(nfcoil+jj)
                factor=-factor/sidif
                rdjdz(kk)=rdjdz(kk)+factor
                do m=1,nsilop
                  gsildz(m)=gsildz(m) + gsilpc(m,kk)*factor
                enddo
                do m=1,magpri
                  gmp2dz(m)=gmp2dz(m) + gmp2pc(m,kk)*factor
                enddo
                if (kstark.gt.0) then
                  do m=1,nstark
                    gbrdz(m)=gbrdz(m) + gbrpc(m,kk)*factor
                    gbzdz(m)=gbzdz(m) + gbzpc(m,kk)*factor
                  enddo
                endif
                if (kece.gt.0) then
                  do m=1,nece
                    recedz(m)=recedz(m) + gecepc(m,kk)*factor
                  enddo
                endif
                if (kecebz.gt.0) then
                  recebzdz=recebzdz + gecebzpc(kk)*factor
                endif
!---------------------------------------------------------------------------
!--             Boundary constraints                                      --
!---------------------------------------------------------------------------
                if (nbdry.gt.0) then
                  do m=1,nbdry
                    gbdrdz(m)=gbdrdz(m) + rbdrpc(m,kk)*factor
                  enddo
                endif
                fgowdz=fgowdz + factor
              enddo
            endif
!----------------------------------------------------------------------
!--         toroidal rotation contributions                      --
!----------------------------------------------------------------------
            if ((kvtor.gt.0).and.(xpsi(kk).ge.0.0).and.(xpsi(kk).le.1.0)) &
              then
              call setpwpp(xpsi(kk),xpsii)
              do jj=kpcurn+1,kwcurn
                jjii=jj-kpcurn
                factor=xpsii(jjii)*rgsvt(i)*www(kk)*pds(3)*brsp(nfcoil+jj)
                factor=-factor/sidif
                rdjdz(kk)=rdjdz(kk)+factor
                do m=1,nsilop
                  gsildz(m)=gsildz(m) + gsilpc(m,kk)*factor
                enddo
                do m=1,magpri
                  gmp2dz(m)=gmp2dz(m) + gmp2pc(m,kk)*factor
                enddo
                if (kstark.gt.0) then
                  do m=1,nstark
                    gbrdz(m)=gbrdz(m) + gbrpc(m,kk)*factor
                    gbzdz(m)=gbzdz(m) + gbzpc(m,kk)*factor
                  enddo
                endif
                if (kece.gt.0) then
                  do m=1,nece
                    recedz(m)=recedz(m) + gecepc(m,kk)*factor
                  enddo
                endif
                if (kecebz.gt.0) then
                  recebzdz=recebzdz + gecebzpc(kk)*factor
                endif
!---------------------------------------------------------------------------
!--             Boundary constraints                                      --
!---------------------------------------------------------------------------
                if (nbdry.gt.0) then
                  do m=1,nbdry
                    gbdrdz(m)=gbdrdz(m) + rbdrpc(m,kk)*factor
                  enddo
                endif
                fgowdz=fgowdz + factor
              enddo
            endif
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
