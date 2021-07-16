!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          currnt computes the current density on the r-z mesh.    **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      subroutine currnt(iter,jtime,ixn,nitett,kerror)
      use commonblocks,only: c,wk,copy,bkx,bky
      use set_kinds
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      dimension pds(6)
      dimension alipc(npcur3,nwcurn),xpspp(nppcur),xpsfp(nffcur)
      dimension wlipc(nwcurn),work(nwcur2),xrsp(npcur3),xpspwp(nwwcur)
      dimension crsp(4*(npcurn-2)+6+npcurn*npcurn,nrsmat)
      dimension b(nrsmat),z(4*(npcurn-2)+6+npcurn*npcurn)
      real*8 :: tcurrt, tcurrtpp, tcurrtffp
      character(len=128) tmpstr
      data initc/0/
      data ten24/1.e4_dp/

      kerror = 0

      initc=initc+1
      if (ivacum.gt.0) return
      if ((nitett.le.1).and.(icinit.eq.1)) return
      if (icinit.gt.0) then
        if ((iconvr.ne.3).and.(iter.le.1)) go to 3100
      endif
      select case (icurrt)
        case (1)
          go to 100
        case (2)
          go to 1100
        case (3)
          go to 2100
        case (4)
          go to 3100
        case (5)
          go to 5100
      end select
  100 continue
!------------------------------------------------------------------------------
!--  uniform current for Solove equilibrium                                  --
!------------------------------------------------------------------------------
      tcurrt=0.0
      do 130 i=1,nw
      do 130 j=1,nh
        kk=(i-1)*nh+j
        pcurrt(kk)=0.0
        pcurrw(kk)=0.0
        if ((xpsi(kk).lt.0.0).or.(xpsi(kk).gt.1.0)) go to 130
        rdiml=rgrid(i)/srma
        pcurrt(kk)=sbeta*rdiml+2.*salpha/rdiml
        if (kvtor.gt.0) then
          pcurrw(kk)=sbetaw*rdiml*(rdiml**2-1.)
          pcurrt(kk)=pcurrt(kk)+pcurrw(kk)
        endif
        pcurrt(kk)=pcurrt(kk)*www(kk)
        pcurrw(kk)=pcurrw(kk)*www(kk)
        tcurrt=tcurrt+pcurrt(kk)
  130 continue
      cratio=cpasma(jtime)/tcurrt
      do kk=1,nwnh
        pcurrt(kk)=pcurrt(kk)*cratio
        pcurrw(kk)=pcurrw(kk)*cratio
      enddo
      cj0=cratio/darea/2.
      fcentr=twopi*2.*saaa**2*cj0
      fbrdy=fcentr*sqrt(1.-4.*salpha)
      bcentr(jtime)=fbrdy*tmu/rcentr
      return
!
 1100 continue
!----------------------------------------------------------------------
!--  polynomial current profile                                      --
!----------------------------------------------------------------------
      if ((nitett.le.1).and.(icinit.gt.0)) go to 3100
      nnn=1
      call green(nnn,jtime,nitett)
      if ((nitett.le.1).and.(icinit.lt.0)) go to 1800
      if (iconvr.eq.3) then
        if (kcgama.gt.0.or.kcalpa.gt.0) go to 1200
      endif
      if ((qenp.le.0.0).or.(nitett.le.1)) go to 1800
      if ((qvfit.le.0.0).or.(fwtqa.le.0.0)) go to 1800
      if ((kedgep.gt.0).or.(kedgef.gt.0)) go to 1800
      tz = 0.0
      aqax=1.0/qvfit/(rqajtor*ppcurr(tz,kppcur) &
                   +rqaftor*fpcurr(tz,kffcur))
      aqax = abs(aqax)
      do i=1,kppcur
        brsp(nfcoil+i)=aqax*brsp(nfcoil+i)
      enddo
      do i=1,kffcur
        brsp(nbase+i)=aqax*brsp(nbase+i)
      enddo
      cwant0=cpasma(jtime)-fgowpc(1)*brsp(nfcoil+1)-fgowpc(kppcur+1) &
             *brsp(nbase+1)
      cwant1=0.0
      do 1150 i=2,kpcurn
        if (i.eq.kppcur+1) go to 1150
        cwant1=cwant1+fgowpc(i)*brsp(nfcoil+i)
 1150 continue
      if (abs(cwant1).le.1.0e-10_dp) then
        kerror=1
        call errctrl_msg('currnt','abs(cwant1) <= 1.0e-10')
        return
      end if

      cwant1=cwant0/cwant1
      do 1170 i=2,kpcurn
        if (i.eq.kppcur+1) go to 1170
        brsp(nfcoil+i)=cwant1*brsp(nfcoil+i)
 1170 continue
      go to 1800
!----------------------------------------------------------------------
!--  Adjust current profile to keep q(0), I, J(1), and others fixed  --
!----------------------------------------------------------------------
 1200 continue
      nj=0
      if (fwtqa.gt.0.0) then
        nj=nj+1
        do 1210 j=1,kppcur
          alipc(nj,j)=rqajx(j)*fwtqa
 1210   continue
        do j=1,kffcur
          alipc(nj,kppcur+j)=rqafx(j)*fwtqa
        enddo
        nnow=kpcurn
        if (kedgep.gt.0) then
          nnow=nnow+1
          alipc(nj,nnow)=fwtqa*rqape
        endif
        if (kedgef.gt.0) then
          nnow=nnow+1
          alipc(nj,nnow)=fwtqa*rqafe
        endif
        xrsp(nj)=fwtqa/qvfit*pasmat(jtime)/abs(pasmat(jtime))
      endif
!
      if (fwtcur.gt.0.0) then
      fwtcux=fwtcur
        nj=nj+1
        do 1220 j=1,kpcurn
          alipc(nj,j)=fwtcux*fgowpc(j)
 1220   continue
        nnow=kpcurn
        if (kedgep.gt.0) then
          nnow=nnow+1
          alipc(nj,nnow)=fwtcux*fgowpe
        endif
        if (kedgef.gt.0) then
          nnow=nnow+1
          alipc(nj,nnow)=fwtcux*fgowfe
        endif
        xrsp(nj)=fwtcux*pasmat(jtime)
      endif
!-----------------------------------------------------------------------
!--  constraints on q at psiwant by successive iterations             --
!-----------------------------------------------------------------------
      if (nqwant.gt.0)   then
        do 13226 i=1,nqwant
        nj=nj+1
        if (initc.ge.jwantm) then
          fwtqqq=fwtxxq/0.001_dp/pasmsw(i)
          do 13220 j=1,kpcurn
            alipc(nj,j)=fwtqqq*fgowsw(j,i)
13220     continue
          xrsp(nj)=fwtqqq/fqsiw(i)*qsiw(i)*pasmsw(i)
        endif
13226   continue
      endif
!-----------------------------------------------------------------------
!--  J at PSIWANT constraint                                          --
!-----------------------------------------------------------------------
      brspmin=max(ten24,abs(brsp(nfcoil+1)))
      fwtxxx=fwtxxj*1000./brspmin
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        nj=nj+1
!----------------------------------------------------------------------
!--  local or flux surface averaged J constraint ?                   --
!----------------------------------------------------------------------
        if (rzeroj(i).gt.0.0) rxxx=rzeroj(i)
        if (rzeroj(i).lt.0.0) rxxx=rseps(1,jtime)/100.
        rxxxf=rxxx
        if (rzeroj(i).eq.0.0) then
             rxxx=1./r1sdry(i)
             rxxxf=r1sdry(i)/r2sdry(i)
        endif
        if (rxxx.le.0.0) then
            rxxx=rcentr
            rxxxf=rxxx
        endif
!
        ysiwant=sizeroj(i)
        call setfp(ysiwant,xpsfp)
        call setpp(ysiwant,xpspp)
        do 1230 j=1,kpcurn
          if (j.le.kppcur) then
            xjj=xpspp(j)
            alipc(nj,j)=rxxx*fwtxxx*xjj
          else
            xjj=xpsfp(j-kppcur)
            alipc(nj,j)=fwtxxx/rxxxf*xjj
          endif
 1230   continue
        nnow=kpcurn
        if (kedgep.gt.0) then
          nnow=nnow+1
          siedge=(sizeroj(i)-pe_psin)/pe_width
          alipc(nj,nnow)=fwtxxx*rxxx/cosh(siedge)**2/pe_width/sidif
        endif
        if (kedgef.gt.0) then
          nnow=nnow+1
          siedge=(sizeroj(i)-fe_psin)/fe_width
          alipc(nj,nnow)=fwtxxx/rxxxf/cosh(siedge)**2/fe_width/sidif
        endif
        xrsp(nj)=fwtxxx*vzeroj(i)*darea*pasmat(jtime)/carea
       enddo
      endif
!-----------------------------------------------------------------------
!--  constraint on betan by successive iterations                     --
!-----------------------------------------------------------------------
      if (fbetan.gt.0.0) then
        brspmin=max(ten24,abs(brsp(nfcoil+jbeta)))
        fwtxxn=fwtxxb*1000./brspmin
        nj=nj+1
        do 13230 j=1,kpcurn
          alipc(nj,j)=0.0
13230   continue
        calpao=brsp(nfcoil+jbeta)
        alipc(nj,jbeta)=fwtxxn
        if (initc.ge.jwantm) then
          xrsp(nj)=fwtxxn/abs(betatn)*fbetan*calpao
        else
          xrsp(nj)=fwtxxn*calpao
        endif
      endif
!-----------------------------------------------------------------------
!--  constraint on li successive iterations                           --
!-----------------------------------------------------------------------
      if (fli.gt.0.0) then
        brspmin=max(ten24,abs(brsp(nfcoil+jli)))
        fwtxxi=fwtxli*1000./brspmin
        nj=nj+1
        do 13240 j=1,kpcurn
          alipc(nj,j)=0.0
13240   continue
        cgamao=brsp(nbase+jli)
        alipc(nj,kppcur+jli)=fwtxxi
        if (initc.ge.jwantm) then
         if (cgamao.lt.0.0) then
          xrsp(nj)=fwtxxi*ali(jtime)/fli*cgamao
         else
          xrsp(nj)=fwtxxi/ali(jtime)*fli*cgamao
         endif
        else
          xrsp(nj)=fwtxxi*cgamao
        endif
      endif
!-----------------------------------------------------------------------
!-- constraints on P' and FF'                                         --
!-----------------------------------------------------------------------
      if (kcalpa.gt.0) then
        do 1240 j=1,kcalpa
        nj=nj+1
        do 1236 i=1,kppcur
          alipc(nj,i)=calpa(i,j)*fwtxxx
 1236   continue
        do 1238 i=kppcur+1,kpcurn
          alipc(nj,i)=0.
 1238   continue
        nnow=kpcurn
        if (kedgep.gt.0) then
          nnow=nnow+1
          alipc(nj,nnow)=calpa(kppcur+1,j)*fwtxxx
        endif
        if (kedgef.gt.0) then
          nnow=nnow+1
          alipc(nj,nnow)=0.0
        endif
        xrsp(nj)=xalpa(j)*darea*fwtxxx
 1240   continue
      endif
!
      if (kcgama.gt.0) then
        do 1250 j=1,kcgama
        nj=nj+1
        do 1244 i=1,kppcur
          alipc(nj,i)=0.0
 1244   continue
        do 1246 i=kppcur+1,kpcurn
          alipc(nj,i)=cgama(i-kppcur,j)*fwtxxx
 1246   continue
        nnow=kpcurn
        if (kedgep.gt.0) then
          nnow=nnow+1
          alipc(nj,nnow)=0.0
        endif
        if (kedgef.gt.0) then
          nnow=nnow+1
          alipc(nj,nnow)=cgama(kffcur+1,j)*fwtxxx
        endif
        xrsp(nj)=xgama(j)*darea*fwtxxx
        if (kfffnc.eq.6) then
           xrsp(nj)=xrsp(nj)/twopi/tmu
        endif
 1250   continue
      endif
!
      nnn=1
      kknow=kcalpa+kcgama
      if (fwtcur.gt.0.0) kknow=kknow+1
      if (fwtqa.gt.0.0) kknow=kknow+1
      if (kzeroj.gt.0) kknow=kknow+kzeroj
      nownow=kpcurn
      if (kedgep.gt.0) nownow=nownow+1
      if (kedgef.gt.0) nownow=nownow+1
      ncrsp = 0
      if (kknow.lt.nownow) then
        nzzzz = 0
        call ppcnst(ncrsp,crsp,z,nzzzz)
        call ffcnst(ncrsp,crsp,z,nzzzz)
      endif
      if (ncrsp .le. 0) then
        call sdecm(alipc,npcur3,nj,nownow,xrsp,npcur3,nnn,wlipc,work,ier)
        if (ier.eq.129) then
          kerror = 1
          call errctrl_msg('currnt','sdecm failed to converge (location 1)')
          return
        end if
        cond=ier
        toler=1.0e-06_dp*wlipc(1)
        do i=1,nownow
          t=0.0
          if (wlipc(i).gt.toler) t=xrsp(i)/wlipc(i)
          work(i)=t
        end do
        do i=1,nownow
          brsp(nfcoil+i)=0.0
          do j=1,nownow
            brsp(nfcoil+i)=brsp(nfcoil+i)+alipc(i,j)*work(j)
          end do
        end do

      else
        do j=1,nj
            b(j) = xrsp(j)
        enddo
        call dgglse(nj,nownow,ncrsp,alipc,npcur3,crsp,4*(npcurn-2)+6+ &
                   npcurn*npcurn,b,z,xrsp,work,nrsma2,info,condno)
        if (info.gt.0) then ! special hack to info in dgglse
          kerror = 1
          write(tmpstr,'(a,i4,a,i4,a)') 'A(',info,',',info,')=0 in dgglse, divide by zero.'
          call errctrl_msg('currnt',tmpstr)
          return
        else if (info.lt.0) then
          kerror = 1
          call errctrl_msg('currnt','calling argument in dgglse was bad')
          return
        endif
        do i=1,nownow
          brsp(nfcoil+i)=xrsp(i)
        enddo
      endif
      nownow=kpcurn
      if (kedgep.gt.0) then
         nownow=nownow+1
         pedge=brsp(nfcoil+nownow)
      endif
      if (kedgef.gt.0) then
         nownow=nownow+1
         f2edge=brsp(nfcoil+nownow)
      endif
!-----------------------------------------------------------------------
!-- update total plasma current if needed to                          --
!-----------------------------------------------------------------------
      if (abs(fwtcur).le.1.e-30_dp .and. nqwant.gt.0) then
        cm=0.0
        do 4999 n=nfcoil+1,nfnpcr
         cm=cm+brsp(n)*fgowpc(n-nfcoil)
 4999   continue
        if (kedgep.gt.0) then
           cm=cm+pedge*fgowpe
        endif
        if (kedgef.gt.0) then
           cm=cm+f2edge*fgowfe
        endif
        cpasma(jtime)=cm
        pasmat(jtime)=cm
      endif
!
 1800 continue
      tcurrt=0.0
      tcurrp=0.0
      tcurrtpp=0.0
      do 2000 i=1,nw
      do 2000 j=1,nh
        kk=(i-1)*nh+j
        pcurrt(kk)=0.0
        pcurrtpp(kk)=pcurrt(kk)
        if (icutfp.eq.0) then
          if ((xpsi(kk).lt.0.0).or.(xpsi(kk).gt.1.0)) go to 2000
          pcurrt(kk)=rgrid(i)*ppcurr(xpsi(kk),kppcur)
          pcurrtpp(kk)=pcurrt(kk)
          pcurrt(kk)=pcurrt(kk) &
                     +fpcurr(xpsi(kk),kffcur)/rgrid(i)
          pcurrt(kk)=pcurrt(kk)*www(kk)
          pcurrtpp(kk)=pcurrtpp(kk)*www(kk)
        else
          if ((xpsi(kk).ge.0.0).and.(xpsi(kk).le.1.0)) &
            pcurrt(kk)=rgrid(i)*ppcurr(xpsi(kk),kppcur)
          pcurrtpp(kk)=pcurrt(kk)
          upsi=xpsi(kk)*xpsimin
          if ((upsi.ge.0.0).and.(upsi.le.1.0)) &
           pcurrt(kk)=pcurrt(kk)+fpcurr(xpsi(kk),kffcur)/rgrid(i)
          pcurrt(kk)=pcurrt(kk)*zero(kk)
          pcurrtpp(kk)=pcurrtpp(kk)*zero(kk)
        endif
        tcurrt=tcurrt+pcurrt(kk)
        tcurrtpp=tcurrtpp+pcurrtpp(kk)
 2000 continue
      tcurrtffp=tcurrt-tcurrtpp
      if ((nitett.le.1).and.(icinit.lt.0)) then
        cratio=1.0
        return
      endif
      if (.not.fixpp) then
      cratio=cpasma(jtime)/tcurrt
      if (abs(cpasma(jtime)).le.1.e-3_dp) cratio=1.0
      cratio_ext = cratio * cratio_ext
      cratiop_ext = cratio_ext
      cratiof_ext = cratio_ext
      do kk=1,nwnh
        pcurrt(kk)=pcurrt(kk)*cratio
        pcurrtpp(kk)=pcurrtpp(kk)*cratio
        if ((xpsi(kk).ge.0.0).and.(xpsi(kk).le.1.0)) then
           tcurrp=tcurrp+pcurrt(kk)
        endif
      enddo
      do i=nfcoil+1,nfnpcr
        brsp(i)=cratio*brsp(i)
      enddo
      if (npsi_ext > 0) then
        prbdry=prbdry*cratio*cratio
      endif
      if (kedgep.gt.0) then
         pedge=pedge*cratio
      endif
      if (kedgef.gt.0) then
         f2edge=f2edge*cratio
      endif
      else
      cratio=1.0
      cratiop_ext = 1.0
      cratiof = (cpasma(jtime)-tcurrtpp)/tcurrtffp
      cratiof_ext = cratiof * cratiof_ext
      pcurrt(1:nwnh)=pcurrtpp(1:nwnh)+(pcurrt(1:nwnh)-pcurrtpp(1:nwnh))*cratiof
      do i=nbase+kppcur+1,nfnpcr
        brsp(i)=cratiof*brsp(i)
      enddo
      if (kedgef.gt.0) then
         f2edge=f2edge*cratiof
      endif
      endif
!----------------------------------------------------------------------------
!-- rigid vertical shift correction ?                                      --
!----------------------------------------------------------------------------
      if (ifitdelz.eq.3) then
      if (fitdelz.and.nitett.ge.ndelzon) then
      call sets2d(psi,c,rgrid,nw,bkx,lkx,zgrid,nh,bky,lky,wk,ier)
      cdelznow=cdelz(nitett-1)/100.
      cdeljsum=0.0
      do i=1,nw
      do j=1,nh
        kk=(i-1)*nh+j
        call seva2d(bkx,lkx,bky,lky,c,rgrid(i),zgrid(j),pds,ier,n333)
!sri-feb09
!        cdeljnow=cdelznow*pds(3)*rdjdz(kk)
        cdeljnow=cdelznow*rdjdz(kk)
        pcurrt(kk)=pcurrt(kk)+cdeljnow
        cdeljsum=cdeljsum+abs(cdeljnow)
      enddo
      enddo
      cdeljsum=abs(cdeljsum/tcurrp)
      endif
      endif
      return
!
 2100 continue
      return
!
 3100 continue
!------------------------------------------------------------------------
!--  GAQ type current profile                                          --
!------------------------------------------------------------------------
      if (kvtor.eq.11) then
        n1set=1
        ypsi=0.5_dp
        pres0=prcur4(n1set,ypsi,kppcur)
        prew0=pwcur4(n1set,ypsi,kwwcur)
        n1set=0
      endif
      tcurrt=0.0
      do 3300 i=1,nw
      do 3300 j=1,nh
        kk=(i-1)*nh+j
        pcurrt(kk)=0.0
        if ((xpsi(kk).lt.0.0).or.(xpsi(kk).gt.1.0)) go to 3300
        rdiml=rgrid(i)/rzero
        pp0       = (1.-xpsi(kk)**enp)**emp*(1.-gammap)+gammap
        pcurrt(kk)=  rbetap/rdiml*((1.-xpsi(kk)**enf)**emf &
                *(1.-gammaf)+gammaf)
!-----------------------------------------------------------------
!--  toroidal rotation ?                                        --
!-----------------------------------------------------------------
        if (kvtor.eq.0) then
         pcurrt(kk)=pcurrt(kk)+pp0*rdiml
        elseif (kvtor.ge.1 .and. kvtor.le.3) then
         ppw= (1.-xpsi(kk)**enw)**emw*(1.-gammaw)+gammaw
         ppw=rbetaw*ppw*rgrvt(i)
         pcurrt(kk)=pcurrt(kk)+(pp0+ppw)*rdiml
        elseif (kvtor.eq.11) then
         ypsi=xpsi(kk)
         pres0=prcur4(n1set,ypsi,kppcur)
         prew0=pwcur4(n1set,ypsi,kwwcur)
         if (abs(pres0).gt.1.e-10_dp) then
          pwop0=prew0/pres0
          ptop0=exp(pwop0*rgrvt(i))
         else
          ptop0=1.0
          pwop0=0.0
         endif
         pp0=pp0*(1.-pwop0*rgrvt(i))
         ppw= (1.-xpsi(kk)**enw)**emw*(1.-gammaw)+gammaw
         ppw=rbetaw*ppw*rgrvt(i)
         pcurrt(kk)=pcurrt(kk)+(pp0+ppw)*rdiml*ptop0
        endif
        pcurrt(kk)=pcurrt(kk)*www(kk)
        tcurrt=tcurrt+pcurrt(kk)
 3300 continue
      if (abs(tcurrt).le.1.0e-10_dp) then
        kerror=1
        call errctrl_msg('currnt','abs(tcurrt) <= 1.0e-10')
        return
      endif
      cratio=cpasma(jtime)/tcurrt
      do kk=1,nwnh
        pcurrt(kk)=pcurrt(kk)*cratio
      enddo
      dfsqe=0.5_dp
      ddpsi=1./real(nw-1,dp)
      sinow=0.0
      do i=2,nw-1
        sinow=sinow+ddpsi
        dfsqe=dfsqe+((1.-sinow**enf)**emf*(1.-gammaf)+gammaf)
      enddo
      dfsqe=dfsqe*2.*twopi/tmu*rzero*cratio*rbetap/darea*ddpsi
      if (nitett.gt.1) go to 4100
      if (icurrt.ne.2.and.icurrt.ne.5) go to 4100
      if (fwtbp.le.0.0) go to 4100
      do i=1,kppcur
        brsp(nfcoil+i)=0.0
      enddo
      do i=1,kffcur
        brsp(nbase+i)=0.0
      enddo
      brsp(nfcoil+1)=cratio/rzero
      brsp(nbase+1)=cratio*rbetap*rzero
      brsp(nfcoil+2)=-brsp(nfcoil+1)
      brsp(nbase+2)=-brsp(nbase+1)
 4100 continue
      return
!---------------------------------------------------------------------
!--  toroidal rotation, node points, bases :  ICURRT=5              --
!---------------------------------------------------------------------
 5100 continue
      if ((nitett.le.1).and.(icinit.gt.0)) go to 3100
      nnn=1
      call green(nnn,jtime,nitett)
      if ((nitett.le.1).and.(icinit.lt.0)) go to 5800
      if (iconvr.ne.3) go to 5800
!----------------------------------------------------------------------
!--  Adjust current profile to keep q(0), I, J(1), and others fixed  --
!----------------------------------------------------------------------
      nj=0
      if (fwtqa.gt.0.0) then
        nj=nj+1
        do 5210 j=1,kppcur
          alipc(nj,j)=rqajx(j)*fwtqa
 5210   continue
        do j=1,kffcur
          alipc(nj,kppcur+j)=rqafx(j)*fwtqa
        enddo
        if (kvtor.gt.0) then
          do j=1,kwwcur
            alipc(nj,kpcurn+j)=rqawx(j)*fwtqa
          enddo
        endif
        xrsp(nj)=fwtqa/qvfit*pasmat(jtime)/abs(pasmat(jtime))
      endif
!
      if (fwtcur.gt.0.0) then
      fwtcux=fwtcur
        nj=nj+1
        do 5220 j=1,kpcurn
          alipc(nj,j)=fwtcux*fgowpc(j)
 5220   continue
        if (kvtor.gt.0) then
          do j=kpcurn+1,kwcurn
            alipc(nj,j)=fwtcux*fgowpc(j)
          enddo
        endif
        xrsp(nj)=fwtcux*pasmat(jtime)
      endif
!-----------------------------------------------------------------------
!--  constraints on q at psiwant by successive iterations             --
!-----------------------------------------------------------------------
      if (nqwant.gt.0)   then
        do 15226 i=1,nqwant
        nj=nj+1
        if (initc.ge.jwantm) then
          fwtqqq=fwtxxq/0.001_dp/pasmsw(i)
          do 15220 j=1,kpcurn
            alipc(nj,j)=fwtqqq*fgowsw(j,i)
15220     continue
          if (kvtor.gt.0) then
            do j=kpcurn+1,kwcurn
              alipc(nj,j)=fwtqqq*fgowsw(j,i)
            enddo
          endif
          xrsp(nj)=fwtqqq/fqsiw(i)*qsiw(i)*pasmsw(i)
        endif
15226   continue
      endif
!-----------------------------------------------------------------------
!--  J at PSIWANT constraint                                          --
!-----------------------------------------------------------------------
      brspmin=max(ten24,abs(brsp(nfcoil+1)))
      fwtxxx=fwtxxj*1000./brspmin
      if (kzeroj.gt.0) then
       do i=1,kzeroj
        nj=nj+1
!----------------------------------------------------------------------
!--  local or flux surface averaged J constraint ?                   --
!----------------------------------------------------------------------
        if (rzeroj(i).gt.0.0) rxxx=rzeroj(i)
        if (rzeroj(i).lt.0.0) rxxx=rseps(1,jtime)/100.
        rxxxf=rxxx
        if (rzeroj(i).eq.0.0) then
             rxxx=1./r1sdry(i)
             rxxxf=r1sdry(i)/r2sdry(i)
        endif
        if (rxxx.le.0.0) then
            rxxx=rcentr
            rxxxf=rxxx
        endif
!
        ysiwant=sizeroj(i)
        call setfp(ysiwant,xpsfp)
        call setpp(ysiwant,xpspp)
!---------------------------------------------------------------------------
!--  rotational term                                                      --
!---------------------------------------------------------------------------
        if (kvtor.gt.0) then
          call setpwp(ysiwant,xpspwp)
          if (kvtor.eq.2.or.kvtor.eq.3) then
            prew0=pwcurr(ysiwant,kwwcur)
            pres0=prcurr(ysiwant,kppcur)
            if (abs(pres0).gt.1.e-10_dp) then
               pwop0=prew0/pres0
            else
               pwop0=0.0
            endif
          endif
          if (rzeroj(i).gt.0.0) then
              rxx2=(rzeroj(i)/rvtor)**2-1.
              rxxw=rxx2*rzeroj(i)
              if (kvtor.eq.2) then
                rxxw=rxxw*(1.+pwop0*rxx2)
                rxxx=rxxx*(1.-0.5_dp*(pwop0*rxx2)**2)
              endif
              if (kvtor.eq.3) then
                pwp0r2=pwop0*rxx2
                ptop0=exp(pwp0r2)
                rxxw=rxxw*ptop0
                rxxx=rxxx*ptop0*(1.-pwp0r2)
              endif
          endif
          if (rzeroj(i).lt.0.0) then
              rxxw=rseps(1,jtime)/100.
              rxx2=(rxxw/rvtor)**2-1.
              rxxw=rxx2*rxxw
              if (kvtor.eq.2) then
                rxxw=rxxw*(1.+pwop0*rxx2)
                rxxx=rxxx*(1.-0.5_dp*(pwop0*rxx2)**2)
              endif
              if (kvtor.eq.3) then
                pwp0r2=pwop0*rxx2
                ptop0=exp(pwp0r2)
                rxxw=rxxw*ptop0
                rxxx=rxxx*ptop0*(1.-pwp0r2)
              endif
          endif
          if (rzeroj(i).eq.0.0) then
              rxx2=  r2wdry/rvtor**2-1.
              rxxw=  rxx2/r1sdry(i)
              if (kvtor.eq.2) then
                rxxw=rxxw+pwop0*r4wdry/r1sdry(i)
                rxxx=rxxx-0.5_dp*pwop0**2*r4wdry/r1sdry(i)
              endif
              if (kvtor.eq.3) then
                rxxx=(rpwdry-pwop0*rp2wdry)/r1sdry(i)
                rxxw=rp2wdry/r1sdry(i)
              endif
          endif
        endif
!
        do 5230 j=1,kwcurn
          if (j.le.kppcur) then
            xjj=xpspp(j)
            alipc(nj,j)=rxxx*fwtxxx*xjj
          elseif (j.le.kpcurn) then
            xjj=xpsfp(j-kppcur)
            alipc(nj,j)=fwtxxx/rxxxf*xjj
          elseif (kvtor.gt.0) then
            xjj=xpspwp(j-kpcurn)
            alipc(nj,j)=rxxw*fwtxxx*xjj
          endif
 5230   continue
        xrsp(nj)=fwtxxx*vzeroj(i)*darea*pasmat(jtime)/carea
       enddo
      endif
!-----------------------------------------------------------------------
!--  constraint on betan by successive iterations                     --
!-----------------------------------------------------------------------
      if (fbetan.gt.0.0) then
        brspmin=max(ten24,abs(brsp(nfcoil+jbeta)))
        fwtxxn=fwtxxb*1000./brspmin
        nj=nj+1
        do 15230 j=1,kpcurn
          alipc(nj,j)=0.0
15230   continue
        calpao=brsp(nfcoil+jbeta)
        alipc(nj,jbeta)=fwtxxn
        if (initc.ge.jwantm) then
          xrsp(nj)=fwtxxn/abs(betatn)*fbetan*calpao
        else
          xrsp(nj)=fwtxxn*calpao
        endif
      endif
!-----------------------------------------------------------------------
!--  constraint on li successive iterations                           --
!-----------------------------------------------------------------------
      if (fli.gt.0.0) then
        brspmin=max(ten24,abs(brsp(nfcoil+jli)))
        fwtxxi=fwtxli*1000./brspmin
        nj=nj+1
        do 15240 j=1,kpcurn
          alipc(nj,j)=0.0
15240   continue
        cgamao=brsp(nbase+jli)
        alipc(nj,kppcur+jli)=fwtxxi
        if (initc.ge.jwantm) then
         if (cgamao.lt.0.0) then
          xrsp(nj)=fwtxxi*ali(jtime)/fli*cgamao
         else
          xrsp(nj)=fwtxxi/ali(jtime)*fli*cgamao
         endif
        else
          xrsp(nj)=fwtxxi*cgamao
        endif
      endif
!-----------------------------------------------------------------------
!-- constraints on P' and FF'                                         --
!-----------------------------------------------------------------------
      if (kcalpa.gt.0) then
        do 5240 j=1,kcalpa
        nj=nj+1
        do 5236 i=1,kppcur
          alipc(nj,i)=calpa(i,j)*fwtxxx
 5236   continue
        do 5238 i=kppcur+1,kpcurn
          alipc(nj,i)=0.
 5238   continue
        xrsp(nj)=xalpa(j)*darea*fwtxxx
 5240   continue
      endif
!
      if (kcgama.gt.0) then
        do 5250 j=1,kcgama
        nj=nj+1
        do 5244 i=1,kppcur
          alipc(nj,i)=0.0
 5244   continue
        do 5246 i=kppcur+1,kpcurn
          alipc(nj,i)=cgama(i-kppcur,j)*fwtxxx
 5246   continue
        xrsp(nj)=xgama(j)*darea*fwtxxx
 5250   continue
      endif
!
      if (kcomega.gt.0) then
        do j=1,kcomega
          nj=nj+1
          do i=1,kwcurn
            alipc(nj,i)=comega(i,j)*fwtxxx
          enddo
          xrsp(nj)=xomega(j)*darea*fwtxxx
        enddo
      endif
!
      if (nj.le.0) go to 5800
      nnn=1
      call sdecm(alipc,npcur3,nj,kwcurn,xrsp,npcur3,nnn,wlipc,work,ier)
      if (ier.eq.129) then
        kerror = 1
        call errctrl_msg('currnt','sdecm failed to converge (location 2)')
        return
      end if
      cond=ier
      toler=1.0e-06_dp*wlipc(1)
      do 5570 i=1,kwcurn
        t=0.0
        if (wlipc(i).gt.toler) t=xrsp(i)/wlipc(i)
        work(i)=t
 5570 continue
      do 5575 i=1,kwcurn
        brsp(nfcoil+i)=0.0
        do 5575 j=1,kwcurn
          brsp(nfcoil+i)=brsp(nfcoil+i)+alipc(i,j)*work(j)
 5575 continue
!-----------------------------------------------------------------------
!-- update total plasma current if needed to                          --
!-----------------------------------------------------------------------
      if (abs(fwtcur).le.1.e-30_dp.and.nqwant.gt.0) then
        cm=0.0
        do 5999 n=nfcoil+1,nfnwcr
         cm=cm+brsp(n)*fgowpc(n-nfcoil)
 5999   continue
        cpasma(jtime)=cm
        pasmat(jtime)=cm
      endif
!
 5800 continue
      tcurrt=0.0
      tcurrp=0.0
      do 6000 i=1,nw
      do 6000 j=1,nh
        kk=(i-1)*nh+j
        pcurrt(kk)=0.0
        pcurrw(kk)=0.0
        if (icutfp.eq.0) then
!----------------------------------------------------------------------
!--  no attached current                                             --
!----------------------------------------------------------------------
          if ((xpsi(kk).lt.0.0).or.(xpsi(kk).gt.1.0)) go to 6000
          pp0=ppcurr(xpsi(kk),kppcur)
          pcurrt(kk)=fpcurr(xpsi(kk),kffcur)/rgrid(i)
          if (kvtor.eq.0) then
             pcurrt(kk)=pcurrt(kk)+pp0*rgrid(i)
          else
!-----------------------------------------------------------------------
!--      rotation                                                     --
!-----------------------------------------------------------------------
             rdimw=rgrid(i)/rvtor
             ppw=pwpcur(xpsi(kk),kwwcur)
             ppw=ppw*(rdimw**2-1.)
             pcurrw(kk)=ppw*rgrid(i)
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
               pcurrw(kk)=pcurrw(kk)*(1.+pwp0r2)
               pp0=pp0*(1.-0.5_dp*pwp0r2**2)
             endif
             if (kvtor.eq.3) then
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
               pcurrw(kk)=pcurrw(kk)*ptop0
               pp0=pp0*ptop0*(1.-pwp0r2)
             endif
             pcurrt(kk)=pcurrt(kk)+pp0*rgrid(i)+pcurrw(kk)
          endif
          pcurrt(kk)=pcurrt(kk)*www(kk)
          pcurrw(kk)=pcurrw(kk)*www(kk)
        else
!------------------------------------------------------------------------
!--  attached current                                                  --
!------------------------------------------------------------------------
          if ((xpsi(kk).ge.0.0).and.(xpsi(kk).le.1.0)) then
             pp0=ppcurr(xpsi(kk),kppcur)
             if (kvtor.eq.0) then
               pcurrt(kk)=rgrid(i)*pp0
             else
!------------------------------------------------------------------------
!--        rotation                                                    --
!------------------------------------------------------------------------
               rdimw=rgrid(i)/rvtor
               ppw=pwpcur(xpsi(kk),kwwcur)
               ppw=ppw*(rdimw**2-1.)
               pcurrw(kk)=ppw*rgrid(i)
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
                 pcurrw(kk)=pcurrw(kk)*(1.+pwp0r2)
                 pp0=pp0*(1.-0.5_dp*pwp0r2**2)
               endif
               if (kvtor.eq.3) then
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
                 pcurrw(kk)=pcurrw(kk)*ptop0
                 pp0=pp0*ptop0*(1.-pwp0r2)
               endif
               pcurrt(kk)=pp0*rgrid(i)+pcurrw(kk)
             endif
          endif
          upsi=xpsi(kk)*xpsimin
          if ((upsi.ge.0.0).and.(upsi.le.1.0)) &
            pcurrt(kk)=pcurrt(kk)+fpcurr(xpsi(kk),kffcur)/rgrid(i)
          pcurrt(kk)=pcurrt(kk)*zero(kk)
        endif
        tcurrt=tcurrt+pcurrt(kk)
 6000 continue
      if ((nitett.le.1).and.(icinit.lt.0)) then
        cratio=1.0
        return
      endif
!-------------------------------------------------------------------
!--  adjust to match total current                                --
!-------------------------------------------------------------------
      cratio=cpasma(jtime)/tcurrt
      do kk=1,nwnh
        pcurrt(kk)=pcurrt(kk)*cratio
        pcurrw(kk)=pcurrw(kk)*cratio
        if ((xpsi(kk).ge.0.0).and.(xpsi(kk).le.1.0)) then
           tcurrp=tcurrp+pcurrt(kk)
        endif
      enddo
      do i=nfcoil+1,nfcoil+kwcurn
        brsp(i)=cratio*brsp(i)
      enddo
      return

      end
