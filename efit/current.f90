!**********************************************************************
!!
!>    currnt computes the current density on the r-z mesh.
!!
!!    @param iter : current profile (outer) loop iteration index
!!    @param jtime : time index
!!    @param nitett : total iteration index (current+equilibirum loops)
!!    @param kerror : error flag
!!
!*********************************************************************
      subroutine currnt(iter,jtime,nitett,kerror)
      use commonblocks,only: c,wk,bkx,bky
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      real*8 fpcurr,ppcurr,prcurr,pwcurr,prcur4,pwcur4,pwpcur

      integer*4, intent(in) :: iter,jtime,nitett
      integer*4, intent(out) :: kerror
      integer*4 i,j,kk,n,nj,ier,info,initc,kknow,n1set,ncrsp,nnow, &
                nownow,npcur2
      real*8 aqax,brspmin,calpao,cdeljnow,cdelznow,cgamao,cj0,cm,cwant0, &
             cwant1,ddpsi,fwtcux,fwtqqq,fwtxxi,fwtxxn,fwtxxx,pp0,ppw, &
             pres0,prew0,ptop0,pwop0,pwp0r2,rdiml,rdimw,rxx2,rxxw,rxxx, &
             rxxxf,siedge,sinow,t,tcurrt,tcurrtpp,tcurrtffp,toler,tz, &
             upsi,xjj,ypsi,ysiwant
      character(len=128) tmpstr
      real*8 pds(6)
      real*8 alipc(npcurn*2,nwcurn),xpspp(nppcur),xpsfp(nffcur)
      real*8 wlipc(nwcurn),work(nwcurn*2),xrsp(npcurn*2),xpspwp(nwwcur)
      real*8 crsp(4*(npcurn-2)+6+npcurn*npcurn,nrsmat)
      real*8 b(nrsmat),z(4*(npcurn-2)+6+npcurn*npcurn)
      integer*4, parameter :: nzzzz=0,nnn=1
      real*8, parameter :: ten24=1.e4_dp
      data initc/0/

      kerror = 0
      npcur2=npcurn*2
      initc=initc+1

      if(ivacum.eq.1) return
      if((nitett.le.1).and.(icinit.eq.1)) return
      GAQ: if (((icinit.gt.0).and.(iconvr.ne.3).and.(iter.le.1)).or.(icurrt.eq.4).or. &
        (((icurrt.eq.2).or.(icurrt.eq.5)).and.(nitett.le.1).and.(icinit.gt.0))) then
!-----------------------------------------------------------------------
!--    GAQ type current profile                                       --
!-----------------------------------------------------------------------
       if (kvtor.eq.11) then
         n1set=1
         ypsi=0.5_dp
         pres0=prcur4(n1set,ypsi)
         prew0=pwcur4(n1set,ypsi)
         n1set=0
       endif
       tcurrt=0.0
       do i=1,nw
        do j=1,nh
         kk=(i-1)*nh+j
         pcurrt(kk)=0.0
         if((xpsi(kk).lt.0.0).or.(xpsi(kk).gt.1.0)) cycle
         rdiml=rgrid(i)/rzero
         pp0=(1.-xpsi(kk)**enp)**emp*(1.-gammap)+gammap
         pcurrt(kk)=rbetap/rdiml*((1.-xpsi(kk)**enf)**emf &
                    *(1.-gammaf)+gammaf)
!-----------------------------------------------------------------
!--      toroidal rotation                                      --
!-----------------------------------------------------------------
         if (kvtor.eq.0) then
           pcurrt(kk)=pcurrt(kk)+pp0*rdiml
         elseif (kvtor.ge.1 .and. kvtor.le.3) then
           ppw=(1.-xpsi(kk)**enw)**emw*(1.-gammaw)+gammaw
           ppw=rbetaw*ppw*rgrvt(i)
           pcurrt(kk)=pcurrt(kk)+(pp0+ppw)*rdiml
         elseif (kvtor.eq.11) then
           ypsi=xpsi(kk)
           pres0=prcur4(n1set,ypsi)
           prew0=pwcur4(n1set,ypsi)
           if (abs(pres0).gt.1.e-10_dp) then
             pwop0=prew0/pres0
             ptop0=exp(pwop0*rgrvt(i))
           else
             ptop0=1.0
             pwop0=0.0
           endif
           pp0=pp0*(1.-pwop0*rgrvt(i))
           ppw=(1.-xpsi(kk)**enw)**emw*(1.-gammaw)+gammaw
           ppw=rbetaw*ppw*rgrvt(i)
           pcurrt(kk)=pcurrt(kk)+(pp0+ppw)*rdiml*ptop0
         endif
         pcurrt(kk)=pcurrt(kk)*www(kk)
         tcurrt=tcurrt+pcurrt(kk)
        enddo
       enddo
       if (abs(tcurrt).le.1.0e-10_dp) then
         ! there is likely a problem with the weights (www)
         kerror=1
         call errctrl_msg('currnt','abs(tcurrt) <= 1.0e-10')
         return
       endif
       cratio=ipmhd(jtime)/tcurrt
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
       if (nitett.gt.1) return
       if (icurrt.ne.2.and.icurrt.ne.5) return
       if (fwtbp.le.0.0) return
       do i=1,kppcur
         brsp(nfsum+i)=0.0
       enddo
       do i=1,kffcur
         brsp(nbase+i)=0.0
       enddo
       brsp(nfsum+1)=cratio/rzero
       brsp(nbase+1)=cratio*rbetap*rzero
       brsp(nfsum+2)=-brsp(nfsum+1)
       brsp(nbase+2)=-brsp(nbase+1)
       return
      endif GAQ
!
      select case (icurrt)
      case default ! 1
!------------------------------------------------------------------------------
!--    uniform current for Solovev equilibrium                               --
!------------------------------------------------------------------------------
       tcurrt=0.0
       do i=1,nw
         do j=1,nh
           kk=(i-1)*nh+j
           pcurrt(kk)=0.0
           pcurrw(kk)=0.0
           if ((xpsi(kk).lt.0.0).or.(xpsi(kk).gt.1.0)) cycle
           rdiml=rgrid(i)/srma
           pcurrt(kk)=sbeta*rdiml+2.*salpha/rdiml
           if (kvtor.gt.0) then
             pcurrw(kk)=sbetaw*rdiml*(rdiml**2-1.)
             pcurrt(kk)=pcurrt(kk)+pcurrw(kk)
           endif
           pcurrt(kk)=pcurrt(kk)*www(kk)
           pcurrw(kk)=pcurrw(kk)*www(kk)
           tcurrt=tcurrt+pcurrt(kk)
         enddo
       enddo
       cratio=ipmhd(jtime)/tcurrt
       do kk=1,nwnh
         pcurrt(kk)=pcurrt(kk)*cratio
         pcurrw(kk)=pcurrw(kk)*cratio
       enddo
       cj0=cratio/darea/2.
       fcentr=twopi*2.*saaa**2*cj0
       fbrdy=fcentr*sqrt(1.-4.*salpha)
       bcentr(jtime)=fbrdy*tmu/rcentr
!
      case (2)
!----------------------------------------------------------------------
!--    polynomial or spline current profile                          --
!----------------------------------------------------------------------
       init_current: if ((nitett.gt.1).or.(icinit.ge.0)) then
       eq_mode: if (iconvr.eq.3) then
         constrain_profs: if (kcgama.gt.0.or.kcalpa.gt.0) then
!----------------------------------------------------------------------
!--      Adjust current profile to keep q(0), I, J(1), and others fixed
!----------------------------------------------------------------------
           if(fwtqa.gt.0.0 .or. fwtcur.gt.0.0 .or. nqwant.gt.0) &
             call green(nnn,jtime,nitett)
           nj=0
           if (fwtqa.gt.0.0) then
             nj=nj+1
             do j=1,kppcur
               alipc(nj,j)=rqajx(j)*fwtqa
             enddo
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
             xrsp(nj)=fwtqa/qvfit*ipmeas(jtime)/abs(ipmeas(jtime))
           endif
!
           if (fwtcur.gt.0.0) then
             fwtcux=fwtcur
             nj=nj+1
             do j=1,kpcurn
               alipc(nj,j)=fwtcux*fgowpc(j)
             enddo
             nnow=kpcurn
             if (kedgep.gt.0) then
               nnow=nnow+1
               alipc(nj,nnow)=fwtcux*fgowpe
             endif
             if (kedgef.gt.0) then
               nnow=nnow+1
               alipc(nj,nnow)=fwtcux*fgowfe
             endif
             xrsp(nj)=fwtcux*ipmeas(jtime)
           endif
!-----------------------------------------------------------------------
!--        constraints on q at psiwant by successive iterations       --
!-----------------------------------------------------------------------
           if (nqwant.gt.0) then
             do i=1,nqwant
               nj=nj+1
               if (initc.ge.jwantm) then
                 fwtqqq=fwtxxq/0.001_dp/pasmsw(i)
                 do j=1,kpcurn
                   alipc(nj,j)=fwtqqq*fgowsw(j,i)
                 enddo
                 xrsp(nj)=fwtqqq/fqsiw(i)*qsiw(i)*pasmsw(i)
               endif
             enddo
           endif
!-----------------------------------------------------------------------
!--        J at PSIWANT constraint                                    --
!-----------------------------------------------------------------------
           brspmin=max(ten24,abs(brsp(nfsum+1)))
           fwtxxx=fwtxxj*1000./brspmin
           if (kzeroj.gt.0) then
            do i=1,kzeroj
             nj=nj+1
!----------------------------------------------------------------------
!--          local or flux surface averaged J constraint ?           --
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
             do j=1,kpcurn
               if (j.le.kppcur) then
                 xjj=xpspp(j)
                 alipc(nj,j)=rxxx*fwtxxx*xjj
               else
                 xjj=xpsfp(j-kppcur)
                 alipc(nj,j)=fwtxxx/rxxxf*xjj
               endif
             enddo
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
             xrsp(nj)=fwtxxx*vzeroj(i)*darea*ipmeas(jtime)/carea
            enddo
           endif
!-----------------------------------------------------------------------
!--        constraint on betan by successive iterations               --
!-----------------------------------------------------------------------
           if (fbetan.gt.0.0) then
            brspmin=max(ten24,abs(brsp(nfsum+jbeta)))
             fwtxxn=fwtxxb*1000./brspmin
             nj=nj+1
             do j=1,kpcurn
               alipc(nj,j)=0.0
             enddo
             calpao=brsp(nfsum+jbeta)
             alipc(nj,jbeta)=fwtxxn
             if (initc.ge.jwantm) then
               xrsp(nj)=fwtxxn/abs(betatn)*fbetan*calpao
             else
               xrsp(nj)=fwtxxn*calpao
             endif
           endif
!-----------------------------------------------------------------------
!--        constraint on li successive iterations                     --
!-----------------------------------------------------------------------
           if (fli.gt.0.0) then
             brspmin=max(ten24,abs(brsp(nfsum+jli)))
             fwtxxi=fwtxli*1000./brspmin
             nj=nj+1
             do j=1,kpcurn
               alipc(nj,j)=0.0
             enddo
             cgamao=brsp(nbase+jli)
             alipc(nj,kppcur+jli)=fwtxxi
             if (initc.ge.jwantm) then
              if (cgamao.lt.0.0) then
               xrsp(nj)=fwtxxi*li(jtime)/fli*cgamao
              else
               xrsp(nj)=fwtxxi/li(jtime)*fli*cgamao
              endif
             else
              xrsp(nj)=fwtxxi*cgamao
             endif
           endif
!-----------------------------------------------------------------------
!--        constraints on P' and FF'                                  --
!-----------------------------------------------------------------------
           if (kcalpa.gt.0) then
             do j=1,kcalpa
               nj=nj+1
               do i=1,kppcur
                 alipc(nj,i)=calpa(i,j)*fwtxxx
               enddo
               do i=kppcur+1,kpcurn
                 alipc(nj,i)=0.
               enddo
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
             enddo
           endif
!
           if (kcgama.gt.0) then
             do j=1,kcgama
               nj=nj+1
               do i=1,kppcur
                 alipc(nj,i)=0.0
               enddo
               do i=kppcur+1,kpcurn
                 alipc(nj,i)=cgama(i-kppcur,j)*fwtxxx
               enddo
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
             enddo
           endif
!
           kknow=kcalpa+kcgama
           if (fwtcur.gt.0.0) kknow=kknow+1
           if (fwtqa.gt.0.0) kknow=kknow+1
           if (kzeroj.gt.0) kknow=kknow+kzeroj
           nownow=kpcurn
           if (kedgep.gt.0) nownow=nownow+1
           if (kedgef.gt.0) nownow=nownow+1
           ncrsp = 0
           if (kknow.lt.nownow) then
             call ppcnst(ncrsp,crsp,z,nzzzz)
             call ffcnst(ncrsp,crsp,z,nzzzz)
           endif
           if (ncrsp .le. 0) then
             call sdecm(alipc,npcur2,nj,nownow,xrsp,npcur2,nnn,wlipc, &
                        work,ier)
             if (ier.eq.129) then
               kerror=1
               call errctrl_msg('currnt', &
                                'sdecm failed to converge (location 1)')
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
               brsp(nfsum+i)=0.0
               do j=1,nownow
                 brsp(nfsum+i)=brsp(nfsum+i)+alipc(i,j)*work(j)
               end do
             end do
           else
             do j=1,nj
               b(j) = xrsp(j)
             enddo
             info=0
             call dgglse(int(nj,8),int(nownow,8),int(ncrsp,8),alipc, &
                         int(npcur2,8),crsp, &
                         int(4*(npcurn-2)+6+npcurn*npcurn,8), &
                         b,z,xrsp,work,int(nwcurn*2,8),int(info,8),condno)
             if (info.gt.0) then ! special hack to info in dgglse
               kerror=1
               write(tmpstr,'(a,i4,a,i4,a)') &
                 'A(',info,',',info,')=0 in dgglse, divide by zero.'
               call errctrl_msg('currnt',tmpstr)
               return
             else if (info.lt.0) then
               kerror=1
               call errctrl_msg('currnt', &
                                'calling argument in dgglse was bad')
               return
             endif
             do i=1,nownow
               brsp(nfsum+i)=xrsp(i)
             enddo
           endif
           nownow=kpcurn
           if (kedgep.gt.0) then
             nownow=nownow+1
             pedge=brsp(nfsum+nownow)
           endif
           if (kedgef.gt.0) then
             nownow=nownow+1
             f2edge=brsp(nfsum+nownow)
           endif
!-----------------------------------------------------------------------
!--        update total plasma current if needed to                   --
!-----------------------------------------------------------------------
           if (abs(fwtcur).le.1.e-30_dp .and. nqwant.gt.0) then
             cm=0.0
             do n=nfsum+1,nfnpcr
               cm=cm+brsp(n)*fgowpc(n-nfsum)
             enddo
             if (kedgep.gt.0) then
               cm=cm+pedge*fgowpe
             endif
             if (kedgef.gt.0) then
               cm=cm+f2edge*fgowfe
             endif
             ipmhd(jtime)=cm
             ipmeas(jtime)=cm
           endif
         endif constrain_profs
       endif eq_mode
       ! TODO: why does qenp value matter?
       if ((qenp.gt.0.0).and.(nitett.gt.1).and.(qvfit.gt.0.0).and. &
           (fwtqa.gt.0.0).and.(kedgep.le.0).and.(kedgef.le.0)) then
         tz = 0.0
         aqax=1.0/qvfit/(rqajtor*ppcurr(tz,kppcur) &
                        +rqaftor*fpcurr(tz,kffcur))
         aqax = abs(aqax)
         do i=1,kppcur
           brsp(nfsum+i)=aqax*brsp(nfsum+i)
         enddo
         do i=1,kffcur
           brsp(nbase+i)=aqax*brsp(nbase+i)
         enddo
         cwant0=ipmhd(jtime)-fgowpc(1)*brsp(nfsum+1)-fgowpc(kppcur+1) &
               *brsp(nbase+1)
         cwant1=0.0
         do i=2,kpcurn
           if(i.eq.kppcur+1) cycle
           cwant1=cwant1+fgowpc(i)*brsp(nfsum+i)
         enddo
         if (abs(cwant1).le.1.0e-10_dp) then
           kerror=1
           call errctrl_msg('currnt','abs(cwant1) <= 1.0e-10')
           return
         end if

         cwant1=cwant0/cwant1
         do i=2,kpcurn
           if(i.eq.kppcur+1) cycle
           brsp(nfsum+i)=cwant1*brsp(nfsum+i)
         enddo
!
       endif
       endif init_current
       tcurrt=0.0
       tcurrp=0.0
       tcurrtpp=0.0
       do i=1,nw
         do j=1,nh
           kk=(i-1)*nh+j
           pcurrt(kk)=0.0
           pcurrtpp(kk)=pcurrt(kk)
           if (icutfp.eq.0) then
             if ((xpsi(kk).lt.0.0).or.(xpsi(kk).gt.1.0)) cycle
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
         enddo
       enddo
       tcurrtffp=tcurrt-tcurrtpp
       if ((nitett.le.1).and.(icinit.lt.0)) then
         cratio=1.0
         return
       endif
       if (.not.fixpp) then
        cratio=ipmhd(jtime)/tcurrt
        if(abs(ipmhd(jtime)).le.1.e-3_dp) cratio=1.0
        cratio_ext = cratio * cratio_ext
        cratiop_ext = cratio_ext
        cratiof_ext = cratio_ext
        do kk=1,nwnh
          pcurrt(kk)=pcurrt(kk)*cratio
          pcurrtpp(kk)=pcurrtpp(kk)*cratio
          if((xpsi(kk).ge.0.0).and.(xpsi(kk).le.1.0)) &
            tcurrp=tcurrp+pcurrt(kk)
        enddo
        do i=nfsum+1,nfnpcr
          brsp(i)=cratio*brsp(i)
        enddo
        if(npsi_ext > 0) &
          prbdry=prbdry*cratio*cratio
        if(kedgep.gt.0) &
          pedge=pedge*cratio
        if(kedgef.gt.0) &
          f2edge=f2edge*cratio
       else
        cratio=1.0
        cratiop_ext = 1.0
        cratiof = (ipmhd(jtime)-tcurrtpp)/tcurrtffp
        cratiof_ext = cratiof * cratiof_ext
        pcurrt(1:nwnh)=pcurrtpp(1:nwnh)+(pcurrt(1:nwnh)-pcurrtpp(1:nwnh))*cratiof
        do i=nbase+kppcur+1,nfnpcr
          brsp(i)=cratiof*brsp(i)
        enddo
        if(kedgef.gt.0) &
          f2edge=f2edge*cratiof
       endif
!----------------------------------------------------------------------------
!--    rigid vertical shift correction ?                                   --
!----------------------------------------------------------------------------
       if (ifitdelz.eq.3) then
        if (fitdelz.and.nitett.ge.ndelzon) then
         !TODO: this is not normally needed here, should only be added
         !when necessary (need to figure out when that is...)
         !call sets2d(psi,c,rgrid,nw,bkx,lkx,zgrid,nh,bky,lky,wk,ier)
         cdelznow=cdelz(nitett-1)/100.
         cdeljsum=0.0
         do i=1,nw
          do j=1,nh
           kk=(i-1)*nh+j
!sri-feb09
!           call seva2d(bkx,lkx,bky,lky,c,rgrid(i),zgrid(j),pds,ier,n333)
!           cdeljnow=cdelznow*pds(3)*rdjdz(kk)
           cdeljnow=cdelznow*rdjdz(kk)
           pcurrt(kk)=pcurrt(kk)+cdeljnow
           cdeljsum=cdeljsum+abs(cdeljnow)
          enddo
         enddo
         cdeljsum=abs(cdeljsum/tcurrp)
        endif
       endif
!
      case (3)
       ! continue
!
      ! case (4) handled by preceeding if statement
      case (5)
!---------------------------------------------------------------------
!--    toroidal rotation, node points, bases :  ICURRT=5            --
!---------------------------------------------------------------------
       init_curr: if (((nitett.gt.1).or.(icinit.ge.0)).and.(iconvr.eq.3)) then
!----------------------------------------------------------------------
!--    Adjust current profile to keep q(0), I, J(1), and others fixed--
!----------------------------------------------------------------------
       if(fwtqa.gt.0.0 .or. fwtcur.gt.0.0 .or. nqwant.gt.0) &
         call green(nnn,jtime,nitett)
       nj=0
       if (fwtqa.gt.0.0) then
         nj=nj+1
         do j=1,kppcur
           alipc(nj,j)=rqajx(j)*fwtqa
         enddo
         do j=1,kffcur
           alipc(nj,kppcur+j)=rqafx(j)*fwtqa
         enddo
         if (kvtor.gt.0) then
           do j=1,kwwcur
             alipc(nj,kpcurn+j)=rqawx(j)*fwtqa
           enddo
         endif
         xrsp(nj)=fwtqa/qvfit*ipmeas(jtime)/abs(ipmeas(jtime))
       endif
!
       if (fwtcur.gt.0.0) then
         fwtcux=fwtcur
         nj=nj+1
         do j=1,kpcurn
           alipc(nj,j)=fwtcux*fgowpc(j)
         enddo
         if (kvtor.gt.0) then
           do j=kpcurn+1,kwcurn
             alipc(nj,j)=fwtcux*fgowpc(j)
           enddo
         endif
         xrsp(nj)=fwtcux*ipmeas(jtime)
       endif
!-----------------------------------------------------------------------
!--    constraints on q at psiwant by successive iterations           --
!-----------------------------------------------------------------------
       if (nqwant.gt.0) then
         do i=1,nqwant
           nj=nj+1
           if (initc.ge.jwantm) then
             fwtqqq=fwtxxq/0.001_dp/pasmsw(i)
             do j=1,kpcurn
               alipc(nj,j)=fwtqqq*fgowsw(j,i)
             enddo
             if (kvtor.gt.0) then
               do j=kpcurn+1,kwcurn
                 alipc(nj,j)=fwtqqq*fgowsw(j,i)
               enddo
             endif
             xrsp(nj)=fwtqqq/fqsiw(i)*qsiw(i)*pasmsw(i)
           endif
         enddo
       endif
!-----------------------------------------------------------------------
!--    J at PSIWANT constraint                                        --
!-----------------------------------------------------------------------
       brspmin=max(ten24,abs(brsp(nfsum+1)))
       fwtxxx=fwtxxj*1000./brspmin
       if (kzeroj.gt.0) then
        do i=1,kzeroj
         nj=nj+1
!----------------------------------------------------------------------
!--      local or flux surface averaged J constraint ?               --
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
!--      rotational term                                                  --
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
         do j=1,kwcurn
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
         enddo
         xrsp(nj)=fwtxxx*vzeroj(i)*darea*ipmeas(jtime)/carea
        enddo
       endif
!-----------------------------------------------------------------------
!--    constraint on betan by successive iterations                   --
!-----------------------------------------------------------------------
       if (fbetan.gt.0.0) then
         brspmin=max(ten24,abs(brsp(nfsum+jbeta)))
         fwtxxn=fwtxxb*1000./brspmin
         nj=nj+1
         do j=1,kpcurn
           alipc(nj,j)=0.0
         enddo
         calpao=brsp(nfsum+jbeta)
         alipc(nj,jbeta)=fwtxxn
         if (initc.ge.jwantm) then
           xrsp(nj)=fwtxxn/abs(betatn)*fbetan*calpao
         else
           xrsp(nj)=fwtxxn*calpao
         endif
       endif
!-----------------------------------------------------------------------
!--    constraint on li successive iterations                         --
!-----------------------------------------------------------------------
       if (fli.gt.0.0) then
         brspmin=max(ten24,abs(brsp(nfsum+jli)))
         fwtxxi=fwtxli*1000./brspmin
         nj=nj+1
         do j=1,kpcurn
           alipc(nj,j)=0.0
         enddo
         cgamao=brsp(nbase+jli)
         alipc(nj,kppcur+jli)=fwtxxi
         if (initc.ge.jwantm) then
          if (cgamao.lt.0.0) then
           xrsp(nj)=fwtxxi*li(jtime)/fli*cgamao
          else
           xrsp(nj)=fwtxxi/li(jtime)*fli*cgamao
          endif
         else
           xrsp(nj)=fwtxxi*cgamao
         endif
       endif
!-----------------------------------------------------------------------
!--    constraints on P' and FF'                                      --
!-----------------------------------------------------------------------
       if (kcalpa.gt.0) then
         do j=1,kcalpa
           nj=nj+1
           do i=1,kppcur
             alipc(nj,i)=calpa(i,j)*fwtxxx
           enddo
           do i=kppcur+1,kpcurn
             alipc(nj,i)=0.
           enddo
           xrsp(nj)=xalpa(j)*darea*fwtxxx
         enddo
       endif
!
       if (kcgama.gt.0) then
         do j=1,kcgama
           nj=nj+1
           do i=1,kppcur
             alipc(nj,i)=0.0
           enddo
           do i=kppcur+1,kpcurn
             alipc(nj,i)=cgama(i-kppcur,j)*fwtxxx
           enddo
           xrsp(nj)=xgama(j)*darea*fwtxxx
         enddo
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
       if (nj.gt.0) then
         call sdecm(alipc,npcur2,nj,kwcurn,xrsp,npcur2,nnn,wlipc,work,ier)
         if (ier.eq.129) then
           kerror=1
           call errctrl_msg('currnt', &
                            'sdecm failed to converge (location 2)')
           return
         endif
         cond=ier
         toler=1.0e-06_dp*wlipc(1)
         do i=1,kwcurn
           t=0.0
           if(wlipc(i).gt.toler) t=xrsp(i)/wlipc(i)
           work(i)=t
         enddo
         do i=1,kwcurn
           brsp(nfsum+i)=0.0
           do j=1,kwcurn
             brsp(nfsum+i)=brsp(nfsum+i)+alipc(i,j)*work(j)
           enddo
         enddo
!-----------------------------------------------------------------------
!--      update total plasma current if needed to                     --
!-----------------------------------------------------------------------
         if (abs(fwtcur).le.1.e-30_dp.and.nqwant.gt.0) then
           cm=0.0
           do n=nfsum+1,nfnwcr
            cm=cm+brsp(n)*fgowpc(n-nfsum)
           enddo
           ipmhd(jtime)=cm
           ipmeas(jtime)=cm
         endif
!
       endif
       endif init_curr
       tcurrt=0.0
       tcurrp=0.0
       do i=1,nw
         do j=1,nh
           kk=(i-1)*nh+j
           pcurrt(kk)=0.0
           pcurrw(kk)=0.0
           SOL_curr: if (icutfp.eq.0) then
!----------------------------------------------------------------------
!--          no attached current                                     --
!----------------------------------------------------------------------
             if ((xpsi(kk).lt.0.0).or.(xpsi(kk).gt.1.0)) cycle
             pp0=ppcurr(xpsi(kk),kppcur)
             pcurrt(kk)=fpcurr(xpsi(kk),kffcur)/rgrid(i)
             if (kvtor.eq.0) then
               pcurrt(kk)=pcurrt(kk)+pp0*rgrid(i)
             else
!-----------------------------------------------------------------------
!--            rotation                                               --
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
           else SOL_curr
!-----------------------------------------------------------------------
!--          attached current                                         --
!-----------------------------------------------------------------------
             if ((xpsi(kk).ge.0.0).and.(xpsi(kk).le.1.0)) then
                pp0=ppcurr(xpsi(kk),kppcur)
                if (kvtor.eq.0) then
                  pcurrt(kk)=rgrid(i)*pp0
                else
!-----------------------------------------------------------------------
!--             rotation                                              --
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
                  pcurrt(kk)=pp0*rgrid(i)+pcurrw(kk)
                endif
             endif
             upsi=xpsi(kk)*xpsimin
             if ((upsi.ge.0.0).and.(upsi.le.1.0)) &
               pcurrt(kk)=pcurrt(kk)+fpcurr(xpsi(kk),kffcur)/rgrid(i)
             pcurrt(kk)=pcurrt(kk)*zero(kk)
           endif SOL_curr
           tcurrt=tcurrt+pcurrt(kk)
         enddo
       enddo
       if ((nitett.le.1).and.(icinit.lt.0)) then
         cratio=1.0
         return
       endif
!-------------------------------------------------------------------
!--    adjust to match total current                              --
!-------------------------------------------------------------------
       cratio=ipmhd(jtime)/tcurrt
       do kk=1,nwnh
         pcurrt(kk)=pcurrt(kk)*cratio
         pcurrw(kk)=pcurrw(kk)*cratio
         if ((xpsi(kk).ge.0.0).and.(xpsi(kk).le.1.0)) then
            tcurrp=tcurrp+pcurrt(kk)
         endif
       enddo
       do i=nfsum+1,nfsum+kwcurn
         brsp(i)=cratio*brsp(i)
       enddo
!
      end select
      return
      end subroutine currnt

!**********************************************************************
!>
!!    fpcurr computes the radial derivative
!!    of the poloidal current ff. ffcurr computes
!!    the poloidal current F=twopi RBt/mu0
!!
!!    @param upsi : flux location to produce result
!!    @param nnn : number of coefficients used in representation
!!
!**********************************************************************
      real*8 function fpcurr(upsi,nnn)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit none
      real*8 fpecrr,ffcurr,seval

      integer*4, intent(in) :: nnn
      real*8, intent(in) :: upsi
      integer*4 i,iiij,iijj,nn
      real*8 f0back,f0edge,fb22,ff22,siedge,xsidif,ypsi
      real*8 xpsii(nffcur)

      if (icutfp.gt.0) then
        ypsi=upsi*xpsimin
      else
        ypsi=upsi
      endif
      if (abs(ypsi).gt.1.0) then
        fpcurr=0.0
        return
      endif
      if (npsi_ext > 0) then
        fpcurr = seval(npsi_ext,ypsi,psin_ext,ffprim_ext,bfp_ext,cfp_ext,dfp_ext)
        fpcurr = fpcurr * cratiof_ext
        return
      endif
      fpcurr=0.0
      call setfp(ypsi,xpsii)
      do iiij=nbase+1,nbase+nnn
        iijj=iiij-nbase
        fpcurr=fpcurr+brsp(iiij)*xpsii(iijj)
      enddo
!----------------------------------------------------------------------
!--   edge hyperbolic tangent component                              --
!----------------------------------------------------------------------
      if (kedgef.eq.0) return
      siedge=(ypsi-fe_psin)/fe_width
      f0edge=f2edge
      f0back=f0edge/fe_width/sidif
      fpcurr=fpcurr+f0back/cosh(siedge)**2
      return
!
      entry fpecrr(upsi,nnn)
      if (icutfp.gt.0) then
        ypsi=upsi*xpsimin
      else
        ypsi=upsi
      endif
      if (abs(ypsi).gt.1.0) then
        fpecrr=0.0
        return
      endif
      fpecrr=0.0
      call setfp(ypsi,xpsii)
      do iiij=nbase+nnn,nbase+nnn
        iijj=iiij-nbase
        fpecrr=fpecrr+brsp(iiij)*xpsii(iijj)
      enddo
      return
!
      entry ffcurr(upsi,nnn)
      if (icutfp.gt.0) then
        ypsi=upsi*xpsimin
        xsidif=-sidif/xpsimin
      else
        ypsi=upsi
        xsidif=-sidif
      endif
      if (abs(ypsi).ge.1.0) then
        ffcurr=fbrdy
        return
      endif
      ffcurr=0.0
      call setff(ypsi,xpsii)
      do i=nbase+1,nbase+nnn
        nn=i-nbase
        ffcurr=ffcurr+brsp(i)*xpsii(nn)
      enddo
      fb22=fbrdy**2
      ff22=fb22+xsidif*ffcurr/constf2
      if(ff22.lt.0.0) ff22=fb22
      ffcurr=sqrt(ff22)*fbrdy/abs(fbrdy)
!----------------------------------------------------------------------
!--   edge hyperbolic tangent component                              --
!----------------------------------------------------------------------
      if(kedgef.eq.0) return
      siedge=(ypsi-fe_psin)/fe_width
      f0edge=f2edge/constf2
      ff22=ff22+f0edge*(tfedge-tanh(siedge))
      if(ffcurr.lt.0.0) ffcurr=fb22
      ffcurr=sqrt(ff22)*fbrdy/abs(fbrdy)
      return
      end function fpcurr
