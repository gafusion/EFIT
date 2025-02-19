#include "config.f"
!**********************************************************************
!>                                                                  
!!     This subroutine computes betas, li, and other physics quantities
!!       that are written to output files (called by shapesurf)
!!
!!     @param jtime : time index
!!
!!     @param rgrid : radius of grid points
!!
!!     @param zgrid : Z of grid points
!!
!!     @param kerror : error flag
!**********************************************************************
      subroutine betaliplus(jtime,rgrid,zgrid,kerror)
      use commonblocks,only: c,wk,copy,bkx,bky,cw,wkw,copyw,bwx, &
                  bwy,sifprw,bwprw,cwprw,dwprw,sfprw,sprwp
      use var_cww
      include 'eparm.inc'
      include 'modules1.inc'
      implicit none

      real*8, external :: fpcurr,ppcurr,prcurr,pwcurr,pwpcur,prcur4,pwcur4,seval
      integer*4, intent(in) :: jtime
      real*8, intent(in) :: rgrid(nw),zgrid(nh)
      integer*4, intent(out) :: kerror
      integer*4 i,ii,ip1,iiqmin,izzmax,j,k,kk,km1,kz,mkk,mx,my,iautoc, &
                ier,nfind
      real*8 aasi,abpol,axout,ayout,bbsi,bp2flx,bpolsq,bpolzs,circum, &
             const,dlbpi,dli,dmui,dsi,dx,dyww,eesi,enr,enz,erhor,erhoz, &
             f222,ff222,fnow,fpnow, &
             pp0,ppnow,ppw,pres0,prew0,prewx,ptop0,pwop0,pwp0r2, &
             rcccc,rcurrm,rdiml,rdimw,rgmvt,rho,rlnr, &
             sbp,sbpi,sdlbp,sdlbpol,sdnrho,sicut,siii,sisi,siwant, &
             sixxx,ssdlbp,sumbp,sumbp2,sumbzz,sumpr2,sumpre,sumprt, &
             sumprw,sumr2,sumvar,sumy2,sumz,xcmax,xcmin,xww,xycmax, &
             xycmin,xym,xyp,ycmax,ycmin,yoxm,yoxp,ypsz,yww,yxcmax, &
             yxcmin,zaaa,zbbb,zccc
      real*8 pds(6),xsier(nercur),worksi(nw),workrm(nw),bwork(nw), &
             cwork(nw),dwork(nw),x(nw),y(nh),dpleng(npoint)
      integer*4, parameter :: licalc=1,idovol=1,inorm=3,ibtcal=2 ! hardcoded options

      kerror = 0

      select case (licalc)
      case (1) ! standard
        call sets2d(psi,c,rgrid,nw,bkx,lkx,zgrid,nh,bky,lky,wk,ier)
        sumbp2=0.0
        do i=1,nw
          do j=1,nh
            kk=(i-1)*nh+j
            call seva2d(bkx,lkx,bky,lky,c,rgrid(i),zgrid(j),pds,ier,n333)
            bpolsq=(pds(2)**2+pds(3)**2)/rgrid(i)
            sumbp2=sumbp2+bpolsq*www(kk)
          enddo
        enddo
        sumbp2=sumbp2*twopi*darea
      case (2) ! rarely used
        sumbp2=sum((psi(1:nwnh)-psibry)*pcurrt(1:nwnh)*www(1:nwnh))*twopi**2*tmu
      end select
!
      psurfa(jtime)=0.0
      plengt(1)=0.0
      do i=1,nfound-1
        ip1=i+1
        dli=sqrt((xout(ip1)-xout(i))**2+(yout(ip1)-yout(i))**2)
        plengt(ip1)=plengt(i)+dli
        dpleng(i)=dli
        psurfa(jtime)=psurfa(jtime)+dli*0.5_dp*(xout(i)+xout(ip1))
      enddo
      psurfa(jtime)=psurfa(jtime)*twopi
!
      s1(jtime)=0.0
      s2(jtime)=0.0
      s3(jtime)=0.0
      qout(jtime)=0.0
      sbp=0.0
      sbpi=0.0
      sumbzz=0.0
      sumr2=0.0
      sumz=0.0
      do i=1,nfound-1
        ip1=i+1
        axout=(xout(i)+xout(ip1))/2.
        ayout=(yout(i)+yout(ip1))/2.
        erhor=axout-rcentr
        erhoz=ayout
        rho=sqrt(erhor**2+erhoz**2)
        enz=(xout(ip1)-xout(i))/2.
        enr=-(yout(ip1)-yout(i))/2.
        erhor=erhor/rho
        erhoz=erhoz/rho
        rlnr=sqrt(enr**2+enz**2)
        enr=enr/rlnr
        enz=enz/rlnr
        sdnrho=enr*erhor+enz*erhoz
        abpol=(bpol(ip1)+bpol(i))/2.
        sumbp=dpleng(i)*abpol
        sumbzz=sumbzz+(bpolz(ip1)+bpolz(i))/2.*(yout(ip1)-yout(i))
        select case (kbetapr)
        case (0)
          sumvar=dpleng(i)*abpol**2*axout
        case (1)
          sumvar=dpleng(i)*(abpol**2+4.0*pi*tmu*prbdry)*axout
        case (2)
          sumvar=dpleng(i)*(abpol**2+4.0*pi*tmu*pressb)*axout
        end select
        s1(jtime)=s1(jtime)+sumvar*sdnrho*rho
        s2(jtime)=s2(jtime)+sumvar*enr
        s3(jtime)=s3(jtime)+sumvar*ayout*enz
        sbp=sbp+sumbp
        sumr2=sumr2+sumbp*axout**2
        sumz=sumz+sumbp*ayout
        dlbpi=dpleng(i)/abpol
        sbpi=sbpi+dlbpi
        qout(jtime)=qout(jtime)+dlbpi/axout**2
      enddo
      select case (inorm)
      case (2) 
        bp2flx=2.*rcentr/volume(jtime)*(pi*tmu*ipmhd(jtime))**2 &
                 *1.0E+06_dp
      case (3)
        bp2flx=(tmu*2.0*pi*ipmhd(jtime)/plengt(nfound))**2
      case (4)
        bp2flx=2.*(tmu*ipmhd(jtime)/aminor(jtime))**2 &
                 /(elong(jtime)**2+1.)*1.0e+04_dp
      case default
        bp2flx=sbp/sbpi
      end select
      bpolav(jtime)=sqrt(bp2flx)
      rcurrt(jtime)=sqrt(sumr2/twopi/tmu/abs(ipmhd(jtime)))*100.
      zcurrt(jtime)=sumz/twopi/tmu/abs(ipmhd(jtime))*100.
      const=twopi/volume(jtime)*1.0e+06_dp/bp2flx
      s1(jtime)=const*s1(jtime)
      s2(jtime)=const*rcentr*s2(jtime)
      s3(jtime)=const*s3(jtime)
      sumbzz=abs(sumbzz)/tmu/abs(ipmhd(jtime))/twopi
      rcurrm=rcentr
      li(jtime)=sumbp2/volume(jtime)/bp2flx*1.0e+06_dp
      li3(jtime)=(tmu*2.0*pi*ipmhd(jtime))**2*rcntr(jtime)/200.0
      li3(jtime)=sumbp2/li3(jtime)
      betap(jtime)=s1(jtime)/4.+s2(jtime)/4.*(1.+rcurrm/rcentr) &
                   -li(jtime)/2.
      betat(jtime)=betap(jtime)*bp2flx/bcentr(jtime)**2*100.
      betat(jtime)=betat(jtime)*(rcntr(jtime)/100./rcentr)**2
      qout(jtime)=qout(jtime)*abs(bcentr(jtime))*rcentr/twopi
      wmhd(jtime)=1.5_dp*betap(jtime)*bp2flx/2./tmu/2./pi*volume(jtime) &
                    /1.0e+06_dp
!---------------------------------------------------------------------
!--   calculations of current moment y2                             --
!---------------------------------------------------------------------
      sumy2=0.0
      rcccc=rcurrt(jtime)/100.
      rmajz0(nw)=0.0
      do i=1,nfound-1
        ip1=i+1
        axout=(xout(i)+xout(ip1))/2./rcccc
        ayout=(yout(i)+yout(ip1))/2./rcccc
        abpol=(bpol(ip1)+bpol(i))/2.
        sumbp=dpleng(i)*abpol
        axout=axout**2
        ayout=ayout**2
        f222=axout*(axout-4.*ayout)
        ff222=f222-2.*axout+1.
        sumy2=sumy2+sumbp*ff222
        if(yout(i)*yout(ip1).le.0.0.and.xout(i).gt.rcentr) &
          rmajz0(nw)=xout(i)-yout(i)*(xout(ip1)-xout(i)) &
                             /(yout(ip1)-yout(i))
      enddo
      yyy2(jtime)=sumy2/(tmu*2.0*pi*ipmhd(jtime))/4. &
                  *(rcurrt(jtime)/aminor(jtime))**2
!
      dsi=(psibry-simag)/(nw-1)
      mx=(rmaxis-rgrid(1))/drgrid+1
      my=(zmaxis-zgrid(1))/dzgrid+1
      mkk=(mx-1)*nh+my+1
      sicut=psi(mkk)
      volp(nw)=volume(jtime)/1.0e+06_dp
      volp(1)=0.0
      if (abs(zmaxis).le.0.001_dp) then
        rmajz0(1)=rmaxis
      else
        rmajz0(1)=0.0
      endif
      r1surf(1)=1./rmaxis
      r2surf(1)=1./rmaxis**2
      r1surf(nw)=r1bdry
      r2surf(nw)=r2bdry
      bpolss(1)=0.0
      bpolss(nw)=bpolav(jtime)
!----------------------------------------------------------------------
!--   rotational terms                                               --
!----------------------------------------------------------------------
      rotation: if (kvtor.gt.0) then
        sumprt=0.0
        sumprw=0.0
        prew0=pwcur4(1,0.5_dp)
        pres0=prcur4(1,0.5_dp)
        pwprim=sprwp
        pressw=sfprw
        do i=1,nw
          do j=1,nh
            kk=(i-1)*nh+j
            pres0=prcur4(0,xpsi(kk))-prbdry
            prew0=pwcur4(0,xpsi(kk))-preswb
            select case (kvtor)
            case (1)
              presst(kk)=pres0+prew0*rgrvt(i)
            case (2)
              if (abs(pres0).gt.1.e-10_dp) then
                pwop0=prew0/pres0
                pwp0r2=pwop0*rgrvt(i)
              else
                pwop0=0.0
                pwp0r2=0.0
              endif
              presst(kk)=pres0*(1.+0.5_dp*pwp0r2**2)+prew0*rgrvt(i)
            case (3,11)
              if (abs(pres0).gt.1.e-10_dp) then
                pwop0=prew0/pres0
                ptop0=exp(pwop0*rgrvt(i))
              else
                ptop0=1.0
              endif
              presst(kk)=pres0*ptop0
            end select
            sumprt=sumprt+rgrid(i)*presst(kk)*www(kk)
            prewx=prew0*(rgrid(i)/rvtor)**2
            sumprw=sumprw+rgrid(i)*prewx*www(kk)
          enddo
        enddo
        sumprt=sumprt*darea*twopi
        sumprw=sumprw*darea*twopi
        call sets2d(presst,cw,rgrid,nw,bwx,lwx,zgrid,nh,bwy,lwy,wkw,ier)
      endif rotation
!
      sumpre=0.0
      rzzmax(1)=rmaxis
      zzmax(1)=zmaxis
      rhovn(1)=0.0
      do i=2,nw-1
        ii=nw-i+1
        siii=(i-1)*dsi
        siwant=psibry-siii
        siii=1.0_dp-1.0_dp/(nw-1)*(i-1)
        single_volume: if (idovol.le.1) then
        rzzmax(ii)=-99.0
        call surfac(siwant,psi,nw,nh,rgrid,zgrid,bpol,bpolz,nfind, &
                    npoint,drgrid,dzgrid,xmin,xmax,ymin,ymax,1, &
                    rmaxis,zmaxis,negcur,kerror,2)
        if (nfind.le.40.and.icntour.eq.0) then
#ifdef DEBUG_LEVEL2
          write (6,*) ' SHAPE/BETALI kerror,i,nfind = ',kerror,i,nfind
#endif
          call cntour(rmaxis,zmaxis,siwant,xcmin,xcmax,ycmin,ycmax, &
                      yxcmin,yxcmax,xycmin,xycmax,30.,drgrid,0.03_dp, &
                      0.01_dp,0.01_dp,xmin,xmax,ymin,ymax,0,iautoc, &
                      bpol,bpolz,nfind,rgrid,nw,zgrid,nh, &
                      c,n222,nh2,nttyo,npoint, &
                      negcur,bkx,lkx,bky,lky,kerror)
#ifdef DEBUG_LEVEL2
          write (6,*) ' BETALI/CNTOUR kerror,nfind = ',kerror,nfind
#endif
          if(kerror /= 0) return
        endif
        found: if (nfind.ge.10) then
        r1surf(ii)=0.0
        r2surf(ii)=0.0
        volp(ii)=0.0
        rhovn(ii)=0.0
        xym=bpol(1)*bpolz(1)
!------------------------------------------------------------------
!--     integration over z from 0 to bpolz                       --
!------------------------------------------------------------------
        xww=bpol(1)
        dyww=bpolz(1)/(nh-1)
        bpolzs=0.5*fpol(ii)
        do kz=2,nh-1
          yww=dyww*(kz-1)
          call seva2d(bkx,lkx,bky,lky,c,xww,yww,pds,ier,n111)
          ypsz=(simag-pds(1))/sidif
          fnow=seval(nw,ypsz,sigrid,fpol,bbfpol,ccfpol,ddfpol)
          bpolzs=bpolzs+fnow
        enddo
        yww=0.0
        call seva2d(bkx,lkx,bky,lky,c,xww,yww,pds,ier,n111)
        ypsz=(simag-pds(1))/sidif
        fnow=seval(nw,ypsz,sigrid,fpol,bbfpol,ccfpol,ddfpol)
        bpolzs=bpolzs+fnow*0.5
        yoxm=bpolzs/bpol(1)*dyww
!
        izzmax=2
        zzmax(ii)=bpolz(1)
        rmajz0(ii)=0.0
        ssdlbp=0.0
        sdlbpol=0.0
        circum=0.0
        do k=2,nfind
          km1=k-1
          xyp=bpol(k)*bpolz(k)
!------------------------------------------------------------------
!--       integration over z from 0 to bpolz                     --
!------------------------------------------------------------------
          xww=bpol(k)
          dyww=bpolz(k)/(nh-1)
          bpolzs=0.5*fpol(ii)
          do kz=2,nh-1
            yww=dyww*(kz-1)
            call seva2d(bkx,lkx,bky,lky,c,xww,yww,pds,ier,n111)
            ypsz=(simag-pds(1))/sidif
            fnow=seval(nw,ypsz,sigrid,fpol,bbfpol,ccfpol,ddfpol)
            bpolzs=bpolzs+fnow
          enddo
          yww=0.0
          call seva2d(bkx,lkx,bky,lky,c,xww,yww,pds,ier,n111)
          ypsz=(simag-pds(1))/sidif
          fnow=seval(nw,ypsz,sigrid,fpol,bbfpol,ccfpol,ddfpol)
          bpolzs=bpolzs+fnow*0.5
          yoxp=bpolzs/bpol(k)*dyww
!
          dx=bpol(k)-bpol(km1)
          volp(ii)=volp(ii)+(xyp+xym)/2.0*dx
          rhovn(ii)=rhovn(ii)+(yoxp+yoxm)/2.0*dx
          xym=xyp
          yoxm=yoxp
          dli=sqrt((bpol(k)-bpol(km1))**2+(bpolz(k)-bpolz(km1))**2)
          xww=(bpol(k)+bpol(km1))/2.
          yww=(bpolz(k)+bpolz(km1))/2.
          call seva2d(bkx,lkx,bky,lky,c,xww,yww,pds,ier,n333)
          bpolsq=sqrt(pds(2)**2+pds(3)**2)/xww
          sdlbp=dli/bpolsq
          sdlbpol=sdlbpol+dli*bpolsq
          circum=circum+dli
          r1surf(ii)=r1surf(ii)+sdlbp/xww
          r2surf(ii)=r2surf(ii)+sdlbp/xww**2
          rr2bpsurf(ii)=r2surf(ii)
          ssdlbp=ssdlbp+sdlbp
          if (bpolz(k).gt.zzmax(ii)) then
            izzmax=k
            zzmax(ii)=bpolz(k)
          endif
          if(bpolz(km1)*bpolz(k).le.0.0.and.bpol(km1).gt.rcentr) &
            rmajz0(ii)=bpol(km1)-bpolz(km1)*(bpol(k)-bpol(km1)) &
                                           /(bpolz(k)-bpolz(km1))
        enddo
        r1surf(ii)=r1surf(ii)/ssdlbp
        r2surf(ii)=r2surf(ii)/ssdlbp
        bpolss(ii)=sdlbpol/circum
        call qfit(n333,bpol(izzmax-1),bpol(izzmax),bpol(izzmax+1), &
                  bpolz(izzmax-1),bpolz(izzmax),bpolz(izzmax+1), &
                  zaaa,zbbb,zccc,ierr)
        if (ierr.ne.0) then
          kerror = 1
          return
        endif
        rzzmax(ii)=-zbbb/2./zaaa
        zzmax(ii)=zaaa*rzzmax(ii)**2+zbbb*rzzmax(ii)+zccc
        volp(ii)=abs(volp(ii))*twopi
        rhovn(ii)=rhovn(ii)/rhovn(nw)
        else found
        bpolss(ii)=0.0
        rmajz0(ii)=0.0
        sisi=simag-siwant
        bbsi=sqrt(sisi/siaz)
        aasi=sqrt(sisi/siar)
        volp(ii)=twopi*rmaxis*pi*aasi*bbsi
        rhovn(ii)=pi*aasi*bbsi/rmaxis/rhovn(nw)*fpol(1)
!
        eesi=bbsi/aasi
        select case (icurrt)
        case (1)
          rdiml=rmaxis/srma
          cjmaxi=cratio*(sbeta*rdiml+2.*salpha/rdiml)/darea
          if(kvtor.gt.0) &
            cjmaxi=cjmaxi+cratio/darea*sbetaw*rdiml*(rdiml**2-1.)
        case (2,5)
          cjmaxi=(rmaxis*ppcurr(0.0,kppcur) &
                        +fpcurr(0.0,kffcur)/rmaxis)*cratio/darea
          if(kvtor.gt.0) &
            cjmaxi=sum(rjjjx(1:kppcur)*brsp(nfsum+1:nfsum+kppcur) &
                      +rjjfx(1:kffcur)*brsp(nfsum+kppcur+1:nfsum+kppcur+kffcur) &
                      +rjjwx(1:kwwcur)*brsp(nfsum+kpcurn+1:nfsum+kpcurn+kwwcur))/darea
        case (4)
          rdiml=rmaxis/rzero
          cjmaxi=cratio/darea*(rdiml+rbetap/rdiml)
          cjmaxi=cratio/darea*(rdiml+rbetap/rdiml)
          if (kvtor.eq.1) then
            rgmvt=(rmaxis/rvtor)**2-1.
            cjmaxi=cjmaxi+cratio/darea*rdiml*rbetaw*rgmvt
          elseif (kvtor.eq.11) then
            pres0=prcur4(1,0.0)
            prew0=pwcur4(1,0.0)
            rgmvt=(rmaxis/rvtor)**2-1.
            pwop0=prew0/pres0
            ptop0=exp(pwop0*rgmvt)
            pp0= 1.-pwop0*rgmvt
            ppw=rbetaw*rgmvt
            cjmaxi=cjmaxi+(pp0+ppw)*rdiml*ptop0
          endif
        end select
        rr2bpsurf(ii)=(eesi**2+1.)/rmaxis**2/tmu/cjmaxi
        r1surf(ii)=1./rmaxis
        endif found
        endif single_volume
        sumpre=sumpre+volp(ii)*pprime(ii)
      enddo
      sumpre=sumpre+volp(nw)*pprime(nw)/2.
      if(ibtcal.le.1) return
      select case (kbetapr)
      case (0)
        sumpre=-sumpre*dsi/volp(nw)
      case (1)
        sumpre=-sumpre*dsi/volp(nw)+prbdry
      case (2)
        sumpre=-sumpre*dsi/volp(nw)+pressb
      end select
      betap(jtime)=sumpre*2.0*twopi*tmu/bp2flx
      betat(jtime)=sumpre*2.0*twopi*tmu/bcentr(jtime)**2
      betat(jtime)=100.*betat(jtime)*(rcntr(jtime)/100./rcentr)**2
      wmhd(jtime)=1.5_dp*betap(jtime)*bp2flx/2./tmu/2./pi*volume(jtime) &
                    /1.0e+06_dp
      pasman=ipmhd(jtime)/1.e4_dp/aminor(jtime)/abs(bcentr(jtime))
      pasman=pasman*rcntr(jtime)/100./rcentr
      betatn=betat(jtime)/pasman
      dmui=1.0e+06_dp*diamag(jtime)*4.*pi*bcentr(jtime)*rcentr &
           /bp2flx/volume(jtime)
      betapd(jtime)=s1(jtime)/2.+s2(jtime)/2.*(1.-rcurrm/rcentr)-dmui
      betatd(jtime)=betapd(jtime)*bp2flx/bcentr(jtime)**2*100.
      betatd(jtime)=betatd(jtime)*(rcntr(jtime)/100./rcentr)**2
      wdia(jtime)=1.5_dp*betapd(jtime)*bp2flx/2./tmu/2./pi*volume(jtime) &
                     /1.0e+06_dp
!-----------------------------------------------------------------------
!--   rotational terms                                                --
!-----------------------------------------------------------------------
      if (kvtor.gt.0) then
        sumprt=sumprt/volp(nw)
        sumprw=sumprw/volp(nw)
        betap(jtime)=sumprt/sumpre*betap(jtime)
        betat(jtime)=sumprt/sumpre*betat(jtime)
        wmhd(jtime)=sumprt/sumpre*wmhd(jtime)
        betapw(jtime)=sumprw/sumprt*betap(jtime)
        betatw(jtime)=sumprw/sumprt*betat(jtime)
        wplasw(jtime)=betapw(jtime)*bp2flx/2./tmu/2./pi &
                      *volume(jtime)/1.0e+06_dp
      endif
!----------------------------------------------------------------------
!--   get normalized radial coordinate square root of toroidal flux  --
!--     at uniform poloidal flux grid sigrid                         --
!----------------------------------------------------------------------
      rhovn(nw)=1.0
      rhovn(2:(nw-1))=sqrt(abs(rhovn(2:(nw-1))))
      call zpline(nw,sigrid,rhovn,brhovn,crhovn,drhovn)
      if (kstark.gt.0.or.kdomse.gt.0) then
        do i=1,nstark
          sixxx=sigam(i)
          if (sixxx.le.1.0) then
            rhogam(i)=seval(nw,sixxx,sigrid,rhovn,brhovn,crhovn,drhovn)
          else
            rhogam(i)=1.1_dp
          endif
        enddo
      endif
!----------------------------------------------------------------------
!--   get electrostatic potential                                    --
!----------------------------------------------------------------------
      if (keecur.gt.0) then
        do i=1,nw
          call seter(sigrid(i),xsier)
          epoten(i)=sum(cerer(1:keecur)*xsier(1:keecur))
        enddo
      endif
      if(idovol.gt.1) return
!-----------------------------------------------------------------------
!--   compute beta*, taking P(1)=0.0    10/25/90                      --
!-----------------------------------------------------------------------
      sumpr2=-2.*sum(volp(2:(nw-1))*pprime(2:(nw-1))*pres(2:(nw-1))) &
                *dsi/volp(nw)
      if (sumpr2.ge.0.0) then
        sumpr2=sqrt(sumpr2)
      else
        sumpr2=0.0
      endif
      betat2=sumpr2*2.0*twopi*tmu/bcentr(jtime)**2
      betat2=100.*betat2*(rcntr(jtime)/100./rcentr)**2
!-----------------------------------------------------------------------
!--   compute the safety factor profile
!-----------------------------------------------------------------------
      qpsi(1:(nw-1))=abs(fpol(1:(nw-1)))/twopi*rr2bpsurf(1:(nw-1))
      qpsi(nw)=qout(jtime)
      qpsi(1)=qmaxis
      qmin=qpsi(1)
      iiqmin=1
      do i=2,nw
        if (qpsi(i).lt.qmin) then
          qmin=qpsi(i)
          iiqmin=i
        endif
      enddo
      rhoqmin=sqrt(volp(iiqmin)/volp(nw))
!
      btaxp(jtime)=fpol(1)/rmaxis
      btaxv(jtime)=fpol(nw)/rmaxis
      vbeta0=pres(1)/sumpre*betat(jtime)
!-----------------------------------------------------------------------
!--   get R at evenly spaced psi along outboard midplane
!-----------------------------------------------------------------------
      do i=1,nw
        workrm(i)=rmaxis+(xmax-rmaxis)/(nw-1)*(i-1)
        call seva2d(bkx,lkx,bky,lky,c,workrm(i),zmaxis,pds,ier,n111)
        worksi(i)=(pds(1)-simag)/(psibry-simag)
      enddo
      call zpline(nw,worksi,workrm,bwork,cwork,dwork)
      do i=1,nw
        sixxx=1.0_dp/(nw-1)*(i-1)
        rpres(i)=seval(nw,sixxx,worksi,workrm,bwork,cwork,dwork)
      enddo
!-----------------------------------------------------------------------
!--   compute <jt/R>/<1/R> profile (ignoring rotation)
!-----------------------------------------------------------------------
      do i=1,nw
        ppnow=pprime(i)
        fpnow=ffprim(i)/twopi/tmu
        cjor(i)=(ppnow+fpnow*r2surf(i))/r1surf(i)
        fpnow=ffprec(i)/twopi/tmu
        cjorec(i)=fpnow*r2surf(i)/r1surf(i)
      enddo
!
      return
 1980 format (1x,i6)
 2000 format (1x,6e12.5)
      end subroutine betaliplus

!**********************************************************************
!!     This subroutine computes only betas and li (part of the fitting
!!       loop)
!!                                                                  
!!     @param jtime : time index
!!
!!     @param rgrid : radius of grid points
!!
!!     @param zgrid : Z of grid points
!!
!!     @param kerror : error flag
!**********************************************************************
      subroutine betali(jtime,rgrid,zgrid,kerror)
      use commonblocks,only: c,wk,copy,bkx,bky
      include 'eparm.inc'
      include 'modules1.inc'
      implicit none

      real*8, external :: ppcurr
      integer*4, intent(in) :: jtime
      real*8, intent(in) :: rgrid(nw),zgrid(nh)
      integer*4, intent(out) :: kerror
      integer*4 i,ii,ip1,j,k,kk,km1,iautoc,ier,nfind
      real*8 aasi,bbsi,bp2flx,bpolsq,const,dli,dsi,dx,siii,sisi,siwant, &
             sumbp2,sumpre,xcmax,xcmin,xycmax,xycmin,xym,xyp,ycmax,ycmin, &
             yxcmax,yxcmin
      real*8 pds(6),x(nw),y(nh),dpleng(npoint),xxs(npoint),yys(npoint)

      kerror = 0
!
      if(ivacum.eq.1) return
      sumbp2=0.0
      call sets2d(psi,c,rgrid,nw,bkx,lkx,zgrid,nh,bky,lky,wk,ier)
      do i=1,nw
        do j=1,nh
          kk=(i-1)*nh+j
          call seva2d(bkx,lkx,bky,lky,c,rgrid(i),zgrid(j),pds,ier,n333)
          bpolsq=(pds(2)**2+pds(3)**2)/rgrid(i)
          sumbp2=sumbp2+bpolsq*www(kk)
        enddo
      enddo
      sumbp2=sumbp2*twopi*darea
!
      plengt(1)=0.0
      do i=1,nfound-1
        ip1=i+1
        dli=sqrt((xout(ip1)-xout(i))**2+(yout(ip1)-yout(i))**2)
        plengt(ip1)=plengt(i)+dli
        dpleng(i)=dli
      enddo
!
      bp2flx=(tmu*2.0*pi*ipmhd(jtime) &
              /plengt(nfound))**2
      const=twopi/volume(jtime)*1.0e+06_dp/bp2flx
      li(jtime)=sumbp2/volume(jtime)/bp2flx*1.0e+06_dp
      li3(jtime)=(tmu*2.0*pi*ipmhd(jtime))**2*rcntr(jtime)/200.0
      li3(jtime)=sumbp2/li3(jtime)
!
      dsi=(psibry-simag)/(nw-1)
      volp(nw)=volume(jtime)/1.0e+06_dp
      volp(1)=0.0
      select case (icurrt)
      case (1)
        pprime(1)=cratio*sbeta/darea/srma
        pprime(nw)=pprime(1)
      case (2,5)
        pprime(nw)=ppcurr(x111,kppcur)/darea
        pprime(1)=ppcurr(x000,kppcur)/darea
      case (4)
        call currnt(n222,jtime,n222,kerror)
        if(kerror.gt.0) return
        pprime(1)=cratio/darea/rzero
        pprime(nw)=pprime(1)*gammap
      end select
      sumpre=0.0
      do i=2,nw-1
        ii=nw-i+1
        siii=(i-1)*dsi
        siwant=psibry-siii
        siii=1.0_dp-1.0_dp/(nw-1)*(i-1)
        call surfac(siwant,psi,nw,nh,rgrid,zgrid,xxs,yys,nfind, &
                    npoint,drgrid,dzgrid,xmin,xmax,ymin,ymax,1, &
                    rmaxis,zmaxis,negcur,kerror,2)
        if (nfind.le.40.and.icntour.eq.0) then
          call cntour(rmaxis,zmaxis,siwant,xcmin,xcmax,ycmin,ycmax, &
                      yxcmin,yxcmax,xycmin,xycmax,30.,drgrid,0.03_dp, &
                      0.01_dp,0.01_dp,xmin,xmax,ymin,ymax,0,iautoc, &
                      xxs,yys,nfind,rgrid,nw,zgrid,nh, &
                      c,n222,nh2,nttyo,npoint, &
                      negcur,bkx,lkx,bky,lky,kerror)
          if(kerror.gt.0) return
        endif
        if (nfind.ge.10) then
          volp(ii)=0.0
          xym=xxs(1)*yys(1)
          do k=2,nfind
            km1=k-1
            xyp=xxs(k)*yys(k)
            dx=xxs(k)-xxs(km1)
            volp(ii)=volp(ii)+(xyp+xym)/2.0*dx
            xym=xyp
          enddo
          volp(ii)=abs(volp(ii))*twopi
        else
          sisi=simag-siwant
          bbsi=sqrt(sisi/siaz)
          aasi=sqrt(sisi/siar)
          volp(ii)=twopi*rmaxis*pi*aasi*bbsi
        endif
!
        select case (icurrt)
        case (1)
          pprime(ii)=pprime(1)
        case (2,5)
          pprime(ii)=ppcurr(siii,kppcur)/darea
        case (4)
          pprime(ii)=(1.-siii**enp)**emp*(1.-gammap)+gammap
          pprime(ii)=pprime(1)*pprime(ii)
        end select
        sumpre=sumpre+volp(ii)*pprime(ii)
      enddo
      sumpre=sumpre+volp(nw)*pprime(nw)/2.
      select case (kbetapr)
      case (0)
        sumpre=-sumpre*dsi/volp(nw)
      case (1)
        sumpre=-sumpre*dsi/volp(nw)+prbdry
      case (2)
        sumpre=-sumpre*dsi/volp(nw)+pressb
      end select
      betap(jtime)=sumpre*2.0*twopi*tmu/bp2flx
      betat(jtime)=sumpre*2.0*twopi*tmu/bcentr(jtime)**2
      betat(jtime)=100.*betat(jtime)*(rcntr(jtime)/100./rcentr)**2
      pasman=ipmhd(jtime)/1.e4_dp/aminor(jtime)/abs(bcentr(jtime))
      pasman=pasman*rcntr(jtime)/100./rcentr
      betatn=betat(jtime)/pasman
!
      return
      end subroutine betali
