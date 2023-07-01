#include "config.f"
!**********************************************************************
!>
!!    computes the dimensionless poloidal fluxes for the r-z grid and
!!      updates additional parameters used in the constraints
!!
!!    @param ix : equilibrium (inner) loop iteration index
!!    @param ixt : total iteration index (current+equilibirum loops)
!!    @param ixout : current profile (outer) loop iteration index
!!    @param jtime : time index
!!    @param kerror : error flag
!!
!**********************************************************************
      subroutine update_params(ix,ixt,ixout,jtime,kerror)
      use commonblocks,only: c,bkx,bky,zeros,xouts,youts, &
           rsplt,zsplt,csplt
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      real*8 ffcurr,fpcurr,ppcurr,prcur4,pwcur4
      integer*4, intent(in) :: ix,ixt,ixout,jtime
      integer*4, intent(out) :: kerror
      integer*4 i,iii,j,k,kk,nn,ier,itot,nfounc,nnerr,nqend,nzz
      real*8 aouter,avebp,bpmin,bpnow,dels,delx,dely,det,dx,epssep, &
             fpnow,pp0,ppw,pres0,prew0,ptop0,pwop0, &
             r1qdry,r2qdry,r2surq,r2surs,rdiml,relsi,rgmvt, &
             sdlbp,sdlobp,sibpmin,sisinow,siwant,xerr, &
             xpsimins,xs,xvsmax,xvsmin,xww,xym,xyma,xyp,xypa,yerr,ys,yww
      real*8 pds(6)
      integer*4, parameter :: nnn=1,nqiter=10,isplit=8
      real*8, parameter :: psitol=1.0e-04_dp,cdum=1.0
      real*8 radbou,xguess,xltrac,yguess,zmaxis_last
      save radbou,xltrac,xguess,yguess
      data zmaxis_last/0.0/
!
      ! initialize variables
      zmaxis=0.0
      simag=0.0
      volume=0.0
      rout=0.0
      aminor=0.0
!
      if(ivacum.eq.1) return
      if (ixt.le.1) then
        xguess=(rgrid(1)+rgrid(nw))/2.
        yguess=(zgrid(1)+zgrid(nh))/2.
        if(zbound.ne.0.0) yguess=zbound
        if(rbound.ne.0.0) xguess=rbound
        xltrac=xlmin
        if(ibound.eq.-1) xltrac=xlmax ! hardcoded option only...
        radbou=(xguess+xltrac)/2.
      endif
!----------------------------------------------------------------------
!--   first set up bi-cubic spline interpolation in findax           --
!----------------------------------------------------------------------
#ifdef DEBUG_LEVEL2
      write (6,*) 'Entering findax'
#endif
      !print *, 'nw,nh,rgrid,zgrid',nw,nh,rgrid,zgrid
      !print *, 'rmaxis,zmaxis,simag', rmaxis,zmaxis,simag
      !print *, 'psibry,rseps(1,jtime),zseps(1,jtime)', psibry,rseps(1,jtime),zseps(1,jtime)
      !print *, 'xout,yout,nfound,psi', xout,yout,nfound,psi
      !print *, 'xmin,xmax,ymin,ymax',xmin,xmax,ymin,ymax
      !print *,  'zxmin,zxmax,rymin,rymax' , zxmin,zxmax,rymin,rymax
      !print *,  'dpsi,bpol,bpolz' , dpsi,bpol,bpolz
      !print *, 'limitr,xlim,ylim,limfag', limitr,xlim,ylim,limfag
      !print *, 'ixt,jtime,kerror', ixt,jtime,kerror
      call findax(nw,nh,rgrid,zgrid,rmaxis,zmaxis,simag, &
                  psibry,rseps(:,jtime),zseps(:,jtime),10, &
                  xout,yout,nfound,psi,xmin,xmax,ymin,ymax, &
                  zxmin,zxmax,rymin,rymax,dpsi,bpol,bpolz, &
                  limitr,xlim,ylim,limfag,ixt,jtime,kerror)

      if(kerror.gt.0) return
#ifdef DEBUG_LEVEL2
      if (nsol.gt.0) then
        ier=0 ! set only for consistent output? (not useful...) 
        write (6,*) 'UPDATE R,Z,Si,Err = ', rsol(1),zsol(1),wsisol,ier
        call seva2d(bkx,lkx,bky,lky,c,rbdry(1),zbdry(1),pds,ier,n111)
        write (6,*) 'UPDATE R,Z,Si,Err = ', rbdry(1),zbdry(1),pds(1),ier
        call seva2d(bkx,lkx,bky,lky,c,rbdry(nbdry),zbdry(nbdry),pds,ier &
                   ,n111)
        write (6,*) 'UPDATE R,Z,Si,Err = ', rbdry(nbdry),zbdry(nbdry) &
                    ,pds(1),ier

        write(6,*) 'rsplt(kk),zsplt(kk)',rbdry(nbdry),zbdry(nbdry)
        write(6,*) 'lkx, lky',lkx,lky
        write(6,*) 'pds,ier,n111', pds,ier,n111
        call seva2d(bkx,lkx,bky,lky,c,xout(1),yout(1),pds,ier,n111)
        write (6,*) 'UPDATE simag,psibry,n1,si = ',simag,psibry,n111,pds(1)
        write (6,*) 'UPDATE Z,R = ', zgrid(33),(rgrid(i),i=45,45)
        write (6,*) 'UPDATE si = ',(psi((i-1)*65+33),i=45,45)
        call seva2d(bkx,lkx,bky,lky,c,rgrid(45),zgrid(33),pds,ier,n111)
        write (6,*) 'UPDATE R,Z,si = ', rgrid(45),zgrid(33),pds(1)
      endif
#endif
!-----------------------------------------------------------------------
!--   Trace boundary, first check for counter beam injection          --
!-----------------------------------------------------------------------
      if (ipmeas(jtime).lt.-1.e3_dp) then
        nnerr=10000
      else
        nnerr=0
      endif
      call bound(psi,nw,nh,nwnh,psibry,xmin,xmax,ymin,ymax, &
                 zero,rgrid,zgrid,xguess,yguess,ixt,limitr,xlim,ylim, &
                 xout,yout,nfound,xltrac,npoint,rymin,rymax,dpsi, &
                 zxmin,zxmax,nnerr,ishot,itime, &
                 limfag,radbou,kbound,tolbndpsi)
      if (nnerr.gt.0) then
        kerror=1
        return
      endif
!----------------------------------------------------------------------
!--   find magnetic axis and poloidal flux at axis simag             --
!----------------------------------------------------------------------
#ifdef DEBUG_LEVEL2
      write (6,*) 'Entering findax at simag'
#endif
      call findax(nw,nh,rgrid,zgrid,rmaxis,zmaxis,simag, &
                  psibry,rseps(:,jtime),zseps(:,jtime),20, &
                  xout,yout,nfound,psi,xmin,xmax,ymin,ymax, &
                  zxmin,zxmax,rymin,rymax,dpsi,bpol,bpolz, &
                  limitr,xlim,ylim,limfag,ixt,jtime,kerror)
      if(kerror.gt.0) return
      sidif=simag-psibry
      eouter=(ymax-ymin)/(xmax-xmin)
      zplasm=(ymin+ymax)/2.
      aouter=(xmax-xmin)/2.
!-----------------------------------------------------------------------
!--   force free current in the scrape-off layer                      --
!-----------------------------------------------------------------------
      SOL_curr: if (icutfp.eq.2) then
        xvsmin=1.e10_dp
        xvsmax=-1.e10_dp
        if (limvs.eq.0) then
          itot=isplit*isplit
          do k=1,nvesel
            call splitc(isplit,rsplt,zsplt,csplt, &
                  rvs(k),zvs(k),wvs(k),hvs(k),avs(k),avs2(k),cdum)
            do kk=2,itot
              call seva2d(bkx,lkx,bky,lky,c,rsplt(kk),zsplt(kk), &
                          pds,ier,n111)
              write(6,*) 'rgrid(i),zgrid(j)',rsplt(kk),zsplt(kk)
              write(6,*) 'lkx, lky',lkx,lky
              write(6,*) 'pds,ier,n111', pds,ier,n111
              xvsmin=min(xvsmin,pds(1))
              xvsmax=max(xvsmax,pds(1))
            enddo
          enddo
        else
          do k=1,limitr-1
            delx=xlim(k+1)-xlim(k)
            dely=ylim(k+1)-ylim(k)
            dels=sqrt(delx**2+dely**2)
            nn=dels/0.002_dp
            nn=max(5,nn)
            delx=delx/(nn-1)
            dely=dely/(nn-1)
            do kk=2,nn
              xww=xlim(k)+delx *(kk-1)
              yww=ylim(k)+dely *(kk-1)
              call seva2d(bkx,lkx,bky,lky,c,xww,yww,pds,ier,n111)
              xvsmin=min(xvsmin,pds(1))
              xvsmax=max(xvsmax,pds(1))
            enddo
          enddo
        endif
!--------------------------------------------------------------
!--     exclude private region flux                          --
!--------------------------------------------------------------
        xvsmax=psibry
!--------------------------------------------------------------
!--     possible second separatrix                           --
!--------------------------------------------------------------
        rsepex=-999.
        yvs2=1000.
        skipvs: if (kskipvs.ne.0) then
        avebp=ipmhd(jtime)*tmu/aouter
        bpmin=avebp
        ! TODO: sibpmin has not yet been defined...
!        sibpold=sibpmin
        do i=1,nw
          do j=1,nh
            kk=(i-1)*nh+j
            if (zero(kk).gt.0.0005_dp.and.www(kk).lt.0.1_dp) then
              call seva2d(bkx,lkx,bky,lky,c,rgrid(i),zgrid(j),pds,ier,n333)
              write(6,*) 'rgrid(i),zgrid(j)',rgrid(i),zgrid(j)
              write(6,*) 'lkx, lky',lkx,lky
              write(6,*) 'pds,ier,n111', pds,ier,n333
              bpnow=sqrt(pds(2)**2+pds(3)**2)/rgrid(i)
              if (bpnow.le.bpmin) then
                if ((abs(dpsi).le.psitol).or.((abs(dpsi).gt.psitol).and. &
                    (zgrid(j)*zseps(1,jtime).lt.0.0))) then
                  bpmin=bpnow
                  xs=rgrid(i)
                  ys=zgrid(j)
                  sibpmin=pds(1)
                endif
              endif
            endif
          enddo
        enddo

        bmin_ave: if (bpmin.ne.avebp) then
          relsi=abs((sibpmin-psibry)/sidif)
          if (bpmin.le.0.10_dp*avebp.and.relsi.gt.0.005_dp) then
!------------------------------------------------------------------
!--         find second separatrix                               --
!------------------------------------------------------------------
            do j=1,40
              call seva2d(bkx,lkx,bky,lky,c,xs,ys,pds,ier,n666)
              write(6,*) 'xs,ys',xs, ys
              write(6,*) 'lkx, lky',lkx,lky
              write(6,*) 'pds,ier,n111', pds,ier,n111

              det=pds(5)*pds(6)-pds(4)*pds(4)
              if (abs(det).lt.1.0e-15_dp) exit
              xerr=(-pds(2)*pds(6)+pds(4)*pds(3))/det
              yerr=(-pds(5)*pds(3)+pds(2)*pds(4))/det
              xs=xs+xerr
              ys=ys+yerr
              if (xerr*xerr+yerr*yerr.lt.1.0e-12_dp) exit
            enddo
            if (xerr*xerr+yerr*yerr.ge.1.0e-12_dp) then
              epssep=xerr*xerr+yerr*yerr
              write (nttyo,11001) epssep,ixt
              if (iand(iout,1).ne.0) write (nout,11001) epssep,ixt
            else
              epssep=-1.0
            endif
            if (epssep.lt.1.0e-10_dp) then
              sibpmin=pds(1)
              yvs2=ys
              rsepex=xs
              relsi=abs((sibpmin-psibry)/sidif)
              if (relsi.gt.0.005_dp) then
                ! TODO: what is the intention here? sibpold is unset...
!                if (ixt.gt.1) sibpmin=sibpmin*(1.-vsdamp)+sibpold*vsdamp
                xvsmin=max(xvsmin,sibpmin)
              endif
            endif
          endif
        endif bmin_ave
        endif skipvs
        if (alphafp.ge.0.0) then
          xpsimin=xvsmin+alphafp*(xvsmax-xvsmin)
        else
          ! TODO: xpsipmins has not yet been defined...
!          xpsimino=xpsimins
          xpsimin=abs(alphafp)*xvsmax
          xpsimins=xpsimin
          ! TODO: what is the intention here? xpsimino is unset...
!          if (ixt.gt.1) xpsimin=xpsimin*(1.-vsdamp)+xpsimino*vsdamp
        endif
        xpsialp=xpsimin
        xpsimin=sidif/(simag-xpsimin)
      endif SOL_curr 
!-----------------------------------------------------------------------
!--   get normalized flux function XPSI                               --
!-----------------------------------------------------------------------
      do i=1,nw
       do j=1,nh
        kk=(i-1)*nh+j
        if (icutfp.eq.0) then
          xpsi(kk)=1.1_dp
          if((rgrid(i).lt.xmin).or.(rgrid(i).gt.xmax)) cycle
          if((zgrid(j).lt.ymin).or.(zgrid(j).gt.ymax)) cycle
          xpsi(kk)=(simag-psi(kk))/sidif
        else
          if (zero(kk).gt.0.0005_dp) then
            xpsi(kk)=(simag-psi(kk))/sidif
            if (xpsi(kk)*xpsimin.le.1.0.and.xpsi(kk)*xpsimin.ge.0.0) then
              if((rgrid(i).lt.rminvs).or.(rgrid(i).gt.rmaxvs)) &
                   xpsi(kk)=1000.
              if ((zgrid(j).lt.zminvs).or.(zgrid(j).gt.zmaxvs)) then
                   xpsi(kk)=1000.
              endif
              if(xpsi(kk).lt.1.0.and.zgrid(j).lt.ymin) xpsi(kk)=1000.
              if(xpsi(kk).lt.1.0.and.zgrid(j).gt.ymax) xpsi(kk)=1000.
              if(abs(zgrid(j)).gt.abs(yvs2).and. &
                zgrid(j)*yvs2.gt.0.0) xpsi(kk)=1000.
            endif
          else
            xpsi(kk)=1000.
          endif
        endif
       enddo
      enddo
!-----------------------------------------------------------------------
!--   get SOL flux if needed                                          --
!-----------------------------------------------------------------------
      if (nsol.gt.0) then
        call seva2d(bkx,lkx,bky,lky,c,rsol(1),zsol(1), &
                     pds,ier,n111)
        wsisol=pds(1)
#ifdef DEBUG_LEVEL2
        write (6,*) 'UPDATE R,Z,Si,Err = ', rsol(1),zsol(1),wsisol,ier 
        call seva2d(bkx,lkx,bky,lky,c,rbdry(1),zbdry(1),pds,ier,n111)
        write (6,*) 'UPDATE R,Z,Si,Err = ', rbdry(1),zbdry(1),pds(1),ier
        call seva2d(bkx,lkx,bky,lky,c,rbdry(nbdry),zbdry(nbdry),pds,ier &
             ,n111)
        write (6,*) 'UPDATE R,Z,Si,Err = ', rbdry(nbdry),zbdry(nbdry) &
            ,pds(1),ier
#endif
      endif
#ifdef DEBUG_LEVEL2
      call seva2d(bkx,lkx,bky,lky,c,xout(1),yout(1),pds,ier,n111)
      write (6,*) 'UPDATE simag,psibry,n1,si = ',simag,psibry,n111,pds(1)
      write (6,*) 'UPDATE Z,R = ', zgrid(33),(rgrid(i),i=45,45)
      write (6,*) 'UPDATE si = ',(psi((i-1)*nw+33),i=45,45)
      call seva2d(bkx,lkx,bky,lky,c,rgrid(45),zgrid(33),pds,ier,n111)
      write (6,*) 'UPDATE R,Z,si = ', rgrid(45),zgrid(33),pds(1)
      write (6,*) 'UPDATE lkx,lky = ',lkx,lky
#endif
!-----------------------------------------------------------------------
!--   get weighting function                                          --
!-----------------------------------------------------------------------
      call weight(rgrid,zgrid)
!-----------------------------------------------------------------------
!--   get response functions for MSE                                  --
!-----------------------------------------------------------------------
      if (kstark.gt.0.or.kdomse.gt.0) then
        do k=1,nstark
          if (rrgam(jtime,k).le.0.0) cycle
          call seva2d(bkx,lkx,bky,lky,c,rrgam(jtime,k) &
                      ,zzgam(jtime,k),pds,ier,n111)
          sistark(k)=pds(1)
          sisinow=(simag-pds(1))/sidif
          sigam(k)=sisinow
          fpnow=ffcurr(sisinow,kffcur)
          btgam(k)=fpnow*tmu/rrgam(jtime,k)
        enddo
      endif
!-----------------------------------------------------------------------
!--   get response functions for MSE-LS                               --
!-----------------------------------------------------------------------
      if (mmbmsels.gt.0.or.kdomsels.gt.0) then
        do k=1,nmsels
#ifdef DEBUG_LEVEL2
          write (6,*) 'UPDATE MSE-LS k,rrmselt= ', k,rrmselt(jtime,k)
#endif
          if (rrmselt(jtime,k).le.0.0) cycle
          call seva2d(bkx,lkx,bky,lky,c,rrmselt(jtime,k) &
                      ,zzmselt(jtime,k),pds,ier,n333)
          simls(k)=pds(1)
          sisinow=(simag-pds(1))/sidif
          sinmls(k)=sisinow
          fpnow=ffcurr(sisinow,kffcur)
          btmls(k)=fpnow*tmu/rrmselt(jtime,k)
          brmls(k)=-pds(3)/rrmselt(jtime,k)
          bzmls(k)=pds(2)/rrmselt(jtime,k)
#ifdef DEBUG_LEVEL2
          write (6,*) 'UPDATE MSE-LS k,rrmselt,br,bz,bt= ',  &
                      k,rrmselt(jtime,k),brmls(k),bzmls(k),btmls(k)
#endif
        enddo
      endif
!
      xouts(1:nfound)=1./xout(1:nfound)**2
      nzz=0
      call fluxav(xouts,xout,yout,nfound,psi,rgrid,nw,zgrid,nh, &
                  r2bdry,nzz,sdlobp,sdlbp)
      xouts(1:nfound)=1./xout(1:nfound)
      call fluxav(xouts,xout,yout,nfound,psi,rgrid,nw,zgrid,nh, &
                  r1bdry,nzz,sdlobp,sdlbp)
!-----------------------------------------------------------------------
!--   get metric elements at PSIWANT for edge constraint if needed    --
!-----------------------------------------------------------------------
      r1sdry(1)=r1bdry
      r2sdry(1)=r2bdry
      if (abs(sizeroj(1)-1.0).gt.1.e-05_dp.or.kzeroj.ne.1) then
        if (kzeroj.gt.0) then
          do i=1,kzeroj
            if(sizeroj(i).ge.1.0) sizeroj(i)=0.99999_dp
            siwant=simag+sizeroj(i)*(psibry-simag)
            call surfac(siwant,psi,nw,nh,rgrid,zgrid,rsplt,zsplt,nfounc, &
                        npoint,drgrid,dzgrid,xmin,xmax,ymin,ymax,nnn, &
                        rmaxis,zmaxis,negcur,kerror,1)
            if(kerror.gt.0) return
            csplt(1:nfounc)=1./rsplt(1:nfounc)**2
            call fluxav(csplt,rsplt,zsplt,nfounc,psi,rgrid,nw,zgrid,nh, &
                        r2sdry(i),nzz,sdlobp,sdlbp)
            csplt(1:nfounc)=1./rsplt(1:nfounc)
            call fluxav(csplt,rsplt,zsplt,nfounc,psi,rgrid,nw,zgrid,nh, &
                        r1sdry(i),nzz,sdlobp,sdlbp)
            r2surs = r2sdry(i)*sdlobp
            fpnow = ffcurr(psiwant,kffcur)
            fpnow = fpnow*tmu
          enddo
        endif
      endif
!-----------------------------------------------------------------------
!--   get metric elements at PSIWANT for q constraint if needed       --
!-----------------------------------------------------------------------
      if (nqwant.gt.0) then
        do i=1,nqwant
          siwant=simag+siwantq(i)*(psibry-simag)
          call surfac(siwant,psi,nw,nh,rgrid,zgrid,rsplt,zsplt,nfounc, &
                      npoint,drgrid,dzgrid,xmin,xmax,ymin,ymax,nnn, &
                      rmaxis,zmaxis,negcur,kerror,1)
          if(kerror.gt.0) return
          csplt(1:nfounc)=1./rsplt(1:nfounc)**2
          call fluxav(csplt,rsplt,zsplt,nfounc,psi,rgrid,nw,zgrid,nh, &
                      r2qdry,nzz,sdlobp,sdlbp)
          csplt(1:nfounc)=1./rsplt(1:nfounc)
          call fluxav(csplt,rsplt,zsplt,nfounc,psi,rgrid,nw,zgrid,nh, &
                      r1qdry,nzz,sdlobp,sdlbp)
          r2surq = r2qdry*sdlobp
          fpnow = ffcurr(siwantq(i),kffcur)
          fpnow = fpnow*tmu
          qsiw(i)= abs(fpnow)/twopi*r2surq
          pasmsw(i)=sdlbp/tmu/twopi
        enddo
      endif
!
      cvolp(ixt)=0.0
      carea=0.0
      xym=xout(1)*yout(1)
      xyma=yout(1)
      do i=2,nfound
        xyp=xout(i)*yout(i)
        xypa=yout(i)
        dx=xout(i)-xout(i-1)
        cvolp(ixt)=cvolp(ixt)+twopi*(xyp+xym)/2.0*dx
        carea=carea+(xyma+xypa)/2.0*dx
        xym=xyp
        xyma=xypa
      enddo
      carea=abs(carea)
      if (iconvr.ne.3) then
        if((ix.gt.1).or.(ixout.le.1)) call chisqr(jtime)
      endif
      cvolp(ixt)=abs(cvolp(ixt))*1.0e+06_dp
      volume(jtime)=cvolp(ixt)
      rout(jtime)=(xmax+xmin)/2.*100.
      aminor(jtime)=100.*(xmax-xmin)/2.0
      csimag(ixt)=simag
      csibry(ixt)=psibry
      crmaxi(ixt)=rmaxis*100.
      czmaxi(ixt)=zmaxis*100.
      cemaxi(ixt)=emaxis
      cchisq(ixt)=chisq(jtime)
      csumip(ixt)=sumip
      tratio(ixt)=cratio
!---------------------------------------------------------------------
!--   get beta and li for constraints                               --
!---------------------------------------------------------------------
      if (fli.gt.0.0.or.fbetan.gt.0.0) then
        call betsli(jtime,rgrid,zgrid,kerror)
        if(kerror.gt.0) return
      endif
!----------------------------------------------------------------------
!--   vertical stabilization information                             --
!----------------------------------------------------------------------
      if (itell.gt.0.and.isetfb.ge.0) then
        delzmm = zmaxis - zmaxis_last
        zmaxis_last=zmaxis
      endif
      if (isetfb.ne.0) then
        fb_plasma(jtime)=abs(brfb(1)*nfbcoil/ipmhd(jtime))
        delzmm = zmaxis - zmaxis_last
        zmaxis_last=zmaxis
      elseif (eelip.gt.2.25_dp .and. itell.eq.0) then
        delzmm = zmaxis - zmaxis_last
        zmaxis_last=zmaxis
      endif
!----------------------------------------------------------------------
!--   magnetic axis parameters if needed                             --
!----------------------------------------------------------------------
      if (((icinit.gt.0).and.(iconvr.ne.3).and.(ixout.le.1)).or.(icurrt.eq.4).or. &
          (((icurrt.eq.2).or.(icurrt.eq.5)).and.(ixt.le.1).and.(icinit.gt.0))) then
        nqend=1
        if(errorm.lt.0.1_dp.and.icurrt.eq.4) nqend=nqiter
        do i=1,nqend
          if (i.gt.1) then
            call currnt(2,jtime,2,kerror)
            if(kerror.gt.0) return
          endif

          fcentr=fbrdy**2+sidif*dfsqe
          if (fcentr.lt.0.0) fcentr=fbrdy**2
          fcentr=sqrt(fcentr)*fbrdy/abs(fbrdy)
          rdiml=rmaxis/rzero
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
          cqmaxi(ixt)=(emaxis**2+1.)*abs(fcentr)/twopi/emaxis &
                      /rmaxis**2/abs(cjmaxi)
          qmaxis=cqmaxi(ixt)
          if (icurrt.eq.4) then
            if(qenp.gt.0.0) enp=enp*qenp/qmaxis
            if(qemp.gt.0.0) emp=emp*qmaxis/qemp
            enf=enp
            emf=emp
          endif
        enddo
        return
      endif
      select case (icurrt)
      case default ! 1
        rdiml=rmaxis/srma
        cjmaxi=cratio*(sbeta*rdiml+2.*salpha/rdiml)/darea
        if(kvtor.gt.0) &
          cjmaxi=cjmaxi+cratio/darea*sbetaw*rdiml*(rdiml**2-1.)
      case (2,5)
        fcentr=ffcurr(x000,kffcur)
        if (kvtor.eq.0) then
          cjmaxi=(rmaxis*ppcurr(x000,kppcur) &
                +fpcurr(x000,kffcur)/rmaxis)*cratio/darea
        else
          cjmaxi=(sum(rjjjx(1:kppcur)*brsp(nfsum+1:nfsum+kppcur)) &
                 +sum(rjjfx(1:kffcur)*brsp(nfsum+kppcur+1:nfsum+kppcur+kffcur)) &
                 +sum(rjjwx(1:kwwcur)*brsp(nfnpcr+1:nfnpcr+kwwcur))) &
                 /darea
        endif
      case (3)
        ! continue
      ! case (4) handled by preceeding if statement
      end select
!
      cqmaxi(ixt)=(emaxis**2+1.)*abs(fcentr)/twopi/emaxis &
                   /rmaxis**2/abs(cjmaxi)
      qmaxis=cqmaxi(ixt)
      return
 2100 format(/,1x,'shot',i6,' at ',i6,' ms ','** Problem in BOUND **')
11001 format(/,1x,'** 2nd seperatrix **',2x,e10.3,2x,i4)
      end subroutine update_params

!**********************************************************************
!>
!!    weight computes the weighting function www
!!    
!!    @param x : R axis of the grid (nw points)
!!    @param y : Z axis of the grid (nh points)
!!
!**********************************************************************
      subroutine weight(x,y)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit none
      real*8, intent(in) :: x(nw),y(nh)
      integer*4 i,ileft,ins,iright,j,jbot,jtop,kk
      real*8 a,a1,a2,a3,a4,b,c,d,p1,p2,p3,p4,psil,xxw,yyw
!
      if (iweigh.le.0) then
        do kk=1,nwnh
          if ((xpsi(kk).lt.0.0).or.(xpsi(kk).gt.1.0)) then
            www(kk)=0.0
          else
            www(kk)=zero(kk)
          endif
        enddo
        return
      endif
!----------------------------------------------------------------------
!--   find a rectangle to search based on output from subroutine bound
!----------------------------------------------------------------------
      www=0.0
      do i=1,nw-1
        ileft=i
        if((x(i).le.xmin).and.(x(i+1).gt.xmin)) exit
      enddo
      do i=ileft,nw-1
        iright=i
        if((x(i).lt.xmax).and.(x(i+1).ge.xmax)) exit
      enddo
      do j=1,nh-1
        jbot=j
        if((y(j).le.ymin).and.(y(j+1).gt.ymin)) exit
      enddo
      do j=jbot,nh-1
        jtop=j
        if((y(j).lt.ymax).and.(y(j+1).ge.ymax)) exit
      enddo
      jtop=min(jtop,nh)
      jbot=max(jbot,1)
      ileft=max(ileft,1)
      iright=min(iright,nw)
      do i = ileft,iright
        do j = jbot,jtop
          kk = j+nh*(i-1)
          a = psi(kk-nh)
          b = psi(kk+nh)
          c = psi(kk-1)
          d = psi(kk+1)
          if(i.eq.1) a = 0.
          if(i.eq.nw) b = 0.
          if(j.eq.1) c = 0.
          if(j.eq.nh) d = 0.
          psil = psibry
          ins = 0
          p1 = 0.
          p2 = 0.
          p3 = 0.
          p4 = 0.
          a1 = 0.
          a2 = 0.
          a3 = 0.
          a4 = 0.
          if (a.ge.psil) then
            p1 = a-psil
          else
            ins = ins+1
            a1 = a-psil
          endif
          if (b.ge.psil) then
            p2 = b-psil
          else
            ins = ins+1
            a2 = b-psil
          endif
          if (c .ge. psil) then
            p3 = c-psil
          else
            ins = ins+1
            a3 = c-psil
          endif
          if (d.ge.psil) then
            p4 = d-psil
          else
            ins = ins+1
            a4 = d-psil
          endif
          ins = ins+1
          select case (ins)
          case (1)
            www(kk) = 1.
          case (2)
            xxw = (p1+p2+p3+p4)/4.
            yyw = (a1+a2+a3+a4)
            yyw = (yyw/(xxw-yyw))**2
            www(kk) = 1.-0.5_dp*yyw
          case (3)
            xxw = (p1+p2+p3+p4)
            yyw = (a1+a2+a3+a4)
            www(kk) = xxw/(xxw-yyw)
          case (4)
            xxw = (p1+p2+p3+p4)
            yyw = (a1+a2+a3+a4)/4.
            xxw = (xxw/(xxw-yyw))**2
            www(kk) = 0.5_dp*xxw
          case (5)
            www(kk) = 0.
          end select
          www(kk) = www(kk)*zero(kk)
        enddo
      enddo
!      www = www*zero
      return
      end subroutine weight

!**********************************************************************
!>
!!    update the computed measurment values and calculate the chisq
!!    figure of merit for fitting
!!
!!    @param jtime : time index
!!
!**********************************************************************
      subroutine chisqr(jtime)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none

      integer*4, intent(in) :: jtime
      integer*4 i,m
      real*8 cm,cmbr,cmbz,ce1rbz,ce2rbz,ce3rbr,saisq
      integer*4, parameter :: nnn=0,nsq=2
!
#ifdef DEBUG_LEVEL2
      write (6,*) 'CHISQR, jtime= ',jtime
#endif
!
      saisq=0.0
      chi2rm(jtime)=0.0
      do m=1,nsilop
        cm=sum(rsilfc(m,:)*brsp(1:nfsum))+sum(gsilpc(m,:)*pcurrt)
        if(ivesel.gt.0) cm=cm+sum(rsilvs(m,:)*vcurrt)
        if (iecurr.eq.1) then
          cm=cm+sum(rsilec(m,:)*ecurrt)
        elseif (iecurr.eq.2) then
          cm=cm+sum(rsilec(m,:)*cecurr)
        endif
        if(iacoil.gt.0) cm=cm+sum(rsilac(m,:)*caccurt(jtime,:))
        if(fitsiref)    cm=cm-csiref
        if (swtsi(m).ne.0.0) then
          saisil(m)=(fwtsi(m)/swtsi(m))**nsq*(silopt(jtime,m)-cm)**2
        else
          saisil(m)=0.0
        endif
        saisq=saisq+saisil(m)
        csilop(m,jtime)=cm
        chi2rm(jtime)=max(chi2rm(jtime),saisil(m))
      enddo
!
      do m=1,magpri
          cm=sum(rmp2fc(m,:)*brsp(1:nfsum))+sum(gmp2pc(m,:)*pcurrt)
        if(ivesel.gt.0) cm=cm+sum(rmp2vs(m,:)*vcurrt)
        if (iecurr.eq.1) then
          cm=cm+sum(rmp2ec(m,:)*ecurrt)
        elseif (iecurr.eq.2) then
          cm=cm+sum(rmp2ec(m,:)*cecurr)
        endif
        if(iacoil.gt.0) cm=cm+sum(rmp2ac(m,:)*caccurt(jtime,:))
        if (swtmp2(m).ne.0.0) then
          saimpi(m)=(fwtmp2(m)/swtmp2(m))**nsq*(expmpi(jtime,m)-cm)**2
        else
          saimpi(m)=0.0
        endif
        saisq=saisq+saimpi(m)
        cmpr2(m,jtime)=cm
        chi2rm(jtime)=max(chi2rm(jtime),saimpi(m))
      enddo
!-------------------------------------------------------------------
!--   calculate ECE chisqr (chiece, tchiece for psi(R+)=psi(R-))
!-------------------------------------------------------------------
      tchiece=0.0
      do m=1,nece
        cm=sum(recefc(m,:)*brsp(1:nfsum)) &
          +sum(recepc(m,1:kwcurn)*brsp(nfsum+1:nfsum+kwcurn))
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
!-------------------------------------------------------------------
!--   calculate ECE chisqr (chiecebz for Bz(receo)=0)
!-------------------------------------------------------------------
      cm=sum(recebzfc*brsp(1:nfsum)) &
        +sum(recebzpc(1:kwcurn)*brsp(nfsum+1:nfsum+kwcurn))
      if(ivesel.gt.0) cm=cm+sum(recevs(m,:)*vcurrt) ! TODO: should this be recebzvs?
      if (iecurr.eq.1) then
        cm=cm+sum(recebzec*ecurrt)
      elseif (iecurr.eq.2) then
        cm=cm+sum(recebzec*cecurr)
      endif
      if(iacoil.gt.0) cm=cm+sum(receac(m,:)*caccurt(jtime,:)) ! TODO: should this be recebzac?
      if (swtecebz.ne.0.0) then
        chiecebz=fwtecebz**nsq*(brspecebz(jtime)-cm)**2
        chiecebz=chiecebz/swtecebz**nsq
      else
        chiecebz=0.0
      endif
      cmecebz(jtime)=cm
!
      cm=sum(pcurrt)
      ipmhd(jtime)=cm
      if(ivesel.gt.0) cm=cm+sum(vcurrt)
      if (swtcur.ne.0.0) then
        chipasma=(fwtcur/swtcur)**nsq*(ipmeas(jtime)-cm)**2
      else
        chipasma=0.0
      endif
      saisq=saisq+chipasma
      chi2rm(jtime)=max(chi2rm(jtime),chipasma)
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
        chi2rm(jtime)=max(chi2rm(jtime),chifcc(i))
      enddo
      if (iecurr.eq.2) then
        do i=1,nesum
          chiecc(i)=0.0
          if (fwtec(i).gt.0.0) then
            chiecc(i)=fwtec(i)**nsq*(cecurr(i)-ecurrt(i))**2
            chiecc(i)=chiecc(i)/swtec(i)**nsq
          endif
          saisq=saisq+chiecc(i)
          chi2rm(jtime)=max(chi2rm(jtime),chiecc(i))
        enddo
      endif
      if (fitsiref) then
        saisref=0.0
        if (fwtref.gt.0.0) then
          saisref=fwtref**nsq*(psiref(jtime)-csiref)**2
          saisref=saisref/swtsi(nslref)**nsq
        endif
        saisq=saisq+saisref
        chi2rm(jtime)=max(chi2rm(jtime),saisref)
      endif
!
      ccbrsp(:,jtime)=brsp(1:nfsum)
!
      chisq(jtime)=saisq
      if (iand(iout,1).ne.0) then
        write(nout,7400) time(jtime),chisq(jtime),ipmhd(jtime)
        write(nout,7420)
        write(nout,7450) (saisil(m),m=1,nsilop)
        write(nout,7430)
        write(nout,7450) (saimpi(m),m=1,magpri)
        write(nout,7460) chipasma
        write(nout,7480) tchifcc
        write(nout,7450) (chifcc(m),m=1,nfsum)
        write(nout,7482) saisref
        write(nout,7485)
        write(nout,7450) (chiecc(m),m=1,nesum)
        if(kprfit.gt.0) write(nout,7470) chipre
        if(kprfit.gt.0) write(nout,7450) (saipre(m),m=1,npress)
        if(kecebz.gt.0) write(nout,7486) chiecebz
        if(kece.gt.0) write(nout,7487) tchiece
        if(kece.gt.0) then
          write(nout,7488)
          write(nout,7450) (chiece(m),m=1,nece)
        endif
      endif
!-------------------------------------------------------------------------
!--   compute signals at MSE locations if requested                     --
!-------------------------------------------------------------------------
      if (kdomse.gt.0) then
        call green(nnn,jtime,n222)
        do m=1,nstark
          cmbr=sum(rbrfc(m,:)*brsp(1:nfsum)) &
              +sum(rbrpc(m,1:kwcurn)*brsp(nfsum+1:nfsum+kwcurn))
          cmbz=sum(rbzfc(m,:)*brsp(1:nfsum)) &
              +sum(rbzpc(m,1:kwcurn)*brsp(nfsum+1:nfsum+kwcurn))
          if (kedgep.gt.0) then
            cmbr=cmbr+rbrpe(m)*pedge
            cmbz=cmbz+rbzpe(m)*pedge
          endif
          if (kedgef.gt.0) then
            cmbr=cmbr+rbrfe(m)*f2edge
            cmbz=cmbz+rbzfe(m)*f2edge
          endif
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
          cm=a2gam(jtime,m)*btgam(m)+a3gam(jtime,m)*cmbr &
                                    +a4gam(jtime,m)*cmbz
          bzmsec(m)=cmbz
          if (keecur.le.0) then
            cm=a1gam(jtime,m)*cmbz/cm
          else
            ce1rbz=sum(e1rbz(m,1:keecur)*cerer(1:keecur))
            ce2rbz=sum(e2rbz(m,1:keecur)*cerer(1:keecur))
            ce3rbr=sum(e3rbr(m,1:keecur)*cerer(1:keecur))
            cm=cm-ce2rbz-ce3rbr
            cm=(a1gam(jtime,m)*cmbz-ce1rbz)/cm
          endif
          cmgam(m,jtime)=cm
        enddo
      endif
!-------------------------------------------------------------------------
!--   compute signals at MSE-LS locations if requested                  --
!-------------------------------------------------------------------------
      if (kdomsels.gt.0) then
        call green(nnn,jtime,n222)
        do m=1,nmsels
          cmbr=sum(rmlspc(m,1:kpcurn)*brsp(nfsum+1:nfsum+kpcurn))
          cmbr=cmbr-rhsmls(jtime,m)
          cmmls(jtime,m)=sqrt(cmbr)
          cmmls2(jtime,m)=l1mselt(jtime,m)*btmls(m)**2+l2mselt(jtime,m)* &
             brmls(m)*btmls(m)+l3mselt(jtime,m)*brmls(m)**2+             &
            (l1mselt(jtime,m)+l3mselt(jtime,m))*bzmls(m)**2
          cmmls2(jtime,m)=sqrt(cmmls2(jtime,m))
        enddo
        write(nout,7585)
        write(nout,7450) (cmmls(jtime,m),m=1,nmsels)
        write(nout,7588)
        write(nout,7450) (cmmls2(jtime,m),m=1,nmsels)
      endif
!
      if (mmbmsels.gt.0.or.kdomsels.gt.0) then
        do m=1,nmsels
          cmmlsv(jtime,m)=0.0
          if (rrmselt(jtime,m).gt.0.0) then
            cmmlsv(jtime,m)=abs(bcentr(jtime)*rcentr/rrmselt(jtime,m))
            cmmlsv(jtime,m)=cmmlsv(jtime,m)*sqrt(l1mselt(jtime,m))
          endif
        enddo 
        write(nout,7590)
        write(nout,7450) (cmmlsv(jtime,m),m=1,nmsels)
      endif
!
      return
 7400 format (/,2x,'time = ',e12.5,2x,'chisq = ',e12.5,2x,'current = ',e12.5)
 7420 format (10x,'chi psi loops:')
 7430 format (10x,'chi inner magnetic probes:')
 7450 format (8(1x,e12.5:,1x))
 7460 format (10x,'chi ip:',/,15x,e12.5)
 7470 format (10x,'chi pressure:         ',/,1x,e12.5)
 7480 format (10x,'chi F-coils:          ',/,10x,e12.5)
 7482 format (10x,'chi psiref:',/,15x,e12.5)
 7485 format (10x,'chi E-coils:          ')
 7486 format (10x,'chi ecebz:            ',/,1x,e12.5)
 7487 format (10x,'chi total eceR+R-:    ',/,1x,e12.5)
 7488 format (10x,'chi eceR+R-:          ')
 7585 format (10x,'Simulated MSE-LS (T): ',/)
 7588 format (10x,'Simulated MSE-LS2 (T): ',/)
 7590 format (10x,'Simulated MSE-LSV (T): ',/)
      end subroutine chisqr
