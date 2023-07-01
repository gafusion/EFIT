#include "config.f"
!**********************************************************************
!>
!!    shape finds the outermost plasma surface
!!    and computes various global plasma parameters.
!!    
!!
!!    @param iges : time index
!!
!!    @param igmax : number of time slices
!!
!!    @param kerror : error flag
!!
!**********************************************************************
      subroutine shapesurf(iges,igmax,kerror)
      use commonblocks,only: c,wk,copy,bkx,bky,psiold,psipold, &
                worka,zeros,byringr,byringz,xouts,youts,bpoo,bpooz, &
                bpooc,bfpol,cfpol,dfpol,xxtra,yxtra,bpxtra,flxtra,fpxtra
      use efit_bdata,only: xlims,ylims,xlmins
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none

      integer*4, intent(in) :: iges,igmax
      integer*4, intent(out) :: kerror
      real*8 ppcurr,fpcurr,fpecrr,speval,pwcurr,prcurr,ffcurr,seval, &
             esradial
      integer*4 i,ii,i2,j,jj,jjj,k,kk,kkk,ks,m,n,kacerr,floorq, &
                iautoc,ichisq,ichisq2,idoqn,iend,ier,inow,ip1,iring, &
                iskip,ixl,ixtpls,ixyz,iznow,jb,jges,jstart,kim,kip, &
                kjm,kjp,kz,m20,mcentral,n22,nbabs,nerr,nfind,nfounc, &
                nfouns,njtwagap,nm1,nmax,nmaxfs,nmin,nminfs,nnn,nzz
      real*8 floorz,abar,aream,areart,aspect,bincp,bp2flx, &
             bpnow,bpolsq,bpolzs,bpzsq,btnow,btttt2,btuse,btvac2, &
             cj1now,cjorka,cons1,cons2,cons3,cons4,cons5,cons6,cons7, &
             cons8,cons9,cons10,cons11,cons12,cons13,cons14,cons15, &
             cons16,cons17,cons21,cons22,consa1,consa2,consa3,cosalp, &
             crymax,crymin,curnow,d11,d22,d33,d2sidr2,d2sidz2,dbott, &
             delarea,delerrb,delerrx,delfp,delli,dells,delr,delrnow,delsbu, &
             delsi,delssi,deltaa,dilnow,dis2p,dismin,disnow,dleft,dli, &
             dlll,dlpol,dltol,dmaxfs,dmaxfs0,dminfs,dminfs0,dminow,dmui, &
             dolnow,atri,dpsis,dright,dsidr,dsiin,dsilim,dsimm,dsiout, &
             dsmin,dtop,dttt,dx,dxtra,dyww,dzring,dzzz1,dzzz2,enrgy, &
             exmui,f_0,fkdel,fmanow,fminow,fnow,fpnow,fsrmax,fszmax, &
             fvsmax,fvsnow,fxmax,fxmin,fxrmax,fxrmin,fxzmax,fxzmin, &
             gap1,gap2,olamda,partic,pasnow,pcnow,pin,pleng,pohm,ppnow, &
             pres0,presped,press0,presss,prettt,prew0,preww0,prewww, &
             psiin,psimm,psinow,psiots,psiout,psiwan,ptotal,pwop0,pwp0r2, &
             qmmm,qppp,qpsic,qpwant,qwant,r2surs,radbou,radp,radum, &
             ravssa,rbar,rcur,rcurrm,resist,rexpmx,rinvs,rkkk,rlinnose,&
             rmer,rmm,rnow,rolnow,routvs,rqmax,rqmin,rrin,rrout,rsepnose, &
             rsnow,rsnow0,rsymax,rsymin,rtemp,rval,rvsnow,rx,rxp,rxxrry, &
             rymaxs,rymins,sbli,sbpli,sbrrs,sdlbp,sdlobp,sepnow,siavej, &
             sigamnow,signr,signz,siii,silimm,silimp,silop_change,sinalp, &
             siwant,siwwww,slope,ssimax,ssimfs,ssinow,ssitra,sssiie, &
             sumbp2,sumbz2,sumf,sumfzp,tevolt,vbpli,volbt,volrt,x11, &
             xdum,xguess,xlam,xlimxs,xlimxx,xmaxs,xmaxx,xminn,xmins, &
             xmnow,xmui,xoutm,xoutp,xpsikk,xpsivs,xrmax,xrmin,xrpres, &
             xsepsl,xsi01,xsi95,xsiww,xsiwww,xww,xxtraa,xxtras,xxxx, &
             xycmax,xycmin,xym,xyma,xyp,xypa,ycmax,ycmin,ycut,ycutm,ydum, &
             yguess,ylimys,ylimyy,ymaxs,ymins,ymnow,yoxm,yoxp,ypsz, &
             yww,yxcmax,yxcmin,yxtraa,yxtras,zavssa,zbar,zcur,zerold, &
             zeta,zexpmx,zhp,zilnow,zinvs,zkkk,znow,zoutvs,zqmax, &
             zringmax,zringmin,zrmin,zsnow,zsnow0,ztemp,zval,zvsnow, &
             zxmins,zxmaxs,zxp,zxx,zzm,zzp
      integer*4 imer(2)
      real*8 pds(6),amer(2,2),bmer(2),wmer(2),temp(ntime)
      real*8 rmid2(2),zerovs(1),ravs(1),zavs(1)
      real*8 sigams(nstark)
      real*8 ringr(6),ringz(6),ringap
      real*8, dimension(:), allocatable :: xsisii,bpres,cpres, &
                    dpres,sjtli,sjtlir,sjtliz,rjtli,bpresw, &
                    cpresw,dpresw,copyn,cjtli,x,y
      character(30) sfname,ofname
      Character(28) xxtitle,yytitle,zztitle
      character(20) zzztitle
      Character(8) jchisq
      character(1) jchisq2
      logical byring,double,onedone
      integer*4, parameter :: idiart=1,limtrs=5
      real*8, parameter :: psitol=1.0e-04_dp,czero=0.0, &
                           rubaf=1.372_dp,zubaf=1.310_dp,&
                           rudom=1.0420_dp,zudom=1.1624_dp
      data floorz/-1.366_dp/
      data double/.false./,onedone/.false./
!---D3D-----------------------D3D----------------------------D3D-----
      data ringr/1.766,1.680,1.674,2*1.671,1.681/
      data ringz/2*-1.255,-1.258,-1.264,-1.327,-1.335/
!---D3D-----------------------D3D----------------------------D3D-----
!
      ALLOCATE(xsisii(nw),bpres(nw),cpres(nw),dpres(nw), &
         sjtli(nw),sjtlir(nw),sjtliz(nw),rjtli(nw), &
         bpresw(nw),cpresw(nw),dpresw(nw),copyn(nwnh), &
         cjtli(nw),x(nw),y(nh))
!-----------------------------------------------------------------------
!-- ringr and ringz define plasma facing surfaces where strike point  --
!-- may contact surface                                               --
!-- byring is true if r>rsep, z<zsep and z>zvs                        --
!-- double refers to double-valued q-profile                          --
!-- onedone=true if the location (psi) of one q=1 surface has been    --
!-- found                                                             --
!-----------------------------------------------------------------------
      kerror = 0
      oring(iges)=999 ! initialize at ridiculous value
      ringap=999
      zavs(1)=0.0
      rhovn=0.0
      betapw(iges)=0.0
      betatw(iges)=0.0
      wplasw(iges)=0.0
      tave(iges)=0.0
!
      xdum=0.0
      ydum=0.0
#ifdef DEBUG_LEVEL1
      write (6,*) 'Enter SHAPE iges, ivacum = ', iges, ivacum
#endif
      is_vacuum: if (ivacum.eq.0) then
      if (iges.gt.1) then
        silop_change=maxval(abs(silopt(iges,:)-silopt(iges-1,:)))
      else
        silop_change=HUGE(siloplim)
      endif
      if (silop_change.ge.siloplim) then
        xguess=(rgrid(1)+rgrid(nw))/2.
        yguess=(zgrid(1)+zgrid(nh))/2.
        xlims(1)=rgrid(2)+0.2_dp*(rgrid(2)-rgrid(1))
        xlims(2)=xlims(1)
        xlims(3)=rgrid(nw-1)-0.2_dp*(rgrid(nw)-rgrid(nw-1))
        xlims(4)=xlims(3)
        xlims(5)=xlims(1)
        ylims(1)=zgrid(2)+0.2_dp*(zgrid(2)-zgrid(1))
        ylims(2)=zgrid(nh-1)-0.2_dp*(zgrid(nh)-zgrid(nh-1))
        ylims(3)=ylims(2)
        ylims(4)=ylims(1)
        ylims(5)=ylims(1)
        xlmins=xlims(1)
        call zlim(zeros,nw,nh,limtrs,xlims,ylims,rgrid,zgrid,limfag)
      endif
!----------------------------------------------------------------------
!--   get outermost flux surface and its shape parameters            --
!----------------------------------------------------------------------
      jges=iges
      itime=time(iges)
      sibdry(iges)=psibry
      psim(iges)=simag
      rout(iges)=(xmin+xmax)/2.0
      zout(iges)=(ymin+ymax)/2.0
      elong(iges)=(ymax-ymin)/(xmax-xmin)
      aminor(iges)=100.*(xmax-xmin)/2.0
      rexpmx=xmax+rexpan
      zexpmx=zxmax
      zzmax(nw)=ymax
      rzzmax(nw)=rymax
      xminn=xmin
      xmaxx=xmax
      xoutp=xout(2)-xout(1)
      xout(nfound+1)=xout(2)
      yout(nfound+1)=yout(2)
      do i=2,nfound
        xoutm=xoutp
        xoutp=xout(i+1)-xout(i)
        if (xoutp*xoutm.ge.0.0) cycle
        if (xoutp.le.0.) then
          if (abs(xout(i)-xmax).le.1.0e-04_dp) cycle
          if (abs(yout(i)-zout(iges)).ge.abs(zxmin-zout(iges))) cycle
          if (xout(i).ge.rout(iges)) cycle
          xminn=xout(i)
        else
          if (abs(xout(i)-xmin).le.1.0e-04_dp) cycle
          if (abs(yout(i)-zout(iges)).ge.abs(zxmax-zout(iges))) cycle
          if (xout(i).le.rout(iges)) cycle
          xmaxx=xout(i)
        endif
      enddo
      indent(iges)=(xminn-xmin+xmax-xmaxx)/2./aminor(iges)*100.
      rout(iges)=100.*rout(iges)
      zout(iges)=100.*zout(iges)
!-----------------------------------------------------------------------
!--   the distance to the top limiter is found only for values of yout(i)
!--   which are greater than yulim below.
!-----------------------------------------------------------------------
      crymin=100.*rymin
      crymax=100.*rymax
      utri(iges)=(rout(iges)-crymax)/aminor(iges)
      ltri(iges)=(rout(iges)-crymin)/aminor(iges)
!---------------------------------------------------------------------
!--   set up P' and FF', then integration                           --
!--   ffprim = (RBt) * d/dpsi(RBt)                                  --
!---------------------------------------------------------------------
      select case (icurrt)
      case (1)
        pprime(1)=cratio*sbeta/darea/srma
        ffprim(1)=cratio*srma*2.*salpha/darea*twopi*tmu
        pprime(nw)=pprime(1)
        ffprim(nw)=ffprim(1)
      case (2,5)
        pprime(nw)=ppcurr(x111,kppcur)/darea
        ffprim(nw)=fpcurr(x111,kffcur)/darea*twopi*tmu
        pprime(1)=ppcurr(x000,kppcur)/darea
        ffprim(1)=fpcurr(x000,kffcur)/darea*twopi*tmu
        if (kfffnc.eq.8) then
          ffprec(nw)=fpecrr(x111,kffcur)/darea*twopi*tmu
          ffprec(1)=fpecrr(x000,kffcur)/darea*twopi*tmu
        else
          ffprec(nw)=0.0
          ffprec(1)=0.0
        endif
      case (4)
        call currnt(n222,iges,n222,kerror)
        if (kerror.gt.0) return
        pprime(1)=cratio/darea/rzero
        ffprim(1)=rbetap*cratio*rzero*twopi*tmu/darea
        ffprim(nw)=ffprim(1)*gammaf
        pprime(nw)=pprime(1)*gammap
      end select
!
      do i=2,nw-1
        ii=nw-i+1
        siii=1.0_dp-1.0_dp/(nw-1)*(i-1)
        sigrid(ii)=siii
        select case (icurrt)
        case (1)
          pprime(ii)=pprime(1)
          ffprim(ii)=ffprim(1)
        case (2,5)
          pprime(ii)=ppcurr(siii,kppcur)/darea
          ffprim(ii)=fpcurr(siii,kffcur)/darea*twopi*tmu
          if (kfffnc.eq.8) then
            ffprec(ii)=fpecrr(siii,kffcur)/darea*twopi*tmu
          else
            ffprec(ii)=0.0
          endif
        case (4)
          pprime(ii)=(1.-siii**enp)**emp*(1.-gammap)+gammap
          ffprim(ii)=ffprim(1)*pprime(ii)
          pprime(ii)=pprime(1)*pprime(ii)
        end select
      enddo
!
      sigrid(1)=0.0
      sigrid(nw)=1.0
      pres(nw)=prbdry
      fpol(nw)=fbrdy*tmu
      sumf=fpol(nw)**2/2.
      delsi=-(psibry-simag)/(nw-1)
      do i=1,nw-1
        pres(nw-i)=pres(nw-i+1)+0.5_dp*(pprime(nw-i+1)+pprime(nw-i))*delsi
        sumf=sumf+0.5_dp*delsi*(ffprim(nw-i+1)+ffprim(nw-i))
        if (sumf .ge. 0.0) then
          fpol(nw-i)=sqrt(2.*sumf)*fpol(nw)/abs(fpol(nw))
        else
          fpol(nw-i)=fpol(nw)
        endif
      enddo
      call zpline(nw,sigrid,fpol,bbfpol,ccfpol,ddfpol)
!
      volume(iges)=0.0
      rhovn(nw)=0.0
      xym=xout(1)*yout(1)
      yoxm=yout(1)/xout(1)
!------------------------------------------------------------------
!--   integration over z from 0 to bpolz                         --
!------------------------------------------------------------------
      xww=xout(1)
      dyww=yout(1)/(nh-1)
      bpolzs=0.5_dp*fpol(nw)
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
      bpolzs=bpolzs+fnow*0.5_dp
      yoxm=bpolzs/xout(1)*dyww
!
      aream=0.0
      xyma=yout(1)
      zzm=zmaxis-yout(1)
      do i=2,nfound
        xyp=xout(i)*yout(i)
!------------------------------------------------------------------
!--     integration over z from 0 to bpolz                       --
!------------------------------------------------------------------
        xww=xout(i)
        dyww=yout(i)/(nh-1)
        bpolzs=0.5_dp*fpol(nw)
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
        bpolzs=bpolzs+fnow*0.5_dp
        yoxp=bpolzs/xout(i)*dyww
!
        xypa=yout(i)
        dx=xout(i)-xout(i-1)
        volume(iges)=volume(iges)+(xyp+xym)/2.0*dx
        aream=aream+(xyma+xypa)/2.0*dx
        rhovn(nw)=rhovn(nw)+(yoxp+yoxm)/2.0*dx
        xym=xyp
        xyma=xypa
        yoxm=yoxp
        zzp=zmaxis-yout(i)
        if (zzp*zzm.le.0.0) then
          slope=(xout(i)-xout(i-1))/(yout(i)-yout(i-1))
          if(xout(i).lt.rmaxis) rminzm=xout(i)+zzp*slope
          if(xout(i).gt.rmaxis) rmaxzm=xout(i)+zzp*slope
        endif
        zzm=zzp
      enddo
      volume(iges)=abs(volume(iges))*1.0e+06_dp*twopi
      aream=abs(aream)
      area(iges)=aream*1.0e+04_dp
      rm(iges)=rmaxis*100.0
      zm(iges)=zmaxis*100.0
      elongm(iges)=emaxis
      qm(iges)=qmaxis
      tflux(iges)=rhovn(nw)
!---------------------------------------------------------------------
!--   gap calculation                                               --
!---------------------------------------------------------------------
      dleft=1.0e+10_dp
      dright=1.0e+10_dp
      dtop=1.0e+10_dp
      dbott=1.0e+10_dp
      do j=1,limitr-1
        call dslant(xout,yout,nfound,xmin,xmax,ymin,ymax, &
                    xlim(j),ylim(j),xlim(j+1),ylim(j+1),disnow)
        if(xlim(j).lt.1.02_dp .and. xlim(j+1).lt.1.02_dp) &
          dleft = min(dleft,disnow)
        if(ylim(j).gt.1.20_dp .and. ylim(j+1).gt.1.20_dp) &
          dtop = min(dtop,disnow)
        if(ylim(j).lt.-1.20_dp .and. ylim(j+1).lt.-1.20_dp) &
          dbott = min(dbott,disnow)
        if(xlim(j).gt.1.70_dp .and. xlim(j+1).gt.1.70_dp) &
          dright = min(dright,disnow)
      enddo
      dismin=min(dleft,dright,dtop,dbott)
      gapin(iges)=dleft
      gapout(iges)=dright
      gaptop(iges)=dtop
      gapbot(iges)=dbott
      dsep(iges)=dismin
      drsep(iges)=40.
      if (abs(dismin).le.0.100_dp) then
        if(dleft.eq.dismin) limloc(iges)='IN '
        if(dright.eq.dismin) limloc(iges)='OUT'
        if(dtop.eq.dismin) limloc(iges)='TOP'
        if(dbott.eq.dismin) limloc(iges)='BOT'
      else
!--------------------------------------------------------------------
!--     diverted configuration                                     --
!--------------------------------------------------------------------
        deltaa=0.0025_dp
        if (delrmax1.lt.deltaa.and.delrmax2.lt.deltaa) then
           limloc(iges)='DN '
           drsep(iges)=max(delrmax1,delrmax2)*100.
           if (zseps(1,iges).lt.zout(iges)) then
             drsep(iges)=-drsep(iges)
           endif
        else
           if (delrmax1.lt.deltaa) then
            if (zseps(1,iges).gt.zout(iges)) then
             limloc(iges)='SNT'
             drsep(iges)=abs(delrmax2)*100.
            else
             limloc(iges)='SNB'
             drsep(iges)=-abs(delrmax2)*100.
            endif
           else
             drsep(iges)=delrmax1*100.
             limloc(iges)='MAR'
           endif
        endif
      endif
!---------------------------------------------------------------------
!--   Helicon gap calculation                                       --
!---------------------------------------------------------------------
      twagap(iges)=1.0e+10_dp
      if (ishot.ge.139282) then
        jtwagap = 59
        njtwagap = 0
        if (ishot.ge.181292) then
          jtwagap = 47
          njtwagap = 1
        elseif (ishot.ge.187873) then
          jtwagap = 48
          njtwagap = 1
        endif
        do j=jtwagap,jtwagap+njtwagap
          call dslant(xout,yout,nfound,xmin,xmax,ymin,ymax, &
                      xlim(j),ylim(j),xlim(j+1),ylim(j+1),disnow)
          twagap(iges) = min(twagap(iges),disnow)
        enddo
      endif
!
      xlimxs=0.0
      ylimys=0.0
      if (dismin.ge.0.500_dp) then
        xsepsl=100.
        if(zseps(1,iges).lt.0.0) xsepsl=rseps(1,iges)/100.
        if(zseps(2,iges).lt.0.0) xsepsl=rseps(2,iges)/100.
        call seva2d(bkx,lkx,bky,lky,c,xlim(1),ylim(1),pds,ier,n111)
        silimp=pds(1)-psibry
        do i=2,limitr
          call seva2d(bkx,lkx,bky,lky,c,xlim(i),ylim(i),pds,ier,n111)
          silimm=silimp
          silimp=pds(1)-psibry
          if(silimp*silimm.gt.0.0) cycle
          dsilim=silimp-silimm
          xlimxx=xlim(i-1)-(xlim(i)-xlim(i-1))/dsilim*silimm
          ylimyy=ylim(i-1)-(ylim(i)-ylim(i-1))/dsilim*silimm
          if(ylimyy.ge.0.0) cycle
          if(xlimxx.lt.xsepsl) cycle
          xlimxs=xlimxx
          ylimys=ylimyy
          exit
        enddo
      endif
!----------------------------------------------------------------------
!--   find separatrix outside a limited plasma if one exists.        --
!--   if the plasma is diverted skip this calculation.               --
!--   the distances are defaulted to 99.0 cm.                        --
!--   DPSI > PSITOL  diverted plasma                                 --
!--          1.e10   well diverted                                   --
!----------------------------------------------------------------------
      m20=-20
      sepin(iges)=-50.0
      sepout(iges)=-50.0
      septop(iges)=-50.0
      sepbot(iges)=-50.0
      if (itrace.le.0) go to 1085
      if(abs(dpsi).gt.psitol.or.dismin.gt.0.1_dp) go to 1085
      radbou=1.10_dp*xmin
      if (ipmeas(iges).lt.-1.e3_dp) then
         nerr=10000
      else
         nerr=0
      endif
      call bound(psi,nw,nh,nwnh,psiots,xmins,xmaxs,ymins,ymaxs,zeros, &
                 rgrid,zgrid,xguess,yguess,jges,limtrs,xlims,ylims, &
                 xouts,youts,nfouns,xlmins,npoint,rymins,rymaxs,dpsis, &
                 zxmins,zxmaxs,nerr,ishot,itime,limfag,radbou,kbound, &
                 tolbndpsi)
      if (nerr.gt.0) then
        kerror = 1
        return
      end if
      !if (nerr.gt.0) then
      !  sepin(iges)=-89.0
      !  sepout(iges)=-89.0
      !  septop(iges)=-89.0
      !  sepbot(iges)=-89.0
      !  dsep(iges)=-89.0
      !  go to 1085
      !endif
      if (abs(dpsis).le.psitol) then
        sepin(iges)=-45.0
        sepout(iges)=-45.0
        septop(iges)=-45.0
        sepbot(iges)=-45.0
        dsep(iges)=-45.0
        go to 1085
      endif
      call findax(nw,nh,rgrid,zgrid,rmaxis,zmaxis,simag, &
                  psiots,rseps(:,iges),zseps(:,iges),m20, &
                  xouts,youts,nfouns,psi,xmins,xmaxs,ymins,ymaxs, &
                  zxmins,zxmaxs,rymins,rymaxs,dpsis,bpoo,bpooz, &
                  limtrs,xlims,ylims,limfag,0,0,kerror)
      if(kerror.gt.0) return
!---------------------------------------------------------------------
!--   gap calculation                                               --
!---------------------------------------------------------------------
      dsep(iges)=1000.
      do j=1,nfouns-1
        call dslant(xout,yout,nfound,xmin,xmax,ymin,ymax, &
                    xouts(j),youts(j),xouts(j+1),youts(j+1),disnow)
        dsep(iges)=-min(abs(dsep(iges)),disnow)
      enddo
!----------------------------------------------------------------------
!--   start distance calculation                                     --
!----------------------------------------------------------------------
      ! TODO: xleft, zleft, xright, zright, xztop, ytop, xzbot, and ybot
      !       are never defined in efit... this affects the computation
      !       of sepin, sepout, septop, and sepbot
      zhp=zout(iges)*0.01_dp
      radp=rout(iges)*0.01_dp
      do j=1,nfouns-1
       if (xouts(j).le.radp) then
!        dxll=(youts(j)-zleft)*(youts(j+1)-zleft)
!        if (dxll.le.0.0) then
!        if (youts(j).eq.zleft) sepin(iges)=(xouts(j)-xleft)*100.0
!        if (youts(j+1).eq.zleft) sepin(iges)=(xouts(j+1)-xleft)*100.0
!        if (xouts(j).eq.xouts(j+1)) sepin(iges)=(xouts(j)-xleft)*100.0
        if (xouts(j).ne.xouts(j+1)) then
         slope=(youts(j+1)-youts(j))/(xouts(j+1)-xouts(j))
         if (slope.ne.0.0) then
          bincp=youts(j)-slope*xouts(j)
!          xl=(zleft-bincp)/slope
!          sepin(iges)=(xl-xleft)*100.0
         endif
        endif
!        endif
       endif
       if (xouts(j).ge.radp) then
!        dxrr=(youts(j+1)-zright)*(youts(j)-zright)
!        if (dxrr.le.0.0) then
!        if (youts(j).eq.zright) sepout(iges)=(xright-xouts(j))*100.0
!        if (youts(j+1).eq.zright) sepout(iges)=(xright-xouts(j+1))*100.0
!        if (xouts(j).eq.xouts(j+1)) sepout(iges)=(xright-xouts(j))*100.0
        if (xouts(j).ne.xouts(j+1)) then
         slope=(youts(j+1)-youts(j))/(xouts(j+1)-xouts(j))
         if (slope.ne.0.0) then
          bincp=youts(j)-slope*xouts(j)
!          xr=(zright-bincp)/slope
!          sepout(iges)=(xright-xr)*100.0
         endif
        endif
!        endif
       endif
       if (youts(j).ge.zhp) then
!        dytt=(xouts(j)-xztop)*(xouts(j+1)-xztop)
!        if (dytt.le.0.0) then
!        if (xouts(j).eq.xouts(j+1)) septop(iges)=(ytop-youts(j))*100.0
        if (xouts(j).ne.xouts(j+1)) then
         slope=(youts(j+1)-youts(j))/(xouts(j+1)-xouts(j))
         bincp=youts(j)-slope*xouts(j)
!         yt=slope*xztop+bincp
!         septop(iges)=(ytop-yt)*100.0
        endif
!        endif
       endif
       if(youts(j).gt.zhp) cycle
!       dybb=(xouts(j)-xzbot)*(xouts(j+1)-xzbot)
!       if(dybb.gt.0.0) cycle
!       if(xouts(j).eq.xouts(j+1)) sepbot(iges)=(youts(j)-ybot)*100.0
       if(xouts(j).eq.xouts(j+1)) cycle
       slope=(youts(j+1)-youts(j))/(xouts(j+1)-xouts(j))
       bincp=youts(j)-slope*xouts(j)
!       yb=slope*xzbot+bincp
!       sepbot(iges)=(yb-ybot)*100.0
      enddo
!
 1085 continue
      call chisqr(iges)
      chifin=chisq(iges)
      nnn=1
      call betali(iges,rgrid,zgrid,nnn,kerror)
      if(kerror.gt.0) return
      peak(iges)=pres(1)/(.667_dp*wmhd(iges)/(volume(iges)/1.e6_dp))
      do i=2,nw
        if(rzzmax(i).gt.0.0) exit
      enddo
      amer(1,1)=2.*(rzzmax(1)-rzzmax(i))
      amer(1,2)=2.*(zzmax(1)-zzmax(i))
      amer(2,1)=2.*(rzzmax(1)-rzzmax(i))
      amer(2,2)=2.*(zzmax(1)-(zzmax(1)-zzmax(i)))
      bmer(1)=rzzmax(1)**2+zzmax(1)**2-rzzmax(i)**2-(zzmax(1)- &
                                                zzmax(i))**2
      bmer(2)=rzzmax(1)**2+zzmax(1)**2-rzzmax(i)**2-(zzmax(1)- &
                                                zzmax(i))**2
      n22=2
      x11=-1.0
      call decomp(n22,n22,amer,x11 ,imer,wmer)
      call solve(n22,n22,amer,bmer,imer)
      rmer=sqrt((bmer(1)-rzzmax(1))**2+(bmer(2)-zzmax(1))**2)
      i=i+1
      if (rzzmax(i).gt.0.0.and.rzzmax(i).lt.rzzmax(i-1)) then
        amer(1,1)=2.*(rzzmax(1)-rzzmax(i))
        amer(1,2)=2.*(zzmax(1)-zzmax(i))
        amer(2,1)=2.*(rzzmax(1)-rzzmax(i))
        amer(2,2)=2.*(zzmax(1)-(zzmax(1)-zzmax(i)))
        bmer(1)=rzzmax(1)**2+zzmax(1)**2-rzzmax(i)**2-(zzmax(1)- &
                                                  zzmax(i))**2
        bmer(2)=rzzmax(1)**2+zzmax(1)**2-rzzmax(i)**2-(zzmax(1)- &
                                                  zzmax(i))**2
        call decomp(n22,n22,amer,x11,imer,wmer)
        call solve(n22,n22,amer,bmer,imer)
!-----------------------------------------------------------------------
!--     need check for RMER 10/91 llao                                  --
!-----------------------------------------------------------------------
        rmer=(sqrt((bmer(1)-rzzmax(1))**2+(bmer(2)-zzmax(1))**2) &
              +rmer)/2.
      endif
      qmerci(iges)=2./(1.+elongm(iges)**2)-2.*(elongm(iges)-1.) &
                 *betap(iges)/elongm(iges)**2/(1.+elongm(iges)) &
                 +(elongm(iges)**2-1.)/(elongm(iges)**2+1.) &
                 *rmaxis/rmer
      if (qmerci(iges).gt.1.e-10_dp) then
        qmerci(iges)=sqrt(1./qmerci(iges))
      else
        qmerci(iges)=-99.
      endif
!
      do i=1,nw-1
        xsisii(i)=real(i-1,dp)/(nw-1)
      enddo
      xsisii(nw)=1.
!-----------------------------------------------------------------------
!--   write out S(shot).(time)_X files in flux space                  --
!-----------------------------------------------------------------------
      kwripre_s: if (kwripre.eq.2) then
        call setfnmd('s',ishot,itime,sfname)
        call setfnmd('o',ishot,itime,ofname)
        if (npress.gt.0) then
          sfname=sfname(1:13)//'_presd'
          open(unit=74,status='old',file=sfname,iostat=ier)
          if(ier.eq.0) close(unit=74,status='delete')
          open(unit=74,status='new',file=sfname)
          do i=1,npress
            xrpres=-rpress(i)
            write (74,*) xrpres,pressr(i),xdum,xdum
          enddo
          close(unit=74)
        endif
        if (nbeam.gt.0) then
          sfname=sfname(1:13)//'_pbeam'
          open(unit=74,status='old',file=sfname,iostat=ier)
          if(ier.eq.0) close(unit=74,status='delete')
          open(unit=74,status='new',file=sfname)
          do i=1,nbeam
            write (74,*) sibeam(i),pbeam(i),xdum,xdum
          enddo
          close(unit=74)
        endif
        sfname=sfname(1:13)//'_qpsi'
        open(unit=74,status='old',file=sfname,iostat=ier)
        if(ier.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=1,nw
          write (74,*) xsisii(i),qpsi(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_jor'
        open(unit=74,status='old',file=sfname,iostat=ier)
        if(ier.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=1,nw
          cjorka=cjor(i)/1000.
          write (74,*) xsisii(i),cjorka,xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_jorec'
        open(unit=74,status='old',file=sfname,iostat=ier)
        if(ier.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=1,nw
          cjorka=cjorec(i)/1000.
          write (74,*) xsisii(i),cjorka,xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_cjmse'
        open(unit=74,status='old',file=sfname,iostat=ier)
        if(ier.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        if (ishot.le.97400) then
          mcentral=15
        else
          mcentral=10
        endif
        do i=1,mcentral+1
          cjorka=cjmse(i)/1000.
          write (74,*) sigam(i),cjorka,xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_cjmsec'
        open(unit=74,status='old',file=sfname,iostat=ier)
        if(ier.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        if (ishot.le.97400) then
          mcentral=15
        else
          mcentral=10
        endif
        do i=1,mcentral+1
          cjorka=cjmsec(i)/1000.
          write (74,*) sigam(i),cjorka,xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_bzmse'
        open(unit=74,status='old',file=sfname,iostat=ier)
        if(ier.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=1,nstark
          cjorka=bzmse(i)
          write (74,*)  sigam(i),cjorka,xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_bzmsec'
        open(unit=74,status='old',file=sfname,iostat=ier)
        if(ier.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=1,nstark
          cjorka=bzmsec(i)
          write (74,*)  sigam(i),cjorka,xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_chigam'
        open(unit=74,status='old',file=sfname,iostat=ier)
        if(ier.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        sigams(1:nstark)=sigam(1:nstark)
        do j=1,nstark
          sigamnow=1.e9_dp
          do i=1,nstark
            if (sigams(i).lt.sigamnow.and.fwtgam(i).gt.0.0) then
              if (rrgam(iges,i).gt.0.0) then
                sigamnow=sigams(i)
                cjorka=chigam(i)
                inow=i
              endif
            endif
          enddo
          write (74,*)  sigamnow,cjorka,xdum,xdum
          sigams(inow)=1.e10_dp
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_presf'
        open(unit=74,status='old',file=sfname,iostat=ier)
        if(ier.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=1,nw
          write (74,*) xsisii(i),pres(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_prespf'
        open(unit=74,status='old',file=sfname,iostat=ier)
        if(ier.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=1,nw
          write (74,*) xsisii(i),pprime(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_prespfn'
        open(unit=74,status='old',file=sfname,iostat=ier)
        if(ier.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=1,nw
          ppnow=pprime(i)/pprime(1)
          write (74,*) xsisii(i),ppnow,xdum,xdum
        enddo
        close(unit=74)
!-----------------------------------------------------------------------
!--     write out MSE-LS files in normalized poloidal flux space      --
!-----------------------------------------------------------------------
        if (mmbmsels.gt.0.or.kdomsels.gt.0) then
          sfname=sfname(1:13)//'_cmls'
          open(unit=74,status='old',file=sfname,iostat=ier)
          if(ier.eq.0) close(unit=74,status='delete')
          open(unit=74,status='new',file=sfname)
          do i=1,nmsels
           if(sinmls(i).gt.0.0) &
             write (74,92924) sinmls(i),cmmls(iges,i),xdum,xdum
          enddo
          close(unit=74)
!
          ofname=ofname(1:13)//'_cmls'
          open(unit=74,status='old',file=ofname,iostat=ier)
          if(ier.eq.0) close(unit=74,status='delete')
          open(unit=74,status='new',file=ofname)
          xxtitle='Normalized Poloidal Flux'
          yytitle='Computed B MSE-LS (T)'
          zztitle='EFIT MSE-LS'
          write (74,93024) xxtitle
          write (74,93024) yytitle
          write (74,93024) zztitle
          do i=1,nmsels
            if(sinmls(i).gt.0.0) &
              write (74,92924) sinmls(i),cmmls(iges,i),xdum,xdum
          enddo
          close(unit=74)
!
          sfname=sfname(1:13)//'_cmls2'
          open(unit=74,status='old',file=sfname,iostat=ier)
          if(ier.eq.0) close(unit=74,status='delete')
          open(unit=74,status='new',file=sfname)
          do i=1,nmsels
           if(sinmls(i).gt.0.0) &
             write (74,92924) sinmls(i),cmmls2(iges,i),xdum,xdum
          enddo
          close(unit=74)
!
          sfname=sfname(1:13)//'_cmlsv'
          open(unit=74,status='old',file=sfname,iostat=ier)
          if(ier.eq.0) close(unit=74,status='delete')
          open(unit=74,status='new',file=sfname)
          do i=1,nmsels
           if(sinmls(i).gt.0.0) &
             write (74,92924) sinmls(i),cmmlsv(iges,i),xdum,xdum
          enddo
          close(unit=74)
!
          ofname=ofname(1:13)//'_cmlsv'
          open(unit=74,status='old',file=ofname,iostat=ier)
          if(ier.eq.0) close(unit=74,status='delete')
          open(unit=74,status='new',file=ofname)
          xxtitle='Normalized Poloidal Flux'
          yytitle='Computed Vacuum B MSE-LS (T)'
          zztitle='EFIT MSE-LS'
          write (74,93024) xxtitle
          write (74,93024) yytitle
          write (74,93024) zztitle
          do i=1,nmsels
           if(sinmls(i).gt.0.0) &
             write (74,92924) sinmls(i),cmmlsv(iges,i),xdum,xdum
          enddo
          close(unit=74)
!
          sfname=sfname(1:13)//'_emls'
          open(unit=74,status='old',file=sfname,iostat=ier)
          if(ier.eq.0) close(unit=74,status='delete')
          open(unit=74,status='new',file=sfname)
          do i=1,nmsels
           if (sinmls(i).gt.0.0) then
             ydum=sbmselt(iges,i)
             if(swtbmselt(iges,i).gt.1.e-06_dp) ydum=ydum/swtbmselt(iges,i)
             write (74,92924) sinmls(i),bmselt(iges,i),xdum,ydum
           endif
          enddo
          close(unit=74)
!
          ofname=ofname(1:13)//'_emls'
          open(unit=74,status='old',file=ofname,iostat=ier)
          if(ier.eq.0) close(unit=74,status='delete')
          open(unit=74,status='new',file=ofname)
          xxtitle='Normalized Poloidal Flux'
          yytitle='Measured B MSE-LS (T)'
          zztitle='EFIT MSE-LS'
          write (74,93024) xxtitle
          write (74,93024) yytitle
          write (74,93024) zztitle
          do i=1,nmsels
           if (sinmls(i).gt.0.0) then
             ydum=sbmselt(iges,i)
             if(swtbmselt(iges,i).gt.1.e-06_dp) ydum=ydum/swtbmselt(iges,i)
             write (74,92924) sinmls(i),bmselt(iges,i),xdum,ydum
           endif
          enddo
          close(unit=74)
        endif
!
        if (mmbmsels.gt.0) then
          open(unit=74,status='old',file=ofname,iostat=ier)
          if(ier.eq.0) close(unit=74,status='delete')
          open(unit=74,status='new',file=ofname)
          xxtitle='Normalized Poloidal Flux'
          yytitle='Chi Square'
          zzztitle='EFIT MSE-LS Chi2  = '
!
          ichisq=tchimls
          ichisq2=(tchimls-ichisq)*10.0
          write (jchisq2,94010) ichisq2
94010     format (i1.1)
          write (jchisq,94020) ichisq,jchisq2
94020     format (i6.6,'.',a1)
          zztitle=zzztitle//jchisq
!
          write (74,93024) xxtitle
          write (74,93024) yytitle
          write (74,93024) zztitle
          do i=1,nmsels
           if(sinmls(i).gt.0.0) &
             write (74,92924) sinmls(i),chimls(i),xdum,xdum
          enddo
          close(unit=74)
!
92924   format (4(1pe12.5,1x))
93024   format (a28)
        endif
!-----------------------------------------------------------------------
!--     write out MSE O files in normalized poloidal flux space       --
!-----------------------------------------------------------------------
        if (kstark.gt.0.or.kdomse.gt.0) then
          ofname=ofname(1:13)//'_cmse'
          open(unit=74,status='old',file=ofname,iostat=ier)
          if(ier.eq.0) close(unit=74,status='delete')
          open(unit=74,status='new',file=ofname)
          xxtitle='Normalized Poloidal Flux'
          yytitle='Computed MSE Pitch Angle'
          zztitle='EFIT MSE'
          write (74,93024) xxtitle
          write (74,93024) yytitle
          write (74,93024) zztitle
          do i=1,nstark
           if(swtgam(i).gt.1.e-6_dp) &
             write (74,92924) sigam(i),cmgam(i,iges),xdum,xdum
          enddo
          close(unit=74)
!
          ofname=ofname(1:13)//'_emse'
          open(unit=74,status='old',file=ofname,iostat=ier)
          if(ier.eq.0) close(unit=74,status='delete')
          open(unit=74,status='new',file=ofname)
          xxtitle='Normalized Poloidal Flux'
          yytitle='Measured MSE Pitch Angle'
          zztitle='EFIT MSE'
          write (74,93024) xxtitle
          write (74,93024) yytitle
          write (74,93024) zztitle
          do i=1,nstark
           if (swtgam(i).gt.1.e-6_dp) then
             ydum=siggam(iges,i)/swtgam(i)
             write (74,92924) sigam(i),tangam(iges,i),xdum,ydum
           endif
          enddo
          close(unit=74)
        endif
!
        if (kstark.gt.0) then
          ofname=ofname(1:13)//'_chi2mse'
          open(unit=74,status='old',file=ofname,iostat=ier)
          if(ier.eq.0) close(unit=74,status='delete')
          open(unit=74,status='new',file=ofname)
          xxtitle='Normalized Poloidal Flux'
          yytitle='Chi Square'
          zzztitle='EFIT MSE Chi2  = '
!
          ichisq=chigamt
          ichisq2=(chigamt-ichisq)*10.0
          write (jchisq2,94010) ichisq2
          write (jchisq,94020) ichisq,jchisq2
          zztitle=zzztitle//jchisq
!
          write (74,93024) xxtitle
          write (74,93024) yytitle
          write (74,93024) zztitle
          do i=1,nstark
           if(fwtgam(i).gt.0.0) &
             write (74,92924) sigam(i),chigam(i),xdum,xdum
          enddo
          close(unit=74)
        endif
      endif kwripre_s
!
!     if two q=1 surfaces, both psi values will be encoded in psiq1
!     one to left of decimal and one to the right
      psiq1=-1000.
!
!     determine if the slope of qpsi is negative at any grid points (idoqn=2)
      idoqn=1
      call zpline(nw,xsisii,qpsi,bfpol,cfpol,dfpol)
      do i=1,nw
         qpwant=speval(nw,xsisii(i),xsisii,qpsi,bfpol,cfpol,dfpol)
         if (qpwant.le.0.0) then
           idoqn=2
           exit
         endif
      enddo
!
#ifdef DEBUG_LEVEL2
      write (6,*) 'SHAPE q=1,2,3 q= ',idoqn,qpsi(1),qpsi(nw)
#endif
      nnn=1
      d11=30.
      d22=0.03_dp
      d33=0.01_dp
      nzz=0
      n22=2
      zxx=0.0
      do i=1,3
        pds(i)=100.
      enddo
      iend=3
      ! splining psi(q) is problematic if it is non-monotonic (has negative slope)
      if(idoqn.eq.1) call zpline(nw,qpsi,xsisii,bfpol,cfpol,dfpol)
      if(qpsi(nw).gt.1.0) iend=1
      if(qpsi(nw).gt.2.0) iend=2
      if(qpsi(nw).gt.3.0) iend=3
      jstart=2
      double=(idoqn.eq.2).and.(qpsi(1).gt.1.)
!---------------------------------------------------------------------
!--   Contour all surfaces with integer safety factors
!---------------------------------------------------------------------
      do i=1,iend
        ! determine the xsisii values of integer q surfaces (siwant)
        ! if q is non-monotonic, then use linear local interpolation
        qwant=i
        if(idoqn.eq.1 .and. qwant.lt.qpsi(1)+0.001_dp) cycle
        if (idoqn.eq.1 .and. i.ge.2) then
         siwant=seval(nw,qwant,qpsi,xsisii,bfpol,cfpol,dfpol)
        else
         do jjj=jstart,nw
          jj=jjj
          qppp=qwant-qpsi(jj)
          qmmm=qwant-qpsi(jj-1)
          if (qppp*qmmm.le.0.0) then
            siwant=xsisii(jj-1)+(xsisii(jj)-xsisii(jj-1))/ &
                   (qpsi(jj)-qpsi(jj-1))*qmmm
            if(jstart.eq.2) onedone=.true.
            jstart=jj+1
            go to 1105
          endif
         enddo
         cycle
        endif
 1105   continue
#ifdef DEBUG_LEVEL2
        write (6,*) 'i, iend, sw, qw = ',i, iend,siwant, qwant
#endif
        if(siwant.le.1.e-5_dp) cycle
        siwant=simag-siwant*(simag-psibry)
        kacerr=0
        call surfac(siwant,psi,nw,nh,rgrid,zgrid,xxtra(1,1),yxtra(1,1), &
                    nfind,npoint,drgrid,dzgrid,xmin,xmax,ymin,ymax,nzz, &
                    rmaxis,zmaxis,negcur,kacerr,2)
        if (nfind.le.40.and.icntour.eq.0) then
#ifdef DEBUG_LEVEL2
          write (6,*) ' SHAPE/SURFAC kacerr,i,nfind,qp,qm,si = ', &
                       kacerr,i,nfind,qppp,qmmm,siwant
#endif
          call cntour(rmaxis,zmaxis,siwant,rqmin,rqmax,ycmin,ycmax, &
                      yxcmin,yxcmax,xycmin,xycmax,d11,drgrid,d22, &
                      d33 ,d33 ,xmin,xmax,ymin,ymax,nzz,iautoc, &
                      xxtra(1,1),yxtra(1,1),nfind,rgrid,nw,zgrid,nh, &
                      c,n22,nh2,nttyo,npoint, &
                      negcur,bkx,lkx,bky,lky,kerror)
#ifdef DEBUG_LEVEL2
          write (6,*) ' SHAPE/CNTOUR kerror,i,nfind = ',kerror,i,nfind
#endif
          if(kerror.gt.0) return

        else
          rqmax=xxtra(1,1)
          rqmin=rqmax
          do k=2,nfind
            rqmax=max(xxtra(k,1),rqmax)
            rqmin=min(xxtra(k,1),rqmin)
          enddo
        endif
        pds(i)=50.*(rqmax-rqmin)
        if (i.eq.1) then
          if(.not.double) psiq1=siwant
!
!     This code will not compile on the HP or SuperCard.  The logic
!     needs to be fixed so that it will compile and then this code
!     can be put back in.
!
          if(double.and.(.not.onedone)) psiq1=psiq1+siwant  ! second value
        endif
      enddo
      aaq1(iges)=pds(1)
      aaq2(iges)=pds(2)
      aaq3(iges)=pds(3)
!---------------------------------------------------------------------
!--   Get 3/2 and 2/1 surface information        2002Jan25          --
!---------------------------------------------------------------------
#ifdef DEBUG_LEVEL2
      write (6,*) 'SHAPE q=3/2, 2/1'
#endif
      iend=2
      psin32(iges)=-99.0
      rq32in(iges)=-99.0
      psin21(iges)=-99.0
      rq21top(iges)=-99.0
      do i=1,iend
        if(i.eq.1) qwant=1.5_dp
        if(i.eq.2) qwant=2.0 ! TODO: already found in the previous loop??
        do j=1,nw-1
          jj=nw-j+1
          qppp=qwant-qpsi(jj)
          qmmm=qwant-qpsi(jj-1)
          if (qppp*qmmm.le.0.0) then
            siwant=xsisii(jj-1)+ (xsisii(jj)-xsisii(jj-1))/ &
              (qpsi(jj)-qpsi(jj-1))*qmmm
            psiwan=simag-siwant*(simag-psibry)
            kacerr=0
            call surfac(psiwan,psi,nw,nh,rgrid,zgrid,xxtra(1,1),yxtra(1,1), &
                        nfind,npoint,drgrid,dzgrid,xmin,xmax,ymin,ymax,nzz, &
                        rmaxis,zmaxis,negcur,kacerr,2)
            if (nfind.le.40.and.icntour.eq.0) then
#ifdef DEBUG_LEVEL2
              write (6,*) ' SHAPE/SURFAC kerror,i,nfind,qp,qm,si = ', &
                           kacerr,i,nfind,qppp,qmmm,psiwan
#endif
              call cntour(rmaxis,zmaxis,psiwan,rqmin,rqmax,ycmin,ycmax, &
                          yxcmin,yxcmax,xycmin,xycmax,d11,drgrid,d22, &
                          d33,d33,xmin,xmax,ymin,ymax,nzz,iautoc, &
                          xxtra(1,1),yxtra(1,1),nfind,rgrid,nw,zgrid,nh, &
                          c,n22,nh2,nttyo,npoint, &
                          negcur,bkx,lkx,bky,lky,kerror)
#ifdef DEBUG_LEVEL2
              write (6,*) ' SHAPE/CNTOUR kerror,i,nfind = ',kerror,i,nfind
#endif
              if(kerror.gt.0) return
            endif
            if (i.eq.1) then
              psin32(iges)=siwant
              rqmin=xxtra(1,1)
              do k=2,nfind
                rqmin=min(xxtra(k,1),rqmin)
              enddo
              rq32in(iges)=100.0*rqmin
            endif
            if (i.eq.2) then
              psin21(iges)=siwant
              zqmax=yxtra(1,1)
              do k=2,nfind
                if (zqmax.lt.yxtra(k,1)) then
                  rqmax=xxtra(k,1)
                  zqmax=yxtra(k,1)
                endif
              enddo
              rq21top(iges)=100.0*rqmax
            endif
            exit
          endif
        enddo
      enddo
!---------------------------------------------------------------------
!--   Get information from the psi=0.95 and 0.01 surfaces (no contour)
!---------------------------------------------------------------------
      if(idoqn.eq.1) call zpline(nw,xsisii,qpsi,bfpol,cfpol,dfpol)
      siwant=0.95_dp
      q95(iges)=seval(nw,siwant,xsisii,qpsi,bfpol,cfpol,dfpol)
      shearb(iges)=speval(nw,siwant,xsisii,qpsi,bfpol,cfpol,dfpol)
      shearb(iges)=shearb(iges)/q95(iges)
      siwant=0.01_dp
      qpsic=seval(nw,siwant,xsisii,qpsi,bfpol,cfpol,dfpol)
      shearc=speval(nw,siwant,xsisii,qpsi,bfpol,cfpol,dfpol)
      shearc=shearc/qpsic
!----------------------------------------------------------------------
!--   if not diverted, get q at PSIWANT by interpolation             --
!----------------------------------------------------------------------
      if ((abs(psiwant-1.0).le.1.e-05_dp).or.(abs(dismin).le.0.100_dp) &
        .or.(psiwant.le.1.e-5_dp)) then
        siwwww=psiwant
        if (siwwww.le.1.e-5_dp) siwwww=1.e-5_dp
        qsiwant(iges)=seval(nw,siwwww,xsisii,qpsi,bfpol,cfpol,dfpol)
      else
       siwant=simag+psiwant*(psibry-simag)
       call surfac(siwant,psi,nw,nh,rgrid,zgrid,bfpol,dfpol,nfounc, &
                   npoint,drgrid,dzgrid,xmin,xmax,ymin,ymax,nnn, &
                   rmaxis,zmaxis,negcur,kerror,1)
       if(kerror.gt.0) return
       do k=1,nfounc
        cfpol(k)=1./bfpol(k)**2
       enddo
       call fluxav(cfpol,bfpol,dfpol,nfounc,psi,rgrid,nw,zgrid,nh, &
                   r2sdry(i),nzz,sdlobp,sdlbp)
       do k=1,nfounc
        cfpol(k)=1./bfpol(k)
       enddo
       call fluxav(cfpol,bfpol,dfpol,nfounc,psi,rgrid,nw,zgrid,nh, &
                   r1sdry(i),nzz,sdlobp,sdlbp)
       rotation_1: if (kvtor.gt.0) then
         do k=1,nfounc
           cfpol(k)=bfpol(k)**2
         enddo
         call fluxav(cfpol,bfpol,dfpol,nfounc,psi,rgrid,nw,zgrid, &
                     nh,r2wdry,nzz,sdlobp,sdlbp)
         do k=1,nfounc
           cfpol(k)=((bfpol(k)/rvtor)**2-1.)**2
         enddo
         call fluxav(cfpol,bfpol,dfpol,nfounc,psi,rgrid,nw,zgrid, &
                     nh,r4wdry,nzz,sdlobp,sdlbp)
         if (kvtor.eq.3) then
           prew0=pwcurr(psiwant,kwwcur)
           pres0=prcurr(psiwant,kppcur)
           if (abs(pres0).gt.1.e-10_dp) then
             pwop0=prew0/pres0
           else
             pwop0=0.0
           endif
           do k=1,nfounc
             cfpol(k)=((bfpol(k)/rvtor)**2-1.)
             cfpol(k)=cfpol(k)*exp(pwop0*cfpol(k))
           enddo
           call fluxav(cfpol,bfpol,dfpol,nfounc,psi,rgrid,nw,zgrid, &
                       nh,rp2wdry,nzz,sdlobp,sdlbp)
           do k=1,nfounc
             cfpol(k)=((bfpol(k)/rvtor)**2-1.)
             cfpol(k)=exp(pwop0*cfpol(k))
           enddo
           call fluxav(cfpol,bfpol,dfpol,nfounc,psi,rgrid,nw,zgrid, &
                       nh,rpwdry,nzz ,sdlobp,sdlbp)
         endif
       endif rotation_1
       r2surs=r2sdry(i)*sdlobp
       fpnow=ffcurr(psiwant,kffcur)
       fpnow=fpnow*tmu
       qsiwant(iges)=abs(fpnow)/twopi*r2surs
      endif
!
      call zpline(nw,xsisii,cjor,bfpol,cfpol,dfpol)
      siwant=0.95_dp
      cjor95(iges)=seval(nw,siwant,xsisii,cjor,bfpol,cfpol,dfpol)
!-----------------------------------------------------------------------
!--   shear=dq/d(sqrtVn)/q and J at psiwant                           --
!--   normalize all J's to I/area                         96/04/11    --
!-----------------------------------------------------------------------
      cjorsw(iges)=seval(nw,psiwant,xsisii,cjor,bfpol,cfpol,dfpol)
      if (abs(ipmhd(iges)).gt.1.e-1_dp) then
        cjor0(iges)=cjor(1)*area(iges)/ipmhd(iges)/1.0e4_dp
        cjor95(iges)=cjor95(iges)*area(iges)/ipmhd(iges)/1.0e4_dp
        cjorsw(iges)=cjorsw(iges)*area(iges)/ipmhd(iges)/1.0e4_dp
        siwant=0.995
        cjor99(iges)=seval(nw,siwant,xsisii,cjor,bfpol,cfpol,dfpol)
        cjor99(iges)=cjor99(iges)*area(iges)/ipmhd(iges)/1.0e4_dp
      endif
      do i=1,nw
        bpres(i)=sqrt(volp(i)/volp(nw))
      enddo
      call zpline(nw,xsisii,bpres,bfpol,cfpol,dfpol)
      siwant=0.95_dp
      xsi95=seval(nw,siwant,xsisii,bpres,bfpol,cfpol,dfpol)
      siwant=0.01_dp
      xsi01=seval(nw,siwant,xsisii,bpres,bfpol,cfpol,dfpol)
      xsiwww=psiwant
      if(xsiwww.le.1.e-5_dp) xsiwww=1.e-5_dp
      xsiww=seval(nw,xsiwww,xsisii,bpres,bfpol,cfpol,dfpol)
      call zpline(nw,bpres,qpsi,bfpol,cfpol,dfpol)
      ssiwant(iges)=speval(nw,xsiww,bpres,qpsi,bfpol,cfpol,dfpol) &
                    /qsiwant(iges)
      ssi95(iges)=speval(nw,xsi95,bpres,qpsi,bfpol,cfpol,dfpol) &
                    /q95(iges)
      ssi01=speval(nw,xsi01,bpres,qpsi,bfpol,cfpol,dfpol) &
                    /qpsic
!-------------------------------------------------------------------
!--    get average current density in outer at 95% flux           --
!-------------------------------------------------------------------
       siavej=0.95_dp
       siwant=siavej
       siwant=simag+siwant*(psibry-simag)
       call surfac(siwant,psi,nw,nh,rgrid,zgrid,bfpol,dfpol,nfounc, &
                   npoint,drgrid,dzgrid,xmin,xmax,ymin,ymax,nnn, &
                   rmaxis,zmaxis,negcur,kerror,1)
       if(kerror.gt.0) return
       call fluxav(bfpol,bfpol,dfpol,nfounc,psi,rgrid,nw,zgrid,nh, &
                   rxxrry,nzz ,sdlobp,sdlbp)
      aream=0.0
      xyma=dfpol(1)
      do i=2,nfounc
        xypa=dfpol(i)
        dx=bfpol(i)-bfpol(i-1)
        aream=aream+(xyma+xypa)/2.0*dx
        xyma=xypa
      enddo
      aream=abs(aream)*1.e4_dp
      delarea=area(iges)-aream
      pasnow=sdlbp/tmu/twopi
      cj1now=(ipmhd(iges)-pasnow)/delarea
      cj1ave(iges)=cj1now/ipmhd(iges)*area(iges)
      cj1now=cj1now*10.
!     if (iand(iout,1).ne.0) then
!     write (nout,*) ipmhd(iges),pasnow,area(iges),aream,cj1ave(iges) &
!                     ,cj1now,xsi95
!     endif
!
      call zpline(nw,xsisii,pprime,bfpol,cfpol,dfpol)
      siwant=0.95_dp
      pp95(iges)=seval(nw,siwant,xsisii,pprime,bfpol,cfpol,dfpol)
      pp95(iges)=pp95(iges)/pres(1)*sidif
!-----------------------------------------------------------------------
!--   switch order of separatrix, Z of first < 0                      --
!-----------------------------------------------------------------------
      if (zseps(1,iges).gt.0.0) then
        !call errctrl_msg('shapesurf','Switching separatrix 1 and 2',3)
        rtemp=rseps(1,iges)
        rseps(1,iges)=rseps(2,iges)
        rseps(2,iges)=rtemp
        rtemp=zseps(1,iges)
        zseps(1,iges)=zseps(2,iges)
        zseps(2,iges)=rtemp
      endif
      aspect=rout(iges)/aminor(iges)
      qstar(iges)=rcentr*abs(bcentr(iges))/tmu/abs(ipmhd(iges))* &
                 (elong(iges)**2+1.)/2./aspect**2
      olamda=betap(iges)+li(iges)/2.
      olamda=1.+0.5_dp*olamda**2
      atri=(utri(iges)+ltri(iges))/2.

      fkdel=1.24_dp-0.54_dp*elong(iges)+0.3_dp*(elong(iges)**2+atri**2) &
           +0.13_dp*atri
      qstar(iges)=qstar(iges)*fkdel*(1.+olamda/aspect**2)
      sepexp(iges)=-100.
      if (kdata.ne.4) then
        pohm=vloopt(iges)*ipmhd(iges)
        pohm=max(czero,pohm)
        pin=pohm+pbinj(iges)
        if (pin.gt.1.e-06_dp) then
          taumhd(iges)=wmhd(iges)/pin*1000.
          taudia(iges)=wdia(iges)/pin*1000.
        endif
      endif
!
      call lenco2(xout,yout,nfound,iges)
!----------------------------------------------------------------------
!--   compute fast and thermal stored energy                         --
!--   using Heibrink and Duong's approximate formula   L.Lao 10/91   --
!----------------------------------------------------------------------
      wtherm(iges)=wmhd(iges)
      wfbeam(iges)=0.0
      taujd3(iges)=0.0
      tauthn(iges)=0.0
      if ((dco2v(iges,2).ge.1.e+10_dp) .and. (pbinj(iges).gt.100.) &
          .and. (abs(ipmhd(iges)).gt.10000.) .and. (wmhd(iges).gt.1000.)) then
        cons1=dco2v(iges,2)*1.e-13_dp                      !normalized density
        cons2=wmhd(iges)/cons1/volume(iges)/3.36_dp/1.e-6_dp  !electron temperature (ev)
        cons3=75000.                                    !w_beam (ev)
        cons4=22.*cons2                                  !w_crit (ev)
        cons5=log(1.+(75000./cons4)**1.5_dp)
        cons6=8.37e-2_dp*(cons2/1000)**1.5_dp/cons1*cons5  !thermalization tau_s (s)
        cons7=sqrt(7.5e4_dp)                               !normalized v_beam
        cons8=sqrt(cons4)                               !normalized v_crit
        cons9=cons7**2-cons7*cons8+cons8**2
        cons10=log((cons7+cons8)**2/cons9)
        cons11=atan((2*cons7-cons8)/1.732_dp/cons8)
        cons12=2._dp/3._dp*(cons8/cons7)**2
        cons13=1.+cons12*(0.5_dp*cons10-1.732_dp*(pi/6+cons11))
        cons14=cons6*3/cons5                            !spitzer tau_se (sec)
        cons15=2.75e-1_dp*pbinj(iges)*cons14*cons13          !stored beam energy (j)
        cons16=wmhd(iges)-cons15                         !thermal energy (j)
        cons17=cons16/wmhd(iges)
        wtherm(iges)=cons16
        wfbeam(iges)=cons15
!-----------------------------------------------------------------------
!--     get thermal and jet-diiid confinement times                   --
!-----------------------------------------------------------------------
        cons21=pohm/1.e6_dp                              !ohmic power in mw
        cons22=taumhd(iges)*cons17
        ptotal = pin/1.e6_dp                              !toal power in mw
        if (ptotal.gt.0.01_dp) then
          consa1=ptotal**0.49_dp
          consa2=(abs(ipmhd(iges))/1.e6_dp)**1.08_dp
          consa3=(rout(iges)/100.)**1.45_dp
          taujd3(iges)=110.*consa2*consa3/consa1
        else
          taujd3(iges)=0.0
        endif
        if  (taujd3(iges).gt.0.0001_dp) then
          tauthn(iges)=cons22/taujd3(iges)
        else
          tauthn(iges)=0.0
        endif
      endif
!----------------------------------------------------------------------------
!--   find the vessel intersecting points if diverted plasmas              --
!----------------------------------------------------------------------------
      dolubaf(iges)=-89.0
      diludom(iges)=-89.0
      dolubafm(iges)=-89.0
      diludomm(iges)=-89.0
      dminux(iges)=-89.0
      dminlx(iges)=-89.0
      ratsol(iges)=-89.0
      rvsiu(iges)=-89.0
      zvsiu(iges)=-89.0
      rvsid(iges)=-89.0
      zvsid(iges)=-89.0
      rvsou(iges)=-89.0
      zvsou(iges)=-89.0
      rvsod(iges)=-89.0
      zvsod(iges)=-89.0
      if ((dpsi.gt.psitol) .and. (limloc(iges).ne.'IN')) then
!----------------------------------------------------------------------
!--   get 1 cm SOL width ratio                                       --
!----------------------------------------------------------------------
      rnow=xmax+0.01_dp
      znow=zxmax
      call seva2d(bkx,lkx,bky,lky,c,rnow,znow,pds,ier,n111)
      psinow=pds(1)
      nmax=(xmin-rgrid(1))/(rgrid(2)-rgrid(1))+1
      znow=zxmin
      rrout=xmin
      psiout=psibry
      dsiout=psinow-psiout
      do i=1,nmax
       k= nmax-i +1
       rrin=rgrid(k)
       call seva2d(bkx,lkx,bky,lky,c,rrin,znow,pds,ier,n111)
       psiin =pds(1)
       dsiin =psinow-psiin
       if (dsiin*dsiout.le.0.0) then
         do n=1,20
          rmm=0.5_dp*(rrin+rrout)
          call seva2d(bkx,lkx,bky,lky,c,rmm ,znow,pds,ier,n111)
          psimm=pds(1)
          dsimm=psinow-psimm
          if (dsimm*dsiin.le.0.0) then
            rrout=rmm
          else
            rrin =rmm
          endif
         enddo
         exit
       endif
      enddo
      ratsol(iges)=100.*(xmin-rrin)
!----------------------------------------------------------------------
!--   get distance from x point to limiter                           --
!----------------------------------------------------------------------
      ycut=0.5_dp*ymax
      if (zfsep.gt.ycut) then
        dminux(iges)=1.e10_dp
        do i=1,limitr
          if (ylim(i).gt.ycut) then
            dminow=(rfsep-xlim(i))**2+(zfsep-ylim(i))**2
            dminux(iges)=min(dminux(iges),dminow)
          endif
        enddo
        dminux(iges)=100.0*sqrt(dminux(iges))
      else
        dminlx(iges)=1.e10_dp
        ycut=0.5_dp*ymin
        do i=1,limitr
          if (ylim(i).lt.ycut) then
            dminow=(rfsep-xlim(i))**2+(zfsep-ylim(i))**2
            dminlx(iges)=min(dminlx(iges),dminow)
          endif
        enddo
        dminlx(iges)=100.0*sqrt(dminlx(iges))
      endif
!---------------------------------------------------------------
!--   Get distance from second X point to limiter             --
!---------------------------------------------------------------
      ycut=0.5_dp*ymax
      if (sissep.gt.-1.0e10_dp) then
        if (zssep.gt.ycut) then
          dminux(iges)=1.e10_dp
          do i=1,limitr
            if (ylim(i).gt.ycut) then
              dminow=(rssep-xlim(i))**2+(zssep-ylim(i))**2
              dminux(iges)=min(dminux(iges),dminow)
            endif
          enddo
          dminux(iges)=100.0*sqrt(dminux(iges))
        else
          dminlx(iges)=1.e10_dp
          ycut=0.5_dp*ymin
          do i=1,limitr
            if (ylim(i).lt.ycut) then
              dminow=(rssep-xlim(i))**2+(zssep-ylim(i))**2
              dminlx(iges)=min(dminlx(iges),dminow)
            endif
          enddo
          dminlx(iges)=100.0*sqrt(dminlx(iges))
        endif
      endif
!
      rkkk=xmax
      zkkk=zxmax
      call seva2d(bkx,lkx,bky,lky,c,rkkk,zkkk,pds,ier,n111)
      ssimax=pds(1)
      rkkk=rmaxfs
      zkkk=zmaxfs
      call seva2d(bkx,lkx,bky,lky,c,rkkk,zkkk,pds,ier,n111)
      ssimfs=pds(1)
      ixyz=-2
73723 dmaxfs=0.00001_dp
      dmaxfs0=dmaxfs
      nmaxfs=1
77723 xxtras=rmaxfs+dmaxfs
      yxtras=zmaxfs
      rkkk=xxtras
      zkkk=yxtras
      call seva2d(bkx,lkx,bky,lky,c,rkkk,zkkk,pds,ier,n111)
      ssitra=pds(1)
      ixl=1
      if (ipmeas(iges).lt.-1.e3_dp) then
         nerr=10000
      else
         nerr=0
      endif
      call bound(psi,nw,nh,nwnh,psiout,xmin,xmax,ymin,ymax,zeros, &
                 rgrid,zgrid,xxtras,yxtras,ixyz,limtrs,xlims,ylims, &
                 xxtra(1,ixl),yxtra(1,ixl),npxtra(ixl),xlims(1),nxtrap, &
                 rymin,rymax,dpsi,zxmin,zxmax,nerr,ishot,itime, &
                 limfag,radum,kbound,tolbndpsi)
      if (nerr.gt.0) then
        kerror = 1
        return
      end if
      dis2p=(xxtra(1,ixl)-xxtra(npxtra(ixl),ixl))**2
      dis2p=sqrt(dis2p+(yxtra(1,ixl)-yxtra(npxtra(ixl),ixl))**2)
      if (dis2p.lt.0.1_dp*drgrid) then
        nmaxfs=nmaxfs+1
        dmaxfs=dmaxfs0*nmaxfs
        if (nmaxfs.le.20) go to 77723
      endif
      rsepnose=-999.
      if (dis2p.ge.0.1_dp*drgrid) then
        do i=1,npxtra(ixl)-1
          if ((yxtra(i,ixl)-znose)*(yxtra(i+1,ixl)-znose).le.0.0) then
            rsepnose=xxtra(i,ixl)+(xxtra(i+1,ixl)-xxtra(i,ixl))/ &
              (yxtra(i+1,ixl)-yxtra(i,ixl))*(znose-yxtra(i,ixl))
            exit
          endif
        enddo
        zerovs(1)=1.0
        do i=1,npxtra(ixl)
          zerold=zerovs(1)
          call zlim(zerovs,nnn,nnn,limitr,xlim,ylim,xxtra(i,ixl), &
                    yxtra(i,ixl),limfag)
          if (zerold.gt.0.01_dp.and.zerovs(1).lt.0.01_dp) exit
        enddo
        if (zerold.gt.0.01_dp.and.zerovs(1).lt.0.01_dp) then
          rinvs=xxtra(i-1,ixl)
          zinvs=yxtra(i-1,ixl)
          routvs=xxtra(i,ixl)
          zoutvs=yxtra(i,ixl)
          do i=1,20
            ravs(1)=0.5_dp*(rinvs+routvs)
            zavs(1)=0.5_dp*(zinvs+zoutvs)
            call zlim(zerovs,nnn,nnn,limitr,xlim,ylim,ravs,zavs,limfag)
            if (zerovs(1).gt.0.01_dp) then
              rinvs=ravs(1)
              zinvs=zavs(1)
            else
              routvs=ravs(1)
              zoutvs=zavs(1)
            endif
          enddo
          rvsout(iges)=ravs(1)*100.
          zvsout(iges)=zavs(1)*100.
        else
          rvsout(iges)=0.
          zvsout(iges)=0.
        endif
      endif ! dis2p.ge.0.1_dp*drgrid
      if ((zavs(1)*zfsep.le.0.0).and.(ixyz.eq.-2)) then
        ixyz=-1
        go to 73723
      endif
      if (dis2p.ge.0.1_dp*drgrid) then
        call seva2d(bkx,lkx,bky,lky,c,ravs(1),zavs(1),pds,ier,n111)
        ssinow=pds(1)
!------------------------------------------------------------------
!--     get distance between inner leg and upper dome
!------------------------------------------------------------------
        ycut=ymax*0.5_dp
        ycutm=ymin*0.5_dp
        if (zavs(1).gt.ycut.and.ravs(1).lt.rymax) then
          rvsiu(iges)=rvsout(iges)
          zvsiu(iges)=zvsout(iges)
          if (ishot.ge.100771) then
            diludom(iges)=1.e6_dp
            do i=1,npxtra(ixl)
              if (yxtra(i,ixl).ge.ycut) then
                dilnow=(rudom - xxtra(i,ixl))**2 + &
                            (zudom - yxtra(i,ixl))**2
                if (dilnow.lt.diludom(iges)) then
                  diludom(iges)=dilnow
                  zilnow=yxtra(i,ixl)
                endif
              endif
            enddo
            diludom(iges)=100.0*sqrt(diludom(iges))
            if(zilnow.lt.zudom) diludom(iges)=-diludom(iges)
          endif
        elseif (zavs(1).gt.ycut.and.ravs(1).gt.rymax) then
!------------------------------------------------------------------
!--       get distance between outer leg and upper baffle
!------------------------------------------------------------------
          rvsou(iges)=rvsout(iges)
          zvsou(iges)=zvsout(iges)
          if (ishot.ge.100771) then
            dolubaf(iges)=1.e6_dp
            do i=1,npxtra(ixl)
              if (yxtra(i,ixl).ge.ycut) then
                dolnow=(rubaf - xxtra(i,ixl))**2 + &
                            (zubaf - yxtra(i,ixl))**2
                if (dolnow.lt.dolubaf(iges)) then
                  dolubaf(iges)=dolnow
                  rolnow=xxtra(i,ixl)
                endif
              endif
            enddo
            dolubaf(iges)=100.0*sqrt(dolubaf(iges))
            if(rolnow.gt.rubaf) dolubaf(iges)=-dolubaf(iges)
          endif
        elseif (zavs(1).lt.ycutm.and.ravs(1).gt.rymin) then
!------------------------------------------------------------------
!--       get lower strike points
!------------------------------------------------------------------
          rvsod(iges)=rvsout(iges)
          zvsod(iges)=zvsout(iges)
        elseif (zavs(1).lt.ycutm.and.ravs(1).lt.rymin) then
          rvsid(iges)=rvsout(iges)
          zvsid(iges)=zvsout(iges)
        endif
      endif ! dis2p.ge.0.1_dp*drgrid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ring calculation here!!!
      if ((lring.ne.0).and.(zvsout(iges).ne.0.)) then
        iring=1
        rxp=rseps(1,iges)/100.
        zxp=zseps(1,iges)/100.
        zrmin=floorz
        byringr(iring)=rxp  ! make first point the xpoint
        byringz(iring)=zxp
        do i=1,npxtra(ixl)
          zrmin=min(zrmin,yxtra(i,ixl))
        enddo
        do i=1,npxtra(ixl)
          byring=((xxtra(i,ixl).gt.rxp).and. &
                  (yxtra(i,ixl).lt.zxp).and. &
                  (yxtra(i,ixl).ge.zrmin))
          if (byring) then ! byring >>> points are close to the ring
            iring=iring+1
            byringr(iring)=xxtra(i,ixl)       ! save points in byring arrays
            byringz(iring)=yxtra(i,ixl)
          endif
        enddo
        if (iring.ge.3) then
          call order(byringr,byringz,iring) ! put point in increasing order in z
          call zpline(iring,byringz,byringr,bfpol,cfpol,dfpol)  ! spline r(z)
          zringmin=zrmin
          zringmax=max(ringz(1)+.005,zxp)
          dzring=(zringmax-zringmin)/nh2
          do i=1,nh2
            zval=zringmin+(i-1)*dzring
            yxtra(i,ixl)=zval
            rval=seval(iring,zval,byringz,byringr,bfpol,cfpol,dfpol)
            xxtra(i,ixl)=rval
          enddo
          do i=1,nh2
            byringr(i)=xxtra(i,ixl)  !overwrite byring arrays with splined values
            byringz(i)=yxtra(i,ixl)
          enddo
          xrmin=1e20
          xrmax=0
          rbar=0
          zbar=0
          do i=1,6
            xrmax=max(xrmax,ringr(i))
            xrmin=min(xrmin,ringr(i))
            rbar=rbar+ringr(i)/6.
            zbar=zbar+ringz(i)/6.
          enddo
          dismin=1e20
          dsmin=1e20
          do i=1,6
            dsmin=min(dsmin,sqrt((rbar-ringr(i))**2+(zbar-ringz(i))**2))
          enddo
          call dslant(byringr,byringz,nh2,rxp,xrmax,floorz,zxp, &
                      ringr(4),ringz(4),ringr(5),ringz(5),gap1) ! vertical ring face
          call dslant(byringr,byringz,nh2,rxp,xrmax,floorz,zxp, &
                      ringr(5),ringz(5),ringr(6),ringz(6),gap2) ! slanted ring face
          ringap=gap2
          if(gap1.lt.gap2) ringap=gap1 ! closest approach, separatrix to ring
          do i=1,nh2
            dismin=min(dismin, &
                       sqrt((rbar-byringr(i))**2+(zbar-byringz(i))**2))
          enddo
          if (ringap.eq.gap2) then
            if(dismin.lt.dsmin) gap2=-gap2 ! separatrix intersects ring
          endif
          if (ringap.eq.gap1) then ! check for intersection more carefully
            do i=1,nh2
              if(dismin.eq.sqrt((rbar-byringr(i))**2+(zbar-byringz(i))**2)) &
                iring=i
            enddo
            do i=1,6
              if ((ringz(i)-byringz(iring).ne.0.).and.(ringap.eq.gap1)) then
                if(atan2(ringr(i)-byringr(iring),ringr(i)-byringz(iring)).lt. &
                    0.) gap1=-gap1   ! separatrix intersects ring
              endif
            enddo
          endif
          if(abs(gap1).gt.1.e-6_dp.and.abs(gap1).lt.1.e6_dp) &
            oring(iges)=gap2*gap1/abs(gap1) ! closest approach to slanted face
          iring=0 !start cleaning up byrung arrays for plots in subr expand
        endif
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ring calculation end!!!!
      ixyz=-2
26199 dminfs=0.0001_dp
      dminfs0=dminfs
      nminfs=1
16199 xxtras=rminfs-dminfs
      yxtras=zminfs
      ixl=1
      if (ipmeas(iges).lt.-1.e3_dp) then
        nerr=10000
      else
        nerr=0
      endif
      call bound(psi,nw,nh,nwnh,psiout,xmin,xmax,ymin,ymax,zeros, &
                 rgrid,zgrid,xxtras,yxtras,ixyz,limtrs,xlims,ylims, &
                 xxtra(1,ixl),yxtra(1,ixl),npxtra(ixl),xlims(1),nxtrap, &
                 rymin,rymax,dpsi,zxmin,zxmax,nerr,ishot,itime, &
                 limfag,radum,kbound,tolbndpsi)
      if (nerr.gt.0) then
        kerror = 1
        return
      end if
      dis2p=(xxtra(1,ixl)-xxtra(npxtra(ixl),ixl))**2
      dis2p=sqrt(dis2p+(yxtra(1,ixl)-yxtra(npxtra(ixl),ixl))**2)
      if (dis2p.lt.0.1_dp*drgrid)then
        nminfs=nminfs+1
        dminfs=dminfs0*nminfs
        if (nminfs.le.20) go to 16199
      endif
      if (dis2p.ge.0.1_dp*drgrid) then
        zerovs(1)=1.0
        do i=1,npxtra(ixl)
          zerold=zerovs(1)
          call zlim(zerovs,nnn,nnn,limitr,xlim,ylim,xxtra(i,ixl), &
                    yxtra(i,ixl),limfag)
          if(zerold.gt.0.01_dp.and.zerovs(1).lt.0.01_dp) exit
        enddo
        if (zerold.gt.0.01_dp.and.zerovs(1).lt.0.01_dp) then
          rinvs=xxtra(i-1,ixl)
          zinvs=yxtra(i-1,ixl)
          routvs=xxtra(i,ixl)
          zoutvs=yxtra(i,ixl)
          do i=1,20
            ravs(1)=0.5_dp*(rinvs+routvs)
            zavs(1)=0.5_dp*(zinvs+zoutvs)
            call zlim(zerovs,nnn,nnn,limitr,xlim,ylim,ravs,zavs,limfag)
            if (zerovs(1).gt.0.01_dp) then
              rinvs=ravs(1)
              zinvs=zavs(1)
            else
              routvs=ravs(1)
              zoutvs=zavs(1)
            endif
          enddo
          rvsin(iges)=ravs(1)*100.
          zvsin(iges)=zavs(1)*100.
        else
          rvsin(iges)=0.
          zvsin(iges)=0.
        endif
      endif
      if ((zavs(1)*zfsep.le.0.0).and.(ixyz.eq.-2)) then
        ixyz=-1
        go to 26199
      endif
      if (dis2p.ge.0.1_dp*drgrid) then
        call seva2d(bkx,lkx,bky,lky,c,ravs(1),zavs(1),pds,ier,n111)
        ssinow=pds(1)
!------------------------------------------------------------------
!--     get distance between inner leg and upper dome            --
!------------------------------------------------------------------
        ycut=ymax*0.5_dp
        ycutm=ymin*0.5_dp
        if (zavs(1).gt.ycut.and.ravs(1).lt.rymax) then
          rvsiu(iges)=rvsin(iges)
          zvsiu(iges)=zvsin(iges)
          if (ishot.ge.100771) then
            diludom(iges)=1.e6_dp
            do i=1,npxtra(ixl)
              if (yxtra(i,ixl).ge.ycut) then
                dilnow=(rudom - xxtra(i,ixl))**2 + &
                          (zudom - yxtra(i,ixl))**2
                if (dilnow.lt.diludom(iges)) then
                  diludom(iges)=dilnow
                  zilnow=yxtra(i,ixl)
                endif
              endif
            enddo
            diludom(iges)=100.0*sqrt(diludom(iges))
            if(zilnow.lt.zudom) diludom(iges)=-diludom(iges)
          endif
        elseif (zavs(1).gt.ycut.and.ravs(1).gt.rymax) then
!------------------------------------------------------------------
!--       get distance between outer leg and upper baffle        --
!------------------------------------------------------------------
          rvsou(iges)=rvsin(iges)
          zvsou(iges)=zvsin(iges)
          if (ishot.ge.100771) then
            dolubaf(iges)=1.e6_dp
            do i=1,npxtra(ixl)
              if (yxtra(i,ixl).ge.ycut) then
                dolnow=(rubaf - xxtra(i,ixl))**2 + &
                          (zubaf - yxtra(i,ixl))**2
                if (dolnow.lt.dolubaf(iges)) then
                  dolubaf(iges)=dolnow
                  rolnow=xxtra(i,ixl)
                endif
              endif
            enddo
            dolubaf(iges)=100.0*sqrt(dolubaf(iges))
            if(rolnow.gt.rubaf) dolubaf(iges)=-dolubaf(iges)
          endif
        elseif (zavs(1).lt.ycutm.and.ravs(1).gt.rymin) then
!------------------------------------------------------------------
!--       get lower strike points                                --
!------------------------------------------------------------------
          rvsod(iges)=rvsin(iges)
          zvsod(iges)=zvsin(iges)
        elseif (zavs(1).lt.ycutm.and.ravs(1).lt.rymin) then
          rvsid(iges)=rvsin(iges)
          zvsid(iges)=zvsin(iges)
        endif
      endif
      if (sissep.gt.-1.e10_dp) then
!------------------------------------------------------------------
!--   second separatrix distances from outboard side             --
!------------------------------------------------------------------
      rnow=rssep
      znow=zssep
      psinow=sissep
      nmin=(xmax-rgrid(1))/(rgrid(2)-rgrid(1))+2
      iznow=(zxmax-zgrid(1)+1.0e-6_dp)/(zgrid(2)-zgrid(1))+1
      znow=zgrid(iznow)
      rrin=xmax
      psiin=psibry
      dsiin =psinow-psiin
      do k=nmin,nw
        rrout=rgrid(k)
        call seva2d(bkx,lkx,bky,lky,c,rrout,znow,pds,ier,n111)
        psiout=pds(1)
        dsiout=psinow-psiout
        if (dsiin*dsiout.le.0.0) then
          do n=1,20
            rmm=0.5_dp*(rrin+rrout)
            call seva2d(bkx,lkx,bky,lky,c,rmm ,znow,pds,ier,n111)
            psimm=pds(1)
            dsimm=psinow-psimm
            if (dsimm*dsiin.le.0.0) then
              rrout=rmm
            else
              rrin=rmm
            endif
          enddo
          exit
        endif
      enddo
      rmaxss=rrout
      zmaxss=znow
      ixyz=-2
66501 xxtras=rmaxss+0.0001_dp
      yxtras=zmaxss
      ixl=1
      if (ipmeas(iges).lt.-1.e3_dp) then
        nerr=10000
      else
        nerr=0
      endif
      call bound(psi,nw,nh,nwnh,psiout,xmin,xmax,ymin,ymax,zeros, &
                 rgrid,zgrid,xxtras,yxtras,ixyz,limtrs,xlims,ylims, &
                 xxtra(1,ixl),yxtra(1,ixl),npxtra(ixl),xlims(1),nxtrap, &
                 rymin,rymax,dpsi,zxmin,zxmax,nerr,ishot,itime, &
                 limfag,radum,kbound,tolbndpsi)
      if (nerr.gt.0) then
        kerror = 1
        return
      end if
      zerovs(1)=1.0
      do i=1,npxtra(ixl)
        zerold=zerovs(1)
        call zlim(zerovs,nnn,nnn,limitr,xlim,ylim,xxtra(i,ixl), &
                  yxtra(i,ixl),limfag)
        if(zerold.gt.0.01_dp.and.zerovs(1).lt.0.01_dp) exit
      enddo
      if (zerold.gt.0.01_dp.and.zerovs(1).lt.0.01_dp.and.i>1) then
        rinvs=xxtra(i-1,ixl)
        zinvs=yxtra(i-1,ixl)
        routvs=xxtra(i,ixl)
        zoutvs=yxtra(i,ixl)
        do i=1,20
          ravs(1)=0.5_dp*(rinvs+routvs)
          zavs(1)=0.5_dp*(zinvs+zoutvs)
          call zlim(zerovs,nnn,nnn,limitr,xlim,ylim,ravs,zavs,limfag)
          if (zerovs(1).gt.0.01_dp) then
            rinvs=ravs(1)
            zinvs=zavs(1)
          else
            routvs=ravs(1)
            zoutvs=zavs(1)
          endif
        enddo
        if ((zavs(1)*zssep.lt.0.0).and.(ixyz.eq.-2)) then
          ixyz=-1
          go to 66501
        endif
        rvsnow=ravs(1)*100.
        zvsnow=zavs(1)*100.
      else
        rvsnow=0.
        zvsnow=0.
      endif
      call seva2d(bkx,lkx,bky,lky,c,ravs(1),zavs(1),pds,ier,n111)
      ssinow=pds(1)
!------------------------------------------------------------------
!--   get distance between inner leg and upper dome              --
!------------------------------------------------------------------
      ycut=ymax*0.5_dp
      ycutm=ymin*0.5_dp
      rsymax=rymax
      rsymin=rymin
      if(zssep.lt.ymin) rsymin=rssep
      if(zssep.gt.ymax) rsymax=rssep
      if (zavs(1)*zssep.lt.0.0) then
        ravssa=0.0
        zavssa=0.0
      else
        ravssa=ravs(1)
        zavssa=zavs(1)
        if (zavs(1).gt.ycut.and.ravs(1).lt.rsymax) then
          rvsiu(iges)=rvsnow
          zvsiu(iges)=zvsnow
          if (ishot.lt.100771) go to 6500
          diludom(iges)=1.e6_dp
          do i=1,npxtra(ixl)
            if (yxtra(i,ixl).ge.ycut) then
              dilnow=(rudom - xxtra(i,ixl))**2 + &
                          (zudom - yxtra(i,ixl))**2
              if (dilnow.lt.diludom(iges)) then
                diludom(iges)=dilnow
                zilnow=yxtra(i,ixl)
              endif
            endif
          enddo
          diludom(iges)=100.0*sqrt(diludom(iges))
          if(zilnow.lt.zudom) diludom(iges)=-diludom(iges)
        elseif (zavs(1).gt.ycut.and.ravs(1).gt.rsymax) then
!------------------------------------------------------------------
!--       get distance between outer leg and upper baffle        --
!------------------------------------------------------------------
          rvsou(iges)=rvsnow
          zvsou(iges)=zvsnow
          if (ishot.lt.100771) go to  6500
          dolubaf(iges)=1.e6_dp
          do i=1,npxtra(ixl)
            if (yxtra(i,ixl).ge.ycut) then
              dolnow=(rubaf - xxtra(i,ixl))**2 + &
                          (zubaf - yxtra(i,ixl))**2
              if (dolnow.lt.dolubaf(iges)) then
                dolubaf(iges)=dolnow
                rolnow=xxtra(i,ixl)
              endif
            endif
          enddo
          dolubaf(iges)=100.0*sqrt(dolubaf(iges))
          if(rolnow.gt.rubaf) dolubaf(iges)=-dolubaf(iges)
        elseif (zavs(1).lt.ycutm.and.ravs(1).gt.rsymin) then
!------------------------------------------------------------------
!--       get lower strike points                                --
!------------------------------------------------------------------
          rvsod(iges)=rvsnow
          zvsod(iges)=zvsnow
        elseif (zavs(1).lt.ycutm.and.ravs(1).lt.rsymin) then
          rvsid(iges)=rvsnow
          zvsid(iges)=zvsnow
        endif
      endif
!------------------------------------------------------------------
!--   second separatrix distances from inboard side              --
!------------------------------------------------------------------
      ixyz=-1
      rnow=rssep
      znow=zssep
      psinow=sissep
      nmax=(xmin-rgrid(1))/(rgrid(2)-rgrid(1))+1
      iznow=(zxmin-zgrid(1)+1.0e-6_dp)/(zgrid(2)-zgrid(1))+1
      znow=zgrid(iznow)
      rrout=xmin
      psiout=psibry
      dsiout=psinow-psiout
      do i=1,nmax
       k=nmax-i+1
       rrin=rgrid(k)
       call seva2d(bkx,lkx,bky,lky,c,rrin,znow,pds,ier,n111)
       psiin=pds(1)
       dsiin=psinow-psiin
       if (dsiin*dsiout.le.0.0) then
         do n=1,20
          rmm=0.5_dp*(rrin+rrout)
          call seva2d(bkx,lkx,bky,lky,c,rmm ,znow,pds,ier,n111)
          psimm=pds(1)
          dsimm=psinow-psimm
          if (dsimm*dsiin.le.0.0) then
            rrout=rmm
          else
            rrin=rmm
          endif
         enddo
         exit
       endif
      enddo
      rminss=rrout
      zminss=znow
66199 xxtras=rminss-0.0001_dp
      yxtras=zminss
      ixl=1
      if (ipmeas(iges).lt.-1.e3_dp) then
        nerr=10000
      else
        nerr=0
      endif
      call bound(psi,nw,nh,nwnh,psiout,xmin,xmax,ymin,ymax,zeros, &
                 rgrid,zgrid,xxtras,yxtras,ixyz,limtrs,xlims,ylims, &
                 xxtra(1,ixl),yxtra(1,ixl),npxtra(ixl),xlims(1),nxtrap, &
                 rymin,rymax,dpsi,zxmin,zxmax,nerr,ishot,itime, &
                 limfag,radum,kbound,tolbndpsi)
      if (nerr.gt.0) then
        kerror = 1
        return
      endif
      zerovs(1)=1.0
      do i=1,npxtra(ixl)
        zerold=zerovs(1)
        call zlim(zerovs,nnn,nnn,limitr,xlim,ylim,xxtra(i,ixl), &
                  yxtra(i,ixl),limfag)
        if(zerold.gt.0.01_dp.and.zerovs(1).lt.0.01_dp) exit
      enddo
      if (zerold.gt.0.01_dp.and.zerovs(1).lt.0.01_dp.and.i.gt.1) then
        rinvs=xxtra(i-1,ixl)
        zinvs=yxtra(i-1,ixl)
        routvs=xxtra(i,ixl)
        zoutvs=yxtra(i,ixl)
        do i=1,20
          ravs(1)=0.5_dp*(rinvs+routvs)
          zavs(1)=0.5_dp*(zinvs+zoutvs)
          call zlim(zerovs,nnn,nnn,limitr,xlim,ylim,ravs,zavs,limfag)
          if (zerovs(1).gt.0.01_dp) then
            rinvs=ravs(1)
            zinvs=zavs(1)
          else
            routvs=ravs(1)
            zoutvs=zavs(1)
          endif
        enddo
        if (((abs(ravssa-ravs(1)).le.1.e-04_dp).and. &
          (abs(zavssa-zavs(1)).le.1.e-04_dp)).and.(ixyz.eq.-1)) then
          ixyz=-2
          go to 66199
        endif
        if ((zavs(1)*zssep.lt.0.0).and.(ixyz.eq.-1)) then
          ixyz=-2
          go to 66199
        endif
        rvsnow=ravs(1)*100.
        zvsnow=zavs(1)*100.
      else
        rvsnow=0.
        zvsnow=0.
      endif
!------------------------------------------------------------------
!--   get distance between inner leg and upper dome              --
!------------------------------------------------------------------
      if (zavs(1)*zssep.ge.0.0) then
        ycut=ymax*0.5_dp
        ycutm=ymin*0.5_dp
        if (zavs(1).gt.ycut.and.ravs(1).lt.rsymax) then
          rvsiu(iges)=rvsnow
          zvsiu(iges)=zvsnow
          if (ishot.ge.100771) then
            diludom(iges)=1.e6_dp
            do i=1,npxtra(ixl)
              if (yxtra(i,ixl).ge.ycut) then
                dilnow=(rudom - xxtra(i,ixl))**2 + &
                          (zudom - yxtra(i,ixl))**2
                if (dilnow.lt.diludom(iges)) then
                  diludom(iges)=dilnow
                  zilnow=yxtra(i,ixl)
                endif
              endif
            enddo
            diludom(iges)=100.0*sqrt(diludom(iges))
            if(zilnow.lt.zudom) diludom(iges)=-diludom(iges)
          endif
        elseif (zavs(1).gt.ycut.and.ravs(1).gt.rsymax) then
!------------------------------------------------------------------
!--       get distance between outer leg and upper baffle          --
!------------------------------------------------------------------
          rvsou(iges)=rvsnow
          zvsou(iges)=zvsnow
          if (ishot.ge.100771) then
            dolubaf(iges)=1.e6_dp
            do i=1,npxtra(ixl)
              if (yxtra(i,ixl).ge.ycut) then
                dolnow=(rubaf - xxtra(i,ixl))**2 + &
                          (zubaf - yxtra(i,ixl))**2
                if (dolnow.lt.dolubaf(iges)) then
                  dolubaf(iges)=dolnow
                  rolnow=xxtra(i,ixl)
                endif
              endif
            enddo
            dolubaf(iges)=100.0*sqrt(dolubaf(iges))
            if(rolnow.gt.rubaf) dolubaf(iges)=-dolubaf(iges)
          endif
        elseif (zavs(1).lt.ycutm.and.ravs(1).gt.rsymin) then
!------------------------------------------------------------------
!--       get lower strike points                                  --
!------------------------------------------------------------------
          rvsod(iges)=rvsnow
          zvsod(iges)=zvsnow
        elseif (zavs(1).lt.ycutm.and.ravs(1).lt.rsymin) then
          rvsid(iges)=rvsnow
          zvsid(iges)=zvsnow
        endif
      endif
 6500 continue
      endif ! sissep.gt.-1.e10_dp
      endif ! (dpsi.gt.psitol) .and. (limloc(iges).ne.'IN')
!
      if (rvsin(iges).gt.rvsout(iges)) then
        rtemp=rvsin(iges)
        ztemp=zvsin(iges)
        rvsin(iges)=rvsout(iges)
        zvsin(iges)=zvsout(iges)
        rvsout(iges)=rtemp
        zvsout(iges)=ztemp
      endif
!---------------------------------------------------------------------------
!--   compute distances at outboard and inboard sides                     --
!---------------------------------------------------------------------------
      if (diludom(iges).gt.0.0) then
        rnow=rudom
        znow=zudom
        call seva2d(bkx,lkx,bky,lky,c,rnow,znow,pds,ier,n111)
        psinow=pds(1)
        nmax=(xmin-rgrid(1))/(rgrid(2)-rgrid(1))+1
        znow=zxmin
        rrout=xmin
        psiout=psibry
        dsiout=psinow-psiout
        do i=1,nmax
          k= nmax-i +1
          rrin=rgrid(k)
          call seva2d(bkx,lkx,bky,lky,c,rrin,znow,pds,ier,n111)
          psiin =pds(1)
          dsiin =psinow-psiin
          if (dsiin*dsiout.le.0.0) then
            do n=1,20
              rmm=0.5_dp*(rrin+rrout)
              call seva2d(bkx,lkx,bky,lky,c,rmm ,znow,pds,ier,n111)
              psimm=pds(1)
              dsimm=psinow-psimm
              if (dsimm*dsiin.le.0.0) then
                rrout=rmm
              else
                rrin=rmm
              endif
            enddo
            exit
          endif
        enddo
        diludomm(iges)=100.0*(xmin-rrin)
      endif
!
      if (dolubaf(iges).gt.0.0) then
        rnow=rubaf
        znow=zubaf
        call seva2d(bkx,lkx,bky,lky,c,rnow,znow,pds,ier,n111)
        psinow=pds(1)
        nmin=(xmax-rgrid(1))/(rgrid(2)-rgrid(1))+2
        znow=zxmax
        rrin=xmax
        psiin=psibry
        dsiin =psinow-psiin
        do k=nmin,nw
          rrout=rgrid(k)
          call seva2d(bkx,lkx,bky,lky,c,rrout,znow,pds,ier,n111)
          psiout=pds(1)
          dsiout=psinow-psiout
          if (dsiin*dsiout.le.0.0) then
            do n=1,20
              rmm=0.5_dp*(rrin+rrout)
              call seva2d(bkx,lkx,bky,lky,c,rmm ,znow,pds,ier,n111)
              psimm=pds(1)
              dsimm=psinow-psimm
              if (dsimm*dsiin.le.0.0) then
                rrout=rmm
              else
                rrin=rmm
              endif
            enddo
            exit
          endif
        enddo
        dolubafm(iges)=100.*(rrout-xmax)
      endif
!---------------------------------------------------------------------------
!--   compute flux expansion parameter                                    --
!---------------------------------------------------------------------------
      fexpan=-10.
      fexpvs=-10.
      if (limloc(iges).eq.'DN '.or.limloc(iges).eq.'SNB'.or. &
          limloc(iges).eq.'SNT') then
        fxrmin=rexpmx-rexpan
        fxzmin=zexpmx
        if (limloc(iges).eq.'DN '.or.limloc(iges).eq.'SNB') then
          ixyz=-2
          zexpmx=zexpmx+0.00001_dp
          fxrmax=rseps(1,iges)/100.
          fxzmax=zseps(1,iges)/100.
        else
          ixyz=-1
          zexpmx=zexpmx-0.0001_dp
          fxrmax=rseps(2,iges)/100.
          fxzmax=zseps(2,iges)/100.
        endif
        if (ipmeas(iges).lt.-1.e3_dp) then
          nerr=10000
        else
          nerr=0
        endif
        ixl=1
        call bound(psi,nw,nh,nwnh,psiout,xmin,xmax,ymin,ymax,zeros, &
                   rgrid,zgrid,rexpmx,zexpmx,ixyz,limtrs,xlims,ylims, &
                   xxtra(1,ixl),yxtra(1,ixl),npxtra(ixl),xlims(1),nxtrap, &
                   rymin,rymax,dpsi,zxmin,zxmax,nerr,ishot,itime, &
                   limfag,radum,kbound,tolbndpsi)
        if (nerr.gt.0) then
          kerror = 1
          return
        end if
!-------------------------------------------------------------------------
!--     get engineering slot parameter in cm, SEPNOSE                   --
!-------------------------------------------------------------------------
        sepnose=-999.
        do i=1,npxtra(ixl)-1
          if ((yxtra(i,ixl)-znose)*(yxtra(i+1,ixl)-znose).le.0.0) then
            rlinnose=xxtra(i,ixl)+(xxtra(i+1,ixl)-xxtra(i,ixl))/ &
             (yxtra(i+1,ixl)-yxtra(i,ixl))*(znose-yxtra(i,ixl))
            sepnose=(rlinnose-rsepnose)*100.
            if(rsepnose.le.0.0) sepnose=-999.
            exit
          endif
        enddo
        fxmin=1.e30_dp
        fxmax=1.e30_dp
        fvsmax=1.e30_dp
        fsrmax=rvsout(iges)/100.
        fszmax=zvsout(iges)/100.
        do i=1,npxtra(ixl)
          fminow=sqrt((fxrmin-xxtra(i,ixl))**2+(fxzmin-yxtra(i,ixl))**2)
          fmanow=sqrt((fxrmax-xxtra(i,ixl))**2+(fxzmax-yxtra(i,ixl))**2)
          fvsnow=sqrt((fsrmax-xxtra(i,ixl))**2+(fszmax-yxtra(i,ixl))**2)
          fxmin=min(fminow,fxmin)
          fxmax=min(fmanow,fxmax)
          fvsmax=min(fvsnow,fvsmax)
        enddo
        if (fxmin.gt.1.e-30_dp) then
          fexpan=fxmax/fxmin
          fexpvs=fvsmax/fxmin
          if(fszmax.ge.0.0) fexpvs=-999.
        endif
      endif
!---------------------------------------------------------------------------
!--   trace external field lines                                          --
!---------------------------------------------------------------------------
      extra_surfs: if (nextra.ne.0) then
      if (ixstrt.ne.-1) then
        xxtraa=xmax
        yxtraa=zxmax
      else
        xxtraa=xmin
        yxtraa=zxmin
      endif
      ixtpls=0
      ixyz=-1
  420 continue
      dxtra=scrape/iabs(nextra)
      do kkk=1,iabs(ixstrt)
       do i=1,iabs(nextra)
        ixl=i+ixtpls*iabs(ixstrt)+(kkk-1)*iabs(nextra)
        scraps(ixl)=1000.*i*dxtra
        if (ixstrt.eq.2) then
         if (kkk.eq.2) then
          xxtraa=xmin
          yxtraa=zxmin
         else
          xxtraa=xmax
          yxtraa=zxmax
         endif
        endif
        xxtras=xxtraa+dxtra*i
        yxtras=yxtraa
        if (kkk.eq.2.or.ixstrt.eq.-1) then
         xxtras=xxtraa-dxtra*i
         yxtras=yxtraa
        endif
        if(xxtras.le.rgrid(1).or.xxtras.ge.rgrid(nw)) cycle
        if (ipmeas(iges).lt.-1.e3_dp) then
         nerr=10000
        else
         nerr=0
        endif
        call bound(psi,nw,nh,nwnh,psiout,xmin,xmax,ymin,ymax,zeros, &
                   rgrid,zgrid,xxtras,yxtras,ixyz,limtrs,xlims,ylims, &
                   xxtra(1,ixl),yxtra(1,ixl),npxtra(ixl),xlims(1),nxtrap, &
                   rymin,rymax,dpsi,zxmin,zxmax,nerr,ishot,itime, &
                   limfag,radum,kbound,tolbndpsi)
        if (nerr.gt.0) then
         kerror = 1
         return
        endif
        if ((i.le.1).and.(ixyz.eq.-2).and.(xlimxs.ge.0.0)) then
         do n=1,npxtra(ixl)
          sepnow=sqrt((xlimxs-xxtra(n,ixl))**2+(ylimys- &
                       yxtra(n,ixl))**2)
          sepexp(iges)=min(abs(sepexp(iges)),sepnow)
         enddo
         sepexp(iges)=sepexp(iges)*100.
        endif
        if(nextra.gt.0) cycle
        do n=1,npxtra(ixl)
         call seva2d(bkx,lkx,bky,lky,c,xxtra(n,ixl),yxtra(n,ixl), &
                     pds,ier,n333)
         bpxtra(n,ixl)=sqrt(pds(2)**2+pds(3)**2)/xxtra(n,ixl)
        enddo
        fpxtra(1,ixl)=0.0
        flxtra(1,ixl)=0.0
        do n=2,npxtra(ixl)
         nm1=n-1
         dlpol=sqrt((xxtra(n,ixl)-xxtra(nm1,ixl))**2 &
                    +(yxtra(n,ixl)-yxtra(nm1,ixl))**2)
         btnow=abs(fbrdy)*tmu/2.*(1./xxtra(n,ixl)+1./xxtra(nm1,ixl))
         bpnow=(bpxtra(n,ixl)+bpxtra(nm1,ixl))/2.
         dltol=dlpol*btnow/bpnow
         dlll=sqrt(dlpol**2+dltol**2)
         fpxtra(n,ixl)=fpxtra(nm1,ixl)+dlpol
         flxtra(n,ixl)=flxtra(nm1,ixl)+dlll
        enddo
       enddo
      enddo
      if (ixyz.ne.-2) then
       ixyz=-2
       ixtpls=iabs(nextra)
       go to 420
      endif
      endif extra_surfs
!-----------------------------------------------------------------------------
!--   compute the diamagnetic flux                                          --
!-----------------------------------------------------------------------------
      vbtot2=0.0
      vbtvac2=0.0
      vbtvac=0.0
      vbtor2=0.0
      volbt=0.0
      cdflux(iges)=0.0
      edflux(iges)=0.0
      btvtot2=0.0
      btvvac2=0.0
      btvtor2=0.0
      volrt=0.0
      areart=0.0
      sumbz2=0.0
      sumbp2=0.0
!
      call zpline(nw,xsisii,pres,bpres,cpres,dpres)
      if (kdopre.gt.0) then
        ii=0
        do i=1,nw,kdopre
          ii=ii+1
          rpress(ii)=-xsisii(i)
          pressr(ii)=pres(i)
          sigpre(ii)=0.1_dp*pres(i)
        enddo
        if (abs(rpress(ii)).ne.xsisii(nw)) then
          ii=ii+1
          i=nw
          rpress(ii)=-xsisii(i)
          pressr(ii)=pres(i)
          sigpre(ii)=0.1_dp*pres(i)
        endif
        npress=ii
      endif
!
      if(kvtor.gt.0) &
        call zpline(nw,xsisii,pressw,bpresw,cpresw,dpresw)
      delerr=0.0
      delerb=0.0
      do i=1,nw
       do j=1,nh
        kk=(i-1)*nh+j
        fnow=fbrdy
        presss=0.0
        prewww=0.0
        prettt=presss
        if (xpsi(kk).le.1.0) then
         fnow=seval(nw,xpsi(kk),sigrid,fpol,bbfpol,ccfpol,ddfpol)
         presss=seval(nw,xpsi(kk),xsisii,pres,bpres,cpres,dpres) &
                -prbdry
         fnow=fnow/tmu
         prettt=presss
         rotation_2: if (kvtor.gt.0) then
           preww0=seval(nw,xpsi(kk),xsisii,pressw,bpresw,cpresw,dpresw) &
                  -preswb
           press0=presss
           presss=press0+preww0*rgrvt(i)
           prewww=preww0*(rgrid(i)/rvtor)**2*2.
           prew0=preww0
           pres0=press0
           if (abs(pres0).gt.1.e-10_dp) then
             pwop0=prew0/pres0
             pwp0r2=pwop0*rgrvt(i)
           else
             pwop0=0.0
             pwp0r2=0.0
           endif
           if (kvtor.eq.2) then
             presss=press0*(1.+0.5_dp*pwp0r2**2)+preww0*rgrvt(i)
           elseif (kvtor.eq.3) then
             presss=press0*exp(pwp0r2)
           endif
           prettt=prewww+presss
         endif rotation_2
        endif
        cdflux(iges)=cdflux(iges)+(fbrdy-fnow)/rgrid(i)*www(kk)
        edflux(iges)=edflux(iges)+(fbrdy**2-fnow**2)/rgrid(i)*www(kk)
        call seva2d(bkx,lkx,bky,lky,c,rgrid(i),zgrid(j),pds,ier,n666)
        bpzsq=pds(2)**2/rgrid(i)
        bpolsq=bpzsq+pds(3)**2/rgrid(i)
        sumbz2=sumbz2+bpzsq*www(kk)
        sumbp2=sumbp2+bpolsq*www(kk)
        bpolsq=bpolsq/rgrid(i)
        btvac2=(tmu*fbrdy/rgrid(i))**2
        btttt2=(tmu*fnow/rgrid(i))**2
        enrgy=(2.*twopi*tmu*prettt+bpolsq+btvac2-btttt2)*www(kk)
        volrt=volrt+rgrid(i)*enrgy
        areart=areart+enrgy
        vbtot2=vbtot2+rgrid(i)*www(kk)*(bpolsq+btttt2)
        vbtvac2=vbtvac2+rgrid(i)*www(kk)*btvac2
        vbtvac=vbtvac+rgrid(i)*www(kk)*sqrt(btvac2)
        vbtor2=vbtor2+rgrid(i)*www(kk)*btttt2
        volbt=volbt+rgrid(i)*www(kk)
        btvvac2=btvvac2+2.*twopi*tmu*presss*www(kk)/btvac2*rgrid(i)
        btvtor2=btvtor2+2.*twopi*tmu*presss*www(kk)/btttt2*rgrid(i)
        btvtot2=btvtot2+2.*twopi*tmu*presss*www(kk)/(bpolsq+btttt2)* &
                rgrid(i)
!--------------------------------------------------------------------------
!--     evaluate -del*psi/R/mu0                                          --
!--------------------------------------------------------------------------
        if(i.eq.1.or.i.eq.nw) cycle
        if(j.eq.1.or.j.eq.nh) cycle
        kip=i*nh+j
        kim=(i-2)*nh+j
        kjp=(i-1)*nh+j+1
        kjm=(i-1)*nh+j-1
        d2sidr2=(psipold(kip)-2.*psipold(kk)+psipold(kim))/drgrid**2
        d2sidz2=(psipold(kjp)-2.*psipold(kk)+psipold(kjm))/dzgrid**2
        dsidr=(psipold(kip)-psipold(kim))/2./drgrid
        delsbu=d2sidr2-dsidr/rgrid(i)+d2sidz2
        delsbu=-delsbu/rgrid(i)/tmu0
        pcnow=pcurrt(kk)/darea
        delerrb=abs(delsbu-pcnow)
        delerb=max(delerrb,delerb)
!--------------------------------------------------------------------------
!--     total psi for inside plasma only                                 --
!--------------------------------------------------------------------------
        if(xpsi(kk).gt.0.999_dp) cycle
        d2sidr2=(psiold(kip)-2.*psiold(kk)+psiold(kim))/drgrid**2
        d2sidz2=(psiold(kjp)-2.*psiold(kk)+psiold(kjm))/dzgrid**2
        dsidr=(psiold(kip)-psiold(kim))/2./drgrid
        delssi=d2sidr2-dsidr/rgrid(i)+d2sidz2
        delssi=-delssi/rgrid(i)/tmu0
        delerrx=abs(delssi-pcnow)
        delerr=max(delerrx,delerr)
       enddo
      enddo
      curnow=abs(ipmhd(iges))/area(iges)*1.e4_dp
      delerr=delerr/curnow
      delerb=delerb/curnow
      alpha(iges)=2*sumbz2/sumbp2
      rttt(iges)=volrt/areart*100.0
      vbtot2=vbtot2/volbt
      vbtvac2=vbtvac2/volbt
      vbtor2=vbtor2/volbt
      vbtvac=vbtvac/volbt
      btuse=bcentr(iges)*rcentr*100./rout(iges)
      btuse=btuse**2
      vbtot2=betat(iges)*btuse/vbtot2
      vbtvac2=betat(iges)*btuse/vbtvac2
      vbtor2=betat(iges)*btuse/vbtor2
      vbtvac=betat(iges)*btuse/vbtvac**2
      btvvac2=btvvac2*100./volbt
      btvtor2=btvtor2*100./volbt
      btvtot2=btvtot2*100./volbt
      vbtmag=betat(iges)*(rmaxis*100./rout(iges))**2
      betped=0.0
      betnped=0.0
      if (kedgep.gt.0) then
        sssiie=pe_psin-pe_width
        presped=seval(nw,sssiie,xsisii,pres,bpres,cpres,dpres) &
                -prbdry
        betped=presped/btuse*tmu*twopi*2.*100.
        betnped=betped/betat(iges)*betatn
      endif
!--------------------------------------------------------------------
!--   change sign convention of computed diamagnetic flux          --
!--   to be consistent with DIII-D measurements.  The sign now     --
!--   is opposite to that in Nuc Fusion paper      06/23/87        --
!--------------------------------------------------------------------
      cdflux(iges)=-cdflux(iges)*darea*tmu
      edflux(iges)=-edflux(iges)*darea*tmu**2
      sbpli=(s1(iges)+s2(iges)*(1.+rttt(iges)*0.01_dp/rcentr))/4.0
      vbpli=betap(iges)+li(iges)/2.
      if (kvtor.gt.0) vbpli=vbpli+betapw(iges)
      dbpli(iges)=abs((sbpli-vbpli)/sbpli)
      chidflux=0.0
      if (fwtdlc.gt.0.0) then
        chidflux=((diamag(iges)-cdflux(iges))/sigdia(iges))**2
      endif
      chitot=chipre+chifin+chidflux
      if (idiart.gt.0) then
        bp2flx=bpolav(iges)**2
        dmui=1.0e+06_dp*diamag(iges)*4.*pi*bcentr(iges)*rcentr &
              /bp2flx/volume(iges)
        rcurrm=rttt(iges)*0.01_dp
        betapd(iges)=s1(iges)/2.+s2(iges)/2.*(1.-rcurrm/rcentr)-dmui
        betatd(iges)=betapd(iges)*bp2flx/bcentr(iges)**2*100.
        betatd(iges)=betatd(iges)*(rout(iges)/100./rcentr)**2
        wdia(iges)=1.5_dp*betapd(iges)*bp2flx/2./tmu/2./pi*volume(iges) &
                      /1.0e+06_dp
        if (kdata.ne.4) then
          pohm=vloopt(iges)*ipmhd(iges)
          pohm=max(czero,pohm)
          pin=pohm+pbinj(iges)
          if (pin.gt.1.e-06_dp) then
            taudia(iges)=wdia(iges)/pin*1000.
          endif
        endif
      endif
!
      xmui=2.*twopi*bcentr(iges)*rcentr*cdflux(iges)/volume(iges) &
           *1.0e+06_dp/bpolav(iges)**2
      exmui=twopi*edflux(iges)/volume(iges) &
            *1.0e+06_dp/bpolav(iges)**2
      sbppa=(s1(iges)+s2(iges)*(1.-rttt(iges)*0.01_dp/rcentr))/2.-xmui
      sbpp=(s1(iges)+s2(iges)*(1.-rttt(iges)*0.01_dp/rcentr))/2.-exmui
      delbp(iges)=abs((sbpp-betap(iges))/sbpp)
      if (elong(iges).gt.1.2_dp) then
        sbli=(s1(iges)/2.+s2(iges)/2.*(1.-rttt(iges)*0.01_dp/rcentr)- &
              s3(iges))/(alpha(iges)-1.)
        delli=abs((sbli-li(iges))/li(iges))
        delbp(iges)=max(delbp(iges),delli)
      endif
!
      if (nbskip.le.0) then
        nbabs=1
      else
        nbabs=max(nfound/min(mbdry,nbdrymx)+1,nbskip)
      endif
      jb=0
      iskip=nbabs-1
!-------------------------------------------------------------------
!--   skip points if needed but make sure X point and last point  --
!--   are included                                                --
!-------------------------------------------------------------------
      do i=1,nfound
        iskip=iskip+1
        dzzz1=abs(zseps(1,iges)/100.-yout(i))
        dzzz2=abs(zseps(2,iges)/100.-yout(i))
        dzzz1=min(dzzz1,dzzz2)
        if ((iskip.eq.nbabs.or.dzzz1.le.1.e-4_dp).or.(i.eq.nfound)) &
          then
          jb=jb+1
          rbbbs(jb)=xout(i)
          zbbbs(jb)=yout(i)
          iskip=0
        endif
      enddo
      nbbbs=jb
!--------------------------------------------------------------------
!--   rmidin and rmidout are boundary R (rbbbs) at Z (zbbbs) = 0.0 --
!--------------------------------------------------------------------
      i2=1
      do i=1,nbbbs-1
        if (zbbbs(i).eq.0.0) then
          rmid2(i2)= rbbbs(i)
          i2=i2+1
        elseif (zbbbs(i)*zbbbs(i+1).lt.0.0) then
          rmid2(i2)=(rbbbs(i)+rbbbs(i+1))/2.0
          i2=i2+1
        endif
        if(i2.gt.2) exit
      enddo
      rmidin(iges) = min(rmid2(1),rmid2(2))
      rmidout(iges)= max(rmid2(1),rmid2(2))
!
      if (kprfit.gt.0) then
        do i=1,npress
          if(rpress(i).le.0.0) cycle
          call seva2d(bkx,lkx,bky,lky,c,rpress(i),zpress(i),pds,ier,n111)
          rpress(i)=-(simag-pds(1))/sidif
        enddo
      endif
!----------------------------------------------------------------
!--   Current at Z=Z_Libeam                                    --
!----------------------------------------------------------------
      kwripre_li: if (kwripre.gt.0) then
       if (nstark.gt.nmselp) then
        znow=zzgam(iges,nmselp+1)
        delrnow=(xmax-rmaxis)/(nw-1)
        do i=1,nw
           rnow=rmaxis+(i-1)*delrnow
           call seva2d(bkx,lkx,bky,lky,c,rnow,znow,pds,ier,n666)
           rjtli(i)=rnow
           sjtliz(i)=(pds(2)/rnow-pds(5))/rnow/twopi/tmu/1000.
           sjtlir(i)=-pds(6)/rnow/twopi/tmu/1000.
           sjtli(i)=sjtlir(i)+sjtliz(i)
           xpsikk=(pds(1)-simag)/(psibry-simag)
           if(xpsikk.lt.0.0) xpsikk=0.0
           cjtli(i)=rnow*ppcurr(xpsikk,kppcur) &
                      +fpcurr(xpsikk,kffcur)/rnow
           cjtli(i)=cjtli(i)/darea/1000.
        enddo
        call setfnmd('q',ishot,itime,sfname)
        sfname=sfname(1:13)//'_jtorli'
        open(unit=74,status='old',file=sfname,iostat=ier)
        if(ier.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        xdum=0.0
        do i=1,nw
          write (74,*) rjtli(i),cjtli(i),xdum,xdum
        enddo
        close(unit=74)
        call setfnmd('q',ishot,itime,sfname)
        sfname=sfname(1:13)//'_jtsli'
        open(unit=74,status='old',file=sfname,iostat=ier)
        if(ier.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        xdum=0.0
        do i=1,nw
          write (74,*) rjtli(i),sjtli(i),xdum,xdum
        enddo
        close(unit=74)
        call setfnmd('q',ishot,itime,sfname)
        sfname=sfname(1:13)//'_jtslir'
        open(unit=74,status='old',file=sfname,iostat=ier)
        if(ier.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=1,nw
          write (74,*) rjtli(i),sjtlir(i),xdum,xdum
        enddo
        close(unit=74)
        call setfnmd('q',ishot,itime,sfname)
        sfname=sfname(1:13)//'_jtsliz'
        open(unit=74,status='old',file=sfname,iostat=ier)
        if(ier.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=1,nw
          write (74,*) rjtli(i),sjtliz(i),xdum,xdum
        enddo
        close(unit=74)
       endif
      endif kwripre_li
!
      if (dco2v(iges,2).gt.1.e+10_dp) then
        partic=dco2v(iges,2)*volume(iges)
        if(partic.gt.1.0e+10_dp) tave(iges)=wmhd(iges)/partic/3. &
                                                         /1.602e-16_dp
        resist=vloopt(iges)/ipmhd(iges)
        zeta=resist*area(iges)/twopi/rout(iges)
        xlam=20.
        if (tave(iges).gt.0.001_dp) then
          tevolt=sqrt((tave(iges)*1000.)**3)
        endif
        zeffr(iges)=zeta*tevolt/xlam/1.03e-02_dp*2.
        if (iges.eq.igmax) then
          do m=1,igmax
            temp(m)=dco2v(m,2)
          enddo
          if (size(time,1).gt.1) then
            dttt=time(2)-time(1)
          else
            dttt=0.0
          endif
          if(igmax.eq.1) dttt=5.
          ! DEPRECATED (no longer implemented)
          !call getzeff(ishot,igmax,time(1),dttt,temp,tave,zeff,ier)
        endif
      endif
!------------------------------------------------------------------
!--   compute vessel forces                                      --
!------------------------------------------------------------------
      vessel_force: if (ivesel.eq.3.or.icutfp.eq.2) then
        call sets2d(psi,c,rgrid,nw,bkx,lkx,zgrid,nh,bky,lky,wk,ier)
        fztor=0.0
        fzpol=0.0
        do i=1,nvesel
          signz=zvs(i)/abs(zvs(i))
          signr=(rvs(i)-rcentr)/abs(rvs(i)-rcentr)
          if (avs(i).eq.0.0.and.avs2(i).eq.0.0) then
            if (wvs(i).lt.hvs(i)) then
              cosalp=0.0
              dells=hvs(i)
            else
              cosalp=1.0*signz
              dells=wvs(i)
            endif
          endif
          if (avs(i).ne.0.) then
            cosalp=abs(cos(avs(i))*180./pi)*signz
            dells=sqrt(hvs(i)**2+wvs(i)**2)
          endif
          if (avs2(i).ne.0.) then
            cosalp=abs(cos(avs2(i))*180./pi)*signz
            dells=sqrt(hvs(i)**2+wvs(i)**2)
          endif
          sinalp=signr*sqrt(1.-cosalp**2)
          rsnow0=rvs(i)-(wvs(i)/2.+wvs(i)/40.)*signz
          zsnow0=zvs(i)+(hvs(i)/2.+hvs(i)/40.)*signr
          sbrrs=0.0
          sumfzp=0.0
          do k=1,20
            rsnow=rsnow0+wvs(i)*k/20.*signz
            zsnow=zsnow0-hvs(i)*k/20.*signr
            call seva2d(bkx,lkx,bky,lky,c,rsnow,zsnow,pds,ier,n333)
            xpsivs=(simag-pds(1))/sidif
            sbrrs=sbrrs-pds(3)
            if (pds(1).lt.psibry.and. &
               (abs(zsnow).lt.abs(yvs2).or.zsnow*yvs2.lt.0.0)) then
              fpnow=ffcurr(xpsivs,kffcur)
            else
              fpnow=fbrdy
            endif
            delfp=fpnow -fbrdy
            btnow=fpnow*tmu/rsnow
            sumfzp=sumfzp+btnow*delfp
          enddo
          vforcet(i)=-sbrrs/20.*twopi*vcurrt(i)
          fztor=fztor+vforcet(i)
          vforcep(i)= cosalp*sumfzp*dells/20.
          fzpol=fzpol+vforcep(i)
        enddo
      endif vessel_force
      if (icutfp.eq.2) then
        xxxx=1./xpsimin
        fpolvs=ffcurr(xxxx,kffcur)-ffcurr(x111,kffcur)
      endif
!
! --- prepare for vertical stability parameters
! --- calculation is done after pltout calls
!
      pleng=0.0
      do i=1,nfound-1
        ip1=i+1
        dli=sqrt((xout(ip1)-xout(i))**2+(yout(ip1)-yout(i))**2)
        pleng=pleng+dli
      enddo
      abar=100.*pleng/2./pi
!
#ifdef DEBUG_LEVEL1
      write (6,*) 'Call SHAPE/PLTOUT kerror = ', kerror
      write (6,*) 'ivesel, icutfp, fpolvs(ka) = ', ivesel,icutfp,fpolvs/1000.
#endif
      if (itek.gt.0) then
        if (idplace.eq.0) then
#ifdef DEBUG_LEVEL1
          write (6,*) 'Before SHAPE/PLTOUT'
#endif
          call pltout(xout,yout,nfound,iges,nnn, &
                      xmin,xmax,ymin,ymax,igmax,kerror)
          if(kerror.gt.0) return
#ifdef DEBUG_LEVEL1
          write (6,*) 'After SHAPE/PLTOUT'
#endif
        endif
      endif
#ifdef DEBUG_LEVEL1
      write (6,*) 'exit SHAPE/PLTOUT kerror = ', kerror
#endif
      if ((itrace.gt.1) .and. (abs(dpsi).le.psitol)) then
        if (itek.gt.0) then
           if (idplace.eq.0) then
             call pltout(xouts,youts,nfouns,jges,n22, &
                         xmins,xmaxs,ymins,ymaxs,igmax,kerror)
             if(kerror.gt.0) return
           endif
        endif
      endif
      if((itek.ge.5).and.(iges.eq.igmax)) close(unit=35)
!-----------------------------------------------------------------------
!--   vertical stability parameter,  reference Nuc Fusion  18(1978)1331
!--   move out of pltout so that vertn, xnnc are indepent of itek value
!-----------------------------------------------------------------------
      do i=1,nw
        do j=1,nh
          kk=(i-1)*nh+j
          copy(i,j)=0.0
          do m=1,nfsum
            copy(i,j)=copy(i,j)+gridfc(kk,m)*brsp(m)
          enddo
          if (ivesel.gt.0) then
            do m=1,nvesel
              copy(i,j)=copy(i,j)+gridvs(kk,m)*vcurrt(m)
            enddo
          endif
          if (iecurr.gt.0) then
            do m=1,nesum
              copy(i,j)=copy(i,j)+gridec(kk,m)*ecurrt(m)
            enddo
          endif
        enddo
      enddo
      do i=1,nw
        do j=1,nh
          k=(i-1)*nh+j
          copyn(k)=copy(i,j)
        enddo
      enddo
      call sets2d(copyn,c,rgrid,nw,bkx,lkx,zgrid,nh,bky,lky,wk,ier)
      rcur=rcurrt(jges)/100
      zcur=zcurrt(jges)/100
      call seva2d(bkx,lkx,bky,lky,c,rcur,zcur,pds,ier,n555)
      vertn(jges)=1.-pds(5)/pds(2)*rcur
      rx=rm(jges)/100.
      f_0=log(8*rout(jges)/abar)-2+betap(jges)+li(jges)/2+.5_dp
      delr=rout(jges)/100.-1.67_dp
!-----------------------------------------------------------------------
!--   metal wall                                                    --
!-----------------------------------------------------------------------
      xnnc(jges)=vertn(jges)/((10.77_dp*delr**2+8.08_dp*delr+2.54_dp)/f_0)
!
      else is_vacuum
      limloc(iges)='VAC'
      call chisqr(iges)
      if (itek.gt.0) then
        call pltout(xout,yout,nzz,iges,nnn,zxx,zxx,zxx,zxx,igmax,kerror)
        if (kerror.gt.0) return
      endif
      if ((itek.ge.5).and.(iges.eq.igmax)) close(unit=35)
      ! initialize output variables that weren't set
      zout=0.0
      elong=0.0
      utri=0.0
      ltri=0.0
      qstar=0.0
      gapin=0.0
      gapout=0.0
      gaptop=0.0
      gapbot=0.0
      sepin=0.0
      sepout=0.0
      septop=0.0
      sepbot=0.0
      q95=0.0
      shearb=0.0
      vertn=0.0
      rco2v=0.0
      dco2v=0.0
      rco2r=0.0
      dco2r=0.0
      sibdry=0.0
      area=0.0
      terror=0.0
      elongm=0.0
      qmerci=0.0
      qm=0.0
      cdflux=0.0
      alpha=0.0
      rttt=0.0
      indent=0.0
      rseps=0.0
      zseps=0.0
      sepexp=0.0
      dsep=0.0
      aaq1=0.0
      aaq2=0.0
      aaq3=0.0
      rm=0.0
      zm=0.0
      psim=0.0
      rcurrt=0.0
      zcurrt=0.0
      li=0.0
      li3=0.0
      betap=0.0
      betat=0.0
      wmhd=0.0
      betapd=0.0
      betatd=0.0
      wdia=0.0
      bpolav=0.0
      s1=0.0
      s2=0.0
      s3=0.0
      qout=0.0
      btaxp=0.0
      btaxv=0.0
      zuperts=0.0
      cjor0=0.0
      cjor95=0.0
      cjor99=0.0
      cj1ave=0.0
      pp95=0.0
      drsep=0.0
      yyy2=0.0
      xnnc=0.0
      fexpan=0.0
      qmin=0.0
      rhoqmin=0.0
      ssi01=0.0
      ssi95=0.0
      fexpvs=0.0
      sepnose=0.0
      rmidin=0.0
      rmidout=0.0
      psurfa=0.0
      peak=0.0
      dminux=0.0
      dminlx=0.0
      dolubaf=0.0
      dolubafm=0.0
      diludom=0.0
      diludomm=0.0
      ratsol=0.0
      rvsiu=0.0
      zvsiu=0.0
      rvsid=0.0
      zvsid=0.0
      rvsou=0.0
      zvsou=0.0
      rvsod=0.0
      zvsod=0.0
      psin21=0.0
      psin32=0.0
      rq21top=0.0
      rq32in=0.0
      tflux=0.0
      twagap=0.0
!
      endif is_vacuum
!
!------------------------------------------------------------------
!--   compute shearing rate eshear                               --
!------------------------------------------------------------------
      if (keecur.gt.0) then
        do i=1,nw
          eshear(i)=esradial(sipmid(i),keecur,rpmid(i),zmaxis)
        enddo
      endif
!-----------------------------------------------------------------------
!--   write out r(shot).(time)_X files in rho space                   --
!-----------------------------------------------------------------------
      kwripre_r: if (kwripre.eq.3) then
        call setfnmd('r',ishot,itime,sfname)
        sfname=sfname(1:13)//'_qpsi'
        open(unit=74,status='old',file=sfname,iostat=ier)
        if(ier.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=1,nw
          write (74,*) rhovn(i),qpsi(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_jor'
        open(unit=74,status='old',file=sfname,iostat=ier)
        if(ier.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=1,nw
          cjorka=cjor(i)/1000.
          write (74,*) rhovn(i),cjorka,xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_jorec'
        open(unit=74,status='old',file=sfname,iostat=ier)
        if(ier.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=1,nw
          cjorka=cjorec(i)/1000.
          write (74,*) rhovn(i),cjorka,xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_cjmse'
        open(unit=74,status='old',file=sfname,iostat=ier)
        if(ier.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        if (ishot.le.97400) then
          mcentral=15
        else
          mcentral=10
        endif
        do i=1,mcentral+1
          cjorka=cjmse(i)/1000.
          write (74,*) rhogam(i),cjorka,xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_cjmsec'
        open(unit=74,status='old',file=sfname,iostat=ier)
        if(ier.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=1,mcentral+1
          cjorka=cjmsec(i)/1000.
          write (74,*) rhogam(i),cjorka,xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_bzmse'
        open(unit=74,status='old',file=sfname,iostat=ier)
        if(ier.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=1,nstark
          write (74,*) rhogam(i),bzmse(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_bzmsec'
        open(unit=74,status='old',file=sfname,iostat=ier)
        if(ier.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=1,nstark
          write (74,*) rhogam(i),bzmsec(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_presf'
        open(unit=74,status='old',file=sfname,iostat=ier)
        if(ier.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=1,nw
          write (74,*) rhovn(i),pres(i),xdum,xdum
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_presd'
        open(unit=74,status='old',file=sfname,iostat=ier)
        if(ier.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=1,npress
          xmnow=-rpress(i)
          ymnow=seval(nw,xmnow,sigrid,rhovn,brhovn,crhovn,drhovn)
          write (74,*) ymnow,pressr(i),xdum,sigpre(i)
        enddo
        close(unit=74)
        sfname=sfname(1:13)//'_pbeam'
        open(unit=74,status='old',file=sfname,iostat=ier)
        if(ier.eq.0) close(unit=74,status='delete')
        open(unit=74,status='new',file=sfname)
        do i=1,nbeam
          xmnow=sibeam(i)
          ymnow=seval(nw,xmnow,sigrid,rhovn,brhovn,crhovn,drhovn)
          write (74,*) ymnow,pbeam(i),xdum,xdum
        enddo
        close(unit=74)
        if (keecur.gt.0) then
          sfname=sfname(1:13)//'_wexb'
          open(unit=74,status='old',file=sfname,iostat=ier)
          if(ier.eq.0) close(unit=74,status='delete')
          open(unit=74,status='new',file=sfname)
          do i=1,nw
            if(rpmid(i).ge.rmaxis) &
              write (74,*) rhopmid(i),eshear(i),xdum,xdum
          enddo
          close(unit=74)
          sfname=sfname(1:13)//'_errho'
          open(unit=74,status='old',file=sfname,iostat=ier)
          if(ier.eq.0) close(unit=74,status='delete')
          open(unit=74,status='new',file=sfname)
          do i=1,nw
            if(rpmid(i).ge.rmaxis) &
              write (74,*) rhopmid(i),ermid(i),xdum,xdum
          enddo
          close(unit=74)
        endif
      endif kwripre_r
!
      DEALLOCATE(xsisii,bpres,cpres,dpres,sjtli,sjtlir,sjtliz, &
                 rjtli,bpresw,cpresw,dpresw,copyn,cjtli,x,y)
!
      return
      end subroutine shapesurf


!**********************************************************************
!>
!!    dslant finds the minimum distance between a curve
!!    represented by (x,y) and the line given by (x1,y1)
!!    and (x2,y2).
!!    
!!
!!    @param x : array of R positions along curve
!!
!!    @param y : array of Z positions along curve
!!
!!    @param np : number of points along curve
!!
!!    @param xmin :
!!
!!    @param xmax :
!!
!!    @param ymin :
!!
!!    @param ymax :
!!
!!    @param x1 : endpoint R of first segment 
!!
!!    @param y1 : endpoint Z of first segment
!!
!!    @param x2 : endpoint R of second segment
!!
!!    @param y2 : endpoint Z of second segment
!!
!!    @param dismin : minimum distance
!!
!**********************************************************************
      subroutine dslant(x,y,np,xmin,xmax,ymin,ymax,x1,y1,x2,y2,dismin)
      use set_kinds, only: dp
      use eparm, only: npoint
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      data nn/30/
      real*8, intent(in) :: xmin, xmax, ymin, ymax, x1, y1, x2, y2
      real*8, intent(inout) :: dismin
      real*8, intent(in) :: x(npoint), y(npoint)

      dismin=1.0e+20_dp
      delx=x2-x1
      dely=y2-y1
      dels=sqrt(delx**2+dely**2)
      nn=dels/0.002_dp
      nn=max(5,nn)
      delx=delx/(nn-1)
      dely=dely/(nn-1)
      do i=1,nn
        xw=x1+delx *(i-1)
        yw=y1+dely *(i-1)
        do m=1,np
          if ((x(m).lt.xmin) .or. (x(m).gt.xmax) .or. (y(m).lt.ymin) .or. (y(m).gt.ymax)) cycle
          disw=sqrt((xw-x(m))**2+(yw-y(m))**2)
          dismin=min(dismin,disw)
        end do
      enddo
      dismin=dismin*100.
      return
      end subroutine dslant
