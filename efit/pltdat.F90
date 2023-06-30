#include "config.f"
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          pltout traces out the plasma outermost surface.  This   **
!**          version implements the plotting using DISSPLA.          **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**********************************************************************
      subroutine pltout(xplt,yplt,nplt,jtime,ipass, &
                        rmin,rmax,zmin,zmax,ktime,kerror)
      use commonblocks, only: worka,c,wk,copy,bkx,bky,cw,wkw,copyw,bwx, &
                              bwy,cj,wkj,copyj,bjx,bjy,cv,wkv,copyv,bvx, &
                              bvy,byringr,byringz, &
                              xxtra,yxtra,bpxtra,flxtra,fpxtra
      use efit_bdata, only: iunit,m_write
      use set_kinds, only: i4
      use curve2d_mod
      use var_cww
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      dimension rrgams(nstark),rrgaml(nstark)
      character(10) :: uday, clocktime
      character(5)  :: zone
      integer*4, dimension(8) :: values
      real*8, dimension(2) :: rwstrip1,zwstrip1,rwstrip2,zwstrip2 
      character(50) dataname,plotname
      character iname,let
      character(2) let2
      character(1000) :: lline
      dimension ae(necoil),ae2(necoil),si(nwf),rmpi(magpri)
      dimension anglem(nstark),anglec(nstark)
      dimension xplt(npoint),yplt(npoint),workb(nwf),workc(nwf),workd(nwf), &
                worke(nwf),bworkb(nwf),cworkb(nwf),dworkb(nwf),workk(nwf), &
                workf(nwf),workg(nwf),workh(nwf),rpressv(mpress), &
                bpressv(mpress),cpressv(mpress),dpressv(mpress), &
                qthom(mpress),zqthom(mpress),vnbeam(mpress),kct(4), &
                xiter(kxiter),czmcm(kxiter),workaw(nwf),ptherm(mpress), &
                presv(nwf),bpresv(nwf),cpresv(nwf),dpresv(nwf)
      character(72) text
      dimension ratray(10),expsi(nsilop),expmp(magpri),pds(6)
      dimension scrat(ntime),bscra(ntime),cscra(ntime), &
                dscra(ntime),bwork(ndata),cwork(ndata),dwork(ndata)
      dimension psecem(nece),psecep(nece), &
                xece(ndim_crv,3+nece),yece(ndim,3+nece)
      dimension thcece(3+nece),sclece(3+nece)
      dimension zz(ix,iy)
      real*8, allocatable :: bfield(:,:), &
                             sivol(:),voln(:),bvoln(:),cvoln(:), &
                             dvoln(:),rscrap(:),curscr(:),pbimf(:), &
                             pmid(:),pmidw(:),bpmid(:),cpmid(:),dpmid(:), &
                             workj(:),copyn(:),copy1(:),xbnow(:), &
                             ybnow(:),rat(:)
      character*20 clrece(3+nece)
      integer*4 dshece(3+nece),dotece(3+nece),cdhece(3+nece), &
                cdtece(3+nece),mrkece(3+nece),nplece(3+nece), &
                cntece(3+nece)
!      equivalence (copy(1,1),copy1(1))
      data ratray/2.,1.,1.,2.,2.,1.,1.,1.,2.,1./
      data istrpl/0/
      integer*4, parameter :: lfile=36,n00=0,n11=1,idsep=1
      integer*4, parameter :: magpri67=29,magpri322=31,magprirdp=8
      real*8, parameter :: one=1.

      ifcoil=1
      kerror = 0
      kgrid=1
      kthkcrv=0
      pleng = 0
      tlen = 0
      ! DIII-D default data (will break if nangle or ntangle change)
      xangle = (/120.,124.,128.,132.,136.,140.,144.,148., &
                 152.,156.,160.,164.,168.,172.,176.,180.,180., &
                 184.,188.,192.,196.,200.,204.,208.,212.,216., &
                 220.,224.,228.,232.,236.,240., &
                 294.5,290.5,286.5,282.5,278.5,274.5,270.5,266.5, &
                 262.5,258.5,254.5,250.5,246.5,242.5,238.5,234.5, &
                 234.5,230.5,226.5,222.5,218.5,214.5,210.5,206.5, &
                 202.5,198.5,194.5,190.5,186.5,182.5,178.5,174.5, &
                 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)
      zxray = (/-10.7,-10.7,-10.7,-10.7,-10.7,-10.7, &
                -10.7,-10.7,-10.7,-10.7,-10.7,-10.7,-10.7,-10.7, &
                -10.7,-10.7,-14.7,-14.7,-14.7,-14.7,-14.7,-14.7, &
                -14.7,-14.7,-14.7,-14.7,-14.7,-14.7,-14.7,-14.7, &
                -14.7,-14.7, &
                130.1,130.1,130.1,130.1,130.1,130.1,130.1,130.1, &
                130.1,130.1,130.1,130.1,130.1,130.1,130.1,130.1, &
                132.6,132.6,132.6,132.6,132.6,132.6,132.6,132.6, &
                132.6,132.6,132.6,132.6,132.6,132.6,132.6,132.6, &
                0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)
      rxray = (/248.9,248.9,248.9,248.9,248.9,248.9,248.9, &
                248.9,248.9,248.9,248.9,248.9,248.9,248.9,248.9,248.9, &
                248.9,248.9,248.9,248.9,248.9,248.9,248.9,248.9, &
                248.9,248.9,248.9,248.9,248.9,248.9,248.9,248.9, &
                197.6,197.6,197.6,197.6,197.6,197.6,197.6,197.6, &
                197.6,197.6,197.6,197.6,197.6,197.6,197.6,197.6, &
                194.5,194.5,194.5,194.5,194.5,194.5,194.5,194.5, &
                194.5,194.5,194.5,194.5,194.5,194.5,194.5,194.5, &
                0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)
      rwstrip1(1)=1.33
      zwstrip1(1)=-1.363
      rwstrip1(2)=1.38
      zwstrip1(2)=-1.363
      rwstrip2(1)=1.4075
      zwstrip2(1)=-1.250
      rwstrip2(2)=1.4575
      zwstrip2(2)=-1.250
!
      if (ndimc.ne.ndim) then
        call errctrl_msg('pltout', &
                         'ndimc is not consistent with ndim')
        stop
      endif
!
      allocate(bfield(nw,nh),  &
         sivol(nw),voln(nw),bvoln(nw), &
         cvoln(nw),dvoln(nw),rscrap(nw),curscr(nw), &
         pbimf(nw),pmid(nw),pmidw(nw),bpmid(nw),cpmid(nw), &
         dpmid(nw),workj(nh),copyn(nwnh),copy1(nwnh),xbnow(npoint), &
         ybnow(npoint),rat(ndim))
!
#ifdef DEBUG_LEVEL1
      write (6,*) 'Enter PLTOUT'
#endif      
      if (kdot.gt.0.and.jtime.ne.kdot+1) return
      if (kdata.eq.4) then
         dismin=min(gapin(jtime),gapout(jtime),gaptop(jtime),gapbot(jtime))
         if (dismin.eq.gapin(jtime)) gapdis=sepin(jtime)
         if (dismin.eq.gapout(jtime)) gapdis=sepout(jtime)
         if (dismin.eq.gapbot(jtime)) gapdis=sepbot(jtime)
         if (dismin.eq.gaptop(jtime)) gapdis=septop(jtime)
         if (dismin.le.0.1_dp) dismin=gapdis
         if (idsep.eq.0) dsep(jtime)=dismin
         if (jtime.lt.ktime) return
      endif
      if (jerror(jtime).gt.0) return
!------------------------------------------------------------------------
!--   Set up flag to write plot file                                   --
!------------------------------------------------------------------------
      idotek = 1
      if (itek.ge.5) then
        mtek = itek - 10
        if (mtek.lt.0) then
           idotek=0
        elseif (mtek.eq.0) then
           idotek=1
        else
           mdotek = ktime/mtek
!          idotek = mod (jtime+ktime, mdotek)
           idotek = mod (jtime, mtek)
!          write (6,*) itek, ktime, mdotek, jtime, idotek
           if (istrpl.eq.0) idotek=0
        endif
      endif
!
      istrpl0: if (istrpl.le.0) then
      almin=xlmin
      almax=xlmax
      blmin=ylmin
      blmax=ylmax
!------------------------------------------------------------------------
!--   read in new SXR detectors geometry for shot > 80744              --
!--      and/or Toroidal x-ray geometry for shot > 91000               --
!------------------------------------------------------------------------
      if (ixray.gt.0) call getsxr(ishot,ixray)      
      if (ixray.gt.0) then       ! ixray could be changed by getsxr
         if (ixray.eq.1) then           ! poroidal xray only
            ixraystart = 1
            ixrayend = nangle
         elseif (ixray.eq.2) then      ! toroidal xray only
            ixraystart = nangle+1
            ixrayend = nangle+ntangle
         else                           ! both
            ixraystart = 1
            ixrayend = nangle+ntangle
         endif
         do i=ixraystart,ixrayend
            rxray(i)=rxray(i)/100.
            zxray(i)=zxray(i)/100.
            almin=min(almin,rxray(i))
            almax=max(almax,rxray(i))
            blmin=min(blmin,zxray(i))
            blmax=max(blmax,zxray(i))
         enddo
      endif
      if (iprobe.ne.0) then
         do i=1,nsilop
            almin=min(almin,rsi(i))
            almax=max(almax,rsi(i))
            blmin=min(blmin,zsi(i))
            blmax=max(blmax,zsi(i))
         enddo
         ibm=1
         iem=magpri67
         idelm=iprobe
         if (iprobe.ge.322) then
            ibm=magpri67+1
            iem=magpri67+magpri322+magprirdp
            idelm=iprobe-322+1
         endif
         idelm=iprobe
         if (iprobe.le.-322) then
            ibm=magpri67+magpri322+1
            iem=magpri67+magpri322+magprirdp
            idelm=1
         endif
         do i=ibm,iem
            almin=min(almin,xmp2(i))
            almax=max(almax,xmp2(i))
            blmin=min(blmin,ymp2(i))
            blmax=max(blmax,ymp2(i))
         enddo
      endif
      do i=1,nfcoil
         almin=min(almin,rf(i))
         almax=max(almax,rf(i))
         blmin=min(blmin,zf(i))
         blmax=max(blmax,zf(i))
      enddo
      if (iecoil.gt.0) then
         do i=1,necoil
            almin=min(almin,re(i))
            almax=max(almax,re(i))
            blmin=min(blmin,ze(i))
            blmax=max(blmax,ze(i))
         enddo
      endif
      if (kframe.gt.0) then
         almin=almin-0.05_dp
         almax=almax+0.05_dp
         blmin=blmin-0.05_dp
         blmax=blmax+0.05_dp
      endif
      elim=(blmax-blmin)/(almax-almin)
      uday=' '
      call date_and_time(uday, clocktime, zone, values)
      yll=7.2_dp
      if (iecoil.le.0) then
         if (klabel.ge.0) then
            yll=6.2_dp
         else
            yll=7.5_dp
         endif
      endif
      xll=yll/elim
!-----------------------------------------------------------------------
!     Open file PLTOUT.OUT to output plot parameters
!-----------------------------------------------------------------------
      if (itek.ge.5.and.idotek.eq.0) then
      iunit = 59
      if (kgraph.eq.0) then
        plotname='pltout.out'
      else
        if (ifitvs.eq.1.or.icutfp.eq.2) then
           write (6,*) 'itimeu = ',itimeu
        endif
        let2 = 'pl'
        call setfnmpl(itimeu,let2,ishot,itime,plotname)
        if (istore .eq. 1)  &
             plotname = store_dir(1:lstdir)//plotname
      endif
      if (m_write .eq. 1) then
        open(unit=iunit,file = plotname, status = 'old',iostat=ioerr)
        if (ioerr.eq.0) close(unit=iunit,status='delete')
        open (unit=iunit, file = plotname, status = 'new')
      elseif (m_write .eq. 0) then
        open (unit=iunit, file = plotname, form = 'unformatted', &
              status = 'old',iostat=ioerr)
        if (ioerr.eq.0) close(unit=iunit,status='delete')
        open (unit=iunit, file = plotname, form = 'unformatted', &
              status = 'new')
      endif
      endif
!
      endif istrpl0
      xmm=2.0
      if (itek.ge.5.and.idotek.eq.0) then
        if (kgraph.eq.1.and.istrpl.gt.0) then
!--------------------------------------------------------------------
!--       ITEK > 100, specific pltout.out name                     --
!--------------------------------------------------------------------
          let2 = 'pl'
          plotname = ' '
          call setfnmpl(itimeu,let2,ishot,itime,plotname)
          if (istore .eq. 1)  &
               plotname = store_dir(1:lstdir)//plotname
          iunit = 35
          close(unit=iunit)
          if (m_write .eq. 1) then
            open(unit=iunit,file = plotname, status = 'old',iostat=ioerr)
            if (ioerr.eq.0) close(unit=iunit,status='delete')
            open (unit=iunit, file = plotname, status = 'new')
          elseif (m_write .eq. 0) then
            open (unit=iunit, file = plotname, form = 'unformatted', &
                  status = 'old',iostat=ioerr)
            if (ioerr.eq.0) close(unit=iunit,status='delete')
            open (unit=iunit, file = plotname, form = 'unformatted', &
                  status = 'new')
          endif
        endif
      endif
!--------------------------------------------------------------------
!--   Set q-files dataname                                         --
!--   If ISTORE = 0 Then dataname is prefixed by qshot.time        --
!--   ELSE dataname is prefixed by /link/efit/qshot.time           --
!--   lprx = 13 prefix is qshot.time                               --
!--   lprx = 25  prefix is /link/efit/qshot.time                   --
!--------------------------------------------------------------------
      let = 'q'
      call setfnmq(let,ishot,itime,istore,dataname)
      lprx = 13
      if (istore .eq. 1) lprx = lstdir+lprx
!
      ! TODO: nfound never defined... using nplt instead?
!      do i=1,nfound-1
      do i=1,nplt-1
         ip1=i+1
         dli=sqrt((xout(ip1)-xout(i))**2+(yout(ip1)-yout(i))**2)
         pleng=pleng+dli
      enddo
      abar=100.*pleng/2./pi
!
      not_time_snap: if ((kdata.ne.4).and.(ivacum.eq.0)) then
      is_vacuum_1: if (ivacum.eq.0) then
      ibrdr = 1
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      if (klabel.ge.0) then
         xphy = 7.0
         yphy = 1.0
      else
         xphy = 3.0
         yphy = 1.6_dp
         bngle = 90.0
         bshft(1) = 6.10_dp
         bshft(2) = 0.95_dp
      endif
      if (kframe.eq.0) then
         xtitle = 'R(m)$'
         nxlen = 0
         ytitle = 'Z(m)$'
         nylen = 0
         xstp = almax - almin
         ystp = blmax - blmin
         igridx = -1
      else
         xtitle = 'R(m)$'
         nxlen = 100
         ytitle = 'Z(m)$'
         nylen = 100
         xstp = 0.0
         ystp = 0.0
      endif
      x11x = xll
      y11y = yll
      ptitle = '$'

      nn = 1
      nxy(nn) = nplt
      clearx(nn)='PINK'
      if (kthkcrv.eq.1) then
         thcrv(nn) = 0.015_dp
         ndshme(nn) = 1
      endif
      do i = 1, nplt
         xx(i,nn) = xplt(i)
         yy(i,nn) = yplt(i)
      enddo
      if (kplotp.ne.0) then
         do i=1,nplt
           xbnow(i)=xplt(i)
           ybnow(i)=yplt(i)
         enddo
        nbnow=nplt
      endif
      do i=1,nplt
         flxtra(i,1)=xplt(i)     ! for transfer to subroutine expand
         fpxtra(i,1)=yplt(i)
      enddo
      npltbdry=nplt

      if (kwripre.gt.0) then
        dataname=dataname(1:lprx)//'_surfb'
        open(unit=62,file=dataname,status='old',iostat=ioerr)
        if (ioerr.eq.0) close(unit=62,status='delete')
        open(unit=62,file=dataname,status='new')
        do i=1,nplt
          write (62,*) xplt(i),yplt(i)
        enddo
        close(unit=62)
      endif
      single_pass: if (ipass.le.1) then
      dismin=min(gapin(jtime),gapout(jtime),gaptop(jtime),gapbot(jtime))
      dis_min: if (dismin.gt.0.1_dp) then
      if (kwripre.gt.0) then
        dataname=dataname(1:lprx)//'_sep'
        open(unit=62,file=dataname,status='old',iostat=ioerr)
        if (ioerr.eq.0) close(unit=62,status='delete')
        open(unit=62,file=dataname,status='new')
      endif
      call surfac(psibry,psi,nw,nh,rgrid,zgrid,xplt,yplt,nplt, &
        npoint,drgrid,dzgrid,rgrid(1),rgrid(nw),zgrid(1), &
        zgrid(nh),n00,rmaxis,zmaxis,negcur,kerror,1)
      if (kerror.gt.0) return
      yminmm=zmin-0.001_dp
      ymaxmm=zmax+0.001_dp
      nn = nn + 1
      clearx(nn) = 'PINK'
      thcrv(nn) = 0.015_dp
      sclpc(nn) = 0.2_dp
      markme(nn) = 15
      ncnct(nn) = -1
      nxy(nn)   =  0
      ji = 0
      do i=1,nplt
         if ((yplt(i).ge.yminmm).and.(yplt(i).le.ymaxmm)) cycle
         nxy(nn) = 1 + nxy(nn)
         xx(nxy(nn) , nn) = xplt(i)
         yy(nxy(nn) , nn) = yplt(i)
         ji=ji+1
         flxtra(ji,2)=xplt(i)       ! for transfer to subroutine expand
         fpxtra(ji,2)=yplt(i)
         npltxpt=ji
         if (kwripre.gt.0) then
           call zlim(zeron,n11,n11,limitr,xlim,ylim,xplt(i),yplt(i),limfag)
           if (zeron.gt.0.001_dp) write (62,*) xplt(i),yplt(i)
         endif
      enddo
      if (kwripre.gt.0) close(unit=62)
      endif dis_min
      extra_surfs: if (nextra.ne.0) then
      nxcurv=0
      nexexx=2*iabs(nextra)*iabs(ixstrt)
      do i=1,nexexx
        if(npxtra(i).le.2) cycle
        nxcurv=nxcurv+1
        nn = nn + 1
        nxy(nn) = npxtra(i)
        clearx(nn) = 'GREE'
        kk = i
        do ii = 1, npxtra(i)
          xx(ii,nn) = xxtra(ii,i)
          yy(ii,nn) = yxtra(ii,i)
        enddo

        if (kwripre.gt.0) then
          if (i.lt.9) then
            write(iname,40023) i
            dataname=dataname(1:lprx)//'_fl'//iname
            open(unit=62,file=dataname,status='old',iostat=ioerr)
            if (ioerr.eq.0) close(unit=62,status='delete')
            open(unit=62,file=dataname,status='new')
            do iii=1,npxtra(i)
              call zlim(zeron,n11,n11,limitr,xlim,ylim,xxtra(iii,i),yxtra(iii,i),limfag)
              if (zeron.gt.0.001_dp) write (62,*) xxtra(iii,i),yxtra(iii,i)
            enddo
            close(unit=62)
          endif
        endif
      enddo
      endif extra_surfs
      endif single_pass
      fixed_bdry: if (nbdry.gt.0) then
      bdry: if (nbdry.le.20) then
      ! TODO: it is impossible to reach this code...
      nn = nn + 1
      markme(nn) = 16
      ncnct(nn) = -1
      nxy(nn) = nbdry
      do ii = 1, nbdry
         xx(ii,nn) = rbdry(ii)
         yy(ii,nn) = zbdry(ii)
      enddo
      else bdry
      bdryp: if (nbdryp <= 0) then
      nbdrz=nbdry+1
      rbdry(nbdrz)=rbdry(1)
      zbdry(nbdrz)=zbdry(1)
      nn = nn + 1
      nxy(nn) = nbdrz
      ndotme(nn) = 1
      do ii = 1, nbdrz
         xx(ii,nn) = rbdry(ii)
         yy(ii,nn) = zbdry(ii)
      enddo
      else bdryp
      nn = nn + 1
      nxy(nn) = nbdryp + 1
      ndotme(nn) = 1
      do ii = 1, nbdryp
         xx(ii,nn) = rbdry(ii)
         yy(ii,nn) = zbdry(ii)
      enddo
      ii=nbdryp+1
      xx(ii,nn) = rbdry(1)
      yy(ii,nn) = zbdry(1)
!  
      nn = nn + 1
      markme(nn) = 16
      ncnct(nn) = -1
      nxy(nn) = nbdry-nbdryp
      do ii = 1,nxy(nn) 
         iijj=ii+nbdryp
         xx(ii,nn) = rbdry(iijj)
         yy(ii,nn) = zbdry(iijj)
      enddo
      endif bdryp
      endif bdry
      endif fixed_bdry
      if (nsol > 0) then
         nn = nn + 1
         markme(nn) = 2
         ncnct(nn) = -1
         clearx(nn) = 'GREE'
         nxy(nn) = nsol
         do ii = 1,nsol
            xx(ii,nn) = rsol(ii)
            yy(ii,nn) = zsol(ii)
         enddo
      endif
      one_pass: if (ipass.le.1) then
      do i=1,nw-1
         worke(i)=real(i-1,dp)/(nw-1)
      enddo
      worke(nw)=1.
      call zpline(nw,worke,qpsi,bworkb,cworkb,dworkb)
      delsi=(psibry-simag)/(nsplot+1)
      nsplot0=nsplot
      idqplot=iqplot
      if (iqplot.ge.10) idqplot=iqplot-10
      if ((iqplot.gt.1).and.(psiq1.gt.0.)) then
         nsplot=nsplot0+1
         if(psiq1.gt.10.0)then
            psiq11=real(int(psiq1),dp)/1000.
            psiq12=psiq1-int(psiq1)
            nsplot=nsplot0+2
         else
            psiq11=psiq1
         endif
      endif
      do i=1,nsplot
         nn = nn + 1
         j=nsplot0-i+1
         if (i.gt.nsplot0.and.iqplot.ge.10) then
            siwant=psiq11
            psiq11=psiq12
         else
            siwant=simag+j*delsi
            if (kvtor.eq.0.or.kplotp.eq.0) ndshme(nn) = 1
         endif
         call surfac(siwant,psi,nw,nh,rgrid,zgrid,xplt,yplt,nplt, &
           npoint,drgrid,dzgrid,rmin,rmax,zmin, &
           zmax,n11,rmaxis,zmaxis,negcur,kerror)
         if (kerror.gt.0) return
         if (kthkcrv.gt.0) thcrv(nn) = 0.010_dp
         clearx(nn) = 'BLUE'
         nxy(nn) = nplt
         do ii = 1,nplt
            xx(ii,nn) = xplt(ii)
            yy(ii,nn) = yplt(ii)
         enddo
!----------------------------------------------------------------------
!--      regular surfaces                                            --
!----------------------------------------------------------------------
         if (kwripre.gt.0) then
            if (i.lt.nsplot.or.iqplot.lt.10) then
               write(iname,40023) i
               dataname=dataname(1:lprx)//'_surf'//iname
               open(unit=62,file=dataname,status='old',iostat=ioerr)
               if (ioerr.eq.0) close(unit=62,status='delete')
               open(unit=62,file=dataname,status='new')
               do iii=1,nplt
                  write (62,*) xplt(iii),yplt(iii)
               enddo
               close(unit=62)
            endif
!
            if (i.eq.nsplot.and.iqplot.ge.10) then
              write(iname,40023) i-nsplot0
              dataname=dataname(1:lprx)//'_surfq1'//iname
              open(unit=62,file=dataname,status='old',iostat=ioerr)
              if (ioerr.eq.0) close(unit=62,status='delete')
              open(unit=62,file=dataname,status='new')
              do iii=1,nplt
                  write (62,*) xplt(iii),yplt(iii)
              enddo
              close(unit=62)
            endif
         endif
!----------------------------------------------------------------------
!--      plot pressure surface for rotational equilibrium if requested
!----------------------------------------------------------------------
         if (kvtor.gt.0.and.kplotp.ne.0) then
           signp=1.
           if (kplotp.lt.0) signp=-1.
           nn = nn + 1
           xplt(nplt+1)=xplt(1)
           yplt(nplt+1)=yplt(1)
           do ii=1,nplt
             if (signp*xplt(ii).gt.signp*rmaxis) then
               if (yplt(ii)*yplt(ii+1).le.0.0) exit
             endif
           enddo
           call seva2d(bwx,lwx,bwy,lwy,cw,xplt(ii),yplt(ii), &
                    pds,ier,n111)
           spwant=pds(1)
           call surfac(spwant,presst,nw,nh,rgrid,zgrid,xplt,yplt, &
                       nplt,npoint,drgrid,dzgrid,rmin,rmax,zmin, &
                       zmax,n111,rmaxis,zmaxis,negcur,kerror)
           if (kerror.gt.0) return
           clearx(nn) = 'YELL'
           nxy(nn) = nplt
           ndshme(nn) = 1
           do ii = 1,nplt
              xx(ii,nn) = xplt(ii)
              yy(ii,nn) = yplt(ii)
           enddo
           if (kwripre.gt.0) then
             if (i.lt.nsplot.or.iqplot.lt.10) then
               write(iname,40023) i
               dataname=dataname(1:lprx)//'_surp'//iname
               open(unit=62,file=dataname,status='old',iostat=ioerr)
               if (ioerr.eq.0) close(unit=62,status='delete')
               open(unit=62,file=dataname,status='new')
               do iii=1,nplt
                 write (62,*) xplt(iii),yplt(iii)
               enddo
               close(unit=62)
             endif
           endif
         endif
!
         if (idqplot.gt.0) then
            if (iqplot.lt.10.or.i.lt.nsplot0) then
               sinow=j/real(nsplot0+1,dp)
               qnow=seval(nw,sinow,worke,qpsi,bworkb,cworkb,dworkb)
               msg = msg + 1
               note(msg) = 4
               anum(msg) = qnow
               iplce(msg) = idqplot
               xpos(msg) = xplt(nplt/6)
               ypos(msg) = yplt(nplt/6)
               ht(msg) = 0.10_dp
            endif
         endif
      enddo
40023 format (i1)
!-----------------------------------------------------------------------
!--   plot MSE locations                                              --
!-----------------------------------------------------------------------
      if (iplots.gt.0.and.kstark.gt.0) then
         nn = nn + 1
         markme(nn) = 3
         sclpc(nn) = 1.0
         innn    = 0
         do i=1,nstark
            if (abs(fwtgam(i)).gt.1.0e-06_dp) then
               innn=innn+1
               xx(innn,nn) = rrgam(jtime,i)
               yy(innn,nn) = zzgam(jtime,i)
            endif
         enddo
         nxy(nn) =  innn
         ncnct(nn)=-1
      endif  
      if (iplots.gt.0.and.mmbmsels.gt.0) then
         nn = nn + 1
         markme(nn) = 4
         sclpc(nn) = 1.0
         innn    = 0
         do i=1,nmsels
            if (abs(fwtbmselt(jtime,i)).gt.1.0e-06_dp) then
               innn=innn+1
               xx(innn,nn) = rrmselt(jtime,i)
               yy(innn,nn) = zzmselt(jtime,i)
            endif
         enddo
         nxy(nn) =  innn
         ncnct(nn)=-1
      endif
      kframe2: if (kframe.eq.2) then
         idelx=1
         idely=1
         mmw=nw
         mmh=nh
         if (mmw.gt.33) idelx=2
         if (mmh.gt.33) idely=2
         innn=0
         nn = nn + 1
         sclpc(nn) = 0.25_dp
         markme(nn) = 3
         ncnct(nn) = -1
         do i=1,nw,idelx
            do j=1,nh,idely
               kk=(i-1)*nh+j
               if (icutfp.eq.0) then
                  if (zero(kk)*www(kk).gt.0.001_dp) then
                     innn=innn+1
                     xx(innn,nn) = rgrid(i)
                     yy(innn,nn) = zgrid(j)
                  endif
               else
                  if (zero(kk)*xpsi(kk).gt.0.001_dp.and. &
                  xpsi(kk)*xpsimin.le.1.) then
                     innn=innn+ 1
                     xx(innn,nn) = rgrid(i)
                     yy(innn,nn) = zgrid(j)
                  endif
               endif
            enddo
         enddo
         nxy(nn) = innn
      endif kframe2
      plot_SOL_curr: if (icutfp.eq.2.and.kframe.ne.1) then
         innn=0
         nn = nn + 1
         sclpc(nn) = 0.3_dp
         markme(nn) = 15
         ncnct(nn) = -1
         rvsmin=rvs(1)
         rvsmax=rvs(1)
         zvsmin=zvs(1)
         zvsmax=zvs(1)
         do i=1,nvesel
            rvsmin=min(rvs(i),rvsmin)
            rvsmax=max(rvs(i),rvsmax)
            zvsmin=min(zvs(i),zvsmin)
            zvsmax=max(zvs(i),zvsmax)
         enddo
         call surfac(xpsialp,psi,nw,nh,rgrid,zgrid,xplt,yplt,nplt, &
           npoint,drgrid,dzgrid,rvsmin,rvsmax,zvsmin, &
           zvsmax,n00,rmaxis,zmaxis,negcur,kerror)
         if (kerror.gt.0) return
         if (iprobe.ne.90) then
           do i=1,nplt
             call zlim(zeron,n11,n11,limitr,xlim,ylim,xplt(i),yplt(i),limfag)
             if (zeron.gt.0.001_dp) then
                  innn=innn+1
                  xx(innn,nn) = xplt(i)
                  yy(innn,nn) = yplt(i)
             endif
           enddo
           nxy(nn) = innn
         endif
!         hgt = 0.10_dp
!         msg = msg + 1
!         note(msg) = 2
!         lmes(msg) = 'alpha =   '
!         imes(msg) = 10
!         xpos(msg) = rvsmax-0.5_dp
!         ypos(msg) = zvsmin
!         ht(msg) = hgt
!         call rlreal(alphafp,2,rvsmax,zvsmin)
!         msg = msg + 1
!         note(msg) = 4
!         anum(msg) = alphafp
!         iplce(msg) = 2
!         xpos(msg) = rvsmax
!         ypos(msg) = zvsmin
!         ht(msg) = 0.10_dp
         write(text,96251) alphafp
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = text
         imes(msg) = 25
         xpos(msg) = rvsmax-0.5_dp
         ypos(msg) = zvsmin
         ht(msg) = 0.10_dp
!
         sumif=0.0
         do i=1,nvesel
            sumif=sumif+vcurrt(i)
         enddo
         sumif=sumif/1000.
         write(text,96252) sumif
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = text
         imes(msg) = 25
         xpos(msg) = rvsmax-0.5_dp
         ypos(msg) = zvsmin-0.1_dp
         ht(msg) = 0.10_dp
         write (6,*) 'Ivt(ka) = ',sumif
!         msg = msg + 1
!         note(msg) = 2
!         lmes(msg) = 'Iv(ka)=   '
!         imes(msg) = 10
!         xpos(msg) = rvsmax-0.5_dp
!         ypos(msg) = zvsmin-0.1_dp
!         ht(msg) = hgt
!         call rlreal(sumif,0,rvsmax,zvsmin-0.1_dp)
!         msg = msg + 1
!         note(msg) = 4
!         anum(msg) = sumif
!         iplce(msg) = 0
!         xpos(msg) = rvsmax
!         ypos(msg) = zvsmin-0.1_dp
!         ht(msg) = 0.10_dp
!         msg = msg + 1
!         note(msg) = 2
!         lmes(msg) = 'Fp(kn)=   '
!         imes(msg) = 10
!         xpos(msg) = rvsmax-0.5_dp
!         ypos(msg) = zvsmin-0.2_dp
!         ht(msg) = hgt
!         sumif=fzpol/1000.
!         call rlreal(sumif,0,rvsmax,zvsmin-0.2_dp)
!         msg = msg + 1
!         note(msg) = 4
!         anum(msg) = sumif
!         iplce(msg) = 0
!         xpos(msg) = rvsmax
!         ypos(msg) = zvsmin-0.2_dp
!         ht(msg) = 0.10_dp
!         msg = msg + 1
!         note(msg) = 2
!         lmes(msg) = 'Ft(kn)=   '
!         imes(msg) = 10
!         xpos(msg) = rvsmax-0.5_dp
!         ypos(msg) = zvsmin-0.3_dp
!         ht(msg) = hgt
!         sumif=fztor/1000.
!         call rlreal(sumif,0,rvsmax,zvsmin-0.3_dp)
!         msg = msg + 1
!         note(msg) = 4
!         anum(msg) = sumif
!         iplce(msg) = 0
!         xpos(msg) = rvsmax
!         ypos(msg) = zvsmin-0.3_dp
!         ht(msg) = 0.10_dp
!
!         msg = msg + 1
!         note(msg) = 2
!         lmes(msg) = 'Ih(ka)=   '
!         imes(msg) = 10
!         xpos(msg) = rvsmax-0.5_dp
!         ypos(msg) = zvsmin-0.4_dp
!         ht(msg) = hgt
         sumif=tcurrp/1000.
         write(text,96253) sumif
         write (6,*) 'Ipt(ka) = ',sumif
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = text
         imes(msg) = 25
         xpos(msg) = rvsmax
         ypos(msg) = zvsmin-0.2_dp
         ht(msg) = 0.10_dp
!         call rlreal(sumif,0,rvsmax,zvsmin-0.4_dp)
!         msg = msg + 1
!         note(msg) = 4
!         anum(msg) = sumif
!         iplce(msg) = 0
!         xpos(msg) = rvsmax
!         ypos(msg) = zvsmin-0.4_dp
!         ht(msg) = 0.10_dp
!         msg = msg + 1
!         note(msg) = 2
!         lmes(msg) = 'Ig(ka)=   '
!         imes(msg) = 10
!         xpos(msg) = rvsmax-0.5_dp
!         ypos(msg) = zvsmin-0.5_dp
!         ht(msg) = hgt
         sumif=fpolvs/1000.
!         call rlreal(sumif,0,rvsmax,zvsmin-0.5_dp)
         write(text,96254) sumif
         write (6,*) 'Ivp(ka) = ',sumif
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = text
         imes(msg) = 25
         xpos(msg) = rvsmax
         ypos(msg) = zvsmin-0.3_dp
         ht(msg) = 0.10_dp
      endif plot_SOL_curr
      endif one_pass
      plot_limiter: if ((iplim.le.0).or.((iplim.ne.10) &
          .and.((limitr.ne.(limid+1)).or.(xlim(limitr+1).eq.0.0)))) then
      nn = nn + 1
      nxy(nn) = limitr
      do ii = 1, limitr
         xx(ii,nn) = xlim(ii)
         yy(ii,nn) = ylim(ii)
      enddo
      if (kwripre.gt.0) then
         dataname=dataname(1:lprx)//'_lim'
         open(unit=62,file=dataname,status='old',iostat=ioerr)
         if (ioerr.eq.0) close(unit=62,status='delete')
         open(unit=62,file=dataname,status='new')
         do i=1,limitr
            write (62,*) xlim(i),ylim(i)
         enddo
         close(unit=62)
      endif
!--------------------------------------------------------------------
!--   Plot W strips                                                --
!--------------------------------------------------------------------
      nn = nn + 1
      nxy(nn) = 2
      clearx(nn)='PINK'
      thcrv(nn) = 0.030_dp
      do ii = 1, 2
         xx(ii,nn) = rwstrip1(ii)
         yy(ii,nn) = zwstrip1(ii)
      enddo
!
      nn = nn + 1
      nxy(nn) = 2
      clearx(nn)='PINK'
      thcrv(nn) = 0.030_dp
      do ii = 1, 2
         xx(ii,nn) = rwstrip2(ii)
         yy(ii,nn) = zwstrip2(ii)
      enddo
      else plot_limiter
!-------------------------------------------------------------------
!--   plot vessel                                                 --
!-------------------------------------------------------------------
      lim_opt: if (iplim.eq.10) then
      call pltcol(nvesel,rvs,zvs,wvs,hvs,avs,avs2,n00, &
      nn, xx, yy, nxy, msg, note, inum, xpos, ypos, ht, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps)
      nn = nn + 1
      thcrv(nn) = 0.02_dp
      nxy(nn) = limbot-limup + 1
      do ii = limup, limup+nxy(nn)-1
         xx(ii,nn) = xlim(ii)
         yy(ii,nn) = ylim(ii)
      enddo
      xplt(1)=xlim((limup+limbot)/2)
      yplt(1)=ylim((limup+limbot)/2)
      xplt(2)=rvs(12)-wvs(12)/2.
      yplt(2)=yplt(1)
      nn = nn + 1
      nxy(nn) = 2
      do ii = 1, nxy(nn)
         xx(ii,nn) = xplt(ii)
         yy(ii,nn) = yplt(ii)
      enddo
      else lim_opt
!
      iupper=0
      ilower=0
      do i=1,limitr
         if ((xlim(i).eq.xlim(limitr+1)).and.(ylim(i).eq.ylim(limitr+1))) &
           iupper=i
         if ((xlim(i).eq.xlim(limitr+12)).and.(ylim(i).eq. &
           ylim(limitr+12))) ilower=i
      enddo
      nn = nn + 1
      nxy(nn) = iupper
      do ii = 1, nxy(nn)
         xx(ii,nn) = xlim(ii)
         yy(ii,nn) = ylim(ii)
      enddo
      nn = nn + 1
      nxy(nn) = 12
      do ii = limitr+1, limitr+nxy(nn)
         xx(ii,nn) = xlim(ii)
         yy(ii,nn) = ylim(ii)
      enddo
      nn = nn + 1
      nxy(nn) = limitr-ilower+1
      do ii = ilower, ilower+nxy(nn)-1
         xx(ii,nn) = xlim(ii)
         yy(ii,nn) = ylim(ii)
      enddo
      open(unit=lfile,status='old',file=filimt)
      read (lfile,14980) npts
      read (lfile,15000) (workc(1),workd(1),i=1,npts) &
                 ,(workc(1),workd(1),i=npts+2,npts+13)
      read (lfile,14980) npts
      read (lfile,15000) (workc(i),workd(i),i=1,npts)
      nn = nn + 1
      nxy(nn) = npts
      do ii = 1, nxy(nn)
         xx(ii,nn) = workc(ii)
         yy(ii,nn) = workd(ii)
      enddo
!
      read (lfile,14980) npts
      read (lfile,15000) (workc(i),workd(i),i=1,npts)
!      call shade(workc,workd,npts,30.,0.1_dp,1,0,0)
      nshd = nshd + 1
      nsxy(nshd) = npts
      sangle(nshd) = 30.0
      sgap(nshd) = 0.1_dp
      ngaps(nshd) = 1
      do ii = 1, nsxy(nshd)
         sxx(ii,nshd) = workc(ii)
         syy(ii,nshd) = workd(ii)
      enddo
      nn = nn + 1
      nxy(nn) = 2
      clearx(nn) = 'GREE'
      do ii = 4, 4+nxy(nn)-1
         xx(ii,nn) = xlim(ii)
         yy(ii,nn) = ylim(ii)
      enddo
      nn = nn + 1
      nxy(nn) = 2
      clearx(nn) = 'GREE'
      do ii = 8, 8+nxy(nn)-1
         xx(ii,nn) = xlim(ii)
         yy(ii,nn) = ylim(ii)
      enddo

      read (lfile,14980) npts
      read (lfile,15000) (workc(i),workd(i),i=1,npts)
!      call shade(workc,workd,npts,30.,0.1_dp,1,0,0)
      nshd = nshd + 1
      nsxy(nshd) = npts
      sangle(nshd) = 30.0
      sgap(nshd) = 0.1_dp
      ngaps(nshd) = 1
      do ii = 1, nsxy(nshd)
         sxx(ii,nshd) = workc(ii)
         syy(ii,nshd) = workd(ii)
      enddo
      nn = nn + 1
      nxy(nn) = 2
      clearx(nn) = 'GREE'
      do ii = 5, 6
         xx(ii,nn) = workc(ii)
         yy(ii,nn) = workd(ii)
      enddo
      nn = nn + 1
      nxy(nn) = 2
      clearx(nn) = 'GREE'
      do ii = 10, 11
         xx(ii,nn) = workc(ii)
         yy(ii,nn) = workd(ii)
      enddo
!
      read (lfile,14980) npts
      read (lfile,15000) (workc(i),workd(i),i=1,npts)
      nn = nn + 1
      nxy(nn) = 3
      do ii =1, nxy(nn)
         xx(ii,nn) = workc(ii)
         yy(ii,nn) = workd(ii)
      enddo
!      call shade(workc,workd,npts,90.,0.01_dp,1,0,0)
      nshd = nshd + 1
      nsxy(nshd) = npts
      sangle(nshd) = 90.0
      sgap(nshd) = 0.01_dp
      ngaps(nshd) = 1
      do ii = 1, nsxy(nshd)
         sxx(ii,nshd) = workc(ii)
         syy(ii,nshd) = workd(ii)
      enddo
!
      do j=1,2
         read (lfile,14980) npts
         read (lfile,15000) (workc(i),workd(i),i=1,npts)
         nn = nn + 1
         nxy(nn) = 2
         do ii = 1, nxy(nn)
            xx(ii,nn) = workc(ii)
            yy(ii,nn) = workd(ii)
         enddo
!         call shade(workc,workd,npts,45.,0.05_dp,1,0,0)
         nshd = nshd + 1
         nsxy(nshd) = npts
         sangle(nshd) = 45.0
         sgap(nshd) = 0.05_dp
         ngaps(nshd) = 1
         do ii = 1, nsxy(nshd)
            sxx(ii,nshd) = workc(ii)
            syy(ii,nshd) = workd(ii)
         enddo
!         call shade(workc,workd,npts,135.,0.05_dp,1,0,0)
         nshd = nshd + 1
         nsxy(nshd) = npts
         sangle(nshd) = 135.0
         sgap(nshd) = 0.05_dp
         ngaps(nshd) = 1
         do ii = 1, nsxy(nshd)
            sxx(ii,nshd) = workc(ii)
            syy(ii,nshd) = workd(ii)
         enddo
      enddo
!
      read (lfile,14980) npts
      read (lfile,15000) (workc(i),workd(i),i=1,npts)
      nn = nn + 1
      nxy(nn) = 2
      do ii = 1, nxy(nn)
         xx(ii,nn) = workc(ii)
         yy(ii,nn) = workd(ii)
      enddo
!     call shade(workc,workd,npts,90.,0.01_dp,1,0,0)
      nshd = nshd + 1
      nsxy(nshd) = npts
      sangle(nshd) = 90.0
      sgap(nshd) = 0.01_dp
      ngaps(nshd) = 1
      do ii = 1, nsxy(nshd)
         sxx(ii,nshd) = workc(ii)
         syy(ii,nshd) = workd(ii)
      enddo
      close(unit=lfile)
      endif lim_opt
      endif plot_limiter
!--------------------------------------------------------------------
!--   vertical feedback ?                                          --
!--------------------------------------------------------------------
      if (isetfb.ne.0) then
         kct(1)=kct1
         kct(2)=kct2
         kct(3)=kct3
         kct(4)=kct4
         do ms=1,4
            if(kct(ms).ne.0)then
               do i=1,nw
                  do j=1,nh
                     kk=(i-1)*nh+j
                     if(kk.eq.kct(ms))then
                        nn = nn + 1
                        markme(nn) = 11
                        sclpc(nn) = 0.75_dp
                        nxy(nn) = 1
                        ncnct(nn) = -1
                        xx(1,nn) = rgrid(i)
                        yy(1,nn) = zgrid(j)
                     endif
                  enddo
               enddo
            endif
         enddo
      endif
      plot_probes: if (iprobe.ne.0) then
      nn = nn + 1
      markme(nn) = 16
      sclpc(nn) = 0.5_dp
      if (kthkcrv.eq.1) sclpc(nn) = 1.2_dp
      nxy(nn) = nsilop
      ncnct(nn) = -1
      do ii = 1, nxy(nn)
         xx(ii,nn) = rsi(ii)
         yy(ii,nn) = zsi(ii)
      enddo
      if (klabel.eq.0) then
         do i=1,nsilop,idelm
!           call rlint(i,rsi(i),zsi(i))
            msg = msg + 1
            note(msg) = 6
            inum(msg) = i
            xpos(msg) = rsi(i)
            ypos(msg) = zsi(i)
            ht(msg) = 0.04_dp
         enddo
      else
         hgt = 0.13_dp
!        if (klabel.lt.0) call hgt = 0.16_dp
         if (klabel.lt.0) hgt = 0.16_dp
         rf16=rf(16)-0.10_dp
         zf16=zf(16)+0.30_dp
         msg = msg + 1
         note(msg) = 2
         lmes(msg) = 'flux loops $'
         imes(msg) = 100
         xpos(msg) = rf16
         ypos(msg) = zf16
         ht(msg) = hgt
!        call rlvec(rf16,zf16,rsi(16),zsi(16),1)
         nvec = nvec + 1
         xfm(nvec) = rf16
         yfm(nvec) = zf16
         xto(nvec) = rsi(16)
         yto(nvec) = zsi(16)
         ivec(nvec) = 1
      endif
      nn = nn + 1
      markme(nn) = 15
      sclpc(nn) = 0.5_dp
      if (kthkcrv.eq.1) sclpc(nn) = 1.2_dp
      ncnct(nn) = -1
      nxy(nn) = iem-ibm+1
      do ii = ibm, ibm+nxy(nn)-1
         xx(ii,nn) = xmp2(ii)
         yy(ii,nn) = ymp2(ii)
      enddo
      klabel0: if (klabel.eq.0) then
         do i=ibm,iem,idelm
!           call rlint(i,xmp2(i),ymp2(i))
            msg = msg + 1
            note(msg) = 6
            inum(msg) = i
            xpos(msg) = xmp2(i)
            ypos(msg) = ymp2(i)
            ht(msg) = 0.04_dp
         enddo
      else klabel0
         hgt = 0.13_dp
         if (klabel.lt.0) hgt = 0.16_dp
         rrff=(rf(8)+rf(9))*0.5_dp
         zzff=(zf(8)+zf(9))*0.5_dp
         msg = msg + 1
         note(msg) = 2
         lmes(msg) = 'magnetic probes $'
         imes(msg) = 100
         xpos(msg) = rrff
         ypos(msg) = zzff
         ht(msg) = hgt
!        call rlvec(rrff,zzff,xmp2(8),ymp2(8),1)
         nvec = nvec + 1
         xfm(nvec) = rrff
         yfm(nvec) = zzff
         xto(nvec) = xmp2(8)
         yto(nvec) = ymp2(8)
         ivec(nvec) = 1
         ifcoil=0
         workd(1)=0.0
         workd(2)=0.0
         workc(1)=almin+0.15_dp
         workc(2)=almax-0.25_dp
         nn = nn + 1
         ncdtme(nn) = 1
         thcrv(nn) = 0.012_dp
         nxy(nn) = 2
         do ii = 1, nxy(nn)
            xx(ii,nn) = workc(ii)
            yy(ii,nn) = workd(ii)
         enddo
!        call reset('chndot')
         workc(3)=almax
         workd(3)=0.20_dp
         msg = msg + 1
         note(msg) = 2
         lmes(msg) = 'CER, $'
         imes(msg) = 100
         xpos(msg) = workc(3)
         ypos(msg) = workd(3)
         ht(msg) = hgt
!        call mixalf('INSTRU')
         mxalf(msg) = 'INSTRU'
         workd(3)=0.10_dp
         msg = msg + 1
         note(msg) = 2
         lmes(msg) = 'CO(LH.8)2(LXHX) chord $'
         imes(msg) = 100
         xpos(msg) = workc(3)
         ypos(msg) = workd(3)
         ht(msg) = hgt
!        call rlvec(workc(3),workd(3),workc(2),workd(2),1)
         nvec = nvec + 1
         xfm(nvec) = workc(3)
         yfm(nvec) = workd(3)
         xto(nvec) = workc(2)
         yto(nvec) = workd(2)
         ivec(nvec) = 1
         workd(1)=zuperts(jtime)*0.012_dp
         workd(2)=zlowerts*0.01_dp-0.400_dp
         workc(1)=rmajts
         workc(2)=rmajts
         nn = nn + 1
         ncdtme(nn) = 1
         thcrv(nn) = 0.012_dp
         nxy(nn) = 2
         do ii = 1, nxy(nn)
            xx(ii,nn) = workc(ii)
            yy(ii,nn) = workd(ii)
         enddo
         workc(3)=rmajts+0.20_dp
         workd(3)=zuperts(jtime)*0.0150_dp
!        call rlvec(workc(3),workd(3),workc(1),workd(1),1)
         nvec = nvec + 1
         xfm(nvec) = workc(3)
         yfm(nvec) = workd(3)
         xto(nvec) = workc(1)
         yto(nvec) = workd(1)
         ivec(nvec) = 1
         msg = msg + 1
         note(msg) = 2
         lmes(msg) = 'Thomson, $'
         imes(msg) = 100
         xpos(msg) = workc(3)
         ypos(msg) = workd(3)
         ht(msg) = hgt
         workd(3)=workd(3)-0.1_dp
         msg = msg + 1
         note(msg) = 2
         lmes(msg) = 'CO(LH.8)2(LXHX) chord $'
         imes(msg) = 100
         xpos(msg) = workc(3)
         ypos(msg) = workd(3)
         ht(msg) = hgt
         workc(1)=chordv(1)
         workc(2)=chordv(1)
         nn = nn + 1
         ncdtme(nn) = 1
         thcrv(nn) = 0.012_dp
         nxy(nn) = 2
         do ii = 1, nxy(nn)
            xx(ii,nn) = workc(ii)
            yy(ii,nn) = workd(ii)
         enddo
         workc(3)=workc(1)-0.4_dp
         workd(3)=workd(1)+0.59_dp
!        call rlvec(workc(3),workd(3),workc(1),workd(1),1)
         nvec = nvec + 1
         xfm(nvec) = workc(3)
         yfm(nvec) = workd(3)
         xto(nvec) = workc(1)
         yto(nvec) = workd(1)
         ivec(nvec) = 1
         msg = msg + 1
         note(msg) = 2
         lmes(msg) = 'CO(LH.8)2(LXHX) chord $'
         imes(msg) = 100
         xpos(msg) = workc(3)
         ypos(msg) = workd(3)
         ht(msg) = hgt
         workc(1)=chordv(3)
         workc(2)=chordv(3)
         workd(1)=workd(1)-0.2_dp
         workd(2)=workd(2)+0.2_dp
         nn = nn + 1
         ncdtme(nn) = 1
         thcrv(nn) = 0.012_dp
         nxy(nn) = 2
         do ii = 1, nxy(nn)
            xx(ii,nn) = workc(ii)
            yy(ii,nn) = workd(ii)
         enddo
         workc(3)=workc(2)+0.2_dp
         workd(3)=workd(2)-0.3_dp
!        call rlvec(workc(3),workd(3),workc(2),workd(2),1)
         nvec = nvec + 1
         xfm(nvec) = workc(3)
         yfm(nvec) = workd(3)
         xto(nvec) = workc(2)
         yto(nvec) = workd(2)
         ivec(nvec) = 1
         msg = msg + 1
         note(msg) = 2
         lmes(msg) = 'CO(LH.8)2(LXHX) chord $'
         imes(msg) = 100
         xpos(msg) = workc(3)
         ypos(msg) = workd(3)
         ht(msg) = hgt
!-----------------------------------------------------------------------
!--      motional Stark                                               --
!-----------------------------------------------------------------------
         workc(1)=1.90_dp
         workd(1)=0.0
         nn = nn + 1
         markme(nn) = 4
         sclpc(nn) = 2.0
         thcrv(nn) = 0.012_dp
         ncnct(nn) = -1
         nxy(nn) = 1
         xx(1,nn) = workc(1)
         yy(1,nn) = workd(1)
         workc(2)=workc(1)-0.32_dp
         workd(2)=workd(1)+0.17_dp
         msg = msg + 1
         note(msg) = 2
         lmes(msg) = 'Stark $'
         imes(msg) = 100
         xpos(msg) = workc(2)
         ypos(msg) = workd(2)
         ht(msg) = hgt
         workc(3)=workc(2)+0.08_dp
         workd(3)=workd(2)-0.03_dp
         workc(4)=workc(1)-0.02_dp
         workd(4)=workd(1)+0.02_dp
!        call rlvec(workc(3),workd(3),workc(4),workd(4),1)
         nvec = nvec + 1
         xfm(nvec) = workc(3)
         yfm(nvec) = workd(3)
         xto(nvec) = workc(4)
         yto(nvec) = workd(4)
         ivec(nvec) = 1
         if (klabel.lt.0) then
           if (klabel.gt.-100) then
             klabela=iabs(klabel)
             write(text, 9950) klabela
             xabs=1.0
             yabs=-0.5_dp
             dyabs = 0.22_dp
             msg = msg + 1
             note(msg) = 1
             lmes(msg) = text
             imes(msg) = 27
             xpos(msg) = xabs
             ypos(msg) = yabs
             yabs = yabs - dyabs
             ht(msg) = 0.16_dp
           endif
           ! call physor(7.0,1.0)
           xphy = 7.0
           yphy = 1.0
           yll=6.2_dp
           xll=yll/elim
           xtitle = 'R(m)$'
           ytitle = 'Z(m)$'
           nxlen = 0
           nylen = 0
           xlen = xll
           ylen = yll
           ! call title('$',-100,'R(m)$',0,'Z(m)$',0,xll,yll)
         endif
      endif klabel0
      endif plot_probes
!
      SXR: if (ixray.ne.0) then
      do i=ixraystart,ixrayend
         xplt(1)=rxray(i)
         yplt(1)=zxray(i)
         theta=xangle(i)*radeg
         mdot = 1
         iddsxr=0
         if (i.le.nangle/2) then
            do ii=1,10
               if (i.eq.ksxr0(ii)) then
                  mdot = 0
                  iddsxr=ii
               endif
            enddo
         else
            j=i-nangle/2
            do ii=1,10
               if (j.eq.ksxr2(ii)) then
                  mdot = 0
                  iddsxr=ii
               endif
            enddo
         endif
         do k=2,npoint
            xplt(k)=k*0.08_dp*cos(theta)+rxray(i)
            yplt(k)=k*0.08_dp*sin(theta)+zxray(i)
            if ((xplt(k).lt.almin).or.(xplt(k).gt.almax)) exit
            if ((yplt(k).lt.blmin).or.(yplt(k).gt.blmax)) exit
         enddo
         if (iddsxr.gt.0.or.idosxr.eq.1) then
            nn = nn + 1
            ndotme(nn) = mdot
            nxy(nn) = k
            do ii = 1, nxy(nn)
               xx(ii,nn) = xplt(ii)
               yy(ii,nn) = yplt(ii)
            enddo
         endif
!----------------------------------------------------------------------
!--      write out SXR files                                         --
!----------------------------------------------------------------------
         if (kwripre.gt.0.and.iddsxr.gt.0) then
            write(iname,40023) iddsxr
            if (i.le.nangle/2) then
               dataname=dataname(1:lprx)//'_sxr0'//iname
            else
               dataname=dataname(1:lprx)//'_sxr2'//iname
            endif
            open(unit=62,file=dataname,status='old',iostat=ioerr)
            if (ioerr.eq.0) close(unit=62,status='delete')
            open(unit=62,file=dataname,status='new')
            do iii=1,k
               call zlim(zeron,n11,n11,limitr,xlim,ylim,xplt(iii), &
               yplt(iii),limfag)
               if (zeron.gt.0.001_dp) write (62,*) xplt(iii),yplt(iii)
            enddo
            close(unit=62)
         endif
      enddo
!
      endif SXR
      if (ifcoil.gt.0) then
        call pltcol(nfcoil,rf,zf,wf,hf,af,af2,n11, &
          nn, xx, yy, nxy, msg, note, inum, xpos, ypos, ht, &
          nshd, sxx, syy, nsxy, sangle, sgap, ngaps)
      endif
      if (iecoil.gt.0) then
        do i=1,necoil
          ae(i)=0.0
          ae2(i)=0.0
          ishade=-ecid(i)
          call pltcol(n11,re(i),ze(i),we(i),he(i),ae(i),ae2(i),ishade, &
            nn, xx, yy, nxy, msg, note, inum, xpos, ypos, ht, &
            nshd, sxx, syy, nsxy, sangle, sgap, ngaps)
          if (ecid(i).eq.1) clearx(nn)='PINK'
          if (ecid(i).eq.2) clearx(nn)='FORE'
          if (ecid(i).eq.3) clearx(nn)='GREE'
          if (ecid(i).eq.4) clearx(nn)='CYAN'
          if (ecid(i).eq.5) clearx(nn)='YELL'
          if (ecid(i).eq.6) clearx(nn)='ORAN'
        enddo
      endif
!
      pasmac=ipmhd(jtime)/1000.
      pasmap=ipmeas(jtime)/1000.
!
      klabel_1: if (klabel.ge.0) then
      xabs=-3.0
      yabs=7.0
      dyabs = 0.22_dp
      if ((kdata.eq.1).or.(kdata.eq.2)) then
         write(text, 8948) ifname(jtime)
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = text
         imes(msg) = 25
         xpos(msg) = xabs
         ypos(msg) = yabs
         yabs = yabs - dyabs
         ht(msg) = 0.14_dp
      endif
      xabs=-6.5_dp
      yabs=7.0
      dyabs = 0.22_dp
      write(text,8950) trim(ch1),trim(ch2),efitvers
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 27
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,8960) uday
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,9000) ishot
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,9020) itime,itimeu
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,9040) chisq(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,9060) rout(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,9080) zout(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,9100) aminor(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,9120) elong(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,9140) utri(jtime),ltri(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,9155) indent(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp

      vm3=volume(jtime)/1.e6_dp
      am2=area(jtime)/1.e4_dp
      write(text,9160) vm3,am2
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,9165) wmhd(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,9180) betat(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,9200) betap(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,9150) betatn,pasman
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,9220) li(jtime),li3(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp

      if (ipass.le.1) then
        write(text,9230) errorm,nitera
        msg = msg + 1
        note(msg) = 1
        lmes(msg) = text
        imes(msg) = 26
        xpos(msg) = xabs
        ypos(msg) = yabs
        yabs = yabs - dyabs
        ht(msg) = 0.14_dp
        write(text,9240) delerb,delerr
        msg = msg + 1
        note(msg) = 1
        lmes(msg) = text
        imes(msg) = 29
        xpos(msg) = xabs
        ypos(msg) = yabs
        yabs = yabs - dyabs
        ht(msg) = 0.14_dp
      endif
      fpolss=fpolvs/1000.
      write(text,9262) dminux(jtime),dminlx(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
!
      edgej=ipmhd(jtime)/area(jtime)*1.e4_dp
      if (icurrt.eq.2.or.icurrt.eq.5) then
        if (rseps(1,jtime).gt.0.0) then
          rsepcm=rseps(1,jtime)/100.
          cjsep=rsepcm*ppcurr(x111,kppcur) &
                  +fpcurr(x111,kffcur)/rsepcm
          cjsep=cjsep/darea
          if (abs(edgej).gt.1.e-6_dp) then
            if (kzeroj.gt.0.and.rzeroj(1).lt.0.0) then
              edgej=cjsep/edgej
            else
              edgej=cjor(nw)/edgej
            endif
          endif
        else
          edgej=cjor(nw)/edgej
        endif
      endif
!
      if (icurrt.eq.4.or.icurrt.eq.1) edgej=cjor(nw)/edgej
      write (text,9265) cjor0(jtime),edgej
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,9268) qout(jtime),q95(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp

      dismin=min(gapin(jtime),gapout(jtime),gaptop(jtime),gapbot(jtime))
      if (dismin.eq.gapin(jtime)) gapdis=sepin(jtime)
      if (dismin.eq.gapout(jtime)) gapdis=sepout(jtime)
      if (dismin.eq.gapbot(jtime)) gapdis=sepbot(jtime)
      if (dismin.eq.gaptop(jtime)) gapdis=septop(jtime)
      if (dismin.le.0.1_dp) dismin=gapdis
      if (idsep.eq.0) dsep(jtime)=dismin
      write(text,9272) dsep(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,9280) rm(jtime),rcurrt(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,9290) zm(jtime),zcurrt(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,9542) betapd(jtime),betapw(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,9544) betatd(jtime),betatw(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,9546) wdia(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,9547) wplasw(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,9300)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,9390) pasmap
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,9400) bcentr(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      iece=kece+kecebz
      write(text,9320) ipsi(jtime),iece
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,9340) imag2(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,9380) iplasm(jtime),ifc(jtime),iec(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      if (mmbmsels.eq.0) then
        write(text,9385) idlopc(jtime),kmtark,klibim
      else
        write(text,99385) idlopc(jtime),kmtark,mmbmsels
      endif
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      if(nextra.ne.0) then
        write(text,9395)(scraps(i),i=1,iabs(nextra))
        msg = msg + 1
        note(msg) = 1
        lmes(msg) = text
        imes(msg) = 40
        xpos(msg) = xabs
        ypos(msg) = yabs
        yabs = yabs - dyabs
        ht(msg) = 0.14_dp
      endif
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = ' '
      imes(msg) = 1
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
!
      if (icurrt.eq.1) then
         xabs=-3.5_dp
         yabs=0.8_dp
         write(text,19603) srm,salpha
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = text
         imes(msg) = 27
         xpos(msg) = xabs
         ypos(msg) = yabs
         yabs = yabs - dyabs
         ht(msg) = 0.14_dp
         write(text,19605) sbeta,sbetaw
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = text
         imes(msg) = 27
         xpos(msg) = xabs
         ypos(msg) = yabs
         yabs = yabs - dyabs
         ht(msg) = 0.14_dp
         saaacm=saaa*100.
         write(text,19610) saaacm
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = text
         imes(msg) = 27
         xpos(msg) = xabs
         ypos(msg) = yabs
         yabs = yabs - dyabs
         ht(msg) = 0.14_dp
      endif
!
      if (icurrt.eq.2.or.icurrt.eq.5) then
         xabs=-3.5_dp
         yabs=0.8_dp
         write(text,9603) icprof,ifitvs
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = text
         imes(msg) = 25
         xpos(msg) = xabs
         ypos(msg) = yabs
         yabs = yabs - dyabs
         ht(msg) = 0.14_dp
         write(text,9605) kffcur,kppcur,keecur
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = text
         imes(msg) = 25
         xpos(msg) = xabs
         ypos(msg) = yabs
         yabs = yabs - dyabs
         ht(msg) = 0.14_dp
         write(text,9610) fcurbd,pcurbd,ecurbd
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = text
         imes(msg) = 25
         xpos(msg) = xabs
         ypos(msg) = yabs
         yabs = yabs - dyabs
         ht(msg) = 0.14_dp
         if (kvtor.ge.1.and.kvtor.le.3) then
           write(text,9612) kwwcur,wcurbd
           msg = msg + 1
           note(msg) = 1
           lmes(msg) = text
           imes(msg) = 25
           xpos(msg) = xabs
           ypos(msg) = yabs
           yabs = yabs - dyabs
           ht(msg) = 0.14_dp
         elseif (kvtor.eq.11) then
           write(text,9614) kvtor
           msg = msg + 1
           note(msg) = 1
           lmes(msg) = text
           imes(msg) = 25
           xpos(msg) = xabs
           ypos(msg) = yabs
           yabs = yabs - dyabs
           ht(msg) = 0.14_dp
         endif
!
         wwtbp=0.0
         wwtqa=0.0
         if (abs(fwtbp).gt.1.e-8_dp) wwtbp=abs(fwtbp/fwtbp)
         if (abs(fwtqa).gt.1.e-8_dp) wwtqa=abs(fwtqa/fwtqa)
         write(text,9615) wwtbp,wwtqa
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = text
         imes(msg) = 25
         xpos(msg) = xabs
         ypos(msg) = yabs
         yabs = yabs - dyabs
         ht(msg) = 0.14_dp

         zelipcm=zelip*100.
         write(text,9635) zelipcm
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = text
         imes(msg) = 25
         xpos(msg) = xabs
         ypos(msg) = yabs
         yabs = yabs - dyabs
         ht(msg) = 0.14_dp

         if (abs(relax-1.0).gt.1.e-4_dp) then
            write(text,9637) relax
            msg = msg + 1
            note(msg) = 1
            lmes(msg) = text
            imes(msg) = 25
            xpos(msg) = xabs
            ypos(msg) = yabs
            yabs = yabs - dyabs
            ht(msg) = 0.14_dp
         endif
         if (kcalpa.gt.0.or.kcgama.gt.0) then
            write(text,9620) kcgama,kcalpa
            msg = msg + 1
            note(msg) = 1
            lmes(msg) = text
            imes(msg) = 25
            xpos(msg) = xabs
            ypos(msg) = yabs
            yabs = yabs - dyabs
         endif
         if (kcomega.gt.0) then
            write(text,9623) kcomega
            msg = msg + 1
            note(msg) = 1
            lmes(msg) = text
            imes(msg) = 25
            xpos(msg) = xabs
            ypos(msg) = yabs
            yabs = yabs - dyabs
         endif
         if (kzeroj.gt.0) then
            write(text,9639) rzeroj(1)
            msg = msg + 1
            note(msg) = 1
            lmes(msg) = text
            imes(msg) = 25
            xpos(msg) = xabs
            ypos(msg) = yabs
            yabs = yabs - dyabs
            ht(msg) = 0.14_dp
         endif
         if (icutfp.eq.2) then
            if (alphafp.ge.0.0) alphamu=alphafp
            write(text,9625) alphamu
            msg = msg + 1
            note(msg) = 1
            lmes(msg) = text
            imes(msg) = 25
            xpos(msg) = xabs
            ypos(msg) = yabs
            yabs = yabs - dyabs
            ht(msg) = 0.14_dp
         endif
         if (kprfit.gt.0) then
            write(text, 9630) kprfit
            msg = msg + 1
            note(msg) = 1
            lmes(msg) = text
            imes(msg) = 25
            xpos(msg) = xabs
            ypos(msg) = yabs
            yabs = yabs - dyabs
            ht(msg) = 0.14_dp
         endif
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = ' '
         imes(msg) = 1
         xpos(msg) = xabs
         ypos(msg) = yabs
         yabs = yabs - dyabs
         ht(msg) = 0.14_dp
      endif
!
      endif klabel_1
      if (ipass.gt.1) then
         return
      endif
!
      if (itek.ge.5.and.idotek.eq.0) then
        ncurve = nn
        ipag = 0
        ibrdr = 1
        grce = -1.0
        iexit = 1
        nplen = 100
        hight = 0.14_dp
!-----------------------------------------------------------------------
!       Write Plot Parameters
!-----------------------------------------------------------------------
        call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy,  &
        iorel, xorl, yorl, hight, bngle, bshft, &
        ptitle, nplen, xtitle, nxlen, ytitle, nylen, x11x, y11y, &
        almin, xstp, almax, blmin, ystp, blmax, &
        iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
        isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
        igridx, igridy, idash, idot, ichdsh, ichdot, &
        thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
        markme, clearx, mrc, tlen, nmrk, rat, &
        xx, yy, nxy, ncnct, &
        icont, nword, zz, ix, iy, zinc, line, mode, &
        lbflg, ithk, ipri, nline, draw, &
        nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
        nvec, xfm, yfm, xto, yto, ivec, &
        msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif
!
      else is_vacuum_1
      call sets2d(psi,c,rgrid,nw,bkx,lkx,zgrid,nh,bky,lky,wk,ier)
      endif is_vacuum_1
      do i=1,nh
         call seva2d(bkx,lkx,bky,lky,c,rmaxis,zgrid(i),pds,ier,n333)
         workj(i)=-pds(3)/rmaxis
      enddo
      dxnow=(rmaxzm-rminzm)/(nw-1)
      drnow=(xlmax-xlmin)/(nw-1)
      do i=1,nw
         xnow=rminzm+(i-1)*dxnow
         call seva2d(bkx,lkx,bky,lky,c,xnow,zmaxis,pds,ier,n111)
         workc(i)=xnow
         workd(i)=pds(1)
         xnow=xlmin+(i-1)*drnow
         call seva2d(bkx,lkx,bky,lky,c,xnow,zmaxis,pds,ier,n333)
         workf(i)=pds(1)+psiref(jtime)
         workg(i)=pds(2)/xnow
         workk(i)=xnow
      enddo
      is_vacuum_2: if (ivacum.eq.0) then
      if (kwripre.gt.0) then
         delz=(zuperts(jtime)-zlowerts)/(mpress-1)
         do i=1,mpress
            zqthom(i)=(zlowerts+(i-1)*delz)/100.
         enddo
      endif
      do i=1,nw
         do j=1,nh
            kk=(i-1)*nh+j
            copyn(kk)=pcurrt(kk)/darea/1.0e+03_dp
         enddo
      enddo
      call sets2d(copyn,cj,rgrid,nw,bjx,ljx,zgrid,nh,bjy,ljy,wkj,ier)
      if (kvtor.gt.0) then
        n1set=1
        ypsi=0.5_dp
        pres0=prcur4(n1set,ypsi)
        prew0=pwcur4(n1set,ypsi)
        n1set=0
        if (icurrt.eq.1) then
          do i=1,nw
            do j=1,nh
              kk=(i-1)*nh+j
              copyn(kk)=pcurrw(kk)/darea/1.0e+03_dp
            enddo
          enddo
          call sets2d(copyn,cv,rgrid,nw,bvx,lvx,zgrid,nh,bvy, &
                      lvy,wkv,ier)
        endif
      endif
      call zpline(nw,worke,pres,bpmid,cpmid,dpmid)
      if (kwripre.gt.0) then
         do i=1,mpress
            call seva2d(bkx,lkx,bky,lky,c,rmajts,zqthom(i), &
                        pds,ier,n111)
            sinow=(pds(1)-simag)/(psibry-simag)
            qthom(i)=seval(nw,sinow,worke,qpsi,bworkb,cworkb,dworkb)
         enddo
      endif
      delscr=(rmaxzm-rminzm)/2./(nw-1)*0.6_dp
      do i=1,nw
         rpmid(i)=workc(i)
         rscrap(i)=rmaxzm+(i-1)*delscr
         call seva2d(bkx,lkx,bky,lky,c,workc(i),zmaxis,pds,ier,n111)
         xpsikk=(pds(1)-simag)/(psibry-simag)
         sipmid(i)=xpsikk
         rhopmid(i)=seval(nw,xpsikk,sigrid,rhovn,brhovn,crhovn,drhovn)
         if (keecur.gt.0) then
            if (xpsikk.lt.0.0.or.xpsikk.gt.1.0) then
               ermid(i)=0.0
            else
               ermid(i)=eradial(xpsikk,keecur,workc(i),zmaxis)/1000.
            endif
         endif
         if (xpsikk.lt.0.0.or.xpsikk.gt.1.0) then
            pmid(i)=0.0
         else
            pmid(i)=seval(nw,xpsikk,worke,pres,bpmid,cpmid,dpmid)
            xmid(i)=xpsikk
         endif
         pmidw(i)=0.0
         rot: if (kvtor.gt.0.and.icurrt.ne.2) then
           call seva2d(bwx,lwx,bwy,lwy,cw,workc(i),zmaxis, &
                       pds,ier,n111)
           pmid(i)=pds(1)
           ypsi=xpsikk
           pres0=prcur4(n1set,ypsi)
           prew0=pwcur4(n1set,ypsi)
           pmidw(i)=prew0
           if (kvtor.eq.11) then
             rdiml=workc(i)/rzero
             if (abs(pres0).gt.1.e-10_dp) then
               pwop0=prew0/pres0
               ptop0=exp(pwop0*rgrvt(i))
             else
               ptop0=1.0
               pwop0=0.0
             endif
             rgrvn=(workc(i)/rvtor)**2-1.
             if (xpsikk.ge.0.0.and.xpsikk.lt.1.0) then
               ppw= (1.-xpsikk**enw)**emw*(1.-gammaw)+gammaw
             else
               ppw=0.0
             endif
             ppw=rbetaw*ppw*rgrvn*rdiml*ptop0*cratio/darea/1.e3_dp
             workaw(i)=ppw
           else
             if (icurrt.eq.4) then
               rdiml=workc(i)/rzero
               rgrvn=(workc(i)/rvtor)**2-1.
               if (xpsikk.ge.0.0.and.xpsikk.lt.1.0) then
                 ppw= (1.-xpsikk**enw)**emw*(1.-gammaw)+gammaw
               else
                 ppw=0.0
               endif
               ppw=rbetaw*ppw*rgrvn*rdiml*cratio/darea/1.e3_dp
               workaw(i)=ppw
             elseif (icurrt.eq.5) then
               rgrvn=(workc(i)/rvtor)**2-1.
               ppw=pwpcur(xpsikk,kwwcur)
               if (kvtor.eq.1) then
                 ppw=ppw*rgrvn*workc(i)
                 workaw(i)=ppw/darea/1.e3_dp
               elseif (kvtor.eq.2) then
                 if (abs(pres0).gt.1.e-10_dp) then
                    pwop0=prew0/pres0
                    ptop0=exp(pwop0) !*rgrvv) ! TODO: rgrvv undefined
                 else
                    ptop0=1.0
                    pwop0=0.0
                 endif
                 ppw=ppw*(1.+pwop0*rgrvn)*rgrvn*workc(i)
                 workaw(i)=ppw/darea/1.e3_dp
               elseif (kvtor.eq.3) then
                 if (abs(pres0).gt.1.e-10_dp) then
                    pwop0=prew0/pres0
                    ptop0=exp(pwop0*rgrvn)
                 else
                    ptop0=1.0
                    pwop0=0.0
                 endif
                 ppw=ppw*ptop0*rgrvn*workc(i)
                 workaw(i)=ppw/darea/1.e3_dp
               endif
             elseif (icurrt.eq.1) then
               call seva2d(bvx,lvx,bvy,lvy,cv,workc(i),zmaxis, &
                           pds,ier,n111)
               workaw(i)=pds(1)
             endif
           endif
         endif rot
!-------------------------------------------------------------------
!--      total current density                                    --
!-------------------------------------------------------------------
         call seva2d(bkx,lkx,bky,lky,c,rscrap(i),zmaxis,pds,ier,n111)
         xpsiss=(pds(1)-simag)/(psibry-simag)
         if (icurrt.eq.2) then
            if (xpsikk.lt.0.0) xpsikk=0.0
            if (xpsikk.gt.1.0001_dp.and.icutfp.eq.0) then
               worka(i)=0.0
            else
               worka(i)=workc(i)*ppcurr(xpsikk,kppcur) &
                      +fpcurr(xpsikk,kffcur)/workc(i)
               worka(i)=worka(i)/darea/1.0e+03_dp
            endif
            if (xpsiss.lt.0.0) xpsiss=0.0
            if (xpsiss.gt.1.0001_dp.and.icutfp.eq.0) then
               curscr(i)=0.0
               cycle
            endif
            curscr(i)=rscrap(i)*ppcurr(xpsiss,kppcur) &
                         +fpcurr(xpsiss,kffcur)/rscrap(i)
            curscr(i)=curscr(i)/darea/1.0e+03_dp
         else
           call seva2d(bjx,ljx,bjy,ljy,cj,workc(i),zmaxis, &
                       pds,ier,n111)
           worka(i)=pds(1)
         endif
      enddo
!
      curmin=worka(1)
      curmax=worka(1)
      do i=1,nw
         curmin=min(curmin,worka(i))
         curmax=max(curmax,worka(i))
         if (kvtor.gt.0) then
           curmin=min(curmin,workaw(i))
           curmax=max(curmax,workaw(i))
         endif
         curmid(i)=worka(i)
      enddo
      drrrr=(workc(nw)-workc(1))/2.
      drpmid=drrrr
      dcurn=curmax-curmin
      curmax=curmax+0.05_dp*dcurn
      dcurn=(curmax-curmin)/2.
      klabel_2: if (klabel.ge.0) then
!
      pl_files_1: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!        Initialize plot parameters
!-----------------------------------------------------------------------
         call init2d
         ibrdr = 1
         xphy = 4.1_dp
         yphy = 5.5_dp
         hight = 0.10_dp
         nplen = 100
         nxlen = -100
         nylen = 100
         iytck = 1
         ixnon = 1
         igridx = 1
         igridy = 1
         idot = 1
         xtitle = '$'
         ytitle = 'J(Zmag, KA/m2)$'
         iexit = 1

         nn=1
         nxy(nn) = nw
         do ii = 1,nw
           xx(ii,nn) = workc(ii)
           yy(ii,nn) = worka(ii)
         enddo
         plot_mid: if (kwripre.gt.0) then
           xdum=0.
           dataname=dataname(1:lprx)//'_jmid'
           open(unit=62,file=dataname,status='old',iostat=ioerr)
           if (ioerr.eq.0) close(unit=62,status='delete')
           open(unit=62,file=dataname,status='new')
           do i=1,nw
             write (62,*) workc(i),worka(i),xdum,xdum
           enddo
           close(unit=62)
           if (keecur.gt.0) then
             dataname=dataname(1:lprx)//'_ermid'
             open(unit=62,file=dataname,status='old',iostat=ioerr)
             if (ioerr.eq.0) close(unit=62,status='delete')
             open(unit=62,file=dataname,status='new')
             do i=1,nw
               write (62,*) rpmid(i),ermid(i),xdum,xdum
             enddo
             close(unit=62)
           endif
           if (kvtor.gt.0) then
             dataname=dataname(1:lprx)//'_jwmid'
             open(unit=62,file=dataname,status='old',iostat=ioerr)
             if (ioerr.eq.0) close(unit=62,status='delete')
             open(unit=62,file=dataname,status='new')
             do i=1,nw
               write (62,*) workc(i),workaw(i),xdum,xdum
             enddo
             close(unit=62)
           endif
         endif plot_mid
!--------------------------------------------------------------------
!--      zero line                                                 --
!--------------------------------------------------------------------
         nn=nn+1
         nxy(nn) = 2
         ndotme(nn) = 1
         xx(1,nn) = 1.0001_dp*workc(1)
         yy(1,nn) = 0.0
         xx(2,nn) = 0.9999_dp*workc(nw)
         yy(2,nn) = 0.0
         if (kvtor.gt.0) then
           nn=nn+1
           nxy(nn) = nw
           ndshme(nn) = 1
           do ii = 1,nw
             xx(ii,nn) = workc(ii)
             yy(ii,nn) = workaw(ii)
           enddo
         endif
         ncurve=nn
!-----------------------------------------------------------------------
!        Write Plot Parameters
!-----------------------------------------------------------------------
         call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy,  &
         iorel, xorl, yorl, hight, bngle, bshft, &
         ptitle, nplen, xtitle, nxlen, ytitle, nylen, x11x, y11y, &
         almin, xstp, almax, blmin, ystp, blmax, &
         iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
         isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
         igridx, igridy, idash, idot, ichdsh, ichdot, &
         thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
         markme, clearx, mrc, tlen, nmrk, rat, &
         xx, yy, nxy, ncnct, &
         icont, nword, zz, ix, iy, zinc, line, mode, &
         lbflg, ithk, ipri, nline, draw, &
         nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
         nvec, xfm, yfm, xto, yto, ivec, &
         msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_1
      endif klabel_2
!
      do i=1,nw
         sinow=(workd(i)-simag)/(psibry-simag)
         if (sinow.le.0.0) sinow=0.
         if (sinow.gt.1.) sinow=1.
         worka(i)=seval(nw,sinow,worke,qpsi,bworkb,cworkb,dworkb)
      enddo
      curmin=worka(1)
      curmax=worka(1)
      do i=1,nw
         curmin=min(curmin,worka(i))
         curmax=max(curmax,worka(i))
      enddo
      drrrr=(workc(nw)-workc(1))/2.
      dcurn=curmax-curmin
      curmax=curmax+0.05_dp*dcurn
      dcurn=(curmax-curmin)/2.
      klabel_3: if (klabel.ge.0) then
         ibrdr = 1
!-----------------------------------------------------------------------
!        Initialize plot parameters
!-----------------------------------------------------------------------
         call init2d
         nn = nn + 1
         nxy(1) = nw
         do ii = 1, nw
           xx(ii,1) = workc(ii)
           yy(ii,1) = worka(ii)
         enddo
         nn=nn+1
         nxy(nn)=2
         ndotme(nn)=1
         xx(1,nn)=rmaxis
         yy(1,nn)=curmin
         xx(2,nn)=rmaxis
         yy(2,nn)=curmax
         if (nqwant.gt.0) then
            do i=1,nqwant
               sinow=siwantq(i)
               worka(i)=seval(nw,sinow,worke,qpsi,bworkb,cworkb,dworkb)
            enddo
            nn = nn + 1
            nxy(nn) = nqwant
            markme(nn) = 3
            ncnct(nn) = 1
            do ii = 1, nxy(nn)
               xx(ii,nn) = qsiw(ii)
               yy(ii,nn) = worka(ii)
            enddo
         endif
         if (kwripre.gt.0) then
            dataname=dataname(1:lprx)//'_qmid'
            open(unit=62,file=dataname,status='old',iostat=ioerr)
            if (ioerr.eq.0) close(unit=62,status='delete')
            open(unit=62,file=dataname,status='new')
            do ip=1,nw
               write(62,*)workc(ip),worka(ip)
            enddo
            close(62)
         endif
         MSE: if (kstark.gt.0.or.kdomse.gt.0) then
          do i=1,nstark
            sinow=(sistark(i)-simag)/(psibry-simag)
            if (sinow.gt.0.0.and.sinow.lt.1.0) then
              qstark(i)=seval(nw,sinow,worke,qpsi,bworkb,cworkb,dworkb)
            else
              qstark(i)=0.0
            endif
          enddo
          plot_mse: if (kwripre.gt.0) then
            dataname=dataname(1:lprx)//'_rzmse'
            open(unit=62,file=dataname,status='old',iostat=ioerr)
            if(ioerr.eq.0) close(unit=62,status='delete')
            open(unit=62,file=dataname,status='new')
            do ip=1,nstark
              write(62,*)rrgam(jtime,ip),zzgam(jtime,ip)
            enddo
            close(62)
            dataname=dataname(1:lprx)//'_bzmse'
            open(unit=62,file=dataname,status='old',iostat=ioerr)
            if (ioerr.eq.0) close(unit=62,status='delete')
            open(unit=62,file=dataname,status='new')
            zerono = 0.0
            do ip=1,nmselp
              if (kstark.gt.0) then
                sigmabz = 0.0
                if(fwtgam(ip).gt.0.0) sigmabz=bzmse(ip)* &
                       siggam(jtime,ip)/tangam(jtime,ip)
                write(62,*) rrgam(jtime,ip),bzmse(ip), &
                            zerono,abs(sigmabz)
              else
                write(62,*) rrgam(jtime,ip),bzmsec(ip), &
                            zerono,zerono
              endif
              rrgams(ip)=rrgam(jtime,ip)
            enddo
            close(62)
            dataname=dataname(1:lprx)//'_bzlim'
            open(unit=62,file=dataname,status='old',iostat=ioerr)
            if(ioerr.eq.0) close(unit=62,status='delete')
            open(unit=62,file=dataname,status='new')
            do ip=nmselp+1,nstark
              if (nstark.gt.nmselp) then
                if (kstark.gt.0) then
                  sigmabz = 0.0
                  if(fwtgam(ip).gt.0.0) sigmabz=bzmse(ip)* &
                         siggam(jtime,ip)/tangam(jtime,ip)
                  write(62,*) rrgam(jtime,ip),bzmse(ip), &
                              zerono,abs(sigmabz)
                else
                  write(62,*) rrgam(jtime,ip),bzmsec(ip), &
                              zerono,zerono
                endif
                rrgaml(ip-nmselp)=rrgam(jtime,ip)
              endif
            enddo
            close(62)
            dataname=dataname(1:lprx)//'_bzmsec'
            open(unit=62,file=dataname,status='old',iostat=ioerr)
            if (ioerr.eq.0) close(unit=62,status='delete')
            open(unit=62,file=dataname,status='new')
            do jp=1,nmselp
              rrgamn=rrgams(1)
              bzgamn=bzmsec(1)
              igamn=1
              do ip=1,nmselp
                if (rrgams(ip).lt.rrgamn) then
                  rrgamn=rrgams(ip)
                  bzgamn=bzmsec(ip)
                  igamn=ip
                endif
              enddo
              write(62,*)rrgamn,bzgamn
              rrgams(igamn)=1.e10_dp
            enddo
            close(62)
            dataname=dataname(1:lprx)//'_bzlimc'
            open(unit=62,file=dataname,status='old',iostat=ioerr)
            if (ioerr.eq.0) close(unit=62,status='delete')
            open(unit=62,file=dataname,status='new')
            do jp=nmselp+1,nstark
              rrgamn=rrgaml(1)
              bzgamn=bzmsec(1+nmselp)
              igamn=1
              do ip=1,libim
                if (rrgaml(ip).lt.rrgamn) then
                  rrgamn=rrgaml(ip)
                  bzgamn=bzmsec(ip+nmselp)
                  igamn=ip
                endif
              enddo
              write(62,*)rrgamn,bzgamn
              rrgaml(igamn)=1.e10_dp
            enddo
            close(62)
          endif plot_mse
         endif MSE
         if (kstark.eq.0.and.mmbmsels.gt.0) then
          do i=1,nmsels
            sinow=sinmls(i)
            if (sinow.gt.0.0.and.sinow.lt.1.0) then
              qstark(i)=seval(nw,sinow,worke,qpsi,bworkb,cworkb,dworkb)
            else
              qstark(i)=0.0
            endif
          enddo
         endif
!
         pl_files_2: if (itek.ge.5.and.idotek.eq.0) then
           ibrdr = 1
           xphy = 4.1_dp
           yphy = 2.7_dp
           hight = 0.10_dp
           nplen = 100
           nxlen = 100
           nylen = 100
           iytck = 1
           if (keecur.gt.0) then
             igridx = 0
           else
             igridx = 1
           endif
           igridy = 1
           idot = 1
           xtitle = 'R(m)$'
           ytitle = 'q$'
           iexit = 2
!-----------------------------------------------------------------------
!--        overlay ER plot                                            --
!-----------------------------------------------------------------------
           if (keecur.gt.0) then
             eurmin=1.0e+10_dp
             eurmax=-1.0e+10_dp
             do i=1,nw
               eurmin=min(eurmin,ermid(i))
               eurmax=max(eurmax,ermid(i))
             enddo
             ecurn=eurmax-eurmin
             eurmax=eurmax+0.05_dp*ecurn
             eurmin=eurmin-0.05_dp*ecurn
             eurmax=max(abs(eurmin),abs(eurmax))
             eurmin=-eurmax
             ecurn=(eurmax-eurmin)
             sorg = eurmin
             stp = ecurn
             smax = eurmax
             slen = xmm
             sname = 'ER(kV/m)$'
             nslen = -100
             xps = xmm
             yps = 0.0
             nn = nn + 1
             isaxs = nn
             nxy(nn) = nw
             ndshme(nn) = 1
             do ii = 1, nxy(nn)
               xx(ii,nn) = rpmid(ii)
               yy(ii,nn) = ermid(ii)
             enddo
             ncurve = nn
           endif
!-----------------------------------------------------------------------
!          Write Plot Parameters
!-----------------------------------------------------------------------
           call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy,  &
           iorel, xorl, yorl, hight, bngle, bshft, &
           ptitle, nplen, xtitle, nxlen, ytitle, nylen, x11x, y11y, &
           almin, xstp, almax, blmin, ystp, blmax, &
           iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
           isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
           igridx, igridy, idash, idot, ichdsh, ichdot, &
           thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
           markme, clearx, mrc, tlen, nmrk, rat, &
           xx, yy, nxy, ncnct, &
           icont, nword, zz, ix, iy, zinc, line, mode, &
           lbflg, ithk, ipri, nline, draw, &
           nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
           nvec, xfm, yfm, xto, yto, ivec, &
           msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
         endif pl_files_2
      endif klabel_3
!
      endif is_vacuum_2
      istrpl=1
      if (itek.eq.1) call tekall(4010,960,0,0,0)
      call grace(0.0)
      xmm=1.6_dp
!
      iexpand=1
      if (oring(jtime).lt.100.*scrape) iexpand=0
      iexp: if (iexpand.eq.0) then
         call expand(npltbdry,npltxpt,nexexx,xmm,jtime)
      else iexp
         curmin=workg(1)
         curmax=workg(1)
         do i=1,nw
            curmin=min(curmin,workg(i))
            curmax=max(curmax,workg(i))
         enddo
         dcurn=curmax-curmin
         curmax=curmax+0.05_dp*dcurn
         curmin=curmin-0.05_dp*dcurn
         drrrr=(workk(nw)-workk(1))
         dcurn=(curmax-curmin)
         xmm=1.6_dp
!
         pl_files_3: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!        Initialize plot parameters
!-----------------------------------------------------------------------
         call init2d
         ibrdr = 1
         xphy = 9.0
         yphy = 1.4_dp
         hight = 0.12_dp
         nplen = 100
         nxlen = 100
         nylen = 100
         ixtck = 2
         iytck = 2
         igridx = 2
         igridy = 2
         idot = 1
         ipag = 0
         xtitle = 'R(m)$'
         ytitle = 'Bz(T)$'
         iexit = 1
!-----------------------------------------------------------------------
!        Write Plot Parameters
!-----------------------------------------------------------------------
         call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy,  &
         iorel, xorl, yorl, hight, bngle, bshft, &
         ptitle, nplen, xtitle, nxlen, ytitle, nylen, x11x, y11y, &
         almin, xstp, almax, blmin, ystp, blmax, &
         iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
         isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
         igridx, igridy, idash, idot, ichdsh, ichdot, &
         thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
         markme, clearx, mrc, tlen, nmrk, rat, &
         xx, yy, nxy, ncnct, &
         icont, nword, zz, ix, iy, zinc, line, mode, &
         lbflg, ithk, ipri, nline, draw, &
         nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
         nvec, xfm, yfm, xto, yto, ivec, &
         msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_3
      endif iexp
!
      is_vacuum_3: if (ivacum.eq.0) then
      curmin=workf(1)
      curmax=workf(1)
      do i=1,nw
         curmin=min(curmin,workf(i))
         curmax=max(curmax,workf(i))
      enddo
      drrrr=(workk(nw)-workk(1))
      dcurn=curmax-curmin
      curmax=curmax+0.05_dp*dcurn
      curmin=curmin-0.05_dp*dcurn
      dcurn=(curmax-curmin)
      xmm=1.6_dp
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      ibrdr = 1
      xphy = 9.0
      yphy = 3.7_dp
      hight = 0.12_dp
      nplen = 100
      nxlen = 100
      nylen = 100
      ixtck = 2
      iytck = 2
      igridx = 2
      igridy = 2
      idot = 1
      ipag = 0
      if (iexpand .eq. 0) then
         xtitle = 'R(m)$'
         ytitle = 'ABS PSI(V-S/rad)$'
      else
         xtitle = '$'
         ytitle = 'ABS PSI(V-S/rad)$'
         ixnon = 1
      endif
      iexit = 1
!
      if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!        Write Plot Parameters
!-----------------------------------------------------------------------
         call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy,  &
         iorel, xorl, yorl, hight, bngle, bshft, &
         ptitle, nplen, xtitle, nxlen, ytitle, nylen, x11x, y11y, &
         almin, xstp, almax, blmin, ystp, blmax, &
         iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
         isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
         igridx, igridy, idash, idot, ichdsh, ichdot, &
         thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
         markme, clearx, mrc, tlen, nmrk, rat, &
         xx, yy, nxy, ncnct, &
         icont, nword, zz, ix, iy, zinc, line, mode, &
         lbflg, ithk, ipri, nline, draw, &
         nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
         nvec, xfm, yfm, xto, yto, ivec, &
         msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif
!
      do i=1,nw
         do j=1,nh
            kk=(i-1)*nh+j
            copyj(i,j)=0.0
            do m=1,nfsum
               copyj(i,j)=copyj(i,j)+gridfc(kk,m)*brsp(m)
            enddo
            if (ivesel.gt.0) then
               do m=1,nvesel
                  copyj(i,j)=copyj(i,j)+gridvs(kk,m)*vcurrt(m)
               enddo
            endif
            if (iecurr.gt.0) then
               do m=1,nesum
                  copyj(i,j)=copyj(i,j)+gridec(kk,m)*ecurrt(m)
               enddo
            endif
         enddo
      enddo
      do i=1,nw
        do j=1,nh
          k=(i-1)*nh+j
          copyn(k)=copyj(i,j)
        enddo
      enddo
      call sets2d(copyn,cj,rgrid,nw,bjx,ljx,zgrid,nh,bjy,ljy,wkj,ier)
      do i=1,nw
         call seva2d(bjx,ljx,bjy,ljy,cj,workk(i),zmaxis,pds,ier,n555)
         workh(i)=1.-pds(5)/pds(2)*workk(i)
      enddo
      rcur=rcurrt(jtime)/100! m
      zcur=zcurrt(jtime)/100! m
      call seva2d(bjx,ljx,bjy,ljy,cj,rcur,zcur,pds,ier,n555)
      vertn(jtime)=1.-pds(5)/pds(2)*rcur
      curmin=workh(1)
      curmax=workh(1)
      do i=1,nw
         curmin=min(curmin,workh(i))
         curmax=max(curmax,workh(i))
      enddo
      dcurn=(curmax-curmin)/2.
      curmax=curmax+dcurn*0.05_dp
      curmin=curmin-dcurn*0.05_dp
      dcurn=(curmax-curmin)
      ibrdr = 1
      workg(1)=workk(1)
      workg(2)=workk(nw)
      workh(1)=0.
      workh(2)=0.
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      ncurve = 2
      nxy(1) = nw
      do i = 1,nw
         xx(i,1) = workk(i)
         yy(i,1) = workh(i)
      enddo
      nxy(2) = 2
      do i = 1, 2
         xx(i,2) = workg(i)
         yy(i,2) = workh(i)
      enddo
      ndshme(2) = 1
!
      pl_files_4: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      ibrdr = 1
      xphy = 9.0
      yphy = 6.0
      hight = 0.12_dp
      nplen = 100
      nxlen = 100
      nylen = 100
      ixtck = 2
      iytck = 2
      ixnon = 1
      igridx = 2
      igridy = 2
      idot = 1
      ipag = 0
      xtitle = '$'
      ytitle = 'n$'
      iexit = 1
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy, &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, xmm, xmm, &
      workk(1), drrrr, workk(nw), curmin, dcurn, curmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      igridx, igridy, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx, yy, nxy, ncnct, &
      icont, nword, zz, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_4
!
      curmin=workj(1)
      curmax=workj(1)
      do i=1,nh
         curmin=min(curmin,workj(i))
         curmax=max(curmax,workj(i))
      enddo
      drrrr=(zgrid(nh)-zgrid(1))
      dcurn=curmax-curmin
      curmax=curmax+0.05_dp*dcurn
      curmin=curmin-0.05_dp*dcurn
      dcurn=(curmax-curmin)
!
      pl_files_5: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      ipag = 0
      ibrdr = 1
      xphy = 6.6_dp
      yphy = 1.4_dp
      hight = 0.12_dp
      nplen = -100
      nxlen = 100
      nylen = 100
      ixtck = 2
      iytck = 2
      igridx = 2
      igridy = 2
      idot = 1
      xtitle = 'Z(m)$'
      ytitle = 'Br(T)$'
      iexit = 1
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy,  &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, x11x, y11y, &
      almin, xstp, almax, blmin, ystp, blmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      igridx, igridy, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx, yy, nxy, ncnct, &
      icont, nword, zz, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_5
!
      curmin=1.0e+10_dp
      curmax=-1.0e+10_dp
      do i=1,nw
         sinow=worke(i)
         worka(i)=speval(nw,sinow,worke,qpsi,bworkb,cworkb,dworkb)
         worka(i)=worka(i)/qpsi(i)
         curmin=min(curmin,worka(i))
         curmax=max(curmax,worka(i))
      enddo
      dcurn=curmax-curmin
      curmax=curmax+0.05_dp*dcurn
      curmin=curmin-0.05_dp*dcurn
      dcurn=(curmax-curmin)
!
      pl_files_6: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      ibrdr = 1
      xphy = 6.5_dp
      yphy = 3.7_dp
      hight = 0.12_dp
      nplen = -100
      nxlen = 100
      nylen = 100
      ixtck = 1
      iytck = 2
      igridx = 1
      igridy = 2
      idot = 1
      xnmax = 0.5_dp
      ipag = 0
      xtitle = 'PSI$'
      ytitle = 'SHEAR$'
      iexit = 1
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy,  &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, x11x, y11y, &
      almin, xstp, almax, blmin, ystp, blmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      igridx, igridy, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx, yy, nxy, ncnct, &
      icont, nword, zz, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_6
      curmin=1.0e+10_dp
      curmax=-1.0e+10_dp
      do i=1,nw
         curmin=min(curmin,pmid(i))
         curmax=max(curmax,pmid(i))
      enddo
      dcurn=curmax-curmin
      curmax=curmax+0.05_dp*dcurn
      curmin=curmin
      dcurn=(curmax-curmin)
      drpmid=drpmid*2.
!
      pl_files_7: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      xtitle = 'R(m)$'
      ytitle = 'P(n/m2)$'
      ibrdr = 1
      xphy = 4.0
      yphy = 1.4_dp
      hight = 0.12_dp
      nplen = -100
      nxlen = 100
      nylen = 100
      ixtck = 2
      iytck = 2
      igridx = 2
      igridy = 2
      idot = 1
      iexit = 1
!
      if (kwripre.gt.0) then
         dataname=dataname(1:lprx)//'_pmid'
         open(unit=62,file=dataname,status='old',err=17940)
         close(unit=62,status='delete')
17940    continue
         open(unit=62,file=dataname,status='new')
         do i=1,nw
            write (62,*) rpmid(i),pmid(i)
         enddo
         close(unit=62)
         dataname=dataname(1:lprx)//'_pwmid'
         open(unit=62,file=dataname,status='old',err=17950)
         close(unit=62,status='delete')
17950    continue
         open(unit=62,file=dataname,status='new')
         do i=1,nw
            write (62,*) rpmid(i),pmidw(i)
         enddo
         close(unit=62)
      endif
!
      nn=1
      nxy(nn) = nw
      do ii = 1,nw
        xx(ii,nn) = rpmid(ii)
        yy(ii,nn) = pmid(ii)
      enddo
      nn=nn+1
      nxy(nn) = nw
      ndshme(nn) = 1
      do ii = 1,nw
        xx(ii,nn) = rpmid(ii)
        yy(ii,nn) = pmidw(ii)
      enddo
      ncurve=nn
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy, &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, xmm, xmm, &
      rpmid(1), drpmid, rpmid(nw), curmin, dcurn, curmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      igridx, igridy, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx   , yy  , nxy, ncnct, &
      icont, nword, zz, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_7
      curmin=1.0e+10_dp
      curmax=-1.0e+10_dp
      do i=1,nitera
         curmin=min(curmin,cerror(i))
         curmax=max(curmax,cerror(i))
         xiter(i)=i
      enddo
      xdel=(xiter(nitera)-xiter(1))/2.
      ibrdr = 1
      if (abs(curmin).gt.1.e-10_dp) then
         ycycle=-log10(curmin)
      else
         ycycle=3.
      endif
      kcycle=ycycle+1
      curmin=10.**(-kcycle)
      ycycle=xmm/kcycle
      if (ycycle.lt.xmm/3.) ycycle=xmm/3.
      idel=nitera/xmm
      xdel=idel+1
!
      pl_files_8: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      ibrdr = 1
      xphy = 4.0
      yphy = 3.7_dp
      hight = 0.12_dp
      nplen = -100
      nxlen = 100
      nylen = 100
      ixtck = 2
      idot = 1
      xstp = 0.0
      ystp = 0.0
      iaxis = 2
      intax = 1
      xtitle = 'ITERS$'
      ytitle = 'ERR$'
      iexit = 1
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy,  &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, x11x, y11y, &
      almin, xstp, almax, blmin, ystp, blmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      igridx, igridy, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx, yy, nxy, ncnct, &
      icont, nword, zz, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_8
      not_eq_mode: if (iconvr.ne.3) then
         curmin=1.0e+10_dp
         curmax=-1.0e+10_dp
         do i=1,nitera
            curmin=min(curmin,cchisq(i))
            curmax=max(curmax,cchisq(i))
         enddo
         ibrdr = 1
         ycycle= log10(curmax)
         ycycle=max(one,ycycle)
         kcycle=ycycle
         ycycle=xmm/kcycle
         if (ycycle.lt.xmm/3.) ycycle=xmm/3.
         curmin=   (log10(curmin))
         kcycle=curmin
         curmin=10.**(kcycle)
!
         pl_files_9: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!          Initialize plot parameters
!-----------------------------------------------------------------------
           call init2d
           ibrdr = 1
           xphy = 4.0
           yphy = 6.0
           hight = 0.12_dp
           nplen = 100
           nxlen = 100
           nylen = 100
           iaxis = 2         
           ixtck = 2
           intax = 1
           xstp = 0.0
           ystp = 0.0
           xtitle = 'ITERS$'
           ytitle = 'CHI2$'
           iexit = 1
!-----------------------------------------------------------------------
!         Write Plot Parameters
!-----------------------------------------------------------------------
          call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy,  &
          iorel, xorl, yorl, hight, bngle, bshft, &
          ptitle, nplen, xtitle, nxlen, ytitle, nylen, x11x, y11y, &
          almin, xstp, almax, blmin, ystp, blmax, &
          iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
          isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
          igridx, igridy, idash, idot, ichdsh, ichdot, &
          thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
          markme, clearx, mrc, tlen, nmrk, rat, &
          xx, yy, nxy, ncnct, &
          icont, nword, zz, ix, iy, zinc, line, mode, &
          lbflg, ithk, ipri, nline, draw, &
          nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
          nvec, xfm, yfm, xto, yto, ivec, &
          msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
        endif pl_files_9
      endif not_eq_mode
      endif is_vacuum_3
!-----------------------------------------------------------------
!--   F coil current limits from Hyatt                          --
!-----------------------------------------------------------------
      do i=1,4
         workc(i)=348.
         workc(i+9)=348.
      enddo
      workc(5)=319.
      workc(14)=319.
      workc(6)=687.
      workc(15)=687.
      workc(7)=687.
      workc(16)=687.
      workc(8)=348.
      workc(17)=348.
      workc(9)=495.
      workc(18)=495.
      do i=1,nfsum
         workb(i)=brsp(i)/1000.
         worka(i)=i
         workd(i)=-workc(i)
         if (fwtfc(1).ne.0.0.or.ivacum.eq.1.or.imag2(jtime).gt.0) then
            workc(i)=fccurt(jtime,i)/1000.
            workd(i)=workc(i)
         endif
      enddo
      curmin=workb(1)
      curmax=workb(1)
      do i=1,nfsum
         curmin=min(curmin,workb(i),workc(i),workd(i))
         curmax=max(curmax,workb(i),workc(i),workd(i))
      enddo
      dcurn=(curmax-curmin)/10.
      curmax=curmax+dcurn
      curmin=curmin-dcurn
      dcurn=(curmax-curmin)
      ibrdr = 1
      iline=1
      if (fwtfc(1).ne.0.0) iline=-1
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
!-----------------------------------------------------------------------
!     Computed F-coil currents, solid line, open circles, cyan
!-----------------------------------------------------------------------
      nn = nn + 1
      nxy(nn) = nfsum
      ncnct(nn) = 1
      markme(nn) = 16
      sclpc(nn) = 0.7_dp
      clearx(nn)='CYAN'
      do i = 1, nfsum
         xx(i,nn) = worka(i)
         yy(i,nn) = workb(i)
      enddo
      F_coils: if (ifcurr.le.0) then
      if (imag2(jtime).eq.0) then
         nn = nn + 1
         nxy(nn) = nfsum
         ncnct(nn) = 1
         markme(nn) = 15
         sclpc(nn) = 0.7_dp
         ndotme(nn) = 1
         do i = 1, nfsum
            xx(i,nn) = worka(i)
            yy(i,nn) = workc(i)
         enddo
         nn = nn + 1
         nxy(nn) = nfsum
         ncnct(nn) = 1
         markme(nn) = 15
         sclpc(nn) = 0.7_dp
         ndotme(nn) = 1
         do i = 1, nfsum
            xx(i,nn) = worka(i)
            yy(i,nn) = workd(i)
         enddo
      else
!------------------------------------------------------------------------
!        Measured F-coil currents, dot curve, solid pink symbols
!------------------------------------------------------------------------
         nn = nn + 1
         nxy(nn) = nfsum
         ncnct(nn) = 1
         markme(nn) = 15
         sclpc(nn) = 0.7_dp
         clearx(nn)='PINK'
         ndotme(nn) = 1
         do i = 1, nfsum
            xx(i,nn) = worka(i)
            yy(i,nn) = workc(i)
         enddo
      endif
      endif F_coils
!-----------------------------------------------------------------------
!--   write out plasma parameters                                     --
!-----------------------------------------------------------------------
      xabs=-6.0
      yabs=1.8_dp
      dyabs = 0.22_dp
      write(text,8950) trim(ch1),trim(ch2),efitvers
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 27
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,8960) uday
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9000) ishot
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9020) itime,itimeu
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9440) chi2rm(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9270) pasmac
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9420) gapin(jtime),gapout(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9460) gaptop(jtime),gapbot(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9441) fexpan,fexpvs
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9467) zuperts(jtime),rlibim(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9520) rvsiu(jtime),rvsou(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9521) rvsid(jtime),rvsod(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9530) rseps(1,jtime),rseps(2,jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9540) zseps(1,jtime),zseps(2,jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9470) sibdry(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9480) elongm(jtime),qm(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9482) psiref(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9484) vertn(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
!-----------------------------------------------------------------------
!--   vertical stability parameter, reference Nuc Fusion 18(1978)1331 --
!-----------------------------------------------------------------------
      if (ivacum.eq.0) then
         rx=rm(jtime)/100.
         pleng=0.0
         f_0=log(8*rout(jtime)/abar)-2+betap(jtime)+li(jtime)/2+.5_dp
         delr=rout(jtime)/100.-1.67_dp
!-----------------------------------------------------------------------
!--      metal wall                                                   --
!-----------------------------------------------------------------------
         xnnc(jtime)=vertn(jtime)/((10.77_dp*delr**2+8.08_dp*delr+2.54_dp)/f_0)
      endif
      write (text,9399) xnnc(jtime),qmerci(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9492) shearb(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      edgesw=cjor95(jtime)*cjor0(jtime)
      write (text,9490) ssi95(jtime),edgesw
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9496) psiwant,qsiwant(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      edgesw=cjorsw(jtime)*cjor0(jtime)
      write (text,9498) ssiwant(jtime),edgesw
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      dnev2=dco2v(jtime,2)
      write (text,9502) dnev2
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      dner1=dco2r(jtime,1)
      erreq=max(dbpli(jtime),delbp(jtime))
      write (text,9510) erreq
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
!----------------------------------------------------------------------------
!--   dummy write statement added for HP version 9.01 optimization problems--
!----------------------------------------------------------------------------
!      if (iand(iout,1).ne.0) write (nout,*) jtime
      write (text,9494) aaq1(jtime),aaq2(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9512) taumhd(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9548) taudia(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      if (q95(jtime).gt.0.1_dp) then
         sq295=ssi95(jtime)/q95(jtime)**2
      else
         sq295=0.0
      endif
      pbinjmw=pbinj(jtime)/1.e+06_dp
      pasmma=abs(ipmeas(jtime))/1.e+06_dp
      routm=rout(jtime)/100.
      tauenn=tauthn(jtime)
!     if (pbinjmw.gt.1.e-03_dp.and.abs(pasmma).gt.1.e-03_dp) then
!     taud3jet=106._dp*pbinjmw**(-0.46_dp)*pasmma**1.03_dp*routm**1.48_dp
!     tauenn=taumhd(jtime)/taud3jet
!     else
!     tauenn=0.0
!     endif
      if (wmhd(jtime).gt.1000.) then
         wfn_bim=wfbeam(jtime)/wmhd(jtime)
      else
         wfn_bim=0.0
      endif
      write (text,9465) sq295,tauenn
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9550) vloopt(jtime),cj1ave(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9934) tave(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9936) pbinjmw,wfn_bim
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      if (isetfb.ne.0) then
         brfbtka=brfb(1)*nfbcoil/1000.
         write (text,9938) brfbtka
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = text
         imes(msg) = 25
         xpos(msg) = xabs
         ypos(msg) = yabs
         yabs = yabs - dyabs
         ht(msg) = 0.14_dp
      elseif (symmetrize) then
         write (text,9937)
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = text
         imes(msg) = 25
         xpos(msg) = xabs
         ypos(msg) = yabs
         yabs = yabs - dyabs
         ht(msg) = 0.14_dp
      endif
      if (kstark.gt.0.or.mmbmsels.gt.0) then
         !kstnow=mse315/2   ! TODO: mse315 is undefined
         if (mmbmsels.eq.0) then
           write (text,9930) chigamt,chilibt
         else
           write (text,99930) chigamt,tchimls
         endif
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = text
         imes(msg) = 30
         xpos(msg) = xabs
         ypos(msg) = yabs
         yabs = yabs - dyabs
         ht(msg) = 0.14_dp
         if (kstark.gt.0) then
         !  write (text,9940) (qstark(i),i=1,kstnow)
         elseif (mmbmsels.gt.0) then
           write (text,99940) (qstark(i),i=1,nmsels/2)
         endif
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = text
         imes(msg) = 70
         xpos(msg) = xabs
         ypos(msg) = yabs
         yabs = yabs - dyabs
         ht(msg) = 0.14_dp
         if (kstark.gt.0) then
         !  write (text,9940) (qstark(i),i=kstnow+1,mse315)
         elseif (mmbmsels.gt.0) then
           write (text,99940) (qstark(i),i=nmsels/2+1,nmsels)
         endif
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = text
         imes(msg) = 70
         xpos(msg) = xabs
         ypos(msg) = yabs
         yabs = yabs - dyabs
         ht(msg) = 0.14_dp
      endif
      if (fwtdlc.gt.0.0) then
         write (text,9495) chidflux
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = text
         imes(msg) = 25
         xpos(msg) = xabs
         ypos(msg) = yabs
         yabs = yabs - dyabs
         ht(msg) = 0.14_dp
      endif
      if (nqwant.gt.0) then
         kstnow=min(nqwant,8)
         write (text,9945) (qsiw(i),i=1,kstnow)
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = text
         imes(msg) = 70
         xpos(msg) = xabs
         ypos(msg) = yabs
         yabs = yabs - dyabs
         ht(msg) = 0.14_dp
      endif
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = ' '
      imes(msg) = 1
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
!
      pl_files_10: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      ncurve = nn
      ibrdr = 1
      xphy = 6.5_dp
      yphy = 6.0
      hight = 0.12_dp
      nplen = 100
      nxlen = 100
      nylen = 100
      ixtck = 2
      iytck = 2
      igridx = 2
      igridy = 2
      idot = 1
      intax = 1
      xtitle = 'F coils$'
      ytitle = 'I(kA)$'
      iexit = 2
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      xfcoil=nfsum-1
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy, &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, xmm, xmm, &
      worka(1), xfcoil, worka(nfsum), curmin, dcurn, curmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      igridx, igridy, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx, yy, nxy, ncnct, &
      icont, nword, zz, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_10
!----------------------------------------------------------------------
!--   plot flux loop chi squares                                     --
!----------------------------------------------------------------------
      eq_file_mode: if (iconvr.eq.3.and.((kdata.eq.1).or.(kdata.eq.2))) then
      if (itek.eq.1) call tekall(4010,960,0,0,0)
      xmm=1.6_dp*2
      else
      chsmin=1.0e30_dp
      chsmax=-1.0e30_dp
      do i=1,nsilop
         chsmin=min(chsmin,saisil(i))
         chsmax=max(chsmax,saisil(i))
         si(i)=i
      enddo
!
      pl_files_11: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      if (kstark.le.1) then
         yxxx=3.0
         xphy = 7.7_dp
         yphy = 1.0
      else
         yxxx=1.8_dp
         xphy = 7.5_dp
         yphy = 0.9_dp
      endif
      ibrdr = 1
      hight = 0.12_dp
      nplen = -100
      nxlen = 100
      nylen = 100
      xlen = 3.0
      if (ivacum.eq.1) then
         xlen = 6.0
         xphy = 4.7_dp
      endif
      xstp = si(nsilop) - si(1)
      ystp = chsmax - chsmin
      intax = 1
      igridy = 2
      igridx = 2
      ixtck = 2
      iytck = 2
      idot = 1
      ipag = 1
      xtitle = 'PSI LOOPS$'
      ytitle = 'CHI**2$'
      sclpc(1) = 0.5_dp
      ncnct(1) = 1
      iexit = 1
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy,  &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, x11x, y11y, &
      almin, xstp, almax, blmin, ystp, blmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      igridx, igridy, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx, yy, nxy, ncnct, &
      icont, nword, zz, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_11
!----------------------------------------------------------------------
!--   plot the probes chi squares                                     --
!----------------------------------------------------------------------
      chsmin=1.0e30_dp
      chsmax=-1.0e30_dp
      do i=1,magpri
         chsmin=min(saimpi(i),chsmin)
         chsmax=max(saimpi(i),chsmax)
         rmpi(i)=i
      enddo
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      if (kstark.le.1) then
         xphy = 7.7_dp
         yphy = 5.0
      else
         xphy = 7.5_dp
         yphy = 3.5_dp
      endif
      xtitle = 'MAGNETIC PROBES$'
      ytitle = 'CHI**2$'
      ipag = 1
      iexit = 1
      if (ivacum.eq.1) iexit=2
      sclpc(1) = 0.5_dp
      ncnct(1) = 1

      xabs=-7.3_dp
      if (ivacum.eq.1) xabs=-4.3_dp
      if (kstark.le.1) then
         yabs=3.0
      else
         yabs=4.5_dp
      endif
      dyabs = 0.22_dp
      write(text,8950) trim(ch1),trim(ch2),efitvers
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 27
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,8960) uday
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9000) ishot
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9020) itime,itimeu
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
!-----------------------------------------------------------------------
!--   print out pointnames                                            --
!-----------------------------------------------------------------------
      if (kstark.le.1) then
         yabs=1.5_dp
      else
         yabs=3.0
      endif
      yabs0=yabs
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = '       probe identification$'
      imes(msg) = 100
      xpos(msg) = xabs
      ypos(msg) = yabs
      ht(msg) = 0.10_dp
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = 'no.   name        no.   name$'
      imes(msg) = 100
      xpos(msg) = xabs
      ypos(msg) = yabs - 0.25_dp
      ht(msg) = 0.10_dp
      yabs=yabs0-0.5_dp
      do i=magpri67+1,magpri67+16
!        call intno(i,xabs,yabs)
         msg = msg + 1
         note(msg) = 5
         inum(msg) = i
         xpos(msg) = xabs
         ypos(msg) = yabs
         ht(msg) = 0.07_dp
!        call messag(%ref(mpnam2(i)),10,-6.8_dp,yabs)
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = mpnam2(i)
         imes(msg) = 10
         xpos(msg) = -6.8_dp
         if (ivacum.eq.1) xpos(msg) = -3.8_dp
         ypos(msg) = yabs
         ht(msg) = 0.07_dp
         yabs=yabs-.13_dp
      enddo
      xabs=-5.8_dp
      if (ivacum.eq.1) xabs=-2.8_dp
      yabs=yabs0-0.5_dp
      do i=magpri67+17,magpri67+magpri322
!        call intno(i,xabs,yabs)
         msg = msg + 1
         note(msg) = 5
         inum(msg) = i
         xpos(msg) = xabs
         ypos(msg) = yabs
         ht(msg) = 0.07_dp
!        call messag(%ref(mpnam2(i)),10,-5.2,yabs)
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = mpnam2(i)
         imes(msg) = 10
         xpos(msg) = -5.2_dp
         if (ivacum.eq.1) xpos(msg) = -2.2_dp
         ypos(msg) = yabs
         ht(msg) = 0.07_dp
         yabs=yabs-.13_dp
      enddo
!
      xabs=-7.3_dp
      if (ivacum.eq.1) xabs=-4.3_dp
      yabs=yabs-.5_dp
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = '     psi loop identification$'
      imes(msg) = 100
      xpos(msg) = xabs
      ypos(msg) = yabs
      ht(msg) = 0.10_dp
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = 'no.   name        no.   name$'
      imes(msg) = 100
      xpos(msg) = xabs
      ypos(msg) = yabs - 0.25_dp
      ht(msg) = 0.10_dp
!     call height(.07)
      yabs=yabs-.37_dp
      yyabs=yabs
      do i=1,21
         yabs=yabs-.13_dp
!        call intno(i,xabs,yabs)
         msg = msg + 1
         note(msg) = 5
         inum(msg) = i
         xpos(msg) = xabs
         ypos(msg) = yabs
         ht(msg) = 0.07_dp
!        call messag(%ref(lpname(i)),10,-6.8,yabs)
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = lpname(i)
         imes(msg) = 10
         xpos(msg) = -6.8_dp
         if (ivacum.eq.1) xpos(msg) = -3.8_dp
         ypos(msg) = yabs
         ht(msg) = 0.07_dp
      enddo
      yabs=yyabs
      do i=22,nsilop
         yabs=yabs-.13
!        call intno(i,-5.8,yabs)
         msg = msg + 1
         note(msg) = 5
         inum(msg) = i
         xpos(msg) = -5.8_dp
         ypos(msg) = yabs
         ht(msg) = 0.07_dp
!        call messag(%ref(lpname(i)),10,-5.2,yabs)
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = lpname(i)
         imes(msg) = 10
         xpos(msg) = -5.2_dp
         if (ivacum.eq.1) xpos(msg) = -2.2_dp
         ypos(msg) = yabs
         ht(msg) = 0.07_dp
      enddo
!
      pl_files_12: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      ibrdr = 1
      hight = 0.12_dp
      nplen = 100
      nxlen = 100
      nylen = 100
      xlen = 3.0
      if (ivacum.eq.1) then
         xlen = 6.0
         xphy = 4.7_dp
      endif
      xstp = rmpi(magpri) - rmpi(1)
      ystp = chsmax - chsmin
      intax = 1
      igridx = 2
      igridy = 2
      ixtck = 2
      iytck = 2
      idot = 1
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy,  &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, x11x, y11y, &
      almin, xstp, almax, blmin, ystp, blmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      igridx, igridy, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx, yy, nxy, ncnct, &
      icont, nword, zz, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_12
!--------------------------------------------------------------------
!--   plot chi2 for MSE                                            --
!--------------------------------------------------------------------
      MSE_chi: if (kstark.gt.1) then
       chsmin=1.0e30_dp
       chsmax=-1.0e30_dp
       do i=1,nstark
         chsmin=min(chigam(i),chsmin)
         chsmax=max(chigam(i),chsmax)
         rmpi(i)=i
       enddo
!
       pl_files_13: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!        Initialize plot parameters
!-----------------------------------------------------------------------
         call init2d
         yorg = chsmin
         ynmax = chsmax
         if (ynmax.le.1.) ynmax=1.
         nn = nn + 1
         ndshme(nn) = 1
         nxy(nn) = nstark
         do ii = 1, nxy(nn)
            xx(ii,nn) = ii
            yy(ii,nn) = chigam(ii)
         enddo
!
         chsmin=1.0e30_dp
         chsmax=-1.0e30_dp
         do i=1,nstark
            anglem(i)=atan(tangam(jtime,i)/180.*pi)
            anglec(i)=atan(cmgam(i,jtime)/180.*pi)
            chsmin=min(anglem(i),chsmin)
            chsmax=max(anglem(i),chsmax)
            chsmin=min(anglec(i),chsmin)
            chsmax=max(anglec(i),chsmax)
         enddo
         dangle=chsmax-chsmin
         chsmax=chsmax+0.1_dp*dangle
         chsmin=chsmin-0.1_dp*dangle
         dangle=(chsmax-chsmin)/3.
         isaxs = 2
         sorg = chsmin
         stp = chsmax - chsmin
         smax = chsmax
         slen = yxxx
         sname = 'TH(d)$'
         nslen = -100
         xps = 3.0
         yps = 0.0
         nn = nn + 1
         ncnct(nn) = -1
         nxy(nn) = nstark
         sclpc(nn) = 0.7_dp
         do ii = 1, nxy(nn)
            xx(ii,nn) = ii
            yy(ii,nn) = anglem(ii)
         enddo
         nn = nn + 1
         nxy(nn) = nstark
         do ii = 1, nxy(nn)
            xx(ii,nn) = ii
            yy(ii,nn) = anglec(ii)
         enddo
         ibrdr = 1
         xphy = 7.5_dp
         yphy = 6.1_dp
         hight = 0.12_dp
         nplen = 100
         nxlen = 100
         nylen = 100
         xlen = 3.0
         xstp = xx(nstark,1) - xx(1,1)
         ystp = ynmax - yorg
         intax = 1
         igridx = 2
         igridy = 2
         ixtck = 2
         iytck = 2
         idot = 1
         xtitle = 'POLARIMETRY$'
         ytitle = 'CHI**2$'
         ncurve = nn
         iexit = 1
!-----------------------------------------------------------------------
!        Write Plot Parameters
!-----------------------------------------------------------------------
         call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy,  &
         iorel, xorl, yorl, hight, bngle, bshft, &
         ptitle, nplen, xtitle, nxlen, ytitle, nylen, x11x, y11y, &
         almin, xstp, almax, blmin, ystp, blmax, &
         iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
         isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
         igridx, igridy, idash, idot, ichdsh, ichdot, &
         thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
         markme, clearx, mrc, tlen, nmrk, rat, &
         xx, yy, nxy, ncnct, &
         icont, nword, zz, ix, iy, zinc, line, mode, &
         lbflg, ithk, ipri, nline, draw, &
         nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
         nvec, xfm, yfm, xto, yto, ivec, &
         msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
       endif pl_files_13
      endif MSE_chi
!-------------------------------------------------------------------------
!--   plot P', FF', and Zm                                              --
!-------------------------------------------------------------------------
      is_vacuum_4: if (ivacum.eq.0) then
      curmin=1.0e+10_dp
      curmax=-1.0e+10_dp
      do i=1,nitera
         curmin=min(curmin,czmaxi(i))
         curmax=max(curmax,czmaxi(i))
         xiter(i)=i
      enddo
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      xtitle = 'ITERS$'
      ytitle = 'Zm(cm)$'
      ipag = 1
      iexit = 1
      nn = 1
      nxy(nn) = nitera
      do ii = 1, nxy(nn)
         xx(ii,nn) = xiter(ii)
         yy(ii,nn) = czmaxi(ii)
      enddo
      xdel=(xiter(nitera)-xiter(1))
      dcurn=curmax-curmin
      curmax=curmax+0.06_dp*dcurn
      curmin=curmin-0.06_dp*dcurn
      dcurn=(curmax-curmin)
      xmm=2.0
      yorg = curmin
      ynmax = curmax
      ystp = ynmax - yorg
!---------------------------------------------------------------------
!--   plot feedback currents                                        --
!---------------------------------------------------------------------
      if (.not.fitdelz.and.(abs(vcurfb(1)).gt.1.e-6_dp)) then
         curmin=1.0e+10_dp
         curmax=-1.0e+10_dp
         do i=1,nitera
            czmcm(i)=tvfbrt(i)
            curmin=min(curmin,czmcm(i))
            curmax=max(curmax,czmcm(i))
         enddo
         dcurn=curmax-curmin
         curmax=curmax+0.05_dp*dcurn
         curmin=curmin-0.05_dp*dcurn
         dcurn=(curmax-curmin)
         isaxs = 2
         sorg = curmin
         stp = dcurn
         smax = curmax
         slen = xmm
         sname = 'Ifb(A)$'
         nslen = -100
         xps = xmm
         yps = 0.0
         nn = nn + 1
         nxy(nn) = nitera
         ndshme(nn) = 1
         do ii = 1, nxy(nn)
            xx(ii,nn) = xiter(ii)
            yy(ii,nn) = czmcm(ii)
         enddo
      endif
!-------------------------------------------------------------
!--   add DELZ iterations                                   --
!-------------------------------------------------------------
      if (fitdelz) then
         curmin=1.0e+10_dp
         curmax=-1.0e+10_dp
         do i=1,nitera
            czmcm(i)=cdelz(i)
            curmin=min(curmin,czmcm(i))
            curmax=max(curmax,czmcm(i))
         enddo
         dcurn=curmax-curmin
         curmax=curmax+0.05_dp*dcurn
         curmin=curmin-0.05_dp*dcurn
         dcurn=(curmax-curmin)
         isaxs = 2
         sorg = curmin
         stp = dcurn
         smax = curmax
         slen = xmm
         sname = 'dZ(cm)$'
         nslen = -100
         xps = xmm
         yps = 0.0
         nn = nn + 1
         nxy(nn) = nitera
         ndshme(nn) = 1
         do ii = 1, nxy(nn)
            xx(ii,nn) = xiter(ii)
            yy(ii,nn) = czmcm(ii)
         enddo
      endif
      ncurve = nn
!
      if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!       Initialize plot parameters
!-----------------------------------------------------------------------
        ibrdr = 1
        xphy = 4.2_dp
        yphy = 5.7_dp
        hight = 0.12_dp
        nplen = 100
        nxlen = 100
        nylen = 100
        ixtck = 2
        iytck = 2
        intax = 1
        igridx = 2
        igridy = 2
        idot = 1
!-----------------------------------------------------------------------
!       Write Plot Parameters
!-----------------------------------------------------------------------
        call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy, &
        iorel, xorl, yorl, hight, bngle, bshft, &
        ptitle, nplen, xtitle, nxlen, ytitle, nylen, xmm, xmm, &
        xiter(1), xdel, xiter(nitera), yorg  , ystp   ,ynmax , &
        iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
        isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
        igridx, igridy, idash, idot, ichdsh, ichdot, &
        thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
        markme, clearx, mrc, tlen, nmrk, rat, &
        xx   , yy    , nxy   , ncnct, &
        icont, nword, zz, ix, iy, zinc, line, mode, &
        lbflg, ithk, ipri, nline, draw, &
        nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
        nvec, xfm, yfm, xto, yto, ivec, &
        msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif
!
      curmin=1.0e+10_dp
      curmax=-1.0e+10_dp
      do i=1,nw
         curmin=min(curmin,pprime(i))
         curmax=max(curmax,pprime(i))
         xiter(i)=real(i-1,dp)/(nw-1)
      enddo
      xdel=0.5_dp
      dcurn=curmax-curmin
      if (abs(dcurn).le.1.e-4_dp) dcurn=5.*curmax
      curmax=curmax+0.05_dp*dcurn
      curmin=curmin-0.05_dp*dcurn
      dcurn=(curmax-curmin)
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      yorg = curmin
      ystp = dcurn
      ynmax = curmax
      nn = nn + 1
      nxy(1) = nw
      do ii = 1, nw
         xx(ii,1) = xiter(ii)
         yy(ii,1) = pprime(ii)
      enddo
      curmin=1.0e+10_dp
      curmax=-1.0e+10_dp
      do i=1,nw
         czmcm(i)=ffprim(i)/tmu/twopi
         curmin=min(curmin,czmcm(i))
         curmax=max(curmax,czmcm(i))
      enddo
      dcurn=curmax-curmin
      if (abs(dcurn).le.1.e-4_dp) dcurn=5.*curmax
      curmax=curmax+0.05_dp*dcurn
      curmin=curmin-0.05_dp*dcurn
      dcurn=(curmax-curmin)
      isaxs = 2
      sorg = curmin
      stp = dcurn
      smax = curmax
      slen = xmm
      sname = 'FFPRIM$'
      nslen = -100
      xps = xmm
      yps = 0.0
      nn = nn + 1
      nxy(2) = nw
      ndshme(2) = 1
      do ii = 1, nxy(2)
         xx(ii,2) = xiter(ii)
         yy(ii,2) = czmcm(ii)
      enddo
      xtitle = 'XSI$'
      ytitle = 'PPRIME$'
      ncurve = nn
      ipag = 1
      iexit = 2
      xabs= -0.5_dp
      xabs= -0.7_dp
      yabs= -0.8_dp
      if ((icurrt.eq.2.or.icurrt.eq.5) &
               .and.abs(brsp(nfsum+1)).gt.1.e-10_dp) then
         do i=nfsum+1,nfsum+kppcur
            xxnorm=brsp(i)/brsp(nfsum+1)
            if (i.eq.nfsum+1) then
               write (text,18950) xxnorm
               msg = msg + 1
               note(msg) = 1
               lmes(msg) = text
               imes(msg) = 18
               xpos(msg) = xabs
               ypos(msg) = yabs
               yabs = yabs - dyabs
               ht(msg) = 0.13_dp*0.7_dp
            else
               write (text,18960) xxnorm
               msg = msg + 1
               note(msg) = 1
               lmes(msg) = text
               imes(msg) = 17
               xpos(msg) = xabs
               ypos(msg) = yabs
               yabs = yabs - dyabs
               ht(msg) = 0.13_dp*0.7_dp
            endif
         enddo
         xxnorm=brsp(nfsum+1)/darea
         write (text,18980) xxnorm
      else
         xxnorm=0.0
      endif
!
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 17
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.13_dp*0.7_dp
      write (text,18983) condno
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 17
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.13_dp*0.7_dp
!
      pr=pres(1)/(.667_dp*wmhd(jtime)/(volume(jtime)/1.e6_dp))
      write (text,18981) pr
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 18
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.13_dp*0.7_dp
!
      if (abs(vcurfb(1)).gt.1.e-06_dp) then
         write (text,18985)  (vcurfb(i),i=1,3)
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = text
         imes(msg) = 28
         xpos(msg) = xabs
         ypos(msg) = yabs
         yabs = yabs - dyabs
         ht(msg) = 0.13_dp*0.7_dp
      endif
      rot_form: if (kvtor.ge.1.and.kvtor.le.3) then
      rot_exists: if (icurrt.eq.5.and. abs(brsp(nfsum+1)).gt.1.e-10_dp) then
         do i=nfnpcr+1,nfnpcr+kwwcur,2
            xxnorm=brsp(i)/brsp(nfsum+1)
            if (i+1.gt.nfnpcr+kwwcur) then
               if (i.eq.nfnpcr+1) then
                  write (text,18971) xxnorm
               else
                  write (text,18973) xxnorm
               endif
            else
               xynorm=brsp(i+1)/brsp(nfsum+1)
               if (i.eq.nfnpcr+1) then
                  write (text,28971) xxnorm,xynorm
               else
                  write (text,28973) xxnorm,xynorm
               endif
            endif
            msg = msg + 1
            note(msg) = 1
            lmes(msg) = text
            imes(msg) = 31
            xpos(msg) = xabs
            ypos(msg) = yabs
            yabs = yabs - dyabs
            ht(msg) = 0.13_dp*0.7_dp
         enddo
      endif rot_exists
      endif rot_form
!
      xabs= 1.2_dp
      yabs= -0.8_dp
      if ((icurrt.eq.2.or.icurrt.eq.5).and. &
                   abs(brsp(nfsum+1)).gt.1.e-10_dp) then
         do i=nfsum+1+kppcur,nfsum+kppcur+kffcur
            xxnorm=brsp(i)/brsp(nfsum+1)
            write (text,18970) xxnorm
            msg = msg + 1
            note(msg) = 1
            lmes(msg) = text
            imes(msg) = 14
            xpos(msg) = xabs
            ypos(msg) = yabs
            yabs = yabs - dyabs
            ht(msg) = 0.13_dp*0.7_dp
         enddo
      endif
!
      pl_files_14: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      ibrdr = 1
      xphy = 4.2_dp
      yphy = 2.80_dp
!     yphy = 3.8_dp
      hight = 0.12_dp
      nplen = 100
      nxlen = 100
      nylen = 100
      ixtck = 1
      iytck = 2
      intax = 1
      igridx = 1
      igridy = 2
      idot = 1
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy, &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, xmm, xmm, &
      xiter(1), xdel, xiter(nw), yorg, ystp, ynmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      igridx, igridy, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx, yy, nxy, ncnct, &
      icont, nword, zz, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_14
      endif is_vacuum_4
!      call endpl(0)

      xmm=1.6_dp
      endif eq_file_mode
!-----------------------------------------------------------------------
!-    plot experimental vs calculated magnetic signals                --
!-----------------------------------------------------------------------
      B_signals: if (iexcal.gt.0) then
      do kk=1,nsilop
         expsi(kk)=silopt(jtime,kk)
         si(kk)=kk
      enddo
      do kk=1,magpri
         expmp(kk)=expmpi(jtime,kk)
         rmpi(kk)=kk
      enddo
!-----------------------------------------------------------------------
!--   limits of psi loops & probes                                    --
!-----------------------------------------------------------------------
      cslmin=1.0e30_dp
      cslmax=-1.0e30_dp
      silmin=1.0e30_dp
      silmax=-1.0e30_dp
      z000=0.0
      do i=1,nsilop
         cslmin=min(cslmin,csilop(i,jtime))
         cslmax=max(cslmax,csilop(i,jtime))
         if (fwtsi(i).gt.0.0) then
           silmin=min(silmin,expsi(i))
           silmax=max(silmax,expsi(i))
         else
           silmin=min(silmin,z000)
           silmax=max(silmax,z000)
         endif
      enddo
      cslmin=min(cslmin,silmin)
      cslmax=max(cslmax,silmax)
!
      cmpmin=1.0e30_dp
      cmpmax=-1.0e30_dp
      empmin=1.0e30_dp
      empmax=-1.0e30_dp
      do i=1,magpri
         cmpmin=min(cmpmin,cmpr2(i,jtime))
         cmpmax=max(cmpmax,cmpr2(i,jtime))
         if (fwtmp2(i).gt.0.0) then
           empmin=min(empmin,expmp(i))
           empmax=max(empmax,expmp(i))
         else
           empmin=min(empmin,z000)
           empmax=max(empmax,z000)
         endif
      enddo
      cmpmin=min(cmpmin,empmin)
      cmpmax=max(cmpmax,empmax)
!
      pl_files_15: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      nn = nn + 1
      sclpc(nn) = 0.5_dp
      nxy(nn) = nsilop
      clearx(nn)='CYAN'
      do ii = 1, nxy(nn)
         xx(ii,nn) = si(ii)
         yy(ii,nn) = csilop(ii,jtime)
      enddo
      do i=1,nsilop
         write (96,*) si(i),csilop(i,jtime)
      enddo
      nn = nn + 1
      markme(nn) = 18
      ncnct(nn) = -1
      sclpc(nn) = 0.5_dp
      nxy(nn) = nsilop
      clearx(nn)='PINK'
      do ii = 1, nxy(nn)
         xx(ii,nn) = si(ii)
         if (fwtsi(ii).gt.0.0) then
           yy(ii,nn) = expsi(ii)
         else
           yy(ii,nn) = 0.0
         endif
      enddo
      do i=1,nsilop
         write (96,*) si(i),expsi(i)
      enddo
      ibrdr = 1
      xphy = 4.7_dp
      yphy = 1.0
      hight = 0.12_dp
      nplen = 100
      nxlen = 100
      nylen = 100
      xlen = 6.0
      ylen = 3.0
      xstp = 0.0
      ystp = 0.0
      intax = 1
      igridx = 1
      igridy = 1
      idot = 1
      xtitle = 'psi loops$'
      ytitle = 'exp vs calc (v-s/r)$'
      ncurve = nn
      iexit = 1
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy, &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, xlen, ylen, &
      si(1), xstp, si(nsilop), cslmin, ystp, cslmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      igridx, igridy, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx, yy, nxy, ncnct, &
      icont, nword, zz, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
!-----------------------------------------------------------------------
!     plot exp vs calc magnetic probes
!
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      nn = nn + 1
      sclpc(nn) = 0.5_dp
      nxy(nn) = magpri
      clearx(nn)='CYAN'
      do ii = 1, nxy(nn)
         xx(ii,nn) = rmpi(ii)
         yy(ii,nn) = cmpr2(ii,jtime)
      enddo

      nn = nn + 1
      markme(nn) = 18
      ncnct(nn) = -1
      sclpc(nn) = 0.5_dp
      clearx(nn)='PINK'
      nxy(nn) = magpri
      do ii = 1, nxy(nn)
         xx(ii,nn) = rmpi(ii)
         if (fwtmp2(ii).gt.0.0) then
           yy(ii,nn) = expmp(ii)
         else
           yy(ii,nn) = 0.0
         endif
      enddo

      do i=1,magpri
         write (96,*) rmpi(i),cmpr2(i,jtime)
      enddo
      do i=1,magpri
         write (96,*) rmpi(i),expmp(i)
      enddo
      xabs=-4.3_dp
      yabs=3.0
      dyabs = 0.22_dp
      write(text,8950) trim(ch1),trim(ch2),efitvers
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 27
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.12_dp
      write (text,8960) uday
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.12_dp
      write (text,9000) ishot
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.12_dp
      write (text,9020) itime,itimeu
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.12_dp
!-----------------------------------------------------------------------
!--   print out pointnames                                            --
!-----------------------------------------------------------------------
      yabs=1.5_dp
!     call messag('       probe identification$',100,xabs,yabs)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = '       probe identification$'
      imes(msg) = 100
      xpos(msg) = xabs
      ypos(msg) = yabs
      ht(msg) = 0.10_dp
!     call messag('no.      name        no.     name$',100,xabs,
!     .yabs-.25)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = 'no.      name        no.     name$'
      imes(msg) = 100
      xpos(msg) = xabs
      ypos(msg) = yabs-.25_dp
      ht(msg) = 0.10_dp
      yabs=1.0
      do i=magpri67+1,magpri67+16
!        call intno(i,xabs,yabs)
         msg = msg + 1
         note(msg) = 5
         inum(msg) = i
         xpos(msg) = xabs
         ypos(msg) = yabs
         ht(msg) = 0.07_dp

!        call messag(%ref(mpnam2(i)),10,-3.8_dp,yabs)
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = (mpnam2(i))
         imes(msg) = 10
         xpos(msg) = -3.8_dp
         ypos(msg) = yabs
         ht(msg) = 0.07_dp
         yabs=yabs-.13_dp
      enddo
      xabs=-2.6_dp
      yabs=1.0
      do i=magpri67+17,magpri67+magpri322
!        call intno(i,xabs,yabs)
         msg = msg + 1
         note(msg) = 5
         inum(msg) = i
         xpos(msg) = xabs
         ypos(msg) = yabs
         ht(msg) = 0.07_dp
!        call messag(%ref(mpnam2(i)),10,-2.0,yabs)
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = (mpnam2(i))
         imes(msg) = 10
         xpos(msg) = -2.0
         ypos(msg) = yabs
         ht(msg) = 0.07_dp
         yabs=yabs-.13_dp
      enddo
!
      xabs=-4.3_dp
      yabs=yabs-.5_dp
!     call height(.10_dp)
!     call messag('     psi loop identification$',100,xabs,yabs)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = '     psi loop identification$'
      imes(msg) = 100
      xpos(msg) = xabs
      ypos(msg) = yabs
      ht(msg) = 0.10_dp
!     call messag('no       name        no.     name$',100,xabs,
!     .yabs-.25_dp)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = 'no       name        no.     name$'
      imes(msg) = 100
      xpos(msg) = xabs
      ypos(msg) = yabs - 0.25_dp
      ht(msg) = 0.10_dp
!     call height(.07_dp)
      yabs=yabs-.37_dp
      yyabs=yabs
      do i=1,21
         yabs=yabs-.13_dp
!        call intno(i,xabs,yabs)
         msg = msg + 1
         note(msg) = 5
         inum(msg) = i
         xpos(msg) = xabs
         ypos(msg) = yabs
         ht(msg) = 0.07_dp
!        call messag(%ref(lpname(i)),10,-3.8_dp,yabs)
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = (lpname(i))
         imes(msg) = 10
         xpos(msg) = -3.8_dp
         ypos(msg) = yabs
         ht(msg) = 0.07_dp
      enddo
      yabs=yyabs
      do i=22,nsilop
         yabs=yabs-.13_dp
!        call intno(i,-2.6,yabs)
         msg = msg + 1
         note(msg) = 5
         inum(msg) = i
         xpos(msg) = -2.6_dp
         ypos(msg) = yabs
         ht(msg) = 0.07_dp
!        call messag(%ref(lpname(i)),10,-2.0,yabs)
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = (lpname(i))
         imes(msg) = 10
         xpos(msg) = -2.0
         ypos(msg) = yabs
         ht(msg) = 0.07_dp
      enddo
      ibrdr = 1
      xphy = 4.7_dp
      yphy = 5.0
      hight = 0.12_dp
      nplen = 100
      nxlen = 100
      nylen = 100
      xlen = 6.0
      ylen = 3.0
      xstp = 0.0
      ystp = 0.0
      intax = 1
      igridx = 1
      igridy = 1
      idot=1
      xtitle = 'magnetic probes$'
      ytitle = 'exp vs calc (T)$'
      ncurve = nn
      iexit = 2
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy, &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, xlen, ylen, &
      rmpi(1), xstp, rmpi(magpri), cmpmin, ystp, cmpmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      igridx, igridy, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx, yy, nxy, ncnct, &
      icont, nword, zz, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_15
      endif B_signals
!
      is_vacuum_5: if (ivacum.eq.0) then
      plot_extra: if (kwripre.gt.0) then
         dataname=dataname(1:lprx)//'_presf'
         open(unit=62,file=dataname,status='old',iostat=ioerr)
         if (ioerr.eq.0) close(unit=62,status='delete')
         open(unit=62,file=dataname,status='new')
         nnww=nw
         xdum=0.
         delvn=1.0_dp/(nw-1)
         do i=1,nw
            workb(i)=sqrt(volp(i)/volp(nw))
            voln(i)=(i-1)*delvn
         enddo
         call zpline(nw,workb,pres,bvoln,cvoln,dvoln)
         voln(1)=0.0
         voln(nw)=1.0
         do i=1,nw
            presnow=seval(nw,voln(i),workb,pres,bvoln,cvoln,dvoln)
            write (62,*) voln(i),presnow,xdum,xdum
            presv(i)=presnow
         enddo
         close(unit=62)
         if (nomegat.gt.0.and.kprfit.eq.3) then
           call zpline(nw,voln,workb,bvoln,cvoln,dvoln)
           do i=1,nomegat
            yxn=rpresws(i)
            rpreswv(i)=seval(nw,yxn,voln,workb,bvoln,cvoln,dvoln)
           enddo
         endif
         mmpress=mpress
         dataname=dataname(1:lprx)//'_qthom'
         open(unit=62,file=dataname,status='old',iostat=ioerr)
         if (ioerr.eq.0) close(unit=62,status='delete')
         open(unit=62,file=dataname,status='new')
         do i=1,mmpress
            write (62,*) zqthom(i),qthom(i),xdum,xdum
         enddo
         close(unit=62)
         dataname=dataname(1:lprx)//'_jscrap'
         open(unit=62,file=dataname,status='old',iostat=ioerr)
         if (ioerr.eq.0) close(unit=62,status='delete')
         open(unit=62,file=dataname,status='new')
         do i=1,nw
            write (62,*) rscrap(i),curscr(i),xdum,xdum
         enddo
         close(unit=62)
         dataname=dataname(1:lprx)//'_qpsi'
         open(unit=62,file=dataname,status='old',iostat=ioerr)
         if (ioerr.eq.0) close(unit=62,status='delete')
         open(unit=62,file=dataname,status='new')
         do i=1,nw
            write (62,*) workb(i),qpsi(i),xdum,xdum
         enddo
         close(unit=62)
         dataname=dataname(1:lprx)//'_jor'
         open(unit=62,file=dataname,status='old',iostat=ioerr)
         if (ioerr.eq.0) close(unit=62,status='delete')
         open(unit=62,file=dataname,status='new')
         do i=1,nw
            cjorka=cjor(i)/1000.
            write (62,*) workb(i),cjorka,xdum,xdum
         enddo
         close(unit=62)
         dataname=dataname(1:lprx)//'_jorec'
         open(unit=62,file=dataname,status='old',iostat=ioerr)
         if (ioerr.eq.0) close(unit=62,status='delete')
         open(unit=62,file=dataname,status='new')
         do i=1,nw
            cjorka=cjorec(i)/1000.
            write (62,*) workb(i),cjorka,xdum,xdum
         enddo
         close(unit=62)
!-----------------------------------------------------------------------
!--      shear   dqdx/q                                               --
!-----------------------------------------------------------------------
         dataname=dataname(1:lprx)//'_shear'
         open(unit=62,file=dataname,status='old',iostat=ioerr)
         if (ioerr.eq.0) close(unit=62,status='delete')
         open(unit=62,file=dataname,status='new')
         call zpline(nw,workb,qpsi,bvoln,cvoln,dvoln)
         do i=1,nw
            qqqnow=seval(nw,voln(i),workb,qpsi,bvoln,cvoln,dvoln)
            qshear=speval(nw,voln(i),workb,qpsi,bvoln,cvoln,dvoln)/ &
                   qqqnow
            write (62,*) voln(i),qshear,xdum,xdum
         enddo
         close(unit=62)
!-----------------------------------------------------------------------
!--      shear ballooning   x*dqdx/q                                  --
!-----------------------------------------------------------------------
         dataname=dataname(1:lprx)//'_shearbal'
         open(unit=62,file=dataname,status='old',iostat=ioerr)
         if (ioerr.eq.0) close(unit=62,status='delete')
         open(unit=62,file=dataname,status='new')
         call zpline(nw,workb,qpsi,bvoln,cvoln,dvoln)
         do i=1,nw
            qqqnow=seval(nw,voln(i),workb,qpsi,bvoln,cvoln,dvoln)
            qshear=speval(nw,voln(i),workb,qpsi,bvoln,cvoln,dvoln)/ &
                   qqqnow*voln(i)
            write (62,*) voln(i),qshear,xdum,xdum
         enddo
         close(unit=62)
!-----------------------------------------------------------------------
!--      avearage poloidal magnetic field in T                        --
!-----------------------------------------------------------------------
         dataname=dataname(1:lprx)//'_bpol'
         open(unit=62,file=dataname,status='old',iostat=ioerr)
         if (ioerr.eq.0) close(unit=62,status='delete')
         open(unit=62,file=dataname,status='new')
         do i=1,nw
            qshear=sqrt(volp(i)/volp(nw))
            write (62,*) qshear,bpolss(i),xdum,xdum
         enddo
         close(unit=62)
!-----------------------------------------------------------------------
!--      avearage poloidal magnetic field in T versus R at Z=0        --
!-----------------------------------------------------------------------
         dataname=dataname(1:lprx)//'_bpolrmaj'
         open(unit=62,file=dataname,status='old',err=31430)
         close(unit=62,status='delete')
31430    continue
         open(unit=62,file=dataname,status='new')
         do i=1,nw
            write (62,*) rmajz0(i),bpolss(i),xdum,xdum
         enddo
         close(unit=62)
!-----------------------------------------------------------------------
!--      R at Z=0                                                     --
!-----------------------------------------------------------------------
         dataname=dataname(1:lprx)//'_rmajz0'
         open(unit=62,file=dataname,status='old',iostat=ioerr)
         if (ioerr.eq.0) close(unit=62,status='delete')
         open(unit=62,file=dataname,status='new')
         do i=1,nw
            qshear=sqrt(volp(i)/volp(nw))
            write (62,*) qshear,rmajz0(i),xdum,xdum
         enddo
         close(unit=62)
!-----------------------------------------------------------------------
!--      shearbal/q2  11/01/92  Lao                                   --
!-----------------------------------------------------------------------
         dataname=dataname(1:lprx)//'_sq2'
         open(unit=62,file=dataname,status='old',iostat=ioerr)
         if (ioerr.eq.0) close(unit=62,status='delete')
         open(unit=62,file=dataname,status='new')
         do i=1,nw
            qqqnow=seval(nw,voln(i),workb,qpsi,bvoln,cvoln,dvoln)
            qshear=speval(nw,voln(i),workb,qpsi,bvoln,cvoln,dvoln)/ &
                   qqqnow*voln(i)
            qshear=qshear/qqqnow**2
            write (62,*) voln(i),qshear,xdum,xdum
         enddo
         close(unit=62)
!---------------------------------------------------------------------
!--      angular speed                                              --
!---------------------------------------------------------------------
         if (nomegat.gt.0.and.kprfit.eq.3) then
            dataname=dataname(1:lprx)//'_omegafit'
            open(unit=62,file=dataname,status='old',iostat=ioerr)
            if (ioerr.eq.0) close(unit=62,status='delete')
            open(unit=62,file=dataname,status='new')
            do i=1,nomegat
               yxn=rpreswv(i)
               write (62,*) yxn,omegat(i)
            enddo
            close(unit=62)
          endif
      endif plot_extra
!-----------------------------------------------------------------------
!--   plots for kinetic fitting option                                --
!-----------------------------------------------------------------------
      kinetic: if ((kprfit.gt.0).and.(npress.gt.0)) then
      if (itek.eq.1) call tekall(4010,960,0,0,0)
      chsmin=1.0e+10_dp
      chsmax=-1.0e+10_dp
      pressbb=0.0
      dsi=1.0_dp/(nw-1)
      do i=1,nw
         workb(i)=sqrt(volp(i)/volp(nw))
         workc(i)=(pres(i)-pressbb)
         worke(i)=(i-1)*dsi
      enddo
      do i=1,npress
         workd(i)=(pressr(i)-pressbb)
         chsmin=min(chsmin,saipre(i))
         chsmax=max(chsmax,saipre(i))
         si(i)=i
      enddo
      call zpline(nw,worke,workb,bworkb,cworkb,dworkb)
!
      pl_files_16: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      nn = nn + 1
      sclpc(nn) = 0.5_dp
      nxy(nn) = npress
      ncnct(nn) = 1
      do ii = 1, nxy(nn)
         xx(ii,nn) = si(ii)
         yy(ii,nn) = saipre(ii)
      enddo
      ibrdr = 1
      xphy = 4.7_dp
      yphy = 1.0
      hight = 0.12_dp
      nplen = 100
      nxlen = 100
      nylen = 100
      xlen = 6.0
      ylen = 3.0
      xstp = (si(npress) - si(1))/4.
      ystp = (chsmax - chsmin)/2.
      ixtck = 1
      iytck = 2
      intax = 1
      intay = 0
      idot = 1
      xtitle = 'PRESSURE$'
      ytitle = 'CHI**2$'
      ncurve = nn
      iexit = 1
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy, &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, xlen, ylen, &
      si(1), xstp, si(npress), chsmin, ystp, chsmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      kgrid, kgrid, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx, yy, nxy, ncnct, &
      icont, nword, zz, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_16
!----------------------------------------------------------------------
!--   plot the pressure profile                                      --
!----------------------------------------------------------------------
      chsmin=1.e10_dp
      chsmax=-1.e10_dp
      do i=1,npress
         chsmin=min(workd(i),chsmin)
         chsmax=max(workd(i)+sigpre(i),chsmax)
         chsmax=max(chsmax,workc(i))
         saipre(i)=seval(nw,-rpress(i),worke,workb,bworkb,cworkb,dworkb)
      enddo
      chsmin=0.0
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      nn = nn + 1
      sclpc(nn) = 0.5_dp
      nxy(nn) = nw
      do ii = 1, nxy(nn)
         xx(ii,nn) = workb(ii)
         yy(ii,nn) = workc(ii)
      enddo
      nn = nn + 1
      sclpc(nn) = 0.5_dp
      markme(nn) = 15
      ncnct(nn) = -1
      nxy(nn) = npress
      do ii = 1, nxy(nn)
         xx(ii,nn) = saipre(ii)
         yy(ii,nn) = workd(ii)
      enddo
!-----------------------------------------------------------------------
!--   add error bars to pressure plot                                 --
!-----------------------------------------------------------------------
      do i=1,npress
         workd(1)=(pressr(i)-pressbb-sigpre(i))
         workd(2)=(pressr(i)-pressbb+sigpre(i))
         workc(1)=saipre(i)
         workc(2)=saipre(i)
         nn = nn + 1
         sclpc(nn) = 0.5_dp
         markme(nn) = 3
         nxy(nn) = 2
         ncnct(nn) = 1
         do ii = 1, nxy(nn)
            xx(ii,nn) = workc(ii)
            yy(ii,nn) = workd(ii)
         enddo
      enddo
      bimbf=0.0
      bimbe=0.0
!-----------------------------------------------------------------------
!--   write out QPLOT data files                                      --
!-----------------------------------------------------------------------
      plot_pres: if (kwripre.gt.0) then
      nnww=nw
      xdum=0.0
      kin: if (kprfit.eq.2.or.kprfit.eq.3) then
         dataname=dataname(1:lprx)//'_presd'
         open(unit=62,file=dataname,status='old',iostat=ioerr)
         if (ioerr.eq.0) close(unit=62,status='delete')
         open(unit=62,file=dataname,status='new')
         do i=1,npress
            write (62,*) saipre(i),pressr(i),xdum,sigpre(i)
         enddo
         close(unit=62)
      elseif (kprfit.eq.1) then kin
!-----------------------------------------------------------------------
!--      pressure input in flux space                                 --
!-----------------------------------------------------------------------
         do i=1,npress
            rpressv(i)=seval(nw,-rpress(i),worke,workb,bworkb, &
            cworkb,dworkb)
         enddo
         if (nbeam.gt.0) then
            do i=1,nbeam
               vnbeam(i)=seval(nw,sibeam(i),worke,workb,bworkb, &
               cworkb,dworkb)
            enddo
         endif
         if (npteth.gt.0) then
            do i=1,npteth
               call seva2d(bkx,lkx,bky,lky,c,rteth(i),zteth(i), &
                           pds,ier,n111)
               workc(i)=(pds(1)-simag)/(psibry-simag)
               workc(i)=seval(nw,workc(i),worke, &
                              workb,bworkb,cworkb,dworkb)
            enddo
         endif
         if (nption.gt.0) then
            do i=npteth+1,npteth+nption
               j=i-npteth
               call seva2d(bkx,lkx,bky,lky,c,rion(j),zion(j), &
                           pds,ier,n111)
                workc(i)=(pds(1)-simag)/(psibry-simag)
                workc(i)=seval(nw,workc(i),worke, &
                              workb,bworkb,cworkb,dworkb)
            enddo
         endif
         call zpline(npress,saipre,pressr,bworkb,cworkb,dworkb)
         dataname=dataname(1:lprx)//'_presd'
         open(unit=62,file=dataname,status='old',iostat=ioerr)
         if (ioerr.eq.0) close(unit=62,status='delete')
         open(unit=62,file=dataname,status='new')
         if (npteth+nption.gt.0) then
            ydumin=.03*seval(npress,workc(npteth+1), &
                             saipre,pressr,bworkb,cworkb,dworkb)
            do i=1,npteth+nption
               workd(i)=seval(npress,workc(i),saipre, &
               pressr,bworkb,cworkb,dworkb)
               ydum=.10_dp*workd(i)
               ydum=max(ydum,ydumin)
               write (62,*) workc(i),workd(i),xdum,ydum
            enddo
            close(unit=62)
         endif
!-----------------------------------------------------------------------
!--      write out pressure data points in flux space for KPRFIT=1    --
!-----------------------------------------------------------------------
         dataname=dataname(1:lprx)//'_presd'//'_psi'
         open(unit=62,file=dataname,status='old',iostat=ioerr)
         if (ioerr.eq.0) close(unit=62,status='delete')
         open(unit=62,file=dataname,status='new')
         do i=1,npress
            write (62,*) rpressv(i),pressr(i),xdum,sigpre(i)
         enddo
         close(unit=62)
!-----------------------------------------------------------------------
!--      now dpdr                                                     --
!-----------------------------------------------------------------------
         dataname=dataname(1:lprx)//'_prespf'//'_psi'
         call zpline(nw,voln,presv,bpresv,cpresv,dpresv)
         open(unit=62,file=dataname,status='old',iostat=ioerr)
         if (ioerr.eq.0) close(unit=62,status='delete')
         open(unit=62,file=dataname,status='new')
         do i=1,npress
            prespv=- &
            speval(nw,rpressv(i),voln,presv,bpresv,cpresv,dpresv)
            write (62,*) rpressv(i),prespv,xdum,xdum
         enddo
         close(unit=62)
!-----------------------------------------------------------------------
!--      now dpdr thermal                                             --
!-----------------------------------------------------------------------
         beams: if (nbeam.gt.0) then
            dataname=dataname(1:lprx)//'_pbeamf'
            open(unit=62,file=dataname,status='old',iostat=ioerr)
            if (ioerr.eq.0) close(unit=62,status='delete')
            open(unit=62,file=dataname,status='new')
            do i=1,npress
               ptherm(i)=pressr(i)-pbeam(i)
            enddo
            call zpline(npress,rpressv,ptherm,bpressv,cpressv,dpressv)
            do i=1,nw
              pbimf(i)=seval(npress,voln(i),rpressv, &
                             ptherm,bpressv,cpressv,dpressv)
              pbimf(i)=presv(i)-pbimf(i)
              write (62,*) voln(i),pbimf(i),xdum,xdum
            enddo
            close(unit=62)
!-----------------------------------------------------------------------
!--         get fitted beam beta                                      --
!-----------------------------------------------------------------------
            do i=1,nw-1
               delvol=(voln(i+1)**2-voln(i)**2)
               bimbf=bimbf+(pbimf(i)+pbimf(i+1))*delvol
            enddo
            bimbf=bimbf*0.5_dp*volume(jtime)/1.e6_dp*1.5_dp
            bimbf=bimbf*betat(jtime)/wmhd(jtime)
!-----------------------------------------------------------------------
!--         now beam data                                             --
!-----------------------------------------------------------------------
            dataname=dataname(1:lprx)//'_pbeam'
            open(unit=62,file=dataname,status='old',iostat=ioerr)
            if (ioerr.eq.0) close(unit=62,status='delete')
            open(unit=62,file=dataname,status='new')
            call zpline(nbeam,sibeam,pbeam,bpressv,cpressv,dpressv)
            call zpline(nw,workb,voln,bvoln,cvoln,dvoln)
            do i=1,nw
              sivol(i)=seval(nw,voln(i),workb, &
                             voln,bvoln,cvoln,dvoln)
              xn=sivol(i)
              workc(i)=seval(nbeam,xn,sibeam,pbeam, &
                            bpressv,cpressv,dpressv)
            enddo
            xdum=0.0
            do i=1,nw
               write (62,*) voln(i),workc(i),xdum,xdum
            enddo
            close(unit=62)
!-----------------------------------------------------------------------
!--         plot beam pressures                                       --
!-----------------------------------------------------------------------
            nn = nn + 1
            ndotme(nn) = 1
            nxy(nn) = nw
            do ii = 1, nxy(nn)
               xx(ii,nn) = voln(ii)
               yy(ii,nn) = pbimf(ii)
            enddo
            nn = nn + 1
            ndshme(nn) = 1
            nxy(nn) = nw
            do ii = 1, nxy(nn)
               xx(ii,nn) = voln(ii)
               yy(ii,nn) = workc(ii)
            enddo
!-----------------------------------------------------------------------
!--         get input beam beta                                       --
!-----------------------------------------------------------------------
            do i=1,nw-1
               delvol=(voln(i+1)**2-voln(i)**2)
               bimbe=bimbe+(workc(i)+workc(i+1))*delvol
            enddo
            bimbe=bimbe*0.5_dp*volume(jtime)/1.e6_dp*1.5_dp
            bimbe=bimbe*betat(jtime)/wmhd(jtime)
         endif beams
      endif kin
      endif plot_pres
      xabs=-4.5_dp
      yabs=3.0
      dyabs = 0.22_dp
      write(text,8950) trim(ch1),trim(ch2),efitvers
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 27
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,8960) uday
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9000) ishot
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9020) itime,itimeu
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9497) chipre
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9555) chitot
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      pressbn=prbdry/pres(1)
      write (text,9560) pressbn
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      pressbn=pressb/pres(1)
      write (text,9562) pressbn
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9565) kppcur
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9570) kffcur
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9575) pcurbd
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9580) fcurbd
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9832) pres(1)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      dp1dxf=ppcurr(x111,kppcur)/darea*(psibry-simag)
      write (text,9835) dp1dxf
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9834) prespb
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9836) kpressb
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9838) cjor95(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9840) pp95(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9842) bimbf,bimbe
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
!
      pr=pres(1)/(.667_dp*wmhd(jtime)/(volume(jtime)/1.e6_dp))
      write (text,18981) pr
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 18
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
!
      pl_files_17: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      ibrdr = 1
      xphy = 4.7_dp
      yphy = 5.0
      hight = 0.12_dp
      nplen = 100
      nxlen = 100
      nylen = 100
      xlen = 6.0
      ylen = 3.0
      xstp = (workb(nw) - workb(1))/4.
      ystp = (chsmax - chsmin)/2.
      idot = 1
      ixtck = 1
      iytck = 2
      xtitle = 'SQRT(V/V(1))$'
      ytitle = 'P(J/m3)$'
      ncurve = nn
      iexit = 2
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy, &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, xlen, ylen, &
      workb(1), xstp, workb(nw), chsmin, ystp, chsmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      kgrid, kgrid, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx, yy, nxy, ncnct, &
      icont, nword, zz, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_17
      endif kinetic
!----------------------------------------------------------------------
!--   rotational page, first plot total pressure contour             --
!----------------------------------------------------------------------
      rotation: if (kvtor.gt.0) then
      if (itek.eq.1) call tekall(4010,960,0,0,0)
!-----------------------------------------------------------------------
!--   make Pt array into 2-d array                                    --
!-----------------------------------------------------------------------
      do ih=1,nh
         do iw=1,nw
            icnt=(iw-1)*nh+ih
            bfield(iw,ih)=presst(icnt)
         enddo
      enddo
      almin=rgrid(1)
      almax=rgrid(nw)
      blmin=zgrid(1)
      blmax=zgrid(nh)
      xll = 2.8_dp
      yll = xll*(blmax-blmin)/(almax-almin)
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      npplot=15
      if (iabs(kplotp).ge.5) then
        npplot=iabs(kplotp)
      endif
      siinc=(rmax-rmaxis)/npplot
!-----------------------------------------------------------------------
!--   plot pressure contours                                          --
!-----------------------------------------------------------------------
      do i=1,npplot-1
           nn=nn+1
           rppnow=rmax-siinc*i
           call seva2d(bwx,lwx,bwy,lwy,cw,rppnow,zmaxis, &
                       pds,ier,n111)
           spwant=pds(1)
           call surfac(spwant,presst,nw,nh,rgrid,zgrid,xplt,yplt, &
                       nplt,npoint,drgrid,dzgrid,rmin,rmax,zmin, &
                       zmax,n11,rmaxis,zmaxis,negcur,kerror)
           if (kerror.gt.0) return
           nxy(nn) = nplt
           ndshme(nn) = 1
           do ii = 1,nplt
              xx(ii,nn) = xplt(ii)
              yy(ii,nn) = yplt(ii)
           enddo
      enddo
!-------------------------------------------------------------------
!--   plot boundary                                               --
!-------------------------------------------------------------------
      nn = nn + 1
      nxy(nn) = nbnow
      clearx(nn) = 'PINK'
      do ii = 1,nxy(nn)
            xx(ii,nn) = xbnow(ii)
            yy(ii,nn) = ybnow(ii)
      enddo
!
      nn = nn + 1
      thcrv(nn) = .02
      nxy(nn) = limitr
      do ii = 1,nxy(nn)
            xx(ii,nn) = xlim(ii)
            yy(ii,nn) = ylim(ii)
      enddo
      call pltcol(nfcoil,rf,zf,wf,hf,af,af2,n11, &
         nn, xx, yy, nxy, msg, note, inum, xpos, ypos, ht, &
         nshd, sxx, syy, nsxy, sangle, sgap, ngaps)
      xphy = 7.625_dp
      yphy = 0.5_dp
      xabs=-7.0
      yabs=7.85_dp
      dyabs = 0.22_dp
      write(text,8950) trim(ch1),trim(ch2),efitvers
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 27
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,8948) ifname(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,8960) uday
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9000) ishot
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9020) itime,itimeu
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9022) chiprw
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9024) rvtor
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9026) kvtor
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
!
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = 'Pt contours (nt/m2)$'
      imes(msg) = 100
      xpos(msg) = 0.0
      ypos(msg) = -.35_dp
      ht(msg) = 0.16_dp
!
      if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!       Initialize plot parameters
!-----------------------------------------------------------------------
        ncurve = nn
        ibrdr = 1
        hight = 0.10_dp
        nplen = 100
        nxlen = 0
        nylen = 0
        xtitle = 'r(m)$'
        ytitle = 'z(m)$'
        igridx = -1
        grce  = -1.0
        iexit = 1
!-----------------------------------------------------------------------
!       Write Plot Parameters
!-----------------------------------------------------------------------
        call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy, &
        iorel, xorl, yorl, hight, bngle, bshft, &
        ptitle, nplen, xtitle, nxlen, ytitle, nylen, xll, yll, &
        almin, almax-almin, almax, blmin, blmax-blmin, blmax, &
        iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
        isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
        igridx, igridy, idash, idot, ichdsh, ichdot, &
        thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
        markme, clearx, mrc, tlen, nmrk, rat, &
        xx, yy, nxy, ncnct, &
        icont, nword, bfield , ix, iy, zinc, line, mode, &
        lbflg, ithk, ipri, nline, draw, &
        nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
        nvec, xfm, yfm, xto, yto, ivec, &
        msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif
!
      do i=1,nw
         worka(i)=bfield(i,nh/2+1)
      enddo
      cslmin=1.0e+30_dp
      cslmax=1.0e-30_dp
      do i=1,nw
         cslmin=min(worka(i),cslmin)
         cslmax=max(worka(i),cslmax)
      enddo
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      xphy = 7.625_dp
      yphy = 0.5_dp
      iorel = 1
      xorl = 0.0
      yorl = yll + 0.5_dp
      nn = nn + 1
      nxy(nn) = nw
      do ii = 1, nxy(nn)
         xx(ii,nn) = rgrid(ii)
         yy(ii,nn) = worka(ii)
      enddo
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = 'R(m)$'
      imes(msg) = 100
      xpos(msg) = 1.0
      ypos(msg) = -.25_dp
      ht(msg) = 0.10_dp
!
      if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!       Initialize plot parameters
!-----------------------------------------------------------------------
        ncurve = nn
        ibrdr = 1
        hight = 0.10_dp
        nplen = -100
        nxlen = 100
        nylen = 100
        ylen = xll*.67_dp
        xtitle = '$'
        ytitle = 'Pt(nt/m2)$'
        iexit = 1
        if (kprfit.ne.3) iexit=2
!-----------------------------------------------------------------------
!       Write Plot Parameters
!-----------------------------------------------------------------------
        call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy, &
        iorel, xorl, yorl, hight, bngle, bshft, &
        ptitle, nplen, xtitle, nxlen, ytitle, nylen, xll, ylen, &
        almin, almax-almin, almax, cslmin, cslmax-cslmin, cslmax, &
        iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
        isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
        igridx, igridy, idash, idot, ichdsh, ichdot, &
        thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
        markme, clearx, mrc, tlen, nmrk, rat, &
        xx, yy, nxy, ncnct, &
        icont, nword, bfield, ix, iy, zinc, line, mode, &
        lbflg, ithk, ipri, nline, draw, &
        nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
        nvec, xfm, yfm, xto, yto, ivec, &
        msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif
!----------------------------------------------------------------------
!--   plot Pw and fitted values                                      --
!----------------------------------------------------------------------
      rotational_pressure: if (kprfit.eq.3) then
      chsmin=1.0e+10_dp
      chsmax=-1.0e+10_dp
      pressbb=0.0
      dsi=1.0_dp/(nw-1)
      do i=1,nw
         workb(i)=sqrt(volp(i)/volp(nw))
         workc(i)=(pressw(i)-preswb)
         worke(i)=(i-1)*dsi
      enddo
      do i=1,npresw
         workd(i)=(presw(i)-preswb)
         chsmin=min(chsmin,saiprw(i))
         chsmax=max(chsmax,saiprw(i))
         si(i)=i
      enddo
      chsmax=max(x111,chsmax)
      call zpline(nw,worke,workb,bworkb,cworkb,dworkb)
!
      pl_files_18: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      nn = nn + 1
      sclpc(nn) = 0.5_dp
      nxy(nn) = npresw
      ncnct(nn) = 1
      do ii = 1, nxy(nn)
         xx(ii,nn) = si(ii)
         yy(ii,nn) = saiprw(ii)
      enddo
      ibrdr = 1
      xphy = 4.0
      yphy = 1.0
      hight = 0.12_dp
      nplen = 100
      nxlen = 100
      nylen = 100
      xlen = 2.8_dp
      ylen = 2.8_dp
      xstp = 0.0
      ystp = 0.0
      intax = 1
      idot = 1
      xtitle = 'PRESSW$'
      ytitle = 'CHI**2$'
      ncurve = nn
      iexit = 1
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy, &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, xlen, ylen, &
      si(1), xstp, si(npresw), chsmin, ystp, chsmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      kgrid, kgrid, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx, yy, nxy, ncnct, &
      icont, nword, zz, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_18
!----------------------------------------------------------------------
!--   plot the pressure profile                                      --
!----------------------------------------------------------------------
      chsmin=1.e10_dp
      chsmax=-1.e10_dp
      do i=1,npresw
         chsmin=min(workd(i),chsmin)
         chsmax=max(workd(i)+sigprw(i),chsmax)
         chsmax=max(chsmax,workc(i))
         saiprw(i)=seval(nw,-rpresw(i),worke,workb,bworkb,cworkb,dworkb)
      enddo
      chsmin=0.0
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      nn = nn + 1
      sclpc(nn) = 0.5_dp
      nxy(nn) = nw
      do ii = 1, nxy(nn)
         xx(ii,nn) = workb(ii)
         yy(ii,nn) = workc(ii)
      enddo
      nn = nn + 1
      sclpc(nn) = 0.5_dp
      markme(nn) = 15
      ncnct(nn) = -1
      nxy(nn) = npresw
      do ii = 1, nxy(nn)
         xx(ii,nn) = saiprw(ii)
         yy(ii,nn) = workd(ii)
      enddo
!-----------------------------------------------------------------------
!--   add error bars to pressure plot                                 --
!-----------------------------------------------------------------------
      do i=1,npresw
         workd(1)=(presw(i)-preswb-sigprw(i))
         workd(2)=(presw(i)-preswb+sigprw(i))
         workc(1)=saiprw(i)
         workc(2)=saiprw(i)
         nn = nn + 1
         sclpc(nn) = 0.5_dp
         markme(nn) = 3
         nxy(nn) = 2
         ncnct(nn) = 1
         do ii = 1, nxy(nn)
            xx(ii,nn) = workc(ii)
            yy(ii,nn) = workd(ii)
         enddo
      enddo
!-----------------------------------------------------------------------
!--   write out QPLOT data files                                      --
!-----------------------------------------------------------------------
      if (kwripre.gt.0) then
        nnww=nw
        xdum=0.0
        dataname=dataname(1:lprx)//'_preswd'
        open(unit=62,file=dataname,status='old',iostat=ioerr)
        if (ioerr.eq.0) close(unit=62,status='delete')
        open(unit=62,file=dataname,status='new')
        do i=1,npresw
          write (62,*) saiprw(i),presw(i),xdum,sigprw(i)
        enddo
        close(unit=62)
        dataname=dataname(1:lprx)//'_presw'
        open(unit=62,file=dataname,status='old',iostat=ioerr)
        if (ioerr.eq.0) close(unit=62,status='delete')
        open(unit=62,file=dataname,status='new')
        do i=1,nw
          write (62,*) workb(i),pressw(i),xdum,xdum
        enddo
        close(unit=62)
      endif
!
      pl_files_19: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      ibrdr = 1
      xphy = 4.0
      yphy = 5.0
      hight = 0.12_dp
      nplen = 100
      nxlen = 100
      nylen = 100
      xlen = 2.8_dp
      ylen = 2.8_dp
      xstp = 0.0
      ystp = 0.0
      idot = 1
      xtitle = 'SQRT(V/V(1))$'
      ytitle = 'Pw(J/m3)$'
      ncurve = nn
      iexit = 2
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy, &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, xlen, ylen, &
      workb(1), xstp, workb(nw), chsmin, ystp, chsmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      kgrid, kgrid, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx, yy, nxy, ncnct, &
      icont, nword, zz, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_19
      endif rotational_pressure
      endif rotation
!----------------------------------------------------------------------
!--   end rotational page                                            --
!----------------------------------------------------------------------
      raw_kinetic: if (kprfit.eq.2) then
      if (itek.eq.1) call tekall(4010,960,0,0,0)
      chsmin=1.0e+10_dp
      chsmax=-1.0e+10_dp
      if (nptef.gt.0) then
         do i=1,nw
            xn=worke(i)
            workc(i)=tefit(1)
            do j=2,nptef
               workc(i)=workc(i)+tefit(j)*xn**(j-1)
            enddo
            chsmin=min(chsmin,workc(i))
            chsmax=max(chsmax,workc(i))
         enddo
      endif
      rsimin=saipre(1)
      do i=1,npress
         rsimin=min(rsimin,saipre(i))
         chsmin=min(chsmin,tethom(i))
         chsmax=max(chsmax,tethom(i)+sgteth(i))
      enddo
      do i=1,nw
         if (workb(i).ge.rsimin) exit
      enddo
      iiim=i
!
      ibrdr = 1
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      if (nptef.gt.0) then
         nn = nn + 1
         sclpc(nn) = 0.5_dp
         nxy(nn) = nw
         do ii = 1, nxy(nn)
            xx(ii,nn) = workb(ii)
            yy(ii,nn) = workc(ii)
         enddo
      endif
      ndodo=npress
      if (kpressb.eq.1) ndodo=npress-1
      nn = nn + 1
      ncnct(nn) = -1
      sclpc(nn) = 0.5_dp
      nxy(nn) = ndodo
      do ii = 1, nxy(nn)
         xx(ii,nn) = saipre(ii)
         yy(ii,nn) = tethom(ii)
      enddo
      if (kwripre.eq.1) then
         xdum=0.0
         write (87,*) npress
         do i=1,npress
            write (87,*) saipre(i),tethom(i),xdum,sgteth(i)
         enddo
         write (87,*) nnww
         do i=1,nw
            write (87,*) workb(i),workc(i),xdum,xdum
         enddo
      endif
      do i=1,ndodo
         workd(1)=(tethom(i)-sgteth(i))
         workd(2)=(tethom(i)+sgteth(i))
         workc(1)=saipre(i)
         workc(2)=saipre(i)
         nn = nn + 1
         sclpc(nn) = 0.5_dp
         markme(nn) = 3
         ncnct(nn) = 1
         nxy(nn) = 2
         do ii = 1, nxy(nn)
            xx(ii,nn) = workc(ii)
            yy(ii,nn) = workd(ii)
         enddo
      enddo
!
      pl_files_20: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      ibrdr = 1
      xphy = 4.7_dp
      yphy = 1.0
      hight = 0.12_dp
      nplen = 100
      nxlen = 100
      nylen = 100
      xlen = 6.0
      ylen = 3.0
      xstp = 0.0
      ystp = 0.0
      xtitle = 'SQRT(V/V(1))$'
      ytitle = 'Te(KeV)$'
      ncurve = nn
      iexit = 1
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy, &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, xlen, ylen, &
      workb(1), xstp, workb(nw), chsmin, ystp, chsmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      kgrid, kgrid, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx, yy, nxy, ncnct, &
      icont, nword, zz, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_20
!----------------------------------------------------------------------
!--   plot the ne profile                                            --
!----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      chsmin=0.0
      chsmax=1.0
      if (npnef.gt.0) then
         do i=1,nw
            xn=worke(i)
            workc(i)=defit(1)
            do j=2,npnef
               workc(i)=workc(i)+defit(j)*xn**(j-1)
            enddo
            chsmin=min(workc(i),chsmin)
            chsmax=max(workc(i),chsmax)
         enddo
      endif
      do i=1,ndodo
         chsmax=max(dnethom(i)+sgneth(i),chsmax)
      enddo
      chsmin=0.0
      if (npnef.gt.0) then
         nn = nn + 1
         sclpc(nn) = 0.5_dp
         nxy(nn) = nw
         do ii = 1, nxy(nn)
            xx(ii,nn) = workb(ii)
            yy(ii,nn) = workc(ii)
         enddo
      endif

      nn = nn + 1
      sclpc(nn) = 0.5_dp
      ncnct(nn) = -1
      nxy(nn) = ndodo
      do ii = 1, nxy(nn)
         xx(ii,nn) = saipre(ii)
         yy(ii,nn) = dnethom(ii)
      enddo

      if (kwripre.eq.1) then
         xdum=0.0
         write (87,*) npress
         do i=1,npress
            write (87,*) saipre(i),dnethom(i),xdum,sgneth(i)
         enddo
         write (87,*) nnww
         do i=1,nw
            write (87,*) workb(i),workc(i),xdum,xdum
         enddo
      endif
      do i=1,ndodo
         workd(1)=(dnethom(i)-sgneth(i))
         workd(2)=(dnethom(i)+sgneth(i))
         workc(1)=saipre(i)
         workc(2)=saipre(i)
         nn = nn + 1
         sclpc(nn) = 0.5_dp
         markme(nn) = 3
         ncnct(nn) = 1
         nxy(nn) = 2
         do ii = 1, nxy(nn)
            xx(ii,nn) = workc(ii)
            yy(ii,nn) = workd(ii)
         enddo
      enddo
!
      xabs=-4.5_dp
      yabs=3.0
      dyabs = 0.22_dp
      write(text,8950) trim(ch1),trim(ch2),efitvers
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 27
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,8960) uday
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9000) ishot
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9020) itime,itimeu
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9800) chisqte
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9805) chisqne
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9810) fco2ne
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9815) nptef
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9820) npnef
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
!
      pl_files_21: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      ibrdr = 1
      xphy = 4.7_dp
      yphy = 5.0
      hight = 0.12_dp
      nplen = 100
      nxlen = 100
      nylen = 100
      xlen = 6.0
      ylen = 3.0
      xstp = 0.0
      ystp = 0.0
      xtitle = 'SQRT(V/V(1))$'
      ytitle = 'ne(cm-3)$'
      ncurve = nn
      iexit = 2
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy, &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, xlen, ylen, &
      workb(1), xstp, workb(nw), chsmin, ystp, chsmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      kgrid, kgrid, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx, yy, nxy, ncnct, &
      icont, nword, zz, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_21
!
      if (itek.eq.1) call tekall(4010,960,0,0,0)
      chsmin=1.0e+10_dp
      chsmax=-1.0e+10_dp
      if (nptionf.ge.50) &
         call zpline(nption,xsiion,tionex,bwork,cwork,dwork)
      do i=1,nw
         xn=worke(i)
         if (nptionf.lt.50) then
            workc(i)=tifit(1)
            do j=2,nptionf
               workc(i)=workc(i)+tifit(j)*xn**(j-1)
            enddo
         else
            workc(i)=seval(nption,xn,xsiion,tionex,bwork,cwork,dwork)
         endif
         chsmin=min(chsmin,workc(i))
         chsmax=max(chsmax,workc(i))
      enddo
      do i=1,nption+1
         chsmin=min(chsmin,tionex(i))
         chsmax=max(chsmax,tionex(i)+sigti(i))
      enddo
      chsmin=0.0
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      ibrdr = 1
      nn = nn + 1
      nxy(nn) = nw
      do ii = 1, nxy(nn)
         xx(ii,nn) = workb(ii)
         yy(ii,nn) = workc(ii)
      enddo
      do i=1,nption+1
         workj(i)=seval(nw,xsiion(i),worke,workb,bworkb,cworkb,dworkb)
      enddo
      nn = nn + 1
      sclpc(nn) = 0.5_dp
      ncnct(nn) = -1
      nxy(nn) = nption+1
      do ii = 1, nxy(nn)
         xx(ii,nn) = workj(ii)
         yy(ii,nn) = tionex(ii)
      enddo
      do i=1,nption+1
         workd(1)=(tionex(i)-sigti(i))
         workd(2)=(tionex(i)+sigti(i))
         workc(1)=workj(i)
         workc(2)=workj(i)
         nn = nn + 1
         sclpc(nn) = 0.5_dp
         markme(nn) = 3
         ncnct(nn) = 1
         nxy(nn) = 2
         do ii = 1, nxy(nn)
            xx(ii,nn) = workc(ii)
            yy(ii,nn) = workd(ii)
         enddo
      enddo
!
      pl_files_22: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      ibrdr = 1
      xphy = 4.7_dp
      yphy = 1.0
      hight = 0.14_dp
      nplen = 100
      nxlen = 100
      nylen = 100
      xlen = 6.0
      ylen = 3.0
      xstp = 0.0
      ystp = 0.0
      xtitle = 'SQRT(V/V(1))$'
      ytitle = 'Ti(KeV)$'
      ncurve = nn
      iexit = 1
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy, &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, xlen, ylen, &
      workb(1), xstp, workb(nw), chsmin, ystp, chsmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      kgrid, kgrid, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx, yy, nxy, ncnct, &
      icont, nword, zz, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_22
      if (kwripre.eq.1) then
         xdum=0.0
         write (87,*) nption
         do i=1,nption
            write (87,*) workj(i),tionex(i),xdum,sigti(i)
         enddo
         write (87,*) nnww
         do i=1,nw
            write (87,*) workb(i),workc(i),xdum,xdum
         enddo
      endif
!----------------------------------------------------------------------
!--   plot the beam profile                                          --
!----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      chsmin=0.0
      chsmax=0.0
      do i=1,npress
         workc(i)=pbimth(i)/(pres(1)-pressbb)
         chsmax=max(workc(i),chsmax)
      enddo
      chsmin=0.0
      nn = nn + 1
      sclpc(nn) = 0.5_dp
      ncnct(nn) = -1
      nxy(nn) = npress
      do ii = 1, nxy(nn)
         xx(ii,nn) = saipre(ii)
         yy(ii,nn) = workc(ii)
      enddo
!
      plot_beam: if (kwripre.eq.1) then
         if (nbeam.gt.0) then
            call zpline(nbeam,sibeam,pbeam,bwork,cwork,dwork)
            call zpline(nw,workb,voln,bvoln,cvoln,dvoln)
            do i=1,nw
               sivol(i)=seval(nw,voln(i),workb,voln,bvoln,cvoln,dvoln)
               xn=sivol(i)
               workc(i)=seval(nbeam,xn,sibeam,pbeam,bwork,cwork,dwork)
            enddo
         endif
         xdum=0.0
         dataname=dataname(1:lprx)//'_pbeam'
         open(unit=62,file=dataname,status='old',iostat=ioerr)
         if (ioerr.eq.0) close(unit=62,status='delete')
         open(unit=62,file=dataname,status='new')
         do i=1,nw
            write (62,*) voln(i),workc(i),xdum,xdum
         enddo
         close(unit=62)
!---------------------------------------------------------------------
!--      beam ion density                                                --
!---------------------------------------------------------------------
         if (nbeam.gt.0) then
            call zpline(nbeam,sibeam,dnbeam,bwork,cwork,dwork)
            do i=1,nw
               xn=sivol(i)
               workc(i)=seval(nbeam,xn,sibeam,dnbeam,bwork,cwork,dwork)
               workc(i)=workc(i)*1.e19_dp
            enddo
         endif
         xdum=0.0
         dataname=dataname(1:lprx)//'_nbeam'
         open(unit=62,file=dataname,status='old',iostat=ioerr)
         if (ioerr.eq.0) close(unit=62,status='delete')
         open(unit=62,file=dataname,status='new')
         do i=1,nw
            write (62,*) voln(i),workc(i),xdum,xdum
         enddo
         close(unit=62)
      endif plot_beam
!
      xabs=-4.5_dp
      yabs=3.0
      dyabs = 0.22_dp
      write(text,8950) trim(ch1),trim(ch2),efitvers
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 27
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,8960) uday
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9000) ishot
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9020) itime,itimeu
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9825) chisqti
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9830) nptionf
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
!
      pl_files_23: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      ibrdr = 1
      xphy = 4.7_dp
      yphy = 5.0
      hight = 0.14_dp
      nplen = 100
      nxlen = 100
      nylen = 100
      xlen = 6.0
      ylen = 3.0
      xstp = 0.0
      ystp = 0.0
      xtitle = 'SQRT(V/V(1))$'
      ytitle = 'rel Pb$'
      ncurve = nn
      iexit = 2
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy, &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, xlen, ylen, &
      workb(1), xstp, workb(nw), chsmin, ystp, chsmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      kgrid, kgrid, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx, yy, nxy, ncnct, &
      icont, nword, zz, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_23
!
      plotpr: if (kplotpr.ne.0) then
      if (itek.eq.1) call tekall(4010,960,0,0,0)
      chsmin=1.0e+10_dp
      chsmax=-1.0e+10_dp
      delz=(zuperts(jtime)-zlowerts)/(nw-1)*0.01_dp
      do i=1,nw
         worke(i)=zlowerts*0.01_dp+(i-1)*delz
         call seva2d(bkx,lkx,bky,lky,c,rmajts,worke(i),pds,ier,n111)
         xn=(pds(1)-simag)/(psibry-simag)
         workb(i)=xn
         workc(i)=tefit(1)
         do j=2,nptef
            workc(i)=workc(i)+tefit(j)*xn**(j-1)
         enddo
         chsmin=min(chsmin,workc(i))
         chsmax=max(chsmax,workc(i))
      enddo
      xxxmin=worke(1)
      xxxmax=worke(nw)
      do i=1,npteth
         chsmin=min(chsmin,tethom(i))
         chsmax=max(chsmax,tethom(i)+sgteth(i))
         xxxmin=min(xxxmin,zteth(i))
         xxxmax=max(xxxmax,zteth(i))
      enddo
      chsmin=0.0
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      ibrdr = 1
      nn = nn + 1
      nxy(nn) = nw
      do ii = 1, nxy(nn)
         xx(ii,nn) = worke(ii)
         yy(ii,nn) = workc(ii)
      enddo
      nn = nn + 1
      sclpc(nn) = 0.5_dp
      ncnct(nn) = -1
      nxy(nn) = npteth
      do ii = 1, nxy(nn)
         xx(ii,nn) = zteth(ii)
         yy(ii,nn) = tethom(ii)
      enddo
      do i=1,npteth
         workd(1)=(tethom(i)-sgteth(i))
         workd(2)=(tethom(i)+sgteth(i))
         workc(1)=zteth(i)
         workc(2)=zteth(i)
         nn = nn + 1
         sclpc(nn) = 0.5_dp
         markme(nn) = 3
         ncnct(nn) = 1
         nxy(nn) = 2
         do ii = 1, nxy(nn)
            xx(ii,nn) = workc(ii)
            yy(ii,nn) = workd(ii)
         enddo
      enddo
      workc(1)=zuperts(jtime)*0.01_dp
      workc(2)=workc(1)
      workd(1)=chsmin
      workd(2)=chsmax
      nn = nn + 1
      ndshme(nn) = 1
      nxy(nn) = 2
      do ii = 1, nxy(nn)
         xx(ii,nn) = workc(ii)
         yy(ii,nn) = workd(ii)
      enddo
!
      pl_files_24: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      ibrdr = 1
      xphy = 4.7_dp
      yphy = 1.0
      hight = 0.14_dp
      nplen = 100
      nxlen = 100
      nylen = 100
      xlen = 6.0
      ylen = 3.0
      xstp = 0.0
      ystp = 0.0
      xtitle = 'Z(M)$'
      ytitle = 'Te(KeV)$'
      ncurve = nn
      iexit = 1
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy, &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, xlen, ylen, &
      xxxmin, xstp, xxxmax, chsmin, ystp, chsmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      kgrid, kgrid, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx, yy, nxy, ncnct, &
      icont, nword, zz, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_24
!----------------------------------------------------------------------
!--   plot the ne profile                                            --
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      chsmin=0.0
      chsmax=1.0
      do i=1,nw
         xn=workb(i)
         workc(i)=defit(1)
         do j=2,npnef
            workc(i)=workc(i)+defit(j)*xn**(j-1)
         enddo
         chsmin=min(workc(i),chsmin)
         chsmax=max(workc(i),chsmax)
      enddo
      do i=1,npneth
         chsmax=max(dnethom(i)+sgneth(i),chsmax)
      enddo
      chsmin=0.0
      nn = nn + 1
      nxy(nn) = nw
      do ii = 1, nxy(nn)
         xx(ii,nn) = worke(ii)
         yy(ii,nn) = workc(ii)
      enddo
      nn = nn + 1
      sclpc(nn) = 0.5_dp
      ncnct(nn) = -1
      nxy(nn) = npneth
      do ii = 1, nxy(nn)
         xx(ii,nn) = zteth(ii)
         yy(ii,nn) = dnethom(ii)
      enddo
      do i=1,npneth
         workd(1)=(dnethom(i)-sgneth(i))
         workd(2)=(dnethom(i)+sgneth(i))
         workc(1)=zteth(i)
         workc(2)=zteth(i)
         nn = nn + 1
         sclpc(nn) = 0.5_dp
         markme(nn) = 3
         ncnct(nn) = 1
         nxy(nn) = 2
         do ii = 1, nxy(nn)
            xx(ii,nn) = workc(ii)
            yy(ii,nn) = workd(ii)
         enddo
      enddo
      workc(1)=zuperts(jtime)*0.01_dp
      workc(2)=workc(1)
      workd(1)=chsmin
      workd(2)=chsmax
      nn = nn + 1
      ndshme(nn) = 3
      nxy(nn) = 2
      do ii = 1, nxy(nn)
         xx(ii,nn) = workc(ii)
         yy(ii,nn) = workd(ii)
      enddo
!
      xabs=-4.5_dp
      yabs=3.0
      write(text,8950) trim(ch1),trim(ch2),efitvers
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 27
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,8960) uday
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9000) ishot
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9020) itime,itimeu
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9800) chisqte
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9805) chisqne
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9810) fco2ne
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9815) nptef
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9820) npnef
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9828) zuperts(jtime),rlibim(jtime)
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
!
      pl_files_25: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      ibrdr = 1
      xphy = 4.7_dp
      yphy = 5.0
      hight = 0.14_dp
      nplen = 100
      nxlen = 100
      nylen = 100
      xlen = 6.0
      ylen = 3.0
      xstp = 0.0
      ystp = 0.0
      xtitle = 'Z(M)$'
      ytitle = 'ne(cm-3)$'
      ncurve = nn
      iexit = 2
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy, &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, xlen, ylen, &
      xxxmin, xstp, xxxmax, chsmin, ystp, chsmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      kgrid, kgrid, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx, yy, nxy, ncnct, &
      icont, nword, zz, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_25
      endif plotpr
      endif raw_kinetic
!
      plot_fieldlines: if (nextra.lt.0) then
      do k=1,2
         fpxmin=fpxtra(1,1)
         fpxmax=fpxmin
         bpxmin=bpxtra(1,1)
         bpxmax=bpxmin
         flxmin=flxtra(1,1)
         flxmax=flxmin
         do j=1,iabs(nextra)
            jj=(k-1)*iabs(nextra)+j
            do i=1,npxtra(jj)
               fpxmin=min(fpxmin,fpxtra(i,jj))
               fpxmax=max(fpxmax,fpxtra(i,jj))
               bpxmin=min(bpxmin,bpxtra(i,jj))
               bpxmax=max(bpxmax,bpxtra(i,jj))
               flxmin=min(flxmin,flxtra(i,jj))
               flxmax=max(flxmax,flxtra(i,jj))
            enddo
         enddo
         dfpx=(fpxmax-fpxmin)/2.
         dbpx=(bpxmax-bpxmin)/2.
         dflx=(flxmax-flxmin)/2.
         if (itek.eq.1) call tekall(4010,960,0,0,0)
         ibrdr = 1
!-----------------------------------------------------------------------
!        Initialize plot parameters
!-----------------------------------------------------------------------
         call init2d
         do j=1,iabs(nextra)
            jj=(k-1)*iabs(nextra)+j
!           call mrscod(0.5_dp,j*2,ratray)
            nn = nn + 1
            nxy(nn) = npxtra(jj)
            do ii = 1, nxy(nn)
               xx(ii,nn) = fpxtra(1,jj+ii-1)
               yy(ii,nn) = bpxtra(1,jj+ii-1)
               rat(ii) = ratray(ii)
            enddo
            mrc(nn) = 1
            tlen(nn) = 0.5_dp
            nmrk(nn) = j*2
!           call rlint(jj,fpxtra(npxtra(jj),jj),bpxtra(npxtra(jj),jj))
            msg = msg + 1
            note(msg) = 6
            inum(msg) = jj
            xpos(msg) = fpxtra(npxtra(jj),jj)
            ypos(msg) = bpxtra(npxtra(jj),jj)
            ht(msg) = 0.14_dp
         enddo
!
         pl_files_26: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!        Initialize plot parameters
!-----------------------------------------------------------------------
         ibrdr = 1
         xphy = 1.5_dp
         yphy = 5.0
         hight = 0.14_dp
         nplen = 100
         nxlen = 100
         nylen = 100
         xlen = 6.0
         ylen = 3.0
         igridx = 2
         igridy = 1
         xtitle = 'Lpol(m)$'
         ytitle = 'bp(T)$'
         ncurve = nn
         iexit = 1
!-----------------------------------------------------------------------
!        Write Plot Parameters
!-----------------------------------------------------------------------
         call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy, &
         iorel, xorl, yorl, hight, bngle, bshft, &
         ptitle, nplen, xtitle, nxlen, ytitle, nylen, xlen, ylen, &
         xxxmin, xstp, xxxmax, chsmin, ystp, chsmax, &
         iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
         isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
         kgrid, kgrid, idash, idot, ichdsh, ichdot, &
         thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
         markme, clearx, mrc, tlen, nmrk, rat, &
         xx, yy, nxy, ncnct, &
         icont, nword, zz, ix, iy, zinc, line, mode, &
         lbflg, ithk, ipri, nline, draw, &
         nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
         nvec, xfm, yfm, xto, yto, ivec, &
         msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
         endif pl_files_26
!-----------------------------------------------------------------------
!        Initialize plot parameters
!-----------------------------------------------------------------------
         call init2d
         do j=1,iabs(nextra)
            jj=(k-1)*iabs(nextra)+j
!           call mrscod(0.5_dp,j*2,ratray)
            nn = nn + 1
            nxy(nn) = npxtra(jj)
            do ii = 1, nxy(nn)
               xx(ii,nn) = fpxtra(1,jj+ii-1)
               yy(ii,nn) = flxtra(1,jj+ii-1)
               rat(ii) = ratray(ii)
            enddo
            mrc(nn) = 1
            tlen(nn) = 0.5_dp
            nmrk(nn) = j*2
!            call rlint(jj,fpxtra(npxtra(jj),jj),flxtra(npxtra(jj),jj))
            msg = msg + 1
            note(msg) = 6
            inum(msg) = jj
            xpos(msg) = fpxtra(npxtra(jj),jj)
            ypos(msg) = flxtra(npxtra(jj),jj)
            ht(msg) = 0.14_dp
         enddo   
!
         pl_files_27: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!        Initialize plot parameters
!-----------------------------------------------------------------------
         ibrdr = 1
         xphy = 1.5_dp
         yphy = 5.0
         hight = 0.14_dp
         nplen = 100
         nxlen = 100
         nylen = 100
         xlen = 6.0
         ylen = 3.0
         igridx = 2
         igridy = 1
         xtitle = 'Lpol(m)$'
         ytitle = 'Ls(m)$'
         ncurve = nn
         iexit = 2
!-----------------------------------------------------------------------
!        Write Plot Parameters
!-----------------------------------------------------------------------
         call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy,  &
         iorel, xorl, yorl, hight, bngle, bshft, &
         ptitle, nplen, xtitle, nxlen, ytitle, nylen, x11x, y11y, &
         almin, xstp, almax, blmin, ystp, blmax, &
         iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
         isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
         igridx, igridy, idash, idot, ichdsh, ichdot, &
         thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
         markme, clearx, mrc, tlen, nmrk, rat, &
         xx, yy, nxy, ncnct, &
         icont, nword, zz, ix, iy, zinc, line, mode, &
         lbflg, ithk, ipri, nline, draw, &
         nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
         nvec, xfm, yfm, xto, yto, ivec, &
         msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
         endif pl_files_27
      enddo
      endif plot_fieldlines
      endif is_vacuum_5
!-----------------------------------------------------------------------
!--   plot psi contours                                               --
!-----------------------------------------------------------------------
      contour_plots: if (iconsi.gt.0) then
      if (itek.eq.1) call tekall(4010,960,0,0,0)
!-----------------------------------------------------------------------
!--   make psi array into 2-d array                                   --
!-----------------------------------------------------------------------
      do ih=1,nh
         do iw=1,nw
            icnt=(iw-1)*nh+ih
            bfield(iw,ih)=psi(icnt)*1000.
         enddo
      enddo
      almin=rgrid(1)
      almax=rgrid(nw)
      blmin=zgrid(1)
      blmax=zgrid(nh)
      xll = 2.8_dp
      yll = xll*(blmax-blmin)/(almax-almin)
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      nn = nn + 1
      icont(nn) =1
      nxy(nn) = 0
      nword = 5000
      line = 0
      mode = 'DASH'
      lbflg = 'NOLABELS'
      ithk = 1
      ipri = 1
      cslmin=1.0e+30_dp
      cslmax=1.0e-30_dp
      do i=1,nw
         do j=1,nh
            cslmin=min(bfield(i,j),cslmin)
            cslmax=max(bfield(i,j),cslmax)
         enddo
      enddo
      siinc=(cslmax-cslmin)/8.
      if (iconsi.ge.5) then
         siinc=(cslmax-cslmin)/iconsi
      endif
      if ((iconsi.eq.2).or.(iconsi.eq.4)) siinc= 0.0
      ix = nw
      iy = nh
      zinc = siinc
      nline = 1
      draw = 'DRAW'
      if (ivesel.eq.0) then
         nn = nn + 1
         thcrv(nn) = .02
         nxy(nn) = limitr
         do ii = 1,nxy(nn)
            xx(ii,nn) = xlim(ii)
            yy(ii,nn) = ylim(ii)
         enddo
!        call reset('thkcrv')
      else
         call pltcol(nvesel,rvs,zvs,wvs,hvs,avs,avs2,n00, &
              nn, xx, yy, nxy, msg, note, inum, xpos, ypos, ht, &
              nshd, sxx, syy, nsxy, sangle, sgap, ngaps)

      endif
      if (ifcoil.gt.0) then
         call pltcol(nfcoil,rf,zf,wf,hf,af,af2,n11, &
            nn, xx, yy, nxy, msg, note, inum, xpos, ypos, ht, &
            nshd, sxx, syy, nsxy, sangle, sgap, ngaps)
      endif
      if (iecoil.gt.0) then
         do i=1,necoil
            ae(i)=0.0
            ae2(i)=0.0
         enddo
         call pltcol(necoil,re,ze,we,he,ae,ae2,n00, &
                  nn, xx, yy, nxy, msg, note, inum, xpos, ypos, ht, &
                  nshd, sxx, syy, nsxy, sangle, sgap, ngaps)
      endif
      xabs=-3.
      yabs=7.85_dp
      dyabs = 0.16_dp
      write(text,8950) trim(ch1),trim(ch2),efitvers
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 27
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.10_dp
      if ((kdata.eq.1).or.(kdata.eq.2)) then
         write (text,8948) ifname(jtime)
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = text
         imes(msg) = 25
         xpos(msg) = xabs
         ypos(msg) = yabs
         yabs = yabs - dyabs
         ht(msg) = 0.10_dp
      endif
      xabs=.5_dp
      yabs=7.8_dp
      write (text,8960)uday
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.10_dp
      xabs=4.05_dp
      yabs=7.85_dp
      write (text,9000)ishot
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.10_dp
      write (text,9020) itime,itimeu
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.10_dp
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = 'psi contours (mv-s/r)$'
      imes(msg) = 100
      xpos(msg) = 0.0
      ypos(msg) = -.35_dp
      ht(msg) = 0.16_dp
      if (iconsi.lt.5) then
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = 'poloidal field (g)$'
         imes(msg) = 100
         xpos(msg) = 3.8_dp
         ypos(msg) = -.35_dp
         ht(msg) = 0.16_dp
      else
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = 'poloidal field (g)$'
         imes(msg) = 100
         xpos(msg) = -3.4_dp
         ypos(msg) = -.35_dp
         ht(msg) = 0.16_dp
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = 'total B field (g)$'
         imes(msg) = 100
         xpos(msg) = 3.8_dp
         ypos(msg) = -.35_dp
         ht(msg) = 0.16_dp
      endif
!
      pl_files_28: if (itek.ge.5.and.idotek.eq.0) then
      ncurve = nn
      ibrdr = 1
      xphy = 4.1_dp
      yphy = 0.5_dp
      hight = 0.10_dp
      nplen = 100
      nxlen = 0
      nylen = 0
      xtitle = 'r(m)$'
      ytitle = 'z(m)$'
      igridx = -1
      grce  = -1.0
      iexit = 1
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy, &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, xll, yll, &
      almin, almax-almin, almax, blmin, blmax-blmin, blmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      igridx, igridy, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx, yy, nxy, ncnct, &
      icont, nword, bfield, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_28
      do i=1,nw
         worka(i)=bfield(i,nh/2+1)
      enddo
      cslmin=1.0e+30_dp
      cslmax=1.0e-30_dp
      do i=1,nw
         cslmin=min(worka(i),cslmin)
         cslmax=max(worka(i),cslmax)
      enddo
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      iorel = 1
      xorl = 0.0
      yorl = yll + 0.5_dp
      nn = nn + 1
      nxy(nn) = nw
      do ii = 1, nxy(nn)
         xx(ii,nn) = rgrid(ii)
         yy(ii,nn) = worka(ii)
      enddo
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = 'R(m)$'
      imes(msg) = 100
      xpos(msg) = 1.0
      ypos(msg) = -.25_dp
      ht(msg) = 0.10_dp
      if (iconsi.lt.5) then
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = 'R(m))$'
         imes(msg) = 100
         xpos(msg) = 4.8_dp
         ypos(msg) = -.25_dp
         ht(msg) = 0.10_dp
      else
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = 'R(m)$'
         imes(msg) = 100
         xpos(msg) = -2.4_dp
         ypos(msg) = -.25_dp
         ht(msg) = 0.10_dp
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = 'R(m)$'
         imes(msg) = 100
         xpos(msg) = 4.8_dp
         ypos(msg) = -.25_dp
         ht(msg) = 0.10_dp
      endif
!
      pl_files_29: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      ncurve = nn
      ibrdr = 1
      xphy = 4.1_dp
      yphy = 0.5_dp
      hight = 0.10_dp
      nplen = -100
      nxlen = 100
      nylen = 100
      ylen = xll*.67_dp
      ixtck = 2
      iytck = 2
      igridx = 2
      igridy = 2
      idot = 1
      xtitle = '$'
      ytitle = 'psi(mv-s/r)$'
      iexit = 1
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy, &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, xll, ylen, &
      almin, almax-almin, almax, cslmin, cslmax-cslmin, cslmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      igridx, igridy, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx, yy, nxy, ncnct, &
      icont, nword, bfield, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_29
      do i=1,nw
         worka(i)=0
      enddo
!-----------------------------------------------------------------------
!--   calculate and plot contours of constant poloidal field          --
!-----------------------------------------------------------------------
      delvn=1.0_dp/(nw-1)
      do i=1,nw
         voln(i)=(i-1)*delvn
      enddo
      call zpline(nw,voln,fpol,bvoln,cvoln,dvoln)
      do ih=1,nh
         do iw=1,nw
            rw=rgrid(iw)
            rh=zgrid(ih)
            kk=(iw-1)*nh+ih
            fnow=fbrdy*tmu
            call seva2d(bkx,lkx,bky,lky,c,rw,rh,pds,ier,n333)
            if (ier.ne.0) then
               write (nttyo,9910) ier,rw,rh
               return
            endif
            dumnow=sqrt(pds(2)**2+pds(3)**2)
            bfield(iw,ih)=(dumnow)/rgrid(iw)*1.e4_dp
            if (ih.eq.nh/2+1) worka(iw)=bfield(iw,ih)
            if (xpsi(kk).le.1.0.and.ivacum.eq.0) &
               fnow=seval(nw,xpsi(kk),voln,fpol,bvoln,cvoln,dvoln)
            btttt=fnow/rgrid(iw)*10000.
            if (iconsi.ge.5) then
               dumnow=sqrt(bfield(iw,ih)**2+btttt**2)
               bfield(iw,ih)=dumnow
               if (ih.eq.nh/2+1) worka(iw)=abs(btttt)
            endif
         enddo
      enddo
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      nn = nn + 1
      icont(nn) =1
      nxy(nn) = 0
      nword = 5000
      line = 0
      mode = 'SOLID'
      lbflg = 'NOLABELS'
      xdv=6.
      cslmin=1.0e+30_dp
      cslmax=1.0e-30_dp
      do i=1,nw
         do j=1,nh
            cslmin=min(bfield(i,j),cslmin)
            cslmax=max(bfield(i,j),cslmax)
         enddo
      enddo
      siinc=(cslmax-cslmin)/5.0
      if (iconsi.ge.3) siinc= 0.0
      if (iconsi.ge.6) siinc=(cslmax-cslmin)/iconsi
      nline = 1
      draw = 'DRAW'
      nn = nn + 1
      icont(nn) = 1
      nxy(nn) = 1
      markme(nn) = 8
      ncnct(nn) = -1
      sclpc(nn) = 1.75_dp
      xx(1,nn) = rmaxis
      yy(1,nn) = zmaxis
      if (ivesel.eq.0) then
         nn = nn + 1
         thcrv(nn) = 0.02_dp
         nxy(nn) = limitr
         do ii = 1, nxy(nn)
            xx(ii,nn) = xlim(ii)
            yy(ii,nn) = ylim(ii)
         enddo
      else
         call pltcol(nvesel,rvs,zvs,wvs,hvs,avs,avs2,n00, &
              nn, xx, yy, nxy, msg, note, inum, xpos, ypos, ht, &
              nshd, sxx, syy, nsxy, sangle, sgap, ngaps)
      endif
      if (ifcoil.gt.0) &
         call pltcol(nfcoil,rf,zf,wf,hf,af,af2,n11, &
            nn, xx, yy, nxy, msg, note, inum, xpos, ypos, ht, &
            nshd, sxx, syy, nsxy, sangle, sgap, ngaps)
      if (iecoil.gt.0) then
         do i=1,necoil
            ae(i)=0.0
            ae2(i)=0.0
         enddo
         call pltcol(necoil,re,ze,we,he,ae,ae2,n00, &
            nn, xx, yy, nxy, msg, note, inum, xpos, ypos, ht, &
            nshd, sxx, syy, nsxy, sangle, sgap, ngaps)
      endif
!
      pl_files_30: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      ncurve = nn
      ibrdr = 1
      xphy = 7.625_dp
      yphy = 0.5_dp
      hight = 0.10_dp
      nplen = -100
      nxlen = 0
      nylen = 0
      ptitle = '$'
      xtitle = 'r(m)$'
      ytitle = 'z(m)$'
      igridx = -1
      grce = -1.0
      iexit = 1
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy,  &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, x11x, y11y, &
      almin, xstp, almax, blmin, ystp, blmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      igridx, igridy, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx, yy, nxy, ncnct, &
      icont, nword, zz, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_30
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      iorel = 1
      xorl = 0.0
      yorl = yll + 0.5_dp
      nn = nn + 1
      nxy(nn) = nw
      do ii = 1, nxy(nn)
         xx(ii,nn) = rgrid(ii)
         yy(ii,nn) = worka(ii)
      enddo
      if (iconsi.ge.5) then
         do i=1,nw
            dumnow=1.e4_dp*bcentr(jtime)*rcentr/rgrid(i)
            worka(i)=abs(dumnow)
         enddo
!
         cslmin=1.0e+30_dp
         cslmax=1.0e-30_dp
         do i=1,nw
            cslmin=min(worka(i) ,cslmin)
            cslmax=max(worka(i) ,cslmax)
         enddo
!
         nn = nn + 1
         ndotme(nn) = 1
         nxy(nn) = nw
         do ii = 1, nxy(nn)
            xx(ii,nn) = rgrid(ii)
            yy(ii,nn) = worka(ii)
         enddo
         do i=1,nw
            cslmin=min(worka(i) ,cslmin)
            cslmax=max(worka(i) ,cslmax)
         enddo
         delcsl=0.10_dp*(cslmax-cslmin)
         cslmax=cslmax+delcsl
         cslmin=cslmin-delcsl
      endif
!
      pl_files_31: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      ncurve = nn
      ibrdr = 1
      xphy = 7.625_dp
      yphy = 0.5_dp
      hight = 0.10_dp
      ylen = xll*.67_dp
      nplen = -100
      nxlen = 100
      nylen = 100
      xtitle = '$'
      ixtck = 2
      iytck = 2
      igridx = 2
      igridy = 2
      idot = 1
      if (iconsi.ge.5) then
        ytitle = 'Bt(gauss)$'
      else
        ytitle = 'Bp(gauss)$'
      endif
      iexit = 1
      if (iconsi.lt.5) iexit=2
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy,  &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, x11x, y11y, &
      almin, xstp, almax, blmin, ystp, blmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      igridx, igridy, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx, yy, nxy, ncnct, &
      icont, nword, zz, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_31
      Bp_contours: if (iconsi.ge.5) then
         do ih=1,nh
            do iw=1,nw
               rw=rgrid(iw)
               rh=zgrid(ih)
               kk=(iw-1)*nh+ih
               call seva2d(bkx,lkx,bky,lky,c,rw,rh,pds,ier,n333)
               if (ier.ne.0) then
                  write (nttyo,9910) ier,rw,rh
                  return
               endif
               dumnow=sqrt(pds(2)**2+pds(3)**2)
               bfield(iw,ih)=(dumnow)/rgrid(iw)*1.e4_dp
            enddo
         enddo
!-----------------------------------------------------------------------
!        Initialize plot parameters
!-----------------------------------------------------------------------
         call init2d
         nn = nn + 1
         icont(nn) =1
         nxy(nn) = 0
         nword = 5000
         line = 0
         mode = 'SOLID'
         lbflg = 'NOLABELS'
         ithk = 1
         ipri = 1
         cslmin=1.0e+30_dp
         cslmax=1.0e-30_dp
         do i=1,nw
            do j=1,nh
               cslmin=min(bfield(i,j),cslmin)
               cslmax=max(bfield(i,j),cslmax)
            enddo
         enddo
         siinc=(cslmax-cslmin)/iconsi
         nline = 1
         draw = 'DRAW'
         nn = nn + 1
         icont(nn)=1
         markme(nn) = 8
         sclpc(nn) = 1.75_dp
         ncnct(nn) = -1
         nxy(nn) = 1
         xx(1,nn) = rmaxis
         yy(1,nn) = zmaxis
         if (ivesel.eq.0) then
            nn = nn + 1
            thcrv(nn) = 0.02_dp
            nxy(nn) = limitr
            do ii = 1, nxy(nn)
               xx(ii,nn) = xlim(ii)
               yy(ii,nn) = ylim(ii)
            enddo
         else
            call pltcol(nvesel,rvs,zvs,wvs,hvs,avs,avs2,n00, &
            nn, xx, yy, nxy, msg, note, inum, xpos, ypos, ht, &
            nshd, sxx, syy, nsxy, sangle, sgap, ngaps)
         endif
         if (ifcoil.gt.0) &
            call pltcol(nfcoil,rf,zf,wf,hf,af,af2,n11, &
            nn, xx, yy, nxy, msg, note, inum, xpos, ypos, ht, &
            nshd, sxx, syy, nsxy, sangle, sgap, ngaps)
         if (iecoil.gt.0) then
            do i=1,necoil
               ae(i)=0.0
               ae2(i)=0.0
            enddo
            call pltcol(necoil,re,ze,we,he,ae,ae2,n00, &
            nn, xx, yy, nxy, msg, note, inum, xpos, ypos, ht, &
            nshd, sxx, syy, nsxy, sangle, sgap, ngaps)
         endif
!
         pl_files_32: if (itek.ge.5.and.idotek.eq.0) then
         ncurve = nn
         ibrdr = 1
         xphy = 0.675_dp
         yphy = 0.5_dp
         hight = 0.10_dp
         nplen = 100
         nxlen = 0
         nylen = 0
         xtitle = 'r(m)$'
         ytitle = 'z(m)$'
         igridx = -1
         grce = -1.0
         iexit = 1
!-----------------------------------------------------------------------
!        Write Plot Parameters
!-----------------------------------------------------------------------
         call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy,  &
         iorel, xorl, yorl, hight, bngle, bshft, &
         ptitle, nplen, xtitle, nxlen, ytitle, nylen, x11x, y11y, &
         almin, xstp, almax, blmin, ystp, blmax, &
         iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
         isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
         igridx, igridy, idash, idot, ichdsh, ichdot, &
         thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
         markme, clearx, mrc, tlen, nmrk, rat, &
         xx, yy, nxy, ncnct, &
         icont, nword, zz, ix, iy, zinc, line, mode, &
         lbflg, ithk, ipri, nline, draw, &
         nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
         nvec, xfm, yfm, xto, yto, ivec, &
         msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
         endif pl_files_32
!-----------------------------------------------------------------------
!        Initialize plot parameters
!-----------------------------------------------------------------------
         call init2d
         iorel = 1
         xorl = 0.0
         yorl = yll + 0.5_dp
         do i=1,nw
            worka(i)=bfield(i,nh/2+1)
         enddo
         nn = nn + 1
         nxy(nn) = nw
         do ii = 1, nxy(nn)
            xx(ii,nn) = rgrid(ii)
            yy(ii,nn) = worka(ii)
         enddo
         cslmin=1.0e+30_dp
         cslmax=1.0e-30_dp
         do i=1,nw
            cslmin=min(worka(i),cslmin)
            cslmax=max(worka(i),cslmax)
         enddo
         z000=0.0
         cslmin=min(cslmin,z000)
         cslmax=max(cslmax,z000)
         delcsl=0.10_dp*(cslmax-cslmin)
         cslmax=cslmax+delcsl
         cslmin=cslmin-delcsl
!
         nn=nn+1
         ndshme(nn) = 1
         nxy(nn)=2
         xx(1,nn)=rgrid(1)+0.0001_dp
         xx(2,nn)=rgrid(nw)-0.0001_dp
         yy(1,nn)=0.0
         yy(2,nn)=0.0
!
         pl_files_33: if (itek.ge.5.and.idotek.eq.0) then
         ncurve = nn
         ibrdr = 1
         xphy = 0.675_dp
         yphy = 0.5_dp
         hight = 0.10_dp
         ylen = xll*.67_dp
         nplen = -100
         nxlen = 100
         nylen = 100
         ixtck = 2
         iytck = 2
         igridx = 2
         igridy = 2
         idot = 1
         xtitle = '$'
         ytitle = 'mod-Bp(gauss)$'
         iexit = 2
!-----------------------------------------------------------------------
!        Write Plot Parameters
!-----------------------------------------------------------------------
         call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy,  &
         iorel, xorl, yorl, hight, bngle, bshft, &
         ptitle, nplen, xtitle, nxlen, ytitle, nylen, x11x, y11y, &
         almin, xstp, almax, blmin, ystp, blmax, &
         iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
         isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
         igridx, igridy, idash, idot, ichdsh, ichdot, &
         thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
         markme, clearx, mrc, tlen, nmrk, rat, &
         xx, yy, nxy, ncnct, &
         icont, nword, zz, ix, iy, zinc, line, mode, &
         lbflg, ithk, ipri, nline, draw, &
         nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
         nvec, xfm, yfm, xto, yto, ivec, &
         msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
         endif pl_files_33
      endif Bp_contours
      endif contour_plots
!-----------------------------------------------------------------------
!     Initialize plot parameters for basis function page
!-----------------------------------------------------------------------
      call ppstore
      call ffstore
      call wwstore
      call eestore
      is_vacuum_6: if (ivacum.eq.0) then
      call init2d
      xmm=1.5_dp
      curmin=1.0e+10_dp
      curmax=-1.0e+10_dp
      do i=1,nw
         curmin=min(curmin,pprime(i))
         curmax=max(curmax,pprime(i))
         xiter(i)=real(i-1,dp)/(nw-1)
      enddo    
      xdel=0.5_dp
      dcurn=curmax-curmin
      if (abs(dcurn).le.1.e-4_dp) dcurn=5.*curmax
      curmax=curmax+0.05_dp*dcurn
      curmin=curmin-0.05_dp*dcurn
      dcurn=(curmax-curmin)
      yorg = curmin
      ystp = dcurn
      ynmax = curmax
      nxy(1) = nw
      do ii = 1, nw
         xx(ii,1) = xiter(ii)
         yy(ii,1) = pprime(ii)
      enddo
      xtitle = 'XSI$'
      ytitle = 'PPRIME$'
      ncurve = 1
      ipag = 0
      iexit = 1
      xabs=-6.0
      yabs=1.5_dp
      dyabs = 0.28_dp
      dyabs = 0.22_dp
      write(text,8950) trim(ch1),trim(ch2),efitvers
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 27
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,8960) uday
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,9000) ishot
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write(text,9020) itime,itimeu
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - 2.0 * dyabs
      ht(msg) = 0.14_dp
!
      msg = msg + 1
      write(text,1001) kppcur,kffcur,kwwcur 
 1001 format(' kppcur = ',i2,' kffcur = ',i2,' kwwcur = ',i2)
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 38
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.10_dp
      msg = msg + 1
!
      pr=pres(1)/(.667_dp*wmhd(jtime)/(volume(jtime)/1.e6_dp))
      write (text,18981) pr
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 18
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.10_dp
      msg = msg + 1
!
      write(text,1002)kppfnc,kppknt,pptens 
 1002 format(' kppfnc = ',i2,' kppknt = ',i2,' pptens = ',f5.2)
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 38
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.10_dp
      write(text,"(' ppknt         ppbdry           pp2bdry ')")
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 45
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.10_dp
      write(text,"('-----------------------------------------')")
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 45
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.10_dp
      do i=1,kppknt
        write(text,"(' ',f4.2,4x,g15.5,2x,g15.5)") ppknt(i),ppbdry(i),pp2bdry(i)
        msg = msg + 1
        note(msg) = 1
        lmes(msg) = text
        imes(msg) = 45
        xpos(msg) = xabs
        ypos(msg) = yabs
        yabs = yabs - dyabs
        ht(msg) = 0.10_dp
      enddo
      yabs = yabs - dyabs
!
      msg = msg + 1
      write(text,1003)kfffnc,kffknt,fftens 
 1003 format(' kfffnc = ',i2,' kffknt = ',i2,' fftens = ',f5.2)
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 38
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.10_dp
      write(text,"(' ffknt         ffbdry           ff2bdry ')")
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 45
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.10_dp
      write(text,"('-----------------------------------------')")
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 45
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.10_dp
      do i=1,kffknt
        write(text,"(' ',f4.2,4x,g15.5,2x,g15.5)") ffknt(i),ffbdry(i),ff2bdry(i)
        msg = msg + 1
        note(msg) = 1
        lmes(msg) = text
        imes(msg) = 45
        xpos(msg) = xabs
        ypos(msg) = yabs
        yabs = yabs - dyabs
        ht(msg) = 0.10_dp
      enddo
      yabs = yabs - dyabs
!
      msg = msg + 1
      write(text,1004) kwwfnc,kwwknt,wwtens
 1004 format(' kwwfnc = ',i2,' kwwknt = ',i2,' wwtens = ',f5.2)
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 38
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.10_dp
      write(text,"(' wwknt         wwbdry           ww2bdry ')")
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 45
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.10_dp
      write(text,"('-----------------------------------------')")
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 45
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.10_dp
      if (kwwknt.gt.0) then
        do i=1,kwwknt
          write(text,"(' ',f4.2,4x,g15.5,2x,g15.5)") wwknt(i),wwbdry(i),ww2bdry(i)
          msg = msg + 1
          note(msg) = 1
          lmes(msg) = text
          imes(msg) = 45
          xpos(msg) = xabs
          ypos(msg) = yabs
          yabs = yabs - dyabs
          ht(msg) = 0.10_dp
        enddo
        yabs = yabs - dyabs
      endif
!
      if (kedgep.gt.0) then
        write(text,1005) 
 1005   format( &
          ' pe_psin       pe_width         pedge       pedgep       ')
        msg = msg + 1
        note(msg) = 1
        lmes(msg) = text
        imes(msg) = 60
        xpos(msg) = xabs
        ypos(msg) = yabs
        yabs = yabs - dyabs
        ht(msg) = 0.10_dp
        write(text,1006) 
 1006   format( &
          '----------------------------------------------------------')
        msg = msg + 1
        note(msg) = 1
        lmes(msg) = text
        imes(msg) = 60
        xpos(msg) = xabs
        ypos(msg) = yabs
        yabs = yabs - dyabs
        ht(msg) = 0.10_dp
        pedgea=pedge/darea
        pedgep=pedgea/pe_width/sidif
        write(text,"(' ',f4.2,4x,g15.5,2x,g15.5,2x,g15.5)") &
                     pe_psin,pe_width,pedgea,pedgep
        msg = msg + 1
        note(msg) = 1
        lmes(msg) = text
        imes(msg) = 60
        xpos(msg) = xabs
        ypos(msg) = yabs
        yabs = yabs - dyabs
        ht(msg) = 0.10_dp
      endif
!
      if (kedgef.gt.0) then
        write(text,"(' fe_psin       fe_width         f2edge  ')")
        msg = msg + 1
        note(msg) = 1
        lmes(msg) = text
        imes(msg) = 45
        xpos(msg) = xabs
        ypos(msg) = yabs
        yabs = yabs - dyabs
        ht(msg) = 0.10_dp
        write(text,"('-----------------------------------------')")
        msg = msg + 1
        note(msg) = 1
        lmes(msg) = text
        imes(msg) = 45
        xpos(msg) = xabs
        ypos(msg) = yabs
        yabs = yabs - dyabs
        ht(msg) = 0.10_dp
        write(text,"(' ',f4.2,4x,g15.5,2x,g15.5)") fe_psin,fe_width,f2edge
        msg = msg + 1
        note(msg) = 1
        lmes(msg) = text
        imes(msg) = 45
        xpos(msg) = xabs
        ypos(msg) = yabs
        yabs = yabs - dyabs
        ht(msg) = 0.10_dp
      endif
      write(text,9852) psiecn, dpsiecn
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 45
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.10_dp
      write(text,9853) cjeccd
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 45
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.10_dp
      write(text,9857) betped, betnped
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 45
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.10_dp
!
      pl_files_34: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      ibrdr = 1
      xphy = 6.5_dp
      yphy = 6.0
      hight = 0.14_dp
      nplen = 100
      nxlen = 100
      nylen = 100
      ixtck = 1
      iytck = 2
      intax = 1
      igridx = 1
      igridy = 2
      idot = 1
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy, &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, 2.0*xmm, xmm, &
      xiter(1), xdel, xiter(nw), yorg, ystp, ynmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      igridx, igridy, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx, yy, nxy, ncnct, &
      icont, nword, zz, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_34
      call init2d
      xmm=1.6_dp
      curmin=1.0e+10_dp
      curmax=-1.0e+10_dp
      do i=1,nw
         czmcm(i)=ffprim(i)/tmu/twopi
         curmin=min(curmin,czmcm(i))
         curmax=max(curmax,czmcm(i))
         xiter(i)=real(i-1,dp)/(nw-1)
      enddo
      xdel=0.5_dp
      dcurn=curmax-curmin
      if (abs(dcurn).le.1.e-4_dp) dcurn=5.*curmax
      curmax=curmax+0.05_dp*dcurn
      curmin=curmin-0.05_dp*dcurn
      dcurn=(curmax-curmin)
      yorg = curmin
      ystp = dcurn
      ynmax = curmax
      nxy(1) = nw
      do ii = 1, nw
         xx(ii,1) = xiter(ii)
         yy(ii,1) = czmcm(ii)
      enddo
      xtitle = 'XSI$'
      ytitle = 'FFPRIME$'
      ncurve = 1
      ipag = 0
      iexit = 1
      if (kvtor.eq.0) iexit = 2
      xabs=-5.5_dp
      yabs=1.5_dp
      dyabs = 0.28_dp
!
      pl_files_35: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      ibrdr = 1
      xphy = 6.5_dp
      yphy = 3.5_dp
      hight = 0.14_dp
      nplen = 100
      nxlen = 100
      nylen = 100
      ixtck = 1
      iytck = 2
      intax = 1
      igridx = 1
      igridy = 2
      idot = 1
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy, &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, 2.0*xmm, xmm, &
      xiter(1), xdel, xiter(nw), yorg, ystp, ynmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      igridx, igridy, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx, yy, nxy, ncnct, &
      icont, nword, zz, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_35
!
      tor_rotation: if (kvtor.gt.0) then
      call init2d
      xmm=1.5_dp
      curmin=1.0e+10_dp
      curmax=-1.0e+10_dp
      do i=1,nw
         czmcm(i)=pwprim(i)
         curmin=min(curmin,czmcm(i))
         curmax=max(curmax,czmcm(i))
         xiter(i)=real(i-1,dp)/(nw-1)
      enddo
      xdel=0.5_dp
      dcurn=curmax-curmin
      if (abs(dcurn).le.1.e-4_dp) dcurn=5.*curmax
      curmax=curmax+0.05_dp*dcurn
      curmin=curmin-0.05_dp*dcurn
      dcurn=(curmax-curmin)
      yorg = curmin
      ystp = dcurn
      ynmax = curmax
      nxy(1) = nw
      do ii = 1, nw
         xx(ii,1) = xiter(ii)
         yy(ii,1) = czmcm(ii)
      enddo
      xtitle = 'XSI$'
      ytitle = 'WPRIME$'
      ncurve = 1
      ipag = 0
      iexit = 2
      xabs=-5.5_dp
      yabs=1.5_dp
      dyabs = 0.28_dp
!
      pl_files_36: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      ibrdr = 1
      xphy = 6.5_dp
      yphy = 1.0
      hight = 0.14_dp
      nplen = 100
      nxlen = 100
      nylen = 100
      ixtck = 1
      iytck = 2
      intax = 1
      igridx = 1
      igridy = 2
      idot = 1
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy, &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, 2.0*xmm, xmm, &
      xiter(1), xdel, xiter(nw), yorg, ystp, ynmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      igridx, igridy, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx, yy, nxy, ncnct, &
      icont, nword, zz, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_36
      endif tor_rotation
      endif is_vacuum_6
#ifdef DEBUG_LEVEL1
      write (6,*) 'End  PLTOUT'
#endif
      endif not_time_snap
      ECE_or_time_snap: if ((kfitece.le.0).and.(kdata.ne.4)) then
#ifdef DEBUG_LEVEL1
      write (6,*) 'Before Exit PLTOUT'
#endif
      deallocate(bfield,sivol,voln, &
         bvoln,cvoln,dvoln,rscrap,curscr,pbimf,pmid,pmidw,bpmid, &
         cpmid,dpmid,workj,copyn,copy1,xbnow,ybnow)
#ifdef DEBUG_LEVEL1
      write (6,*) 'After Exit PLTOUT'
#endif
      return
      elseif (kdata.eq.4) then ECE_or_time_snap
!----------------------------------------------------------------------
!--   plot time history for kdata=4
!----------------------------------------------------------------------
      let = 'qt'
      call setfnmt(let,ishot,itime,dataname)
      dataname=dataname(3:7)//'_efitipmhd.dat '
      call curvec(dataname,jerror,time,ipmhd,ktime,0_i4)
      call zpline(ktime,time,sibdry,bscra,cscra,dscra)
      do i=1,ktime
         if (jerror(i).le.0) cycle
         scrat(i)=-speval(ktime,time(i),time,sibdry,bscra,cscra,dscra)
         scrat(i)=twopi*scrat(i)*1000.+vloopt(i)
      enddo
      dataname=dataname(3:7)//'_efitvsurf.dat '
      call curvec(dataname,jerror,time,scrat,ktime,0_i4)
      dataname=dataname(3:7)//'_efitkappa.dat '
      call curvec(dataname,jerror,time,elong,ktime,0_i4)
      dataname=dataname(3:7)//'_efitaminor.dat'
      call curvec(dataname,jerror,time,aminor,ktime,0_i4)
      dataname=dataname(3:7)//'_efitzsurf.dat '
      call curvec(dataname,jerror,time,zout,ktime,0_i4)
      dataname=dataname(3:7)//'_efitrsurf.dat '
      call curvec(dataname,jerror,time,rout,ktime,0_i4)
      dataname=dataname(3:7)//'_efitbetap.dat '
      call curvec(dataname,jerror,time,betap,ktime,0_i4)
      dataname=dataname(3:7)//'_efitli.dat    '
      call curvec(dataname,jerror,time,li,ktime,0_i4)
      dataname=dataname(3:7)//'_efitq95.dat   '
      call curvec(dataname,jerror,time,qpsi,ktime,0_i4)
      dataname=dataname(3:7)//'_efitbetat.dat '
      call curvec(dataname,jerror,time,betat,ktime,0_i4)
      dataname=dataname(3:7)//'_efitdensv2.dat'
      call curvec(dataname,jerror,time,dco2v(1,2),ktime,0_i4)
      dataname=dataname(3:7)//'_efitwmhd.dat  '
      call curvec(dataname,jerror,time,wmhd,ktime,0_i4)
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      call init2d
      yorg = curmin
      ystp = dcurn
      ynmax = curmax
      nn = nn + 1
      nxy(1) = nw
      do ii = 1, nw
         xx(ii,1) = xiter(ii)
         yy(ii,1) = pprime(ii)
      enddo
      curmin=1.0e+10_dp
      curmax=-1.0e+10_dp
      do i=1,nw
         czmcm(i)=ffprim(i)/tmu/twopi
         curmin=min(curmin,czmcm(i))
         curmax=max(curmax,czmcm(i))
      enddo
      dcurn=curmax-curmin
      if (abs(dcurn).le.1.e-4_dp) dcurn=5.*curmax
      curmax=curmax+0.05_dp*dcurn
      curmin=curmin-0.05_dp*dcurn
      dcurn=(curmax-curmin)
      nslen = -100
      xps = xmm
      yps = 0.0
      nn = nn + 1
      nxy(2) = nw
      ndshme(2) = 1
      do ii = 1, nxy(2)
         xx(ii,2) = xiter(ii)
         yy(ii,2) = czmcm(ii)
      enddo
      xtitle = 'XSI$'
      ytitle = 'PPRIME$'
      ncurve = nn
      ipag = 1
      iexit = 2
      xabs= -0.5_dp
      xabs= -0.6_dp
      yabs= -0.8_dp
!
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 17
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.13_dp*0.7_dp
      write (text,18983) condno
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 17
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.13_dp*0.7_dp
!
      pr=pres(1)/(.667_dp*wmhd(jtime)/(volume(jtime)/1.e6_dp))
      write (text,18981) pr
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 18
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.13_dp*0.7_dp
!
      if (abs(vcurfb(1)).gt.1.e-06_dp) then
         write (text,18985)  (vcurfb(i),i=1,3)
         msg = msg + 1
         note(msg) = 1
         lmes(msg) = text
         imes(msg) = 28
         xpos(msg) = xabs
         ypos(msg) = yabs
         yabs = yabs - dyabs
         ht(msg) = 0.13_dp*0.7_dp
      endif
      poly_rotation: if (kvtor.ge.1.and.kvtor.le.3) then
         if (icurrt.eq.5.and. &
                      abs(brsp(nfsum+1)).gt.1.e-10_dp) then
            do i=nfnpcr+1,nfnpcr+kwwcur,2
               xxnorm=brsp(i)/brsp(nfsum+1)
               if (i+1.gt.nfnpcr+kwwcur) then
                  if (i.eq.nfnpcr+1) then
                     write (text,18971) xxnorm
                  else
                     write (text,18973) xxnorm
                  endif
               else
                  xynorm=brsp(i+1)/brsp(nfsum+1)
                  if (i.eq.nfnpcr+1) then
                     write (text,28971) xxnorm,xynorm
                  else
                     write (text,28973) xxnorm,xynorm
                  endif
               endif
               msg = msg + 1
               note(msg) = 1
               lmes(msg) = text
               imes(msg) = 31
               xpos(msg) = xabs
               ypos(msg) = yabs
               yabs = yabs - dyabs
               ht(msg) = 0.13_dp*0.7_dp
            enddo
         endif
      endif poly_rotation
!
      xabs= 1.2_dp
      yabs= -0.8_dp
      if ((icurrt.eq.2.or.icurrt.eq.5).and. &
                   abs(brsp(nfsum+1)).gt.1.e-10_dp) then
         do i=nfsum+1+kppcur,nfsum+kppcur+kffcur
            xxnorm=brsp(i)/brsp(nfsum+1)
            write (text,18970) xxnorm
            msg = msg + 1
            note(msg) = 1
            lmes(msg) = text
            imes(msg) = 14
            xpos(msg) = xabs
            ypos(msg) = yabs
            yabs = yabs - dyabs
            ht(msg) = 0.13_dp*0.8_dp
         enddo
      endif
!
      pl_files_37: if (itek.ge.5.and.idotek.eq.0) then
!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      ibrdr = 1
      xphy = 4.2_dp
      yphy = 2.80_dp
      hight = 0.12_dp
      nplen = 100
      nxlen = 100
      nylen = 100
      ixtck = 1
      iytck = 2
      intax = 1
      igridx = 1
      igridy = 2
      idot = 1
!-----------------------------------------------------------------------
!     Write Plot Parameters
!-----------------------------------------------------------------------
      call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy, &
      iorel, xorl, yorl, hight, bngle, bshft, &
      ptitle, nplen, xtitle, nxlen, ytitle, nylen, xmm, xmm, &
      xiter(1), xdel, xiter(nw), yorg, ystp, ynmax, &
      iaxis, ixtck, iytck, ixnon, iynon, intax, intay, &
      isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps, &
      igridx, igridy, idash, idot, ichdsh, ichdot, &
      thcrv, sclpc, ndshme, ndotme, ncdhme, ncdtme, &
      markme, clearx, mrc, tlen, nmrk, rat, &
      xx, yy, nxy, ncnct, &
      icont, nword, zz, ix, iy, zinc, line, mode, &
      lbflg, ithk, ipri, nline, draw, &
      nshd, sxx, syy, nsxy, sangle, sgap, ngaps, &
      nvec, xfm, yfm, xto, yto, ivec, &
      msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_37
      endif ECE_or_time_snap
!---------------------------------------------------------------------
!--   plot for ece                                                  --
!--   Assign value to work arrays. workg: Bz. workk: R. workf: Psi  --
!---------------------------------------------------------------------
      if (kfitece.eq.0) return
      dxnow=(rmaxzm-rminzm)/(nw-1)
      drnow=(xlmax-xlmin)/(nw-1)
      do i=1,nw
         xnow=xlmin+(i-1)*drnow
         call seva2d(bkx,lkx,bky,lky,c,xnow,zeceo,pds,ier,n333)
         workf(i)=pds(1)+psiref(jtime)
         workg(i)=pds(2)/xnow
         workk(i)=xnow
      enddo
!---------------------------------------------------------------------
!--   writing to xece and yece, for Bz plot                         --
!---------------------------------------------------------------------
      nplece(1)=nw
      do i=1,nw
        xece(i,1)=workk(i)
        yece(i,1)=workg(i)
      enddo
      nplece(2)=1
      xnow=receo
      xece(1,2)=xnow
      call seva2d(bkx,lkx,bky,lky,c,xnow,zeceo,pds,ier,n333)
      yece(1,2)=pds(2)/xnow
      nplece(3)=2
      xece(1,3)=xlmin
      yece(1,3)=0.0
      xece(2,3)=xlmax
      yece(2,3)=0.0
!---------------------------------------------------------------------
!--   writing to pltout.out, for Bz plot                            --
!---------------------------------------------------------------------
      curmin=workg(1)
      curmax=workg(1)
      do i=1,nw
         curmin=min(curmin,workg(i))
         curmax=max(curmax,workg(i))
      enddo
      dcurn=curmax-curmin
      curmax=curmax+0.05_dp*dcurn
      curmin=curmin-0.05_dp*dcurn
      drrrr=(workk(nw)-workk(1))
      dcurn=(curmax-curmin)
      xmm=1.6_dp
!
      pl_files_38: if (itek.ge.5.and.idotek.eq.0) then
!------------------------------------------------------------------------
!--      Initialize plot parameters                                    --
!------------------------------------------------------------------------
         call init2d
         ibrdr = 1
         xphy = 9.0
         yphy = 1.4_dp
         hight = 0.12_dp
         nplen = 100
         nxlen = 100
         nylen = 100
         ixtck = 2
         iytck = 2
         igridx = 2
         igridy = 2
         idot = 1
         ipag = 0
         xtitle = 'R(m)$'
         ytitle = 'Bz(T)$'
         iexit = 1
         ncurve=3
         thcece(1)=0.0
         sclece(1)=0.0
         dshece(1)=0
         dotece(1)=0
         cdhece(1)=0
         cdtece(1)=0
         mrkece(1)=0
         clrece(1)='FOREGROUND'
         cntece(1)=0
         thcece(2)=0.0
         sclece(2)=0.70_dp
         dshece(2)=0
         dotece(2)=0
         cdhece(2)=0
         cdtece(2)=0
         mrkece(2)=0
         clrece(2)='PINK'
         cntece(2)=-1
         thcece(3)=0.0
         sclece(3)=0.0
         dshece(3)=0
         dotece(3)=1
         cdhece(3)=0
         cdtece(3)=0
         mrkece(3)=0
         clrece(3)='GREE'
         cntece(3)=0
         call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy, &
         iorel, xorl, yorl, hight, bngle, bshft,             &
         ptitle, nplen, xtitle, nxlen, ytitle, nylen, xmm, xmm,&
         workk(1), drrrr, workk(nw), curmin, dcurn, curmax,&
         iaxis, ixtck, iytck, ixnon, iynon, intax, intay,&
         isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps,&
         igridx, igridy, idash, idot, ichdsh, ichdot,&
         thcece, sclece, dshece, dotece, cdhece, cdtece,&
         mrkece, clrece, mrc, tlen, nmrk, rat,&
         xece, yece, nplece,      cntece,&
         icont, nword, zz, ix, iy, zinc, line, mode,&
         lbflg, ithk, ipri, nline, draw,&
         nshd, sxx, syy, nsxy, sangle, sgap, ngaps,&
         nvec, xfm, yfm, xto, yto, ivec,&
         msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_38
!------------------------------------------------------------------------------
!--   writing to xece and yece, for psi plot                                 --
!------------------------------------------------------------------------------
      nplece(1)=nw
      do i=1,nw
        xece(i,1)=workk(i)
        yece(i,1)=workf(i)
      enddo
      nplece(2)=2*nece+1
      do i=1,nece
         xnow=recem(i)
         call seva2d(bkx,lkx,bky,lky,c,xnow,zeceo,pds,ier,n333)
         psecem(i)=pds(1)+psiref(jtime)
         xnow=recep(i)
         call seva2d(bkx,lkx,bky,lky,c,xnow,zeceo,pds,ier,n333)
         psecep(i)=pds(1)+psiref(jtime)
      enddo
      xnow=receo
      call seva2d(bkx,lkx,bky,lky,c,xnow,zeceo,pds,ier,n333)
      pseceo=pds(1)+psiref(jtime)
      do i=1,nece
         xece(nece+1-i,2)=recem(i)
         yece(nece+1-i,2)=psecem(i)
         xece(nece+1+i,2)=recep(i)
         yece(nece+1+i,2)=psecep(i)
      enddo
      xece(nece+1,2)=receo
      yece(nece+1,2)=pseceo
!------------------------------------------------------------------------------
!--   write to coloumn 3 to 2+nece                                           --
!------------------------------------------------------------------------------
      do i=1,nece
        nplece(2+i)=2
        xece(1,2+i)=recem(i)
        yece(1,2+i)=psecem(i)
        xece(2,2+i)=recep(i)
        yece(2,2+i)=psecem(i)
      enddo
!
      curmin=workf(1)
      curmax=workf(1)
      do i=1,nw
         curmin=min(curmin,workf(i))
         curmax=max(curmax,workf(i))
      enddo
      drrrr=(workk(nw)-workk(1))
      dcurn=curmax-curmin
      curmax=curmax+0.05_dp*dcurn
      curmin=curmin-0.05_dp*dcurn
      dcurn=(curmax-curmin)
      xmm=1.6_dp
!------------------------------------------------------------------------------
!--   Initialize plot parameters                                             --
!------------------------------------------------------------------------------
      call init2d
      ibrdr = 1
      xphy = 9.0
      yphy = 3.7_dp
      hight = 0.12_dp
      nplen = 100
      nxlen = 100
      nylen = 100
      ixtck = 2
      iytck = 2
      igridx = 2
      igridy = 2
      idot = 1
      ipag = 0
!      if (iexpand .eq. 0) then
       xtitle = 'R(m)$'
       ytitle = 'PSI(VS/r)$'
!      else
!         xtitle = '$'
!         ytitle = 'ABS PSI(V-S/rad)$'
!         ixnon = 1
!      endif
      iexit = 1
!
      ncurve=2+nece
      thcece(1)=0.0
      sclece(1)=0.0
      dshece(1)=0
      dotece(1)=0
      cdhece(1)=0
      cdtece(1)=0
      mrkece(1)=0
      clrece(1)='FOREGROUND'
      cntece(1)=0
      thcece(2)=0.0
      sclece(2)=0.70_dp
      dshece(2)=0
      dotece(2)=0
      cdhece(2)=0
      cdtece(2)=0
      mrkece(2)=0
      clrece(2)='PINK'
      cntece(2)=-1
      do i=3,ncurve
        thcece(i)=0.0
        sclece(i)=0.0
        dshece(i)=0
        dotece(i)=1
        cdhece(i)=0
        cdtece(i)=0
        mrkece(i)=1
        clrece(i)='GREE'
        cntece(i)=0
      enddo
      pl_files_39: if (itek.ge.5.and.idotek.eq.0) then
        call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy,&
        iorel, xorl, yorl, hight, bngle, bshft,&
        ptitle, nplen, xtitle, nxlen, ytitle, nylen, xmm, xmm,&
        workk(1), drrrr, workk(nw), curmin, dcurn, curmax,&
        iaxis, ixtck, iytck, ixnon, iynon, intax, intay,&
        isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps,&
        igridx, igridy, idash, idot, ichdsh, ichdot,&
        thcece, sclece, dshece, dotece, cdhece, cdtece,&
        mrkece, clrece, mrc, tlen, nmrk, rat,&
        xece, yece, nplece,       cntece,&
        icont, nword, zz, ix, iy, zinc, line, mode,&
        lbflg, ithk, ipri, nline, draw,&
        nshd, sxx, syy, nsxy, sangle, sgap, ngaps,&
        nvec, xfm, yfm, xto, yto, ivec,&
        msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif pl_files_39
!------------------------------------------------------------------------------
!--   psip-psim plot                                                         --
!------------------------------------------------------------------------------
      nplece(1)=nece
      ddpsi=abs(psim(jtime)-sibdry(jtime))
      curmin=-10.
      curmax=10.
      do i=1,nece
        xece(i,1)=i
        yece(i,1)=(psecep(i)-psecem(i))/ddpsi*100.
        curmin=min(curmin,yece(i,1))
        curmax=max(curmax,yece(i,1))
      enddo
      curmin=int(curmin) -1.0
      curmax=int(curmax) +1.0
      dcurn = (curmax - curmin) / 2.
      nplece(2)=2
      xece(1,2)=0
      yece(1,2)=0.0
      xece(2,2)=nece+2
      yece(2,2)=0.0
!
      xorg=0
      drrrr=2.
      xmax=nece+2
      xmm=1.6_dp
      call init2d
      intax = 1
      intay = 1
      ibrdr = 1
      xphy = 9.0
      yphy = 6.0
      hight = 0.12_dp
      nplen = 100
      nxlen = 100
      nylen = 100
      ixtck =0
      iytck =0
      igridx = 0
      igridy = 2
      idot = 1
      ipag = 0
      xtitle = 'ECE$'
      ytitle = 'dPSI(%)$'
      iexit = 1
      ncurve=2
      thcece(1)=0.0
      sclece(1)=0.7_dp
      dshece(1)=0
      dotece(1)=0
      cdhece(1)=0
      cdtece(1)=0
      mrkece(1)=16
      clrece(1)='PINK'
      cntece(1)=-1
      thcece(2)=0.0
      sclece(2)=0.0
      dshece(2)=0
      dotece(2)=1
      cdhece(2)=0
      cdtece(2)=0
      mrkece(2)=0
      clrece(2)='FOREGROUND'
      cntece(2)=0
      if (itek.ge.5.and.idotek.eq.0) then
        call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy,&
        iorel, xorl, yorl, hight, bngle, bshft,&
        ptitle, nplen, xtitle, nxlen, ytitle, nylen, xmm, xmm,&
        xorg, drrrr, xmax, curmin, dcurn, curmax,&
        iaxis, ixtck, iytck, ixnon, iynon, intax, intay,&
        isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps,&
        igridx, igridy, idash, idot, ichdsh, ichdot,&
        thcece, sclece, dshece, dotece, cdhece, cdtece,&
        mrkece, clrece, mrc, tlen, nmrk, rat,&
        xece, yece, nplece,       cntece,&
        icont, nword, zz, ix, iy, zinc, line, mode,&
        lbflg, ithk, ipri, nline, draw,&
        nshd, sxx, syy, nsxy, sangle, sgap, ngaps,&
        nvec, xfm, yfm, xto, yto, ivec,&
        msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif
!------------------------------------------------------------------------------
!--   chi2ece plot                                                           --
!------------------------------------------------------------------------------
      nplece(1)=nece
      do i=1,nece
        xece(i,1)=i
        yece(i,1)=chiece(i)
      enddo
      nplece(2)=1
      xece(1,2)=nece+1
      yece(1,2)=chiecebz
      xorg=0
      drrrr=2
      xmax=nece+2.0
      curmax=chiecebz
      do i=1,nece
         curmax=max(curmax,chiece(i))
      enddo
      curmax=int(curmax)+1
      curmin=0.0
!      dcurn=0.0
      dcurn=(curmax-curmin)/2.
      xmm=1.6_dp
      call init2d
      intax = 1
      ibrdr = 1
      xphy = 6.5_dp
      yphy = 6.0
      hight = 0.12_dp
      nplen = 100
      nxlen = 100
      nylen = 100
      ixtck =0
      iytck = 0
      igridx = 1
      igridy =1
      idot = 1
      ipag = 0
      xtitle = 'ECE$'
      ytitle = 'CHI**2$'
      iexit = 1
      if ((kfixro.gt.0).and.(kfixrece.gt.0)) iexit = 2
      ncurve=2
      thcece(1)=0.0
      sclece(1)=0.7_dp
      dshece(1)=0
      dotece(1)=0
      cdhece(1)=0
      cdtece(1)=0
      mrkece(1)=16
      clrece(1)='PINK'
      cntece(1)=1
      thcece(2)=0.0
      sclece(2)=0.7_dp
      dshece(2)=0
      dotece(2)=1
      cdhece(2)=0
      cdtece(2)=0
      mrkece(2)=15
      clrece(2)='GREEN'
      cntece(2)=-1
!
      xabs=-6.0
      yabs=1.8_dp
      dyabs = 0.22_dp
      write(text,8950) trim(ch1),trim(ch2),efitvers
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 27
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,8960) uday
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9000) ishot
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9020) itime,itimeu
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9960) receo
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9962) zeceo
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9965)
      do k=1,nece
        msg = msg + 1
        note(msg) = 1
        lmes(msg) = text
        imes(msg) = 25
        xpos(msg) = xabs
        ypos(msg) = yabs
        yabs = yabs - dyabs
        ht(msg) = 0.14_dp
        write (text,9968) recem(k),recep(k)
      enddo
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9969) xfit(1)
      do k=2,nfit
        msg = msg + 1
        note(msg) = 1
        lmes(msg) = text
        imes(msg) = 25
        xpos(msg) = xabs
        ypos(msg) = yabs
        yabs = yabs - dyabs
        ht(msg) = 0.14_dp
        write (text,9972) xfit(k)
      enddo
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9976) chisqfit
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9978) tchiece
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9979) chiecebz
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
!
      if (itek.ge.5.and.idotek.eq.0) then
        call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy,&
        iorel, xorl, yorl, hight, bngle, bshft,&
        ptitle, nplen, xtitle, nxlen, ytitle, nylen, xmm, xmm,&
        xorg, drrrr, xmax, curmin, dcurn, curmax,&
        iaxis, ixtck, iytck, ixnon, iynon, intax, intay,&
        isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps,&
        igridx, igridy, idash, idot, ichdsh, ichdot,&
        thcece, sclece, dshece, dotece, cdhece, cdtece,&
        mrkece, clrece, mrc, tlen, nmrk, rat,&
        xece, yece, nplece,       cntece,&
        icont, nword, zz, ix, iy, zinc, line, mode,&
        lbflg, ithk, ipri, nline, draw,&
        nshd, sxx, syy, nsxy, sangle, sgap, ngaps,&
        nvec, xfm, yfm, xto, yto, ivec,&
        msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif
!-----------------------------------------------------------------------
!     writing to xece and yece, for ECE data fitting
!         kfixro kfixrece >0 real space, no mapping
!         kfixro kfixrece =0, fitting Te(R)
!         kfixro kfixrece <0, fitting Te(B)
!-----------------------------------------------------------------------
      ECE: if ((kfixro.eq.0).or.(kfixrece.eq.0)) then
!
      nplece(1)=mecein
      curmax=teeceinr(1)
      do i=1,mecein
        xece(i,1)=recein(i)
        yece(i,1)=teeceinr(i)
        curmax=max(curmax,teeceinr(i))
      enddo
      nplece(2)=(nnnte-1)/4+1
      do i=1,(nnnte-1)/4+1
        ii=4*i-3
        xece(i,2)=rrr(ii)
        yece(i,2)=teecer(ii)
        curmax=max(curmax,teecer(ii))
      enddo
      curmax=int(curmax)+1.
!------------------------------------------------------------------------------
!--   writing to pltout.out. For ECE data fitting                            --
!------------------------------------------------------------------------------
      xorg=1.2_dp
      drrrr=0.4_dp
      xmax=2.4_dp
      curmax=8.0
      curmin=0.0
      dcurn=(curmax - curmin)/2.
      xmm=1.6_dp
      call init2d
      intax = 1
      ibrdr = 1
      xphy = 6.5_dp
      yphy = 3.7_dp
      hight = 0.12_dp
      nplen = 100
      nxlen = 100
      nylen = 100
      ixtck = 1
      iytck = 1
      igridx = 1
      igridy = 1
      idot = 1
      ipag = 0
      xtitle = 'R(m)$'
      ytitle = 'Te(keV)$'
      iexit = 2
      ncurve=2
      thcece(1)=0.0
      sclece(1)=0.7_dp
      dshece(1)=0
      dotece(1)=0
      cdhece(1)=0
      cdtece(1)=0
      mrkece(1)=16
      clrece(1)='PINK'
      cntece(1)=-1
      thcece(2)=0.0
      sclece(2)=0
      dshece(2)=0
      dotece(2)=0
      cdhece(2)=0
      cdtece(2)=0
      mrkece(2)=0
      clrece(2)='FOREGROUND'
      cntece(2)=0
      if (itek.ge.5.and.idotek.eq.0) then
        call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy,&
        iorel, xorl, yorl, hight, bngle, bshft,&
        ptitle, nplen, xtitle, nxlen, ytitle, nylen, xmm, xmm,&
        xorg, drrrr, xmax, curmin, dcurn, curmax,&
        iaxis, ixtck, iytck, ixnon, iynon, intax, intay,&
        isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps,&
        igridx, igridy, idash, idot, ichdsh, ichdot,&
        thcece, sclece, dshece, dotece, cdhece, cdtece,&
        mrkece, clrece, mrc, tlen, nmrk, rat,&
        xece, yece, nplece,       cntece,&
        icont, nword, zz, ix, iy, zinc, line, mode,&
        lbflg, ithk, ipri, nline, draw,&
        nshd, sxx, syy, nsxy, sangle, sgap, ngaps,&
        nvec, xfm, yfm, xto, yto, ivec,&
        msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif
!
      elseif ((kfixro.lt.0).and.(kfixrece.lt.0)) then ECE
!----------------------------------------------------------------------
!--   Te(B)
!----------------------------------------------------------------------
      nplece(1)=necein
      do i=1,necein
        xece(i,1)=becein(i)
        yece(i,1)=teecein0(i)
      enddo
      nplece(2)=(nnnte-1)/4+1
      do i=1,(nnnte-1)/4+1
        ii=4*i-3
        xece(i,2)=bbf(ii)
        yece(i,2)=teeceb(ii)
      enddo
      xorg=1.6_dp
      drrrr=0.4_dp
      xmax=3.4_dp
      curmax=8.0
      curmin=0.0
      dcurn=2.0
      xmm=1.6
      call init2d
      intax = 1
      ibrdr = 1
      xphy = 6.5_dp
      yphy = 3.7_dp
      hight = 0.12_dp
      nplen = 100
      nxlen = 100
      nylen = 100
      ixtck =0
      iytck = 0
      igridx = 0
      igridy =0
      idot = 1
      ipag = 0
      xtitle = 'B(T)$'
      ytitle = 'Te(keV)$'
      iexit = 1
      ncurve=2
      thcece(1)=0.0
      sclece(1)=0.7_dp
      dshece(1)=0
      dotece(1)=0
      cdhece(1)=0
      cdtece(1)=0
      mrkece(1)=16
      clrece(1)='PINK'
      cntece(1)=-1
      thcece(2)=0.0
      sclece(2)=0
      dshece(2)=0
      dotece(2)=0
      cdhece(2)=0
      cdtece(2)=0
      mrkece(2)=0
      clrece(2)='FOREGROUND'
      cntece(2)=0
      if (itek.ge.5.and.idotek.eq.0) then
        call curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy,&
        iorel, xorl, yorl, hight, bngle, bshft,&
        ptitle, nplen, xtitle, nxlen, ytitle, nylen, xmm, xmm,&
        xorg, drrrr, xmax, curmin, dcurn, curmax,&
        iaxis, ixtck, iytck, ixnon, iynon, intax, intay,&
        isaxs, sorg, stp, smax, slen, sname, nslen, xps, yps,&
        igridx, igridy, idash, idot, ichdsh, ichdot,&
        thcece, sclece, dshece, dotece, cdhece, cdtece,&
        mrkece, clrece, mrc, tlen, nmrk, rat,&
        xece, yece, nplece,       cntece,&
        icont, nword, zz, ix, iy, zinc, line, mode,&
        lbflg, ithk, ipri, nline, draw,&
        nshd, sxx, syy, nsxy, sangle, sgap, ngaps,&
        nvec, xfm, yfm, xto, yto, ivec,&
        msg, note, lmes, imes, anum, iplce, inum, xpos, ypos, ht, iexit)
      endif
!
      endif ECE
      xabs=-6.0
      yabs=1.8_dp
      dyabs = 0.22_dp
      write(text,8950) trim(ch1),trim(ch2),efitvers
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 27
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,8960) uday
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9000) ishot
      msg = msg + 1
      note(msg) = 1
      lmes(msg) = text
      imes(msg) = 25
      xpos(msg) = xabs
      ypos(msg) = yabs
      yabs = yabs - dyabs
      ht(msg) = 0.14_dp
      write (text,9020) itime,itimeu
!
#ifdef DEBUG_LEVEL1
      write (6,*) 'Before Exit PLTOUT'
#endif
      deallocate(bfield,sivol,voln, &
         bvoln,cvoln,dvoln,rscrap,curscr,pbimf,pmid,pmidw,bpmid, &
         cpmid,dpmid,workj,copyn,copy1,xbnow,ybnow)
#ifdef DEBUG_LEVEL1
      write (6,*) 'After Exit PLTOUT'
#endif
!
      return
 8948 format (a25)
 8950 format (1x,1('h'),'EFIT-AI ',a4,' x ',a4,' ',a7,1('h'))
 8960 format (' date ran = ',a10)
 9000 format (' shot #   = ',i10)
 9002 format (' shot #   = ',i10,' date ran = ',a10)
 9020 format (' t(ms,us) = ',i6,1x,i6)
 9022 format (' chiprw   = ',1pe10.3)
 9024 format (' Rvt(m)   = ',1pe10.3)
 9026 format (' Kvt(m)   = ',i10)
 9040 format (' chi2(mag)= ',1pe10.3)
 9060 format (' rout(cm) = ',f10.2)
 9080 format (' zout(cm) = ',f10.3)
 9100 format (' a(cm)    = ',f10.2)
 9120 format (' elong    = ',f10.3)
 9140 format (' utri,ltri= ',f6.2,1x,f6.2)
 9150 format (' betan,In = ',f6.2,1x,f6.2)
 9155 format (' indent   = ',f10.3)
 9160 format (' V,A(m3,2)= ',f6.2,1x,f6.3)
 9165 format (' energy(j)= ',1pe10.3)
 9180 format (' betat(%) = ',f10.2)
 9200 format (' betap    = ',f10.2)
 9220 format (' li,li3   = ',f6.2,1x,f6.2)
 9230 format (' error, # = ',1pe10.3,1x,i3)
 9240 format (' delstar  = ',1pe8.1,1x,1pe8.1)
 9260 format (' fpols(kA)= ',f10.1)
 9262 format (' ux,lx(cm)= ',f6.2,1x,f6.2)
 9265 format (' J0n,J1n  = ',f6.2,1x,f6.2)
 9268 format (' q1,q95   = ',f6.2,1x,f6.2)
 9272 format (' dsep(cm) = ',f10.3)
 9270 format (' ipf(ka)  = ',f10.1)
 9280 format (' rm,rc(cm)= ',f6.1,1x,f6.1)
 9290 format (' zm,zc(cm)= ',f6.2,1x,f6.2)
 9300 format (4x,'   data used:   ')
 9320 format (1x,i2,' flux loops, ',i2,' ECE')
 9340 format (1x,i2,' magnetic probes   ')
 9360 format (1x,i2,' partial rogowskis')
 9380 format (1x,i2,' rogow',1x,i2,' fc ',1x,i2,' ec ')
 9385 format (1x,i2,' di, ',i2,' mse ',i2,' li  ')
99385 format (1x,i2,' di, ',i2,' mse ',i2,' mls ')
 9390 format (' ip(ka)   = ',f10.1)
 9395 format (' scrapof(mm)',1x,5f5.1)
 9399 format (' n/nc,qmer= ',f6.3,1x,f6.3)
 9400 format (' bt0(t)   = ',f10.3)
 9408 format ('          = ',f10.3)
 9420 format (' lin,o(cm)= ',f6.2,1x,f6.2)
 9440 format (' chi2max  = ',1pe10.3)
 9441 format (' fexpx,vs = ',f6.2,1x,f6.2)
 9460 format (' lt,b(cm) = ',f6.2,1x,f6.2)
 9465 format (' s/q2,taun= ',f6.2,1x,f6.2)
 9467 format (' Zt(cm),RL= ',f6.2,1x,f6.2)
 9470 format (' sib(vs/r)= ',1pe10.3)
 9480 format (' elongm,qm= ',f6.2,1x,f6.3)
 9482 format (' sir(vs/r)= ',1pe10.3)
 9484 format (' n(Rc,Zc) = ',f10.3)
 9490 format (' shearb,Jb= ',f6.2,1x,f6.3)
 9492 format (' dqdsi/q-b= ',f10.3)
 9494 format (' aq1,2(cm)= ',f6.2,1x,f6.2)
 9496 format (' siw, q   = ',f6.3,1x,f6.2)
 9498 format (' shearw,Jw= ',f6.2,1x,f6.3)
 9495 format (' chidlc   = ',1pe10.3)
 9497 format (' chipre   = ',1pe10.3)
 9500 format (' errbp    = ',1pe10.3)
 9502 format (' nev2(cm3)= ',1pe10.3)
 9504 format (' erbmax   = ',1pe10.3)
 9506 format (' erbave   = ',1pe10.3)
 9507 format (' ner0(cm3)= ',1pe10.3)
 9510 format (' errbpli2 = ',1pe10.3)
 9512 format (' taue(ms) = ',1pe10.3)
 9520 format (' Rvsu(cm) = ',f6.1,1x,f6.1)
 9521 format (' Rvsd(cm) = ',f6.1,1x,f6.1)
 9522 format (' rtch(cm) = ',f10.3)
 9524 format (' ztch(cm) = ',f10.3)
 9530 format (' rsep(cm) = ',f6.1,1x,f6.1)
 9540 format (' zsep(cm) = ',f6.1,1x,f6.1)
 9542 format (' betapd,w = ',f6.2,1x,f6.2)
 9544 format (' betatd,w = ',f6.2,1x,f6.2)
 9546 format (' wdia(J)  = ',1pe10.3)
 9547 format (' wtor(J)  = ',1pe10.3)
 9548 format (' taued(ms)= ',1pe10.3)
 9550 format (' vl(V),Jav= ',f6.2,1x,f6.2)
 9552 format (' qmericer = ',f10.2)
 9555 format (' chitot   = ',1pe10.3)
 9560 format (' p1/p0(f) = ',1pe10.3)
 9562 format (' p1/p0(d) = ',1pe10.3)
 9565 format (' kppcur   = ',i10)
 9570 format (' kffcur   = ',i10)
 9575 format (' pcurbd   = ',f10.3)
 9580 format (' fcurbd   = ',f10.3)
 9603 format (' icp,  ivs= ',2i5)
 9605 format (' kf, kp, E= ',i4,2i3)
 9610 format (' fpd, p, E= ',3f4.1)
 9612 format (' kw ,  wbd= ',i5,f5.1)
 9614 format (' kvt      = ',i5)
 9615 format (' fwbp, fwq= ',2f5.1)
 9620 format (' kcf,  kcp= ',2i5)
 9623 format (' kcw      = ',1i5)
 9625 format (' alfp     = ',f10.3)
96251 format (' alfp     = ',f10.3)
96252 format (' Ivt(ka)  = ',f10.3)
96253 format (' Ipt(ka)  = ',f10.3)
96254 format (' Ivp(ka)  = ',f10.3)
 9630 format (' kprfit   = ',i10)
 9635 format (' Ze(cm)   = ',f10.3)
 9637 format (' relax    = ',f10.3)
 9639 format (' Rj(m)    = ',f10.3)
 9700 format ('solid lines = calculated values')
 9710 format ('symbols = experimental values')
 9800 format (' chisqte  = ',1pe10.3)
 9805 format (' chisqne  = ',1pe10.3)
 9810 format (' fco2ne   = ',f10.3)
 9815 format (' nptef    = ',i10)
 9820 format (' npnef    = ',i10)
 9825 format (' chisqti  = ',1pe10.3)
 9828 format (' zt(cm),RL= ',f6.2,1x,f6.2)
 9830 format (' nptionf  = ',i10)
 9832 format (' p0(n/m2) = ',1pe10.3)
 9834 format (' dp1dx(d) = ',1pe10.3)
 9835 format (' dp1dx(f) = ',1pe10.3)
 9836 format (' kpressb  = ',i10)
 9838 format (' J(0.95)  = ',1pe10.3)
 9840 format (' pp(.95)  = ',1pe10.3)
 9842 format (' bimf,e(%)= ',f6.2,1x,f6.2)
 9852 format (' siec,dsi = ',f8.4,1x,f8.4)
 9853 format (' eccd(kA) = ',1pe10.3)
 9857 format (' betped,bn= ',f8.4,1x,f8.4)
 9900 format ('  error in sets2d = ',i4)
 9910 format ('  error in spline =',i4,' (r,z)= ( ',f5.2,',',f5.2,')')
 9930 format (' chi2ms,li= ',f7.2,1x,f7.2)
99930 format (' chi2ms,ls= ',f7.2,1x,1pe10.3)
 9932 format (' Zeff-res = ',f10.2)
 9934 format (' Tave(KeV)= ',f10.3)
 9936 format (' Pb(MW),Wn= ',f6.2,1x,f6.2)
 9937 format (' SYMMETRIZED SOLUTION')
 9938 format (' Ifb(KA)  = ',1pe10.3)
 9940 format (' qgam     = ',8(1x,f5.2))
99940 format (' qmls     = ',8(1x,f5.2))
 9945 format (' qsiw     = ',8(1x,f5.2))
 9950 format (' Fig. ',i3)
 9960 format (' receo(m) = ',f6.3)
 9962 format (' zeceo(m) = ',f6.3)
 9965 format (' recem(m),  recep(m)')
 9968 format (2x,f6.3,5x,f6.3)
 9969 format (' xfit     = ',f6.3)
 9972 format (12x,f6.3)
 9976 format (' chisqfit = ',f6.3)
 9978 format (' tchiece  = ',1pe10.3)
 9979 format (' chiecebz = ',1pe10.3)
10000 format (6e12.6)
10020 format (5e10.4)
14980 format (i5)
15000 format (2e12.6)
18950 format (' a, g =',1pe10.3)
18960 format ('       ',1pe10.3)
18970 format (' ',1pe10.3)
18971 format (' w    = ',1pe11.4)
18973 format ('        ',1pe11.4)
18980 format (' a0  = ',1pe10.3)
18981 format (' p0/<p>=',1pe9.1)
18983 format (' cno = ',1pe10.3)
18985 format (' Vfb = ',1pe10.3,'    ',0pf3.0,' ',0pf3.0)
19603 format (' rm,al    = ',f7.3,1x,f7.3)
19605 format (' bt,bw    = ',f7.3,1x,f7.3)
19610 format (' saaa(cm) = ',f7.3)
28971 format (' w    = ',1pe11.4,1x,1pe11.4)
28973 format ('        ',1pe11.4,1x,1pe11.4)
      end subroutine pltout

!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          expand.                                                 **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**********************************************************************
      subroutine expand(n1,n2,nexexx,xmm,jtime)
      use commonblocks,only: worka,byringr,byringz,xxtra,yxtra, &
                             bpxtra,flxtra,fpxtra
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)

      real*8,dimension(:),allocatable :: xpltloc,ypltloc,xplxloc,yplxloc
!      equivalence (xpltloc,flxtra(1,1))
!      equivalence (ypltloc,fpxtra(1,1))
!      equivalence (xplxloc,flxtra(1,2))
!      equivalence (yplxloc,fpxtra(1,2))
      data almin,almax,blmin,blmax/1.51_dp,1.79_dp,-1.38_dp,-1.21_dp/
      almin=1.62_dp
      blmax=-1.20_dp
      call nobrdr
      call grace(0.0)
      call xticks(2)
      call yticks(2)
      call physor(9.0,1.4)
!
      allocate(xpltloc(npoint),ypltloc(npoint),xplxloc(npoint),yplxloc(npoint))
!
      ilow=0  
      if(zseps(1,jtime).gt.-900.) ilow=1
      if((zseps(2,jtime).gt.-900.).and. &
       (zseps(2,jtime).lt.zseps(1,jtime))) ilow=2
!      if(ilow.ne.0) then
!         blmax=max(zseps(ilow,jtime)/100.,blmax)
!         almin=min(rseps(ilow,jtime)/100.,almin)
!      endif
      call title('$',-100,'R(m)$',100,'Z(m)$',100,xmm,xmm)
      call graf(almin,'SCAL',almax,blmin,'SCAL',blmax)
      call frame
      !call grid(0,0) !TODO: undefined
      call thkcrv(0.02_dp)
      call curve(xlim,ylim,limitr,0)
      call dash
      call curve(xpltloc,ypltloc,n1,0)
      call dot
      call curve(byringr,byringz,nh2,0)
      call reset('dot')
      call marker(15)
      call sclpic(.5)
      call curve(xplxloc,yplxloc,n2,-1)
      call reset('sclpic')
      call reset('marker')
      do i=1,nexexx
       if(npxtra(i).le.2) cycle
       call curve(xxtra(1,i),yxtra(1,i),npxtra(i),0)
      enddo
      fpxtra(1:nxtrap,1:nxtram)=0
      flxtra(1:nxtrap,1:nxtram)=0
!
      deallocate(xpltloc,ypltloc,xplxloc,yplxloc)
!
      return
      end subroutine expand

!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          pltcol plots the external coils.                        **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          08/05/85..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine pltcol(ncoil,rc,zc,wc,hc,ac,ac2,inum, &
      nn, xx, yy, nxy, nmg, note, num, xpos, ypos, ht, &
      nshd, sx, sy, nsxy, sangle, sgap, ngaps)
      use set_kinds, only: dp
      use global_constants, only: pi
      use eparm, only: ndim
      use curve2d_mod, only: ncrv, mdim
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension rc(ncoil), zc(ncoil), wc(ncoil), &
      hc(ncoil), xx(ndim,ncrv), yy(ndim,ncrv), nxy(ncrv), &
      num(mdim), note(mdim), xpos(mdim), ypos(mdim), ht(mdim), &
      sx(ncrv,ncrv), sy(ncrv,ncrv), nsxy(ncrv), sangle(ncrv), &
      sgap(ncrv), ngaps(ncrv)
      dimension ac(ncoil),ac2(ncoil)
!
      do i=1,ncoil
         x= rc(i)
         y= zc(i)
         dx=wc(i)/2.
         dy=hc(i)/2.
         sn1=sin(ac(i))*180./pi
         cos1=cos(ac(i))*180./pi
         tac=sn1/cos1
         sn2=sin(ac2(i))*180./pi
         cos2=cos(ac2(i))*180./pi
         tac2=sn2/cos2
         if(ac2(i).ne.0) then
            nn = nn + 1
            xx(1,nn)=x-dx-dy/tac2
            xx(3,nn)=x+dx+dy/tac2
            xx(5,nn)=xx(1,nn)
            xx(2,nn)=xx(3,nn)-wc(i)
            xx(4,nn)=xx(1,nn)+wc(i)
!
            yy(1,nn)=y-dy
            yy(2,nn)=y+dy
            yy(3,nn)=y+dy
            yy(4,nn)=y-dy
            yy(5,nn)=y-dy
         else
!
            dac=0.5_dp*wc(i)*tac
            nn = nn + 1
            xx(1,nn)=x-dx
            xx(2,nn)=x-dx
            xx(3,nn)=x+dx
            xx(4,nn)=x+dx
            xx(5,nn)=x-dx
!
            yy(1,nn)=y-dy-dac
            yy(2,nn)=y+dy-dac
            yy(3,nn)=y+dy+dac
            yy(4,nn)=y-dy+dac
            yy(5,nn)=y-dy-dac
         endif
         nxy(nn) = 5
         if (inum.eq.-100) then
            nshd = nshd + 1
            sangle(nshd) = 90.0
            sgap(nshd) = 0.01_dp
            ngaps(nshd) = 1
            nsxy(nshd) = 5
            do ii = 1, nsxy(nshd)
               sx(ii, nshd) = xx(ii,nn)
               sy(ii, nshd) = yy(ii,nn)
            enddo
         endif
         dwh=.015
         if (inum.gt.0) then
            nmg = nmg + 1
            note(nmg) = 6
            num(nmg) = i
            xpos(nmg) = rc(i) - dwh
            ypos(nmg) = zc(i) - dwh
            ht(nmg) = 0.04_dp
         endif
      enddo
      return
      end subroutine pltcol

!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          CURVEC plots xp and yp based on array jerror.           **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          11/04/86..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine curvec(dataname,jerror,xp,yp,np,ns)
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension jerror(1),xp(1),yp(1)
      character*(*) dataname
!
      open(unit=62,file=dataname,status='old',iostat=ioerr)
      if (ioerr.eq.0) close(unit=62,status='delete')
      open(unit=62,file=dataname,status='new')
      do i=1,np
        if (jerror(i).gt.0) cycle
        write (62,6200) xp(i),yp(i)
      enddo
      close(unit=62)
 6200 format (1x,e12.5,1x,e12.5)
      return
      end subroutine curvec

!***********************************************************************
!**                                                                   **
!**   SUBPROGRAM DESCRIPTION:                                         **
!**     Determines whether plot parameters are to be written in       **
!**     ASCII or BINARY format                                        **
!**                                                                   **
!**   CALLING ARGUMENTS:                                              **
!**     ncurve  Number of curves                                      **
!**     ipag page dimension flag                                      **
!**       0 = page dimension 11 x 8.5                                 **
!**       1 = page dimension 8.5 x 11                                 **
!**     ibrdr   Page border flag                                      **
!**       1 = Suppress page border                                    **
!**       0 = Suppress page border                                    **
!**     grce >= 0 Enable grace margin                                 **
!**     xphy X-coordinate of physical origin                          **
!**     yphy Y-coordinate of physical origin                          **
!**     hight   Font size                                             **
!**     bngle   Base rotation angle                                   **
!**     bshft   Base translation                                      **
!**     ptitle  Plot title                                            **
!**     pltlen  Length of plot title                                  **
!**     xtitle  X-axis name                                           **
!**     xnlen   Length of x-axis name                                 **
!**     ytitle  Y-axis name                                           **
!**     ynlen   Length of y-axis name                                 **
!**     xlen Length of x-axis legend                                  **
!**     ylen    Length of y-axis legend                               **
!**     iaxis   Axes flag                                             **
!**       0 = Linear axis                                             **
!**       1 = Log-Linear axis                                         **
!**       2 = Linear-Log axis                                         **
!**       3 = Logarithmic axis                                        **
!**     xtck    X-axis tick marks            **
!**     ytck    Y-axis tick marks            **
!**     ixnon   X-axis tick marks or labels flag     **
!**       1 = Suppress x-axis tick marks or labels **
!**       0 = X-axis tick marks or labels     **
!**     iynon   Y-axis tick marks or labels flag     **
!**       1 = Suppress y-axis tick marks or labels **
!**       0 = Y-axis tick marks or labels     **
!**     intax   Trailing zeroes on x-axis flag     **
!**       1 = Suppress trailing zeroes on x-axis      **
!**       0 = Trailing zeroes on x-axis      **
!**     intay   Trailing zeroes on y-axis flag     **
!**       1 = Suppress trailing zeroes on y-axis      **
!**       0 = Trailing zeroes on y-axis      **
!**     xorg    X-value at the physical origin     **
!**     xstp    X step size in units            **
!**     xmax    Value at x-axis limit            **
!**     yorg    X-value at the physical origin                 **
!**     ystp    X step size in units            **
!**     ymax    Value at y-axis limit            **
!**     iorel   New physical origin relative to current  **
!**     origin flag             **
!**       1 = New physical origin relative to current  **
!**           origin             **
!**       0 = No new physical origin relative to current  **
!**           origin             **
!**     xorl    X-coordinate of relative origin     **
!**     yorl    Y-coordinate of relative origin     **
!**     igridx  Number of grid lines per step on x-axis       **
!**     igridy  Number of grid lines per step on y-axis       **
!**     idash   Dash grid lines flag            **
!**       1 = Dash grid lines            **
!**       0 = No dash grid lines            **
!**     idot    Dot grid lines flag            **
!**       1 = Dot grid lines            **
!**       0 = No dot grid lines            **
!**     ichdot  Chain dot grid lines flag                       **
!**       1 = Chain dot grid lines                        **
!**       0 = No chain dot grid lines                     **
!**     ichdsh  Chain dash grid lines flag                      **
!**       1 = Chain dash grid lines                       **
!**       0 = No chain dash grid lines                    **
!**     thcrv   Curve thickness     dim = ncurve       **
!**     sclpc   Scale curve marker  dim = ncurve            **
!**     dashme  Dash curve flag     dim = ncurve     **
!**       1 = Dash curve             **
!**       0 = No dash curve            **
!**     dotme   Dot curve flag     dim = ncurve     **
!**       1 = Dot curve             **
!**       0 = No dot curve            **
!**     chdhme  Chain dash curve flag dim = ncurve     **
!**       1 = Chain dash curve            **
!**       0 = No chain dash curve            **
!**     chdtme  Chain ot curve flag dim = ncurve     **
!**       1 = Chain dot curve            **
!**       0 = No chain dot curve            **
!**     markme  Curve marker          dim = ncurve            **
!**     clearx  Color     dim = ncurve     **
!**     x       Array of x-coordinates            **
!**                 dim = (nplt, ncurve) **
!**     y       Array of y-coordinates            **
!**                 dim = (nplt, ncurve) **
!**     nplt    Number of (x,y) points to be ploted      **
!**                 dim = ncurve     **
!**     ncnct   Marker specification            **
!**       0 = Points connected with no symbols drawn  **
!**       i = Points connected and a symbol drawn at  **
!**         every ith point            **
!**            -i = Points not connected and a symbol at every  **
!**     mrc     1 = Custom interrupted line style     **
!**     tlen    overall length of the pattern     **
!**     nmrk total number of marks and spaces     **
!**     rat Array of ratios of marks and spaces to overall **
!**             length                 **
!**     icont Contour plotting flag            **
!**       1 = Contour plotting            **
!**       0 = No contour plotting   dim = ncurve     **
!**     nword   Number of words available in common block **
!**     zmat    2-dimensional array containing Z data surface **
!**     values                 **
!**     ix X dimension of zmat            **
!**     iy Y dimension of zmat            **
!**     zinc Increment between z levels     **
!**     line Index number in the group of curve      **
!**     characteristic             **
!**     mode 'DOT' for dotted lines            **
!**     'DASH' for dashed lines       **
!**     lbflg 'NOLABELS' do not use labels      **
!**     'LABELS' use labels            **
!**     ithk Line thickness             **
!**     ipri Line priority             **
!**     draw 'DRAW' draw curves            **
!**     'NODRAW' do not draw curves     **
!**     nshd Number of shades, shade area between 2 curves **
!**     sx Array of x-coordinates     dim = nsxy **
!**     sy Array of y-coordinates     dim = nsxy **
!**     nsxy Number of (sx,sy) pairs            **
!**     sangle angle of shading lines            **
!**     sgap Array of shading gaps            **
!**     ngaps Number of elements in sgaps     **
!**     nvec Number of vectors            **
!**     xfm X value of originating point dim = nvec **
!**     yfm Y value of originating point dim = nvec **
!**     xto X value of point pointed to dim = nvec **
!**     yto Y value of point pointed to dim = nvec **
!**     ivec Describes vector arrowhead dim = nvec **
!**     m Number of text lines to be written     **
!**     note Text string flag     dim = m     **
!**       1 = Text is character string     **
!**       2 = Text is character string drawn relative  **
!**           to the physical origin     **
!**       3 = Text is real number string     **
!**       4 = Text is real number string drawn relative  **
!**           to the physical origin     **
!**       5 = Text is integer string     **
!**       6 = Text is integer string drawn relative  **
!**           to the physical origin     **
!**     lmes Character string text     dim = m     **
!**     imes Number of characters in lmes dim = m     **
!**     anum Real number text string     dim = m     **
!**     iplce Number of decimal place     dim = m     **
!**     inum Integer number text string dim = m     **
!**     xmpos X-position of text string dim = m     **
!**     ympos Y-position of text string dim = m     **
!**     hgt Font size of the text string dim = m     **
!**                                                                   **
!**     REFERENCES:                  **
!**          (1)  CA-DISSPLA manual              **
!**                                                                   **
!**     RECORD OF MODIFICATION:                                       **
!**          04/19/93..........first created                          **
!**                                                                   **
!*************************************************************************
      subroutine curve2d(ncurve, ipag, ibrdr, grce, xphy, yphy, iorel, &
           xorl, yorl, hight, bngle, bshft, ptitle, pltlen, xtitle, &
           xnlen, ytitle, ynlen, xlen, ylen, xorg, xstp, xmax, yorg, &
           ystp, ymax, iaxis, xtck, ytck, ixnon, iynon, intax, intay, &
           isaxs, sorg, stp, smax, slen, sname, nslen, xpos, ypos, &
           igridx, igridy, idash, idot, ichdsh, ichdot,  &
           thcrv, sclpc, dashme, dotme, chdhme, chdtme, markme, &
           clearx, mrc, tlen, nmrk, rat, x, y, nplt, ncnct,  &
           icont, nword, zmat, ix, iy, zinc, line, mode, &
           lbflg, ithk, ipri, nline, draw, &
           nshd, sx, sy, nsxy, &
           sangle, sgap, ngaps, nvec, xfm, yfm, xto, yto, ivec, &
           m, note, lmes, imes, anum, iplce, inum, xmpos, ympos, hgt, &
           iexit)
      include 'curve2d_var.inc'

      if (m_write .eq. 1) then
!-----------------------------------------------------------------------
!       Write curve2d parameters in ASCII format
!-----------------------------------------------------------------------
        call curve2d_asci(ncurve, ipag, ibrdr, grce, xphy, yphy, iorel, &
           xorl, yorl,hight, bngle, bshft, ptitle, pltlen,  xtitle,  &
           xnlen, ytitle, ynlen, xlen, ylen, xorg, xstp, xmax, yorg, &
           ystp, ymax, iaxis, xtck, ytck, ixnon, iynon, intax, intay, &
           isaxs, sorg, stp, smax, slen, sname, nslen,xpos,ypos, &
           igridx, igridy, idash, idot, ichdsh, ichdot,  &
           thcrv, sclpc, dashme, dotme, chdhme, chdtme, markme, &
           clearx, mrc, tlen, nmrk, rat, x, y, nplt, ncnct,  &
           icont, nword, zmat, ix, iy, zinc, line, mode, &
           lbflg, ithk, ipri, nline, draw, &
           nshd, sx, sy, nsxy, &
           sangle, sgap, ngaps, nvec, xfm, yfm, xto, yto, ivec, &
           m, note, lmes, imes, anum, iplce, inum, xmpos, ympos, hgt, &
           iexit)

      elseif (m_write .eq. 0) then
!-----------------------------------------------------------------------
!       Write curve2d parameters in BINARY format
!-----------------------------------------------------------------------
        call curve2d_bin (ncurve, ipag, ibrdr, grce, xphy, yphy, iorel, &
           xorl, yorl,hight, bngle, bshft, ptitle, pltlen,  xtitle, &
           xnlen, ytitle, ynlen, xlen, ylen, xorg, xstp, xmax, yorg, &
           ystp, ymax, iaxis, xtck, ytck, ixnon, iynon, intax, intay, &
           isaxs, sorg, stp, smax, slen, sname, nslen,xpos,ypos,  &
           igridx, igridy, idash, idot, ichdsh, ichdot,  &
           thcrv, sclpc, dashme, dotme, chdhme, chdtme, markme, &
           clearx, mrc, tlen, nmrk, rat, x, y, nplt, ncnct,  &
           icont, nword, zmat, ix, iy, zinc, line, mode, &
           lbflg, ithk, ipri, nline, draw, &
           nshd, sx, sy, nsxy, &
           sangle, sgap, ngaps, nvec, xfm, yfm, xto, yto, ivec, &
            m, note, lmes, imes, anum, iplce, inum, xmpos, ympos, hgt, &
           iexit)
      endif
      return
      end subroutine curve2d

!*************************************************************************
!**                                                                   **
!**   SUBPROGRAM DESCRIPTION:                                       **
!**     Writes plot parameters in ASCII format           **
!**                                                                   **
!**   CALLING ARGUMENTS:                                            **
!**     ncurve  Number of curves            **
!**     ipag page dimension flag            **
!**       0 = page dimension 11 x 8.5     **
!**       1 = page dimension 8.5 x 11     **
!**     ibrdr   Page border flag            **
!**       1 = Suppress page border     **
!**       0 = Suppress page border     **
!**     grce >= 0 Enable grace margin     **
!**     xphy X-coordinate of physical origin     **
!**     yphy Y-coordinate of physical origin     **
!**     hight   Font size             **
!**     bngle   Base rotation angle            **
!**     bshft   Base translation            **
!**     ptitle  Plot title             **
!**     pltlen  Length of plot title                        **
!**     xtitle  X-axis name                              **
!**     xnlen   Length of x-axis name                           **
!**     ytitle  Y-axis name                                **
!**     ynlen   Length of y-axis name            **
!**     xlen Length of x-axis legend            **
!**     ylen    Length of y-axis legend                         **
!**     iaxis   Axes flag             **
!**       0 = Linear axis             **
!**       1 = Log-Linear axis            **
!**       2 = Linear-Log axis            **
!**       3 = Logarithmic axis            **
!**     xtck    X-axis tick marks            **
!**     ytck    Y-axis tick marks            **
!**     ixnon   X-axis tick marks or labels flag     **
!**       1 = Suppress x-axis tick marks or labels **
!**       0 = X-axis tick marks or labels     **
!**     iynon   Y-axis tick marks or labels flag     **
!**       1 = Suppress y-axis tick marks or labels **
!**       0 = Y-axis tick marks or labels     **
!**     intax   Trailing zeroes on x-axis flag     **
!**       1 = Suppress trailing zeroes on x-axis      **
!**       0 = Trailing zeroes on x-axis      **
!**     intay   Trailing zeroes on y-axis flag     **
!**       1 = Suppress trailing zeroes on y-axis      **
!**       0 = Trailing zeroes on y-axis      **
!**     xorg    X-value at the physical origin     **
!**     xstp    X step size in units            **
!**     xmax    Value at x-axis limit            **
!**     yorg    X-value at the physical origin                 **
!**     ystp    X step size in units            **
!**     ymax    Value at y-axis limit            **
!**     iorel   New physical origin relative to current  **
!**     origin flag             **
!**       1 = New physical origin relative to current  **
!**           origin             **
!**       0 = No new physical origin relative to current  **
!**           origin             **
!**     xorl    X-coordinate of relative origin     **
!**     yorl    Y-coordinate of relative origin     **
!**     igridx  Number of grid lines per step on x-axis       **
!**     igridy  Number of grid lines per step on y-axis       **
!**     idash   Dash grid lines flag            **
!**       1 = Dash grid lines            **
!**       0 = No dash grid lines            **
!**     idot    Dot grid lines flag            **
!**     1 = Dot grid lines            **
!**     0 = No dot grid lines            **
!**     ichdot  Chain dot grid lines flag                       **
!**       1 = Chain dot grid lines                        **
!**       0 = No chain dot grid lines                     **
!**     ichdsh  Chain dash grid lines flag                      **
!**       1 = Chain dash grid lines                       **
!**       0 = No chain dash grid lines                    **
!**     thcrv   Curve thickness     dim = ncurve       **
!**     sclpc   Scale curve marker  dim = ncurve            **
!**     dashme  Dash curve flag     dim = ncurve     **
!**       1 = Dash curve             **
!**       0 = No dash curve            **
!**     dotme   Dot curve flag     dim = ncurve     **
!**       1 = Dot curve             **
!**       0 = No dot curve            **
!**     chdhme  Chain dash curve flag dim = ncurve     **
!**       1 = Chain dash curve            **
!**       0 = No chain dash curve            **
!**     chdtme  Chain ot curve flag dim = ncurve     **
!**       1 = Chain dot curve            **
!**       0 = No chain dot curve            **
!**     markme  Curve marker          dim = ncurve            **
!**     clearx  Color     dim = ncurve     **
!**     x       Array of x-coordinates            **
!**             dim = (nplt, ncurve) **
!**     y       Array of y-coordinates            **
!**                 dim = (nplt, ncurve) **
!**     nplt    Number of (x,y) points to be ploted      **
!**                 dim = ncurve     **
!**     ncnct   Marker specification            **
!**       0 = Points connected with no symbols drawn  **
!**       i = Points connected and a symbol drawn at  **
!**         every ith point            **
!**            -i = Points not connected and a symbol at every  **
!**     mrc     1 = Custom interrupted line style     **
!**     tlen    overall length of the pattern     **
!**     nmrk total number of marks and spaces     **
!**     rat Array of ratios of marks and spaces to overall **
!**             length                 **
!**     icont Contour plotting flag            **
!**       1 = Contour plotting            **
!**       0 = No contour plotting            **
!**     nword   Number of words available in common block **
!**     zmat    2-dimensional array containing Z data surface **
!**     values                 **
!**     ix X dimension of zmat            **
!**     iy Y dimension of zmat            **
!**     zinc Increment between z levels     **
!**     line Index number in the group of curve      **
!**     characteristic             **
!**     mode 'DOT' for dotted lines            **
!**     'DASH' for dashed lines       **
!**     lbflg 'NOLABELS' do not use labels      **
!**     'LABELS' use labels            **
!**     ithk Line thickness             **
!**     ipri Line priority             **
!**     draw 'DRAW' draw curves            **
!**     'NODRAW' do not draw curves     **
!**     nshd Number of shades, shade area between 2 curves **
!**     sx Array of x-coordinates     dim = nsxy **
!**     sy Array of y-coordinates     dim = nsxy **
!**     nsxy Number of (sx,sy) pairs            **
!**     sangle angle of shading lines            **
!**     sgap Array of shading gaps            **
!**     ngaps Number of elements in sgaps     **
!**     nvec Number of vectors            **
!**     xfm X value of originating point dim = nvec **
!**     yfm Y value of originating point dim = nvec **
!**     xto X value of point pointed to dim = nvec **
!**     yto Y value of point pointed to dim = nvec **
!**     ivec Describes vector arrowhead dim = nvec **
!**     m Number of text lines to be written     **
!**     note Text string flag     dim = m     **
!**       1 = Text is character string     **
!**       2 = Text is character string drawn relative  **
!**           to the physical origin     **
!**       3 = Text is real number string     **
!**       4 = Text is real number string drawn relative  **
!**           to the physical origin     **
!**       5 = Text is integer string     **
!**       6 = Text is integer string drawn relative  **
!**           to the physical origin     **
!**     lmes Character string text     dim = m     **
!**     imes Number of characters in lmes dim = m     **
!**     anum Real number text string     dim = m     **
!**     iplce Number of decimal place     dim = m     **
!**     inum Integer number text string dim = m     **
!**     xmpos X-position of text string dim = m     **
!**     ympos Y-position of text string dim = m     **
!**     hgt Font size of the text string dim = m     **
!**                                                                   **
!**     REFERENCES:                  **
!**          (1)  CA-DISSPLA manual              **
!**                                                                   **
!**     RECORD OF MODIFICATION:                                       **
!**          04/19/93..........first created                          **
!**                                                                   **
!*************************************************************************
      subroutine curve2d_asci(ncurve, ipag, ibrdr, grce, xphy, yphy, iorel, &
           xorl, yorl, hight, bngle, bshft, ptitle, pltlen, xtitle, &
           xnlen, ytitle, ynlen, xlen, ylen, xorg, xstp, xmax, yorg, &
           ystp, ymax, iaxis, xtck, ytck, ixnon, iynon, intax, intay, &
           isaxs, sorg, stp, smax, slen, sname, nslen, xpos, ypos, &
           igridx, igridy, idash, idot, ichdsh, ichdot,  &
           thcrv, sclpc, dashme, dotme, chdhme, chdtme, markme, &
           clearx, mrc, tlen, nmrk, rat, x, y, nplt, ncnct,  &
           icont, nword, zmat, ix, iy, zinc, line, mode, &
           lbflg, ithk, ipri, nline, draw, &
           nshd, sx, sy, nsxy, &
           sangle, sgap, ngaps, nvec, xfm, yfm, xto, yto, ivec, &
           m, note, lmes, imes, anum, iplce, inum, xmpos, ympos, hgt, &
           iexit)
      include 'curve2d_var.inc'

!-----------------------------------------------------------------------
!     Write page parameters
!-----------------------------------------------------------------------
      write(iunit,100) ncurve, ipag, ibrdr, grce, xphy, yphy
100   format(3i5, 3f9.3)
      write(iunit,101) iorel, xorl, yorl
101   format(i5, 2f9.3)
      write(iunit,102) hight, bngle, bshft
102   format(4f9.3)

!-----------------------------------------------------------------------
!     Write title parameters
!-----------------------------------------------------------------------
      write(iunit,103) ptitle, pltlen,  xtitle, xnlen
103   format(a20, i5, a20, i5)
      write(iunit,104) ytitle, ynlen, xlen, ylen
104   format(a20, i5, 2f9.3)

!-----------------------------------------------------------------------
!     Write axis parameters
!-----------------------------------------------------------------------
      write(iunit,105) xorg, xstp, xmax
      write(iunit,105) yorg, ystp, ymax
105   format(3f15.5)
      write(iunit,106) iaxis, xtck, ytck, ixnon, iynon, intax, intay
106   format(7i5)

!-----------------------------------------------------------------------
!     Write secondary axis parameters
!-----------------------------------------------------------------------
      write(iunit,117) isaxs
      write(iunit,107) sorg, stp, smax, slen, sname
107   format(4f15.5, a20)
      write(iunit,108) nslen,xpos,ypos
108   format(i5, 2f9.3)

!-----------------------------------------------------------------------
!     Write grid parameters
!-----------------------------------------------------------------------
      write(iunit,109) igridx, igridy, idash, idot, ichdsh, ichdot
109   format(6i5)

!-----------------------------------------------------------------------
!     Do for ncurve
!-----------------------------------------------------------------------
      if (ncurve.gt.0) then
        do i = 1, ncurve
!-----------------------------------------------------------------------
!         Write curve parameters
!-----------------------------------------------------------------------
          write(iunit,110) thcrv(i), sclpc(i), dashme(i), dotme(i)
110       format(2f9.3,2i5)
          write(iunit,111) chdhme(i), chdtme(i), markme(i), clearx(i)
111       format(3i5, a20)
          write(iunit,113) mrc(i), tlen(i), nmrk(i)
113       format(i5, f9.3, i5)
          do j = 1, nmrk(i), 4
            write(iunit, 114)  (rat(k), k=j,j+3)
114         format (4f9.3)
          enddo
          write(iunit,115) nplt(i), ncnct(i)
115       format(2i5)
!-----------------------------------------------------------------------
!         Write (X,Y) co-ordinates
!-----------------------------------------------------------------------
          if (nplt(i).gt.0) then
            write(iunit,116) (x(j,i), y(j,i),j=1,nplt(i))
116         format (1x,6e12.5)
          endif    !   (X,Y)

          write(iunit,117) icont(i)
117       format(i5)
!-----------------------------------------------------------------------
!         Write contour parameters
!-----------------------------------------------------------------------
          if (icont(i).gt.0) then
            write(iunit,118) nword, ix, iy, zinc, line, mode
118         format(3i5, f15.5, i5, a20)
            write(iunit,119) ((zmat(k,j),k=1,ix),j=1,iy)
119         format (1x,6e12.5)
            write(iunit,120) lbflg, ithk, ipri, nline, draw
120         format(a20, 3i5, a20)
          endif         !     contour

          write(iunit,121) nshd
121       format(i5)

          if (nshd.gt.0) then
!-----------------------------------------------------------------------
!           Write shade parameters
!-----------------------------------------------------------------------
            do j = 1, nshd
              write(iunit,122) nsxy(j)
122           format (i5)
              write(iunit, 123) (sx(k,j), sy(k,j), k = 1, nsxy(j))
123           format (1x,10e10.3)
              write(iunit,124) sangle(j), sgap(j), ngaps(j)
124           format (1x,2e10.3, i5)
            enddo
          endif      !  shade

          write(iunit,125) nvec
125       format(i5)
!-----------------------------------------------------------------------
!         Write vector parameters
!-----------------------------------------------------------------------
          if (nvec.gt.0) then
            do j = 1,nvec
              write(iunit,126) xfm(j), yfm(j), xto(j), yto(j), ivec(j)
126           format(1x,4e10.3, i5)
            enddo
          endif      ! vector
!-----------------------------------------------------------------------
!       EndDo for ncurve
!-----------------------------------------------------------------------
        enddo
      endif    !  ncurve

      write(iunit,127) m
127   format(i5)
      if (m .gt. 0) then
        do i = 1, m
!-----------------------------------------------------------------------
!         Write annotation parameters
!-----------------------------------------------------------------------
          write(iunit,128) note(i), imes(i)
128       format(i5, i5)
          write(iunit,132) lmes(i)
132       format(a72)
          write(iunit,129) anum(i), iplce(i), inum(i)
129       format(f15.5, 2i5)
          write(iunit,130) xmpos(i), ympos(i), hgt(i)
130       format(3f9.3)
        enddo
      endif

!-----------------------------------------------------------------------
!     Write plot exit parameters
!-----------------------------------------------------------------------
      write(iunit,131) iexit
131   format(i5)

      return
      end subroutine curve2d_asci

!*************************************************************************
!**                                                                   **
!**   SUBPROGRAM DESCRIPTION:                                       **
!**     Writes plot parameters in BINARY format           **
!**                                                                   **
!**   CALLING ARGUMENTS:                                            **
!**     ncurve  Number of curves            **
!**     ipag page dimension flag            **
!**       0 = page dimension 11 x 8.5     **
!**       1 = page dimension 8.5 x 11     **
!**     ibrdr   Page border flag            **
!**       1 = Suppress page border     **
!**       0 = Suppress page border     **
!**     grce >= 0 Enable grace margin     **
!**     xphy X-coordinate of physical origin     **
!**     yphy Y-coordinate of physical origin     **
!**     hight   Font size             **
!**     bngle   Base rotation angle            **
!**     bshft   Base translation            **
!**     ptitle  Plot title             **
!**     pltlen  Length of plot title                        **
!**     xtitle  X-axis name                              **
!**     xnlen   Length of x-axis name                           **
!**     ytitle  Y-axis name                                **
!**     ynlen   Length of y-axis name            **
!**     xlen Length of x-axis legend            **
!**     ylen    Length of y-axis legend                         **
!**     iaxis   Axes flag             **
!**       0 = Linear axis             **
!**       1 = Log-Linear axis            **
!**       2 = Linear-Log axis            **
!**       3 = Logarithmic axis            **
!**     xtck    X-axis tick marks            **
!**     ytck    Y-axis tick marks            **
!**     ixnon   X-axis tick marks or labels flag     **
!**       1 = Suppress x-axis tick marks or labels **
!**       0 = X-axis tick marks or labels     **
!**     iynon   Y-axis tick marks or labels flag     **
!**       1 = Suppress y-axis tick marks or labels **
!**       0 = Y-axis tick marks or labels     **
!**     intax   Trailing zeroes on x-axis flag     **
!**       1 = Suppress trailing zeroes on x-axis      **
!**       0 = Trailing zeroes on x-axis      **
!**     intay   Trailing zeroes on y-axis flag     **
!**       1 = Suppress trailing zeroes on y-axis      **
!**       0 = Trailing zeroes on y-axis      **
!**     xorg    X-value at the physical origin     **
!**     xstp    X step size in units            **
!**     xmax    Value at x-axis limit            **
!**     yorg    X-value at the physical origin                 **
!**     ystp    X step size in units            **
!**     ymax    Value at y-axis limit            **
!**     iorel   New physical origin relative to current  **
!**     origin flag             **
!**       1 = New physical origin relative to current  **
!**           origin             **
!**       0 = No new physical origin relative to current  **
!**           origin             **
!**     xorl    X-coordinate of relative origin     **
!**     yorl    Y-coordinate of relative origin     **
!**     igridx  Number of grid lines per step on x-axis       **
!**     igridy  Number of grid lines per step on y-axis       **
!**     idash   Dash grid lines flag            **
!**       1 = Dash grid lines            **
!**       0 = No dash grid lines            **
!**     idot    Dot grid lines flag            **
!**       1 = Dot grid lines            **
!**       0 = No dot grid lines            **
!**     ichdot  Chain dot grid lines flag                       **
!**       1 = Chain dot grid lines                        **
!**       0 = No chain dot grid lines                     **
!**     ichdsh  Chain dash grid lines flag                      **
!**       1 = Chain dash grid lines                       **
!**       0 = No chain dash grid lines                    **
!**     thcrv   Curve thickness     dim = ncurve       **
!**     sclpc   Scale curve marker  dim = ncurve            **
!**     dashme  Dash curve flag     dim = ncurve     **
!**       1 = Dash curve             **
!**       0 = No dash curve            **
!**     dotme   Dot curve flag     dim = ncurve     **
!**       1 = Dot curve             **
!**       0 = No dot curve            **
!**     chdhme  Chain dash curve flag dim = ncurve     **
!**       1 = Chain dash curve            **
!**       0 = No chain dash curve            **
!**     chdtme  Chain ot curve flag dim = ncurve     **
!**       1 = Chain dot curve            **
!**       0 = No chain dot curve            **
!**     markme  Curve marker          dim = ncurve            **
!**     clearx  Color     dim = ncurve     **
!**     x       Array of x-coordinates            **
!**                dim = (nplt, ncurve) **
!**     y       Array of y-coordinates            **
!**                dim = (nplt, ncurve) **
!**     nplt    Number of (x,y) points to be ploted      **
!**                dim = ncurve     **
!**     ncnct   Marker specification            **
!**       0 = Points connected with no symbols drawn  **
!**       i = Points connected and a symbol drawn at  **
!**         every ith point            **
!**            -i = Points not connected and a symbol at every  **
!**     mrc     1 = Custom interrupted line style     **
!**     tlen    overall length of the pattern     **
!**     nmrk total number of marks and spaces     **
!**     rat Array of ratios of marks and spaces to overall **
!**             length                 **
!**     icont Contour plotting flag            **
!**       1 = Contour plotting            **
!**       0 = No contour plotting            **
!**     nword   Number of words available in common block **
!**     zmat    2-dimensional array containing Z data surface **
!**     values                 **
!**     ix X dimension of zmat            **
!**     iy Y dimension of zmat            **
!**     zinc Increment between z levels     **
!**     line Index number in the group of curve      **
!**     characteristic             **
!**     mode 'DOT' for dotted lines            **
!**     'DASH' for dashed lines       **
!**     lbflg 'NOLABELS' do not use labels      **
!**     'LABELS' use labels            **
!**     ithk Line thickness             **
!**     ipri Line priority             **
!**     draw 'DRAW' draw curves            **
!**     'NODRAW' do not draw curves     **
!**     nshd Number of shades, shade area between 2 curves **
!**     sx Array of x-coordinates     dim = nsxy **
!**     sy Array of y-coordinates     dim = nsxy **
!**     nsxy Number of (sx,sy) pairs            **
!**     sangle angle of shading lines            **
!**     sgap Array of shading gaps            **
!**     ngaps Number of elements in sgaps     **
!**     nvec Number of vectors            **
!**     xfm X value of originating point dim = nvec **
!**     yfm Y value of originating point dim = nvec **
!**     xto X value of point pointed to dim = nvec **
!**     yto Y value of point pointed to dim = nvec **
!**     ivec Describes vector arrowhead dim = nvec **
!**     m Number of text lines to be written     **
!**     note Text string flag     dim = m     **
!**       1 = Text is character string     **
!**       2 = Text is character string drawn relative  **
!**           to the physical origin     **
!**       3 = Text is real number string     **
!**       4 = Text is real number string drawn relative  **
!**           to the physical origin     **
!**       5 = Text is integer string     **
!**       6 = Text is integer string drawn relative  **
!**           to the physical origin     **
!**     lmes Character string text     dim = m     **
!**     imes Number of characters in lmes dim = m     **
!**     anum Real number text string     dim = m     **
!**     iplce Number of decimal place     dim = m     **
!**     inum Integer number text string dim = m     **
!**     xmpos X-position of text string dim = m     **
!**     ympos Y-position of text string dim = m     **
!**     hgt Font size of the text string dim = m     **
!**                                                                   **
!**     REFERENCES:                  **
!**          (1)  CA-DISSPLA manual              **
!**                                                                   **
!**     RECORD OF MODIFICATION:                                       **
!**          03/26/93..........first created                          **
!**                                                                   **
!*************************************************************************
      subroutine curve2d_bin(ncurve, ipag, ibrdr, grce, xphy, yphy, iorel, &
           xorl, yorl, hight, bngle, bshft, ptitle, pltlen, xtitle, &
           xnlen, ytitle, ynlen, xlen, ylen, xorg, xstp, xmax, yorg, &
           ystp, ymax, iaxis, xtck, ytck, ixnon, iynon, intax, intay, &
           isaxs, sorg, stp, smax, slen, sname, nslen, xpos, ypos, &
           igridx, igridy, idash, idot, ichdsh, ichdot,  &
           thcrv, sclpc, dashme, dotme, chdhme, chdtme, markme, &
           clearx, mrc, tlen, nmrk, rat, x, y, nplt, ncnct,  &
           icont, nword, zmat, ix, iy, zinc, line, mode, &
           lbflg, ithk, ipri, nline, draw, &
           nshd, sx, sy, nsxy, &
           sangle, sgap, ngaps, nvec, xfm, yfm, xto, yto, ivec, &
           m, note, lmes, imes, anum, iplce, inum, xmpos, ympos, hgt, &
           iexit)
      include 'curve2d_var.inc'

!-----------------------------------------------------------------------
!     Write page parameters
!-----------------------------------------------------------------------
      write(iunit) ncurve, ipag, ibrdr, grce, xphy, yphy
      write(iunit) iorel, xorl, yorl
      write(iunit) hight, bngle, bshft

!-----------------------------------------------------------------------
!     Write title parameters
!-----------------------------------------------------------------------
      write(iunit) ptitle, pltlen,  xtitle, xnlen
      write(iunit) ytitle, ynlen, xlen, ylen

!-----------------------------------------------------------------------
!     Write axis parameters
!-----------------------------------------------------------------------
      write(iunit) xorg, xstp, xmax
      write(iunit) yorg, ystp, ymax
      write(iunit) iaxis, xtck, ytck, ixnon, iynon, intax, intay

!-----------------------------------------------------------------------
!     Write secondary axis parameters
!-----------------------------------------------------------------------
      write(iunit) isaxs
      write(iunit) sorg, stp, smax, slen, sname
      write(iunit) nslen,xpos,ypos

!-----------------------------------------------------------------------
!     Write grid parameters
!-----------------------------------------------------------------------
      write(iunit) igridx, igridy, idash, idot, ichdsh, ichdot

!-----------------------------------------------------------------------
!     Do for ncurve
!-----------------------------------------------------------------------
      if (ncurve.gt.0) then
        do i = 1, ncurve
!-----------------------------------------------------------------------
!         Write curve parameters
!-----------------------------------------------------------------------
          write(iunit) thcrv(i), sclpc(i), dashme(i), dotme(i)
          write(iunit) chdhme(i), chdtme(i), markme(i), clearx(i)
          write(iunit) mrc(i), tlen(i), nmrk(i)
          do j = 1, nmrk(i), 4
            write(iunit)  (rat(k), k=j,j+3)
          enddo
          write(iunit) nplt(i), ncnct(i)
!-----------------------------------------------------------------------
!         Write (X,Y) coordinates
!-----------------------------------------------------------------------
          if (nplt(i).gt.0) then
            write(iunit) (x(j,i), y(j,i), j=1, nplt(i))
          endif
          write(iunit) icont(i)

!-----------------------------------------------------------------------
!         Write contour parameters
!-----------------------------------------------------------------------
          if (icont(i).gt.0) then
            write (iunit) nword, ix, iy, zinc, line, mode
            write (iunit) ((zmat(k,j) ,k=1,ix),j=1,iy)
            write (iunit) lbflg, ithk, ipri, nline, draw
          endif       !   contour

          write(iunit) nshd
          if (nshd.gt.0) then
            do j = 1, nshd
!-----------------------------------------------------------------------
!             Write shade parameters
!-----------------------------------------------------------------------
              write(iunit) nsxy(j)
              write(iunit) (sx(k,j), sy(k,j), k = 1, nsxy(j))
              write(iunit) sangle(j), sgap(j), ngaps(j)
            enddo
          endif    !   shade

          write(iunit) nvec
          if (nvec.gt.0) then
            do j = 1,nvec
!-----------------------------------------------------------------------
!             Write vector parameters
!-----------------------------------------------------------------------
              write(iunit) xfm(j), yfm(j), xto(j), yto(j), ivec(j)
            enddo
          endif
!-----------------------------------------------------------------------
!       Do for ncurve
!-----------------------------------------------------------------------
        enddo
      endif   !    ncurve

      write(iunit) m
      if (m .gt. 0) then
        do i = 1, m
!-----------------------------------------------------------------------
!         Write annotation parameters
!-----------------------------------------------------------------------
          write(iunit) note(i), imes(i)
          write(iunit) lmes(i)
          write(iunit) anum(i), iplce(i), inum(i)
          write(iunit) xmpos(i), ympos(i), hgt(i)
        enddo
      endif

!-----------------------------------------------------------------------
!     Write plot exit parameters
!-----------------------------------------------------------------------
      write(iunit) iexit

      return
      end subroutine curve2d_bin

!***********************************************************************
!**                                                                   **
!**     SUBPROGRAM DESCRIPTION:                                       **
!**          Initializes plot parameters                              **
!**                                                                   **
!**     RECORD OF MODIFICATION:                                       **
!**          04/28/93..........first created                          **
!**                                                                   **
!***********************************************************************
      subroutine init2d
      use set_kinds, only: dp
      use curve2d_mod
      implicit integer*4 (i-n), real*8 (a-h, o-z)

!-----------------------------------------------------------------------
!     Initialize plot parameters
!-----------------------------------------------------------------------
      nn = 0
      ncurve = 1
!-----------------------------------------------------------------------
!     Initialize page size to 11x8.5
!-----------------------------------------------------------------------
      ipag = 0
!-----------------------------------------------------------------------
!     Initialize grace margin be coincident with the grid frame
!-----------------------------------------------------------------------
      grce = 0.0
!-----------------------------------------------------------------------
!     Initialize page border "do not suppress page border"
!-----------------------------------------------------------------------
      ibrdr = 0
!-----------------------------------------------------------------------
!     Initialize page dimension
!-----------------------------------------------------------------------
      xpag = 0.0
      ypag = 0.0
!-----------------------------------------------------------------------
!     Initialize physical origin parameters
!-----------------------------------------------------------------------
      xphy = 1.0
      yphy = 1.0
!-----------------------------------------------------------------------
!     Initialize relative origin parameters
!-----------------------------------------------------------------------
      iorel = 0
      xorl = 0.0
      yorl = 0.0
!-----------------------------------------------------------------------
!     Initialize font size
!-----------------------------------------------------------------------
      hight = 0.5_dp
!-----------------------------------------------------------------------
!     Initialize base rotation
!-----------------------------------------------------------------------
      bngle = 0.0
!-----------------------------------------------------------------------
!     Initialize translation
!-----------------------------------------------------------------------
      bshft(1) = 0.0
      bshft(2) = 0.0
!-----------------------------------------------------------------------
!     Initialize type of coordinate system (Cartesian coordinate)
!-----------------------------------------------------------------------
      iaxis = 0
!-----------------------------------------------------------------------
!     Initialize axes tick marks off
!-----------------------------------------------------------------------
      ixtck = 0
      iytck = 0
!-----------------------------------------------------------------------
!     Initialize axes labels not to be suppressed
!-----------------------------------------------------------------------
      ixnon = 0
      iynon = 0
!-----------------------------------------------------------------------
!     Initialize trailing zeroes on axes
!-----------------------------------------------------------------------
      intax = 0
      intay = 0
!-----------------------------------------------------------------------
!     Initialize secondary axes parameters
!-----------------------------------------------------------------------
      isaxs = 0
      sorg = 0.0
      stp = 0.0
      smax = 0.0
      slen = 0.0
      sname = ' '
      nslen = 0
      xps = 0.0
      yps = 0.0
!-----------------------------------------------------------------------
!     Initialize plot title legend
!-----------------------------------------------------------------------
      ptitle = '$'
      nplen = -100
!-----------------------------------------------------------------------
!     Initialize x-axis title legend
!-----------------------------------------------------------------------
      xtitle = '$'
      nxlen = 100
!-----------------------------------------------------------------------
!     Initialize y-axis title legend
!-----------------------------------------------------------------------
      ytitle = '$'
      nylen = 100
!-----------------------------------------------------------------------
!     Initialize grids to be turned off
!-----------------------------------------------------------------------
      igridx = 0
      igridy = 0
!-----------------------------------------------------------------------
!     Initialize, grid dashes off
!-----------------------------------------------------------------------
      idash = 0
!-----------------------------------------------------------------------
!     Initialize grid dots off
!-----------------------------------------------------------------------
      idot = 0
!-----------------------------------------------------------------------
!     Initialize grid chain dot off
!-----------------------------------------------------------------------
      ichdot = 0
!-----------------------------------------------------------------------
!     Initialize grid chain dash off
!-----------------------------------------------------------------------
      ichdsh = 0
!-----------------------------------------------------------------------
!     Initialize number of text lines to zero
!-----------------------------------------------------------------------
      msg = 0
!-----------------------------------------------------------------------
!     Initialize number of shades to zero
!-----------------------------------------------------------------------
      nshd = 0
!-----------------------------------------------------------------------
!     Initialize number of vectors to zero
!-----------------------------------------------------------------------
      nvec = 0
!-----------------------------------------------------------------------
!     Initialize contour drawing to none
!-----------------------------------------------------------------------
      do i = 1, ncrv
        icont(i) = 0
      enddo
!-----------------------------------------------------------------------
!     Initialize curve drawing parameters
!-----------------------------------------------------------------------
      do i = 1, ncrv
!-----------------------------------------------------------------------
!       Initialize thickness of curve
!-----------------------------------------------------------------------
        thcrv(i) = 0.0
!-----------------------------------------------------------------------
!       Initialize scaling of marker
!-----------------------------------------------------------------------
        sclpc(i) = 0.0
!-----------------------------------------------------------------------
!       Initialize curve dot off
!-----------------------------------------------------------------------
        ndotme(i) = 0
!-----------------------------------------------------------------------
!       Initialize curve dash off
!-----------------------------------------------------------------------
        ndshme(i) = 0
!-----------------------------------------------------------------------
!       Initialize curve chain dot off
!-----------------------------------------------------------------------
        ncdtme(i) = 0
!-----------------------------------------------------------------------
!       Initialize curve chain dash off
!-----------------------------------------------------------------------
        ncdhme(i) = 0
!-----------------------------------------------------------------------
!       Initialize marker
!-----------------------------------------------------------------------
        markme(i) = -2
!-----------------------------------------------------------------------
!       Initialize curve to be drawn with solid line wit no symbols
!-----------------------------------------------------------------------
        ncnct(i) = 0
!-----------------------------------------------------------------------
!       Initialize color
!-----------------------------------------------------------------------
        clearx(i) = 'FOREGROUND'
!-----------------------------------------------------------------------
!       Initialize shade parameters
!-----------------------------------------------------------------------
!        sangle(i) = 0.0
!        sgap(i) = 0.0
!        ngaps(i) = 0
!        mrc(i) = 0
      enddo
!-----------------------------------------------------------------------
!     Initialize Annotation parameters
!-----------------------------------------------------------------------
      do i = 1, mdim
        note(i) = 0
        lmes(i) = ' '
        imes(i) = 0
        anum(i) = 0.0
        iplce(i) = 0
        inum(i) = 0
        xpos(i) = 0.0
        ypos(i) = 0.0
        ht(i) = 0.14_dp
      enddo
!-----------------------------------------------------------------------
!     Initialize Contour Dimension
!-----------------------------------------------------------------------
      ix = 1
      iy = 1

      return
      end subroutine init2d
