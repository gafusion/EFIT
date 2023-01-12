!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          efun generates the necessary response                   **
!**          functions used by efit for reconstruction of the        **
!**          magnetic surfaces and current density profile.          **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1) d.w. swain and g.h. neilson, nucl. fusion           **
!**              22 (1982) 1015.                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**          28/01/85..........modified for D3D                      **
!**                                                                  **
!**********************************************************************
      program efund

      write(*,*) 'Reading namelist'
      call efund_getsizes
      call efund_getset
      write(*,*) 'Calling matrix subroutine'
      call efund_matrix
      write(*,*) 'Calling grid subroutine'
      call efund_grid

      stop 'GREEN TABLE GENERATED!'
      end program efund
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          e1coef computes the response functions due to           **
!**          the thin flux loops and E coils.                        **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          15/07/83..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine e1coef(coef,nl,ne)
      use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,nesum,&
                      nfsum,nvsum,nvesel,nacoil
      use cecoil
      use coilsp
      use consta
      use nio
      use siloop
      implicit integer*4 (i-n), real*8 (a-h, o-z)
!
      radeg=pi/180.
      isplit=17
      itot=isplit*isplit
      fitot=itot
!
      m=nl
      k=ne
      psict=0
      aaa=0.0
      bbb=0.0
      call splitc(isplit,rsplt,zsplt, &
                  re(k),ze(k),we(k),he(k),aaa,bbb)
      do l=1,itot
         a=rsplt(l)
         r1=rsi(m)
         z1=zsi(m)-zsplt(l)
         psic=psical(a,r1,z1)*tmu
         psict=psict+psic/fitot
      enddo 
      coef=psict
!
      return
      end subroutine e1coef
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          e2coef computes the response functions due to           **
!**          the magnetic probes and E coils.                        **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine e2coef(coef,mp,ne)
      use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,nesum,&
                      nfsum,nvsum,nvesel,nacoil
      use cecoil
      use coilsp
      use consta
      use mprobe
      implicit integer*4 (i-n), real*8 (a-h, o-z)
!
      radeg=pi/180.
      isplit=17
      itot=isplit*isplit
      fitot=itot
!
      m=mp
      if (smp2(m).gt.0.0) then
         sinm=sin(radeg*amp2(m))
         cosm=cos(radeg*amp2(m))
         delsx=smp2(m)/nsmp2*cosm
         delsy=smp2(m)/nsmp2*sinm
      else
!------------------------------------------------------------------------------
!--  perpendicular probes    96/02/04                                        --
!------------------------------------------------------------------------------
         sinm=sin(radeg*amp2(m))
         cosm=cos(radeg*amp2(m))
         sinms=sin(radeg*(amp2(m)+90.))
         cosms=cos(radeg*(amp2(m)+90.))
         delsx=abs(smp2(m))/nsmp2*cosms
         delsy=abs(smp2(m))/nsmp2*sinms
      endif
      xmp20=xmp2(m)-(nsmp2-1)/2.*delsx
      ymp20=ymp2(m)-(nsmp2-1)/2.*delsy
      k=ne
      brct=0
      bzct=0
      aaa=0.0
      bbb=0.0
      call splitc(isplit,rsplt,zsplt, &
                  re(k),ze(k),we(k),he(k),aaa,bbb)
      do l=1,itot
         a=rsplt(l)
         do mmm=1,nsmp2
            r1=xmp20+(mmm-1)*delsx
            z1=ymp20+(mmm-1)*delsy-zsplt(l)
            brc=br(a,r1,z1)*tmu
            bzc=bz(a,r1,z1)*tmu
            brct=brct+brc/fitot
            bzct=bzct+bzc/fitot
         enddo 
      enddo 
      coef=(brct*cosm+bzct*sinm)/nsmp2
!
      return
      end subroutine e2coef
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          egrid computes the response functions due to            **
!**          the grid points and E coils.                            **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          15/07/83..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine egrid(coef,rgrid,nr,zgrid,nz,ne)
      use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,nesum,&
                      nfsum,nvsum,nvesel,nacoil
      use cecoil
      use coilsp
      use consta
      use nio
      use mprobe
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      real*8,dimension(nr) :: rgrid
      real*8,dimension(nz) :: zgrid
!
!
      radeg=pi/180.
      isplit=17
      itot=isplit*isplit
      fitot=itot
!
      k=ne
      psict=0
      aaa=0.0
      bbb=0.0
      call splitc(isplit,rsplt,zsplt, &
                  re(k),ze(k),we(k),he(k),aaa,bbb)
      do l=1,itot
         a=rsplt(l)
         r1=rgrid(nr)
         z1=zgrid(nz)-zsplt(l)
         psic=psical(a,r1,z1)*tmu
         psict=psict+psic/fitot
      enddo 
      coef=psict
      return
      end subroutine egrid
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          flux computes the mutual inductance/2/pi between        **
!**          two circulars of rectangular cross section.             **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**       r1,r2...........radius of first and second coil            **
!**       z1,z2...........elevation                                  **
!**       w1,w2...........width                                      **
!**       t1,t12,t2,t22...angle                                      **
!**                                                                  **
!**********************************************************************
      subroutine flux(r1,z1,w1,h1,t1,t12,r2,z2,w2,h2,t2,t22,fuxx)
      use exparm, only: mgaus1,mgaus2
      implicit none
      real*8, intent(in) :: r1,z1,w1,h1,t1,t12,r2,z2,w2,h2,t2,t22
      real*8, intent(out) :: fuxx
      integer*4 i,j,ner
      real*8 x,zf,zs,hf2,hs2,xl1,xr1,xbm1,xtm1,xl2,xr2,xbm2,xtm2,rf,rs, &
             hfa,hsa,solx
      real*8 post1(mgaus1),wght1(mgaus1),post2(mgaus2),wght2(mgaus2)
!      data init/0/
!
!      init=0
!      if (init.le.0) then
!vas introduced to make it ok in hydra
         call lgauss(post1,wght1,mgaus1,ner)
         call lgauss(post2,wght2,mgaus2,ner)
!         init=1
!      endif
!
      x = 0.
      fuxx = 0.
      xbm1 = 0.
      xtm1 = 0.
      xbm2 = 0.
      xtm2 = 0.
      zf = z1
      zs = z2
      hf2 = h1*.5
      hs2 = h2*.5
!
      if (t12.ne.0) then
         xl1 = r1-0.5*w1-0.5*h1/abs(t12)
         xr1 = r1+0.5*w1+0.5*h1/abs(t12)
         xbm1 = xl1+w1
         xtm1 = xr1-w1
         if(t12.lt.0.) xbm1 = xr1 - w1
         if(t12.lt.0.) xtm1 = xl1 + w1
      endif
!
      if (t22.ne.0) then
         xl2 = r2-0.5*w2-0.5*h2/abs(t22)
         xr2 = r2+0.5*w2+0.5*h2/abs(t22)
         xbm2 = xl2+w2
         xtm2 = xr2-w2
         if(t22.lt.0.) xbm2 = xr2 - w2
         if(t22.lt.0.) xtm2 = xl2 + w2
      endif
!
      do i = 1,mgaus1
         rf = r1+.5*w1*post1(i)
         if(t12.ne.0) rf = r1+(0.5*w1+0.5*h1/abs(t12))*post1(i)
         do j = 1,mgaus2
            rs = r2+0.5*w2*post2(j)
            if(t22.ne.0) rs = r2+(0.5*w2+0.5*h2/abs(t22))*post2(j)
            call soleno(r1,z1,w1,h1,t1,t12,r2,z2,w2,h2,t2,t22,xbm1,xbm2, &
                        xtm1,xtm2,hfa,hsa,rf,rs,solx)
            fuxx = fuxx+wght1(i)*wght2(j)/hfa/hsa*solx
         enddo 
      enddo 
!
      return
      end subroutine flux
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          gacoil gets the Green's functions for the advance       **
!**          divertor coil.                                          **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          16/08/90..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine gacoil(rsilac,rmp2ac,gridac,rgrid,mw, zgrid,mh)
      use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,nesum,&
                      nfsum,nvsum,nvesel,nacoil,nw,nh,nwnh
      use consta
      use cacoil
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension rsilac(nsilop,nacoil),rmp2ac(magpr2,nacoil)
      real*8,dimension(mw) ::  rgrid
      real*8,dimension(mh) ::  zgrid
      real*8,dimension(nwnh,nacoil) :: gridac
      do j=1,nsilop
         jj=j
         do i=1,nacoil
            ii=i
            call a1coef(work,jj,ii)
            rsilac(j,i)=work
         enddo 
      enddo 
!
      do  j=1,magpr2
         jj=j
         do  i=1,nacoil
            ii=i
            call a2coef(work,jj,ii)
            rmp2ac(j,i)=work
         enddo 
      enddo 
!
      do i=1,nw
         nr=i
         do j=1,nh
            nz=j
            kk=(i-1)*nh+j
            do n=1,nacoil
               nn=n
               call agrid(work,rgrid,nr,zgrid,nz,nn)
               gridac(kk,n)=work
            enddo 
         enddo
      enddo 
      return
      end subroutine gacoil
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          gecoil gets the Green's functions for E coils.          **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          30/01/85..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine gecoil(rsilec,rmp2ec,gridec,rgrid,mw,zgrid,mh, &
                        rfcec,recec,rsisec)
      use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,nesum,&
                      nfsum,nvsum,nvesel,nacoil,nw,nh,nwnh
      use consta
      use fcoil
      use cecoil
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension rsilec(nsilop,nesum),rmp2ec(magpr2,nesum), &
                rfcec(nfcoil,nesum),recec(nesum,nesum), &
                rsisec(nesum)
            real*8,dimension(mw) :: rgrid
            real*8,dimension(mh) :: zgrid
            real*8,dimension(nwnh,nesum) :: gridec
      dimension taf(nfcoil),taf2(nfcoil)
            real*8 :: zetaec = 3.5e-08
      do j=1,nsilop
         do i=1,nesum
            rsilec(j,i)=0.0
         enddo 
         jj=j
         do i=1,necoil
            ii=i
            call e1coef(work,jj,ii)
            kkm=ecid(i)
            rsilec(j,kkm)=rsilec(j,kkm)+work*ecturn(i)
         enddo 
      enddo
!
      do j=1,magpr2
         do i=1,nesum
            rmp2ec(j,i)=0.0
         enddo 
         jj=j
         do i=1,necoil
            ii=i
            call e2coef(work,jj,ii)
            kkm=ecid(i)
            rmp2ec(j,kkm)=rmp2ec(j,kkm)+work*ecturn(i)
         enddo 
      enddo
!
      do i=1,nw
         nr=i
         do j=1,nh
            nz=j
            kk=(i-1)*nh+j
            do m=1,nesum
               gridec(kk,m)=0.0
            enddo 
            do n=1,necoil
               nn=n
               call egrid(work,rgrid,nr,zgrid,nz,nn)
               kkkm=ecid(n)
               gridec(kk,kkkm)=gridec(kk,kkkm)+work*ecturn(n)
            enddo 
         enddo
      enddo
!
      aaa=0.0
      bbb=0.0
      do j=1,nfcoil
         do i=1,nesum
            rfcec(j,i)=0.0
         enddo 
         taf(j)=tan(af(j)*pi/180.)
         taf2(j)=tan(af2(j)*pi/180.)
         do i=1,necoil
            call flux(re(i),ze(i),we(i),he(i),aaa,bbb, &
                      rf(j),zf(j),wf(j),hf(j),taf(j),taf2(j),work)
            work=work*0.5/pi
            kkm=ecid(i)
            kk=kkm+1
            rfcec(j,kkm)=rfcec(j,kkm)+work*(kk-ecid(i))
         enddo 
      enddo
!
      do j=1,nesum
         rsisec(j)=0.0
         do i=1,nesum
            recec(j,i)=0.0
         enddo 
      enddo 
      do j=1,necoil
         jjm=ecid(j)
         do i=1,necoil
            call flux(re(i),ze(i),we(i),he(i),aaa,bbb, &
                      re(j),ze(j),we(j),he(j),aaa,bbb,work)
            work=work*0.5/pi
            kkm=ecid(i)
            recec(jjm,kkm)=recec(jjm,kkm)+work
         enddo 
         rsisec(jjm)=rsisec(jjm)+2.*pi*re(j)/we(j)/he(j)*zetaec
      enddo
      return
      end subroutine gecoil
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          grid computes the green's functions at (r,z)            **
!**          due to plasma currents and f coils.                     **
!**                                                                  **
!**********************************************************************
      subroutine efund_grid
      use exparm, only: nfcoil,nsilop,magpr2,nrogow,necoil,nesum, &
                        nfsum,nvsum,nvesel,nacoil,nw,nh,nwnh
      use consta
      use pmodel
      use input
      use fcoil
      use nio
      use var_filech
      implicit none
      real*8 psical
      integer*4 i,ih,ii,iw,j,k,kk,ndr,ndz,ni,nj
      real*8 ab,cmut,dhc,drc,dwc,dzc,fridpc,rtmp,rsum,z1,z2,zdif
      real*8, dimension(:,:), allocatable :: gridfc,gridpc,ggridfc
      real*8 taf(nfcoil),taf2(nfcoil)
      integer*4, parameter :: isplit=17
      real*8, parameter :: aaa=0.0,tole=1.0e-10
!
      if (.not. allocated(gridfc)) then
        allocate(gridfc(nwnh,nfcoil))
        gridfc = 0.0
      endif
      if (.not. allocated(gridpc)) then
        allocate(gridpc(nwnh,nw))
        gridpc = 0.0
      endif
      if (.not. allocated(ggridfc)) then
        allocate(ggridfc(nwnh,nfsum))
        ggridfc = 0.0
      endif
!
      if(igrid.eq.0) return
!----------------------------------------------------------------------
!--  compute the green's functions at (r,z) due to f coils           --
!----------------------------------------------------------------------
      ndr = isplit
      ndz = isplit
      do ii=1,nfcoil
         taf(ii) = tan(af(ii)*pi/180.0)
         taf2(ii) = tan(af2(ii)*pi/180.0)
         do ni=1,nw
            do nj=1,nh
               kk = (ni-1)*nh + nj
               rsum = 0.
               dwc = wf(ii)/ndr
               dhc = hf(ii)/ndz
               if (af2(ii) .eq. 0.) then
                  z1 = zf(ii)-taf(ii)*(wf(ii)-dwc)/2.-.5*hf(ii)+.5*dhc
                  ab = rf(ii)-.5*wf(ii)+.5*dwc
                  do iw = 1,isplit
                     drc = ab+(iw-1)*dwc+iw*tole
                     z2 = z1+(iw-1)*taf(ii)*dwc
                     do ih = 1,isplit
                        dzc = z2+(ih-1)*dhc
                        rtmp = psical(drc,rgrid(ni),zgrid(nj)-dzc)
                        rsum = rsum+rtmp
                     enddo
                  enddo
               else
                  do ih = 1,ndz
                     dzc = zf(ii)-.5*hf(ii)+.5*dhc+dhc*(ih-1)
                     do iw = 1,ndr
                        drc = rf(ii)-.5*wf(ii)-.5*hf(ii)/taf2(ii) &
                              +.5*dwc+.5*dhc/taf2(ii) &
                              +dhc/taf2(ii)*(ih-1)+dwc*(iw-1)
                        rtmp = psical(drc,rgrid(ni),zgrid(nj)-dzc)
                        rsum = rsum+rtmp
                     enddo
                  enddo
               endif
               cmut = rsum*2.e-07/(isplit*isplit)
               gridfc(kk,ii) = cmut
            enddo
         enddo
      enddo
!----------------------------------------------------------------------
!--  compute the green's functions at (r,z) due to itself            --
!----------------------------------------------------------------------
      do i=1,nw
         do j=1,nh
            kk=(i-1)*nh+j
            do ni=1,nw
               if ((j.gt.1).or.(i.ne.ni)) then
                  zdif=(j-1)*dz
                  gridpc(kk,ni)=psical(rgrid(i),rgrid(ni),zdif)*tmu
               else
                  call flux(rgrid(ni),aaa,dr,dz,aaa,aaa,rgrid(ni),aaa, &
                            dr,dz,aaa,aaa,fridpc)
                  gridpc(kk,ni)=fridpc*0.5/pi
               endif
            enddo 
         enddo
      enddo
!----------------------------------------------------------------------
!--  store green's function table                                    --
!----------------------------------------------------------------------
      ggridfc=0.0
!
      do i=1,nfcoil
         k=abs(fcid(i))
         ggridfc(:,k)=ggridfc(:,k)+fcturn(i)*gridfc(:,i)
      enddo
!
      print*,'file name : ','ec'//trim(ch1)//trim(ch2)//'.ddd' 
!
!vasorg      open(unit=ncontr,status='unknown',file='econto.dat', &
      open(unit=ncontr,status='unknown',file='ec'//trim(ch1)// &
           trim(ch2)//'.ddd',form='unformatted')
      write (ncontr) nw,nh
      write (ncontr) rgrid,zgrid
      write (ncontr) ggridfc
      write (ncontr) gridpc
!vas just for testing
!      open(35,file='test-ec1.dat',status='new')
!      write (35,*) nw,nh
!      write (35,1009) rgrid,zgrid
!      close(35)
!      open(35,file='test-ec2.dat',status='new')
!      write (35,1009) ggridfc
!      close(35)
!      open(35,file='test-ec3.dat',status='new')
!      write (35,1009) gridpc
!      close(35)
!1009  format(3(1x,e14.8))
      close(unit=ncontr)
!
      if (allocated(gridfc)) deallocate(gridfc)
      if (allocated(ggridfc)) deallocate(ggridfc)
      if (allocated(gridpc)) deallocate(gridpc)
!
      return
      end subroutine efund_grid
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          gsilop computes the green's FUNCTION at si loops due    **
!**          to filament currents flowing in r(n) and z(n).          **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**       rr..............r coordinates                              **
!**       zz..............z coordinates                              **
!**       nr..............DIMENSION of r                             **
!**       rspfun..........computed green's functions values          **
!**       nz..............DIMENSION of z                             **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine gsilop(rr, nr, zz, nz, rspfun, ns, rsi, zsi, wsi &
           , hsi, as, as2, ndim)
      use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,&
                      nw,nh,nwnh
      use consta
      implicit integer*4 (i-n), real*8 (a-h, o-z)
            real*8,dimension(nr) :: rr
            real*8,dimension(nz) :: zz
            real*8,dimension(ns) :: rsi,zsi,wsi,hsi,as,as2
            real*8,dimension(ndim,nwnh) :: rspfun
!     dimension rsi(1),zsi(1),wsi(1),hsi(1),as(1),as2(1)
!     dimension rr(1),zz(1),rspfun(ndim,1)
      dimension taf(nfcoil),taf2(nfcoil)
      data isplit/17/,tole/1.0e-10/

      ndr = isplit
      ndz = isplit
      in=ns
        taf(in) = tan(as(in)*pi/180.0)
        taf2(in) = tan(as2(in)*pi/180.0)
        do nc=1,nr
          do nd=1,nz
            nk=(nc-1)*nz+nd
            rsum = 0.
            dwc = wsi(in)/ndr
            dhc = hsi(in)/ndz
            if (as2(in) .eq. 0.) then
              z1 = zsi(in) - taf(in)*(wsi(in)-dwc)/2. - .5*hsi(in) + .5*dhc
              ab = rsi(in) - .5*wsi(in) + .5*dwc
              do iw = 1, isplit
                drc = ab + (iw-1)*dwc + iw*tole
                z2 = z1 + (iw-1)*taf(in)*dwc
                do ih = 1, isplit
                  dzc = z2 + (ih-1)*dhc
                  rtmp = psical(drc,rr(nc),zz(nd)-dzc)
                  rsum = rsum + rtmp
                enddo
              enddo
            else
              do ih = 1, ndz
                dzc = zsi(in) - .5*hsi(in) + .5*dhc + dhc*(ih-1)
                do iw = 1, ndr
                  drc = rsi(in) - .5*wsi(in) - .5*hsi(in)/taf2(in)&
                         + .5*dwc + .5*dhc/taf2(in)&
                         + dhc/taf2(in)*(ih-1) + dwc*(iw-1)
                  rtmp = psical(drc,rr(nc),zz(nd)-dzc)
                  rsum = rsum + rtmp
                enddo
              enddo
            endif
            cmut = rsum*2.e-07/(isplit*isplit)
            rspfun(in,nk)=cmut
          enddo
        enddo

      return
      end subroutine gsilop
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          gvesel gets the Green's functions for the vessel        **
!**          segments.                                               **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/11/85..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine gvesel(rsilvs,rmp2vs,gridvs,rgrid,mw, &
                        zgrid,mh,rfcvs,rvsfc,rvsec)
      use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,nesum,&
                      nfsum,nvsum,nvesel,nacoil,nw,nh,nwnh
      use consta
      use fcoil
      use cecoil
      use cvesel
      use input,only:iecoil
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension rsilvs(nsilop,nvesel),rmp2vs(magpr2,nvesel), &
                rgrid(1),zgrid(1),rvsec(nvesel,nesum), &
                rfcvs(nfcoil,nvesel),rvsfc(nvesel,nfcoil)
      dimension gridvs(mw*mh,nvesel)
      dimension taf(nfcoil),taf2(nfcoil)
      dimension tas(nvesel),tas2(nvesel)
      do j=1,nsilop
         jj=j
         do i=1,nvesel
            ii=i
            call v1coef(work,jj,ii)
            rsilvs(j,i)=work
         enddo 
      enddo 
!
      do j=1,magpr2
         jj=j
         do i=1,nvesel
            ii=i
            call v2coef(work,jj,ii)
            rmp2vs(j,i)=work
         enddo 
      enddo 
!
      do i=1,nw
         nr=i
         do j=1,nh
            nz=j
            kk=(i-1)*nh+j
            do n=1,nvesel
               nn=n
               call vgrid(work,rgrid,nr,zgrid,nz,nn)
               gridvs(kk,n)=work
            enddo 
         enddo 
      enddo
!
      do j=1,nfcoil
         taf(j)=tan(af(j)*pi/180.)
         taf2(j)=tan(af2(j)*pi/180.)
         do i=1,nvesel
            tas(i)=tan(avs(i)*pi/180.)
            tas2(i)=tan(avs2(i)*pi/180.)
            call flux(rvs(i),zvs(i),wvs(i),hvs(i),tas(i),tas2(i), &
                      rf(j),zf(j),wf(j),hf(j),taf(j),taf2(j),work)
            work=work*0.5/pi
            rfcvs(j,i)=work
         enddo 
      enddo 
!
      aaa=0.0
      rvsec(:,:)=0.0
      do j=1,nvesel
         tas(j)=tan(avs(j)*pi/180.)
         tas2(j)=tan(avs2(j)*pi/180.)
         if (iecoil.eq.1) then
            do i=1,necoil
               call flux(re(i),ze(i),we(i),he(i),aaa,aaa, &
                         rvs(j),zvs(j),wvs(j),hvs(j),tas(j),tas2(j),work)
               work=work*0.5/pi
               kkm=ecid(i)
               kk=kkm+1
               rvsec(j,kkm)=rvsec(j,kkm)+work*(kk-ecid(i))
            enddo
         endif
      enddo 
!
      do j=1,nvesel
         tas(j)=tan(avs(j)*pi/180.)
         tas2(j)=tan(avs2(j)*pi/180.)
         do i=1,nfcoil
            taf(i)=tan(af(i)*pi/180.)
            taf2(i)=tan(af2(i)*pi/180.)
            call flux(rf(i),zf(i),wf(i),hf(i),taf(i),taf2(i), &
                      rvs(j),zvs(j),wvs(j),hvs(j),tas(j),tas2(j),work)
            work=work*0.5/pi
            rvsfc(j,i)=work
         enddo 
      enddo 
      return
      end subroutine gvesel
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          lgauss computes the zeroes of the legendre polynomial   **
!**          and their associated weights for a gaussian quadrature. **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**       x...............zeroes of legendre polynomial              **
!**       w...............weights                                    **
!**       n...............order of legendre polynomial               **
!**       nn..............error flag                                 **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine lgauss(x,w,n,nn)
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      real*8,dimension(n) ::  x,w
!
      nn = 0
      if (n-1.lt.0) then 
         nn = 1
         return
      elseif (n-1.eq.0) then
!----------------------------------------------------------------------
!--      request for a zero point formula is meaningless             --
!----------------------------------------------------------------------
         x(1) = 0.
         w(1) = 2.
         return
      endif
!----------------------------------------------------------------------
!--   for a one point formula, send back                             --
!--   results without computing.                                     --
!----------------------------------------------------------------------
      r = n
      g = -1.
!----------------------------------------------------------------------
!--   the initial guess for the smallest root                        --
!--   of p(n) is taken as -1.                                        --
!----------------------------------------------------------------------
      do i = 1,n
         test = -2.
         ic = n+1-i
!----------------------------------------------------------------------
!--      whenever we find a root of the                              --
!--      polynomial, its negative is also a root.                    --
!--      the index ic tells where to store the other root            --
!----------------------------------------------------------------------
         if (ic.lt.i) return
   40    s = g
         t = 1.
         u = 1.
         v = 0.
!----------------------------------------------------------------------
!--      evaluation of the n-th legendre polynomial                  --
!--      and its first derivative.                                   --
!--      where   u = ds/dx                                           --
!--              v = dt/dx                                           --
!--              dp=dp/dx                                            --
!----------------------------------------------------------------------
         do k = 2,n
            a = k
            p = ((2.0*a-1.0)*s*g-(a-1.0)*t)/a
             dp = ((2.0*a-1.0)*(s+g*u)-(a-1.0)*v)/a
            v = u
            u = dp
            t = s
            s = p
         enddo 
         if (abs((test-g)/(test+g)).ge.0.0000005) then
            sum = 0.
            if (i.ne.1) then
!----------------------------------------------------------------------
!--            the following computes the reduced                    --
!--            legendre polynomial and its derivative.               --
!----------------------------------------------------------------------
               do k = 2,i
                  sum = sum+1./(g-x(k-1))
               enddo
            endif
            test = g
            g = g-p/(dp-p*sum)
            go to 40
         endif
         x(ic) = -g
         x(i) = g
         w(i) = 2./(r*t*dp)
         w(ic) = w(i)
         g = g-r*t/((r+2.)*g*dp+r*v-2.*r*t*sum)
      enddo
      return
      end subroutine lgauss
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          m1coef computes the response functions due to           **
!**          the thin flux loops.                                    **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          15/07/83..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine m1coef(rr, zz, nr, nz, coef,  nl, nf)
      use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil
      use fcoil
      use coilsp
      use consta
      use nio
      use siloop
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      real*8,dimension(nr) :: rr
      real*8,dimension(nz) :: zz
      real*8,dimension(nsilop,nr*nz) :: coef
!
      radeg=pi/180.
      isplit=17
      itot=isplit*isplit
      fitot=itot
!
      m=nl
      if (nf.le.0) then
         do ii=1,nr
            do jj=1,nz
               kk=(ii-1)*nz+jj
               a=rr(ii)
               r=rsi(m)
               z=zsi(m)-zz(jj)
               cmp2=psical(a,r,z)*tmu
               coef(m,kk)=cmp2
            enddo 
         enddo 
      else
         k=nf
         psict=0
         call splitc(isplit,rsplt,zsplt, &
                     rf(k),zf(k),wf(k),hf(k),af(k),af2(k))
         do l=1,itot
            a=rsplt(l)
            r1=rsi(m)
            z1=zsi(m)-zsplt(l)
            psic=psical(a,r1,z1)*tmu
            psict=psict+psic/fitot
         enddo 
         coef(m,k)=psict
      endif
!
      return
      end subroutine m1coef
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          m2coef computes the response functions due to           **
!**          the magnetic probes.                                    **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine m2coef(rr, nr, zz, nz, coef,  mp, nc)
      use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,&
                      nw,nh,nwnh
      use fcoil
      use coilsp
      use mprobe
      use pmodel
      use consta
      use fshift
      use bfgrid
      implicit integer*4 (i-n), real*8 (a-h, o-z)
            real*8,dimension(nr) :: rr
            real*8,dimension(nz) :: zz
      dimension coef(mp,nc)

      radeg=pi/180.
      isplit=17
      itot=isplit*isplit
      fitot=itot
!
      do m=1,magpr2
         if (smp2(m).gt.0.0) then
            sinm=sin(radeg*amp2(m))
            cosm=cos(radeg*amp2(m))
            delsx=smp2(m)/nsmp2*cosm
            delsy=smp2(m)/nsmp2*sinm
         else
!------------------------------------------------------------------------------
!--      perpendicular probes                                                --
!------------------------------------------------------------------------------
            sinm=sin(radeg*amp2(m))
            cosm=cos(radeg*amp2(m))
            sinms=sin(radeg*(amp2(m)+90.))
            cosms=cos(radeg*(amp2(m)+90.))
            delsx=abs(smp2(m))/nsmp2*cosms
            delsy=abs(smp2(m))/nsmp2*sinms
         endif
         xmp20=xmp2(m)-(nsmp2-1)/2.*delsx
         ymp20=ymp2(m)-(nsmp2-1)/2.*delsy
         if (nz.gt.0) then
            do ii=1,nr
               do jj=1,nz
                  kk=(ii-1)*nz+jj
                  a=rr(ii)
                  brct=0.0
                  bzct=0.0
                  do mmm=1,nsmp2
                     r=xmp20+(mmm-1)*delsx
                     z=ymp20+(mmm-1)*delsy-zz(jj)
                     brtmp=br(a,r,z)*tmu
                     bztmp=bz(a,r,z)*tmu
                     brct=brtmp+brct
                     bzct=bztmp+bzct
                  enddo 
                  cmp2=(brct*cosm+bzct*sinm)/nsmp2
                  coef(m,kk)=cmp2
               enddo
            enddo
         else
            do k=1,nfcoil
!---------------------------------------------------------------
!--         shifted f-coil                                    --
!---------------------------------------------------------------
               if (k.eq.nshiftrz(k)) then
                  pmnow=radeg*pmprobe(m)
                  pfnow=radeg*pshift(k)
               endif
!
               brct=0
               bzct=0
               call splitc(isplit,rsplt,zsplt, &
                           rf(k),zf(k),wf(k),hf(k),af(k),af2(k))
               do l=1,itot
                  a=rsplt(l)
                  do mmm=1,nsmp2
                     r1=xmp20+(mmm-1)*delsx
                     z1=ymp20+(mmm-1)*delsy-zsplt(l)
!---------------------------------------------------------------
!--                  shifted f-coil                           --
!---------------------------------------------------------------
                     if (k.eq.nshiftrz(k)) then
                        rcos=r1*cos(pmnow)-rshift(k)*cos(pfnow)
                        rsin=r1*sin(pmnow)-rshift(k)*sin(pfnow)
                        r1=sqrt(rcos**2+rsin**2)
                        z1=z1-zshift(k)
                     endif
!
                     brc=br(a,r1,z1)*tmu
                     bzc=bz(a,r1,z1)*tmu
                     bzct=bzct+bzc/fitot
!---------------------------------------------------------------
!--                  shifted f-coil ?                         --
!---------------------------------------------------------------
                     if (k.eq.nshiftrz(k)) then
                        cospp=rcos/r1
                        sinpp=rsin/r1
                        cfactor=cos(pmnow)*cospp+sin(pmnow)*sinpp
                        brct=brct+brc/fitot*cfactor
                     else
                        brct=brct+brc/fitot
                     endif
                  enddo
               enddo
!
               coef(m,k)=(brct*cosm+bzct*sinm)/nsmp2
            enddo
         endif
      enddo
!----------------------------------------------------------------
!--   br, bz at grid due to f-coils                            --
!----------------------------------------------------------------
      if (nz.gt.0) then
         return
      else
         do ii=1,nw
            do jj=1,nh
               kk=(ii-1)*nh+jj
               do k=1,nfcoil
                  brct=0
                  bzct=0
                  call splitc(isplit,rsplt,zsplt, &
                              rf(k),zf(k),wf(k),hf(k),af(k),af2(k))
                  do l=1,itot
                     a=rsplt(l)
                     r1=rgrid(ii)
                     z1=zgrid(jj)-zsplt(l)
                     brc=br(a,r1,z1)*tmu
                     bzc=bz(a,r1,z1)*tmu
                     bzct=bzct+bzc/fitot
                     brct=brct+brc/fitot
                  enddo 
                  brgridfc(kk,k)=brct
                  bzgridfc(kk,k)=bzct
               enddo
            enddo
         enddo
      endif
!
      return
      end subroutine m2coef
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          matrix calculates the appropriate response functions.   **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine efund_matrix
      use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,nesum,&
                      nfsum,nvsum,nvesel,nacoil,nw,nh,nwnh
      use consta
      use nio
      use cvesel
      use input
      use cacoil
      use pmodel
      use siloop
      use fcoil
      use fshift
      use bfgrid
!vas
      use var_filech
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension rfcfc(nfcoil,nfcoil)
      dimension rsilfc(nsilop,nfcoil),rmp2fc(magpr2,nfcoil), &
                rgowfc(nrogow,nfcoil)
      dimension gsilfc(nsilop,nfsum),gmp2fc(magpr2,nfsum)
      dimension rsilec(nsilop,nesum),rmp2ec(magpr2,nesum), &
                rfcec(nfcoil,nesum), &
                recec(nesum,nesum),rsisec(nesum)
      dimension rsilvs(nsilop,nvesel),rmp2vs(magpr2,nvesel), &
                rfcvs(nfcoil,nvesel), &
                rvsec(nvesel,nesum),rvsfc(nvesel,nfcoil), &
                rvsvs(nvesel,nvesel),tas(nvesel),tas2(nvesel)
      dimension gsilvs(nsilop,nvsum),gmp2vs(magpr2,nvsum)
      dimension taf(nfcoil),taf2(nfcoil)
      dimension rsilac(nsilop,nacoil),rmp2ac(magpr2,nacoil)
      dimension xdum(1),ydum(1)
      real*8,dimension(:,:),allocatable :: rfcpc,brgrfc,bzgrfc, &
                                         rsilpc,rmp2pc,rgowpc, &
                                         gridec,gridvs,ggridvs, &
                                         gridac
!
      if (.not. allocated(rfcpc)) then
        allocate(rfcpc(nfcoil,nwnh))
        rfcpc(:,:) = 0.0
      endif
      if (.not. allocated(brgrfc)) then
        allocate(brgrfc(nwnh,nfsum))
        brgrfc(:,:) = 0.0
      endif
      if (.not. allocated(bzgrfc)) then
        allocate(bzgrfc(nwnh,nfsum))
        bzgrfc(:,:) = 0.0
      endif
      if (.not. allocated(rsilpc)) then
        allocate(rsilpc(nsilop,nwnh))
        rsilpc(:,:) = 0.0
      endif
      if (.not. allocated(rmp2pc)) then
        allocate(rmp2pc(magpr2,nwnh))
        rmp2pc(:,:) = 0.0
      endif
      if (.not. allocated(rgowpc)) then
        allocate(rgowpc(nrogow,nwnh))
        rgowpc(:,:) = 0.0
      endif
      if (.not. allocated(gridec)) then
        allocate(gridec(nwnh,nesum))
        gridec(:,:) = 0.0
      endif
      if (.not. allocated(gridvs)) then
        allocate(gridvs(nwnh,nvesel))
        gridvs(:,:) = 0.0
      endif
      if (.not. allocated(ggridvs)) then
        allocate(ggridvs(nwnh,nvsum))
        ggridvs(:,:) = 0.0
      endif
      if (.not. allocated(gridac)) then
        allocate(gridac(nwnh,nacoil))
        gridac(:,:) = 0.0
      endif
!
      if (ifcoil.eq.1) then
!----------------------------------------------------------------------
!--      calculate the response function of psi loops due to f coils --
!----------------------------------------------------------------------
         do i=1,nfcoil
            taf(i)=tan(af(i)*pi/180.)
            taf2(i)=tan(af2(i)*pi/180.)
         enddo 
         if (islpfc.eq.1) then
            do i=1,nfcoil
               do j=1,nfcoil
                  call flux(rf(i),zf(i),wf(i),hf(i),taf(i),taf2(i), &
                            rf(j),zf(j),wf(j),hf(j),taf(j),taf2(j), &
                            rfcfc(j,i))
                  rfcfc(j,i)=rfcfc(j,i)*0.5/pi

               enddo
               ii=i
               call gsilop(rgrid,nw,zgrid,nh,rfcpc,ii,rf,zf,wf,hf,af,af2 &
                          ,nfcoil)
            enddo

            print*,'file name : ','fc'//trim(ch1)//trim(ch2)//'.ddd' 
!vasorg      open(unit=nrspfc,status='unknown',file='fcfcpc.dat', &
            open(unit=nrspfc,status='unknown',file='fc'//trim(ch1)// & 
                 trim(ch2)//'.ddd' , &
                 form='unformatted')
            write (nrspfc) rfcfc
            write (nrspfc) rfcpc
            close(unit=nrspfc)
         endif
!---------------------------------------------------------------------
!--      flux loops                                                 --
!---------------------------------------------------------------------
         if (nsilop.gt.1) then
!---------------------------------------------------------------------
!           finite size flux loops                                  --
!---------------------------------------------------------------------
            if (isize.gt.0) then
               do i=1,nfcoil
                  do j=1,isize
                     taz=tan(as(j)*pi/180.)
                     taz2=tan(as2(j)*pi/180.)
                     call flux(rf(i),zf(i),wf(i),hf(i),taf(i),taf2(i), &
                               rsi(j),zsi(j),wsi(j),hsi(j),taz,taz2, &
                               rsilfc(j,i))
                     rsilfc(j,i)=rsilfc(j,i)*0.5/pi
                  enddo 
               enddo
            endif
!---------------------------------------------------------------------
!           thin flux loops                                         --
!---------------------------------------------------------------------
            if (isize.lt.nsilop) then
               do i=1,nfcoil
                  ii=i
                  do j=isize+1,nsilop
                     jj=j
                     call m1coef(xdum,xdum,nsilop,nfcoil,rsilfc,jj,ii)
                  enddo 
               enddo 
            endif 
         endif 

         !
         if (.not. allocated(brgridfc)) then
            allocate(brgridfc(nwnh,nfcoil))
            brgridfc(:,:) = 0.0
         endif
         if (.not. allocated(bzgridfc)) then
            allocate(bzgridfc(nwnh,nfcoil))
            bzgridfc(:,:) = 0.0
         endif
!
!-----------------------------------------------------------------------
!--      compute the response function of magnetic probes due to f coils
!-----------------------------------------------------------------------
         magprr=magpr2
         if (magprr.gt.1) then
            call m2coef(xdum,0,ydum,0,rmp2fc,magpr2,nfcoil)
         endif
!----------------------------------------------------------------------
!--      compute the response function of partial rogowski loops due to
!--      f coils
!----------------------------------------------------------------------
         mrogow=nrogow
         if (mrogow.gt.1) then
            call rogowc(xdum,0,ydum,0,rgowfc,nrogow,nfcoil)
         endif 
!----------------------------------------------------------------------
!--      write f coil response functions                             --
!----------------------------------------------------------------------
         do i=1,nfsum
            do j=1,nsilop
               gsilfc(j,i)=0.0
            enddo
            do j=1,magpr2
               gmp2fc(j,i)=0.0
            enddo
         enddo
!
         do i=1,nfcoil
            k=abs(fcid(i))
            do j=1,nsilop
               gsilfc(j,k)=gsilfc(j,k)+fcturn(i)*rsilfc(j,i)
            enddo
            do j=1,magpr2
               gmp2fc(j,k)=gmp2fc(j,k)+fcturn(i)*rmp2fc(j,i)
            enddo
         enddo
!
         print*,'file name : ','rfcoil.ddd' 
!vasorg      open(unit=nrspfc,status='unknown',file='rfcoil.dat', &
         open(unit=nrspfc,status='unknown',file='rfcoil.ddd', &
              form='unformatted')
         write (nrspfc) gsilfc
         write (nrspfc) gmp2fc
         close(unit=nrspfc)
!
         do i=1,nfsum
            do j=1,nwnh
               brgrfc(j,i)=0.0
               bzgrfc(j,i)=0.0
            enddo
         enddo
         do i=1,nfcoil
            k=abs(fcid(i))
            do j=1,nwnh
               brgrfc(j,k)=brgrfc(j,k)+fcturn(i)*brgridfc(j,i)
               bzgrfc(j,k)=bzgrfc(j,k)+fcturn(i)*bzgridfc(j,i)
            enddo
         enddo
!
         open(unit=nrspfc,status='unknown',file='brzgfc.dat', &
              form='unformatted')
         write (nrspfc) brgrfc
         write (nrspfc) bzgrfc
         close(unit=nrspfc)
      endif
!----------------------------------------------------------------------
!--   plasma response functions                                      --
!----------------------------------------------------------------------
      if (igrid.eq.1) then
         msilop=nsilop
         if (msilop.gt.1) then
!----------------------------------------------------------------------
!--         filament plasma current model                            --
!----------------------------------------------------------------------
            if (isize.gt.0) then
               do j=1,isize
                  jj=j
                  call gsilop(rgrid,nw,zgrid,nh,rsilpc,jj, &
                              rsi,zsi,wsi,hsi,as,as2,nsilop)
               enddo
            endif
         endif
         if (isize.lt.nsilop) then
            do j=isize+1,nsilop
               jj=j
               call m1coef(rgrid,zgrid,nw,nh,rsilpc,jj,0)
            enddo 
         endif
         magprr=magpr2
         if (magprr.gt.1) then
            call m2coef(rgrid,nw,zgrid,nh,rmp2pc,magpr2,nwnh)
         endif
         mrogow=nrogow
         if (mrogow.gt.1) then
            call rogowc(rgrid,nw,zgrid,nh,rgowpc,nrogow,nwnh)
         endif
!----------------------------------------------------------------------
!--      write the plasma response function                          --
!----------------------------------------------------------------------
         print*,'file name : ','ep'//trim(ch1)//trim(ch2)//'.ddd' 
!vasorg      open(unit=nrsppc,status='unknown',file='eplasm.dat', &
         open(unit=nrsppc,status='unknown',file='ep'//trim(ch1)// & 
              trim(ch2)//'.ddd' , &
              form='unformatted')
         write (nrsppc) rsilpc
         write (nrsppc) rmp2pc
         close(unit=nrsppc)
!
      endif
      if (iecoil.eq.1) then
         call gecoil(rsilec,rmp2ec,gridec,rgrid,nw, &
                     zgrid,nh,rfcec,recec,rsisec)
         print*,'file name : ','re'//trim(ch1)//trim(ch2)//'.ddd' 
!vasorg      open(unit=nrsppc,status='unknown',file='recoil.dat', &
         open(unit=nrsppc,status='unknown',file='re'//trim(ch1)// & 
                     trim(ch2)//'.ddd', &
              form='unformatted')
         write (nrsppc) rsilec
         write (nrsppc) rmp2ec
         write (nrsppc) gridec
         close(unit=nrsppc)
      endif
!
      if (ivesel.eq.1) then
         call gvesel(rsilvs,rmp2vs,gridvs,rgrid,nw, &
                     zgrid,nh,rfcvs,rvsfc,rvsec)
         do i=1,nvesel
           tas(i)=tan(avs(i)*pi/180.)
           tas2(i)=tan(avs2(i)*pi/180.)
         enddo 
         do i=1,nvesel
            do j=1,nvesel
               call flux(rvs(i),zvs(i),wvs(i),hvs(i),tas(i),tas2(i), &
                         rvs(j),zvs(j),wvs(j),hvs(j),tas(j),tas2(j), &
                         rvsvs(j,i))
               rvsvs(j,i)=rvsvs(j,i)*0.5/pi
            enddo 
         enddo 
!
         do i=1,nvsum
            do j=1,nsilop
               gsilvs(j,i)=0.0
            enddo
            do j=1,magpr2
               gmp2vs(j,i)=0.0
            enddo
            do j=1,nwnh
               ggridvs(j,i)=0.0
            enddo
         enddo
!
         do i=1,nvesel
            k=abs(vsid(i))
            do j=1,nsilop
               gsilvs(j,k)=gsilvs(j,k)+rsilvs(j,i)
            enddo
            do j=1,magpr2
               gmp2vs(j,k)=gmp2vs(j,k)+rmp2vs(j,i)
            enddo
            do j=1,nwnh
               ggridvs(j,k)=ggridvs(j,k)+gridvs(j,i)
            enddo
         enddo
!
!vas
      print*,'file name : ','rv'//trim(ch1)//trim(ch2)//'.ddd' 
!vasorg         open(unit=nrsppc,status='unknown',file='rvesel.dat', &
         open(unit=nrsppc,status='unknown',file='rv'//trim(ch1)// & 
              trim(ch2)//'.ddd' , &
              form='unformatted')
         write (nrsppc) gsilvs
         write (nrsppc) gmp2vs
         write (nrsppc) ggridvs
         close(unit=nrsppc)
      endif
!---------------------------------------------------------------------
!--   advance divertor coil                                         --
!---------------------------------------------------------------------
      if (iacoil.eq.1) then
         call gacoil(rsilac,rmp2ac,gridac,rgrid,nw, &
                     zgrid,nh)
!vas
         print*,'file name : ','ra'//trim(ch1)//trim(ch2)//'.ddd' 
!vasorg         open(unit=nrsppc,status='unknown',file='racoil.dat', &
         open(unit=nrsppc,status='unknown',file='ra'//trim(ch1)// & 
              trim(ch2)//'.ddd' , &
              form='unformatted')
         write (nrsppc) gridac
         write (nrsppc) rsilac
         write (nrsppc) rmp2ac
         close(unit=nrsppc)
      endif
!
      if (allocated(rfcpc)) deallocate(rfcpc)
      if (allocated(brgrfc)) deallocate(brgrfc)
      if (allocated(bzgrfc)) deallocate(bzgrfc)
      if (allocated(rsilpc)) deallocate(rsilpc)
      if (allocated(rmp2pc)) deallocate(rmp2pc)
      if (allocated(rgowpc)) deallocate(rgowpc)
      if (allocated(gridec)) deallocate(gridec)
      if (allocated(gridvs)) deallocate(gridvs)
      if (allocated(ggridvs)) deallocate(ggridvs)
      if (allocated(gridac)) deallocate(gridac)
!
      return
      end subroutine efund_matrix
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          psical computes mutual inductance/2/pi between two      **
!**          circular filaments of radii a1 and r1 and               **
!**          separation of z1, for mks units multiply returned       **
!**          value by 2.0e-07.                                       **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**       a1..............first filament radius                      **
!**       r1..............second filament radius                     **
!**       z1..............vertical separation                        **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1) f.w. mcclain and b.b. brown, ga technologies        **
!**              report ga-a14490 (1977).                            **
!**                                                                  **
!**********************************************************************
      real*8 function psical(a,r,z)
      implicit none
      real*8 br,bz,xmdelk,xmdele
      real*8, intent(in) :: a,r,z
      integer*4 isw
      real*8 den,xk,x1,cay,ee
!
      isw=1
      go to 10
      entry br(a,r,z)
      isw=2
      go to 10
      entry bz(a,r,z)
      isw=3
!
   10 continue
      den=a*a+r*r+z*z+2.*a*r
      xk=4.*a*r/den
      x1=(a*a+r*r+z*z-2.*a*r)/den
      if(x1.lt.1.0e-10) x1=1.0e-10
      cay=xmdelk(x1)
      ee=xmdele(x1)
      select case (isw)
      case (1)
!----------------------------------------------------------------------
!--      psi computation                                             --
!----------------------------------------------------------------------
         psical=sqrt(den)*((1.e+00-0.5e+00*xk)*cay-ee)
         return
      case (2)
!----------------------------------------------------------------------
!--      br  computation                                             --
!----------------------------------------------------------------------
         br=z/(r*sqrt(den))*(-cay+(a*a+r*r+z*z)/((a-r)*(a-r)+z*z)*ee)
         return
      case (3)
!----------------------------------------------------------------------
!--      bz  computation                                             --
!----------------------------------------------------------------------
         bz=(cay+(a*a-r*r-z*z)/((a-r)*(a-r)+z*z)*ee)/sqrt(den)
         return
      end select
      end function psical
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine rogowc(rr, nrr, zz, nz, coef, nr, nc)
      use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil
      use rogowl
      use coilsp
      use consta
      use fcoil
      implicit integer*4 (i-n), real*8 (a-h, o-z)
!      dimension rogpth(nrogow)
      real*8,dimension(nr) :: rr
      real*8,dimension(nz) :: zz
      dimension coef(nr,nc)
!
      ngrid=25
      isplit=17
      itot=isplit*isplit
      fitot=itot
      dels = 0.
      mm=1
!
      do m=1,nrogow
         if (nz.gt.0) then
            do inn=1,nrr
               do imm=1,nz
                  ikk=(inn-1)*nz+imm
                  coef(m,ikk)=0.0
               enddo 
            enddo 
         else
            do k=1,nfcoil
               coef(m,k)=0.0
            enddo 
         endif
         k=m
         call rogrid(ngrid,mm,k,dels)
         mm=mm+narc(m)+1
         do i=1,ngrid
            iii=i
            if(i.eq.ngrid) then
               zl=zpg(i)-zpg(i-1)
               rl=rpg(i)-rpg(i-1)
            else
               zl=zpg(i+1)-zpg(i)
               rl=rpg(i+1)-rpg(i)
            endif
            hl=sqrt(zl*zl+rl*rl)
            sint=zl/hl
            cost=rl/hl
!
            if (nz.le.0) then
               do k=1,nfcoil
                  call splitc(isplit,rsplt,zsplt, &
                           rf(k),zf(k),wf(k),hf(k),af(k),af2(k))
                  do l=1,itot
                     a=rsplt(l)
                     r1=rpg(i)
                     z1=zpg(i)-zsplt(l)
                     brc=br(a,r1,z1)*tmu/fitot
                     bzc=bz(a,r1,z1)*tmu/fitot
                     part=brc*cost+bzc*sint
                     call simpf(iii,fact)
                     ! todo: rogpth is never defined in efit...
                     coef(m,k)=coef(m,k)+fact*part*dels !/rogpth(m)
                  enddo 
               enddo
            else
               do inn=1,nrr
                  do imm=1,nz
                     ikk=(inn-1)*nz+imm
                     a=rr(inn)
                     r1=rpg(i)
                     z1=zpg(i)-zz(imm)
                     brg=br(a,r1,z1)*tmu
                     bzg=bz(a,r1,z1)*tmu
                     part=brg*cost+bzg*sint
                     call simpf(iii,fact)
                     ! todo: rogpth is never defined in efit...
                     coef(m,ikk)=coef(m,ikk)+fact*part*dels !/rogpth(m)
                  enddo 
               enddo 
            endif
         enddo
      enddo
!
      return
      end subroutine rogowc
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          rogrid calculates grid points along the given arc       **
!**          made up of at most six straight line segments.          **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**       ngrid...........                                           **
!**       mm..............                                           **
!**       m...............                                           **
!**       dels............                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine rogrid(ngrid,mm,m,dels)
      use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil
      use rogowl
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension sl(6)
!
      s = 0.
      mm1 = mm+narc(m)-1
      j = 1
      do i = mm,mm1
         sl(j) = sqrt((rp(i+1)-rp(i))**2 + (zp(i+1)-zp(i))**2)
         s = s+sl(j)
         j = j+1
      enddo 
      dels = s/dble(ngrid-1)
      ds = 0.
      i1 = 1
      do j = 1,narc(m)
         dr = rp(mm+j)-rp(mm+j-1)
         dz = zp(mm+j)-zp(mm+j-1)
         rpg(i1) = rp(mm+j-1)+dr*ds/sl(j)
         zpg(i1) = zp(mm+j-1)+dz*ds/sl(j)
         dd = sl(j)-ds
         n1 = int(dd/dels)+i1
         i2 = i1+1
         dr = dr*dels/sl(j)
         dz = dz*dels/sl(j)
         do i = i2,n1
            rpg(i) = rpg(i-1)+dr
            zpg(i) = zpg(i-1)+dz
         enddo 
         ds = dels-(dd-dble(n1-i1)*dels)
         i1 = n1+1
      enddo
      rpg(ngrid) = rp(mm+narc(m))
      zpg(ngrid) = zp(mm+narc(m))
      return
      end subroutine rogrid
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine simpf(i,f)
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      if (i.eq.1 .or. i.eq.25) then
         f = 1./3.
      elseif ((i/2)*2.eq.i) then
         f = 4./3.
      else
         f = 2./3.
      endif
      return
      end subroutine simpf
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          soleno computes the inductance for a solenoid           **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine soleno(ra,z1,w1,h1,t1,t12,r2,z2,w2,h2,t2,t22,xbm1,xbm2, &
                        xtm1,xtm2,hfa,hsa,rf,rs,sol)
      use consta,only:pi
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension z(2,2)
      real*8 ksq,kpsq,msl,mut
!      data init/0/
!
      rpi=pi
      rpi2=rpi*0.5
      rh=rpi*2.0e-07
      ut=rh/3.0
      err=1.0e-05
!      init=1
!
      dr = rf-rs
      drsq = dr*dr
      tr = rf+rs
      trsq = tr*tr
      fr = 4.*rf*rs
      sol = 0.
      csq = fr/trsq
      cpsq = 1.-csq
!
      if (t12.eq.0.) then
         r = rf-ra+0.5*w1
         zc1 = z1-0.5*h1-0.5*w1*t1
         zb1 = zc1+t1*r
         zt1 = zb1+h1
!
      else
         zb1 = z1-h1/2.
         if (t12.lt.0.) then
            if (rf .lt. xbm1) zb1 = zb1 + t12*(rf-xbm1)
            zt1 = z1 + h1/2.
            if (rf .gt. xtm1) zt1 = zt1 - t12*(xtm1-rf)
!
         else
            if (rf.gt.xbm1) zb1 = zb1+t12*(rf-xbm1)
            zt1 = z1+h1/2.
            if (rf.lt.xtm1) zt1 = zt1-t12*(xtm1-rf)
         endif
      endif
!
      if (t22.eq.0.) then
         r = rs-r2+0.5*w2
         zc2 = z2-0.5*h2-0.5*w2*t2
         zb2 = zc2+t2*r
         zt2 = zb2+h2
!
      else
         zb2 = z2-h2/2.
         if (t22 .lt. 0.) then
            if (rs .lt. xbm2) zb2 = zb2 + t22*(rs-xbm2)
            zt2 = z2 + h2/2.
            if (rs .gt. xtm2) zt2 = zt2 - t22*(xtm2-rs)
!
         else
            if (rs.gt.xbm2) zb2 = zb2+t22*(rs-xbm2)
            zt2 = z2+h2/2.
            if (rs.lt.xtm2) zt2 = zt2-t22*(xtm2-rs)
         endif
      endif
!
      z(1,1) = zb1
      z(2,1) = zb2
      z(1,2) = zt1
      z(2,2) = zt2
      hfa = zt1-zb1
      hsa = zt2-zb2
!
      do i = 1,2
         do j = 1,2
            sign = -.25
            if (i .ne. j) sign = .25
            dz = z(1,i)-z(2,j)
            dzsq = dz*dz
            r2sq = dzsq+drsq
            r1sq = dzsq+trsq
            r1 = sqrt(r1sq)
            ksq = fr/r1sq
            t = 2./ksq-1.
            kpsq = 1.-ksq
            alpha = 1.
!--------------------------------------------------------------------------
!--         to avoid numerical truncation                                --
!--------------------------------------------------------------------------
            if (kpsq .lt. 1.0e-30) kpsq = 1.0e-30
            beta = sqrt(kpsq)
            if (beta .lt. 1.0e-30) beta = 1.0e-10
            if (cpsq .lt. 1.0e-30) cpsq = 1.0e-10
            delta = cpsq/beta
            epsi = csq/cpsq
            zeta = 0.
            sinf = 0.
            sa = .25
!
  100       continue
            sa = 2.*sa
            ambsq = (alpha-beta)*(alpha-beta)
            sinf = sinf+sa*ambsq
            alphat = alpha
            epsit = epsi
            alpha = .5*(alpha+beta)
            beta = sqrt(alphat*beta)
            epsi = (delta*epsi+zeta)/(1.+delta)
            delta = beta/4./alpha*(2.+delta+1./delta)
            zeta = .5*(epsit+zeta)
            if (abs(delta-1.) .gt. err) go to 100
            if (ambsq .gt. 1.e-14) go to 100
            cay = rpi2/alpha
            pik = cay*zeta
            ek = .5*cay*(ksq+sinf)
            msl = rh*dzsq*(r1*ek-drsq*pik/r1)
            if (csq==1.) msl = rh*dzsq*(r1*ek-cay*fr/r1*.5)
!
            mut = msl+ut*fr*r1*(cay-t*ek)
            sol = sol+sign*mut
         enddo
      enddo
!
      return
      end subroutine soleno
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine splitc(is,rs,zs,rc,zc,wc,hc,ac,ac2)
      use consta
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      real*8,dimension(is*is) :: rs,zs
!
      frd=pi/180.
!----------------------------------------------------------------------
!--   rectangle                                                      --
!----------------------------------------------------------------------
      if(ac+ac2.eq.0.) then
          wdelt=wc/is
          hdelt=hc/is
          rstrt=rc-wc/2.+wdelt/2.
          zstrt=zc-hc/2.+hdelt/2.
          zz=zstrt
          ic=0
          do ii=1,is
             rr=rstrt
             do jj=1,is
                ic=ic+1
                zs(ic)=zz
                rs(ic)=rr
                rr=rr+wdelt
             enddo 
             zz=zz+hdelt
          enddo
          return
!----------------------------------------------------------------------
!--   ac .ne. 0                                                      --
!----------------------------------------------------------------------
      elseif(ac.ne.0.) then
          side=tan(frd*ac)*wc
          hdelt=hc/is
          wdelt=wc/is
          zdelt=tan(frd*ac)*wdelt
          rstrt=rc-wc/2.+wdelt/2.
          tsid=hc+side
          zstrt =zc-tsid/2.+tsid/2.*1./is
          rr=rstrt
          ic=0
          do ii=1,is
             zz=zstrt+(ii-1)*zdelt
             do jj=1,is
                ic=ic+1
                zs(ic)=zz
                rs(ic)=rr
                zz=zz+hdelt
             enddo 
             rr=rr+wdelt
          enddo
          return
!----------------------------------------------------------------------
!--   ac2 .ne. 0                                                     --
!----------------------------------------------------------------------
      elseif(ac2.ne.0.) then
          side=hc/tan(frd*ac2)
          hdelt=hc/is
          wdelt=wc/is
          zstrt=zc-hc/2.+hdelt/2.
          rdelt=hdelt/tan(frd*ac2)
          rstrt=rc-side/2.-wc/2.+rdelt/2.+wdelt/2.
          side=hc/tan(frd*ac2)
          wtot=side+wc
          whaf=(side+wc)/2.
          rcorn=rc-whaf
          rcornr=rc+whaf
          rcorn2=rcorn+wtot/is
          rstrt=(rcorn+rcorn2)/2.
          zz=zstrt
          ic=0
          do ii=1,is
             rr=rstrt+(ii-1)*rdelt
             do jj=1,is
                ic=ic+1
                zs(ic)=zz
                rs(ic)=rr
                rr=rr+wdelt
             enddo 
             zz=zz+hdelt
          enddo
      endif
!
      return
      end subroutine splitc
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          e1coef computes the response functions due to           **
!**          the thin flux loops and the vessel segments.            **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          15/07/83..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine v1coef(coef,  nl, ne)
      use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,&
                      nesum,nfsum,nvsum,nvesel,nacoil
      use coilsp
      use consta
      use nio
      use siloop
      use cvesel
      implicit integer*4 (i-n), real*8 (a-h, o-z)
!
      radeg=pi/180.
      isplit=17
      itot=isplit*isplit
      fitot=itot
!
      m=nl
      k=ne
      psict=0
      call splitc(isplit,rsplt,zsplt, &
                  rvs(k),zvs(k),wvs(k),hvs(k),avs(k),avs2(k))
      do l=1,itot
         a=rsplt(l)
         r1=rsi(m)
         z1=zsi(m)-zsplt(l)
         psic=psical(a,r1,z1)*tmu
         psict=psict+psic/fitot
      enddo 
      coef=psict
!
      return
      end subroutine v1coef
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          e2coef computes the response functions due to           **
!**          the magnetic probes and the vessel segments.            **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine v2coef(coef, mp, ne)
      use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,&
                      nesum,nfsum,nvsum,nvesel,nacoil
      use coilsp
      use mprobe
      use consta
      use cvesel
      implicit integer*4 (i-n), real*8 (a-h, o-z)
!
      radeg=pi/180.
      isplit=17
      itot=isplit*isplit
      fitot=itot
!
      m=mp
      if (smp2(m).gt.0.0) then
         sinm=sin(radeg*amp2(m))
         cosm=cos(radeg*amp2(m))
         delsx=smp2(m)/nsmp2*cosm
         delsy=smp2(m)/nsmp2*sinm
      else
!------------------------------------------------------------------------------
!--      perpendicular probes                                                --
!------------------------------------------------------------------------------
         sinm=sin(radeg*amp2(m))
         cosm=cos(radeg*amp2(m))
         sinms=sin(radeg*(amp2(m)+90.))
         cosms=cos(radeg*(amp2(m)+90.))
         delsx=abs(smp2(m))/nsmp2*cosms
         delsy=abs(smp2(m))/nsmp2*sinms
      endif
      xmp20=xmp2(m)-(nsmp2-1)/2.*delsx
      ymp20=ymp2(m)-(nsmp2-1)/2.*delsy
      k=ne
      brct=0
      bzct=0
      call splitc(isplit,rsplt,zsplt, &
                  rvs(k),zvs(k),wvs(k),hvs(k),avs(k),avs2(k))
      do l=1,itot
         a=rsplt(l)
         do mmm=1,nsmp2
            r1=xmp20+(mmm-1)*delsx
            z1=ymp20+(mmm-1)*delsy-zsplt(l)
            brc=br(a,r1,z1)*tmu
            bzc=bz(a,r1,z1)*tmu
            brct=brct+brc/fitot
            bzct=bzct+bzc/fitot
         enddo 
      enddo 
      coef=(brct*cosm+bzct*sinm)/nsmp2
!
      return
      end subroutine v2coef
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          egrid computes the response functions due to            **
!**          the grid points and the vessel segments.                **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/11/85..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      subroutine vgrid(coef, rgrid, nr, zgrid, nz, ne)
      use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,&
                      nesum,nfsum,nvsum,nvesel,nacoil
      use coilsp
      use siloop
      use consta
      use nio
      use cvesel
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      real*8,dimension(nr) :: rgrid
      real*8,dimension(nz) :: zgrid
      data init/0/
!
      radeg=pi/180.
      isplit=17
      itot=isplit*isplit
      fitot=itot
!
      k=ne
      psict=0
      call splitc(isplit,rsplt,zsplt, &
                  rvs(k),zvs(k),wvs(k),hvs(k),avs(k),avs2(k))
      do l=1,itot
         a=rsplt(l)
         r1=rgrid(nr)
         z1=zgrid(nz)-zsplt(l)
         psic=psical(a,r1,z1)*tmu
         psict=psict+psic/fitot
      enddo 
      coef=psict
!
      return
      end subroutine vgrid
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          a1coef computes the response functions due to           **
!**          the thin flux loops and the advance divertor coil.      **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          16/08/90..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine a1coef(coef,  nl, ne)
      use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,&
                      nesum,nfsum,nvsum,nvesel,nacoil
      use coilsp
      use consta
      use cacoil
      use nio
      use siloop
      implicit integer*4 (i-n), real*8 (a-h, o-z)
!
      radeg=pi/180.
      isplit=17
      itot=isplit*isplit
      fitot=itot
!
      m=nl
      k=ne
      psict=0
      aaa=0.0
      bbb=0.0
      call splitc(isplit,rsplt,zsplt, &
                  racoil(k),zacoil(k),wacoil(k),hacoil(k),aaa,bbb)
      do l=1,itot
         a=rsplt(l)
         r1=rsi(m)
         z1=zsi(m)-zsplt(l)
         psic=psical(a,r1,z1)*tmu
         psict=psict+psic/fitot
      enddo 
      coef=psict
!
      return
      end subroutine a1coef
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          a2coef computes the response functions due to           **
!**          the magnetic probes and the advance divertor coil.      **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          16/08/90..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine a2coef(coef, mp, ne)
      use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,&
                      nesum,nfsum,nvsum,nvesel,nacoil
      use coilsp
      use mprobe
      use consta
      use cacoil
      implicit integer*4 (i-n), real*8 (a-h, o-z)
!
      radeg=pi/180.
      isplit=17
      itot=isplit*isplit
      fitot=itot
!
      m=mp
      if (smp2(m).gt.0.0) then
         sinm=sin(radeg*amp2(m))
         cosm=cos(radeg*amp2(m))
         delsx=smp2(m)/nsmp2*cosm
         delsy=smp2(m)/nsmp2*sinm
      else
!------------------------------------------------------------------------------
!--      perpendicular probes                                                --
!------------------------------------------------------------------------------
         sinm=sin(radeg*amp2(m))
         cosm=cos(radeg*amp2(m))
         sinms=sin(radeg*(amp2(m)+90.))
         cosms=cos(radeg*(amp2(m)+90.))
         delsx=abs(smp2(m))/nsmp2*cosms
         delsy=abs(smp2(m))/nsmp2*sinms
      endif
      xmp20=xmp2(m)-(nsmp2-1)/2.*delsx
      ymp20=ymp2(m)-(nsmp2-1)/2.*delsy
      k=ne
      brct=0
      bzct=0
      aaa=0.0
      bbb=0.0
      call splitc(isplit,rsplt,zsplt, &
                  racoil(k),zacoil(k),wacoil(k),hacoil(k),aaa,bbb)
      do l=1,itot
         a=rsplt(l)
         do mmm=1,nsmp2
            r1=xmp20+(mmm-1)*delsx
            z1=ymp20+(mmm-1)*delsy-zsplt(l)
            brc=br(a,r1,z1)*tmu
            bzc=bz(a,r1,z1)*tmu
            brct=brct+brc/fitot
            bzct=bzct+bzc/fitot
         enddo 
      enddo 
      coef=(brct*cosm+bzct*sinm)/nsmp2
!
      return
      end subroutine a2coef
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          agrid computes the response functions due to            **
!**          the grid points and the advance divertor coil.          **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          16/08/90..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine agrid(coef, rgrid, nr, zgrid, nz, ne)
      use exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,&
                      nesum,nfsum,nvsum,nvesel,nacoil
      use coilsp
      use consta
      use cacoil
      use nio
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      real*8,dimension(nr) :: rgrid
      real*8,dimension(nz) :: zgrid
      data init/0/
!
      radeg=pi/180.
      isplit=17
      itot=isplit*isplit
      fitot=itot
!
      k=ne
      psict=0
      aaa=0.0
      bbb=0.0
      call splitc(isplit,rsplt,zsplt, &
                  racoil(k),zacoil(k),wacoil(k),hacoil(k),aaa,bbb)
      do l=1,itot
         a=rsplt(l)
         r1=rgrid(nr)
         z1=zgrid(nz)-zsplt(l)
         psic=psical(a,r1,z1)*tmu
         psict=psict+psic/fitot
      enddo 
      coef=psict
!
      return
      end subroutine agrid
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**       xmdele computes the elliptic integral e.                   **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**       xm1.............argument of elliptic integral e            **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**********************************************************************
      function xmdele(xm1)
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension a(4),b(4)
      real*8 a,b,xm1,xmdele
      data a(1),a(2),a(3),a(4)/.44325141463,.06260601220,&
        .04757383546,.01736506451/
      data b(1),b(2),b(3),b(4)/.24998368310,.09200180037,&
        .04069697526,.00526449639/
!
      xmdele=1.0+xm1*(a(1)+xm1*(a(2)+xm1*(a(3)+xm1*a(4))))&
       +xm1*(b(1)+xm1*(b(2)+xm1*(b(3)+xm1*b(4))))*log(1.0/xm1)
      return
      end function xmdele
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**       xmdelk computes the elliptic integral k.                   **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**       xm1.............argument of elliptic integral k            **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**********************************************************************
      function xmdelk(xm1)
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension a(5),b(5)
      real*8 a,b,xm1,xmdelk
      data a(1),a(2),a(3),a(4),a(5)/1.38629436112,.09666344259,&
        .03590092383,.03742563713,.01451196212/
      data b(1),b(2),b(3),b(4),b(5)/.5,.12498593597,.06880248576,&
        .03328355346,.00441787012/
!
      xmdelk=a(1)+xm1*(a(2)+xm1*(a(3)+xm1*(a(4)+xm1*a(5))))&
       +(b(1)+xm1*(b(2)+xm1*(b(3)+xm1*(b(4)+xm1*b(5)))))&
       *log(1.0/xm1)
      return
      end function xmdelk
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          for SCCS control revision information.                  **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          11/07/95..........first created                         **
!**                                                                  **
!**********************************************************************
      subroutine efundu_rev(i)
      character*100 opt
      character*10 s
      if( i .eq. 0)  &
      s='@(#)efund.for,v 2.3 1996/10/17 15:53:28 lao exp\000'
      return
      end subroutine efundu_rev
