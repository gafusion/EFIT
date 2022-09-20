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
      PROGRAM efund

      write(*,*) 'Reading namelist'
      CALL efund_getsizes
      CALL efund_getset
      write(*,*) 'Calling matrix subroutine'
      CALL efund_matrix
      write(*,*) 'Calling grid subroutine'
      CALL efund_grid

      STOP 'GREEN TABLE GENERATED!'
      END PROGRAM efund
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
      SUBROUTINE e1coef(coef,  nl, ne)
      USE exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,nesum,&
                      nfsum,nvsum,nvesel,nacoil
      USE cecoil
      USE coilsp
      USE consta
      USE nio
      USE siloop
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
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
      CALL splitc(isplit,rsplt,zsplt, &
                  re(k),ze(k),we(k),he(k),aaa,bbb)
      DO l=1,itot
         a=rsplt(l)
         r1=rsi(m)
         z1=zsi(m)-zsplt(l)
         psic=psical(a,r1,z1)*tmu
         psict=psict+psic/fitot
      ENDDO 
      coef=psict
!
      RETURN
      END SUBROUTINE e1coef
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
      SUBROUTINE e2coef(coef, mp, ne)
      USE exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,nesum,&
                      nfsum,nvsum,nvesel,nacoil
      USE cecoil
      USE coilsp
      USE consta
      USE mprobe
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
!
      radeg=pi/180.
      isplit=17
      itot=isplit*isplit
      fitot=itot
!
      m=mp
      IF (smp2(m).gt.0.0) THEN
         sinm=sin(radeg*amp2(m))
         cosm=cos(radeg*amp2(m))
         delsx=smp2(m)/nsmp2*cosm
         delsy=smp2(m)/nsmp2*sinm
      ELSE
!------------------------------------------------------------------------------
!--  perpendicular probes    96/02/04                                        --
!------------------------------------------------------------------------------
         sinm=sin(radeg*amp2(m))
         cosm=cos(radeg*amp2(m))
         sinms=sin(radeg*(amp2(m)+90.))
         cosms=cos(radeg*(amp2(m)+90.))
         delsx=abs(smp2(m))/nsmp2*cosms
         delsy=abs(smp2(m))/nsmp2*sinms
      ENDIF
      xmp20=xmp2(m)-(nsmp2-1)/2.*delsx
      ymp20=ymp2(m)-(nsmp2-1)/2.*delsy
      k=ne
      brct=0
      bzct=0
      aaa=0.0
      bbb=0.0
      CALL splitc(isplit,rsplt,zsplt, &
                  re(k),ze(k),we(k),he(k),aaa,bbb)
      DO l=1,itot
         a=rsplt(l)
         DO mmm=1,nsmp2
            r1=xmp20+(mmm-1)*delsx
            z1=ymp20+(mmm-1)*delsy-zsplt(l)
            brc=br(a,r1,z1)*tmu
            bzc=bz(a,r1,z1)*tmu
            brct=brct+brc/fitot
            bzct=bzct+bzc/fitot
         ENDDO 
      ENDDO 
      coef=(brct*cosm+bzct*sinm)/nsmp2
!
      RETURN
      END SUBROUTINE e2coef
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
      SUBROUTINE egrid(coef, rgrid, nr, zgrid, nz, ne)
      USE exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,nesum,&
                      nfsum,nvsum,nvesel,nacoil
      USE cecoil
      USE coilsp
      USE consta
      USE nio
      USE mprobe
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
      REAL*8,DIMENSION(nr) :: rgrid
      REAL*8,DIMENSION(nz) :: zgrid
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
      CALL splitc(isplit,rsplt,zsplt, &
                  re(k),ze(k),we(k),he(k),aaa,bbb)
      DO l=1,itot
         a=rsplt(l)
         r1=rgrid(nr)
         z1=zgrid(nz)-zsplt(l)
         psic=psical(a,r1,z1)*tmu
         psict=psict+psic/fitot
      ENDDO 
      coef=psict
      RETURN
      END SUBROUTINE egrid
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
!**       h1,h2...........height                                     **
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
      SUBROUTINE gacoil(rsilac,rmp2ac,gridac,rgrid,mw, zgrid,mh)
      USE exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,nesum,&
                      nfsum,nvsum,nvesel,nacoil,nw,nh,nwnh
      USE consta
      USE cacoil
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
      DIMENSION rsilac(nsilop,nacoil),rmp2ac(magpr2,nacoil)
      REAL*8,DIMENSION(mw) ::  rgrid
      REAL*8,DIMENSION(mh) ::  zgrid
      REAL*8,DIMENSION(nwnh,nacoil) :: gridac
      DO j=1,nsilop
         jj=j
         DO i=1,nacoil
            ii=i
            CALL a1coef(work,jj,ii)
            rsilac(j,i)=work
         ENDDO 
      ENDDO 
!
      DO  j=1,magpr2
         jj=j
         DO  i=1,nacoil
            ii=i
            CALL a2coef(work,jj,ii)
            rmp2ac(j,i)=work
         ENDDO 
      ENDDO 
!
      DO i=1,nw
         nr=i
         DO j=1,nh
            nz=j
            kk=(i-1)*nh+j
            DO n=1,nacoil
               nn=n
               CALL agrid(work,rgrid,nr,zgrid,nz,nn)
               gridac(kk,n)=work
            ENDDO 
         ENDDO
      ENDDO 
      RETURN
      END SUBROUTINE gacoil
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
      SUBROUTINE gecoil(rsilec,rmp2ec,gridec,rgrid,mw,zgrid,mh, &
                        rfcec,recec,rsisec)
      USE exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,nesum,&
                      nfsum,nvsum,nvesel,nacoil,nw,nh,nwnh
      USE consta
      USE fcoil
      USE cecoil
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
      DIMENSION rsilec(nsilop,nesum),rmp2ec(magpr2,nesum), &
                rfcec(nfcoil,nesum),recec(nesum,nesum), &
                rsisec(nesum)
            REAL*8,DIMENSION(mw) :: rgrid
            REAL*8,DIMENSION(mh) :: zgrid
            REAL*8,DIMENSION(nwnh,nesum) :: gridec
      DIMENSION taf(nfcoil),taf2(nfcoil)
            REAL*8 :: zetaec = 3.5e-08
      DO j=1,nsilop
         DO i=1,nesum
            rsilec(j,i)=0.0
         ENDDO 
         jj=j
         DO i=1,necoil
            ii=i
            CALL e1coef(work,jj,ii)
            kkm=ecid(i)
            rsilec(j,kkm)=rsilec(j,kkm)+work*ecturn(i)
         ENDDO 
      ENDDO
!
      DO j=1,magpr2
         DO i=1,nesum
            rmp2ec(j,i)=0.0
         ENDDO 
         jj=j
         DO i=1,necoil
            ii=i
            CALL e2coef(work,jj,ii)
            kkm=ecid(i)
            rmp2ec(j,kkm)=rmp2ec(j,kkm)+work*ecturn(i)
         ENDDO 
      ENDDO
!
      DO i=1,nw
         nr=i
         DO j=1,nh
            nz=j
            kk=(i-1)*nh+j
            DO m=1,nesum
               gridec(kk,m)=0.0
            ENDDO 
            DO n=1,necoil
               nn=n
               CALL egrid(work,rgrid,nr,zgrid,nz,nn)
               kkkm=ecid(n)
               gridec(kk,kkkm)=gridec(kk,kkkm)+work*ecturn(n)
            ENDDO 
         ENDDO
      ENDDO
!
      aaa=0.0
      bbb=0.0
      DO j=1,nfcoil
         DO i=1,nesum
            rfcec(j,i)=0.0
         ENDDO 
         taf(j)=tan(af(j)*pi/180.)
         taf2(j)=tan(af2(j)*pi/180.)
         DO i=1,necoil
            CALL flux(re(i),ze(i),we(i),he(i),aaa,bbb, &
                      rf(j),zf(j),wf(j),hf(j),taf(j),taf2(j),work)
            work=work*0.5/pi
            kkm=ecid(i)
            kk=kkm+1
            rfcec(j,kkm)=rfcec(j,kkm)+work*(kk-ecid(i))
         ENDDO 
      ENDDO
!
      DO j=1,nesum
         rsisec(j)=0.0
         DO i=1,nesum
            recec(j,i)=0.0
         ENDDO 
      ENDDO 
      DO j=1,necoil
         jjm=ecid(j)
         DO i=1,necoil
            CALL flux(re(i),ze(i),we(i),he(i),aaa,bbb, &
                      re(j),ze(j),we(j),he(j),aaa,bbb,work)
            work=work*0.5/pi
            kkm=ecid(i)
            recec(jjm,kkm)=recec(jjm,kkm)+work
         ENDDO 
         rsisec(jjm)=rsisec(jjm)+2.*pi*re(j)/we(j)/he(j)*zetaec
      ENDDO
      RETURN
      END SUBROUTINE gecoil
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
      if (igrid.le.0) return
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
      SUBROUTINE gsilop(rr, nr, zz, nz, rspfun, ns, rsi, zsi, wsi &
           , hsi, as, as2, ndim)
      USE exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,&
                      nw,nh,nwnh
      USE consta
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
            REAL*8,DIMENSION(nr) :: rr
            REAL*8,DIMENSION(nz) :: zz
            REAL*8,DIMENSION(ns) :: rsi,zsi,wsi,hsi,as,as2
            REAL*8,DIMENSION(ndim,nwnh) :: rspfun
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
      END SUBROUTINE gsilop
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
      SUBROUTINE gvesel(rsilvs,rmp2vs,gridvs,rgrid,mw, &
                        zgrid,mh,rfcvs,rvsfc,rvsec)
      USE exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,nesum,&
                      nfsum,nvsum,nvesel,nacoil,nw,nh,nwnh
      USE consta
      USE fcoil
      USE cecoil
      USE cvesel
      USE input,only:iecoil
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
      DIMENSION rsilvs(nsilop,nvesel),rmp2vs(magpr2,nvesel), &
                rgrid(1),zgrid(1),rvsec(nvesel,nesum), &
                rfcvs(nfcoil,nvesel),rvsfc(nvesel,nfcoil)
      DIMENSION gridvs(mw*mh,nvesel)
      DIMENSION taf(nfcoil),taf2(nfcoil)
      DIMENSION tas(nvesel),tas2(nvesel)
      DO j=1,nsilop
         jj=j
         DO i=1,nvesel
            ii=i
            CALL v1coef(work,jj,ii)
            rsilvs(j,i)=work
         ENDDO 
      ENDDO 
!
      DO j=1,magpr2
         jj=j
         DO i=1,nvesel
            ii=i
            CALL v2coef(work,jj,ii)
            rmp2vs(j,i)=work
         ENDDO 
      ENDDO 
!
      DO i=1,nw
         nr=i
         DO j=1,nh
            nz=j
            kk=(i-1)*nh+j
            DO n=1,nvesel
               nn=n
               CALL vgrid(work,rgrid,nr,zgrid,nz,nn)
               gridvs(kk,n)=work
            ENDDO 
         ENDDO 
      ENDDO
!
      DO j=1,nfcoil
         taf(j)=tan(af(j)*pi/180.)
         taf2(j)=tan(af2(j)*pi/180.)
         DO i=1,nvesel
            tas(i)=tan(avs(i)*pi/180.)
            tas2(i)=tan(avs2(i)*pi/180.)
            CALL flux(rvs(i),zvs(i),wvs(i),hvs(i),tas(i),tas2(i), &
                      rf(j),zf(j),wf(j),hf(j),taf(j),taf2(j),work)
            work=work*0.5/pi
            rfcvs(j,i)=work
         ENDDO 
      ENDDO 
!
      aaa=0.0
      rvsec(:,:) = 0.0
      DO j=1,nvesel
         tas(j)=tan(avs(j)*pi/180.)
         tas2(j)=tan(avs2(j)*pi/180.)
         DO i=1,nesum
            rvsec(j,i)=0.0
         ENDDO
         IF  (iecoil.eq.1) then
            DO i=1,necoil
               CALL flux(re(i),ze(i),we(i),he(i),aaa,aaa, &
                      rvs(j),zvs(j),wvs(j),hvs(j),tas(j),tas2(j),work)
               work=work*0.5/pi
               kkm=ecid(i)
               kk=kkm+1
               rvsec(j,kkm)=rvsec(j,kkm)+work*(kk-ecid(i))
            ENDDO
         ENDIF
      ENDDO 
!
      DO j=1,nvesel
         tas(j)=tan(avs(j)*pi/180.)
         tas2(j)=tan(avs2(j)*pi/180.)
         DO i=1,nfcoil
            taf(i)=tan(af(i)*pi/180.)
            taf2(i)=tan(af2(i)*pi/180.)
            CALL flux(rf(i),zf(i),wf(i),hf(i),taf(i),taf2(i), &
                      rvs(j),zvs(j),wvs(j),hvs(j),tas(j),tas2(j),work)
            work=work*0.5/pi
            rvsfc(j,i)=work
         ENDDO 
      ENDDO 
      RETURN
      END SUBROUTINE gvesel
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
      SUBROUTINE lgauss(x,w,n,nn)
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
      REAL*8,DIMENSION(n) ::  x,w
!
      nn = 0
      IF (n-1.lt.0) then 
         nn = 1
         RETURN
      ELSEIF (n-1.eq.0) then
!----------------------------------------------------------------------
!--      request for a zero point formula is meaningless             --
!----------------------------------------------------------------------
         x(1) = 0.
         w(1) = 2.
         RETURN
      ENDIF
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
      DO i = 1,n
         test = -2.
         ic = n+1-i
!----------------------------------------------------------------------
!--      whenever we find a root of the                              --
!--      polynomial, its negative is also a root.                    --
!--      the index ic tells WHERE to store the other root            --
!----------------------------------------------------------------------
         IF (ic.lt.i) RETURN
   40    s = g
         t = 1.
         u = 1.
         v = 0.
!----------------------------------------------------------------------
!--      evaluation of the n-th legendre polynomial                  --
!--      and its first derivative.                                   --
!--      WHERE   u = ds/dx                                           --
!--              v = dt/dx                                           --
!--              dp=dp/dx                                            --
!----------------------------------------------------------------------
         DO k = 2,n
            a = k
            p = ((2.0*a-1.0)*s*g-(a-1.0)*t)/a
             dp = ((2.0*a-1.0)*(s+g*u)-(a-1.0)*v)/a
            v = u
            u = dp
            t = s
            s = p
         ENDDO 
         IF (abs((test-g)/(test+g)).ge.0.0000005) THEN
            sum = 0.
            IF (i.ne.1) THEN
!----------------------------------------------------------------------
!--            the following computes the reduced                    --
!--            legendre polynomial and its derivative.               --
!----------------------------------------------------------------------
               DO k = 2,i
                  sum = sum+1./(g-x(k-1))
               ENDDO
            ENDIF
            test = g
            g = g-p/(dp-p*sum)
            go to 40
         ENDIF
         x(ic) = -g
         x(i) = g
         w(i) = 2./(r*t*dp)
         w(ic) = w(i)
         g = g-r*t/((r+2.)*g*dp+r*v-2.*r*t*sum)
      ENDDO
      RETURN
      END SUBROUTINE lgauss
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
      SUBROUTINE m1coef(rr, zz, nr, nz, coef,  nl, nf)
      USE exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil
      USE fcoil
      USE coilsp
      USE consta
      USE nio
      USE siloop
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
      REAL*8,DIMENSION(nr) :: rr
      REAL*8,DIMENSION(nz) :: zz
      REAL*8,DIMENSION(nsilop,nr*nz) :: coef
!
      radeg=pi/180.
      isplit=17
      itot=isplit*isplit
      fitot=itot
!
      m=nl
      IF (nf.le.0) THEN
         DO ii=1,nr
            DO jj=1,nz
               kk=(ii-1)*nz+jj
               a=rr(ii)
               r=rsi(m)
               z=zsi(m)-zz(jj)
               cmp2=psical(a,r,z)*tmu
               coef(m,kk)=cmp2
            ENDDO 
         ENDDO 
      ELSE
         k=nf
         psict=0
         CALL splitc(isplit,rsplt,zsplt, &
                     rf(k),zf(k),wf(k),hf(k),af(k),af2(k))
         DO l=1,itot
            a=rsplt(l)
            r1=rsi(m)
            z1=zsi(m)-zsplt(l)
            psic=psical(a,r1,z1)*tmu
            psict=psict+psic/fitot
         ENDDO 
         coef(m,k)=psict
      ENDIF
!
      RETURN
      END SUBROUTINE m1coef
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
      SUBROUTINE m2coef(rr, nr, zz, nz, coef,  mp, nc)
      USE exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,&
                      nw,nh,nwnh
      USE fcoil
      USE coilsp
      USE mprobe
      USE pmodel
      USE consta
      USE fshift
      USE bfgrid
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
            REAL*8,DIMENSION(nr) :: rr
            REAL*8,DIMENSION(nz) :: zz
      DIMENSION coef(mp,nc)

      radeg=pi/180.
      isplit=17
      itot=isplit*isplit
      fitot=itot
!
      DO m=1,magpr2
         IF (smp2(m).gt.0.0) THEN
            sinm=sin(radeg*amp2(m))
            cosm=cos(radeg*amp2(m))
            delsx=smp2(m)/nsmp2*cosm
            delsy=smp2(m)/nsmp2*sinm
         ELSE
!------------------------------------------------------------------------------
!--      perpendicular probes                                                --
!------------------------------------------------------------------------------
            sinm=sin(radeg*amp2(m))
            cosm=cos(radeg*amp2(m))
            sinms=sin(radeg*(amp2(m)+90.))
            cosms=cos(radeg*(amp2(m)+90.))
            delsx=abs(smp2(m))/nsmp2*cosms
            delsy=abs(smp2(m))/nsmp2*sinms
         ENDIF
         xmp20=xmp2(m)-(nsmp2-1)/2.*delsx
         ymp20=ymp2(m)-(nsmp2-1)/2.*delsy
         IF (nz.gt.0) THEN
            DO ii=1,nr
               DO jj=1,nz
                  kk=(ii-1)*nz+jj
                  a=rr(ii)
                  brct=0.0
                  bzct=0.0
                  DO mmm=1,nsmp2
                     r=xmp20+(mmm-1)*delsx
                     z=ymp20+(mmm-1)*delsy-zz(jj)
                     brtmp=br(a,r,z)*tmu
                     bztmp=bz(a,r,z)*tmu
                     brct=brtmp+brct
                     bzct=bztmp+bzct
                  ENDDO 
                  cmp2=(brct*cosm+bzct*sinm)/nsmp2
                  coef(m,kk)=cmp2
               ENDDO
            ENDDO
         ELSE
            DO k=1,nfcoil
!---------------------------------------------------------------
!--         Shifted F-coil                                    --
!---------------------------------------------------------------
               IF (k.eq.nshiftrz(k)) THEN
                  pmnow=radeg*pmprobe(m)
                  pfnow=radeg*pshift(k)
               ENDIF
!
               brct=0
               bzct=0
               CALL splitc(isplit,rsplt,zsplt, &
                           rf(k),zf(k),wf(k),hf(k),af(k),af2(k))
               DO l=1,itot
                  a=rsplt(l)
                  DO mmm=1,nsmp2
                     r1=xmp20+(mmm-1)*delsx
                     z1=ymp20+(mmm-1)*delsy-zsplt(l)
!---------------------------------------------------------------
!--                  Shifted F-coil                           --
!---------------------------------------------------------------
                     IF (k.eq.nshiftrz(k)) THEN
                        rcos=r1*cos(pmnow)-rshift(k)*cos(pfnow)
                        rsin=r1*sin(pmnow)-rshift(k)*sin(pfnow)
                        r1=sqrt(rcos**2+rsin**2)
                        z1=z1-zshift(k)
                     ENDIF
!
                     brc=br(a,r1,z1)*tmu
                     bzc=bz(a,r1,z1)*tmu
                     bzct=bzct+bzc/fitot
!---------------------------------------------------------------
!--                  Shifted F-coil ?                         --
!---------------------------------------------------------------
                     IF (k.eq.nshiftrz(k)) THEN
                        cospp=rcos/r1
                        sinpp=rsin/r1
                        cfactor=cos(pmnow)*cospp+sin(pmnow)*sinpp
                        brct=brct+brc/fitot*cfactor
                     ELSE
                        brct=brct+brc/fitot
                     ENDIF
                  ENDDO
               ENDDO
!
               coef(m,k)=(brct*cosm+bzct*sinm)/nsmp2
            ENDDO
         ENDIF
      ENDDO
!----------------------------------------------------------------
!--   BR, BZ at grid due to F-coils                            --
!----------------------------------------------------------------
      IF (nz.gt.0) THEN
         return
      ELSE
         DO ii=1,nw
            DO jj=1,nh
               kk=(ii-1)*nh+jj
               DO k=1,nfcoil
                  brct=0
                  bzct=0
                  CALL splitc(isplit,rsplt,zsplt, &
                              rf(k),zf(k),wf(k),hf(k),af(k),af2(k))
                  DO l=1,itot
                     a=rsplt(l)
                     r1=rgrid(ii)
                     z1=zgrid(jj)-zsplt(l)
                     brc=br(a,r1,z1)*tmu
                     bzc=bz(a,r1,z1)*tmu
                     bzct=bzct+bzc/fitot
                     brct=brct+brc/fitot
                  ENDDO 
                  brgridfc(kk,k)=brct
                  bzgridfc(kk,k)=bzct
               ENDDO
            ENDDO
         ENDDO
      ENDIF
!
      RETURN
      END SUBROUTINE m2coef
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
      SUBROUTINE efund_matrix
      USE exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,nesum,&
                      nfsum,nvsum,nvesel,nacoil,nw,nh,nwnh
      USE consta
      USE nio
      USE cvesel
      USE input
      USE cacoil
      USE pmodel
      USE siloop
      USE fcoil
      USE fshift
      USE bfgrid
!vas
      use var_filech
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
      DIMENSION rfcfc(nfcoil,nfcoil)
      DIMENSION rsilfc(nsilop,nfcoil),rmp2fc(magpr2,nfcoil), &
                rgowfc(nrogow,nfcoil)
      DIMENSION gsilfc(nsilop,nfsum),gmp2fc(magpr2,nfsum)
      DIMENSION rsilec(nsilop,nesum),rmp2ec(magpr2,nesum), &
                rfcec(nfcoil,nesum), &
                recec(nesum,nesum),rsisec(nesum)
      DIMENSION rsilvs(nsilop,nvesel),rmp2vs(magpr2,nvesel), &
                rfcvs(nfcoil,nvesel), &
                rvsec(nvesel,nesum),rvsfc(nvesel,nfcoil), &
                rvsvs(nvesel,nvesel),tas(nvesel),tas2(nvesel)
      DIMENSION gsilvs(nsilop,nvsum),gmp2vs(magpr2,nvsum)
      DIMENSION taf(nfcoil),taf2(nfcoil)
      DIMENSION rsilac(nsilop,nacoil),rmp2ac(magpr2,nacoil)
      DIMENSION xdum(1),ydum(1)
      REAL*8,DIMENSION(:,:),ALLOCATABLE :: rfcpc,brgrfc,bzgrfc, &
                                         rsilpc,rmp2pc,rgowpc, &
                                         gridec,gridvs,ggridvs, &
                                         gridac
!
      IF (.NOT. ALLOCATED(rfcpc)) THEN
        ALLOCATE(rfcpc(nfcoil,nwnh))
        rfcpc(:,:) = 0.0
      ENDIF
      IF (.NOT. ALLOCATED(brgrfc)) THEN
        ALLOCATE(brgrfc(nwnh,nfsum))
        brgrfc(:,:) = 0.0
      ENDIF
      IF (.NOT. ALLOCATED(bzgrfc)) THEN
        ALLOCATE(bzgrfc(nwnh,nfsum))
        bzgrfc(:,:) = 0.0
      ENDIF
      IF (.NOT. ALLOCATED(rsilpc)) THEN
        ALLOCATE(rsilpc(nsilop,nwnh))
        rsilpc(:,:) = 0.0
      ENDIF
      IF (.NOT. ALLOCATED(rmp2pc)) THEN
        ALLOCATE(rmp2pc(magpr2,nwnh))
        rmp2pc(:,:) = 0.0
      ENDIF
      IF (.NOT. ALLOCATED(rgowpc)) THEN
        ALLOCATE(rgowpc(nrogow,nwnh))
        rgowpc(:,:) = 0.0
      ENDIF
      IF (.NOT. ALLOCATED(gridec)) THEN
        ALLOCATE(gridec(nwnh,nesum))
        gridec(:,:) = 0.0
      ENDIF
      IF (.NOT. ALLOCATED(gridvs)) THEN
        ALLOCATE(gridvs(nwnh,nvesel))
        gridvs(:,:) = 0.0
      ENDIF
      IF (.NOT. ALLOCATED(ggridvs)) THEN
        ALLOCATE(ggridvs(nwnh,nvsum))
        ggridvs(:,:) = 0.0
      ENDIF
      IF (.NOT. ALLOCATED(gridac)) THEN
        ALLOCATE(gridac(nwnh,nacoil))
        gridac(:,:) = 0.0
      ENDIF
!
      IF (ifcoil.gt.0) THEN
!----------------------------------------------------------------------
!--      calculate the response FUNCTION of psi loops due to f coils --
!----------------------------------------------------------------------
         DO i=1,nfcoil
            taf(i)=tan(af(i)*pi/180.)
            taf2(i)=tan(af2(i)*pi/180.)
         ENDDO 
         IF (islpfc.gt.0) THEN
            DO i=1,nfcoil
               DO j=1,nfcoil
                  CALL flux(rf(i),zf(i),wf(i),hf(i),taf(i),taf2(i), &
                            rf(j),zf(j),wf(j),hf(j),taf(j),taf2(j), &
                            rfcfc(j,i))
                  rfcfc(j,i)=rfcfc(j,i)*0.5/pi

               ENDDO
               ii=i
               CALL gsilop(rgrid,nw,zgrid,nh,rfcpc,ii,rf,zf,wf,hf,af,af2 &
                          ,nfcoil)
            ENDDO

            print*,'file name : ','fc'//trim(ch1)//trim(ch2)//'.ddd' 
!vasorg      OPEN(unit=nrspfc,status='unknown',file='fcfcpc.dat', &
            OPEN(unit=nrspfc,status='unknown',file='fc'//trim(ch1)// & 
                 trim(ch2)//'.ddd' , &
                 form='unformatted')
            WRITE (nrspfc) rfcfc
            WRITE (nrspfc) rfcpc
            CLOSE(unit=nrspfc)
         ENDIF
!---------------------------------------------------------------------
!--      flux loops                                                 --
!---------------------------------------------------------------------
         IF (nsilop.gt.1) THEN
!---------------------------------------------------------------------
!           finite size flux loops                                  --
!---------------------------------------------------------------------
            IF (isize.gt.0) THEN
               DO i=1,nfcoil
                  DO j=1,isize
                     taz=tan(as(j)*pi/180.)
                     taz2=tan(as2(j)*pi/180.)
                     CALL flux(rf(i),zf(i),wf(i),hf(i),taf(i),taf2(i), &
                               rsi(j),zsi(j),wsi(j),hsi(j),taz,taz2, &
                               rsilfc(j,i))
                     rsilfc(j,i)=rsilfc(j,i)*0.5/pi
                  ENDDO 
               ENDDO
            ENDIF
!---------------------------------------------------------------------
!           thin flux loops                                         --
!---------------------------------------------------------------------
            IF (isize.lt.nsilop) then
               DO i=1,nfcoil
                  ii=i
                  DO j=isize+1,nsilop
                     jj=j
                     CALL m1coef(xdum,xdum,nsilop,nfcoil,rsilfc,jj,ii)
                  ENDDO 
               ENDDO 
            ENDIF 
         ENDIF 

         !
         IF (.NOT. ALLOCATED(brgridfc)) THEN
            ALLOCATE(brgridfc(nwnh,nfcoil))
            brgridfc(:,:) = 0.0
         ENDIF
         IF (.NOT. ALLOCATED(bzgridfc)) THEN
            ALLOCATE(bzgridfc(nwnh,nfcoil))
            bzgridfc(:,:) = 0.0
         ENDIF
!
!-----------------------------------------------------------------------
!--      compute the response FUNCTION of magnetic probes due to f coils
!-----------------------------------------------------------------------
         magprr=magpr2
         IF (magprr.gt.1) THEN
            CALL m2coef(xdum,0,ydum,0,rmp2fc,magpr2,nfcoil)
         ENDIF
!----------------------------------------------------------------------
!--      compute the response FUNCTION of partial rogowski loops due to
!--      f coils
!----------------------------------------------------------------------
         mrogow=nrogow
         IF (mrogow.gt.1) THEN
            CALL rogowc(xdum,0,ydum,0,rgowfc,nrogow,nfcoil)
         ENDIF 
!----------------------------------------------------------------------
!--      WRITE f coil response functions                             --
!----------------------------------------------------------------------
         DO i=1,nfsum
            DO j=1,nsilop
               gsilfc(j,i)=0.0
            ENDDO
            DO j=1,magpr2
               gmp2fc(j,i)=0.0
            ENDDO
         ENDDO
!
         DO i=1,nfcoil
            k=abs(fcid(i))
            DO j=1,nsilop
               gsilfc(j,k)=gsilfc(j,k)+fcturn(i)*rsilfc(j,i)
            ENDDO
            DO j=1,magpr2
               gmp2fc(j,k)=gmp2fc(j,k)+fcturn(i)*rmp2fc(j,i)
            ENDDO
         ENDDO
!
         print*,'file name : ','rfcoil.ddd' 
!vasorg      OPEN(unit=nrspfc,status='unknown',file='rfcoil.dat', &
         OPEN(unit=nrspfc,status='unknown',file='rfcoil.ddd', &
              form='unformatted')
         WRITE (nrspfc) gsilfc
         WRITE (nrspfc) gmp2fc
         CLOSE(unit=nrspfc)
!
         DO i=1,nfsum
            DO j=1,nwnh
               brgrfc(j,i)=0.0
               bzgrfc(j,i)=0.0
            ENDDO
         ENDDO
         DO i=1,nfcoil
            k=abs(fcid(i))
            DO j=1,nwnh
               brgrfc(j,k)=brgrfc(j,k)+fcturn(i)*brgridfc(j,i)
               bzgrfc(j,k)=bzgrfc(j,k)+fcturn(i)*bzgridfc(j,i)
            ENDDO
         ENDDO
!
         OPEN(unit=nrspfc,status='unknown',file='brzgfc.dat', &
              form='unformatted')
         WRITE (nrspfc) brgrfc
         WRITE (nrspfc) bzgrfc
         CLOSE(unit=nrspfc)
      ENDIF
!----------------------------------------------------------------------
!--   plasma response functions                                      --
!----------------------------------------------------------------------
      IF (igrid.gt.0) THEN
         msilop=nsilop
         IF (msilop.gt.1) THEN
!----------------------------------------------------------------------
!--         filament plasma current model                            --
!----------------------------------------------------------------------
            IF (isize.gt.0) THEN
               DO j=1,isize
                  jj=j
                  CALL gsilop(rgrid,nw,zgrid,nh,rsilpc,jj, &
                              rsi,zsi,wsi,hsi,as,as2,nsilop)
               ENDDO
            ENDIF
         ENDIF
         IF (isize.lt.nsilop) THEN
            DO j=isize+1,nsilop
               jj=j
               CALL m1coef(rgrid,zgrid,nw,nh,rsilpc,jj,0)
            ENDDO 
         ENDIF
         magprr=magpr2
         IF (magprr.gt.1) THEN
            CALL m2coef(rgrid,nw,zgrid,nh,rmp2pc,magpr2,nwnh)
         ENDIF
         mrogow=nrogow
         IF (mrogow.gt.1) THEN
            CALL rogowc(rgrid,nw,zgrid,nh,rgowpc,nrogow,nwnh)
         ENDIF
!----------------------------------------------------------------------
!--      WRITE the plasma response FUNCTION                          --
!----------------------------------------------------------------------
         print*,'file name : ','ep'//trim(ch1)//trim(ch2)//'.ddd' 
!vasorg      OPEN(unit=nrsppc,status='unknown',file='eplasm.dat', &
         OPEN(unit=nrsppc,status='unknown',file='ep'//trim(ch1)// & 
              trim(ch2)//'.ddd' , &
              form='unformatted')
         WRITE (nrsppc) rsilpc
         WRITE (nrsppc) rmp2pc
         CLOSE(unit=nrsppc)
!
      ENDIF
      IF (iecoil.gt.0) THEN
         CALL gecoil(rsilec,rmp2ec,gridec,rgrid,nw, &
                     zgrid,nh,rfcec,recec,rsisec)
      ENDIF
!
      print*,'file name : ','re'//trim(ch1)//trim(ch2)//'.ddd' 
!vasorg      OPEN(unit=nrsppc,status='unknown',file='recoil.dat', &
      OPEN(unit=nrsppc,status='unknown',file='re'//trim(ch1)// & 
                     trim(ch2)//'.ddd', &
           form='unformatted')
      WRITE (nrsppc) rsilec
      WRITE (nrsppc) rmp2ec
      WRITE (nrsppc) gridec
      CLOSE(unit=nrsppc)
!
      IF (ivesel.gt.0) THEN
         CALL gvesel(rsilvs,rmp2vs,gridvs,rgrid,nw, &
                     zgrid,nh,rfcvs,rvsfc,rvsec)
         DO i=1,nvesel
           tas(i)=tan(avs(i)*pi/180.)
           tas2(i)=tan(avs2(i)*pi/180.)
         ENDDO 
         DO i=1,nvesel
            DO j=1,nvesel
               CALL flux(rvs(i),zvs(i),wvs(i),hvs(i),tas(i),tas2(i), &
                         rvs(j),zvs(j),wvs(j),hvs(j),tas(j),tas2(j), &
                         rvsvs(j,i))
               rvsvs(j,i)=rvsvs(j,i)*0.5/pi
            ENDDO 
         ENDDO 
!
         DO i=1,nvsum
            DO j=1,nsilop
               gsilvs(j,i)=0.0
            ENDDO
            DO j=1,magpr2
               gmp2vs(j,i)=0.0
            ENDDO
            DO j=1,nwnh
               ggridvs(j,i)=0.0
            ENDDO
         ENDDO
!
         DO i=1,nvesel
            k=abs(vsid(i))
            DO j=1,nsilop
               gsilvs(j,k)=gsilvs(j,k)+rsilvs(j,i)
            ENDDO
            DO j=1,magpr2
               gmp2vs(j,k)=gmp2vs(j,k)+rmp2vs(j,i)
            ENDDO
            DO j=1,nwnh
               ggridvs(j,k)=ggridvs(j,k)+gridvs(j,i)
            ENDDO
         ENDDO
!
!vas
      print*,'file name : ','rv'//trim(ch1)//trim(ch2)//'.ddd' 
!vasorg         OPEN(unit=nrsppc,status='unknown',file='rvesel.dat', &
         OPEN(unit=nrsppc,status='unknown',file='rv'//trim(ch1)// & 
              trim(ch2)//'.ddd' , &
              form='unformatted')
         WRITE (nrsppc) gsilvs
         WRITE (nrsppc) gmp2vs
         WRITE (nrsppc) ggridvs
         CLOSE(unit=nrsppc)
      ENDIF
!---------------------------------------------------------------------
!--   advance divertor coil                                         --
!---------------------------------------------------------------------
      IF (iacoil.gt.0) THEN
         CALL gacoil(rsilac,rmp2ac,gridac,rgrid,nw, &
                     zgrid,nh)
!vas
      print*,'file name : ','ra'//trim(ch1)//trim(ch2)//'.ddd' 
!vasorg         OPEN(unit=nrsppc,status='unknown',file='racoil.dat', &
         OPEN(unit=nrsppc,status='unknown',file='ra'//trim(ch1)// & 
              trim(ch2)//'.ddd' , &
              form='unformatted')
         WRITE (nrsppc) gridac
         WRITE (nrsppc) rsilac
         WRITE (nrsppc) rmp2ac
         CLOSE(unit=nrsppc)
      ENDIF
!
      IF (ALLOCATED(rfcpc)) DEALLOCATE(rfcpc)
      IF (ALLOCATED(brgrfc)) DEALLOCATE(brgrfc)
      IF (ALLOCATED(bzgrfc)) DEALLOCATE(bzgrfc)
      IF (ALLOCATED(rsilpc)) DEALLOCATE(rsilpc)
      IF (ALLOCATED(rmp2pc)) DEALLOCATE(rmp2pc)
      IF (ALLOCATED(rgowpc)) DEALLOCATE(rgowpc)
      IF (ALLOCATED(gridec)) DEALLOCATE(gridec)
      IF (ALLOCATED(gridvs)) DEALLOCATE(gridvs)
      IF (ALLOCATED(ggridvs)) DEALLOCATE(ggridvs)
      IF (ALLOCATED(gridac)) DEALLOCATE(gridac)
!
      RETURN
      END SUBROUTINE efund_matrix
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
      SUBROUTINE rogowc(rr, nrr, zz, nz, coef, nr, nc)
      USE exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil
      USE rogowl
      USE coilsp
      USE consta
      USE fcoil
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
!      DIMENSION rogpth(nrogow)
      REAL*8,DIMENSION(nr) :: rr
      REAL*8,DIMENSION(nz) :: zz
      DIMENSION coef(nr,nc)
!
      ngrid=25
      isplit=17
      itot=isplit*isplit
      fitot=itot
      dels = 0.
      mm=1
!
      DO m=1,nrogow
         IF (nz.gt.0) THEN
            DO inn=1,nrr
               DO imm=1,nz
                  ikk=(inn-1)*nz+imm
                  coef(m,ikk)=0.0
               ENDDO 
            ENDDO 
         ELSE
            DO k=1,nfcoil
               coef(m,k)=0.0
            ENDDO 
         ENDIF
         k=m
         CALL rogrid(ngrid,mm,k,dels)
         mm=mm+narc(m)+1
         DO i=1,ngrid
            iii=i
            IF(i.eq.ngrid) THEN
               zl=zpg(i)-zpg(i-1)
               rl=rpg(i)-rpg(i-1)
            ELSE
               zl=zpg(i+1)-zpg(i)
               rl=rpg(i+1)-rpg(i)
            ENDIF
            hl=sqrt(zl*zl+rl*rl)
            sint=zl/hl
            cost=rl/hl
!
            IF (nz.le.0) THEN
               DO k=1,nfcoil
                  CALL splitc(isplit,rsplt,zsplt, &
                           rf(k),zf(k),wf(k),hf(k),af(k),af2(k))
                  DO l=1,itot
                     a=rsplt(l)
                     r1=rpg(i)
                     z1=zpg(i)-zsplt(l)
                     brc=br(a,r1,z1)*tmu/fitot
                     bzc=bz(a,r1,z1)*tmu/fitot
                     part=brc*cost+bzc*sint
                     CALL simpf(iii,fact)
                     ! TODO: rogpth is never defined in efit...
                     coef(m,k)=coef(m,k)+fact*part*dels !/rogpth(m)
                  ENDDO 
               ENDDO
            ELSE
               DO inn=1,nrr
                  DO imm=1,nz
                     ikk=(inn-1)*nz+imm
                     a=rr(inn)
                     r1=rpg(i)
                     z1=zpg(i)-zz(imm)
                     brg=br(a,r1,z1)*tmu
                     bzg=bz(a,r1,z1)*tmu
                     part=brg*cost+bzg*sint
                     CALL simpf(iii,fact)
                     ! TODO: rogpth is never defined in efit...
                     coef(m,ikk)=coef(m,ikk)+fact*part*dels !/rogpth(m)
                  ENDDO 
               ENDDO 
            ENDIF
         ENDDO
      ENDDO
!
      RETURN
      END SUBROUTINE rogowc
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
      SUBROUTINE rogrid(ngrid,mm,m,dels)
      USE exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil
      USE rogowl
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
      DIMENSION sl(6)
!
      s = 0.
      mm1 = mm+narc(m)-1
      j = 1
      DO i = mm,mm1
         sl(j) = sqrt((rp(i+1)-rp(i))**2 + (zp(i+1)-zp(i))**2)
         s = s+sl(j)
         j = j+1
      ENDDO 
      dels = s/dble(ngrid-1)
      ds = 0.
      i1 = 1
      DO j = 1,narc(m)
         dr = rp(mm+j)-rp(mm+j-1)
         dz = zp(mm+j)-zp(mm+j-1)
         rpg(i1) = rp(mm+j-1)+dr*ds/sl(j)
         zpg(i1) = zp(mm+j-1)+dz*ds/sl(j)
         dd = sl(j)-ds
         n1 = int(dd/dels)+i1
         i2 = i1+1
         dr = dr*dels/sl(j)
         dz = dz*dels/sl(j)
         DO i = i2,n1
            rpg(i) = rpg(i-1)+dr
            zpg(i) = zpg(i-1)+dz
         ENDDO 
         ds = dels-(dd-dble(n1-i1)*dels)
         i1 = n1+1
      ENDDO
      rpg(ngrid) = rp(mm+narc(m))
      zpg(ngrid) = zp(mm+narc(m))
      RETURN
      END SUBROUTINE rogrid
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
      SUBROUTINE simpf(i,f)
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
      IF (i.eq.1 .or. i.eq.25) THEN
         f = 1./3.
      ELSEIF ((i/2)*2.eq.i) THEN
         f = 4./3.
      ELSE
         f = 2./3.
      ENDIF
      RETURN
      END SUBROUTINE simpf
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
      SUBROUTINE soleno(ra,z1,w1,h1,t1,t12,r2,z2,w2,h2,t2,t22,xbm1,xbm2, &
                        xtm1,xtm2,hfa,hsa,rf,rs,sol)
      USE consta,only:pi
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
      DIMENSION z(2,2)
      REAL*8 ksq,kpsq,msl,mut
!      DATA init/0/
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
      IF (t12.eq.0.) THEN
         r = rf-ra+0.5*w1
         zc1 = z1-0.5*h1-0.5*w1*t1
         zb1 = zc1+t1*r
         zt1 = zb1+h1
!
      ELSE
         zb1 = z1-h1/2.
         IF (t12.lt.0.) THEN
            IF (rf .lt. xbm1) zb1 = zb1 + t12*(rf-xbm1)
            zt1 = z1 + h1/2.
            IF (rf .gt. xtm1) zt1 = zt1 - t12*(xtm1-rf)
!
         ELSE
            IF (rf.gt.xbm1) zb1 = zb1+t12*(rf-xbm1)
            zt1 = z1+h1/2.
            IF (rf.lt.xtm1) zt1 = zt1-t12*(xtm1-rf)
         ENDIF
      ENDIF
!
      IF (t22.eq.0.) THEN
         r = rs-r2+0.5*w2
         zc2 = z2-0.5*h2-0.5*w2*t2
         zb2 = zc2+t2*r
         zt2 = zb2+h2
!
      ELSE
         zb2 = z2-h2/2.
         IF (t22 .lt. 0.) THEN
            IF (rs .lt. xbm2) zb2 = zb2 + t22*(rs-xbm2)
            zt2 = z2 + h2/2.
            IF (rs .gt. xtm2) zt2 = zt2 - t22*(xtm2-rs)
!
         ELSE
            IF (rs.gt.xbm2) zb2 = zb2+t22*(rs-xbm2)
            zt2 = z2+h2/2.
            IF (rs.lt.xtm2) zt2 = zt2-t22*(xtm2-rs)
         ENDIF
      ENDIF
!
      z(1,1) = zb1
      z(2,1) = zb2
      z(1,2) = zt1
      z(2,2) = zt2
      hfa = zt1-zb1
      hsa = zt2-zb2
!
      DO i = 1,2
         DO j = 1,2
            sign = -.25
            IF (i .ne. j) sign = .25
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
!--         To avoid numerical truncation                                --
!--------------------------------------------------------------------------
            IF (kpsq .lt. 1.0e-30) kpsq = 1.0e-30
            beta = sqrt(kpsq)
            IF (beta .lt. 1.0e-30) beta = 1.0e-10
            IF (cpsq .lt. 1.0e-30) cpsq = 1.0e-10
            delta = cpsq/beta
            epsi = csq/cpsq
            zeta = 0.
            sinf = 0.
            sa = .25
!
  100       CONTINUE
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
            IF (abs(delta-1.) .gt. err) go to 100
            IF (ambsq .gt. 1.e-14) go to 100
            cay = rpi2/alpha
            pik = cay*zeta
            ek = .5*cay*(ksq+sinf)
            msl = rh*dzsq*(r1*ek-drsq*pik/r1)
            IF (csq==1.) msl = rh*dzsq*(r1*ek-cay*fr/r1*.5)
!
            mut = msl+ut*fr*r1*(cay-t*ek)
            sol = sol+sign*mut
         ENDDO
      ENDDO
!
      RETURN
      END SUBROUTINE soleno
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
      SUBROUTINE splitc(is,rs,zs,rc,zc,wc,hc,ac,ac2)
      USE consta
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
      REAL*8,DIMENSION(is*is) :: rs,zs
!
      frd=pi/180.
!----------------------------------------------------------------------
!--   rectangle                                                      --
!----------------------------------------------------------------------
      IF(ac+ac2.eq.0.) THEN
          wdelt=wc/is
          hdelt=hc/is
          rstrt=rc-wc/2.+wdelt/2.
          zstrt=zc-hc/2.+hdelt/2.
          zz=zstrt
          ic=0
          DO ii=1,is
             rr=rstrt
             DO jj=1,is
                ic=ic+1
                zs(ic)=zz
                rs(ic)=rr
                rr=rr+wdelt
             ENDDO 
             zz=zz+hdelt
          ENDDO
          RETURN
!----------------------------------------------------------------------
!--   ac .ne. 0                                                      --
!----------------------------------------------------------------------
      ELSEIF(ac.ne.0.) THEN
          side=tan(frd*ac)*wc
          hdelt=hc/is
          wdelt=wc/is
          zdelt=tan(frd*ac)*wdelt
          rstrt=rc-wc/2.+wdelt/2.
          tsid=hc+side
          zstrt =zc-tsid/2.+tsid/2.*1./is
          rr=rstrt
          ic=0
          DO ii=1,is
             zz=zstrt+(ii-1)*zdelt
             DO jj=1,is
                ic=ic+1
                zs(ic)=zz
                rs(ic)=rr
                zz=zz+hdelt
             ENDDO 
             rr=rr+wdelt
          ENDDO
          RETURN
!----------------------------------------------------------------------
!--   ac2 .ne. 0                                                     --
!----------------------------------------------------------------------
      ELSEIF(ac2.ne.0.) THEN
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
          DO ii=1,is
             rr=rstrt+(ii-1)*rdelt
             DO jj=1,is
                ic=ic+1
                zs(ic)=zz
                rs(ic)=rr
                rr=rr+wdelt
             ENDDO 
             zz=zz+hdelt
          ENDDO
      ENDIF
!
      RETURN
      END SUBROUTINE splitc
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
      SUBROUTINE v1coef(coef,  nl, ne)
      USE exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,&
                      nesum,nfsum,nvsum,nvesel,nacoil
      USE coilsp
      USE consta
      USE nio
      USE siloop
      USE cvesel
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
!
      radeg=pi/180.
      isplit=17
      itot=isplit*isplit
      fitot=itot
!
      m=nl
      k=ne
      psict=0
      CALL splitc(isplit,rsplt,zsplt, &
                  rvs(k),zvs(k),wvs(k),hvs(k),avs(k),avs2(k))
      DO l=1,itot
         a=rsplt(l)
         r1=rsi(m)
         z1=zsi(m)-zsplt(l)
         psic=psical(a,r1,z1)*tmu
         psict=psict+psic/fitot
      ENDDO 
      coef=psict
!
      RETURN
      END SUBROUTINE v1coef
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
      SUBROUTINE v2coef(coef, mp, ne)
      USE exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,&
                      nesum,nfsum,nvsum,nvesel,nacoil
      USE coilsp
      USE mprobe
      USE consta
      USE cvesel
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
!
      radeg=pi/180.
      isplit=17
      itot=isplit*isplit
      fitot=itot
!
      m=mp
      IF (smp2(m).gt.0.0) THEN
         sinm=sin(radeg*amp2(m))
         cosm=cos(radeg*amp2(m))
         delsx=smp2(m)/nsmp2*cosm
         delsy=smp2(m)/nsmp2*sinm
      ELSE
!------------------------------------------------------------------------------
!--      perpendicular probes                                                --
!------------------------------------------------------------------------------
         sinm=sin(radeg*amp2(m))
         cosm=cos(radeg*amp2(m))
         sinms=sin(radeg*(amp2(m)+90.))
         cosms=cos(radeg*(amp2(m)+90.))
         delsx=abs(smp2(m))/nsmp2*cosms
         delsy=abs(smp2(m))/nsmp2*sinms
      ENDIF
      xmp20=xmp2(m)-(nsmp2-1)/2.*delsx
      ymp20=ymp2(m)-(nsmp2-1)/2.*delsy
      k=ne
      brct=0
      bzct=0
      CALL splitc(isplit,rsplt,zsplt, &
                  rvs(k),zvs(k),wvs(k),hvs(k),avs(k),avs2(k))
      DO l=1,itot
         a=rsplt(l)
         DO mmm=1,nsmp2
            r1=xmp20+(mmm-1)*delsx
            z1=ymp20+(mmm-1)*delsy-zsplt(l)
            brc=br(a,r1,z1)*tmu
            bzc=bz(a,r1,z1)*tmu
            brct=brct+brc/fitot
            bzct=bzct+bzc/fitot
         ENDDO 
      ENDDO 
      coef=(brct*cosm+bzct*sinm)/nsmp2
!
      RETURN
      END SUBROUTINE v2coef
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
      SUBROUTINE vgrid(coef, rgrid, nr, zgrid, nz, ne)
      USE exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,&
                      nesum,nfsum,nvsum,nvesel,nacoil
      USE coilsp
      USE siloop
      USE consta
      USE nio
      USE cvesel
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
      REAL*8,DIMENSION(nr) :: rgrid
      REAL*8,DIMENSION(nz) :: zgrid
      DATA init/0/
!
      radeg=pi/180.
      isplit=17
      itot=isplit*isplit
      fitot=itot
!
      k=ne
      psict=0
      CALL splitc(isplit,rsplt,zsplt, &
                  rvs(k),zvs(k),wvs(k),hvs(k),avs(k),avs2(k))
      DO l=1,itot
         a=rsplt(l)
         r1=rgrid(nr)
         z1=zgrid(nz)-zsplt(l)
         psic=psical(a,r1,z1)*tmu
         psict=psict+psic/fitot
      ENDDO 
      coef=psict
!
      RETURN
      END SUBROUTINE vgrid
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
      SUBROUTINE a1coef(coef,  nl, ne)
      USE exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,&
                      nesum,nfsum,nvsum,nvesel,nacoil
      USE coilsp
      USE consta
      USE cacoil
      USE nio
      USE siloop
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
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
      CALL splitc(isplit,rsplt,zsplt, &
                  racoil(k),zacoil(k),wacoil(k),hacoil(k),aaa,bbb)
      DO l=1,itot
         a=rsplt(l)
         r1=rsi(m)
         z1=zsi(m)-zsplt(l)
         psic=psical(a,r1,z1)*tmu
         psict=psict+psic/fitot
      ENDDO 
      coef=psict
!
      RETURN
      END SUBROUTINE a1coef
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
      SUBROUTINE a2coef(coef, mp, ne)
      USE exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,&
                      nesum,nfsum,nvsum,nvesel,nacoil
      USE coilsp
      USE mprobe
      USE consta
      USE cacoil
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
!
      radeg=pi/180.
      isplit=17
      itot=isplit*isplit
      fitot=itot
!
      m=mp
      IF (smp2(m).gt.0.0) THEN
         sinm=sin(radeg*amp2(m))
         cosm=cos(radeg*amp2(m))
         delsx=smp2(m)/nsmp2*cosm
         delsy=smp2(m)/nsmp2*sinm
      ELSE
!------------------------------------------------------------------------------
!--      perpendicular probes                                                --
!------------------------------------------------------------------------------
         sinm=sin(radeg*amp2(m))
         cosm=cos(radeg*amp2(m))
         sinms=sin(radeg*(amp2(m)+90.))
         cosms=cos(radeg*(amp2(m)+90.))
         delsx=abs(smp2(m))/nsmp2*cosms
         delsy=abs(smp2(m))/nsmp2*sinms
      ENDIF
      xmp20=xmp2(m)-(nsmp2-1)/2.*delsx
      ymp20=ymp2(m)-(nsmp2-1)/2.*delsy
      k=ne
      brct=0
      bzct=0
      aaa=0.0
      bbb=0.0
      CALL splitc(isplit,rsplt,zsplt, &
                  racoil(k),zacoil(k),wacoil(k),hacoil(k),aaa,bbb)
      DO l=1,itot
         a=rsplt(l)
         DO mmm=1,nsmp2
            r1=xmp20+(mmm-1)*delsx
            z1=ymp20+(mmm-1)*delsy-zsplt(l)
            brc=br(a,r1,z1)*tmu
            bzc=bz(a,r1,z1)*tmu
            brct=brct+brc/fitot
            bzct=bzct+bzc/fitot
         ENDDO 
      ENDDO 
      coef=(brct*cosm+bzct*sinm)/nsmp2
!
      RETURN
      END SUBROUTINE a2coef
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
      SUBROUTINE agrid(coef, rgrid, nr, zgrid, nz, ne)
      USE exparm,only:nfcoil,nsilop,magpr2,nrogow,necoil,&
                      nesum,nfsum,nvsum,nvesel,nacoil
      USE coilsp
      USE consta
      USE cacoil
      USE nio
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
      REAL*8,DIMENSION(nr) :: rgrid
      REAL*8,DIMENSION(nz) :: zgrid
      DATA init/0/
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
      CALL splitc(isplit,rsplt,zsplt, &
                  racoil(k),zacoil(k),wacoil(k),hacoil(k),aaa,bbb)
      DO l=1,itot
         a=rsplt(l)
         r1=rgrid(nr)
         z1=zgrid(nz)-zsplt(l)
         psic=psical(a,r1,z1)*tmu
         psict=psict+psic/fitot
      ENDDO 
      coef=psict
!
      RETURN
      END SUBROUTINE agrid
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
      FUNCTION xmdele(xm1)
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
      DIMENSION a(4),b(4)
      REAL*8 a,b,xm1,xmdele
      DATA a(1),a(2),a(3),a(4)/.44325141463,.06260601220,&
        .04757383546,.01736506451/
      DATA b(1),b(2),b(3),b(4)/.24998368310,.09200180037,&
        .04069697526,.00526449639/
!
      xmdele=1.0+xm1*(a(1)+xm1*(a(2)+xm1*(a(3)+xm1*a(4))))&
       +xm1*(b(1)+xm1*(b(2)+xm1*(b(3)+xm1*b(4))))*log(1.0/xm1)
      RETURN
      END FUNCTION xmdele
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
      FUNCTION xmdelk(xm1)
      IMPLICIT INTEGER*4 (i-n), REAL*8 (a-h, o-z)
      DIMENSION a(5),b(5)
      REAL*8  a,b,xm1,xmdelk
      DATA a(1),a(2),a(3),a(4),a(5)/1.38629436112,.09666344259,&
        .03590092383,.03742563713,.01451196212/
      DATA b(1),b(2),b(3),b(4),b(5)/.5,.12498593597,.06880248576,&
        .03328355346,.00441787012/
!
      xmdelk=a(1)+xm1*(a(2)+xm1*(a(3)+xm1*(a(4)+xm1*a(5))))&
       +(b(1)+xm1*(b(2)+xm1*(b(3)+xm1*(b(4)+xm1*b(5)))))&
       *log(1.0/xm1)
      RETURN
      END FUNCTION xmdelk
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
      SUBROUTINE efundu_rev(i)
      CHARACTER*100 opt
      CHARACTER*10 s
      IF( i .eq. 0)  &
      s='@(#)efund.for,v 2.3 1996/10/17 15:53:28 lao Exp\000'
      RETURN
      END SUBROUTINE efundu_rev
