      module vessel

      use machine, only: nvesel,nvsum
      implicit none
      public

      integer*4 ivesel
      real*8, dimension (:), allocatable :: rvs,zvs,wvs,hvs,avs,avs2, &
                                            rsisvs
      integer*4, dimension(:), allocatable :: vsid

      public gvesel
      private v1coef,v2coef,vgrid

      contains
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          gvesel gets the Green's functions for the vessel        **
!**          segments.                                               **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**********************************************************************
      subroutine gvesel(rsilvs,rmp2vs,gridvs,rgrid,mw, &
                        zgrid,mh,rfcvs,rvsfc,rvsec)
      use machine, only: nfcoil,nsilop,magpri,necoil,nesum,&
                        nvesel,nw,nh
      use utils, only: pi,flux
      use fcoil
      use ecoil
      implicit none
      integer*4, intent(in) :: mw,mh
      real*8, intent(in) :: rgrid(mw),zgrid(mh)
      real*8, intent(out) :: rsilvs(nsilop,nvesel),rmp2vs(magpri,nvesel), &
                             rfcvs(nfcoil,nvesel),rvsec(nvesel,necoil), &
                             rvsfc(nvesel,nfcoil),gridvs(mw*mh,nvesel)
      integer*4 i,j,kk,n
      real*8 aaa,work
      real*8 taf(nfcoil),taf2(nfcoil),tav(nvesel),tav2(nvesel)
!
      do j=1,nsilop
         do i=1,nvesel
            call v1coef(work,j,i)
            rsilvs(j,i)=work
         enddo 
      enddo 
!
      do j=1,magpri
         do i=1,nvesel
            call v2coef(work,j,i)
            rmp2vs(j,i)=work
         enddo 
      enddo 
!
      do i=1,nw
         do j=1,nh
            kk=(i-1)*nh+j
            do n=1,nvesel
               call vgrid(work,rgrid,i,zgrid,j,n)
               gridvs(kk,n)=work
            enddo 
         enddo 
      enddo
!
      do j=1,nfcoil
         taf(j)=tan(af(j)*pi/180.)
         taf2(j)=tan(af2(j)*pi/180.)
         do i=1,nvesel
            tav(i)=tan(avs(i)*pi/180.)
            tav2(i)=tan(avs2(i)*pi/180.)
            call flux(rvs(i),zvs(i),wvs(i),hvs(i),tav(i),tav2(i), &
                      rf(j),zf(j),wf(j),hf(j),taf(j),taf2(j),work)
            work=work*0.5/pi
            rfcvs(j,i)=work
         enddo 
      enddo 
!
      aaa=0.0
      do j=1,nvesel
         tav(j)=tan(avs(j)*pi/180.)
         tav2(j)=tan(avs2(j)*pi/180.)
         if (iecoil.eq.1) then
            do i=1,necoil
               call flux(re(i),ze(i),we(i),he(i),aaa,aaa, &
                         rvs(j),zvs(j),wvs(j),hvs(j),tav(j),tav2(j),work)
               work=work*0.5/pi
               rvsec(j,i)=work
            enddo
         endif
      enddo 
!
      do j=1,nvesel
         tav(j)=tan(avs(j)*pi/180.)
         tav2(j)=tan(avs2(j)*pi/180.)
         do i=1,nfcoil
            taf(i)=tan(af(i)*pi/180.)
            taf2(i)=tan(af2(i)*pi/180.)
            call flux(rf(i),zf(i),wf(i),hf(i),taf(i),taf2(i), &
                      rvs(j),zvs(j),wvs(j),hvs(j),tav(j),tav2(j),work)
            work=work*0.5/pi
            rvsfc(j,i)=work
         enddo 
      enddo 
      return
      end subroutine gvesel
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          v1coef computes the response functions due to           **
!**          the thin flux loops and the vessel segments.            **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**********************************************************************
      subroutine v1coef(coef,nl,ne)
      use coilsp
      use utils, only: pi,tmu,psical
      use siloop
      implicit none
      integer*4, intent(in) :: nl,ne
      real*8, intent(out) :: coef
      integer*4 l
      real*8 a,psic,psict,r1,z1
      integer*4, parameter :: isplit=17,itot=isplit*isplit
      real*8, parameter :: radeg=pi/180.,fitot=itot
!
      psict=0
      call splitc(isplit,rsplt,zsplt, &
                  rvs(ne),zvs(ne),wvs(ne),hvs(ne),avs(ne),avs2(ne))
      do l=1,itot
         a=rsplt(l)
         r1=rsi(nl)
         z1=zsi(nl)-zsplt(l)
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
!**********************************************************************
      subroutine v2coef(coef,mp,ne)
      use coilsp
      use mprobe
      use utils, only: pi,tmu,br,bz
      implicit none
      integer*4, intent(in) :: mp,ne
      real*8, intent(out) :: coef
      integer*4 l,mmm
      real*8 a,brc,brct,bzc,bzct,cosm,cosms,sinm,sinms,delsx,delsy,r1,z1, &
             xmp20,ymp20
      integer*4, parameter :: isplit=17,itot=isplit*isplit
      real*8, parameter :: fitot=itot,radeg=pi/180.
!
      if (smp2(mp).gt.0.0) then
         sinm=sin(radeg*amp2(mp))
         cosm=cos(radeg*amp2(mp))
         delsx=smp2(mp)/nsmp2*cosm
         delsy=smp2(mp)/nsmp2*sinm
      else
!------------------------------------------------------------------------------
!--      perpendicular probes                                                --
!------------------------------------------------------------------------------
         sinm=sin(radeg*amp2(mp))
         cosm=cos(radeg*amp2(mp))
         sinms=sin(radeg*(amp2(mp)+90.))
         cosms=cos(radeg*(amp2(mp)+90.))
         delsx=abs(smp2(mp))/nsmp2*cosms
         delsy=abs(smp2(mp))/nsmp2*sinms
      endif
      xmp20=xmp2(mp)-(nsmp2-1)/2.*delsx
      ymp20=ymp2(mp)-(nsmp2-1)/2.*delsy
      brct=0
      bzct=0
      call splitc(isplit,rsplt,zsplt, &
                  rvs(ne),zvs(ne),wvs(ne),hvs(ne),avs(ne),avs2(ne))
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
!**********************************************************************
      subroutine vgrid(coef,rgrid,nr,zgrid,nz,ne)
      use coilsp
      use siloop
      use utils, only: pi,tmu,psical
      implicit none
      integer*4, intent(in) :: nr,nz,ne
      real*8, intent(in) :: rgrid(nr),zgrid(nz)
      real*8, intent(out) :: coef
      integer*4 l
      real*8 a,psic,psict,r1,z1
      integer*4, parameter :: isplit=17,itot=isplit*isplit
      real*8, parameter :: fitot=itot,radeg=pi/180.
!
      psict=0
      call splitc(isplit,rsplt,zsplt, &
                  rvs(ne),zvs(ne),wvs(ne),hvs(ne),avs(ne),avs2(ne))
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

      end module vessel
