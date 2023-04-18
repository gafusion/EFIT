      module siloop

      use machine, only: nsilop
      implicit none

      integer*4 islpfc,isize
      real*8, dimension (:), allocatable :: rsi,zsi,wsi,hsi,&
                                            as,as2

      contains
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
!**********************************************************************
      subroutine gsilop(rr,nr,zz,nz,rspfun,ns,rsi,zsi,wsi, &
                        hsi,as,as2,ndim)
      use machine, only: nsilop,nwnh
      use utils, only: pi,psical
      implicit none
      integer*4, intent(in) :: nr,nz,ns,ndim
      real*8, dimension(ns), intent(in) :: rsi,zsi,wsi,hsi,as,as2
      real*8, intent(in) :: rr(nr),zz(nz)
      real*8, intent(out) :: rspfun(ndim,nwnh)
      integer*4 ih,iw,nc,nd,nk
      real*8 ab,cmut,dhc,drc,dwc,dzc,rsum,rtmp,tas(nsilop),tas2(nsilop), &
             z,z1,z2
      integer*4, parameter :: isplit=17
      real*8, parameter :: tole=1.0e-10
!
      tas(ns) = tan(as(ns)*pi/180.0)
      tas2(ns) = tan(as2(ns)*pi/180.0)
      do nc=1,nr
        do nd=1,nz
          nk=(nc-1)*nz+nd
          rsum = 0.
          dwc = wsi(ns)/isplit
          dhc = hsi(ns)/isplit
          if (as(ns) .eq. 0. .and. as2(ns) .eq. 0.) then
            z = zsi(ns) - .5*hsi(ns) + .5*dhc
            ab = rsi(ns) - .5*wsi(ns) + .5*dwc
            do iw = 1, isplit
              drc = ab + (iw-1)*dwc + iw*tole
              do ih = 1, isplit
                dzc = z + (ih-1)*dhc
                rtmp = psical(drc,rr(nc),zz(nd)-dzc)
                rsum = rsum + rtmp
              enddo
            enddo
          elseif (as(ns) .ne. 0.) then
            z1 = zsi(ns) - tas(ns)*(wsi(ns)-dwc)/2. - .5*hsi(ns) + .5*dhc
            ab = rsi(ns) - .5*wsi(ns) + .5*dwc
            do iw = 1, isplit
              drc = ab + (iw-1)*dwc + iw*tole
              z2 = z1 + (iw-1)*tas(ns)*dwc
              do ih = 1, isplit
                dzc = z2 + (ih-1)*dhc
                rtmp = psical(drc,rr(nc),zz(nd)-dzc)
                rsum = rsum + rtmp
              enddo
            enddo
          else
            do ih = 1,isplit
              dzc = zsi(ns) - .5*hsi(ns) + .5*dhc + dhc*(ih-1)
              do iw = 1,isplit
                drc = rsi(ns) - .5*wsi(ns) - .5*hsi(ns)/tas2(ns)&
                       + .5*dwc + .5*dhc/tas2(ns)&
                       + dhc/tas2(ns)*(ih-1) + dwc*(iw-1)
                rtmp = psical(drc,rr(nc),zz(nd)-dzc)
                rsum = rsum + rtmp
              enddo
            enddo
          endif
          cmut = rsum*2.e-07/(isplit*isplit)
          rspfun(ns,nk)=cmut
        enddo
      enddo
      return
      end subroutine gsilop
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          m1coef computes the response functions due to           **
!**          the thin flux loops.                                    **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**********************************************************************
      subroutine m1coef(rr,zz,nr,nz,coef,nl,nf)
      use machine, only: nsilop
      use fcoil
      use coilsp
      use utils, only: pi,tmu,psical
      implicit none
      integer*4, intent(in) :: nr,nz,nl,nf
      real*8, intent(in) :: rr(nr),zz(nz)
      real*8, intent(out) :: coef(nsilop,nr*nz)
      integer*4 ii,jj,l,kk
      real*8 a,cmp2,psic,psict,r,r1,z,z1
      integer*4, parameter :: isplit=17,itot=isplit*isplit
      real*8, parameter :: radeg=pi/180.,fitot=itot
!
      if (nf.le.0) then
         do ii=1,nr
            do jj=1,nz
               kk=(ii-1)*nz+jj
               a=rr(ii)
               r=rsi(nl)
               z=zsi(nl)-zz(jj)
               cmp2=psical(a,r,z)*tmu
               coef(nl,kk)=cmp2
            enddo 
         enddo 
      else
         psict=0
         call splitc(isplit,rsplt,zsplt, &
                     rf(nf),zf(nf),wf(nf),hf(nf),af(nf),af2(nf))
         do l=1,itot
            a=rsplt(l)
            r1=rsi(nl)
            z1=zsi(nl)-zsplt(l)
            psic=psical(a,r1,z1)*tmu
            psict=psict+psic/fitot
         enddo 
         coef(nl,nf)=psict
      endif
!
      return
      end subroutine m1coef

      end module siloop
