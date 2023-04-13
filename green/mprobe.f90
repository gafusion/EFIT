      module mprobe

      use machine, only: magpri
      implicit none
      public

      integer*4 :: nsmp2
      real*8,dimension (:), allocatable :: xmp2,ymp2,amp2,smp2
      ! these could live in the fcoil module but are only used
      ! in the subroutine here
      integer*4, dimension (:), allocatable :: nshiftrz
      real*8, dimension (:), allocatable :: rshift,zshift,pshift
      real*8, dimension (:), allocatable :: pmprobe

      contains
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          m2coef computes the response functions due to           **
!**          the magnetic probes.                                    **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**********************************************************************
      subroutine m2coef(rr,nr,zz,nz,coef,mp,nc)
      use machine, only: nw,nh
      use grid, only: rgrid,zgrid,brgridfc,bzgridfc
      use fcoil
      use coilsp
      use utils, only: pi,tmu,br,bz
      implicit none
      integer*4, intent(in) :: nr,nz,mp,nc
      real*8, intent(in) :: rr(nr),zz(nz)
      real*8, intent(out) :: coef(mp,nc)
      integer*4 ii,jj,l,k,kk,m,mmm
      real*8 a,brc,brct,brtmp,bzc,bzct,bztmp,cfactor,cmp2, &
             cosm,cosms,cospp,sinm,sinms,sinpp, &
             delsx,delsy,pfnow,pmnow,r,r1,rcos,rsin,z,z1,xmp20,ymp20
      integer*4, parameter :: isplit=17,itot=isplit*isplit
      real*8, parameter :: fitot=itot,radeg=pi/180.
!
      do m=1,magpri
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
!--         shifted f-coil and/or probe
!---------------------------------------------------------------
               if (nshiftrz(k).gt.0) then
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
!--                  shifted f-coil and/or probe
!---------------------------------------------------------------
                     if (nshiftrz(k).gt.0) then
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
!--                  shifted magneic probe                    --
!---------------------------------------------------------------
                     if (nshiftrz(k).gt.0) then
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

      end module mprobe
