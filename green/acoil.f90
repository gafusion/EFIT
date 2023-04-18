      module acoil

      use machine, only: nacoil
      implicit none

      integer*4 iacoil
      real*8, dimension (:), allocatable :: racoil,zacoil,wacoil,hacoil

      public gacoil
      private a1coef,a2coef,agrid

      contains
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          gacoil gets the Green's functions for the advance       **
!**          divertor coil.                                          **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**********************************************************************
      subroutine gacoil(rsilac,rmp2ac,gridac,rgrid,mw,zgrid,mh)
      use machine, only: nsilop,magpri,nacoil,nw,nh,nwnh
      implicit none
      integer*4, intent(in) :: mw,mh
      real*8, intent(in) :: rgrid(mw),zgrid(mh)
      real*8, intent(out) :: rsilac(nsilop,nacoil),rmp2ac(magpri,nacoil), &
                             gridac(nwnh,nacoil)
      integer*4 i,j,kk,n
      real*8 work
!
      do j=1,nsilop
         do i=1,nacoil
            call a1coef(work,j,i)
            rsilac(j,i)=work
         enddo 
      enddo 
!
      do  j=1,magpri
         do  i=1,nacoil
            call a2coef(work,j,i)
            rmp2ac(j,i)=work
         enddo 
      enddo 
!
      do i=1,nw
         do j=1,nh
            kk=(i-1)*nh+j
            do n=1,nacoil
               call agrid(work,rgrid,i,zgrid,j,n)
               gridac(kk,n)=work
            enddo 
         enddo
      enddo 
      return
      end subroutine gacoil
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          a1coef computes the response functions due to           **
!**          the thin flux loops and the advance divertor coil.      **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**********************************************************************
      subroutine a1coef(coef,nl,ne)
      use coilsp
      use utils, only: pi,tmu,psical
      use siloop
      implicit none
      integer*4, intent(in) :: nl,ne
      real*8, intent(out) :: coef
      integer*4 l
      real*8 a,aaa,bbb,psic,psict,r1,z1
      integer*4, parameter :: isplit=17,itot=isplit*isplit
      real*8, parameter :: fitot=itot,radeg=pi/180.
!
      psict=0
      aaa=0.0
      bbb=0.0
      call splitc(isplit,rsplt,zsplt, &
                  racoil(ne),zacoil(ne),wacoil(ne),hacoil(ne),aaa,bbb)
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
      end subroutine a1coef
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          a2coef computes the response functions due to           **
!**          the magnetic probes and the advance divertor coil.      **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**********************************************************************
      subroutine a2coef(coef,mp,ne)
      use coilsp
      use mprobe
      use utils, only: pi,tmu,br,bz
      implicit none
      integer*4, intent(in) :: mp,ne
      real*8, intent(out) :: coef
      integer*4 l,mmm
      real*8 a,aaa,bbb,brc,brct,bzc,bzct,cosm,cosms,sinm,sinms, &
             delsx,delsy,r1,z1,xmp20,ymp20
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
      aaa=0.0
      bbb=0.0
      call splitc(isplit,rsplt,zsplt, &
                  racoil(ne),zacoil(ne),wacoil(ne),hacoil(ne),aaa,bbb)
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
!**********************************************************************
      subroutine agrid(coef,rgrid,nr,zgrid,nz,ne)
      use coilsp
      use utils, only: pi,tmu,psical
      implicit none
      integer*4, intent(in) :: nr,nz,ne
      real*8, intent(in) :: rgrid(nr),zgrid(nz)
      real*8, intent(out) :: coef
      integer*4 l
      real*8 a,aaa,bbb,psic,psict,r1,z1
      integer*4, parameter :: isplit=17,itot=isplit*isplit
      real*8, parameter :: fitot=itot,radeg=pi/180.
!
      psict=0
      aaa=0.0
      bbb=0.0
      call splitc(isplit,rsplt,zsplt, &
                  racoil(ne),zacoil(ne),wacoil(ne),hacoil(ne),aaa,bbb)
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

      end module acoil
