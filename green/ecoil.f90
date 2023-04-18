      module ecoil

      use machine, only: necoil,nesum
      implicit none

      integer*4 iecoil
      real*8, dimension (:), allocatable :: re,ze,we,he,ecturn
      integer*4, dimension(:), allocatable :: ecid

      public gecoil
      private e1coef,e2coef,egrid

      contains
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          gecoil gets the Green's functions for E coils.          **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**********************************************************************
      subroutine gecoil(rsilec,rmp2ec,gridec,rgrid,mw,zgrid,mh, &
                        rfcec,recec,rsisec)
      use machine, only: nfcoil,nsilop,magpri,necoil,nesum,&
                      nw,nh,nwnh
      use utils, only: pi,flux
      use fcoil
      implicit none
      integer*4, intent(in) :: mw,mh
      real*8, intent(in) :: rgrid(mw),zgrid(mh)
      real*8, intent(out) :: rsilec(nsilop,nesum),rmp2ec(magpri,nesum), &
                             rfcec(nfcoil,nesum),recec(nesum,nesum), &
                             rsisec(nesum),gridec(nwnh,nesum)
      integer*4 i,j,kk,n
      real*8 aaa,bbb,work
      real*8 taf(nfcoil),taf2(nfcoil)
      real*8, parameter :: zetaec = 3.5e-08
!
      rsilec=0.0
      do j=1,nsilop
         do i=1,necoil
            call e1coef(work,j,i)
            rsilec(j,ecid(i))=rsilec(j,ecid(i))+work*ecturn(i)
         enddo 
      enddo
!
      rmp2ec=0.0
      do j=1,magpri
         do i=1,necoil
            call e2coef(work,j,i)
            rmp2ec(j,ecid(i))=rmp2ec(j,ecid(i))+work*ecturn(i)
         enddo 
      enddo
!
      gridec=0.0
      do i=1,nw
         do j=1,nh
            kk=(i-1)*nh+j
            do n=1,necoil
               call egrid(work,rgrid,i,zgrid,j,n)
               gridec(kk,ecid(n))=gridec(kk,ecid(n))+work*ecturn(n)
            enddo 
         enddo
      enddo
!
      aaa=0.0
      bbb=0.0
      rfcec=0.0
      do j=1,nfcoil
         taf(j)=tan(af(j)*pi/180.)
         taf2(j)=tan(af2(j)*pi/180.)
         do i=1,necoil
            call flux(re(i),ze(i),we(i),he(i),aaa,bbb, &
                      rf(j),zf(j),wf(j),hf(j),taf(j),taf2(j),work)
            work=work*0.5/pi
            rfcec(j,ecid(i))=rfcec(j,ecid(i))+work
         enddo 
      enddo
!
      rsisec=0.0
      recec=0.0
      do j=1,necoil
         do i=1,necoil
            call flux(re(i),ze(i),we(i),he(i),aaa,bbb, &
                      re(j),ze(j),we(j),he(j),aaa,bbb,work)
            work=work*0.5/pi
            recec(ecid(j),ecid(i))=recec(ecid(j),ecid(i))+work
         enddo 
         rsisec(ecid(j))=rsisec(ecid(j))+2.*pi*re(j)/we(j)/he(j)*zetaec
      enddo
      return
      end subroutine gecoil
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          e1coef computes the response functions due to           **
!**          the thin flux loops and E coils.                        **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**********************************************************************
      subroutine e1coef(coef,nl,ne)
      use coilsp
      use utils, only: pi,tmu,psical
      use siloop
      implicit none
      integer*4, intent(in) :: nl,ne
      real*8, intent(out) :: coef
      integer*4 l
      real*8 a,aaa,bbb,psic,psict,r1,z1
      integer*4, parameter :: isplit=17,itot=isplit*isplit
      real*8, parameter :: radeg=pi/180.,fitot=itot
!
      psict=0
      aaa=0.0
      bbb=0.0
      call splitc(isplit,rsplt,zsplt, &
                  re(ne),ze(ne),we(ne),he(ne),aaa,bbb)
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
      end subroutine e1coef
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          e2coef computes the response functions due to           **
!**          the magnetic probes and E coils.                        **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**********************************************************************
      subroutine e2coef(coef,mp,ne)
      use coilsp
      use utils, only: pi,tmu,br,bz
      use mprobe
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
!--  perpendicular probes    96/02/04                                        --
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
                  re(ne),ze(ne),we(ne),he(ne),aaa,bbb)
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
!**********************************************************************
      subroutine egrid(coef,rgrid,nr,zgrid,nz,ne)
      use coilsp
      use utils, only: pi,tmu,psical
      use mprobe
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
                  re(ne),ze(ne),we(ne),he(ne),aaa,bbb)
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

      end module ecoil
