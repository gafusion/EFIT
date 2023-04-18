      module coilsp

      implicit none
      real*8, dimension(300) :: rsplt,zsplt

      contains
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          divides up a coil cross section into smaller elements   **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**********************************************************************
      subroutine splitc(is,rs,zs,rc,zc,wc,hc,ac,ac2)
      use utils, only: pi
      implicit none
      integer*4, intent(in) :: is
      real*8, intent(in) :: rc,zc,wc,hc,ac,ac2
      real*8, dimension(is*is), intent(out) :: rs,zs
      integer*4 ic,ii,jj
      real*8 hdelt,wdelt,rdelt,zdelt,htot,wtot,side,rr,rstrt,zstrt,zz
      real*8, parameter :: frd=pi/180.
!
      wdelt=wc/is
      hdelt=hc/is
!----------------------------------------------------------------------
!--   rectangle                                                      --
!----------------------------------------------------------------------
      if(ac+ac2.eq.0.) then
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
!----------------------------------------------------------------------
!--   ac .ne. 0                                                      --
!----------------------------------------------------------------------
      elseif(ac.ne.0.) then
          side=tan(frd*ac)*wc
          zdelt=side/is
          rstrt=rc-wc/2.+wdelt/2.
          htot=hc+side
          zstrt=zc-htot/2.+htot/is/2.
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
!----------------------------------------------------------------------
!--   ac2 .ne. 0                                                     --
!----------------------------------------------------------------------
      elseif(ac2.ne.0.) then
          side=hc/tan(frd*ac2)
          rdelt=side/is
          zstrt=zc-hc/2.+hdelt/2.
          wtot=wc+side
          rstrt=rc-wtot/2.+wtot/is/2.
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

      end module coilsp
