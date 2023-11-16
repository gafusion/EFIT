!**********************************************************************
!>
!!    This subroutine sets the P' and FF' basis funciton parameters
!!    
!!
!**********************************************************************
      subroutine set_basis_params()
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
!---------------------------------------------------------------------
!--   specific choice of current profile                            --
!--       ICPROF=1  no edge current density allowed                 --
!--       ICPROF=2  free edge current density                       --
!--       ICPROF=3  weak edge current density constraint            --
!---------------------------------------------------------------------
      select case (icprof)
      case (1)
        kffcur=2
        kppcur=2
        fcurbd=1.
        pcurbd=1.
        fwtbp=1.
        fwtqa=0.
        qvfit=0.
      case (2)
        kffcur=2
        kppcur=2
        fcurbd=0.
        pcurbd=0.
        fwtbp=0.
        fwtqa=0.
        qvfit=0.
      case (3)
        kffcur=3
        kppcur=2
        fcurbd=0.
        pcurbd=0.
        fwtbp=0.
        fwtqa=0.
        qvfit=0.
        kcalpa=1
        calpa(1,1)=0.1_dp
        calpa(2,1)=0.1_dp
        calpa(3,1)=0.1_dp
        xalpa(1)=0.0
        kcgama=1
        cgama(1,1)=0.1_dp
        cgama(2,1)=0.1_dp
        cgama(3,1)=0.1_dp
        xgama(1)=0.0
      end select
      if(mse_usecer .eq. 1) keecur = 0
      if (mse_usecer .eq. 2 .and. keecur .eq. 0) then
        keecur = 2
        keefnc = 0
        itek = 5
      endif
      if (imagsigma.gt.0) then
        call errctrl_msg('set_basis_params', &
                         'nonlinear mag sigma calculation deprecated')
      endif
!---------------------------------------------------------------------
!--   adjust fit parameters based on basis function selected        --
!---------------------------------------------------------------------
      select case (kfffnc)
      case (3,4)
        kffcur = 4 * (kffknt - 1)
      case (5)
        kffcur = kffcur * (kffknt - 1)
      case (6)
        kffcur = kffknt * 2
      end select
      select case (kppfnc)
      case (3,4)
        kppcur = 4 * (kppknt - 1)
      case (5)
        kppcur = kppcur * (kppknt - 1)
      case (6)
        kppcur = kppknt * 2
      end select
      select case (kwwfnc)
      case (3,4)
        kwwcur = 4 * (kwwknt - 1)
      case (5)
        kwwcur = kwwcur * (kwwknt - 1)
      case (6)
        kwwcur = kwwknt * 2
      end select
      if (keecur.gt.0) then
        select case (keefnc)
        case (3,4)
          keecur = 4 * (keeknt - 1)
        case (5)
          keecur = keecur * (keeknt - 1)
        case (6)
          keecur = keeknt * 2
        end select
      endif
!
      if (kzeroj.eq.1.and.sizeroj(1).lt.0.0) sizeroj(1)=psiwant
      end subroutine set_basis_params

!**********************************************************************
!>
!!    setff sets up the basis functions for ff-ff(1).
!!    
!!
!!    @param ypsi :
!!
!!    @param xpsii :
!!
!**********************************************************************
      subroutine setff(ypsi,xpsii)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      real*8, intent(in) :: ypsi
      real*8, intent(out) :: xpsii(kffcur)

      do i=1,kffcur
        xpsii(i) = bsffin(kfffnc,i,ypsi)
      enddo
      return
      end subroutine setff

!**********************************************************************
!>
!!    setfp sets up the basis functions for ffp.
!!    
!!
!!    @param ypsi :
!!
!!    @param xpsii :
!!
!**********************************************************************
      subroutine setfp(ypsi,xpsii)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      real*8, intent(in) :: ypsi
      real*8, intent(out) :: xpsii(kffcur)

      do i=1,kffcur
        xpsii(i) = bsffel(kfffnc,i,ypsi)
      enddo
      return
      end subroutine setfp

!**********************************************************************
!>
!!    setfpp computes derivative of ffp.
!!    
!!
!!    @param ypsi :
!!
!!    @param xpsii :
!!
!**********************************************************************
      subroutine setfpp(ypsi,xpsii)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      real*8, intent(in) :: ypsi
      real*8, intent(out) :: xpsii(kffcur)

      do i=1,kffcur
        xpsii(i) = bsffpel(kfffnc,i,ypsi)
      enddo
      return
      end subroutine setfpp

!**********************************************************************
!>
!!    setpp sets up the basis functions for pp.
!!    
!!
!!    @param ypsi :
!!
!!    @param xpsii :
!!
!**********************************************************************
      subroutine setpp(ypsi,xpsii)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      real*8, intent(in) :: ypsi
      real*8, intent(out) :: xpsii(kppcur)

      do i=1,kppcur
        xpsii(i) = bsppel(kppfnc,i,ypsi)
      enddo
      return
      end subroutine setpp

!**********************************************************************
!>
!!    setppp computes derivative of pp.
!!    
!!
!!    @param ypsi :
!!
!!    @param xpsii :
!!
!**********************************************************************
      subroutine setppp(ypsi,xpsii)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      real*8, intent(in) :: ypsi
      real*8, intent(out) :: xpsii(kppcur)

      do i=1,kppcur
        xpsii(i) = bspppel(kppfnc,i,ypsi)
      enddo
      return
      end subroutine setppp

!**********************************************************************
!>
!!    setpr sets up the basis functions for p-p(1).
!!    
!!
!!    @param ypsi :
!!
!!    @param xpsii :
!!
!**********************************************************************
      subroutine setpr(ypsi,xpsii)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      real*8, intent(in) :: ypsi
      real*8, intent(out) :: xpsii(kppcur)

      do i=1,kppcur
        xpsii(i) = bsppin(kppfnc,i,ypsi)
      enddo
      return
      end subroutine setpr

!**********************************************************************
!>
!!    setpw sets up the basis functions for pw-pw(1).
!!    
!!
!!    @param ypsi :
!!
!!    @param xpsii :
!!
!**********************************************************************
      subroutine setpw(ypsi,xpsii)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      real*8, intent(in) :: ypsi
      real*8, intent(out) :: xpsii(kwwcur)

      do i=1,kwwcur
        xpsii(i) = bswwin(kwwfnc,i,ypsi)
      enddo
      end subroutine setpw

!**********************************************************************
!>
!!    setpwp sets up the basis functions for pwp.
!!    
!!
!!    @param ypsi :
!!
!!    @param xpsii :
!!
!**********************************************************************
      subroutine setpwp(ypsi,xpsii)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      real*8, intent(in) :: ypsi
      real*8, intent(out) :: xpsii(kwwcur)

      do i=1,kwwcur
        xpsii(i) = bswwel(kwwfnc,i,ypsi)
      enddo
      return
      end subroutine setpwp

!**********************************************************************
!>
!!    setpwpp computes derivative of pwp.
!!    
!!
!!    @param ypsi :
!!
!!    @param xpsii :
!!
!**********************************************************************
      subroutine setpwpp(ypsi,xpsii)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      real*8, intent(in) :: ypsi
      real*8, intent(out) :: xpsii(kwwcur)

      do i=1,kwwcur
        xpsii(i) = bswwpel(kwwfnc,i,ypsi)
      enddo
      return
      end subroutine setpwpp

!**********************************************************************
!>
!!    seter sets up the basis functions for Er.
!!    
!!
!!    @param ypsi :
!!
!!    @param xpsii :
!!
!**********************************************************************
      subroutine seter(ypsi,xpsii)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      real*8, intent(in) :: ypsi
      real*8, intent(out) :: xpsii(keecur)

      do i=1,keecur
        xpsii(i) = bserel(keefnc,i,ypsi)
      enddo
      return
      end subroutine seter

!**********************************************************************
!>
!!    seterp computes derivative of er.
!!    
!!
!!    @param ypsi :
!!
!!    @param xpsii :
!!
!**********************************************************************
      subroutine seterp(ypsi,xpsii)
      include 'eparm.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      real*8, intent(in) :: ypsi
      real*8, intent(out) :: xpsii(kffcur)

      do i=1,kffcur
        xpsii(i) = bserpel(keefnc,i,ypsi)
      enddo
      return
      end subroutine seterp
