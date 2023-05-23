!**********************************************************************
!>
!!    This function returns the matrix element for the
!!    selected basis function.
!!    
!!
!!    @param ifunc : basis function number
!!
!!    @param iparm : basis function parameter number
!!
!!    @param ypsi : independent variable value
!!
!********************************************************************** 
      real*8 function bserel(ifunc,iparm,ypsi)
      
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      integer*4, intent(in) :: ifunc, iparm
      real*8, intent(in) :: ypsi
      integer*4 iorder,nk
      real*8 eetens2,tpsi,w
      
      bserel = 0.0
      select case (ifunc)
      case (0)
        if (iparm.eq.1) then
          bserel=1.0 - ypsi**keecur*ecurbd
        else
          bserel = ypsi**(iparm - 1) - ypsi**keecur*ecurbd
        endif
      case (1)
        tpsi = ypsi - 1.0
        if (iparm.eq.1) then
          bserel=1.0 - tpsi**keecur*ecurbd
        else
          bserel = tpsi**(iparm - 1) - tpsi**keecur*ecurbd
        endif
      case (2)
        if (iparm.eq.1) then
          bserel=-(1.0 - ypsi**keecur*ecurbd)
        else
          bserel = -(ypsi**(iparm - 1) - ypsi**keecur*ecurbd)
        endif
      case (3)
        nk = (iparm - 1) / 4 + 1
        if(nk .ge. keeknt) nk = keeknt - 1
        if ((nk .lt. (keeknt - 1) .and. ypsi .lt. eeknt(nk+1) &
          .and.  ypsi .ge. eeknt(nk)) &
          .or. (nk .eq. (keeknt - 1) .and. ypsi .le. eeknt(nk+1) &
          .and.  ypsi .ge. eeknt(nk))) then
          w = eeknt(nk+1) - eeknt(nk)
          if(mod(iparm,4) .eq. 1) bserel = 1.0
          if(mod(iparm,4) .eq. 2) bserel = ypsi
          if(mod(iparm,4) .eq. 3) bserel = cos(w*eetens*ypsi)
          if(mod(iparm,4) .eq. 0) bserel = sin(w*eetens*ypsi)
        endif
      case (4)
        nk = (iparm - 1) / 4 + 1
        if(nk .ge. keeknt) nk = keeknt - 1
        if ((nk .lt. (keeknt - 1) .and. ypsi .lt. eeknt(nk+1) &
          .and.  ypsi .ge. eeknt(nk)) &
          .or. (nk .eq. (keeknt - 1) .and. ypsi .le. eeknt(nk+1) &
          .and.  ypsi .ge. eeknt(nk))) then
          if(mod(iparm,4) .eq. 1) bserel = 1.0
          if(mod(iparm,4) .eq. 2) bserel = ypsi
          if(mod(iparm,4) .eq. 3) bserel = cos(eetens*ypsi)
          if(mod(iparm,4) .eq. 0) bserel = sin(eetens*ypsi)
        endif
      case (5)
        iorder = keecur / (keeknt - 1)
        nk = (iparm - 1) / iorder + 1
        if(nk .ge. keeknt) nk = keeknt - 1
        if ((nk .lt. (keeknt - 1) .and. ypsi .lt. eeknt(nk+1) &
          .and.  ypsi .ge. eeknt(nk)) &
          .or. (nk .eq. (keeknt - 1) .and. ypsi .le. eeknt(nk+1) &
          .and.  ypsi .ge. eeknt(nk))) then
          w = eeknt(nk+1) - eeknt(nk)
          tpsi = (ypsi - eeknt(nk)) / w
          if (mod(iparm,iorder) .eq. 0) then
            if (iorder.eq.1) then
              bserel = 1.0
            else
              bserel = tpsi**(iorder-1)
            endif
          else
            bserel = tpsi**(mod(iparm,iorder)-1)
          endif
        endif
      case (6)
        nk = ((iparm - 1) / 2) + 1
        eetens2 = abs(eetens)*(keeknt-1)/(eeknt(keeknt)-eeknt(1))
        if (nk .gt. 1 ) then
          if (ypsi .le. eeknt(nk) .and. &
            ypsi .ge. eeknt(nk-1)) then
            w = eeknt(nk) - eeknt(nk-1)
            if (mod(iparm,2) .eq. 0) then
              bserel = (sinh(eetens2*(ypsi-eeknt(nk-1)))/ &
                sinh(eetens2*w) - (ypsi-eeknt(nk-1))/w) &
                / (eetens2*eetens2)
            else
              bserel = (ypsi-eeknt(nk-1))/w
            endif
          endif
        endif
        if (nk .lt. keeknt) then
          if (ypsi .ge. eeknt(nk) .and. &
            ypsi .le. eeknt(nk+1)) then
            w = eeknt(nk+1) - eeknt(nk)
            if (mod(iparm,2) .eq. 0) then
              bserel = (sinh(eetens2*(eeknt(nk+1)-ypsi))/ &
                sinh(eetens2*w) - (eeknt(nk+1)-ypsi)/w) &
                / (eetens2*eetens2)
            else
              bserel = (eeknt(nk+1) - ypsi)/w
            endif
          endif
        endif
      case (7)
        if (iparm.eq.keecur) then
          bserel = ypsi**(keehord)
        elseif (iparm .eq. 1) then
          bserel = 1.0
        else
          bserel = ypsi**(iparm - 1)
        endif
      end select

      if (ifunc .ne. keefnc) &
        write(6,*)'ifunc .ne. keefnc ',ifunc,keefnc
      return
      end function bserel

!**********************************************************************
!>
!!    This function returns the matrix element for the
!!    first derivative of selected basis function.
!!    
!!
!!    @param ifunc : basis function number
!!
!!    @param iparm : basis function parameter number
!!
!!    @param ypsi : independent variable value
!!
!**********************************************************************  
      real*8 function bserpel(ifunc,iparm,ypsi)

      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      integer*4, intent(in) :: ifunc, iparm
      real*8, intent(in) :: ypsi
      integer*4 iorder,nk
      real*8 jparm,eetens2,tpsi,w

      bserpel = 0.0
      select case (ifunc) 
      case (0)
        if (iparm.eq.1) then
          if (keecur.eq.1) then
            bserpel = -ecurbd
          else
            bserpel = -keecur*ypsi**(keecur-1)*ecurbd
          endif
        elseif (iparm.eq.2) then
          bserpel = 1. - keecur*ypsi**(keecur-1)*ecurbd
        elseif (iparm.gt.2) then
          bserpel = (iparm - 1)*ypsi**(iparm - 2) - &
            keecur*ypsi**(keecur-1)*ecurbd
        endif
      case (1)
        tpsi = ypsi - 1.0
        if (iparm.eq.1) then
          if (keecur.eq.1) then
            bserpel = -ecurbd
          else
            bserpel = -keecur*tpsi**(keecur-1)*ecurbd
          endif
        elseif (iparm.eq.2) then
          bserpel = 1. - keecur*tpsi**(keecur-1)*ecurbd
        elseif (iparm.gt.2) then
          bserpel = (iparm - 1)*tpsi**(iparm - 2) - &
            keecur*tpsi**(keecur-1)*ecurbd
        endif
      case (2)
        if (iparm.eq.1) then
          if (keecur.eq.1) then
            bserpel = -ecurbd
          else
            bserpel = -keecur*ypsi**(keecur-1)*ecurbd
          endif
        elseif (iparm.eq.2) then
          bserpel = 1. - keecur*ypsi**(keecur-1)*ecurbd
        elseif (iparm.gt.2) then
          bserpel = (iparm - 1)*ypsi**(iparm - 2) - &
            keecur*ypsi**(keecur-1)*ecurbd
        endif
        bserpel = - bserpel
      case (3)
        nk = (iparm - 1) / 4 + 1
        if(nk .ge. keeknt) nk = keeknt - 1
        if ((nk .lt. (keeknt - 1) .and. ypsi .lt. eeknt(nk+1) &
          .and.  ypsi .ge. eeknt(nk)) &
          .or. (nk .eq. (keeknt - 1) .and. ypsi .le. eeknt(nk+1) &
          .and.  ypsi .ge. eeknt(nk))) then
          w = eeknt(nk+1) - eeknt(nk)
          if(mod(iparm,4) .eq. 1) bserpel = 0.0
          if(mod(iparm,4) .eq. 2) bserpel = 1.0
          if(mod(iparm,4) .eq. 3) bserpel = &
            -w*eetens*sin(w*eetens*ypsi)
          if(mod(iparm,4) .eq. 0) bserpel = &
            w*eetens*cos(w*eetens*ypsi)
        endif
      case (4)
        nk = (iparm - 1) / 4 + 1
        if(nk .ge. keeknt) nk = keeknt - 1
        if ((nk .lt. (keeknt - 1) .and. ypsi .lt. eeknt(nk+1) &
          .and.  ypsi .ge. eeknt(nk)) &
          .or. (nk .eq. (keeknt - 1) .and. ypsi .le. eeknt(nk+1) &
          .and.  ypsi .ge. eeknt(nk))) then
          if(mod(iparm,4) .eq. 1) bserpel = 0.0
          if(mod(iparm,4) .eq. 2) bserpel = 1.0
          if(mod(iparm,4) .eq. 3) bserpel = -eetens*sin(eetens*ypsi)
          if(mod(iparm,4) .eq. 0) bserpel = eetens*cos(eetens*ypsi)
        endif
      case (5)
        iorder = keecur / (keeknt - 1)
        nk = (iparm - 1) / iorder + 1
        if(nk .ge. keeknt) nk = keeknt - 1
        if ((nk .lt. (keeknt - 1) .and. ypsi .lt. eeknt(nk+1) &
          .and.  ypsi .ge. eeknt(nk)) &
          .or. (nk .eq. (keeknt - 1) .and. ypsi .le. eeknt(nk+1) &
          .and.  ypsi .ge. eeknt(nk))) then
          w = eeknt(nk+1) - eeknt(nk)
          tpsi = (ypsi - eeknt(nk)) / w
          jparm = mod(iparm,iorder)
          if (jparm.eq.1) then
            bserpel = 0.0
          elseif (jparm.eq.2) then
            bserpel = 1./w
          elseif (jparm.gt.2) then
            bserpel = (jparm - 1)/w*tpsi**(jparm - 2)
          endif
        endif
      case (6)
        nk = ((iparm - 1) / 2) + 1
        eetens2 = abs(eetens)*(keeknt-1)/(eeknt(keeknt)-eeknt(1))
        if (nk .gt. 1) then
          if (ypsi .le. eeknt(nk) .and. &
            ypsi .ge. eeknt(nk-1)) then
            w = eeknt(nk) - eeknt(nk-1)
            if (mod(iparm,2) .eq. 0) then
              bserpel = (eetens2*cosh(eetens2* &
                (ypsi-eeknt(nk-1)))/sinh(eetens2*w) - (1.0/w)) &
                / (eetens2*eetens2)
            else
              bserpel = 1.0/w
            endif
          endif
        endif
        if (nk .lt. keeknt) then
          if (ypsi .ge. eeknt(nk) .and. &
            ypsi .le. eeknt(nk+1)) then
            w = eeknt(nk+1) - eeknt(nk)
            if (mod(iparm,2) .eq. 0) then
              bserpel = (-eetens2*cosh(eetens2* &
                (eeknt(nk+1)-ypsi))/sinh(eetens2*w)+(1.0/w)) &
                / (eetens2*eetens2)
            else
              bserpel = -1.0/w
            endif
          endif
        endif
      case (7)
        if (iparm.eq.keecur) then
          bserpel = keehord*ypsi**(keehord-1)
        elseif (iparm.eq.1 ) then
          bserpel = 0.
        elseif (iparm.eq.2 ) then
          bserpel = 1.
        elseif (iparm.gt.2) then
          bserpel = (iparm - 1)*ypsi**(iparm - 2)
        endif
      end select
      return
      end function bserpel

!**********************************************************************
!>
!!    This function returns the matrix element for the
!!     selected basis function.
!!    
!!
!!    @param ifunc : basis function number
!!
!!    @param iparm : basis function parameter number
!!
!!    @param ypsi : independent variable value
!!
!**********************************************************************        
      real*8 function bserin(ifunc,iparm,ypsi)
      
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      integer*4, intent(in) :: ifunc, iparm
      real*8, intent(in) :: ypsi
      integer*4 iorder,nk
      real*8 b1,b2,eetens2,tpsi,tpsi2,w,ypsi1,ypsi2
      
      bserin = 0.0
      ypsi2 = 1.0
      select case (ifunc)
      case (0)
        bserin = (ypsi**iparm)/iparm - &
          (ypsi**(keecur+1))/(keecur+1)*ecurbd
        bserin = bserin - ((ypsi2**iparm)/iparm - &
          (ypsi2**(keecur+1))/(keecur+1)*ecurbd)
      case (1)
        tpsi = ypsi - 1.0
        tpsi2 = ypsi2 - 1.0
        bserin = (tpsi**iparm)/iparm - &
          (tpsi**(keecur+1))/(keecur+1)*ecurbd
        bserin = bserin - ((tpsi2**iparm)/iparm - &
          (tpsi2**(keecur+1))/(keecur+1)*ecurbd)
      case (2)
        bserin = -((ypsi**iparm)/iparm - &
          (ypsi**(keecur+1))/(keecur+1)*ecurbd)
        bserin = bserin - (-((ypsi2**iparm)/iparm - &
          (ypsi2**(keecur+1))/(keecur+1)*ecurbd))
      case (3)
        nk = (iparm - 1) / 4 + 1
        if(nk .ge. keeknt) nk = keeknt - 1
        if (ypsi .ge. eeknt(nk+1)) then
          bserin = 0
          return
        endif
        if (1.0 .le. eeknt(nk)) then
          bserin = 0
          return
        endif
        if (ypsi .ge. eeknt(nk)) then
          ypsi1 = ypsi
        else
          ypsi1 = eeknt(nk)
        endif
        w = eeknt(nk+1) - eeknt(nk)
        if (1.0 .ge. eeknt(nk+1)) then
          ypsi2 = eeknt(nk+1)
        else
          ypsi2 = 1.0
        endif
        if(mod(iparm,4) .eq. 1) b1 = ypsi1
        if(mod(iparm,4) .eq. 2) b1 = (ypsi1**2) / 2.0
        if(mod(iparm,4) .eq. 3) &
          b1 = sin(w*eetens*ypsi1)/w*eetens
        if(mod(iparm,4) .eq. 0) &
          b1 = -cos(w*eetens*ypsi1)/w*eetens
        if(mod(iparm,4) .eq. 1) b2 = ypsi2
        if(mod(iparm,4) .eq. 2) b2 = (ypsi2**2) / 2.0
        if(mod(iparm,4) .eq. 3) &
          b2 = sin(w*eetens*ypsi2)/w*eetens
        if(mod(iparm,4) .eq. 0) &
          b2 = -cos(w*eetens*ypsi2)/w*eetens
        bserin = b1 - b2
        !     write(6,*)'for ypsi=',ypsi,' integrate from ',ypsi1,' to ',
        !     $                       ypsi2,' = ',bserin
      case (4)
        nk = (iparm - 1) / 4 + 1
        if(nk .ge. keeknt) nk = keeknt - 1
        if (ypsi .ge. eeknt(nk+1)) then
          bserin = 0
          return
        endif
        if (1.0 .le. eeknt(nk)) then
          bserin = 0
          return
        endif
        if (ypsi .ge. eeknt(nk)) then
          ypsi1 = ypsi
        else
          ypsi1 = eeknt(nk)
        endif
        if (1.0 .ge. eeknt(nk+1)) then
          ypsi2 = eeknt(nk+1)
        else
          ypsi2 = 1.0
        endif
        if(mod(iparm,4) .eq. 1) b1 = ypsi1
        if(mod(iparm,4) .eq. 2) b1 = (ypsi1**2) / 2.0
        if(mod(iparm,4) .eq. 3) &
          b1 = sin(eetens*ypsi1)/eetens
        if(mod(iparm,4) .eq. 0) &
          b1 = -cos(eetens*ypsi1)/eetens
        if(mod(iparm,4) .eq. 1) b2 = ypsi2
        if(mod(iparm,4) .eq. 2) b2 = (ypsi2**2) / 2.0
        if(mod(iparm,4) .eq. 3) &
          b2 = sin(eetens*ypsi2)/eetens
        if(mod(iparm,4) .eq. 0) &
          b2 = -cos(eetens*ypsi2)/eetens
        bserin = b1 - b2
        !     write(6,*)'for ypsi=',ypsi,' integrate from ',ypsi1,' to ',
        !     $                       ypsi2,' = ',bserin
      case (5)
        iorder = keecur / (keeknt - 1)
        nk = (iparm - 1) / iorder + 1
        if(nk .ge. keeknt) nk = keeknt - 1
        if (ypsi .ge. eeknt(nk+1)) then
          bserin = 0
          return
        endif
        if (1.0 .le. eeknt(nk)) then
          bserin = 0
          return
        endif
        if (ypsi .ge. eeknt(nk)) then
          ypsi1 = ypsi
        else
          ypsi1 = eeknt(nk)
        endif
        if (1.0 .ge. eeknt(nk+1)) then
          ypsi2 = eeknt(nk+1)
        else
          ypsi2 = 1.0
        endif
        w = eeknt(nk+1) - eeknt(nk)
         
        tpsi=(ypsi1**2/2.0 - ypsi1*eeknt(nk))/2
        if (mod(iparm,iorder) .eq. 0) then
          b1 = tpsi**iorder / iorder
        else
          b1 = tpsi**mod(iparm,iorder) / mod(iparm,iorder)
        endif
        tpsi=(ypsi2**2/2.0 - ypsi2*eeknt(nk))/2
        if (mod(iparm,iorder) .eq. 0) then
          b2 = tpsi**iorder / iorder
        else
          b2 = tpsi**mod(iparm,iorder) / mod(iparm,iorder)
        endif
        bserin = b1 - b2
        !     write(6,*)'for ypsi=',ypsi,' integrate from ',ypsi1,' to ',
        !     $                       ypsi2,' = ',bserin
      case (6)
        nk = ((iparm - 1) / 2) + 1
        eetens2 = abs(eetens)*(keeknt-1)/(eeknt(keeknt)-eeknt(1))
        bserin = 0
        if (nk .gt.1) then
          if (ypsi .le. eeknt(nk)) then
            if (ypsi .le. eeknt(nk-1)) then
              ypsi1 = eeknt(nk-1)
            else
              ypsi1 = ypsi
            endif
            if (1.0 .le. eeknt(nk)) then
              ypsi2 = 1.0
            else
              ypsi2 = eeknt(nk)
            endif
            w = eeknt(nk) - eeknt(nk-1)
            if (mod(iparm,2) .eq. 0) then
              b1 = (cosh(eetens2*(ypsi1-eeknt(nk-1)))/ &
                (eetens2*sinh(eetens2*w)) - (ypsi1*ypsi1/2.0- &
                eeknt(nk-1)*ypsi1)/w) &
                / (eetens2*eetens2)
            else
              b1 = (ypsi1*ypsi1/2.0-eeknt(nk-1)*ypsi1)/w
            endif
            if (mod(iparm,2) .eq. 0) then
              b2 = (cosh(eetens2*(ypsi2-eeknt(nk-1)))/ &
                (eetens2*sinh(eetens2*w)) - (ypsi2*ypsi2/2.0- &
                eeknt(nk-1)*ypsi2)/w) &
                / (eetens2*eetens2)
            else
              b2 = (ypsi2*ypsi2/2.0-eeknt(nk-1)*ypsi2)/w
            endif
            bserin = bserin + b1 - b2
          endif
        endif
        if (nk .lt. keeknt) then
          if (ypsi .le. eeknt(nk+1)) then
            if (ypsi .le. eeknt(nk)) then
              ypsi1 = eeknt(nk)
            else
              ypsi1 = ypsi
            endif
            if (1.0 .le. eeknt(nk+1)) then
              ypsi2 = 1.0
            else
              ypsi2 = eeknt(nk+1)
            endif
            w = eeknt(nk+1) - eeknt(nk)
            if (mod(iparm,2) .eq. 0) then
              b1 = (-cosh(eetens2*(eeknt(nk+1)-ypsi1))/ &
                (eetens2*sinh(eetens2*w))-(ypsi1*eeknt(nk+1) &
                -ypsi1*ypsi1/2.0)/w) &
                / (eetens2*eetens2)
            else
              b1 = (eeknt(nk+1)*ypsi1-ypsi1*ypsi1/2.0)/w
            endif
            if (mod(iparm,2) .eq. 0) then
              b2 = (-cosh(eetens2*(eeknt(nk+1)-ypsi2))/ &
                (eetens2*sinh(eetens2*w))-(ypsi2*eeknt(nk+1) &
                -ypsi2*ypsi2/2.0)/w) &
                / (eetens2*eetens2)
            else
              b2 = (eeknt(nk+1)*ypsi2-ypsi2*ypsi2/2.0)/w
            endif
            bserin = bserin + b1 - b2
          endif
        endif
      case (7)
        if (iparm .eq. keecur) then
          bserin = (ypsi**(keehord+1))/(keehord+1)
          bserin = bserin - ((ypsi2**(keehord+1))/(keehord+1))
        else
          bserin = (ypsi**iparm)/iparm
          bserin = bserin - ((ypsi2**iparm)/iparm)
        endif
      end select

      if (ifunc .ne. keefnc) &
        write(6,*) 'ifunc .ne. keefnc ',ifunc,keefnc
      return
      end function bserin 
       
!**********************************************************************
!>
!!    In addition to the least squares constraints that
!!    efit already uses, some basis functions have exact
!!    constraints. Most notable is the spline function
!!    whose continuity constraints are exact, not LSE.
!!    
!!
!!    @param ncrsp : number of constraint equations
!!
!!    @param crsp : constraint matrix
!!
!!    @param z : value vector
!!
!!    @param nffcoi : array index for setting up crsp
!!
!**********************************************************************      
      subroutine eecnst(ncrsp,crsp,z,nffcoi)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      real*8 bserel
      integer*4, intent(in) :: nffcoi
      integer*4, intent(inout) :: ncrsp
      real*8, intent(out) :: crsp(4*(npcurn-2)+6+npcurn*npcurn,nrsmat), &
                             z(4*(npcurn-2)+6+npcurn*npcurn)
      integer*4 i,j,iorder
      real*8 h,w,eetens2

      select case (keefnc)
      case (3)
        if (keeknt .gt. 2) then
!     
!     first set of constraints is that splines must be equal at the knots
!     
          do i = 2,keeknt-1
            ncrsp = ncrsp + 1
            crsp(ncrsp,:) = 0.0
            z(ncrsp) = 0.0
            h = eeknt(i) - eeknt(i-1)
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 1) = 1.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 2) = eeknt(i)
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 3) = &
              cos(h * eetens * eeknt(i))
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 4) =  &
              sin(h * eetens * eeknt(i))
               
            h = eeknt(i+1) - eeknt(i)
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 5) = -1.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 6) = -eeknt(i)
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 7) = &
              -cos(h * eetens * eeknt(i))
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 8) =  &
              -sin(h * eetens * eeknt(i))
          enddo
          !
          !     second set of constraints is that splines have equal first
          !     derivative at the knots
          !
          do i = 2,keeknt-1
            ncrsp = ncrsp + 1
            crsp(ncrsp,:) = 0.0
            z(ncrsp) = 0.0
            h = eeknt(i) - eeknt(i-1)
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 1) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 2) = 1.0
            crsp(ncrsp,nffcoi + kppcur + kffcur + &
              4*(i-2) + 3) = &
              -h * eetens * sin(h * eetens * eeknt(i))
            crsp(ncrsp,nffcoi + kppcur + kffcur + &
              4*(i-2) + 4) =  &
              h * eetens * cos(h * eetens * eeknt(i))
               
            h = eeknt(i+1) - eeknt(i)
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 5) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 6) = -1.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 7) = &
              h * eetens * sin(h * eetens * eeknt(i))
            crsp(ncrsp,nffcoi + kppcur + kffcur + 4*(i-2) + 8) = &
              -h * eetens * cos(h * eetens * eeknt(i))
          enddo
          !
          !     second set of constraints is that splines have equal second
          !     derivative at the knots
          !
          do i = 2,keeknt-1
            ncrsp = ncrsp + 1
            crsp(ncrsp,:) = 0.0
            z(ncrsp) = 0.0
            h = eeknt(i) - eeknt(i-1)
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 1) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 2) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 3) = &
              -h*h*eetens*eetens*cos(h * eetens * eeknt(i))
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 4) =  &
              -h*h*eetens*eetens*sin(h * eetens * eeknt(i))
               
            h = eeknt(i+1) - eeknt(i)
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 5) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 6) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 7) = &
              h*h*eetens*eetens*cos(h * eetens * eeknt(i))
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 8) =  &
              h*h*eetens*eetens*sin(h * eetens * eeknt(i))
          enddo
        endif
      case (4)
        if (keeknt .le. 2) then
          !
          !     first set of constraints is that splines must be equal at the knots
          !
          do i = 2,keeknt-1
            ncrsp = ncrsp + 1
            crsp(ncrsp,:) = 0.0
            z(ncrsp) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 1) = 1.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 2) = eeknt(i)
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 3) = &
              cos(eetens * eeknt(i))
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 4) =  &
              sin(eetens * eeknt(i))
               
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 5) = -1.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 6) = -eeknt(i)
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 7) = &
              -cos(eetens * eeknt(i))
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 8) =  &
              -sin(eetens * eeknt(i))
          enddo
          !
          !     second set of constraints is that splines have equal first
          !     derivative at the knots
          !
          do i = 2,keeknt-1
            ncrsp = ncrsp + 1
            crsp(ncrsp,:) = 0.0
            z(ncrsp) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 1) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 2) = 1.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 3) = &
              -eetens * sin(eetens * eeknt(i))
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 4) =  &
              eetens * cos(eetens * eeknt(i))
               
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 5) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 6) = -1.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 7) = &
              eetens * sin(eetens * eeknt(i))
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 8) =  &
              -eetens * cos(eetens * eeknt(i))
          enddo
          !
          !     second set of constraints is that splines have equal second
          !     derivative at the knots
          !
          do i = 2,keeknt-1
            ncrsp = ncrsp + 1
            crsp(ncrsp,:) = 0.0
            z(ncrsp) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 1) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 2) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 3) = &
              -eetens*eetens*cos(eetens * eeknt(i))
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 4) =  &
              -eetens*eetens*sin(eetens * eeknt(i))
               
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 5) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 6) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 7) = &
              eetens*eetens*cos(eetens * eeknt(i))
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 8) =  &
              eetens*eetens*sin(eetens * eeknt(i))
          enddo
        endif
      case (5)
        iorder = keecur / (keeknt - 1)
        if (keeknt .le. 2) then
          !
          !     first set of constraints is that splines must be equal at the knots
          !
          do i = 2,keeknt-1
            ncrsp = ncrsp + 1
            crsp(ncrsp,:) = 0.0
            z(ncrsp) = 0.0
            do j= 1,iorder
              crsp(ncrsp,nffcoi + kppcur + kffcur &
                + iorder*(i-2) + j)  = 1.0
            enddo
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + iorder*(i-1) + 1)  = -1.0
          enddo
          !
          !     second set of constraints is that splines have equal first
          !     derivative at the knots
          !
          do i = 2,keeknt-1
            ncrsp = ncrsp + 1
            crsp(ncrsp,:) = 0.0
            z(ncrsp) = 0.0
            do j= 2,iorder
              crsp(ncrsp,nffcoi + kppcur + kffcur &
                + iorder*(i-2) + j)  = (j-1)
            enddo
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + iorder*(i-1) + 2)  = -1.0
          enddo
          !
          !     second set of constraints is that splines have equal second
          !     derivative at the knots
          !
          do i = 2,keeknt-1
            ncrsp = ncrsp + 1
            crsp(ncrsp,:) = 0.0
            z(ncrsp) = 0.0
            do j= 3,iorder
              crsp(ncrsp,nffcoi + kppcur + kffcur &
                + iorder*(i-2) + j)  = (j-1)*(j-2)
            enddo
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + iorder*(i-1) + 3)  = -2.0
          enddo
        endif
      case (6)
        !
        !     first set of constraints is that splines have equal first
        !     derivative at the knots
        !
        if (keeknt .gt. 2) then
          eetens2 = abs(eetens)*(keeknt-1)/(eeknt(keeknt)-eeknt(1))
          do i = 2,keeknt-1
            ncrsp = ncrsp + 1
            crsp(ncrsp,:) = 0.0
            z(ncrsp) = 0.0
            w = eeknt(i+1) - eeknt(i)
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 2*(i-1) + 1) = -1.0/w
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 2*(i-1) + 2) = (-eetens2* &
              cosh(eetens2*w)/sinh(eetens2*w) &
              + 1.0/w)/(eetens2*eetens2)
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 2*(i-1) + 3) = 1.0/w
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 2*(i-1) + 4) = (eetens2/ &
              sinh(eetens2*w) - 1.0/w)/(eetens2*eetens2)
               
            w = eeknt(i) - eeknt(i-1)
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 2*(i-1) - 1) = 1.0/w
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 2*(i-1) + 0) = -(-eetens2/ &
              sinh(eetens2*w) + 1.0/w)/(eetens2*eetens2)
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 2*(i-1) + 1) =  &
              crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 2*(i-1) + 1) - 1.0/w
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 2*(i-1) + 2) =  &
              crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 2*(i-1) + 2) - (eetens2* &
              cosh(eetens2*w)/sinh(eetens2*w) &
              - 1.0/w)/(eetens2*eetens2)
          enddo
        endif
        do i = 1,keeknt
          if (keebdry(i) .eq. 1) then
            ncrsp = ncrsp + 1
            crsp(ncrsp,:) = 0.0
            z(ncrsp) = eebdry(i)
            crsp(ncrsp,nffcoi + kppcur + kffcur+2*i - 1) = 1.0
          endif
          if (kee2bdry(i) .eq. 1) then
            ncrsp = ncrsp + 1
            crsp(ncrsp,:) = 0.0
            z(ncrsp) = ee2bdry(i)
            crsp(ncrsp,nffcoi + kppcur + kffcur+2*i) = 1.0
          endif
        enddo
      end select
      select case (keefnc)
      case (3,4,5,6)
        if (ecurbd .eq. 1.0) then
          ncrsp = ncrsp + 1
          crsp(ncrsp,:) = 0.0
          z(ncrsp) = 0.0
          do j = 1,keecur
            crsp(ncrsp,nffcoi+kppcur+kffcur+j) &
              = bserel(keefnc,j,1.0)
          enddo
        endif
      end select
      return
      end subroutine eecnst

!**********************************************************************
!>
!!    Store the solution coefs into eebdry and ee2bdry
!!    
!!
!**********************************************************************
      subroutine eestore()
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      integer*4 i

      if (keefnc .gt. 0 .and. keefnc .le. 2) then
         do i = 1,keecur
           eebdry(i) = cerer(i)
           ee2bdry(i) = 0.0
         enddo
      else if (keefnc .eq. 6)then
         do i = 1,keeknt
           if(keebdry(i) .ne. 1)   eebdry(i) = cerer(2*i - 1)
           if(kee2bdry(i) .ne. 1) ee2bdry(i) = cerer(2*i)
         enddo
      endif
      return
      end subroutine eestore
