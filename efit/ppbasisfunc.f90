!**********************************************************************
!>
!!    This function returns the matrix element for the
!!    selected basis function.
!!
!!    @param ifunc : basis function number
!!
!!    @param iparm : basis function parameter number
!!
!!    @param ypsi : independent variable value
!!
!**********************************************************************
  real*8 function bsppel(ifunc,iparm,ypsi)
      
  include 'eparm.inc'
  include 'modules2.inc'
  include 'modules1.inc'
  implicit none
  integer*4, intent(in) :: ifunc, iparm
  real*8, intent(in) :: ypsi
  integer*4 iorder,nk
  real*8 tpsi,pptens2,w
      
  bsppel = 0.0
  select case (ifunc)
  case (0)
    if (iparm.eq.1) then
      bsppel = 1.0 - ypsi**kppcur*pcurbd
    else
      bsppel = ypsi**(iparm - 1) - ypsi**kppcur*pcurbd
    endif
  case (1)
    tpsi = ypsi - 1.0
    if (iparm.eq.1) then
      bsppel = 1.0 - tpsi**kppcur*pcurbd
    else
      bsppel = tpsi**(iparm - 1) - tpsi**kppcur*pcurbd
    endif
  case (2)
    if (iparm.eq.1) then
      bsppel = -(1.0 - ypsi**kppcur*pcurbd)
    else
      bsppel = -(ypsi**(iparm - 1) - ypsi**kppcur*pcurbd)
    endif
  case (3)
    nk = (iparm - 1) / 4 + 1
    if(nk .ge. kppknt) nk = kppknt - 1
    if ((nk .lt. (kppknt - 1) .and. ypsi .lt. ppknt(nk+1) &
      .and.  ypsi .ge. ppknt(nk)) &
      .or. (nk .eq. (kppknt - 1) .and. ypsi .le. ppknt(nk+1) &
      .and.  ypsi .ge. ppknt(nk))) then
      w = ppknt(nk+1) - ppknt(nk)
      if(mod(iparm,4) .eq. 1) bsppel = 1.0
      if(mod(iparm,4) .eq. 2) bsppel = ypsi
      if(mod(iparm,4) .eq. 3) bsppel = cos(w*pptens*ypsi)
      if(mod(iparm,4) .eq. 0) bsppel = sin(w*pptens*ypsi)
    endif
  case (4)
    nk = (iparm - 1) / 4 + 1
    if(nk .ge. kppknt) nk = kppknt - 1
    if ((nk .lt. (kppknt - 1) .and. ypsi .lt. ppknt(nk+1) &
      .and.  ypsi .ge. ppknt(nk)) &
      .or. (nk .eq. (kppknt - 1) .and. ypsi .le. ppknt(nk+1) &
      .and.  ypsi .ge. ppknt(nk))) then
      if(mod(iparm,4) .eq. 1) bsppel = 1.0
      if(mod(iparm,4) .eq. 2) bsppel = ypsi
      if(mod(iparm,4) .eq. 3) bsppel = cos(pptens*ypsi)
      if(mod(iparm,4) .eq. 0) bsppel = sin(pptens*ypsi)
    endif
  case (5)
    iorder = kppcur / (kppknt - 1)
    nk = (iparm - 1) / iorder + 1
    if(nk .ge. kppknt)nk = kppknt - 1
    if ((nk .lt. (kppknt - 1) .and. ypsi .lt. ppknt(nk+1) &
      .and.  ypsi .ge. ppknt(nk)) &
      .or. (nk .eq. (kppknt - 1) .and. ypsi .le. ppknt(nk+1) &
      .and.  ypsi .ge. ppknt(nk))) then
      w = ppknt(nk+1) - ppknt(nk)
      tpsi = (ypsi - ppknt(nk)) / w
      if (mod(iparm,iorder) .eq. 0) then
        if (iorder.eq.1) then
          bsppel=1.0
        else
          bsppel = tpsi**(iorder-1)
        endif
      else
        bsppel = tpsi**(mod(iparm,iorder)-1)
      endif
    endif
  case (6)
    nk = ((iparm - 1) / 2) + 1
    pptens2 = abs(pptens)*(kppknt-1)/(ppknt(kppknt)-ppknt(1))
    if (nk .gt. 1) then
      if (ypsi .le. ppknt(nk) .and.  &
        ypsi .ge. ppknt(nk-1)) then
        w = ppknt(nk) - ppknt(nk-1)
        if (mod(iparm,2) .eq. 0) then
          bsppel = (sinh(pptens2*(ypsi-ppknt(nk-1)))/ &
            sinh(pptens2*w) - (ypsi-ppknt(nk-1))/w) &
            / (pptens2*pptens2)
        else
          bsppel = (ypsi-ppknt(nk-1))/w
        endif
      endif
    endif
    if (nk .lt. kppknt) then
      if (ypsi .ge. ppknt(nk) .and.  &
        ypsi .le. ppknt(nk+1)) then
        w = ppknt(nk+1) - ppknt(nk)
        if(mod(iparm,2) .eq. 0) then
          bsppel = (sinh(pptens2*(ppknt(nk+1)-ypsi))/ &
            sinh(pptens2*w) - (ppknt(nk+1)-ypsi)/w) &
            / (pptens2*pptens2)
        else
          bsppel = (ppknt(nk+1) - ypsi)/w
        endif 
      endif
    endif
  case (7)
    if (iparm.eq.kppcur) then
      bsppel = ypsi**(kpphord)
    elseif (iparm .eq. 1) then
      bsppel = 1.0
    else
      bsppel = ypsi**(iparm - 1)
    endif
  end select

  if (ifunc .ne. kppfnc)  &
    write(6,*)'ifunc .ne. kppfnc ',ifunc,kppfnc
  return
  end function bsppel

!**********************************************************************
!>    
!!    Function bspppel(ifunc,iparm,ypsi)
!!    
!!    This function returns the matrix element for the
!!    first derivative of selected basis function.
!!    
!!    @param ifunc : basis function number
!!
!!    @param iparm : basis function parameter number
!!
!!    @param ypsi : independent variable value
!!    
!**********************************************************************
  real*8 function bspppel(ifunc,iparm,ypsi)

  include 'eparm.inc'
  include 'modules2.inc'
  include 'modules1.inc'
  implicit none
  integer*4, intent(in) :: ifunc, iparm
  real*8, intent(in) :: ypsi
  integer*4 iorder,nk
  real*8 jparm,pptens2,tpsi,w

  bspppel = 0.0
  select case (ifunc)
  case (0)
    if (iparm.eq.1) then
      if (kppcur.eq.1) then
        bspppel = -pcurbd
      else
        bspppel = -kppcur*ypsi**(kppcur-1)*pcurbd
      endif
    elseif (iparm.eq.2) then
      bspppel = 1. - kppcur*ypsi**(kppcur-1)*pcurbd
    elseif (iparm.gt.2) then
      bspppel = (iparm - 1)*ypsi**(iparm - 2) - &
        kppcur*ypsi**(kppcur-1)*pcurbd
    endif
  case (1)
    tpsi = ypsi - 1.0
    if (iparm.eq.1) then
      if (kppcur.eq.1) then
        bspppel = -pcurbd
      else
        bspppel = -kppcur*tpsi**(kppcur-1)*pcurbd
      endif
    elseif (iparm.eq.2) then
      bspppel = 1. - kppcur*tpsi**(kppcur-1)*pcurbd
    elseif (iparm.gt.2) then
      bspppel = (iparm - 1)*tpsi**(iparm - 2) - &
        kppcur*tpsi**(kppcur-1)*pcurbd
    endif
  case (2)
    if (iparm.eq.1) then
      if (kppcur.eq.1) then
        bspppel = -pcurbd
      else
        bspppel = -kppcur*ypsi**(kppcur-1)*pcurbd
      endif
    elseif (iparm.eq.2) then
      bspppel = 1. - kppcur*ypsi**(kppcur-1)*pcurbd
    elseif (iparm.gt.2) then
      bspppel = (iparm - 1)*ypsi**(iparm - 2) - &
        kppcur*ypsi**(kppcur-1)*pcurbd
    endif
    bspppel = - bspppel
  case (3)
    nk = (iparm - 1) / 4 + 1
    if(nk .ge. kppknt) nk = kppknt - 1
    if ((nk .lt. (kppknt - 1) .and. ypsi .lt. ppknt(nk+1) &
      .and.  ypsi .ge. ppknt(nk)) &
      .or. (nk .eq. (kppknt - 1) .and. ypsi .le. ppknt(nk+1) &
      .and.  ypsi .ge. ppknt(nk))) then
      w = ppknt(nk+1) - ppknt(nk)
      if(mod(iparm,4) .eq. 1) bspppel = 0.0
      if(mod(iparm,4) .eq. 2) bspppel = 1.0
      if(mod(iparm,4) .eq. 3) bspppel =  &
        -w*pptens*sin(w*pptens*ypsi)
      if(mod(iparm,4) .eq. 0) bspppel =  &
        w*pptens*cos(w*pptens*ypsi)
    endif
  case (4)
    nk = (iparm - 1) / 4 + 1
    if(nk .ge. kppknt) nk = kppknt - 1
    if ((nk .lt. (kppknt - 1) .and. ypsi .lt. ppknt(nk+1) &
      .and.  ypsi .ge. ppknt(nk)) &
      .or. (nk .eq. (kppknt - 1) .and. ypsi .le. ppknt(nk+1) &
      .and.  ypsi .ge. ppknt(nk))) then
      if(mod(iparm,4) .eq. 1) bspppel = 0.0
      if(mod(iparm,4) .eq. 2) bspppel = 1.0
      if(mod(iparm,4) .eq. 3) bspppel = -pptens*sin(pptens*ypsi)
      if(mod(iparm,4) .eq. 0) bspppel = pptens*cos(pptens*ypsi)
    endif
  case (5)
    iorder = kppcur / (kppknt - 1)
    nk = (iparm - 1) / iorder + 1
    if(nk .ge. kppknt) nk = kppknt - 1
    if ((nk .lt. (kppknt - 1) .and. ypsi .lt. ppknt(nk+1) &
      .and.  ypsi .ge. ppknt(nk)) &
      .or. (nk .eq. (kppknt - 1) .and. ypsi .le. ppknt(nk+1) &
      .and.  ypsi .ge. ppknt(nk))) then
      w = ppknt(nk+1) - ppknt(nk)
      tpsi = (ypsi - ppknt(nk)) / w
      jparm = mod(iparm,iorder)
      if (jparm.eq.1) then
        bspppel = 0.0
      elseif (jparm.eq.2) then
        bspppel = 1./w
      elseif (jparm.gt.2) then
        bspppel = (jparm - 1)/w*tpsi**(jparm - 2)
      endif
    endif
  case (6)
    nk = ((iparm - 1) / 2) + 1
    pptens2 = abs(pptens)*(kppknt-1)/(ppknt(kppknt)-ppknt(1))
    if (nk .gt. 1) then
      if (ypsi .le. ppknt(nk) .and.  &
        ypsi .ge. ppknt(nk-1)) then
        w = ppknt(nk) - ppknt(nk-1)
        if(mod(iparm,2) .eq. 0) then
          bspppel = (pptens2*cosh(pptens2* &
            (ypsi-ppknt(nk-1)))/sinh(pptens2*w) - (1.0/w)) &
            / (pptens2*pptens2)
        else
          bspppel = 1.0/w
        endif
      endif
    endif
    if (nk .lt. kppknt) then
      if (ypsi .ge. ppknt(nk) .and.  &
        ypsi .le. ppknt(nk+1)) then
        w = ppknt(nk+1) - ppknt(nk)
        if(mod(iparm,2) .eq. 0) then
          bspppel = ((-pptens2*cosh(pptens2* &
            (ppknt(nk+1)-ypsi)))/sinh(pptens2*w) + (1.0/w)) &
            / (pptens2*pptens2)
        else
          bspppel = -1.0/w
        endif
      endif
    endif
  case (7)
    if (iparm.eq.kppcur) then
      bspppel = kpphord*ypsi**(kpphord-1)
    elseif (iparm.eq.1) then
      bspppel = 0.
    elseif (iparm.eq.2) then
      bspppel = 1.
    elseif (iparm.gt.2) then
      bspppel = (iparm - 1)*ypsi**(iparm - 2)
    endif
  end select
  return
  end function bspppel

!**********************************************************************
!>    
!!    Function bsppin(ifunc,iparm,ypsi)
!!    
!!    This function returns the matrix element for the
!!    selected basis function.
!!    
!!    @param ifunc : basis function number
!!
!!    @param iparm : basis function parameter number
!!
!!    @param ypsi : independent variable value
!!    
!**********************************************************************
  real*8 function bsppin(ifunc,iparm,ypsi)
      
  include 'eparm.inc'
  include 'modules2.inc'
  include 'modules1.inc'
  implicit none
  integer*4, intent(in) :: ifunc, iparm
  real*8, intent(in) :: ypsi
  integer*4 iorder,nk
  real*8 b1,b2,pptens2,tpsi,tpsi2,w,ypsi1,ypsi2
      
  bsppin = 0.0
  ypsi2 = 1.0
  select case (ifunc)
  case (0)
    bsppin = (ypsi**iparm)/iparm - &
      (ypsi**(kppcur+1))/(kppcur+1)*pcurbd
    bsppin = bsppin - ((ypsi2**iparm)/iparm - &
      (ypsi2**(kppcur+1))/(kppcur+1)*pcurbd)
  case (1)
    tpsi = ypsi - 1.0
    tpsi2 = ypsi2 - 1.0
    bsppin = (tpsi**iparm)/iparm - &
      (tpsi**(kppcur+1))/(kppcur+1)*pcurbd
    bsppin = bsppin - ((tpsi2**iparm)/iparm - &
      (tpsi2**(kppcur+1))/(kppcur+1)*pcurbd)
  case (2)
    bsppin = -((ypsi**iparm)/iparm - &
      (ypsi**(kppcur+1))/(kppcur+1)*pcurbd)
    bsppin = bsppin - (-((ypsi2**iparm)/iparm - &
      (ypsi2**(kppcur+1))/(kppcur+1)*pcurbd))
  case (3)
    nk = (iparm - 1) / 4 + 1
    if(nk .ge. kppknt) nk = kppknt - 1
    if (ypsi .ge. ppknt(nk+1)) then
      bsppin = 0
      return
    endif
    if (1.0 .le. ppknt(nk)) then
      bsppin = 0
      return
    endif
    if (ypsi .ge. ppknt(nk)) then
      ypsi1 = ypsi
    else
      ypsi1 = ppknt(nk)
    endif
    w = ppknt(nk+1) - ppknt(nk)
    if (1.0 .ge. ppknt(nk+1)) then
      ypsi2 = ppknt(nk+1)
    else
      ypsi2 = 1.0
    endif
    if(mod(iparm,4) .eq. 1) b1 = ypsi1
    if(mod(iparm,4) .eq. 2) b1 = (ypsi1**2) / 2.0
    if(mod(iparm,4) .eq. 3) &
      b1 = sin(w*pptens*ypsi1)/w*pptens
    if(mod(iparm,4) .eq. 0) &
      b1 = -cos(w*pptens*ypsi1)/w*pptens
    if(mod(iparm,4) .eq. 1) b2 = ypsi2
    if(mod(iparm,4) .eq. 2) b2 = (ypsi2**2) / 2.0
    if(mod(iparm,4) .eq. 3) &
      b2 = sin(w*pptens*ypsi2)/w*pptens
    if(mod(iparm,4) .eq. 0) &
      b2 = -cos(w*pptens*ypsi2)/w*pptens
    bsppin = b1 - b2
    !     write(6,*)'for ypsi=',ypsi,' integrate from ',ypsi1,' to ',
    !     $                       ypsi2,' = ',bsppin
  case (4)
    nk = (iparm - 1) / 4 + 1
    if(nk .ge. kppknt) nk = kppknt - 1
    if (ypsi .ge. ppknt(nk+1)) then
      bsppin = 0
      return
    endif
    if (1.0 .le. ppknt(nk)) then
      bsppin = 0
      return
    endif
    if (ypsi .ge. ppknt(nk)) then
      ypsi1 = ypsi
    else
      ypsi1 = ppknt(nk)
    endif
    if (1.0 .ge. ppknt(nk+1)) then
      ypsi2 = ppknt(nk+1)
    else
      ypsi2 = 1.0
    endif
    if(mod(iparm,4) .eq. 1) b1 = ypsi1
    if(mod(iparm,4) .eq. 2) b1 = (ypsi1**2) / 2.0
    if(mod(iparm,4) .eq. 3) &
      b1 = sin(pptens*ypsi1)/pptens
    if(mod(iparm,4) .eq. 0) &
      b1 = -cos(pptens*ypsi1)/pptens
    if(mod(iparm,4) .eq. 1) b2 = ypsi2
    if(mod(iparm,4) .eq. 2) b2 = (ypsi2**2) / 2.0
    if(mod(iparm,4) .eq. 3) &
      b2 = sin(pptens*ypsi2)/pptens
    if(mod(iparm,4) .eq. 0) &
      b2 = -cos(pptens*ypsi2)/pptens
    bsppin = b1 - b2
    !     write(6,*)'for ypsi=',ypsi,' integrate from ',ypsi1,' to ',
    !     $                       ypsi2,' = ',bsppin
  case (5)
    iorder = kppcur / (kppknt - 1)
    nk = (iparm - 1) / iorder + 1
    if(nk .ge. kppknt) nk = kppknt - 1
    if (ypsi .ge. ppknt(nk+1)) then
      bsppin = 0
      return
    endif
    if (1.0 .le. ppknt(nk)) then
      bsppin = 0
      return
    endif
    if (ypsi .ge. ppknt(nk)) then
      ypsi1 = ypsi
    else
      ypsi1 = ppknt(nk)
    endif
    if (1.0 .ge. ppknt(nk+1)) then
      ypsi2 = ppknt(nk+1)
    else
      ypsi2 = 1.0
    endif
    w = ppknt(nk+1) - ppknt(nk)
         
    tpsi=(ypsi1**2/2.0 - ypsi1*ppknt(nk))/2
    if (mod(iparm,iorder) .eq. 0) then
      b1 = tpsi**iorder / iorder
    else
      b1 = tpsi**mod(iparm,iorder) / mod(iparm,iorder)
    endif
    tpsi=(ypsi2**2/2.0 - ypsi2*ppknt(nk))/2
    if (mod(iparm,iorder) .eq. 0) then
      b2 = tpsi**iorder / iorder
    else
      b2 = tpsi**mod(iparm,iorder) / mod(iparm,iorder)
    endif
    bsppin = b1 - b2     
    !     write(6,*)'for ypsi=',ypsi,' integrate from ',ypsi1,' to ',
    !     $                       ypsi2,' = ',bsppin
  case (6)
    nk = ((iparm - 1) / 2) + 1
    pptens2 = abs(pptens)*(kppknt-1)/(ppknt(kppknt)-ppknt(1))
    bsppin = 0
    if (nk .gt.1) then
      if (ypsi .le. ppknt(nk)) then
        if (ypsi .le. ppknt(nk-1)) then
          ypsi1 = ppknt(nk-1)
        else
          ypsi1 = ypsi
        endif
        if (1.0 .le. ppknt(nk)) then
          ypsi2 = 1.0
        else
          ypsi2 = ppknt(nk)
        endif
        w = ppknt(nk) - ppknt(nk-1)
        if (mod(iparm,2) .eq. 0) then
          b1 = (cosh(pptens2*(ypsi1-ppknt(nk-1)))/ &
            (pptens2*sinh(pptens2*w)) - (ypsi1*ypsi1/2.0- &
            ppknt(nk-1)*ypsi1)/w) &
            / (pptens2*pptens2)
        else
          b1 = (ypsi1*ypsi1/2.0-ppknt(nk-1)*ypsi1)/w
        endif
        if (mod(iparm,2) .eq. 0) then
          b2 = (cosh(pptens2*(ypsi2-ppknt(nk-1)))/ &
            (pptens2*sinh(pptens2*w)) - (ypsi2*ypsi2/2.0- &
            ppknt(nk-1)*ypsi2)/w) &
            / (pptens2*pptens2)
        else
          b2 = (ypsi2*ypsi2/2.0-ppknt(nk-1)*ypsi2)/w
        endif
        bsppin = bsppin + b1 - b2
      endif
    endif
    if (nk .lt. kppknt) then
      if (ypsi .le. ppknt(nk+1)) then
        if (ypsi .le. ppknt(nk)) then
          ypsi1 = ppknt(nk)
        else
          ypsi1 = ypsi
        endif
        if (1.0 .le. ppknt(nk+1)) then
          ypsi2 = 1.0
        else
          ypsi2 = ppknt(nk+1)
        endif
        w = ppknt(nk+1) - ppknt(nk)
        if (mod(iparm,2) .eq. 0) then
          b1 = (-cosh(pptens2*(ppknt(nk+1)-ypsi1))/ &
            (pptens2*sinh(pptens2*w))-(ypsi1*ppknt(nk+1) &
            -ypsi1*ypsi1/2.0)/w) &
            / (pptens2*pptens2)
        else
          b1 = (ppknt(nk+1)*ypsi1-ypsi1*ypsi1/2.0)/w
        endif
        if (mod(iparm,2) .eq. 0) then
          b2 = (-cosh(pptens2*(ppknt(nk+1)-ypsi2))/ &
            (pptens2*sinh(pptens2*w))-(ypsi2*ppknt(nk+1) &
            -ypsi2*ypsi2/2.0)/w) &
            / (pptens2*pptens2)
        else
          b2 = (ppknt(nk+1)*ypsi2-ypsi2*ypsi2/2.0)/w
        endif
        bsppin = bsppin + b1 - b2
      endif
    endif
  case (7)
    if (iparm .eq. kppcur) then
      bsppin = (ypsi**(kpphord+1))/(kpphord+1)
      bsppin = bsppin - ((ypsi2**(kpphord+1))/(kpphord+1))
    else
      bsppin = (ypsi**iparm)/iparm
      bsppin = bsppin - ((ypsi2**iparm)/iparm)
    endif
  end select

  if (ifunc .ne. kppfnc)  &
    write(6,*)'ifunc .ne. kppfnc ',ifunc,kppfnc
  return
  end function bsppin

!**********************************************************************
!>
!!    In addition to the least squares constraints that
!!    efit already uses, some basis functions have exact
!!    constraints. Most notable is the spline function
!!    whose continuity constraints are exact, not LSE.
!!
!!    @param ncrsp : number of constraint equations
!!
!!    @param crsp : constraint matrix
!!
!!    @param z : value vector
!!
!!    @param nffcoi : array index for seeting up crsp
!!
!**********************************************************************
  subroutine ppcnst(ncrsp,crsp,z,nffcoi)
  include 'eparm.inc'
  include 'modules2.inc'
  include 'modules1.inc'
  implicit none
  real*8 bsppel
  integer*4, intent(in) :: nffcoi
  integer*4, intent(inout) :: ncrsp
  real*8, intent(out) :: crsp(4*(npcurn-2)+6+npcurn*npcurn,nrsmat), &
                         z(3*(npcurn-2)+6+npcurn*npcurn)
  integer*4 i,j,iorder
  real*8 h,wl,wr,pptens2

  select case (kppfnc)
  case (3)
    if (kppknt .gt. 2) then
      !
      !     first set of constraints is that splines must be equal at the knots
      !
      do i = 2,kppknt-1
        ncrsp = ncrsp + 1
        crsp(ncrsp,:) = 0.0
        z(ncrsp) = 0.0
        h = ppknt(i) - ppknt(i-1)
        crsp(ncrsp,nffcoi + 4*(i-2) + 1) = 1.0
        crsp(ncrsp,nffcoi + 4*(i-2) + 2) = ppknt(i)
        crsp(ncrsp,nffcoi + 4*(i-2) + 3) = &
          cos(h * pptens * ppknt(i))
        crsp(ncrsp,nffcoi + 4*(i-2) + 4) = &
          sin(h * pptens * ppknt(i))
               
        h = ppknt(i+1) - ppknt(i)
        crsp(ncrsp,nffcoi + 4*(i-2) + 5) = -1.0
        crsp(ncrsp,nffcoi + 4*(i-2) + 6) = -ppknt(i)
        crsp(ncrsp,nffcoi + 4*(i-2) + 7) = &
          -cos(h * pptens * ppknt(i))
        crsp(ncrsp,nffcoi + 4*(i-2) + 8) = &
          -sin(h * pptens * ppknt(i))
      enddo
      !
      !     second set of constraints is that splines have equal first
      !     derivative at the knots
      !
      do i = 2,kppknt-1
        ncrsp = ncrsp + 1
        crsp(ncrsp,:) = 0.0
        z(ncrsp) = 0.0
        h = ppknt(i) - ppknt(i-1)
        crsp(ncrsp,nffcoi + 4*(i-2) + 1) = 0.0
        crsp(ncrsp,nffcoi + 4*(i-2) + 2) = 1.0
        crsp(ncrsp,nffcoi + 4*(i-2) + 3) = &
          -h * pptens * sin(h * pptens * ppknt(i))
        crsp(ncrsp,nffcoi + 4*(i-2) + 4) = &
          h * pptens * cos(h * pptens * ppknt(i))
               
        h = ppknt(i+1) - ppknt(i)
        crsp(ncrsp,nffcoi + 4*(i-2) + 5) = 0.0
        crsp(ncrsp,nffcoi + 4*(i-2) + 6) = -1.0
        crsp(ncrsp,nffcoi + 4*(i-2) + 7) = &
          h * pptens * sin(h * pptens * ppknt(i))
        crsp(ncrsp,nffcoi + 4*(i-2) + 8) = &
          -h * pptens * cos(h * pptens * ppknt(i))
      enddo
      !
      !     second set of constraints is that splines have equal second
      !     derivative at the knots
      !
      do i = 2,kppknt-1
        ncrsp = ncrsp + 1
        crsp(ncrsp,:) = 0.0
        z(ncrsp) = 0.0
        h = ppknt(i) - ppknt(i-1)
        crsp(ncrsp,nffcoi + 4*(i-2) + 1) = 0.0
        crsp(ncrsp,nffcoi + 4*(i-2) + 2) = 0.0
        crsp(ncrsp,nffcoi + 4*(i-2) + 3) = &
          -h*h*pptens*pptens*cos(h * pptens * ppknt(i))
        crsp(ncrsp,nffcoi + 4*(i-2) + 4) = &
          -h*h*pptens*pptens*sin(h * pptens * ppknt(i))
               
        h = ppknt(i+1) - ppknt(i)
        crsp(ncrsp,nffcoi + 4*(i-2) + 5) = 0.0
        crsp(ncrsp,nffcoi + 4*(i-2) + 6) = 0.0
        crsp(ncrsp,nffcoi + 4*(i-2) + 7) = &
          h*h*pptens*pptens*cos(h * pptens * ppknt(i))
        crsp(ncrsp,nffcoi + 4*(i-2) + 8) = &
          h*h*pptens*pptens*sin(h * pptens * ppknt(i))
      enddo
    endif
  case (4)
    if (kppknt .le. 2) then
      !
      !     first set of constraints is that splines must be equal at the knots
      !
      do i = 2,kppknt-1
        ncrsp = ncrsp + 1
        crsp(ncrsp,:) = 0.0
        z(ncrsp) = 0.0
        crsp(ncrsp,nffcoi + 4*(i-2) + 1) = 1.0
        crsp(ncrsp,nffcoi + 4*(i-2) + 2) = ppknt(i)
        crsp(ncrsp,nffcoi + 4*(i-2) + 3) = &
          cos(pptens * ppknt(i))
        crsp(ncrsp,nffcoi + 4*(i-2) + 4) = &
          sin(pptens * ppknt(i))
               
        crsp(ncrsp,nffcoi + 4*(i-2) + 5) = -1.0
        crsp(ncrsp,nffcoi + 4*(i-2) + 6) = -ppknt(i)
        crsp(ncrsp,nffcoi + 4*(i-2) + 7) = &
          -cos(pptens * ppknt(i))
        crsp(ncrsp,nffcoi + 4*(i-2) + 8) = &
          -sin(pptens * ppknt(i))
      enddo
      !
      !     second set of constraints is that splines have equal first
      !     derivative at the knots
      !
      do i = 2,kppknt-1
        ncrsp = ncrsp + 1
        crsp(ncrsp,:) = 0.0
        z(ncrsp) = 0.0
        crsp(ncrsp,nffcoi + 4*(i-2) + 1) = 0.0
        crsp(ncrsp,nffcoi + 4*(i-2) + 2) = 1.0
        crsp(ncrsp,nffcoi + 4*(i-2) + 3) = &
          -pptens * sin(pptens * ppknt(i))
        crsp(ncrsp,nffcoi + 4*(i-2) + 4) = &
          pptens * cos(pptens * ppknt(i))
               
        crsp(ncrsp,nffcoi + 4*(i-2) + 5) = 0.0
        crsp(ncrsp,nffcoi + 4*(i-2) + 6) = -1.0
        crsp(ncrsp,nffcoi + 4*(i-2) + 7) = &
          pptens * sin(pptens * ppknt(i))
        crsp(ncrsp,nffcoi + 4*(i-2) + 8) = &
          -pptens * cos(pptens * ppknt(i))
      enddo
      !
      !     second set of constraints is that splines have equal second
      !     derivative at the knots
      !
      do i = 2,kppknt-1
        ncrsp = ncrsp + 1
        crsp(ncrsp,:) = 0.0
        z(ncrsp) = 0.0
        crsp(ncrsp,nffcoi + 4*(i-2) + 1) = 0.0
        crsp(ncrsp,nffcoi + 4*(i-2) + 2) = 0.0
        crsp(ncrsp,nffcoi + 4*(i-2) + 3) = &
          -pptens*pptens*cos(pptens * ppknt(i))
        crsp(ncrsp,nffcoi + 4*(i-2) + 4) = &
          -pptens*pptens*sin(pptens * ppknt(i))
               
        crsp(ncrsp,nffcoi + 4*(i-2) + 5) = 0.0
        crsp(ncrsp,nffcoi + 4*(i-2) + 6) = 0.0
        crsp(ncrsp,nffcoi + 4*(i-2) + 7) = &
          pptens*pptens*cos(pptens * ppknt(i))
        crsp(ncrsp,nffcoi + 4*(i-2) + 8) = &
          pptens*pptens*sin(pptens * ppknt(i))
      enddo
    endif
  case (5)
    iorder = kppcur / (kppknt - 1)
    if (kppknt .le. 2) then
      !
      !     first set of constraints is that splines must be equal at the knots
      !
      do i = 2,kppknt-1
        ncrsp = ncrsp + 1
        crsp(ncrsp,:) = 0.0
        z(ncrsp) = 0.0
        do j= 1,iorder
          crsp(ncrsp,nffcoi + iorder*(i-2) + j) = 1.0
        enddo
        crsp(ncrsp,nffcoi + iorder*(i-1) + 1) = -1.0
      enddo
      !
      !     second set of constraints is that splines have equal first
      !     derivative at the knots
      !
      do i = 2,kppknt-1
        ncrsp = ncrsp + 1
        crsp(ncrsp,:) = 0.0
        z(ncrsp) = 0.0
        do j= 2,iorder
          crsp(ncrsp,nffcoi + iorder*(i-2) + j) = (j-1)
        enddo
        crsp(ncrsp,nffcoi + iorder*(i-1) + 2) = -1.0
      enddo
      !
      !     second set of constraints is that splines have equal second
      !     derivative at the knots
      !
      do i = 2,kppknt-1
        ncrsp = ncrsp + 1
        crsp(ncrsp,:) = 0.0
        z(ncrsp) = 0.0
        do j= 3,iorder
          crsp(ncrsp,nffcoi + iorder*(i-2) + j) = (j-1)*(j-2)
        enddo
        crsp(ncrsp,nffcoi + iorder*(i-1) + 3) = -2.0
      enddo
    endif
  case (6)
    !
    !     first set of constraints is that splines have equal first
    !     derivative at the knots
    !
    if (kppknt .gt. 2) then
      pptens2 = abs(pptens)*(kppknt-1)/(ppknt(kppknt)-ppknt(1))
      do i = 2,kppknt-1
        ncrsp = ncrsp + 1
        crsp(ncrsp,:) = 0.0
        z(ncrsp) = 0.0
        wr = ppknt(i+1) - ppknt(i)
        wl = ppknt(i) - ppknt(i-1)
        crsp(ncrsp,nffcoi + 2*(i-1) - 1) = 1.0/wl
        crsp(ncrsp,nffcoi + 2*(i-1) + 0) = -(-pptens2/ &
          sinh(pptens2*wl) + 1.0/wl)/(pptens2*pptens2)
        crsp(ncrsp,nffcoi + 2*(i-1) + 1) = -1.0/wr - 1.0/wl
        crsp(ncrsp,nffcoi + 2*(i-1) + 2) = (-pptens2* &
          cosh(pptens2*wr)/sinh(pptens2*wr) &
          + 1.0/wr)/(pptens2*pptens2) - (pptens2* &
          cosh(pptens2*wl)/sinh(pptens2*wl) &
          - 1.0/wl)/(pptens2*pptens2)
        crsp(ncrsp,nffcoi + 2*(i-1) + 3) = 1.0/wr
        crsp(ncrsp,nffcoi + 2*(i-1) + 4) = (pptens2/ &
          sinh(pptens2*wr) - 1.0/wr)/(pptens2*pptens2)
      enddo
    endif
    do i = 1,kppknt
      if (kppbdry(i) .eq. 1) then
        ncrsp = ncrsp + 1
        crsp(ncrsp,:) = 0.0
        z(ncrsp) = ppbdry(i)*darea
        crsp(ncrsp,nffcoi+2*i - 1) = 1.0
      endif
      if (kpp2bdry(i) .eq. 1) then
        ncrsp = ncrsp + 1
        crsp(ncrsp,:) = 0.0
        z(ncrsp) = pp2bdry(i)*darea
        crsp(ncrsp,nffcoi+2*i) = 1.0
      endif
    enddo
  end select
  select case (kppfnc)
  case (3,4,5,6,7)
    if (pcurbd .eq. 1.0) then
      ncrsp = ncrsp + 1
      crsp(ncrsp,:) = 0.0
      z(ncrsp) = 0.0
      do j = 1,kppcur
        crsp(ncrsp,nffcoi+j) = bsppel(kppfnc,j,1.0)
      enddo
    endif
  end select

  return
  end subroutine ppcnst

!**********************************************************************
!>
!!     Store the solution coefs into ppbdry and pp2bdry
!!     
!**********************************************************************
  subroutine ppstore()
  include 'eparm.inc'
  include 'modules2.inc'
  include 'modules1.inc'
  implicit none
  integer*4 i

  if (kppfnc .ge. 0 .and. kppfnc .le. 2) then
    do i = 1,kppcur
      ppbdry(i) = brsp(nfsum+i)/darea
      pp2bdry(i) = 0.0
    enddo
  elseif (kppfnc .eq. 6) then
    do i = 1,kppknt
      if(kppbdry(i) .ne. 1)   ppbdry(i) = brsp(nfsum+2*i-1)/darea
      if(kpp2bdry(i) .ne. 1) pp2bdry(i) = brsp(nfsum+2*i)/darea
    enddo
  endif
  return
  end subroutine ppstore
