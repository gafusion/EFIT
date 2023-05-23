!**********************************************************************
!>
!!    Function bswwel(ifunc,iparm,ypsi)
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
      real*8 function bswwel(ifunc,iparm,ypsi)
      
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      integer*4, intent(in) :: ifunc, iparm
      real*8, intent(in) :: ypsi
      integer*4 iorder,nk
      real*8 tpsi,w,wwtens2
      
      bswwel = 0.0
      select case (ifunc)
      case (0)
        if (iparm.eq.1) then
          bswwel=1.0 - ypsi**kwwcur*wcurbd
        else
          bswwel = ypsi**(iparm - 1) - ypsi**kwwcur*wcurbd
        endif
      case (1)
        tpsi = ypsi - 1.0
        if (iparm.eq.1) then
          bswwel=1.0 - tpsi**kwwcur*wcurbd
        else
          bswwel = tpsi**(iparm - 1) - tpsi**kwwcur*wcurbd
        endif
      case (2)
        if (iparm.eq.1) then
          bswwel=-(1.0 - ypsi**kwwcur*wcurbd)
        else
          bswwel = -(ypsi**(iparm - 1) - ypsi**kwwcur*wcurbd)
        endif
      case (3)
        nk = (iparm - 1) / 4 + 1
        if(nk .ge. kwwknt) nk = kwwknt - 1
        if ((nk .lt. (kwwknt - 1) .and. ypsi .lt. wwknt(nk+1) &
          .and.  ypsi .ge. wwknt(nk)) &
          .or. (nk .eq. (kwwknt - 1) .and. ypsi .le. wwknt(nk+1) &
          .and.  ypsi .ge. wwknt(nk))) then
          w = wwknt(nk+1) - wwknt(nk)
          if(mod(iparm,4) .eq. 1) bswwel = 1.0
          if(mod(iparm,4) .eq. 2) bswwel = ypsi
          if(mod(iparm,4) .eq. 3) bswwel = cos(w*wwtens*ypsi)
          if(mod(iparm,4) .eq. 0) bswwel = sin(w*wwtens*ypsi)
        endif
      case (4)
        nk = (iparm - 1) / 4 + 1
        if(nk .ge. kwwknt) nk = kwwknt - 1
        if ((nk .lt. (kwwknt - 1) .and. ypsi .lt. wwknt(nk+1) &
          .and.  ypsi .ge. wwknt(nk)) &
          .or. (nk .eq. (kwwknt - 1) .and. ypsi .le. wwknt(nk+1) &
          .and.  ypsi .ge. wwknt(nk))) then
          if(mod(iparm,4) .eq. 1) bswwel = 1.0
          if(mod(iparm,4) .eq. 2) bswwel = ypsi
          if(mod(iparm,4) .eq. 3) bswwel = cos(wwtens*ypsi)
          if(mod(iparm,4) .eq. 0) bswwel = sin(wwtens*ypsi)
        endif
      case (5)
        iorder = kwwcur / (kwwknt - 1)
        nk = (iparm - 1) / iorder + 1
        if(nk .ge. kwwknt) nk = kwwknt - 1
        if ((nk .lt. (kwwknt - 1) .and. ypsi .lt. wwknt(nk+1) &
          .and.  ypsi .ge. wwknt(nk)) &
          .or. (nk .eq. (kwwknt - 1) .and. ypsi .le. wwknt(nk+1) &
          .and.  ypsi .ge. wwknt(nk))) then
          w = wwknt(nk+1) - wwknt(nk)
          tpsi = (ypsi - wwknt(nk)) / w
          if (mod(iparm,iorder) .eq. 0) then
            if (iorder.eq.1) then
              bswwel = 1.0
            else
              bswwel = tpsi**(iorder-1)
            endif
          else
            bswwel = tpsi**(mod(iparm,iorder)-1)
          endif
        endif
      case (6)
        nk = ((iparm - 1) / 2) + 1
        wwtens2 = abs(wwtens)*(kwwknt-1)/(wwknt(kwwknt)-wwknt(1))
        if (nk .gt. 1 ) then
          if (ypsi .le. wwknt(nk) .and.  &
            ypsi .ge. wwknt(nk-1)) then
            w = wwknt(nk) - wwknt(nk-1)
            if (mod(iparm,2) .eq. 0) then
              bswwel = (sinh(wwtens2*(ypsi-wwknt(nk-1)))/ &
                sinh(wwtens2*w) - (ypsi-wwknt(nk-1))/w) &
                / (wwtens2*wwtens2)
            else
              bswwel = (ypsi-wwknt(nk-1))/w
            endif
          endif
        endif
        if (nk .lt. kwwknt) then
          if (ypsi .ge. wwknt(nk) .and.  &
            ypsi .le. wwknt(nk+1)) then
            w = wwknt(nk+1) - wwknt(nk)
            if(mod(iparm,2) .eq. 0) then
              bswwel = (sinh(wwtens2*(wwknt(nk+1)-ypsi))/ &
                sinh(wwtens2*w) - (wwknt(nk+1)-ypsi)/w) &
                / (wwtens2*wwtens2)
            else
              bswwel = (wwknt(nk+1) - ypsi)/w
            endif
          endif
        endif
      case (7)
        if (iparm.eq.kwwcur) then
          bswwel = ypsi**(kwwhord)
        elseif (iparm .eq. 1) then
          bswwel = 1.0
        else
          bswwel = ypsi**(iparm - 1)
        endif
      end select

      if(ifunc .ne. kwwfnc)  &
        write(6,*)'ifunc .ne. kwwfnc ',ifunc,kwwfnc
      return
      end function bswwel

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
      real*8 function bswwpel(ifunc,iparm,ypsi)

      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      integer*4, intent(in) :: ifunc, iparm
      real*8, intent(in) :: ypsi
      integer*4 iorder,nk
      real*8 jparm,tpsi,w,wwtens2

      bswwpel = 0.0
      select case (ifunc)
      case (0)
        if (iparm.eq.1) then
          if (kwwcur.eq.1) then
            bswwpel = -wcurbd
          else
            bswwpel = -kwwcur*ypsi**(kwwcur-1)*wcurbd
          endif
        elseif (iparm.eq.2) then
          bswwpel = 1. - kwwcur*ypsi**(kwwcur-1)*wcurbd
        elseif (iparm.gt.2) then
          bswwpel = (iparm - 1)*ypsi**(iparm - 2) - &
            kwwcur*ypsi**(kwwcur-1)*wcurbd
        endif
      case (1)
        tpsi = ypsi - 1.0
        if (iparm.eq.1) then
          if (kwwcur.eq.1) then
            bswwpel = -wcurbd
          else
            bswwpel = -kwwcur*tpsi**(kwwcur-1)*wcurbd
          endif
        elseif (iparm.eq.2) then
          bswwpel = 1. - kwwcur*tpsi**(kwwcur-1)*wcurbd
        elseif (iparm.gt.2) then
          bswwpel = (iparm - 1)*tpsi**(iparm - 2) - &
            kwwcur*tpsi**(kwwcur-1)*wcurbd
        endif
      case (2)
        if (iparm.eq.1) then
          if (kwwcur.eq.1) then
            bswwpel = -wcurbd
          else
            bswwpel = -kwwcur*ypsi**(kwwcur-1)*wcurbd
          endif
        elseif (iparm.eq.2) then
          bswwpel = 1. - kwwcur*ypsi**(kwwcur-1)*wcurbd
        elseif (iparm.gt.2) then
          bswwpel = (iparm - 1)*ypsi**(iparm - 2) - &
            kwwcur*ypsi**(kwwcur-1)*wcurbd
        endif
        bswwpel = - bswwpel
      case (3)
        nk = (iparm - 1) / 4 + 1
        if(nk .ge. kwwknt) nk = kwwknt - 1
        if ((nk .lt. (kwwknt - 1) .and. ypsi .lt. wwknt(nk+1) &
          .and.  ypsi .ge. wwknt(nk)) &
          .or. (nk .eq. (kwwknt - 1) .and. ypsi .le. wwknt(nk+1) &
          .and.  ypsi .ge. wwknt(nk))) then
          w = wwknt(nk+1) - wwknt(nk)
          if(mod(iparm,4) .eq. 1) bswwpel = 0.0
          if(mod(iparm,4) .eq. 2) bswwpel = 1.0
          if(mod(iparm,4) .eq. 3) bswwpel =  &
            -w*wwtens*sin(w*wwtens*ypsi)
          if(mod(iparm,4) .eq. 0) bswwpel =  &
            w*wwtens*cos(w*wwtens*ypsi)
        endif
      case (4)
        nk = (iparm - 1) / 4 + 1
        if(nk .ge. kwwknt) nk = kwwknt - 1
        if ((nk .lt. (kwwknt - 1) .and. ypsi .lt. wwknt(nk+1) &
          .and.  ypsi .ge. wwknt(nk)) &
          .or. (nk .eq. (kwwknt - 1) .and. ypsi .le. wwknt(nk+1) &
          .and.  ypsi .ge. wwknt(nk))) then
          if(mod(iparm,4) .eq. 1) bswwpel = 0.0
          if(mod(iparm,4) .eq. 2) bswwpel = 1.0
          if(mod(iparm,4) .eq. 3) bswwpel = -wwtens*sin(wwtens*ypsi)
          if(mod(iparm,4) .eq. 0) bswwpel = wwtens*cos(wwtens*ypsi)
        endif
      case (5)
        iorder = kwwcur / (kwwknt - 1)
        nk = (iparm - 1) / iorder + 1
        if(nk .ge. kwwknt) nk = kwwknt - 1
        if ((nk .lt. (kwwknt - 1) .and. ypsi .lt. wwknt(nk+1) &
          .and.  ypsi .ge. wwknt(nk)) &
          .or. (nk .eq. (kwwknt - 1) .and. ypsi .le. wwknt(nk+1) &
          .and.  ypsi .ge. wwknt(nk))) then
          w = wwknt(nk+1) - wwknt(nk)
          tpsi = (ypsi - wwknt(nk)) / w
          jparm = mod(iparm,iorder)
          if (jparm.eq.1) then
            bswwpel = 0.0
          elseif (jparm.eq.2) then
            bswwpel = 1./w
          elseif (jparm.gt.2) then
            bswwpel = (jparm - 1)/w*tpsi**(jparm - 2)
          endif
        endif
      case (6)
        nk = ((iparm - 1) / 2) + 1
        wwtens2 = abs(wwtens)*(kwwknt-1)/(wwknt(kwwknt)-wwknt(1))
        if (nk .gt. 1) then
          if (ypsi .le. wwknt(nk) .and.  &
            ypsi .ge. wwknt(nk-1)) then
            w = wwknt(nk) - wwknt(nk-1)
            if (mod(iparm,2) .eq. 0) then
              bswwpel = (wwtens2*cosh(wwtens2* &
                (ypsi-wwknt(nk-1)))/sinh(wwtens2*w) - (1.0/w)) &
                / (wwtens2*wwtens2)
            else
              bswwpel = 1.0/w
            endif
          endif
        endif
        if (nk .lt. kwwknt) then
          if (ypsi .ge. wwknt(nk) .and.  &
            ypsi .le. wwknt(nk+1)) then
            w = wwknt(nk+1) - wwknt(nk)
            if (mod(iparm,2) .eq. 0) then
              bswwpel = (-wwtens2*cosh(wwtens2* &
                (wwknt(nk+1)-ypsi))/sinh(wwtens2*w)+(1.0/w)) &
                / (wwtens2*wwtens2)
            else
              bswwpel = -1.0/w
            endif
          endif
        endif
      case (7)
        if (iparm.eq.kwwcur) then
          bswwpel = kwwhord*ypsi**(kwwhord-1)
        elseif (iparm.eq.1 ) then
          bswwpel = 0.
        elseif (iparm.eq.2 ) then
          bswwpel = 1.
        elseif (iparm.gt.2) then
          bswwpel = (iparm - 1)*ypsi**(iparm - 2)
        endif
      end select
      return
      end function bswwpel
 
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
      real*8 function bswwin(ifunc,iparm,ypsi)
      
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      integer*4, intent(in) :: ifunc, iparm
      real*8, intent(in) :: ypsi
      integer*4 iorder,nk
      real*8 b1,b2,wwtens2,tpsi,tpsi2,w,ypsi1,ypsi2

      bswwin = 0.0
      ypsi2 = 1.0
      select case (ifunc)
      case (0)
         bswwin = (ypsi**iparm)/iparm - &
              (ypsi**(kwwcur+1))/(kwwcur+1)*wcurbd                
         bswwin = bswwin - ((ypsi2**iparm)/iparm - &
              (ypsi2**(kwwcur+1))/(kwwcur+1)*wcurbd)
      case (1)
         tpsi = ypsi - 1.0
         tpsi2 = ypsi2 - 1.0
         bswwin = (tpsi**iparm)/iparm - &
              (tpsi**(kwwcur+1))/(kwwcur+1)*wcurbd                
         bswwin = bswwin - ((tpsi2**iparm)/iparm - &
              (tpsi2**(kwwcur+1))/(kwwcur+1)*wcurbd)         
      case (2)
         bswwin = -((ypsi**iparm)/iparm - &
              (ypsi**(kwwcur+1))/(kwwcur+1)*wcurbd)
         bswwin = bswwin - (-((ypsi2**iparm)/iparm - &
              (ypsi2**(kwwcur+1))/(kwwcur+1)*wcurbd))
      case (3)
         nk = (iparm - 1) / 4 + 1
         if(nk .ge. kwwknt) nk = kwwknt - 1
         if (ypsi .ge. wwknt(nk+1)) then
            bswwin = 0
            return
         endif
         if (1.0 .le. wwknt(nk)) then
            bswwin = 0
            return
         endif
         if (ypsi .ge. wwknt(nk)) then
            ypsi1 = ypsi
         else
            ypsi1 = wwknt(nk)
         endif
         w = wwknt(nk+1) - wwknt(nk)
         if (1.0 .ge. wwknt(nk+1)) then
            ypsi2 = wwknt(nk+1)
         else
            ypsi2 = 1.0
         endif
         if(mod(iparm,4) .eq. 1) b1 = ypsi1
         if(mod(iparm,4) .eq. 2) b1 = (ypsi1**2) / 2.0
         if(mod(iparm,4) .eq. 3) &
              b1 = sin(w*wwtens*ypsi1)/w*wwtens
         if(mod(iparm,4) .eq. 0) &
              b1 = -cos(w*wwtens*ypsi1)/w*wwtens
         if(mod(iparm,4) .eq. 1) b2 = ypsi2
         if(mod(iparm,4) .eq. 2) b2 = (ypsi2**2) / 2.0
         if(mod(iparm,4) .eq. 3) &
              b2 = sin(w*wwtens*ypsi2)/w*wwtens
         if(mod(iparm,4) .eq. 0) &
              b2 = -cos(w*wwtens*ypsi2)/w*wwtens
         bswwin = b1 - b2
         !     write(6,*)'for ypsi=',ypsi,' integrate from ',ypsi1,' to ',
         !     $                       ypsi2,' = ',bswwin
      case (4)
         nk = (iparm - 1) / 4 + 1
         if(nk .ge. kwwknt) nk = kwwknt - 1
         if (ypsi .ge. wwknt(nk+1)) then
            bswwin = 0
            return
         endif
         if (1.0 .le. wwknt(nk)) then
            bswwin = 0
            return
         endif
         if (ypsi .ge. wwknt(nk)) then
            ypsi1 = ypsi
         else
            ypsi1 = wwknt(nk)
         endif
         if (1.0 .ge. wwknt(nk+1)) then
            ypsi2 = wwknt(nk+1)
         else
            ypsi2 = 1.0
         endif
         if(mod(iparm,4) .eq. 1) b1 = ypsi1
         if(mod(iparm,4) .eq. 2) b1 = (ypsi1**2) / 2.0
         if(mod(iparm,4) .eq. 3) &
              b1 = sin(wwtens*ypsi1)/wwtens
         if(mod(iparm,4) .eq. 0) &
              b1 = -cos(wwtens*ypsi1)/wwtens
         if(mod(iparm,4) .eq. 1) b2 = ypsi2
         if(mod(iparm,4) .eq. 2) b2 = (ypsi2**2) / 2.0
         if(mod(iparm,4) .eq. 3) &
              b2 = sin(wwtens*ypsi2)/wwtens
         if(mod(iparm,4) .eq. 0) &
              b2 = -cos(wwtens*ypsi2)/wwtens
         bswwin = b1 - b2
         !     write(6,*)'for ypsi=',ypsi,' integrate from ',ypsi1,' to ',
         !     $                       ypsi2,' = ',bswwin
      case (5)
         iorder = kwwcur / (kwwknt - 1)
         nk = (iparm - 1) / iorder + 1
         if(nk .ge. kwwknt) nk = kwwknt - 1
         if (ypsi .ge. wwknt(nk+1)) then
            bswwin = 0
            return
         endif
         if (1.0 .le. wwknt(nk)) then
            bswwin = 0
            return
         endif
         if (ypsi .ge. wwknt(nk)) then
            ypsi1 = ypsi
         else
            ypsi1 = wwknt(nk)
         endif
         if (1.0 .ge. wwknt(nk+1)) then
            ypsi2 = wwknt(nk+1)
         else
            ypsi2 = 1.0
         endif
         w = wwknt(nk+1) - wwknt(nk)
         
         tpsi=(ypsi1**2/2.0 - ypsi1*wwknt(nk))/2
         if (mod(iparm,iorder) .eq. 0) then
            b1 = tpsi**iorder / iorder
         else
            b1 = tpsi**mod(iparm,iorder) / mod(iparm,iorder)
         endif
         tpsi=(ypsi2**2/2.0 - ypsi2*wwknt(nk))/2
         if (mod(iparm,iorder) .eq. 0) then
            b2 = tpsi**iorder / iorder
         else
            b2 = tpsi**mod(iparm,iorder) / mod(iparm,iorder)
         endif
         bswwin = b1 - b2
         !     write(6,*)'for ypsi=',ypsi,' integrate from ',ypsi1,' to ',
         !     $                       ypsi2,' = ',bswwin
      case (6)
         nk = ((iparm - 1) / 2) + 1
         wwtens2 = abs(wwtens)*(kwwknt-1)/(wwknt(kwwknt)-wwknt(1))
         bswwin = 0
         if (nk .gt.1) then
            if (ypsi .le. wwknt(nk)) then
               if (ypsi .le. wwknt(nk-1)) then
                  ypsi1 = wwknt(nk-1)
               else
                  ypsi1 = ypsi
               endif
               if (1.0 .le. wwknt(nk)) then
                  ypsi2 = 1.0
               else
                  ypsi2 = wwknt(nk)
               endif
               w = wwknt(nk) - wwknt(nk-1)
               if (mod(iparm,2) .eq. 0) then
                  b1 = (cosh(wwtens2*(ypsi1-wwknt(nk-1)))/ &
                       (wwtens2*sinh(wwtens2*w)) - (ypsi1*ypsi1/2.0- &
                       wwknt(nk-1)*ypsi1)/w) &
                       / (wwtens2*wwtens2)
               else
                  b1 = (ypsi1*ypsi1/2.0-wwknt(nk-1)*ypsi1)/w
               endif
               if (mod(iparm,2) .eq. 0) then
                  b2 = (cosh(wwtens2*(ypsi2-wwknt(nk-1)))/ &
                       (wwtens2*sinh(wwtens2*w)) - (ypsi2*ypsi2/2.0- &
                       wwknt(nk-1)*ypsi2)/w) &
                       / (wwtens2*wwtens2)
               else
                  b2 = (ypsi2*ypsi2/2.0-wwknt(nk-1)*ypsi2)/w
               endif
               bswwin = bswwin + b1 - b2
            endif
         endif
         if (nk .lt. kwwknt) then
            if (ypsi .le. wwknt(nk+1)) then
               if (ypsi .le. wwknt(nk)) then
                  ypsi1 = wwknt(nk)
               else
                  ypsi1 = ypsi
               endif
               if (1.0 .le. wwknt(nk+1)) then
                  ypsi2 = 1.0
               else
                  ypsi2 = wwknt(nk+1)
               endif
               w = wwknt(nk+1) - wwknt(nk)
               if (mod(iparm,2) .eq. 0) then
                  b1 = (-cosh(wwtens2*(wwknt(nk+1)-ypsi1))/ &
                       (wwtens2*sinh(wwtens2*w))-(ypsi1*wwknt(nk+1) &
                       -ypsi1*ypsi1/2.0)/w) &
                       / (wwtens2*wwtens2)
               else
                  b1 = (wwknt(nk+1)*ypsi1-ypsi1*ypsi1/2.0)/w
               endif
               if (mod(iparm,2) .eq. 0) then
                  b2 = (-cosh(wwtens2*(wwknt(nk+1)-ypsi2))/ &
                       (wwtens2*sinh(wwtens2*w))-(ypsi2*wwknt(nk+1) &
                       -ypsi2*ypsi2/2.0)/w) &
                       / (wwtens2*wwtens2)
               else
                  b2 = (wwknt(nk+1)*ypsi2-ypsi2*ypsi2/2.0)/w
               endif
               bswwin = bswwin + b1 - b2
            endif
         endif
      case (7)
         if (iparm .eq. kwwcur) then
            bswwin = (ypsi**(kwwhord+1))/(kwwhord+1)
            bswwin = bswwin - ((ypsi2**(kwwhord+1))/(kwwhord+1))
         else
            bswwin = (ypsi**iparm)/iparm
            bswwin = bswwin - ((ypsi2**iparm)/iparm)
         endif
      end select

      if(ifunc .ne. kwwfnc)  &
           write(6,*)'ifunc .ne. kwwfnc ',ifunc,kwwfnc
      return
      end function bswwin 
          
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
      subroutine wwcnst(ncrsp,crsp,z,nffcoi)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      real*8 bswwel
      integer*4, intent(in) :: nffcoi
      integer*4, intent(inout) :: ncrsp
      real*8, intent(out) :: crsp(4*(npcurn-2)+6+npcurn*npcurn,nrsmat), &
                             z(4*(npcurn-2)+6+npcurn*npcurn)
      integer*4 i,j,iorder
      real*8 h,w,wwtens2

      select case (kwwfnc)
      case (3)
        if (kwwknt .gt. 2) then
          !
          !     first set of constraints is that splines must be equal at the knots
          !
          do i = 2,kwwknt-1
            ncrsp = ncrsp + 1
            crsp(ncrsp,:) = 0.0
            z(ncrsp) = 0.0
            h = wwknt(i) - wwknt(i-1)
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 1) = 1.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 2) = wwknt(i)
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 3) = &
              cos(h * wwtens * wwknt(i))
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 4) =  &
              sin(h * wwtens * wwknt(i))
               
            h = wwknt(i+1) - wwknt(i)
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 5) = -1.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 6) = -wwknt(i)
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 7) = &
              -cos(h * wwtens * wwknt(i))
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 8) =  &
              -sin(h * wwtens * wwknt(i))
          enddo
          !
          !     second set of constraints is that splines have equal first
          !     derivative at the knots
          !
          do i = 2,kwwknt-1
            ncrsp = ncrsp + 1
            crsp(ncrsp,:) = 0.0
            z(ncrsp) = 0.0
            h = wwknt(i) - wwknt(i-1)
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 1) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 2) = 1.0
            crsp(ncrsp,nffcoi + kppcur + kffcur + &
              4*(i-2) + 3) = &
              -h * wwtens * sin(h * wwtens * wwknt(i))
            crsp(ncrsp,nffcoi + kppcur + kffcur + &
              4*(i-2) + 4) =  &
              h * wwtens * cos(h * wwtens * wwknt(i))
               
            h = wwknt(i+1) - wwknt(i)
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 5) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 6) = -1.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 7) = &
              h * wwtens * sin(h * wwtens * wwknt(i))
            crsp(ncrsp,nffcoi + kppcur + kffcur + 4*(i-2) + 8) = &
              -h * wwtens * cos(h * wwtens * wwknt(i))
          enddo
          !
          !     second set of constraints is that splines have equal second
          !     derivative at the knots
          !
          do i = 2,kwwknt-1
            ncrsp = ncrsp + 1
            crsp(ncrsp,:) = 0.0
            z(ncrsp) = 0.0
            h = wwknt(i) - wwknt(i-1)
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 1) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 2) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 3) = &
              -h*h*wwtens*wwtens*cos(h * wwtens * wwknt(i))
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 4) =  &
              -h*h*wwtens*wwtens*sin(h * wwtens * wwknt(i))
               
            h = wwknt(i+1) - wwknt(i)
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 5) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 6) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 7) = &
              h*h*wwtens*wwtens*cos(h * wwtens * wwknt(i))
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 8) =  &
              h*h*wwtens*wwtens*sin(h * wwtens * wwknt(i))
          enddo
        endif
      case (4)
        if (kwwknt .le. 2) then
          !
          !     first set of constraints is that splines must be equal at the knots
          !
          do i = 2,kwwknt-1
            ncrsp = ncrsp + 1
            crsp(ncrsp,:) = 0.0
            z(ncrsp) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 1) = 1.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 2) = wwknt(i)
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 3) = &
              cos(wwtens * wwknt(i))
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 4) =  &
              sin(wwtens * wwknt(i))
               
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 5) = -1.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 6) = -wwknt(i)
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 7) = &
              -cos(wwtens * wwknt(i))
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 8) =  &
              -sin(wwtens * wwknt(i))
          enddo
          !
          !     second set of constraints is that splines have equal first
          !     derivative at the knots
          !
          do i = 2,kwwknt-1
            ncrsp = ncrsp + 1
            crsp(ncrsp,:) = 0.0
            z(ncrsp) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 1) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 2) = 1.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 3) = &
              -wwtens * sin(wwtens * wwknt(i))
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 4) =  &
              wwtens * cos(wwtens * wwknt(i))
               
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 5) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 6) = -1.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 7) = &
              wwtens * sin(wwtens * wwknt(i))
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 8) =  &
              -wwtens * cos(wwtens * wwknt(i))
          enddo
          !
          !     second set of constraints is that splines have equal second
          !     derivative at the knots
          !
          do i = 2,kwwknt-1
            ncrsp = ncrsp + 1
            crsp(ncrsp,:) = 0.0
            z(ncrsp) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 1) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 2) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 3) = &
              -wwtens*wwtens*cos(wwtens * wwknt(i))
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 4) =  &
              -wwtens*wwtens*sin(wwtens * wwknt(i))
               
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 5) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 6) = 0.0
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 7) = &
              wwtens*wwtens*cos(wwtens * wwknt(i))
            crsp(ncrsp,nffcoi + kppcur + kffcur &
              + 4*(i-2) + 8) =  &
              wwtens*wwtens*sin(wwtens * wwknt(i))
          enddo
        endif
      case (5)
        iorder = kwwcur / (kwwknt - 1)
        if (kwwknt .le. 2) then
          !
          !     first set of constraints is that splines must be equal at the knots
          !
          do i = 2,kwwknt-1
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
          do i = 2,kwwknt-1
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
          do i = 2,kwwknt-1
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
         if (kwwknt .gt. 2) then
            wwtens2 = abs(wwtens)*(kwwknt-1)/(wwknt(kwwknt)-wwknt(1))
            do i = 2,kwwknt-1
               ncrsp = ncrsp + 1
               crsp(ncrsp,:) = 0.0
               z(ncrsp) = 0.0
               w = wwknt(i+1) - wwknt(i)
               crsp(ncrsp,nffcoi + kppcur + kffcur &
                    + 2*(i-1) + 1) = -1.0/w
               crsp(ncrsp,nffcoi + kppcur + kffcur &
                    + 2*(i-1) + 2) = (-wwtens2* &
                    cosh(wwtens2*w)/sinh(wwtens2*w) &
                    + 1.0/w)/(wwtens2*wwtens2)
               crsp(ncrsp,nffcoi + kppcur + kffcur &
                    + 2*(i-1) + 3) = 1.0/w
               crsp(ncrsp,nffcoi + kppcur + kffcur &
                    + 2*(i-1) + 4) = (wwtens2/ &
                    sinh(wwtens2*w) - 1.0/w)/(wwtens2*wwtens2)
               
               w = wwknt(i) - wwknt(i-1)
               crsp(ncrsp,nffcoi + kppcur + kffcur &
                    + 2*(i-1) - 1) = 1.0/w
               crsp(ncrsp,nffcoi + kppcur + kffcur &
                    + 2*(i-1) + 0) = -(-wwtens2/ &
                    sinh(wwtens2*w) + 1.0/w)/(wwtens2*wwtens2)
               crsp(ncrsp,nffcoi + kppcur + kffcur &
                    + 2*(i-1) + 1) =  &
                    crsp(ncrsp,nffcoi + kppcur + kffcur &
                    + 2*(i-1) + 1) - 1.0/w
               crsp(ncrsp,nffcoi + kppcur + kffcur &
                    + 2*(i-1) + 2) =  &
                    crsp(ncrsp,nffcoi + kppcur + kffcur &
                    + 2*(i-1) + 2) - (wwtens2* &
                    cosh(wwtens2*w)/sinh(wwtens2*w) &
                    - 1.0/w)/(wwtens2*wwtens2)
            enddo
         endif
         do i = 1,kwwknt
            if (kwwbdry(i) .eq. 1) then
               ncrsp = ncrsp + 1
               crsp(ncrsp,:) = 0.0
               z(ncrsp) = wwbdry(i)*darea
               crsp(ncrsp,nffcoi + kppcur + kffcur+2*i - 1) = 1.0
            endif
            if (kww2bdry(i) .eq. 1) then
               ncrsp = ncrsp + 1
               crsp(ncrsp,:) = 0.0
               z(ncrsp) = ww2bdry(i)*darea
               crsp(ncrsp,nffcoi + kppcur + kffcur+2*i) = 1.0
            endif
         enddo
      end select
      select case (kwwfnc)
      case (3,4,5,6,7)
         if (wcurbd .eq. 1.0) then
            ncrsp = ncrsp + 1
            crsp(ncrsp,:) = 0.0
            z(ncrsp) = 0.0
            do j = 1,kwwcur
               crsp(ncrsp,nffcoi+kppcur+kffcur+j) = bswwel(kwwfnc,j,1.0)
            enddo
         endif
      end select
      return
      end subroutine wwcnst

!**********************************************************************
!>
!!    Store the solution coefs into wwbdry and ww2bdry
!!    
!**********************************************************************     
      subroutine wwstore()
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      integer*4 i

      if (kwwfnc .ge. 0 .and. kwwfnc .le. 2) then
         do i = 1,kwwcur
            wwbdry(i) = brsp(nfsum+kppcur+kffcur+i)/darea
            ww2bdry(i) = 0.0
         enddo
      else if (kwwfnc .eq. 6) then
         do i = 1,kwwknt
            if(kwwbdry(i) .ne. 1) &
                wwbdry(i) = brsp(nfsum+kppcur+kffcur+2*i-1)/darea
            if(kww2bdry(i) .ne. 1) &
               ww2bdry(i) = brsp(nfsum+kppcur+kffcur+2*i)/darea
         enddo
      endif
      return
      end subroutine wwstore
