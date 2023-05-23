!**********************************************************************
!>
!!    This function returns the matrix element for the
!!    selected basis function.
!!
!!    @param ifunc : basis function number
!!    @param iparm : basis function parameter number
!!    @param ypsi : independent variable value
!**********************************************************************
      real*8 function bsffel(ifunc,iparm,ypsi)
      
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      integer*4, intent(in) :: ifunc, iparm
      real*8, intent(in) :: ypsi
      integer*4 iorder,kffcm1,nk
      real*8 fftens2,tpsi,w
      
      bsffel = 0.0
      select case (ifunc)
      case (0)
         if (iparm.eq.1) then
            bsffel = 1.0 - ypsi**kffcur*fcurbd
         else
            bsffel = ypsi**(iparm - 1) - ypsi**kffcur*fcurbd
         endif
      case (1)
         tpsi = ypsi - 1.0
         if (iparm.eq.1) then
            bsffel = 1.0 - tpsi**kffcur*fcurbd
         else
            bsffel = tpsi**(iparm - 1) - tpsi**kffcur*fcurbd
         endif
!----------------------------------------------------------------------
!-- new local cos2 representation                                    --
!----------------------------------------------------------------------
      case (8)
         kffcm1 = kffcur - 1
         if (iparm.eq.1) then
            bsffel = 1.0 - ypsi**kffcm1*fcurbd
         elseif (iparm.eq.kffcur) then
            tpsi = ypsi - psiecn
            if (abs(tpsi).gt.dpsiecn) then
              bsffel = 0.0
            else
              bsffel = cos(rkec*tpsi)
              bsffel = bsffel**2
           endif
         else
            bsffel = ypsi**(iparm - 1) - ypsi**kffcm1*fcurbd
         endif
      case (2)
         if (iparm.eq.1) then
            bsffel = -(1.0 - ypsi**kffcur*fcurbd)
         else
            bsffel = -(ypsi**(iparm - 1) - ypsi**kffcur*fcurbd)
         endif
      case (3)
         nk = (iparm - 1) / 4 + 1
         if(nk .ge. kffknt) nk = kffknt - 1
         if ((nk .lt. (kffknt - 1) .and. ypsi .lt. ffknt(nk+1) &
              .and.  ypsi .ge. ffknt(nk)) &
              .or. (nk .eq. (kffknt - 1) .and. ypsi .le. ffknt(nk+1) &
              .and.  ypsi .ge. ffknt(nk))) then
            w = ffknt(nk+1) - ffknt(nk)
            if(mod(iparm,4) .eq. 1) bsffel = 1.0
            if(mod(iparm,4) .eq. 2) bsffel = ypsi
            if(mod(iparm,4) .eq. 3) bsffel = cos(w*fftens*ypsi)
            if(mod(iparm,4) .eq. 0) bsffel = sin(w*fftens*ypsi)
         endif
      case (4)
         nk = (iparm - 1) / 4 + 1
         if(nk .ge. kffknt) nk = kffknt - 1
         if ((nk .lt. (kffknt - 1) .and. ypsi .lt. ffknt(nk+1) &
              .and.  ypsi .ge. ffknt(nk)) &
              .or. (nk .eq. (kffknt - 1) .and. ypsi .le. ffknt(nk+1) &
              .and.  ypsi .ge. ffknt(nk))) then
            if(mod(iparm,4) .eq. 1) bsffel = 1.0
            if(mod(iparm,4) .eq. 2) bsffel = ypsi
            if(mod(iparm,4) .eq. 3) bsffel = cos(fftens*ypsi)
            if(mod(iparm,4) .eq. 0) bsffel = sin(fftens*ypsi)
         endif
      case (5)
         iorder = kffcur / (kffknt - 1)
         nk = (iparm - 1) / iorder + 1
         if(nk .ge. kffknt) nk = kffknt - 1
         if ((nk .lt. (kffknt - 1) .and. ypsi .lt. ffknt(nk+1) &
              .and.  ypsi .ge. ffknt(nk)) &
              .or. (nk .eq. (kffknt - 1) .and. ypsi .le. ffknt(nk+1) &
              .and.  ypsi .ge. ffknt(nk))) then
            w = ffknt(nk+1) - ffknt(nk)
            tpsi = (ypsi - ffknt(nk)) / w
            if (mod(iparm,iorder) .eq. 0) then
               if (iorder.eq.1) then
                  bsffel = 1.0
               else
                  bsffel = tpsi**(iorder-1)
               endif
            else
               bsffel = tpsi**(mod(iparm,iorder)-1)
            endif
         endif
      case (6)
         nk = ((iparm - 1) / 2) + 1
         fftens2 = abs(fftens)*(kffknt-1)/ &
              (ffknt(kffknt)-ffknt(1)) 
         if (nk .gt. 1) then
            if (ypsi .le. ffknt(nk) .and.  &
                 ypsi .ge. ffknt(nk-1)) then
               w = ffknt(nk) - ffknt(nk-1)
               if (mod(iparm,2) .eq. 0) then
                  bsffel = (sinh(fftens2*(ypsi-ffknt(nk-1)))/ &
                       sinh(fftens2*w) - (ypsi-ffknt(nk-1))/w) &
                       / (fftens2*fftens2)
               else
                  bsffel = (ypsi-ffknt(nk-1))/w
               endif
            endif
         endif
         if (nk .lt. kffknt) then
            if (ypsi .ge. ffknt(nk) .and.  &
                 ypsi .le. ffknt(nk+1)) then
               w = ffknt(nk+1) - ffknt(nk)
               if (mod(iparm,2) .eq. 0) then
                  bsffel = (sinh(fftens2*(ffknt(nk+1)-ypsi))/ &
                       sinh(fftens2*w) - (ffknt(nk+1)-ypsi)/w) &
                       / (fftens2*fftens2)
               else
                  bsffel = (ffknt(nk+1) - ypsi)/w
               endif
            endif
         endif
      case (7)
         if (iparm.eq.kffcur) then
            bsffel = ypsi**(kffhord)
         elseif (iparm .eq. 1) then
            bsffel = 1.0
         else
            bsffel = ypsi**(iparm - 1) 
         endif
      end select

      if(ifunc .ne. kfffnc)  &
           write(6,*)'ifunc .ne. kfffnc ',ifunc,kfffnc
      return
      end function bsffel

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
!!    @param ypsi : ndependent variable value
!!
!**********************************************************************
      real*8 function bsffpel(ifunc,iparm,ypsi)

      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      integer*4, intent(in) :: ifunc, iparm
      real*8, intent(in) :: ypsi
      integer*4 iorder,nk
      real*8 jparm,kffcm1,fftens2,tpsi,w

      bsffpel = 0.0
      select case (ifunc)
      case (0)
         if (iparm.eq.1) then
            if (kffcur.eq.1) then
               bsffpel = -fcurbd
            else
               bsffpel = -kffcur*ypsi**(kffcur-1)*fcurbd
            endif
         elseif (iparm.eq.2) then
            bsffpel = 1. - kffcur*ypsi**(kffcur-1)*fcurbd
         elseif (iparm.gt.2) then
            bsffpel = (iparm - 1)*ypsi**(iparm - 2) - &
                 kffcur*ypsi**(kffcur-1)*fcurbd
         endif
!----------------------------------------------------------------------
!-- new local cos2 representation                                    --
!----------------------------------------------------------------------
      case (8)
         kffcm1 = kffcur - 1
         if (iparm.eq.1) then
            if (kffcm1.eq.1) then
               bsffpel = -fcurbd
            else
               bsffpel = -kffcm1*ypsi**(kffcm1-1)*fcurbd
            endif
         elseif (iparm.eq.2) then
            bsffpel = 1. - kffcm1*ypsi**(kffcm1-1)*fcurbd
         elseif (iparm.eq.kffcur) then
            tpsi = ypsi - psiecn
            if (abs(tpsi).gt.dpsiecn) then
              bsffpel = 0.0
            else
              bsffpel = 2.0*(rkec*tpsi)
              bsffpel = -rkec*sin(bsffpel)
            endif
         elseif (iparm.gt.2) then
            bsffpel = (iparm - 1)*ypsi**(iparm - 2) - &
                 kffcm1*ypsi**(kffcm1-1)*fcurbd
         endif
      case (1)
         tpsi = ypsi - 1.0
         if (iparm.eq.1) then
            if (kffcur.eq.1) then
               bsffpel = -fcurbd
            else
               bsffpel = -kffcur*tpsi**(kffcur-1)*fcurbd
            endif
         elseif (iparm.eq.2) then
            bsffpel = 1. - kffcur*tpsi**(kffcur-1)*fcurbd
         elseif (iparm.gt.2) then
            bsffpel = (iparm - 1)*tpsi**(iparm - 2) - &
                 kffcur*tpsi**(kffcur-1)*fcurbd
         endif
      case (2)
         if (iparm.eq.1) then
            if (kffcur.eq.1) then
               bsffpel = -fcurbd
            else
               bsffpel = -kffcur*ypsi**(kffcur-1)*fcurbd
            endif
         elseif (iparm.eq.2) then
            bsffpel = 1. - kffcur*ypsi**(kffcur-1)*fcurbd
         elseif (iparm.gt.2) then
            bsffpel = (iparm - 1)*ypsi**(iparm - 2) - &
                 kffcur*ypsi**(kffcur-1)*fcurbd
         endif
         bsffpel = - bsffpel
      case (3)
         nk = (iparm - 1) / 4 + 1
         if(nk .ge. kffknt) nk = kffknt - 1
         if ((nk .lt. (kffknt - 1) .and. ypsi .lt. ffknt(nk+1) &
              .and.  ypsi .ge. ffknt(nk)) &
              .or. (nk .eq. (kffknt - 1) .and. ypsi .le. ffknt(nk+1) &
              .and.  ypsi .ge. ffknt(nk))) then
            w = ffknt(nk+1) - ffknt(nk)
            if(mod(iparm,4) .eq. 1) bsffpel = 0.0
            if(mod(iparm,4) .eq. 2) bsffpel = 1.0
            if(mod(iparm,4) .eq. 3) bsffpel =  &
                 -w*fftens*sin(w*fftens*ypsi)
            if(mod(iparm,4) .eq. 0) bsffpel =  &
                 w*fftens*cos(w*fftens*ypsi)
         endif
      case (4)
         nk = (iparm - 1) / 4 + 1
         if(nk .ge. kffknt) nk = kffknt - 1
         if ((nk .lt. (kffknt - 1) .and. ypsi .lt. ffknt(nk+1) &
              .and.  ypsi .ge. ffknt(nk)) &
              .or. (nk .eq. (kffknt - 1) .and. ypsi .le. ffknt(nk+1) &
              .and.  ypsi .ge. ffknt(nk))) then
            if(mod(iparm,4) .eq. 1) bsffpel = 0.0
            if(mod(iparm,4) .eq. 2) bsffpel = 1.0
            if(mod(iparm,4) .eq. 3) bsffpel = -fftens*sin(fftens*ypsi)
            if(mod(iparm,4) .eq. 0) bsffpel = fftens*cos(fftens*ypsi)
         endif
      case (5)
        iorder = kffcur / (kffknt - 1)
        nk = (iparm - 1) / iorder + 1
        if(nk .ge. kffknt) nk = kffknt - 1
        if ((nk .lt. (kffknt - 1) .and. ypsi .lt. ffknt(nk+1) &
          .and.  ypsi .ge. ffknt(nk)) &
          .or. (nk .eq. (kffknt - 1) .and. ypsi .le. ffknt(nk+1) &
          .and.  ypsi .ge. ffknt(nk))) then
          w = ffknt(nk+1) - ffknt(nk)
          tpsi = (ypsi - ffknt(nk)) / w
          jparm = mod(iparm,iorder)
          if (jparm.eq.1) then
            bsffpel = 0.0
          elseif (jparm.eq.2) then
            bsffpel = 1./w
          elseif (jparm.gt.2) then
            bsffpel = (jparm - 1)/w*tpsi**(jparm - 2)
          endif
        endif
      case (6)
        nk = ((iparm - 1) / 2) + 1
        fftens2 = abs(fftens)*(kffknt-1)/ &
          (ffknt(kffknt)-ffknt(1))
        if (nk .gt. 1) then
          if (ypsi .le. ffknt(nk) .and.  &
            ypsi .ge. ffknt(nk-1)) then
            w = ffknt(nk) - ffknt(nk-1)
            if (mod(iparm,2) .eq. 0) then
              bsffpel = (fftens2*cosh(fftens2* &
                (ypsi-ffknt(nk-1)))/sinh(fftens2*w) - (1.0/w)) &
                / (fftens2*fftens2)
            else
              bsffpel = 1.0/w
            endif
          endif
        endif
        if (nk .lt. kffknt) then
          if (ypsi .ge. ffknt(nk) .and.  &
            ypsi .le. ffknt(nk+1)) then
            w = ffknt(nk+1) - ffknt(nk)
            if (mod(iparm,2) .eq. 0) then
              bsffpel = (-fftens2*cosh(fftens2* &
                (ffknt(nk+1)-ypsi))/sinh(fftens2*w)+(1.0/w)) &
                / (fftens2*fftens2)
            else
              bsffpel = -1.0/w
            endif
          endif
        endif
      case (7)
        if (iparm.eq.kffcur) then
          bsffpel = kffhord*ypsi**(kffhord-1)
        elseif (iparm.eq.1) then
          bsffpel = 0.
        elseif (iparm.eq.2) then
          bsffpel = 1.
        elseif (iparm.gt.2) then
          bsffpel = (iparm - 1)*ypsi**(iparm - 2)
        endif
      end select
      return
      end function bsffpel

!**********************************************************************
!>
!!    In addition to the least squares constraints that
!!    efit already uses, some basis functions have exact
!!    constraints. Most notable is the spline function
!!    whose continuity constraints are exact, not LSE.
!!
!!    @param ncrsp : number of constraint equations
!!
!!    @param crsp :  constraint matrix
!!
!!    @param z : value vector
!!
!********************************************************************** 
      real*8 function bsffin(ifunc,iparm,ypsi)
      
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      integer*4, intent(in) :: ifunc, iparm
      real*8, intent(in) :: ypsi
      integer*4 iorder,nk
      real*8 b1,b2,fftens2,tpsi,tpsi2,w,ypsi1,ypsi2
      
      bsffin = 0.0
      ypsi2 = 1.0
      select case (ifunc)
      case (0)
         bsffin = (ypsi**iparm)/iparm - &
              (ypsi**(kffcur+1))/(kffcur+1)*fcurbd                
         bsffin = bsffin - ((ypsi2**iparm)/iparm - &
              (ypsi2**(kffcur+1))/(kffcur+1)*fcurbd)
      case (8)
         if (iparm.ne.kffcur) then
         bsffin = (ypsi**iparm)/iparm - &
              (ypsi**kffcur)/kffcur*fcurbd
         bsffin = bsffin - ((ypsi2**iparm)/iparm - &
              (ypsi2**kffcur)/kffcur*fcurbd)
         else
            tpsi = ypsi - psiecn
            if (tpsi.lt.-dpsiecn) then
              bsffin = 0.0
            elseif (tpsi.gt.dpsiecn) then
              bsffin = dpsiecn
            else
              bsffin = 2.0*rkec*tpsi
              bsffin = tpsi + dpsiecn + sin(bsffin) /2.0/rkec
              bsffin = bsffin/2.0
            endif
            bsffin = bsffin - dpsiecn
         endif
      case (1)
         tpsi = ypsi - 1.0
         tpsi2 = ypsi2 - 1.0
         bsffin = (tpsi**iparm)/iparm - &
              (tpsi**(kffcur+1))/(kffcur+1)*fcurbd                
         bsffin = bsffin - ((tpsi2**iparm)/iparm - &
              (tpsi2**(kffcur+1))/(kffcur+1)*fcurbd)         
      case (2)
         bsffin = -((ypsi**iparm)/iparm - &
              (ypsi**(kffcur+1))/(kffcur+1)*fcurbd)
         bsffin = bsffin - (-((ypsi2**iparm)/iparm - &
              (ypsi2**(kffcur+1))/(kffcur+1)*fcurbd))
      case (3)
         nk = (iparm - 1) / 4 + 1
         if(nk .ge. kffknt) nk = kffknt - 1
         if (ypsi .ge. ffknt(nk+1)) then
            bsffin = 0
            return
         endif
         if (1.0 .le. ffknt(nk)) then
            bsffin = 0
            return
         endif
         if (ypsi .ge. ffknt(nk)) then
            ypsi1 = ypsi
         else
            ypsi1 = ffknt(nk)
         endif
         w = ffknt(nk+1) - ffknt(nk)
         if (1.0 .ge. ffknt(nk+1)) then
            ypsi2 = ffknt(nk+1)
         else
            ypsi2 = 1.0
         endif
         if(mod(iparm,4) .eq. 1) b1 = ypsi1
         if(mod(iparm,4) .eq. 2) b1 = (ypsi1**2) / 2.0
         if(mod(iparm,4) .eq. 3) &
              b1 = sin(w*fftens*ypsi1)/w*fftens
         if(mod(iparm,4) .eq. 0) &
              b1 = -cos(w*fftens*ypsi1)/w*fftens
         if(mod(iparm,4) .eq. 1) b2 = ypsi2
         if(mod(iparm,4) .eq. 2) b2 = (ypsi2**2) / 2.0
         if(mod(iparm,4) .eq. 3) &
              b2 = sin(w*fftens*ypsi2)/w*fftens
         if(mod(iparm,4) .eq. 0) &
              b2 = -cos(w*fftens*ypsi2)/w*fftens
         bsffin = b1 - b2
         !  write(6,*)'for ypsi=',ypsi,' integrate from ',ypsi1,' to ',
         !  $                       ypsi2,' = ',bsffin
      case (4)
         nk = (iparm - 1) / 4 + 1
         if(nk .ge. kffknt) nk = kffknt - 1
         if (ypsi .ge. ffknt(nk+1)) then
            bsffin = 0
            return
         endif
         if (1.0 .le. ffknt(nk)) then
            bsffin = 0
            return
         endif
         if (ypsi .ge. ffknt(nk)) then
            ypsi1 = ypsi
         else
            ypsi1 = ffknt(nk)
         endif
         if (1.0 .ge. ffknt(nk+1)) then
            ypsi2 = ffknt(nk+1)
         else
            ypsi2 = 1.0
         endif
         if(mod(iparm,4) .eq. 1) b1 = ypsi1
         if(mod(iparm,4) .eq. 2) b1 = (ypsi1**2) / 2.0
         if(mod(iparm,4) .eq. 3) &
              b1 = sin(fftens*ypsi1)/fftens
         if(mod(iparm,4) .eq. 0) &
              b1 = -cos(fftens*ypsi1)/fftens
         if(mod(iparm,4) .eq. 1) b2 = ypsi2
         if(mod(iparm,4) .eq. 2) b2 = (ypsi2**2) / 2.0
         if(mod(iparm,4) .eq. 3) &
              b2 = sin(fftens*ypsi2)/fftens
         if(mod(iparm,4) .eq. 0) &
              b2 = -cos(fftens*ypsi2)/fftens
         bsffin = b1 - b2
         !     write(6,*)'for ypsi=',ypsi,' integrate from ',ypsi1,' to ',
         !     $                       ypsi2,' = ',bsffin
      case (5)
         iorder = kffcur / (kffknt - 1)
         nk = (iparm - 1) / iorder + 1
         if(nk .ge. kffknt) nk = kffknt - 1
         if (ypsi .ge. ffknt(nk+1)) then
            bsffin = 0
            return
         endif
         if (1.0 .le. ffknt(nk)) then
            bsffin = 0
            return
         endif
         if (ypsi .ge. ffknt(nk)) then
            ypsi1 = ypsi
         else
            ypsi1 = ffknt(nk)
         endif
         if (1.0 .ge. ffknt(nk+1)) then
            ypsi2 = ffknt(nk+1)
         else
            ypsi2 = 1.0
         endif
         w = ffknt(nk+1) - ffknt(nk)
         
         tpsi=(ypsi1**2/2.0 - ypsi1*ffknt(nk))/2
         if (mod(iparm,iorder) .eq. 0) then
            b1 = tpsi**iorder / iorder
         else
            b1 = tpsi**mod(iparm,iorder) / mod(iparm,iorder)
         endif
         tpsi=(ypsi2**2/2.0 - ypsi2*ffknt(nk))/2
         if (mod(iparm,iorder) .eq. 0) then
            b2 = tpsi**iorder / iorder
         else
            b2 = tpsi**mod(iparm,iorder) / mod(iparm,iorder)
         endif
         bsffin = b1 - b2
         !     write(6,*)'for ypsi=',ypsi,' integrate from ',ypsi1,' to ',
         !     $                       ypsi2,' = ',bsffin
      case (6)
         nk = ((iparm - 1) / 2) + 1
         fftens2 = abs(fftens)*(kffknt-1)/ &
              (ffknt(kffknt)-ffknt(1))
         bsffin = 0
         if (nk .gt.1) then
            if (ypsi .le. ffknt(nk)) then
               if (ypsi .le. ffknt(nk-1)) then
                  ypsi1 = ffknt(nk-1)
               else
                  ypsi1 = ypsi
               endif
               if (1.0 .le. ffknt(nk)) then
                  ypsi2 = 1.0
               else
                  ypsi2 = ffknt(nk)
               endif
               w = ffknt(nk) - ffknt(nk-1)
               if (mod(iparm,2) .eq. 0) then
                  b1 = (cosh(fftens2*(ypsi1-ffknt(nk-1)))/ &
                       (fftens2*sinh(fftens2*w)) - (ypsi1*ypsi1/2.0- &
                       ffknt(nk-1)*ypsi1)/w) &
                       / (fftens2*fftens2)
               else
                  b1 = (ypsi1*ypsi1/2.0-ffknt(nk-1)*ypsi1)/w
               endif
               if (mod(iparm,2) .eq. 0) then
                  b2 = (cosh(fftens2*(ypsi2-ffknt(nk-1)))/ &
                       (fftens2*sinh(fftens2*w)) - (ypsi2*ypsi2/2.0- &
                       ffknt(nk-1)*ypsi2)/w) &
                       / (fftens2*fftens2)
               else
                  b2 = (ypsi2*ypsi2/2.0-ffknt(nk-1)*ypsi2)/w
               endif
               bsffin = bsffin + b1 - b2
            endif
         endif
         if (nk .lt. kffknt) then
            if (ypsi .le. ffknt(nk+1)) then
               if (ypsi .le. ffknt(nk)) then
                  ypsi1 = ffknt(nk)
               else
                  ypsi1 = ypsi
               endif
               if (1.0 .le. ffknt(nk+1)) then
                  ypsi2 = 1.0
               else
                  ypsi2 = ffknt(nk+1)
               endif
               w = ffknt(nk+1) - ffknt(nk)
               if (mod(iparm,2) .eq. 0) then
                  b1 = (-cosh(fftens2*(ffknt(nk+1)-ypsi1))/ &
                       (fftens2*sinh(fftens2*w))-(ypsi1*ffknt(nk+1) &
                       -ypsi1*ypsi1/2.0)/w) &
                       / (fftens2*fftens2)
               else
                  b1 = (ffknt(nk+1)*ypsi1-ypsi1*ypsi1/2.0)/w
               endif
               if (mod(iparm,2) .eq. 0) then
                  b2 = (-cosh(fftens2*(ffknt(nk+1)-ypsi2))/ &
                       (fftens2*sinh(fftens2*w))-(ypsi2*ffknt(nk+1) &
                       -ypsi2*ypsi2/2.0)/w) &
                       / (fftens2*fftens2)
               else
                  b2 = (ffknt(nk+1)*ypsi2-ypsi2*ypsi2/2.0)/w
               endif
               bsffin = bsffin + b1 - b2
            endif
         endif
      case (7)
        if (iparm .eq. kffcur) then
          bsffin = (ypsi**(kffhord+1))/(kffhord+1)
          bsffin = bsffin - ((ypsi2**(kffhord+1))/(kffhord+1))
        else
          bsffin = (ypsi**iparm)/iparm
          bsffin = bsffin - ((ypsi2**iparm)/iparm)
        endif
      end select

      if(ifunc .ne. kfffnc) write(6,*)'ifunc .ne. kfffnc ',ifunc,kfffnc
      return
      end function bsffin

!**********************************************************************
!>
!!    In addition to the least squares constraints that
!!    efit already uses, some basis functions have exact
!!    constraints. Most notable is the spline function
!!    whose continuity constraints are exact, not LSE.
!!
!!    @param ncrsp : number of constraint equations
!!
!!    @param crsp :  constraint matrix
!!
!!    @param z : value vector
!!
!!    @param nffcoi : array index for setting up crsp
!!
!**********************************************************************  
      subroutine ffcnst(ncrsp,crsp,z,nffcoi)
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      real*8 bsffel
      integer*4, intent(in) :: nffcoi
      integer*4, intent(inout) :: ncrsp
      real*8, intent(out) :: crsp(4*(npcurn-2)+6+npcurn*npcurn,nrsmat), &
                             z(3*(npcurn-2)+6+npcurn*npcurn)
      integer*4 i,j,iorder
      real*8 h,w,fftens2

      select case (kfffnc)
      case (3)
         if (kffknt .gt. 2) then
!     
!     first set of constraints is that splines must be equal at the knots
!     
            do i = 2,kffknt-1
               ncrsp = ncrsp + 1
               crsp(ncrsp,:) = 0.0
               z(ncrsp) = 0.0
               h = ffknt(i) - ffknt(i-1)
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 1) = 1.0
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 2) = ffknt(i)
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 3) = &
                    cos(h * fftens * ffknt(i))
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 4) = &
                    sin(h * fftens * ffknt(i))
               
               h = ffknt(i+1) - ffknt(i)
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 5) = -1.0
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 6) = -ffknt(i)
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 7) = &
                    -cos(h * fftens * ffknt(i))
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 8) = &
                    -sin(h * fftens * ffknt(i))
            enddo
!     
!     second set of constraints is that splines have equal first 
!     derivative at the knots
!     
            do i = 2,kffknt-1
               ncrsp = ncrsp + 1
               crsp(ncrsp,:) = 0.0
               z(ncrsp) = 0.0
               h = ffknt(i) - ffknt(i-1)
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 1) = 0.0
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 2) = 1.0
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 3) = &
                    -h * fftens * sin(h * fftens * ffknt(i))
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 4) = &
                    h * fftens * cos(h * fftens * ffknt(i))
               
               h = ffknt(i+1) - ffknt(i)
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 5) = 0.0
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 6) = -1.0
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 7) = &
                    h * fftens * sin(h * fftens * ffknt(i))
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 8) = &
                    -h * fftens * cos(h * fftens * ffknt(i))
            enddo
!     
!     second set of constraints is that splines have equal second 
!     derivative at the knots
!     
            do i = 2,kffknt-1
               ncrsp = ncrsp + 1
               crsp(ncrsp,:) = 0.0
               z(ncrsp) = 0.0
               h = ffknt(i) - ffknt(i-1)
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 1) = 0.0
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 2) = 0.0
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 3) = &
                    -h*h*fftens*fftens*cos(h * fftens * ffknt(i))
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 4) = &
                    -h*h*fftens*fftens*sin(h * fftens * ffknt(i))
               
               h = ffknt(i+1) - ffknt(i)
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 5) = 0.0
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 6) = 0.0
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 7) = &
                    h*h*fftens*fftens*cos(h * fftens * ffknt(i))
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 8) = &
                    h*h*fftens*fftens*sin(h * fftens * ffknt(i))
            enddo
            
         endif
      case (4)
         if (kffknt .le. 2) then
!     
!     first set of constraints is that splines must be equal at the knots
!     
            do i = 2,kffknt-1
               ncrsp = ncrsp + 1
               crsp(ncrsp,:) = 0.0
               z(ncrsp) = 0.0
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 1) = 1.0
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 2) = ffknt(i)
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 3) = &
                    cos(fftens * ffknt(i))
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 4) = &
                    sin(fftens * ffknt(i))
               
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 5) = -1.0
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 6) = -ffknt(i)
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 7) = &
                    -cos(fftens * ffknt(i))
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 8) = &
                    -sin(fftens * ffknt(i))
            enddo
!     
!     second set of constraints is that splines have equal first 
!     derivative at the knots
!     
            do i = 2,kffknt-1
               ncrsp = ncrsp + 1
               crsp(ncrsp,:) = 0.0
               z(ncrsp) = 0.0
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 1) = 0.0
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 2) = 1.0
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 3) = &
                    -fftens * sin(fftens * ffknt(i))
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 4) = &
                    fftens * cos(fftens * ffknt(i))
               
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 5) = 0.0
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 6) = -1.0
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 7) = &
                    fftens * sin(fftens * ffknt(i))
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 8) = &
                    -fftens * cos(fftens * ffknt(i))
            enddo
!     
!     second set of constraints is that splines have equal second 
!     derivative at the knots
!     
            do i = 2,kffknt-1
               ncrsp = ncrsp + 1
               crsp(ncrsp,:) = 0.0
               z(ncrsp) = 0.0
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 1) = 0.0
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 2) = 0.0
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 3) = &
                    -fftens*fftens*cos(fftens * ffknt(i))
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 4) = &
                    -fftens*fftens*sin(fftens * ffknt(i))
               
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 5) = 0.0
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 6) = 0.0
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 7) = &
                    fftens*fftens*cos(fftens * ffknt(i))
               crsp(ncrsp,nffcoi + kppcur + 4*(i-2) + 8) = &
                    fftens*fftens*sin(fftens * ffknt(i))
            enddo
            
         endif
      case (5)
         iorder = kffcur / (kffknt - 1)
         if (kffknt .le. 2) then
!     
!     first set of constraints is that splines must be equal at the knots
!     
            do i = 2,kffknt-1
               ncrsp = ncrsp + 1
               crsp(ncrsp,:) = 0.0
               z(ncrsp) = 0.0
               do j= 1,iorder
                  crsp(ncrsp,nffcoi + kppcur + iorder*(i-2) + j)  = 1.0
               enddo
               crsp(ncrsp,nffcoi + kppcur + iorder*(i-1) + 1)  = -1.0
            enddo
!     
!     second set of constraints is that splines have equal first 
!     derivative at the knots
!     
            do i = 2,kffknt-1
               ncrsp = ncrsp + 1
               crsp(ncrsp,:) = 0.0
               z(ncrsp) = 0.0
               do j= 2,iorder
                  crsp(ncrsp,nffcoi + kppcur + &
                       iorder*(i-2) + j)  = (j-1)
               enddo
               crsp(ncrsp,nffcoi + kppcur + &
                    iorder*(i-1) + 2)  = -1.0
            enddo
!     
!     second set of constraints is that splines have equal second 
!     derivative at the knots
!     
            do i = 2,kffknt-1
               ncrsp = ncrsp + 1
               crsp(ncrsp,:) = 0.0
               z(ncrsp) = 0.0
               do j= 3,iorder
                  crsp(ncrsp,nffcoi + kppcur + &
                       iorder*(i-2) + j)  = (j-1)*(j-2)
               enddo
               crsp(ncrsp,nffcoi + kppcur + &
                    iorder*(i-1) + 3)  = -2.0
            enddo
            
         endif
      case (6)
        !
        !     first set of constraints is that splines have equal first
        !     derivative at the knots
        !
        if(kffknt .gt. 2)then
          fftens2 = abs(fftens)*(kffknt-1)/ &
            (ffknt(kffknt)-ffknt(1))
          do i = 2,kffknt-1
            ncrsp = ncrsp + 1
            crsp(ncrsp,:) = 0.0
            z(ncrsp) = 0.0
            w = ffknt(i+1) - ffknt(i)
            crsp(ncrsp,nffcoi + kppcur + &
              2*(i-1) + 1) = -1.0/w
            crsp(ncrsp,nffcoi + kppcur + &
              2*(i-1) + 2) = (-fftens2* &
              cosh(fftens2*w)/sinh(fftens2*w) &
              + 1.0/w)/(fftens2*fftens2)
            crsp(ncrsp,nffcoi + kppcur &
              + 2*(i-1) + 3) = 1.0/w
            crsp(ncrsp,nffcoi + kppcur &
              + 2*(i-1) + 4) = (fftens2/ &
              sinh(fftens2*w) - 1.0/w)/(fftens2*fftens2)

            w = ffknt(i) - ffknt(i-1)
            crsp(ncrsp,nffcoi + kppcur &
              + 2*(i-1) - 1) = 1.0/w
            crsp(ncrsp,nffcoi + kppcur &
              + 2*(i-1) + 0) = -(-fftens2/ &
              sinh(fftens2*w) + 1.0/w)/(fftens2*fftens2)
            crsp(ncrsp,nffcoi + kppcur + 2*(i-1) + 1) = &
              crsp(ncrsp,nffcoi + kppcur &
              + 2*(i-1) + 1) - 1.0/w
            crsp(ncrsp,nffcoi + kppcur + 2*(i-1) + 2) = &
              crsp(ncrsp,nffcoi + kppcur &
              + 2*(i-1) + 2) - (fftens2* &
              cosh(fftens2*w)/sinh(fftens2*w) &
              - 1.0/w)/(fftens2*fftens2)
          enddo
        endif
        do i = 1,kffknt
          if (kffbdry(i) .eq. 1) then
            ncrsp = ncrsp + 1
            crsp(ncrsp,:) = 0.0
            z(ncrsp) = ffbdry(i)*darea/twopi/tmu
            crsp(ncrsp,nffcoi + kppcur+2*i - 1) = 1.0
          endif
          if (kff2bdry(i) .eq. 1) then
            ncrsp = ncrsp + 1
            crsp(ncrsp,:) = 0.0
            z(ncrsp) = ff2bdry(i)*darea/twopi/tmu
            crsp(ncrsp,nffcoi + kppcur+2*i) = 1.0
          endif
        enddo
      end select
      select case (kfffnc)
      case (3,4,5,6,7)
        if (fcurbd .eq. 1.0) then
          ncrsp = ncrsp + 1
          crsp(ncrsp,:) = 0.0
          z(ncrsp) = 0.0
          do j = 1,kffcur
            crsp(ncrsp,nffcoi+kppcur+j) = bsffel(kfffnc,j,1.0)
          enddo
        endif
      end select
      return
      end subroutine ffcnst
      
!**********************************************************************
!>
!!    This subroutine stores the solution coefs into ffbdry and ff2bdry
!!
!********************************************************************** 
      subroutine ffstore()
      include 'eparm.inc'
      include 'modules2.inc'
      include 'modules1.inc'
      implicit none
      integer*4 i

      if (kfffnc .ge. 0 .and. kfffnc .le. 2) then
         do i = 1,kffcur
            ffbdry(i) = brsp(nfsum+kppcur+i)*twopi*tmu/darea
            ff2bdry(i) = 0.0
         enddo
      elseif (kfffnc .eq. 6) then
         do i = 1,kffknt
            if(kffbdry(i) .ne. 1) &
               ffbdry(i) = brsp(nfsum+kppcur+2*i - 1)*twopi*tmu/darea
            if(kff2bdry(i) .ne. 1) &
              ff2bdry(i) = brsp(nfsum+kppcur+2*i)*twopi*tmu/darea
         enddo
      endif
      return
      end subroutine ffstore
