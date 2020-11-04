      function psical(a1,r1,z1)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          psical computes mutual inductance/2/pi between two      **
!**          circular filaments of radii a1 and r1 and               **
!**          separation of z1, for mks units multiply returned       **
!**          value by 2.0e-07.                                       **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**       a1..............first filament radius                      **
!**       r1..............second filament radius                     **
!**       z1..............vertical separation                        **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1) f.w. mcclain and b.b. brown, ga technologies        **
!**              report ga-a14490 (1977).                            **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use set_kinds
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      real*8 x1,cay,ee,xmdelk,xmdele
!
      isw=1
      go to 10
      entry br(a1,r1,z1)
      isw=2
      go to 10
      entry bz(a1,r1,z1)
      isw=3
!
   10 continue
      a=a1
      r=r1
      z=z1
      den=a*a+r*r+z*z+2.*a*r
      xk=4.*a*r/den
      x1=(a*a+r*r+z*z-2.*a*r)/den
      if (x1.lt.1.0e-10_dp) x1=1.0e-10_dp
      cay=xmdelk(x1)
      ee=xmdele(x1)
      select case (isw)
      case (1)
        go to 20
      case (2)
        go to 30
      case (3)
        go to 40
      end select
!----------------------------------------------------------------------
!--   psi computation                                                --
!----------------------------------------------------------------------
   20 psical= sqrt(den)*((1.0-0.5_dp*xk)*cay-ee)
      return
!----------------------------------------------------------------------
!--   br  computation                                                --
!----------------------------------------------------------------------
   30 psical=z/(r* sqrt(den))*(-cay+(a*a+r*r+z*z)/((a-r)*(a-r)+z*z)*ee)
      return
!----------------------------------------------------------------------
!--   bz  computation                                                --
!----------------------------------------------------------------------
   40 psical=(cay+(a*a-r*r-z*z)/((a-r)*(a-r)+z*z)*ee)/ sqrt(den)
      return
      end
      function xmdele(xm1)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**       xmdele computes the elliptic integral e.                   **
!**                                                                  **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**       xm1.............argument of elliptic integral e            **
!**                                                                  **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use set_kinds
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension a(4),b(4)
      real*8 a,b,xm1,xmdele
      data a(1),a(2),a(3),a(4)/.44325141463_dp,.06260601220_dp, &
         .04757383546_dp,.01736506451_dp/
      data b(1),b(2),b(3),b(4)/.24998368310_dp,.09200180037_dp, &
         .04069697526_dp,.00526449639_dp/
!
      xmdele=1.0+xm1*(a(1)+xm1*(a(2)+xm1*(a(3)+xm1*a(4)))) &
        +xm1*(b(1)+xm1*(b(2)+xm1*(b(3)+xm1*b(4))))*log(1.0/xm1)
      return
      end
      function xmdelk(xm1)
!**********************************************************************
!**                                                                  **
!**     MAIN PROGRAM:  MHD FITTING CODE                              **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**       xmdelk computes the elliptic integral k.                   **
!**                                                                  **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**       xm1.............argument of elliptic integral k            **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**          (2)                                                     **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          26/04/83..........first created                         **
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
      use set_kinds
      implicit integer*4 (i-n), real*8 (a-h, o-z)
      dimension a(5),b(5)
      real*8  a,b,xm1,xmdelk
      data a(1),a(2),a(3),a(4),a(5)/1.38629436112_dp,.09666344259_dp, &
         .03590092383_dp,.03742563713_dp,.01451196212_dp/
      data b(1),b(2),b(3),b(4),b(5)/.5_dp,.12498593597_dp,.06880248576_dp, &
         .03328355346_dp,.00441787012_dp/
!
      xmdelk=a(1)+xm1*(a(2)+xm1*(a(3)+xm1*(a(4)+xm1*a(5)))) &
        +(b(1)+xm1*(b(2)+xm1*(b(3)+xm1*(b(4)+xm1*b(5))))) &
        *log(1.0/xm1)
      return
      end
