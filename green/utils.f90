      module utils

      implicit none

      integer*4 mgaus1,mgaus2

      real*8, parameter :: pi = 3.1415926535897932, &
                          tmu = 2.0e-07

      public flux,psical,br,bz
      private lgauss,soleno,xmdele,xmdelk

      contains
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          flux computes the mutual inductance/2/pi between        **
!**          two conductors with rectangular cross sections.         **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**       r1,r2...........radius of first and second coil            **
!**       z1,z2...........elevation                                  **
!**       w1,w2...........width                                      **
!**       t1,t12,t2,t22...angle                                      **
!**                                                                  **
!**********************************************************************
      subroutine flux(r1,z1,w1,h1,t1,t12,r2,z2,w2,h2,t2,t22,fuxx)
      implicit none
      real*8, intent(in) :: r1,z1,w1,h1,t1,t12,r2,z2,w2,h2,t2,t22
      real*8, intent(out) :: fuxx
      integer*4 i,j,ner
      real*8 x,zf,zs,hf2,hs2,xl1,xr1,xbm1,xtm1,xl2,xr2,xbm2,xtm2,rf,rs, &
             hfa,hsa,solx
      real*8 post1(mgaus1),wght1(mgaus1),post2(mgaus2),wght2(mgaus2)
!      data init/0/
!
!      init=0
!      if (init.le.0) then
!vas introduced to make it ok in hydra
         call lgauss(post1,wght1,mgaus1,ner)
         call lgauss(post2,wght2,mgaus2,ner)
!         init=1
!      endif
!
      x = 0.
      fuxx = 0.
      xbm1 = 0.
      xtm1 = 0.
      xbm2 = 0.
      xtm2 = 0.
      zf = z1
      zs = z2
      hf2 = h1*.5
      hs2 = h2*.5
!
      if (t1.eq.0 .and. t12.ne.0) then
         xl1 = r1-0.5*w1-0.5*h1/abs(t12)
         xr1 = r1+0.5*w1+0.5*h1/abs(t12)
         xbm1 = xl1+w1
         xtm1 = xr1-w1
         if(t12.lt.0.) xbm1 = xr1 - w1
         if(t12.lt.0.) xtm1 = xl1 + w1
      endif
!
      if (t2.eq.0 .and. t22.ne.0) then
         xl2 = r2-0.5*w2-0.5*h2/abs(t22)
         xr2 = r2+0.5*w2+0.5*h2/abs(t22)
         xbm2 = xl2+w2
         xtm2 = xr2-w2
         if(t22.lt.0.) xbm2 = xr2 - w2
         if(t22.lt.0.) xtm2 = xl2 + w2
      endif
!
      do i = 1,mgaus1
         rf = r1+.5*w1*post1(i)
         if(t1.eq.0 .and. t12.ne.0) &
            rf = r1+(0.5*w1+0.5*h1/abs(t12))*post1(i)
         do j = 1,mgaus2
            rs = r2+0.5*w2*post2(j)
            if(t2.eq.0 .and. t22.ne.0) &
               rs = r2+(0.5*w2+0.5*h2/abs(t22))*post2(j)
            call soleno(r1,z1,w1,h1,t1,t12,r2,z2,w2,h2,t2,t22,xbm1,xbm2, &
                        xtm1,xtm2,hfa,hsa,rf,rs,solx)
            fuxx = fuxx+wght1(i)*wght2(j)/hfa/hsa*solx
         enddo 
      enddo 
!
      return
      end subroutine flux
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          lgauss computes the zeroes of the legendre polynomial   **
!**          and their associated weights for a gaussian quadrature. **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**       x...............zeroes of legendre polynomial              **
!**       w...............weights                                    **
!**       n...............order of legendre polynomial               **
!**       nn..............error flag                                 **
!**                                                                  **
!**********************************************************************
      subroutine lgauss(x,w,n,nn)
      implicit none
      integer*4, intent(in) :: n
      integer*4, intent(out) :: nn
      real*8, dimension(n), intent(out) ::  x,w
      integer*4 i,ic,k
      real*8 a,dp,g,p,r,s,su,t,test,u,v
!
      nn = 0
      if (n-1.lt.0) then 
         nn = 1
         return
      elseif (n-1.eq.0) then
!----------------------------------------------------------------------
!--      request for a zero point formula is meaningless             --
!----------------------------------------------------------------------
         x(1) = 0.
         w(1) = 2.
         return
      endif
!----------------------------------------------------------------------
!--   for a one point formula, send back                             --
!--   results without computing.                                     --
!----------------------------------------------------------------------
      r = n
      g = -1.
!----------------------------------------------------------------------
!--   the initial guess for the smallest root                        --
!--   of p(n) is taken as -1.                                        --
!----------------------------------------------------------------------
      do i = 1,n
         test = -2.
         ic = n+1-i
!----------------------------------------------------------------------
!--      whenever we find a root of the                              --
!--      polynomial, its negative is also a root.                    --
!--      the index ic tells where to store the other root            --
!----------------------------------------------------------------------
         if (ic.lt.i) return
   40    s = g
         t = 1.
         u = 1.
         v = 0.
!----------------------------------------------------------------------
!--      evaluation of the n-th legendre polynomial                  --
!--      and its first derivative.                                   --
!--      where   u = ds/dx                                           --
!--              v = dt/dx                                           --
!--              dp=dp/dx                                            --
!----------------------------------------------------------------------
         do k = 2,n
            a = k
            p = ((2.0*a-1.0)*s*g-(a-1.0)*t)/a
           dp = ((2.0*a-1.0)*(s+g*u)-(a-1.0)*v)/a
            v = u
            u = dp
            t = s
            s = p
         enddo 
         if (abs((test-g)/(test+g)).ge.0.0000005) then
            su = 0.
            if (i.ne.1) then
!----------------------------------------------------------------------
!--            the following computes the reduced                    --
!--            legendre polynomial and its derivative.               --
!----------------------------------------------------------------------
               do k = 2,i
                  su = su+1./(g-x(k-1))
               enddo
            endif
            test = g
            g = g-p/(dp-p*su)
            go to 40
         endif
         x(ic) = -g
         x(i) = g
         w(i) = 2./(r*t*dp)
         w(ic) = w(i)
         g = g-r*t/((r+2.)*g*dp+r*v-2.*r*t*su)
      enddo
      return
      end subroutine lgauss
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          soleno computes the inductance for a solenoid           **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**********************************************************************
      subroutine soleno(ra,z1,w1,h1,t1,t12,r2,z2,w2,h2,t2,t22,xbm1,xbm2, &
                        xtm1,xtm2,hfa,hsa,rf,rs,sol)
      implicit none
      real*8, intent(in) :: ra,z1,w1,h1,t1,t12,r2,z2,w2,h2,t2,t22, &
                            xbm1,xbm2,xtm1,xtm2,rf,rs
      real*8, intent(out) :: hfa,hsa,sol
      integer*4 i,j
      real*8 alpha,alphat,ambsq,beta,cay,csq,cpsq,delta,dr,drsq,dz,dzsq,ek, &
             epsi,epsit,fr,ksq,kpsq,msl,mut,pik,r,r1,r1sq,r2sq,sa,sig,sinf, &
             t,tr,trsq,z(2,2),zb1,zc1,zt1,zb2,zc2,zt2,zeta
      real*8, parameter :: rpi2=pi*0.5,rh=pi*2.0e-07,ut=rh/3.0,er=1.0e-05
!      data init/0/
!
!      init=1
!
      dr = rf-rs
      drsq = dr*dr
      tr = rf+rs
      trsq = tr*tr
      fr = 4.*rf*rs
      sol = 0.
      csq = fr/trsq
      cpsq = 1.-csq
!
      if (t1.ne.0.) then
         r = rf-ra+0.5*w1
         zc1 = z1-0.5*h1-0.5*w1*t1
         zb1 = zc1+t1*r
         zt1 = zb1+h1
!
      else
         zb1 = z1-h1/2.
         if (t12.lt.0.) then
            if (rf .lt. xbm1) zb1 = zb1 + t12*(rf-xbm1)
            zt1 = z1 + h1/2.
            if (rf .gt. xtm1) zt1 = zt1 - t12*(xtm1-rf)
!
         else
            if (rf.gt.xbm1) zb1 = zb1+t12*(rf-xbm1)
            zt1 = z1+h1/2.
            if (rf.lt.xtm1) zt1 = zt1-t12*(xtm1-rf)
         endif
      endif
!
      if (t2.ne.0.) then
         r = rs-r2+0.5*w2
         zc2 = z2-0.5*h2-0.5*w2*t2
         zb2 = zc2+t2*r
         zt2 = zb2+h2
!
      else
         zb2 = z2-h2/2.
         if (t22 .lt. 0.) then
            if (rs .lt. xbm2) zb2 = zb2 + t22*(rs-xbm2)
            zt2 = z2 + h2/2.
            if (rs .gt. xtm2) zt2 = zt2 - t22*(xtm2-rs)
!
         else
            if (rs.gt.xbm2) zb2 = zb2+t22*(rs-xbm2)
            zt2 = z2+h2/2.
            if (rs.lt.xtm2) zt2 = zt2-t22*(xtm2-rs)
         endif
      endif
!
      z(1,1) = zb1
      z(2,1) = zb2
      z(1,2) = zt1
      z(2,2) = zt2
      hfa = zt1-zb1
      hsa = zt2-zb2
!
      do i = 1,2
         do j = 1,2
            sig = -.25
            if(i .ne. j) sig = .25
            dz = z(1,i)-z(2,j)
            dzsq = dz*dz
            r2sq = dzsq+drsq
            r1sq = dzsq+trsq
            r1 = sqrt(r1sq)
            ksq = fr/r1sq
            t = 2./ksq-1.
            kpsq = 1.-ksq
            alpha = 1.
!--------------------------------------------------------------------------
!--         to avoid numerical truncation                                --
!--------------------------------------------------------------------------
            if(kpsq .lt. 1.0e-30) kpsq = 1.0e-30
            beta = sqrt(kpsq)
            if(beta .lt. 1.0e-30) beta = 1.0e-10
            if(cpsq .lt. 1.0e-30) cpsq = 1.0e-10
            delta = cpsq/beta
            epsi = csq/cpsq
            zeta = 0.
            sinf = 0.
            sa = .25
!
            sa = 2.*sa
            ambsq = (alpha-beta)*(alpha-beta)
            sinf = sinf+sa*ambsq
            alphat = alpha
            epsit = epsi
            alpha = .5*(alpha+beta)
            beta = sqrt(alphat*beta)
            epsi = (delta*epsi+zeta)/(1.+delta)
            delta = beta/4./alpha*(2.+delta+1./delta)
            zeta = .5*(epsit+zeta)
            do while (abs(delta-1.) .gt. er .or. ambsq .gt. 1.e-14)
               sa = 2.*sa
               ambsq = (alpha-beta)*(alpha-beta)
               sinf = sinf+sa*ambsq
               alphat = alpha
               epsit = epsi
               alpha = .5*(alpha+beta)
               beta = sqrt(alphat*beta)
               epsi = (delta*epsi+zeta)/(1.+delta)
               delta = beta/4./alpha*(2.+delta+1./delta)
               zeta = .5*(epsit+zeta)
            enddo
            cay = rpi2/alpha
            pik = cay*zeta
            ek = .5*cay*(ksq+sinf)
            msl = rh*dzsq*(r1*ek-drsq*pik/r1)
            if(csq==1.) msl = rh*dzsq*(r1*ek-cay*fr/r1*.5)
!
            mut = msl+ut*fr*r1*(cay-t*ek)
            sol = sol+sig*mut
         enddo
      enddo
!
      return
      end subroutine soleno
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          psical computes mutual inductance/2/pi between two      **
!**          circular filaments of radii a and r and                 **
!**          separation of z, for mks units multiply returned        **
!**          value by 2.0e-07.                                       **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**       a..............first filament radius                       **
!**       r..............second filament radius                      **
!**       z..............vertical separation                         **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1) f.w. mcclain and b.b. brown, ga technologies        **
!**              report ga-a14490 (1977).                            **
!**                                                                  **
!**********************************************************************
      real*8 function psical(a,r,z)
      implicit none
      real*8, intent(in) :: a,r,z
      real*8 den,xk,x1,cay,ee
!
      den=a*a+r*r+z*z+2.*a*r
      xk=4.*a*r/den
      x1=(a*a+r*r+z*z-2.*a*r)/den
      if(x1.lt.1.0e-10) x1=1.0e-10
      cay=xmdelk(x1)
      ee=xmdele(x1)
      psical=sqrt(den)*((1.e+00-0.5e+00*xk)*cay-ee)
      return
      end function psical
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          br computes mutual inductance/2/pi radial field         **
!**          between two circular filaments of radii a and r and     **
!**          separation of z, for mks units multiply returned        **
!**          value by 2.0e-07.                                       **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**       a..............first filament radius                       **
!**       r..............second filament radius                      **
!**       z..............vertical separation                         **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1) f.w. mcclain and b.b. brown, ga technologies        **
!**              report ga-a14490 (1977).                            **
!**                                                                  **
!**********************************************************************
      real*8 function br(a,r,z)
      implicit none
      real*8, intent(in) :: a,r,z
      real*8 den,x1,cay,ee
!
      den=a*a+r*r+z*z+2.*a*r
      x1=(a*a+r*r+z*z-2.*a*r)/den
      if(x1.lt.1.0e-10) x1=1.0e-10
      cay=xmdelk(x1)
      ee=xmdele(x1)
      br=z/(r*sqrt(den))*(-cay+(a*a+r*r+z*z)/((a-r)*(a-r)+z*z)*ee)
      return
      end function br
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          bz computes mutual inductance/2/pi vertical field       **
!**          between two circular filaments of radii a and r and     **
!**          separation of z, for mks units multiply returned        **
!**          value by 2.0e-07.                                       **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**       a..............first filament radius                       **
!**       r..............second filament radius                      **
!**       z..............vertical separation                         **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1) f.w. mcclain and b.b. brown, ga technologies        **
!**              report ga-a14490 (1977).                            **
!**                                                                  **
!**********************************************************************
      real*8 function bz(a,r,z)
      implicit none
      real*8, intent(in) :: a,r,z
      real*8 den,x1,cay,ee
!
      den=a*a+r*r+z*z+2.*a*r
      x1=(a*a+r*r+z*z-2.*a*r)/den
      if(x1.lt.1.0e-10) x1=1.0e-10
      cay=xmdelk(x1)
      ee=xmdele(x1)
      bz=(cay+(a*a-r*r-z*z)/((a-r)*(a-r)+z*z)*ee)/sqrt(den)
      return
      end function bz
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**       xmdele computes the elliptic integral e.                   **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**       xm1.............argument of elliptic integral e            **
!**                                                                  **
!**********************************************************************
      real*8 function xmdele(xm1)
      implicit none
      real*8, intent(in) :: xm1
      real*8, parameter :: a1=.44325141463,a2=.06260601220,a3=.04757383546, &
        a4=.01736506451,b1=.24998368310,b2=.09200180037,b3=.04069697526, &
        b4=.00526449639
!
      xmdele=1.0+xm1*(a1+xm1*(a2+xm1*(a3+xm1*a4))) &
                +xm1*(b1+xm1*(b2+xm1*(b3+xm1*b4)))*log(1.0/xm1)
      return
      end function xmdele
!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**       xmdelk computes the elliptic integral k.                   **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**       xm1.............argument of elliptic integral k            **
!**                                                                  **
!**********************************************************************
      real*8 function xmdelk(xm1)
      implicit none
      real*8, intent(in) :: xm1
      real*8, parameter :: a1=1.38629436112,a2=.09666344259,a3=.03590092383, &
        a4=.03742563713,a5=.01451196212,b1=.5,b2=.12498593597,b3=.06880248576,&
        b4=.03328355346,b5=.00441787012
!
      xmdelk=a1+xm1*(a2+xm1*(a3+xm1*(a4+xm1*a5))) &
                   +(b1+xm1*(b2+xm1*(b3+xm1*(b4+xm1*b5))))*log(1.0/xm1)
      return
      end function xmdelk

      end module utils
