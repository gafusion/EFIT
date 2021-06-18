!**********************************************************************
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          weight computes the weighting function w.               **
!**          boy, you're going to carry that weight                  **
!**          carry that weight a long time                           **
!**                                                                  **
!**     CALLING ARGUMENTS:                                           **
!**                                                                  **
!**     RECORD OF MODIFICATION:                                      **
!**          21/10/83..........first created                         **
!**          24/07/85..........revised                               **
!**                                                                  **
!**********************************************************************
      subroutine weight(x,y)
      include 'eparmdud129.inc'
      include 'modules1.inc'
      implicit integer*4 (i-n), real*8 (a-h,o-z)
      dimension x(1),y(1)
!
      if (iweigh.le.0) go to 5000
!----------------------------------------------------------------------
!-- find a rectangle to search based on output from subroutine bound --
!----------------------------------------------------------------------
      do 100 j=1,nwnh
  100 www(j)=0.0
      do 110 i=1,nw-1
      ileft=i
  110 if((x(i).le.xmin).and.(x(i+1).gt.xmin))go to 120
  120 do 130 i=ileft,nw-1
      iright=i
  130 if((x(i).lt.xmax).and.(x(i+1).ge.xmax))go to 140
  140 do 150 j=1,nh-1
      jbotm=j
  150 if((y(j).le.ymin).and.(y(j+1).gt.ymin))go to 160
  160 do 170 j=jbotm,nh-1
      jtop=j
  170 if((y(j).lt.ymax).and.(y(j+1).ge.ymax))go to 180
  180 continue 
      jtop=min(jtop,nh)
      jbotm=max(jbotm,1)
      ileft=max(ileft,1)
      iright=min(iright,nw)
      do 320 i = ileft,iright
      do 320 j = jbotm,jtop
      kk = j+nh*(i-1)
      a = psi(kk-nh)
      b = psi(kk+nh)
      c = psi(kk-1)
      d = psi(kk+1)
      if (i.eq.1) a = 0.
      if (i.eq.nw) b = 0.
      if (j.eq.1) c = 0.
      if (j.eq.nh) d = 0.
      psil = psibry
      in = 0
      p1 = 0.
      p2 = 0.
      p3 = 0.
      p4 = 0.
      a1 = 0.
      a2 = 0.
      a3 = 0.
      a4 = 0.
      if (a.ge.psil) go to 1
      in = in+1
      a1 = a-psil
      go to 2
    1 p1 = a-psil
    2 if (b.ge.psil) go to 3
      in = in+1
      a2 = b-psil
      go to 4
    3 p2 = b-psil
    4 if (c .ge. psil) go to 5
      in = in+1
      a3 = c-psil
      go to 6
    5 p3 = c-psil
    6 if (d.ge.psil) go to 7
      in = in+1
      a4 = d-psil
      go to 8
    7 p4 = d-psil
    8 in = in+1
      select case (in)
      case (1)
        go to 10
      case (2)
        go to 11
      case (3)
        go to 12
      case (4)
        go to 13
      case (5)
        go to 14
      end select
   10 www(kk) = 1.
      go to 20
   11 xxw = (p1+p2+p3+p4)/4.
      yyw = (a1+a2+a3+a4)
      yyw = (yyw/(xxw-yyw))**2
      www(kk) = 1.-0.5_dp*yyw
      go to 20
   12 xxw = (p1+p2+p3+p4)
      yyw = (a1+a2+a3+a4)
      www(kk) = xxw/(xxw-yyw)
      go to 20
   13 xxw = (p1+p2+p3+p4)
      yyw = (a1+a2+a3+a4)/4.
      xxw = (xxw/(xxw-yyw))**2
      www(kk) = 0.5_dp*xxw
      go to 20
   14 www(kk) = 0.
   20 continue
!     do 15 kk = 1,nwnh
      www(kk) = www(kk)*zero(kk)
   15 continue
  320 continue
      return
!
 5000 continue
      do 5500 kk=1,nwnh
        www(kk)=0.0
        if ((xpsi(kk).lt.0.0).or.(xpsi(kk).gt.1.0)) go to 5500
        www(kk)=zero(kk)
 5500 continue
      return
      end subroutine weight
