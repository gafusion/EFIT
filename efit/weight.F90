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
      if (iweigh.le.0) then
        do kk=1,nwnh
          www(kk)=0.0
          if ((xpsi(kk).lt.0.0).or.(xpsi(kk).gt.1.0)) cycle
          www(kk)=zero(kk)
        enddo
        return
      endif
!----------------------------------------------------------------------
!-- find a rectangle to search based on output from subroutine bound --
!----------------------------------------------------------------------
      do j=1,nwnh
        www(j)=0.0
      enddo
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
      if (a.ge.psil) then
        p1 = a-psil
      else
        in = in+1
        a1 = a-psil
      endif
      if (b.ge.psil) then
        p2 = b-psil
      else
        in = in+1
        a2 = b-psil
      endif
      if (c .ge. psil) then
        p3 = c-psil
      else
        in = in+1
        a3 = c-psil
      endif
      if (d.ge.psil) then
        p4 = d-psil
      else
        in = in+1
        a4 = d-psil
      endif
      in = in+1
      select case (in)
      case (1)
        www(kk) = 1.
      case (2)
        xxw = (p1+p2+p3+p4)/4.
        yyw = (a1+a2+a3+a4)
        yyw = (yyw/(xxw-yyw))**2
        www(kk) = 1.-0.5_dp*yyw
      case (3)
        xxw = (p1+p2+p3+p4)
        yyw = (a1+a2+a3+a4)
        www(kk) = xxw/(xxw-yyw)
      case (4)
        xxw = (p1+p2+p3+p4)
        yyw = (a1+a2+a3+a4)/4.
        xxw = (xxw/(xxw-yyw))**2
        www(kk) = 0.5_dp*xxw
      case (5)
        www(kk) = 0.
      end select
!     do 15 kk = 1,nwnh
      www(kk) = www(kk)*zero(kk)
!   15 continue
  320 continue
      return
      end subroutine weight
