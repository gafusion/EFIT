      module var_cecoil
      use set_kinds
      use eparmdud129,only:nsilop,nesum,magpri,nwnh,nfcoil,mbdry
      integer*4 iecurr
      real*8,dimension(nsilop,nesum) :: rsilec
      real*8,dimension(magpri,nesum) :: rmp2ec
      real*8,dimension(:,:),allocatable :: gridec
      real*8,dimension(nfcoil,nesum) :: rfcec
      real*8,dimension(nesum) :: ecurrt,pecur
      real*8,dimension(mbdry,nesum) :: rbdrec, rsolec
      real*8,dimension(nesum,nesum) :: recec
      data ecurrt(3)/-1.e10_dp/,ecurrt(4)/-1.e10_dp/, &
           ecurrt(5)/-1.e10_dp/,ecurrt(6)/-1.e10_dp/, &
           ecurrt(1)/0.0/,ecurrt(2)/0.0/
      end module var_cecoil

      module var_cvesel
      use eparmdud129,only:nsilop,magpri,nvesel,nesum,nwnh, &
                           nfcoil,mbdry
      real*8,dimension(nvesel) :: vcurrt
      real*8,dimension(nsilop,nvesel) :: rsilvs
      real*8,dimension(magpri,nvesel) :: rmp2vs
      real*8,dimension(:,:),allocatable :: gridvs
      real*8,dimension(nfcoil,nvesel) :: rfcvs
      real*8,dimension(mbdry,nvesel) :: rbdrvs, rsolvs
      real*8,dimension(nvesel,nfcoil) :: rvsfc
      real*8,dimension(nvesel,nesum) :: rvsec
      real*8,dimension(nvesel,nvesel) :: rvsvs
      end module var_cvesel

!vas      common/gtable/gridfc(nwnh,nfcoil),rgrid(nw),zgrid(nh) &
!vas           ,rsilfc(nsilop,nfcoil),rmp2fc(magpri,nfcoil) &
!vas           ,gridpc(nwnh,nw),gsilpc(nsilop,nwnh) &
!vas           ,gmp2pc(magpri,nwnh),gwork(nbwork,nwnh) &
!vas           ,rfcfc(nfcoil,nfcoil) &
!vas           ,cjfgtable(mbdry*nwnh-magpri*nwnh-nbwork*nwnh-nfcoil*nfcoil)

      module var_stable
      use eparmdud129,only:nstark,nwnh
      real*8,dimension(:,:),allocatable :: gbrpc,gbzpc
      end module var_stable

      module var_rtable
      use eparmdud129,only:nsilop,npcurn,magpri,mpress,nppcur
      real*8,dimension(nsilop,npcurn) :: rsilpc
      real*8,dimension(npcurn) :: fgowpc
      real*8,dimension(magpri,npcurn) :: rmp2pc
      real*8,dimension(mpress,nppcur) :: rprepc
      end module var_rtable

      module var_atable
      use eparmdud129,only:nsilop,nacoil,magpri,nwnh
      real*8,dimension(nsilop,nacoil) :: rsilac
      real*8,dimension(magpri,nacoil) :: rmp2ac
      real*8,dimension(:,:),allocatable :: gridac
      end module var_atable

      module var_fourier
      use eparmdud129,only:nvesel,nfourier
      real*8,dimension(nvesel) :: thetav
      real*8,dimension(nfourier,nvesel) :: sinta,costa
      real*8,dimension(2*nfourier+1,nvesel) :: vecta
      end module var_fourier
