#include "config.f"
!**********************************************************************
!>
!!    sets default DIII-D values for eparm module
!!
!**********************************************************************
      subroutine set_eparm_defaults()
      use eparm
      use var_gwork1, only: mfila
      use var_exdata, only: ifitvs
      use var_cecoil, only: iecurr
      use var_vessel, only: ivesel
      use var_input, only: icutfp
      use errlims
      implicit none
      integer*4 i,j,nhpwr

      device="DIII-D"
      nsilds=3
      nsilol=41
      nfcoil=18
      nrogow=1
      nacoil=1
      mfcoil=18
      necoil=122
      nvesel=24
      mpress=201
      nesum=6
      magpri67=29
      magpri322=31
      magprirdp=8
      magudom=5
      maglds=3
      mse315=11
      mse45=15
      mse15=10
      mse1h=4
      mse315_2=5
      mse210=24
      libim=32
      nmsels=16
      nnece=40
      nnecein=80
      neceo=1
      nnnte=801
      ngam_vars=9
      ngam_u=5
      ngam_w=3
      nlimit=160
      nlimbd=6
      nangle=64
      ntangle=12
      nfbcoil=12
      mccoil=6
      micoil=12

      ndata=61
      nwwcur=32
      nffcur=32
      nppcur=32
      nercur=32

      ntime=1001
      ndim=3200
      kxiter=515
      mqwant=30
      mbdry=300
      mbdry1=110
      nxtram=10
      nxtlim=9
      nco2v=3
      nco2r=2
      modef=4
      modep=4
      modew=4
      nfourier=5

      ! set a loop size for cyclic reduction that is compatible with any
      ! grid size (this is not the exact number of loops, but a slight
      ! overestimate - some exact numbers are: 33[123], 65[307],
      ! 129[741], 257[1731])
      i=1
      nhpwr=1
      do j=1,11
        i=i*2
        if (i.eq.(nh-1)) then
          nhpwr=j
          exit
        endif
      enddo
      icycred_loopmax=nh*(nhpwr-1)

      mfila=10

      ! these parameters are part of eparm, but get used when reading
      ! tables, etc. so they need to be initialized here
      icutfp=0
      iecurr=1
      ifitvs=0
      ivesel=0

      ! checks for solution validity (see chkerr.f90) 
      li_max=2.5
      li_min=0.05
      betap_lim=6.0
      plasma_diff=0.08
      aminor_max=75.0
      aminor_min=30.
      elong_max=4.0
      elong_min=0.8
      rout_max=240.
      rout_min=90.0
      zout_max=100.
      zout_min=-100.
      rcurrt_max=240.
      rcurrt_min=90.0
      zcurrt_max=100.
      zcurrt_min=-100.
      qstar_max=200.
      qstar_min=1.
      betat_lim=25.
      gapin_lim=-0.2
      gapout_lim=-0.2
      gaptop_lim=-0.2
      sepin_check=-90.0
      qout_max=200.
      qout_min=1.
      dbpli_lim=0.05
      delbp_lim=0.08
      end subroutine set_eparm_defaults

!**********************************************************************
!>
!!    this subroutine calculates experiment variables that depend on
!!    other experiment variables
!!
!**********************************************************************
      subroutine set_eparm_dependents()
      use eparm
      implicit none

      nmtark=mse315+mse45+mse15+mse1h+mse315_2+mse210
      nstark=nmtark+libim

      magpol=magpri67+magpri322+magprirdp+magudom
      magpri=magpol+maglds

      nsilop=nsilds+nsilol

      nbwork=nsilop
      msbdry=mbdry+nsilop+nfcoil+1
      msbdr2=2*msbdry

      npcurn=nffcur+nppcur
      mfnpcr=nfcoil+npcurn+nvesel+nwwcur+nesum+nfcoil+nercur 
      nrsmat=nsilop+magpri+nrogow+nffcur+1+npcurn+nwwcur+mpress+nfcoil+nstark+nnece+neceo

      nwcurn=nwwcur+npcurn

      ncurrt=nvesel+nesum+nfcoil

      end subroutine set_eparm_dependents
