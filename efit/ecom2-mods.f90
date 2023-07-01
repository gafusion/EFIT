      module var_cecoil
      integer*4 iecurr
      real*8,dimension(:,:),allocatable :: rsilec
      real*8,dimension(:,:),allocatable :: rmp2ec
      real*8,dimension(:,:),allocatable :: gridec
      real*8,dimension(:,:),allocatable :: rfcec
      real*8,dimension(:),allocatable :: ecurrt,pecur
      real*8,dimension(:,:),allocatable :: rbdrec, rsolec
      real*8,dimension(:,:),allocatable :: recec
      end module var_cecoil

      module var_cvesel
      real*8,dimension(:),allocatable :: vcurrt
      real*8,dimension(:,:),allocatable :: rsilvs
      real*8,dimension(:,:),allocatable :: rmp2vs
      real*8,dimension(:,:),allocatable :: gridvs
      real*8,dimension(:,:),allocatable :: rfcvs
      real*8,dimension(:,:),allocatable :: rbdrvs, rsolvs
      real*8,dimension(:,:),allocatable :: rvsfc
      real*8,dimension(:,:),allocatable :: rvsec
      real*8,dimension(:,:),allocatable :: rvsvs
      end module var_cvesel

      module var_stable
      real*8,dimension(:,:),allocatable :: gbrpc,gbzpc
      end module var_stable

      module var_rtable
      real*8,dimension(:,:),allocatable :: rsilpc
      real*8,dimension(:),allocatable :: fgowpc
      real*8,dimension(:,:),allocatable :: rmp2pc
      real*8,dimension(:,:),allocatable :: rprepc
      end module var_rtable

      module var_atable
      real*8,dimension(:,:),allocatable :: rsilac
      real*8,dimension(:,:),allocatable :: rmp2ac
      real*8,dimension(:,:),allocatable :: gridac
      end module var_atable

      module var_fourier
      real*8,dimension(:),allocatable :: thetav
      real*8,dimension(:,:),allocatable :: sinta,costa
      real*8,dimension(:,:),allocatable :: vecta
      end module var_fourier

      module var_bscom
      integer*4 :: kfffnc,kffhord,kffknt,kppfnc,kpphord,kppknt, &
                   kwwfnc,kwwhord,kwwknt
      real*8 :: pptens,fftens,wwtens
      integer*4,dimension(:),allocatable :: kffbdry,kff2bdry, &
                                            kppbdry,kpp2bdry, &
                                            kwwbdry,kww2bdry
      real*8,dimension(:),allocatable :: ffbdry,ff2bdry,ffknt, &
                                         ppbdry,pp2bdry,ppknt, &
                                         wwbdry,ww2bdry,wwknt
      end module var_bscom

      module var_autokknot
      integer*4 :: kautoknt,kakloop,kakiter
      real*8 :: akchiwt,akerrwt,aktol,akgamwt,akprewt
      end module var_autokknot

      module var_autok
      integer*4 :: ks_a,lconvr_a,ktime_a,kerror_a,kadknt, &
                   kappknt,kaffknt,kawwknt, kaeeknt,mxiter_a
      real*8,dimension(:),allocatable :: aeeknt,awwknt,affknt,appknt
      end module var_autok

      module var_fixstark
      integer*4 :: ifirst
      real*8,dimension(:,:),allocatable :: save_gam,save_tangam
      end module var_fixstark

      module var_cwork2 
      real*8,dimension(:),allocatable :: wrsp_cw2,work_cw2,bdata_cw2
      real*8,dimension(:,:),allocatable ::ematrix_cw2,einv_cw2,arsp_cw2

      end module var_cwork2

      module var_cwork3
      integer*4 :: lkx, lky
      real*8,dimension(:),allocatable ::scrat_cw3,bscra_cw3,work_cw3, &
                                        cscra_cw3,dscra_cw3
      end module var_cwork3

      module var_gwork1
      integer*4 mfila
      integer*4,dimension(:),allocatable :: irfila,jzfila
      real*8,dimension(:),allocatable :: rmx,zmx,wsilpc,wmp2pc, &
                                         wfcpc,wecpc,wvspc
      real*8,dimension(:,:),allocatable :: rsilpf, rmp2pf
      real*8 wpcpc
      end module var_gwork1

      module var_cwork4
      real*8,dimension(:),allocatable :: scraps
      integer*4,dimension(:),allocatable ::npxtra
      end module var_cwork4

      module var_jwork4
      real*8,dimension(:),allocatable :: workb_jw4
      end module var_jwork4

      module var_cww
      integer*4 lwx,lwy
      end module var_cww
