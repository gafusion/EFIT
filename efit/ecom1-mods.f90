      module var_acoilrz
      real*8,dimension(:),allocatable  :: racoil,zacoil,wacoil, &
                     hacoil
      end module var_acoilrz

      module var_ccurrn
      use eparm,only:nwnh
      real*8,dimension(:),allocatable :: pcurrt,pcurrw,pcurrtpp
      integer*4 icurrt,icinit,islpfc,icprof
      data islpfc/0/

      real*8 fconst,betap0
      end module var_ccurrn

      module var_czero
      real*8,dimension(:),allocatable :: zero,www
      integer*4 iweigh
      end module var_czero

      module var_parame
      real*8,dimension(:),allocatable :: volp,pprime,pres,ffprim,fpol &
         ,qpsi,r2surf,rpres,curmid,cjor,bpolss &
         ,rmajz0,bprmaj,btrmaj,r1surf,r2surg &
         ,bbfpol,ccfpol,ddfpol
      real*8 r1bdry,r2bdry,carea,jwantm,r2wdry,r4wdry,rpwdry,rp2wdry 
      real*8,dimension(:),allocatable :: r1sdry,r2sdry,vzeroj,sizeroj 
      data jwantm/3/
      end module var_parame

      module var_cqfit
      use set_kinds

      real*8 fwtqa,qvfit,fbrdy,fcentr,cjmaxi &
           ,goli,gocur,goqq,gobeta,cbrsp &
           ,rqajtor,rqaftor,rqapetor,rqafetor
      data cbrsp/3.0e-5_dp/
      integer*4 nqiter,ivacum 
      data nqiter/10/
      real*8,dimension(:),allocatable :: rqajx,rjjjx
      real*8,dimension(:),allocatable :: rqafx,rjjfx
      real*8,dimension(:),allocatable :: rqawx,rjjwx
      end module var_cqfit

      module var_rmatri
      real*8,dimension(:),allocatable :: brsp,brspss
      real*8,dimension(:),allocatable :: tsaisq
      real*8 cond,tsaifc
      integer*4 nparam,nfnpcr,nfnwcr,nbase
      end module var_rmatri
  
      module var_chimag
      real*8,dimension(:),allocatable :: chi2rm
      end module var_chimag

      module var_exdat2
      use set_kinds
      real*8,dimension(:),allocatable :: bcentr
      real*8 rcentr
      data rcentr/1.6955_dp/
      end module var_exdat2

      module var_exdata
      use set_kinds
      integer*4 ishot,itime,ifitvs,iacoil,itimeu,kersil
      data errsil/0.03_dp/
      data kersil/2/

      real*8 serror,fwtcur,elomin,fwtbp,fwtdlc,errsil 
      real*8,dimension(:),allocatable :: fwtsi,rwtsi
      real*8,dimension(:),allocatable :: fwtmp2,rwtmp2
      real*8,dimension(:),allocatable :: fwacoil
      real*8,dimension(:),allocatable :: fwtfc
      real*8,dimension(:),allocatable :: cbrspv
      data elomin/0.90_dp/
      data iacoil/0/
      end module var_exdata

      module var_texdat
      integer*4 ierpla,ico2,ierrdi 
      real*8 bitip
      real*8,dimension(:,:),allocatable ::  silopt
      real*8,dimension(:,:),allocatable ::  expmpi
      real*8,dimension(:,:),allocatable ::  accurt
      real*8,dimension(:,:),allocatable :: fccurt
      real*8,dimension(:,:),allocatable ::eccurt
      real*8,dimension(:,:),allocatable :: denvt
      real*8,dimension(:,:),allocatable :: denrt
      integer*4,dimension(:,:),allocatable :: iccurt
      real*8,dimension(:),allocatable :: pasmat,time,pbinj 
      integer*4,dimension(:),allocatable :: ierpsi
      integer*4,dimension(:),allocatable :: ierpr
      integer*4,dimension(:),allocatable :: iermpi
      real*8,dimension(:),allocatable :: psibit
      real*8,dimension(:),allocatable :: prbit
      real*8,dimension(:),allocatable :: bitmpi
      real*8 :: vbit
           
      real*8,dimension(:),allocatable :: vloopt,psiref,diamag,sigdia,psirefs 
      real*8,dimension(:),allocatable :: bitfc
      integer*4,dimension(:),allocatable :: ierfc
      real*8,dimension(:),allocatable :: bitec
      integer*4,dimension(:),allocatable :: ierec
      real*8,dimension(:),allocatable ::bitic
      integer*4,dimension(:),allocatable ::ieric
      integer*4 ierdia(3)
!sri-mpi
      integer*4 nopbinj
      data vbit/10./
      end module var_texdat

      module var_savfit
      real*8,dimension(:),allocatable ::  swtsi
      real*8,dimension(:),allocatable ::  swtmp2
      real*8,dimension(:),allocatable ::  swtfc
      real*8,dimension(:),allocatable ::  swtec
      real*8 swtcur,swtdlc
      end module var_savfit

      module var_cstark
      use set_kinds
      real*8,dimension(:,:),allocatable :: tangam, siggam,a1gam &
           ,a2gam,a3gam,a4gam,tangam_uncor
      real*8,dimension(:),allocatable  :: fwtgam,chigam,swtgam
      real*8 v30lt,v30rt,v210rt,v210lt
      integer*4,dimension(:),allocatable :: iergam
      integer*4,dimension(:),allocatable :: mseport,mse_spave_on
      integer*4  kstark,iplots,kmtark,klibim,kdomse &
                 ,msebkp,msefitfun,kwaitmse,mse_quiet &
                 ,mse_strict,ok_30rt,ok_210lt,mse_usecer &
                 ,mse_certree,mse_use_cer330,mse_use_cer210
      data iplots/1/
      data  msebkp/0/,msefitfun/1/
      data kwaitmse/0/,dtmsefull/0.0/
      data mse_strict/0/,t_max_beam_off/0.0/,ok_30rt/0/,ok_210lt/0/
      data kdomse/0/
      data mse_usecer/0/,mse_certree/0/
      data mse_use_cer210/0/,mse_use_cer330/0/
      data v30lt/0.0/,v30rt/0.0/
      data v210lt/0.0/,v210rt/0.0/

      real*8 chigamt,chilibt,dtmsefull,t_max_beam_off 
      real*8,dimension(:,:),allocatable :: rbrpc,rbzpc,rgampc
      real*8,dimension(:,:),allocatable :: rbrfc,rbzfc
      real*8,dimension(:,:),allocatable :: rbrec,rbzec,rgamec
      real*8,dimension(:,:),allocatable :: rbrvs,rbzvs,rgamvs
      real*8,dimension(:,:),allocatable :: rgamfc
      real*8,dimension(:,:),allocatable :: rhsgam,rrgam,zzgam &
           ,starkar,starkaz,a5gam,a6gam,a7gam,a8gam
      real*8,dimension(:,:),allocatable :: cmgam,spatial_fix
      real*8,dimension(:,:),allocatable :: rgamac,rbrac,rbzac
      real*8,dimension(:),allocatable :: btgam,sistark,qstark &
           ,rmse_gain,rmse_slope,rmse_scale,rmse_offset,sigam &
           ,bzmse,bzmsec,cjmse,cjmsec,rhogam
      real*8,dimension(:,:,:,:),allocatable :: spatial_avg_gam
      end module var_cstark

      module var_msels
      use set_kinds
      real*8,dimension(:,:),allocatable :: bmselt,sbmselt,fwtbmselt, &
           rrmselt,zzmselt,l1mselt,l2mselt,l4mselt,emselt,semselt, &
           fwtemselt,swtbmselt,swtemselt,cmmls,cmels,rhsmls,rhsels, &
           l3mselt,cmmls2,cmmlsv
      real*8,dimension(:),allocatable :: swtbmsels,swtemsels,chimls,chiels
      real*8 :: avemsels, fmlscut, tchimls, tchiels
      real*8,dimension(:,:),allocatable :: rmlspc
      real*8,dimension(:,:),allocatable :: rmlsec,relsec
      real*8,dimension(:,:),allocatable :: rmlsvs,relsvs
      real*8,dimension(:),allocatable ::  simls,sinmls,btmls,brmls,bzmls
      real*8,dimension(:),allocatable ::  chi2mls, chi2gamt

      integer*4,dimension(:,:),allocatable :: iermselt
      integer*4,dimension(:),allocatable :: kerrot
      integer*4 kdomsels, mmbmsels
      character*3 :: synmsels
      data avemsels/10./, synmsels/'SYN'/, kdomsels/0/, fmlscut/1.e-6_dp/
      end module var_msels

      module var_rpedge
      real*8,dimension(:),allocatable :: rsilpe
      real*8,dimension(:),allocatable :: rmp2pe
      real*8 fgowpe
      real*8,dimension(:),allocatable :: rprepe
      real*8,dimension(:),allocatable ::  gbdrpe
      real*8,dimension(:),allocatable :: rgampe,rbrpe,rbzpe
      real*8,dimension(:),allocatable ::  rmlspe,relspe
      end module var_rpedge

      module var_rfedge
      real*8,dimension(:),allocatable :: rsilfe
      real*8,dimension(:),allocatable :: rmp2fe
      real*8 fgowfe
      real*8,dimension(:),allocatable :: gbdrfe
      real*8,dimension(:),allocatable :: rgamfe,rbrfe,rbzfe
      real*8,dimension(:),allocatable :: rmlsfe,relsfe
      end module var_rfedge

      module var_cece
      use set_kinds
      integer*4 ierecebz,necein,jo,kfixro,nece,kece,kecebz,mecein &
                ,kfitece,kinputece,kcallece,nharm,kdoece,kgeteceb & 
                ,mtxece,ksetece,kfixrece,nfit,kcmin,nconstr

      data ksetece/0/,kgeteceb/0/,kcallece/2/
      data kdoece/0/,nconstr/1/,mtxece/0/

      character*4 eceiter
      data eceiter/'pair'/

      real*8 receo,fwtecebz,fwtecebz0,rteo,zteo,robit,swtecebz,chiecebz &
             ,ecebzfit,ecebzbit,fwtnow,zeceo,chisqfit,xfit(10),tchiece &
             ,recebzdz,gecebzdz,eceerror

      data eceerror/0.03_dp/

      real*8,dimension(:),allocatable :: teecein,feece,errorece &
           ,teecein0,feece0,errorece0 &
           ,becein,recein,teeceinr
      real*8,dimension(:),allocatable :: recem,recep,fwtece,fwtece0 &
                                 ,swtece,chiece,ecefit,ecebit
      real*8,dimension(:),allocatable :: recebzfc
      real*8,dimension(:),allocatable :: gecebzpc
      real*8,dimension(:),allocatable :: recebzec
      real*8,dimension(:,:),allocatable  :: recefc
      real*8,dimension(:,:),allocatable :: gecepc,geceppc,gecempc
      real*8,dimension(:,:),allocatable  :: receec
      integer*4,dimension(:),allocatable :: ierece
      real*8,dimension(:,:),allocatable :: recepc
      real*8,dimension(:,:),allocatable :: brspece
      real*8,dimension(:),allocatable :: recebzpc
      real*8,dimension(:),allocatable :: brspecebz,cmecebz
      real*8,dimension(:,:),allocatable :: cmece
      real*8,dimension(:,:),allocatable :: recevs
      real*8,dimension(:),allocatable :: recebzvs
      real*8,dimension(:),allocatable :: recedz,gecedz,rtep,rtem,rpbit,rmbit
      real*8,dimension(:,:),allocatable :: receac
      real*8,dimension(:),allocatable  :: recebzac
      real*8,dimension(:),allocatable  :: teecer,rrr,bbf,teeceb
      real*8,dimension(:),allocatable :: receoi
      real*8,dimension(:,:),allocatable :: recemi,recepi
      end module var_cece

      module var_calchi
      real*8 chipre,chitot,saiip
      real*8,dimension(:),allocatable :: saisil
      real*8,dimension(:),allocatable :: saimpi
      real*8,dimension(:),allocatable :: saipr
      real*8,dimension(:),allocatable :: saipre
      real*8,dimension(:),allocatable :: saifc
      end module var_calchi

      module var_dlc
      real*8 sigdlc
      real*8,dimension(:),allocatable :: dfluxc,cdflux,edflux
      real*8,dimension(:),allocatable :: rspdlc
      end module var_dlc
    
      module var_comco2
      use set_kinds
      real*8,dimension(:,:),allocatable :: rco2r
      real*8,dimension(:,:),allocatable :: rco2v
      real*8,dimension(:),allocatable :: chordv
      real*8,dimension(:),allocatable :: chordr
      real*8 zcentr
      real*8,dimension(:,:),allocatable :: dco2r
      real*8,dimension(:,:),allocatable :: dco2v
      data zcentr/0./
      end module var_comco2

      module var_check
      use eparm,only:ntime
      integer*4,dimension(:,:),allocatable ::erflag
      integer*4 kflag(30),lflag,ktimeo
      end module var_check

      module var_consum
      use set_kinds
      real*8 condno,condin 
      data condin/1.0e-06_dp/
      real*8,dimension(:),allocatable :: cerror,csibry,csimag, &
                  cvolp,crmaxi,czmaxi,cemaxi,cqmaxi,cchisq &
                  ,brfbc,tvfbrt,cdelz
      end module var_consum

      module var_cxray
      real*8,dimension(:),allocatable :: rxray,zxray,xangle
      integer*4 ksxr0(10),ksxr2(10),idosxr

      end module var_cxray

      module var_mercie
      integer*4 imerci
      real*8,dimension(:),allocatable :: rzzmax,zzmax
      end module var_mercie

      ! OPT_INPUT >>>
      module opt_input
         integer*4 mode_in
         character cmdfile_in*15, shotfile_in*15
         integer*4 shot_in
         real*8 starttime_in
         real*8 deltatime_in
         integer*4 steps_in
         logical use_opt_input
         character snapext_in*82
         character(80),dimension(:),allocatable :: inpfile_in
      end module opt_input
      ! OPT_INPUT <<<

      ! MPI >>>
      module mpi_info
         integer :: rank, nproc, ierr
         integer*4,dimension(:),allocatable :: dist_data
         integer*4,dimension(:),allocatable :: dist_data_displs
         double precision,dimension(:,:),allocatable :: fwtgam_mpi
      end module mpi_info
      ! MPI <<<

!jal 2/23/04 add iplcout=1 print plasma and pf currents to gfile
      module var_input1
      use set_kinds
      logical write_Kfile ,fitfcsum

      integer*4 icondn,itek,kdata,itrace,ierchk,iconvr,ixray,itell &
           ,kprfit,licalc,ibound,ibatch,idite,ilaser, lookfw &
            ,kdot, icutfp, keqdsk,kdofit,kbetapr,kplotpr,kpressb &
            ,kdoqn,kfcurb,kpcurb, kzeroj, ncstne,ncstte,ncstfp,ncstpp &
            ,kgrid,kwripre,negcur, kframe,kskipvs,icntour, iavdpl &
            ,jwake,limvs, kbound,kgraph,istore,iout,kdopre &
            ,iishot,kktime,iplcout,ksigma, kwritime 
      integer*4 iteks, mxiters, zelipss, n1coils
      integer*4 itekt, mxitert, zeliptt, n1coilt

      data kgrid/1/,kwripre/0/,kwritime/0/
      data licalc/1/
      data ksigma/0/
      data kdoqn/0/
      data icntour/0/,cstabne/0.0/,cstabte/0.0/ &
           ,limvs/1/,ncstpp/1/,ncstfp/1/,kzeroj/0/
      data kbound/0/
      data kdofit/0/
      data kdopre/0/

      real*8 cutip, dtdot, xpsimin,fcurbd,pcurbd,prbdry,sgprmin &
           ,prespb,tipbry,tepbry,dipbry,depbry,pbimpb,sigppb,sigpreb &
           ,sigtipb,sigtepb,sigdipb,sigdepb,fwtpd,fwtfd,cstabte,cstabne &
           ,dp1dxf, sgtimin, alphafp,xpsialp, vsdamp &
           ,rminvs,rmaxvs,zminvs,zmaxvs,relbps,zbound,yvs2,saimin &
           ,fztor,fzpol,tcurrp,fpolvs, rbound, dnmin, sigprebi,pressbi &
           ,alphamu,saicon,rsepex, ttimeb,ddtime
      data alphafp/0./,sigppb/1000./
      data kframe/0/,rminvs/0/,rmaxvs/100./,zminvs/-100./,zmaxvs/100./
      data kskipvs/0/,vsdamp/0/,relbps/0.004_dp/,zbound/0.0/,rbound/0.0/
      data dnmin/1.0/
      data saimin/60.0/,saicon/60.0/

      real*8,dimension(:),allocatable :: rzeroj
      real*8,dimension(:),allocatable :: fwtpre
      real*8,dimension(:),allocatable :: vforcep,vforcet
      real*8,dimension(:),allocatable :: fwtfcsum
      character*82  snap_file 
      end module var_input1

      module var_inputc

      character*12 mfitpop
      character*5 mfvers(2)
      data mfvers(1)/'11/23'/,mfvers(2)/'/2020'/
      character(4),dimension(:),allocatable :: limloc
      character(10),dimension(:),allocatable :: vsname
      character(10),dimension(:),allocatable :: mpnam2
      character(10),dimension(:),allocatable :: lpname
      character  filimt*100,cshot*6,jdebug*4
      integer idebug,efitversion
      data idebug/0/,efitversion/20201123/
      data jdebug/'NONE'/
      end module var_inputc

      module var_input4
      character(80),dimension(:),allocatable :: ifname
      end module var_input4

      module var_switch
      integer*4 nqaxis,isumip,jbeta,jli,nqwant
      data nqwant/0/,jbeta/1/,jli/2/
      real*8 sumip,fbetap,fli,fbetat,fbetan 
      real*8,dimension(:),allocatable :: qsiw, pasmsw,fqsiw,siwantq
      real*8,dimension(:,:),allocatable :: fgowsw
      end module var_switch

      module var_siloop
      real*8,dimension(:),allocatable :: rsi,zsi,wsi,hsi
      end module var_siloop

      module var_slname
      real*8,dimension(:),allocatable :: as,as2
      integer*4 nslref
      data nslref/1/
      end module var_slname

      module var_ecoil 
      use eparm,only:necoil
      real*8,dimension(:),allocatable :: re,ze,we,he,ecid
      end module var_ecoil 

      module var_fcoil
      real*8,dimension(:),allocatable :: rf,zf,wf,hf, &
           af,af2,rsisfc,turnfc,fcid,fcturn
 
      end module var_fcoil 

      module var_mprobe
      real*8,dimension(:),allocatable :: xmp2,ymp2,amp2,smp2,patmp2
      integer*4 nsmp2
      end module var_mprobe

      module var_limite
      integer*4 limitr,iplim,limfag,limitr_180
      data limfag/2/
      real*8,dimension(:),allocatable :: xlim,ylim, xlim_180, ylim_180
      real*8,dimension(2) :: rwstrip1, zwstrip1,rwstrip2, zwstrip2 
      real*8 yltype, yltype_180
      logical dowstrip

      data dowstrip/.false./
      end module var_limite

      module var_mimite
      use eparm,only:nlimit
      integer*4 mimitr,mimitr_180
      real*8,dimension(:),allocatable :: xmim,ymim,xmim_180,ymim_180
      end module var_mimite

      module var_udata
      integer*4,dimension(:),allocatable :: ipsi,irogw,imag2,iplasm &
           ,idlopc,ifc,iec
      end module var_udata

      module var_morsum
      real*8,dimension(:),allocatable :: csumip,tratio,aveerr
      integer*4,dimension(:),allocatable :: iermax,jermax
      end module var_morsum

      module var_bdsend
      integer*4 nbbbs,nbskip,nbdrymx, nbdryp
      data nbskip/2/,nbdrymx/110/
      real*8,dimension(:),allocatable :: rbbbs,zbbbs
      end module var_bdsend

      module var_fxbry
      use set_kinds
      integer*4 nbdry,nbdryss,nsol
      data nsol/0/
      logical fitts 
      real*8 wsisol,dselsum/1.e-12_dp/
      real*8,dimension(:),allocatable :: rbdry,zbdry,fwtbdry,fwtbry,sigrbd &
                                 ,sigzbd,rbdry0, zbdry0,rbdryss,zbdryss
      real*8,dimension(:),allocatable ::rsol, zsol, fwtsol, fwtsolw
      real*8,dimension(:,:),allocatable :: rbdrfc, rsolfc
      real*8,dimension(:,:),allocatable :: rbdrac
      real*8,dimension(:,:),allocatable :: rbdrpc,rsolpc
      real*8,dimension(:,:),allocatable :: gbdrpc
      end module var_fxbry

      module var_fwtdz
      use set_kinds
      logical fitdelz
      integer*4 ndelzon,ifitdelz 
      data fitdelz/.false./,ndelzon/999/,relaxdz/1.0/, &
           stabdz/-1.e-4_dp/,scaledz/1.e-03_dp/,ifitdelz/1/
      data errdelz/0.06_dp/
      real*8 errdelz,fgowdz,scaledz,stabdz,relaxdz,cdeljsum
      real*8,dimension(:),allocatable :: gmp2dz
      real*8,dimension(:),allocatable :: gsildz
      real*8,dimension(:),allocatable :: gbrdz,gbzdz,rgamdz
      real*8,dimension(:),allocatable :: rmlsdz,relsdz
      real*8,dimension(:),allocatable :: rpredz,rprwdz
      real*8,dimension(:),allocatable :: gbdrdz
      real*8,dimension(:),allocatable :: rdjdz
      end module var_fwtdz
      
      module var_combry
      use set_kinds
      real*8,dimension(:),allocatable ::erbloc, erbsloc
      real*8 erbmax,erbave,erbcom,erbsmax,erbsave
      data erbcom/1.0e-2_dp/
      end module var_combry

      module var_fbysta
      real*8 cfcoil, psibry0
      integer*4 ifref
      real*8,dimension(:),allocatable :: fcsum,fczero
      end module var_fbysta

      module var_prdata
      integer*4 npress,npteth,nption,npneth,nbeam,ndokin,nbrmcrd &
               ,nptef,npnef,nptionf,nmass 
      data ndokin/1/
      real*8 pbeamb,pressb,zeffvs,zlowimp
      integer*4,dimension(:),allocatable :: ivbcuse

      real*8,dimension(:),allocatable :: pressr,rpress,zpress,tethom,rteth &
                   ,zteth,tionex,rion,zion,dnethom,rneth,zneth &
                   ,pbeam,sibeam,sgteth,sigti,sigpre,sgneth,precal &
                   ,dnbeam,dmass,scalepr,premea,saipre2
      real*8,dimension(:),allocatable :: pbimth,bremin,bremsig,brmrtan, &
                                 brmzelev,dnbthom
      end module var_prdata

      module var_cerfit
      use set_kinds
      integer*4 keecur,keefnc,keeknt,needer,keehord
      integer*4,dimension(:),allocatable :: keebdry,kee2bdry
      real*8 ecurbd,eetens
      real*8,dimension(:,:),allocatable :: rgamer
      real*8,dimension(:,:),allocatable :: rmlser,relser
      real*8,dimension(:),allocatable ::eeknt,eebdry,ee2bdry,cerer
      real*8,dimension(:),allocatable :: ermse
      real*8,dimension(:,:),allocatable :: e1rbz,e2rbz,e3rbr
      real*8,dimension(:),allocatable :: ermid, eshear,epoten,rhovn, &
                    rpmid,xmid,sigrid,sipmid, &
                    brhovn,crhovn,drhovn,rhopmid
      data keecur/0/,ecurbd/0.0/,keefnc/0/,eetens/5.0_dp/
      end module var_cerfit

      module var_ccgama
      use set_kinds
      integer*4 kcalpa,kcgama,kcomega
      real*8 fwtxx,fwtxxj,fwtxxq,fwtxxb,fwtxli 
      data fwtxx/0.2_dp/,fwtxxj/1./
      data  fwtxxq/1./,fwtxxb/1./,fwtxli/1./
      data kcalpa/0/,kcgama/0/,kcomega/0/

      real*8,dimension(:,:),allocatable ::  calpa
      real*8,dimension(:),allocatable :: xalpa, alpax
      real*8,dimension(:,:),allocatable ::  cgama
      real*8,dimension(:),allocatable :: xgama, gamax
      real*8,dimension(:,:),allocatable ::  comega
      real*8,dimension(:),allocatable :: xomega
      end module var_ccgama

      module var_cccoils
      integer*4 kccoils,kcloops
      data kccoils/0/
      real*8,dimension(:,:),allocatable :: ccoils
      real*8,dimension(:),allocatable :: xcoils
      real*8,dimension(:,:),allocatable :: cloops
      real*8,dimension(:),allocatable :: xloops
      end module var_cccoils

      module var_tionfit
      real*8 chisqti, tibdry,stibdry
      real*8,dimension(:),allocatable :: tifit
      real*8,dimension(:),allocatable :: stitho,xsiion,tithom
      end module var_tionfit

      module var_telnfit
      real*8 chisqte,tebdry,stebdry
      real*8,dimension(:),allocatable :: tefit
      end module var_telnfit

      module var_dionfit
      real*8 dibdry,sdibdry
      real*8,dimension(:),allocatable :: snitho,dnitho
      end module var_dionfit

      module var_delnfit
      real*8,dimension(:),allocatable :: defit
      real*8 chisqne,debdry,sdebdry,fco2ne
      
      end module var_delnfit

      module var_tsrz
      use set_kinds
      real*8,dimension(:),allocatable :: zuperts,rlibim
      real*8 zlowerts,rmajts,zlibim
      data rmajts/1.94_dp/,zlibim/-0.127_dp/
      end module var_tsrz

      module var_climxx
      real*8,dimension(:),allocatable :: xlimbd,ylimbd
      end module var_climxx

      module var_fdbkgr
      real*8 gsum
      real*8,dimension(:,:),allocatable :: grdfdb
      end module var_fdbkgr

      module var_fdbkcl
      real*8 brfb(2)
      integer*4 kct1,kct2,kct3,kct4
      real*8,dimension(:),allocatable :: fb_plasma
      end module var_fdbkcl

      module var_coiln1
      use eparm,only:ntime,magpri
      integer*4 n1coil,iern1
      data n1coil/0/
      real*8,dimension(:),allocatable :: signn1
      real*8,dimension(:),allocatable :: curtn1
      end module var_coiln1

!***** !EJS(2014)
!  Need to update values -- wherever they are specified:
! mccoil = 6   (probably is 3 now)
! micoil = 12  (probably is 6 now)
!*****
      module var_coilcc
      use eparm,only:ntime,mccoil,micoil
      integer*4 nccoil,iercc,nicoil
      data nccoil/1/,nicoil/1/
      logical oldccomp,oldcomp 
      data oldccomp/.false./,oldcomp/.false./
!***** !EJS(2014)
! Not sure why these are separate arrays.
! They seem to be the same as the columns of curccoi and curicoi
! Are they only used in magsigma?  If so, the following changes are
! probably not needed until some future time when magsigma is updated.
!*****
      real*8,dimension(:),allocatable :: curc139,curc79,curc199,curiu30 &
                    ,curiu90,curiu150,curil30,curil90,curil150 &
                    ,curc259,curc319,curc19        & !EJS(2014)
                    ,curiu210,curiu270,curiu330    & !EJS(2014)
                    ,curil210,curil270,curil330      !EJS(2014)
      real*8,dimension(:,:),allocatable :: curccoi
      real*8,dimension(:,:),allocatable :: curicoi
      end module var_coilcc

      module var_btcomp
      use eparm,only:ntime,magpri
      integer*4 ibtcomp,ierbtc
      data ibtcomp/1/
      real*8,dimension(:),allocatable :: signbt
      real*8,dimension(:),allocatable :: bti322
      end module var_btcomp

      module var_subic
      use eparm,only:modef,modep,modew
      integer*4 nodef,nodep,nodew
      real*8,dimension(:),allocatable :: xnodef
      real*8,dimension(:),allocatable :: xnodep
      real*8,dimension(:),allocatable :: xnodew
      integer*4,dimension(:),allocatable :: kbasef
      integer*4,dimension(:),allocatable :: kbasep
      integer*4,dimension(:),allocatable :: kbasew
      end module var_subic

      module var_vtor
      use set_kinds
      integer*4 nomegat,kvtor,kwwcur,kwcurn,kplotp,npresw &
                ,nsplot,kdovt
      real*8 betapw0,enw,emw,rvtor,wcurbd,gammaw,rbetaw, & 
             preswb,chiprw
      data kvtor/0/,rvtor/1.70_dp/,preswb/0.0/,wcurbd/0.0/,betapw0/0.0/ &
           ,kwwcur/2/,kplotp/1/,kdovt/0/,nsplot/4/ 
      real*8,dimension(:),allocatable :: omegat,rpresw,zpresw,presw, &
                       sigprw,rpresws,fwtprw,romegat,zomegat &
                       ,sigome,rpreswv,scalepw
      real*8,dimension(:,:),allocatable :: rprwpc
      real*8,dimension(:),allocatable :: betapw,betatw,wplasw
      real*8,dimension(:),allocatable :: rgrvt,pwprim,pressw,prwcal, &
                              saiprw,rgsvt
      real*8,dimension(:),allocatable :: presst
      end module var_vtor

      module var_fitec
      real*8,dimension(:),allocatable :: fwtec,cecurr,saiec
      logical writepc
      data writepc/.false./
      end module var_fitec

      module var_pflocal
      real*8 psiecn,dpsiecn,rkec,cjeccd
      data psiecn/0.0/,dpsiecn/0.0/
      real*8,dimension(:),allocatable :: cjorec,ffprec
      end module var_pflocal

      module var_ctanhts
      character*2 fitzts
      data fitzts/'no'/
      real*8,dimension(:),allocatable :: ztssym,ztswid,ptssym
      logical,dimension(:),allocatable :: ztserr
      end module var_ctanhts

      module var_qsurfac
      real*8,dimension(:),allocatable :: psin32,psin21,rq32in, &
                     rq21top
      end module var_qsurfac
!----------------------------------------------------------------------
!-- New magnetic uncertainties commons                               --
!----------------------------------------------------------------------
      module var_initerror
      real*8 sigmaip0
      real*8,dimension(:),allocatable :: sigmaf0
      real*8,dimension(:),allocatable :: sigmae0
      real*8,dimension(:),allocatable :: sigmafl0
      real*8,dimension(:),allocatable :: sigmamp0
      end module var_initerror

      module var_magerror
      use set_kinds
      integer*4 imagsigma,icountmagsigma
      data imagsigma/0/, errmag/1.0e-3_dp/, errmagb/1.e-2_dp/
      real*8 errmag,errmagb
      real*8,dimension(:,:),allocatable :: sigmaf
      real*8,dimension(:),allocatable :: sigmab,sigmaip
      real*8,dimension(:,:),allocatable :: sigmae
      real*8,dimension(:,:),allocatable :: sigmafl,gradsfl
      real*8,dimension(:,:),allocatable :: sigmamp,gradsmp,bpermp
      end module var_magerror

      module var_psilopdat
      real*8,dimension(:),allocatable :: psircg,psi_k,psi_rc,vrespsi,t0psi
      real*8,dimension(:,:),allocatable :: devpsi,rnavpsi
      integer*4,dimension(:,:),allocatable :: navpsi
      end module var_psilopdat

      module var_plasmacurrdat
      real*8 prcg,p_k,p_rc,vresp,t0p
      real*8,dimension(:),allocatable :: devp,navp,rnavp
      end module var_plasmacurrdat

      module var_ccoilsdat
      real*8,dimension(:),allocatable :: ccrcg,cc_k,cc_rc,vrescc,t0cc
      real*8,dimension(:,:),allocatable :: devcc
      integer*4,dimension(:,:),allocatable :: navcc
      end module var_ccoilsdat

      module var_icoilsdat
      real*8,dimension(:),allocatable :: xicrcg,xic_k,xic_rc,vresxic,t0xic
      real*8,dimension(:,:),allocatable :: devxic
      integer*4,dimension(:,:),allocatable :: navxic
      end module var_icoilsdat

      module var_n1coildat
      real*8 xn1rcg,xn1_k,xn1_rc,vresxn1,t0xn1
      real*8,dimension(:),allocatable :: devxn1
      integer*4,dimension(:),allocatable :: navxn1
      end module var_n1coildat

      module var_vloopdat
      real*8 vlrcg,vl_k,vl_rc, vresvl,t0vl
      real*8,dimension(:),allocatable :: devvl
      integer*4,dimension(:),allocatable :: navvl
      end module var_vloopdat

      module var_diamdat
      real*8 diamrcg,diam_k,diam_rc,vresdiam,t0diam
      real*8,dimension(:),allocatable :: devdiam
      integer*4,dimension(:),allocatable :: navdiam
      end module var_diamdat

      module var_denvdat
      real*8,dimension(:),allocatable :: denvrcg,denv_k,denv_rc,vresdenv,t0denv
      real*8,dimension(:,:),allocatable :: devdenv
      integer*4,dimension(:,:),allocatable :: navdenv
      end module var_denvdat

      module var_denrdat
      real*8,dimension(:),allocatable :: denrrcg,denr_k,denr_rc,vresdenr,t0denr
      real*8,dimension(:,:),allocatable :: devdenr
      integer*4,dimension(:,:),allocatable :: navdenr
      end module var_denrdat

      module var_magprobdat
      real*8,dimension(:),allocatable :: xmprcg,xmp_k,xmp_rc,vresxmp,t0xmp
      real*8,dimension(:,:),allocatable :: devxmp,rnavxmp
      integer*4,dimension(:,:),allocatable :: navxmp
      end module var_magprobdat

      module var_btcompdat
      real*8 btrcg,bt_k,bt_rc,vresbt,t0bt
      real*8,dimension(:),allocatable :: devbt
      integer*4,dimension(:),allocatable :: navbt
      end module var_btcompdat

      module var_btordat
      real*8 bcrcg,bc_k,bc_rc,vresbc,t0bc
      real*8,dimension(:),allocatable :: devbc, rnavbc
      integer*4,dimension(:),allocatable :: navbc
      end module var_btordat

      module var_fcoildat
      real*8,dimension(:),allocatable :: fcrcg,fc_k,fc_rc,vresfc,t0fc
      real*8,dimension(:,:),allocatable :: devfc,rnavfc
      integer*4,dimension(:,:),allocatable :: navfc
      end module var_fcoildat

      module var_ecoildat
      real*8,dimension(:),allocatable :: ercg,e_k,e_rc,vrese,t0e
      real*8,dimension(:,:),allocatable :: deve,rnavec
      integer*4,dimension(:,:),allocatable :: navec
      end module var_ecoildat

      module var_beamdat
      real*8 beamrcg,beam_k,beam_rc,vresbeam,t0beam
      real*8,dimension(:),allocatable :: devbeam
      integer*4,dimension(:),allocatable :: navbeam
      end module var_beamdat
!
      module commonblocks
      real*8,allocatable :: c(:,:,:,:),wk(:),copy(:,:),bkx(:),bky(:), &
           cw(:,:,:,:),wkw(:),copyw(:,:),bwx(:),bwy(:), &
           sifprw(:),bwprw(:),cwprw(:),dwprw(:),psirz(:,:), &
           sfprw(:),sprwp(:),wgridpc(:),rfcpc(:,:),ct(:,:,:,:), &
           wkt(:),bkrt(:),bkzt(:),psiold(:),psipold(:),psipp(:), &
           work(:),sifpre(:),bwpre(:),cwpre(:),dwpre(:),sfpre(:), &
           sprep(:),worka(:),zeros(:),byringr(:),byringz(:), &
           cjrf(:),cj(:,:,:,:),wkj(:),copyj(:,:),bjx(:),bjy(:), &
           cv(:,:,:,:),wkv(:),copyv(:,:),bvx(:),bvy(:),xsisii(:), &
           f_a(:),f_b(:),f_c(:),f_d(:),pp_a(:),pp_b(:),pp_c(:), &
           pp_d(:),chi_c(:,:,:,:),chi_wk(:),chi_copy(:,:),chi_bkx(:), &
           chi_bky(:),wxin(:),wyin(:),wxout(:),wyout(:),xouts(:), &
           youts(:),bpoo(:),bpooz(:),bpooc(:),bfpol(:),cfpol(:), &
           dfpol(:),rsplt(:),zsplt(:),csplt(:),xxtra(:,:),yxtra(:,:), &
           bpxtra(:,:),flxtra(:,:),fpxtra(:,:),worka2(:)
      end module commonblocks

      subroutine set_ecom_mod1_arrays
      use set_kinds
      use eparm
      use var_parame, only: vzeroj,sizeroj
      use var_exdata, only: fwacoil
      use var_comco2, only: chordv,chordr
      use var_cxray, only: ksxr0,ksxr2,idosxr,xangle,zxray,rxray
      use var_fcoil, only: fcid
      use var_limite, only: rwstrip1, zwstrip1,rwstrip2,zwstrip2
      use var_delnfit, only: fco2ne
      use var_rmatri, only: tsaisq
      use var_consum, only: cdelz
      implicit none   

      tsaisq=0.0
      cdelz=0.0

      vzeroj(1)=0.0
      sizeroj(1)=-1.0
      fwacoil=1*0.
      
      chordv=(/1.486_dp,1.945_dp,2.098_dp/)
      chordr=(/0._dp, 0.1524_dp/)
      ksxr0=10*0
      ksxr2=10*0
      idosxr=1
      xangle= (/    120.,124.,128.,132.,136.,140.,144.,148., &
                     152.,156.,160.,164.,168.,172.,176.,180.,180., &
                      184.,188.,192.,196.,200.,204.,208.,212.,216., &
                      220.,224.,228.,232.,236.,240., &
                   294.5,290.5,286.5,282.5,278.5,274.5,270.5,266.5, &
                   262.5,258.5,254.5,250.5,246.5,242.5,238.5,234.5, &
                   234.5,230.5,226.5,222.5,218.5,214.5,210.5,206.5, &
                   202.5,198.5,194.5,190.5,186.5,182.5,178.5,174.5, &
                   0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)
      zxray= (/-10.7,-10.7,-10.7,-10.7,-10.7,-10.7, &
                    -10.7,-10.7,-10.7,-10.7,-10.7,-10.7,-10.7,-10.7, &
                    -10.7,-10.7,-14.7,-14.7,-14.7,-14.7,-14.7,-14.7, &
                    -14.7,-14.7,-14.7,-14.7,-14.7,-14.7,-14.7,-14.7, &
                    -14.7,-14.7, &
                    130.1,130.1,130.1,130.1,130.1,130.1,130.1,130.1, &
                    130.1,130.1,130.1,130.1,130.1,130.1,130.1,130.1, &
                    132.6,132.6,132.6,132.6,132.6,132.6,132.6,132.6, &
                    132.6,132.6,132.6,132.6,132.6,132.6,132.6,132.6, &
                   0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)
      rxray = (/248.9,248.9,248.9,248.9,248.9,248.9,248.9, &
                248.9,248.9,248.9,248.9,248.9,248.9,248.9,248.9,248.9, &
                248.9,248.9,248.9,248.9,248.9,248.9,248.9,248.9, &
                248.9,248.9,248.9,248.9,248.9,248.9,248.9,248.9, &
                197.6,197.6,197.6,197.6,197.6,197.6,197.6,197.6, &
                197.6,197.6,197.6,197.6,197.6,197.6,197.6,197.6, &
                194.5,194.5,194.5,194.5,194.5,194.5,194.5,194.5, &
                194.5,194.5,194.5,194.5,194.5,194.5,194.5,194.5, &
                   0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)

!---D3D-----------------------D3D----------------------------D3D-----

      fcid= (/1. ,2. ,3. ,4. ,5. ,6. ,7. ,8. ,9. , &
                10.,11.,12.,13.,14.,15.,16.,17.,18./)
      rwstrip1(1)=1.33
      zwstrip1(1)=-1.363
      rwstrip1(2)=1.38
      zwstrip1(2)=-1.363
      rwstrip2(1)=1.4075
      zwstrip2(1)=-1.250
      rwstrip2(2)=1.4575
      zwstrip2(2)=-1.250
      fco2ne=1.0
      end subroutine
