ALLOCATE( &
! 
!-------from ecom1-mods.f90---------
   volp(nw),pprime(nw),pres(nw),ffprim(nw), &
   fpol(nw),qpsi(nw),r2surf(nw),rpres(nw),curmid(nw), &
   cjor(nw),bpolss(nw),rmajz0(nw),bprmaj(nw),btrmaj(nw), &
   r1surf(nw),r2surg(nw),bbfpol(nw),ccfpol(nw), &
   ddfpol(nw),rzzmax(nw),zzmax(nw),ermid(nw),eshear(nw), &
   epoten(nw),rhovn(nw),rpmid(nw),xmid(nw),sigrid(nw), &
   sipmid(nw),brhovn(nw),crhovn(nw),drhovn(nw),rhopmid(nw), &
   rgrvt(nw),pwprim(nw),pressw(nw),prwcal(nw),saiprw(nw), &
   rgsvt(nw),cjorec(nw),ffprec(nw),pcurrt(nwnh), &
   pcurrw(nwnh),zero(nwnh),www(nwnh),gecebzpc(nwnh), &
   gecepc(nnece,nwnh),geceppc(nnece,nwnh),gecempc(nnece,nwnh), &
   rbdrpc(mbdry,nwnh),rdjdz(nwnh),grdfdb(nwnh,nfbcoil), &
   rsolpc(mbdry,nwnh), &
   presst(nwnh),rbbbs(npoint),zbbbs(npoint),bpol(npoint), &
   plengt(npoint),bpolz(npoint), pcurrtpp(nwnh), &
! OPT_INPUT >>>
   inpfile_in(ntime), &
! OPT_INPUT <<<
!
!-------from ecom2-mods.f90--------
   gridec(nwnh,nesum),gridvs(nwnh,nvesel),gbrpc(nstark,nwnh), &
   gbzpc(nstark,nwnh),gridac(nwnh,nacoil), &
!
!-------from modules-efitx.f90-------   
   beti(icycred_loopmax,nw-2),abeti(icycred_loopmax,nw-2), &
   wk1(icycred_loopmax,nw-2),rhsdumy1(nwnh),phi(nw),v(nw), &
   wk2(nw),diagl(nw),diagu(nw),psi(nwnh),xpsi(nwnh), &
   vfbrrt(nwnh),psipla(nwnh),xout(npoint),yout(npoint), &
!
!-------from commonblocks module------
   c(kubicx,lubicx,kubicy,lubicy),wk(nwrk), &
   copy(nw,nh),bkx(lubicx+1),bky(lubicy+1), &
   cw(kubicx,lubicx,kubicy,lubicy),wkw(nwrk), &
   copyw(nw,nh),bwx(lubicx+1),bwy(lubicy+1), &
   sifprw(nw),bwprw(nw),cwprw(nw),dwprw(nw), &
   psirz(nw,nh),sfprw(nw),sprwp(nw), &
   wgridpc(nwnh),rfcpc(nfcoil,nwnh), &
   ct(kubicx,lubicx,kubicy,lubicy), &
   wkt(nwrk),bkrt(lubicx+1),bkzt(lubicy+1), &
   psiold(nwnh),psipold(nwnh),psipp(nwnh), &
   work(nwnh),sifpre(nw),bwpre(nw),cwpre(nw), &
   dwpre(nw),sfpre(nw),sprep(nw),worka(nwf),worka2(nw), &
   zeros(nwnh),byringr(nh2),byringz(nh2), &
   cjrf(nwf+5*nxtrap*nxtram + 2*nxtram - 4*npoint), &
   cj(kubicx,lubicx,kubicy,lubicy),wkj(nwrk), &
   copyj(nw,nh),bjx(lubicx+1),bjy(lubicy+1), &
   cv(kubicx,lubicx,kubicy,lubicy),wkv(nwrk), &
   copyv(nw,nh),bvx(lubicx+1),bvy(lubicy+1), &
   xsisii(nw),f_a(nw),f_b(nw),f_c(nw),f_d(nw), &
   pp_a(nw),pp_b(nw),pp_c(nw),pp_d(nw), &
   chi_c(kubicx,lubicx,kubicy,lubicy),chi_wk(nwrk), &
   chi_copy(nw,nh),chi_bkx(lubicx+1),chi_bky(lubicy+1), &
   wxin(npoint),wyin(npoint),wxout(npoint),wyout(npoint), &
   xouts(npoint),youts(npoint),bpoo(npoint),bpooz(npoint), &
   bpooc(npoint),bfpol(npoint),cfpol(npoint),dfpol(npoint), &
   rsplt(npoint),zsplt(npoint),csplt(npoint),xxtra(nxtrap,nxtram), &
   yxtra(nxtrap,nxtram),bpxtra(nxtrap,nxtram),flxtra(nxtrap,nxtram), &
   fpxtra(nxtrap,nxtram) &
)
