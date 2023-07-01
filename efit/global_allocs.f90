 
! 
!-------from ecom1-mods.f90---------
ALLOCATE(volp(nw),pprime(nw),pres(nw),ffprim(nw), &
   fpol(nw),qpsi(nw),r2surf(nw),rpres(nw),curmid(nw), &
   cjor(nw),bpolss(nw),rmajz0(nw),bprmaj(nw),btrmaj(nw), &
   r1surf(nw),r2surg(nw),bbfpol(nw),ccfpol(nw), &
   ddfpol(nw),rzzmax(nw),zzmax(nw),ermid(nw),eshear(nw), &
   epoten(nw),rhovn(nw),rpmid(nw),xmid(nw),sigrid(nw), &
   sipmid(nw),brhovn(nw),crhovn(nw),drhovn(nw),rhopmid(nw), &
   rgrvt(nw),pwprim(nw),pressw(nw),prwcal(nw), &
   rgsvt(nw),cjorec(nw),ffprec(nw),pcurrt(nwnh), &
   pcurrw(nwnh),zero(nwnh),www(nwnh),gecebzpc(nwnh), &
   gecepc(nnece,nwnh),geceppc(nnece,nwnh),gecempc(nnece,nwnh), &
   rbdrpc(mbdry,nwnh),rdjdz(nwnh),grdfdb(nwnh,nfbcoil), &
   rsolpc(mbdry,nwnh), &
   presst(nwnh),rbbbs(npoint),zbbbs(npoint),bpol(npoint), &
   plengt(npoint),bpolz(npoint),pcurrtpp(nwnh))
     
ALLOCATE(racoil(nacoil),zacoil(nacoil), wacoil(nacoil), hacoil(nacoil), &
   r1sdry(ndata),r2sdry(ndata),vzeroj(ndata),sizeroj(ndata), &
   rqajx(nppcur),rjjjx(nppcur),rqafx(nffcur),rjjfx(nffcur), &
   rqawx(nwwcur),rjjwx(nwwcur), &
   brsp(nrsmat),chisq(ntime), &
   chi2rm(ntime),bcentr(ntime),fwtsi(nsilop),rwtsi(nsilop), &
   fwtmp2(magpri),rwtmp2(magpri),fwacoil(nacoil),fwtfc(nfsum), &
   cbrspv(npcurn),silopt(ntime,nsilop),expmpi(ntime,magpri), &
   accurt(ntime,nacoil),fccurt(ntime,nfsum),eccurt(ntime,nesum), &
   denvt(ntime,nco2v),denrt(ntime,nco2r),iccurt(ntime,micoil), &
   ipmeas(ntime),time(ntime),pbinj(ntime),ierpsi(nsilop), &
   ierpr(nrogow),iermpi(magpri),psibit(nsilop),prbit(nrogow), &
   bitmpi(magpri),vloopt(ntime),psiref(ntime),diamag(ntime), &
   sigdia(ntime),swtsi(nsilop),swtmp2(magpri), &
   swtfc(nfsum),swtec(nesum),bitfc(nfsum),ierfc(nfsum),bitec(nesum), &
   ierec(nesum),bitic(micoil),ieric(micoil), &
   tangam(ntime,nstark),siggam(ntime,nstark),a1gam(ntime,nstark), &
   a2gam(ntime,nstark),a3gam(ntime,nstark),a4gam(ntime,nstark),tangam_uncor(ntime,nstark), &
   fwtgam(nstark),chigam(nstark),swtgam(nstark),iergam(nstark),mseport(nstark),mse_spave_on(nstark))
   
ALLOCATE(rsilpe(nsilop),rmp2pe(magpri),rprepe(mpress),gbdrpe(mbdry), &
   rgampe(nstark),rbrpe(nstark),rbzpe(nstark),rmlspe(nmsels),relspe(nmsels), &
   bmselt(ntime,nmsels),sbmselt(ntime,nmsels),fwtbmselt(ntime,nmsels), &
   rrmselt(ntime,nmsels),zzmselt(ntime,nmsels),l1mselt(ntime,nmsels),l2mselt(ntime,nmsels), &
   l4mselt(ntime,nmsels),emselt(ntime,nmsels),semselt(ntime,nmsels), &
   fwtemselt(ntime,nmsels),swtbmselt(ntime,nmsels),swtemselt(ntime,nmsels), &
   cmmls(ntime,nmsels),cmels(ntime,nmsels),rhsmls(ntime,nmsels),rhsels(ntime,nmsels), &
   l3mselt(ntime,nmsels),cmmls2(ntime,nmsels),cmmlsv(ntime,nmsels), &
   swtbmsels(nmsels),swtemsels(nmsels),chimls(nmsels),chiels(nmsels), &
   rmlspc(nmsels,npcurn),rmlsec(nmsels,nesum),relsec(nmsels,nesum),rmlsvs(nmsels,nvesel), &
   relsvs(nmsels,nvesel),simls(nmsels),sinmls(nmsels),btmls(nmsels),brmls(nmsels), &
   bzmls(nmsels),chi2mls(ntime),chi2gamt(ntime),iermselt(ntime,nmsels),kerrot(ntime), &
   rsilfe(nsilop),rmp2fe(magpri),gbdrfe(mbdry), &
   rgamfe(nstark),rbrfe(nstark),rbzfe(nstark),rmlsfe(nmsels),relsfe(nmsels), &
   teecein(nnecein),feece(nnecein),errorece(nnecein), &
   teecein0(nnecein),feece0(nnecein),errorece0(nnecein), &
   becein(nnecein),recein(nnecein),teeceinr(nnecein), &
   recem(nnece),recep(nnece),fwtece(nnece),fwtece0(nnece), &
   swtece(nnece),chiece(nnece),ecefit(nnece),ecebit(nnece), &
   recebzfc(nfsum),recebzec(nesum),recefc(nnece,nfsum),receec(nnece,nesum), &
   ierece(nnece),recepc(nnece,npcurn),brspece(ntime,nnece),recebzpc(npcurn), &
   brspecebz(ntime),cmecebz(ntime),cmece(nnece,ntime),recevs(nnece,nvesel), &
   recebzvs(nvesel),recedz(nnece),gecedz(nnece),rtep(nnece),rtem(nnece),rpbit(nnece),rmbit(nnece), &
   receac(nnece,nacoil),recebzac(nacoil),teecer(nnnte),rrr(nnnte),bbf(nnnte),teeceb(nnnte), &
   receoi(kxiter),recemi(kxiter,nnece),recepi(kxiter,nnece),iwp(nnece),iwm(nnece))
   
ALLOCATE(rbrpc(nstark,npcurn),rbzpc(nstark,npcurn),rgampc(nstark,npcurn),rbrfc(nstark,nfsum), &
   rbzfc(nstark,nfsum),rbrec(nstark,nesum),rbzec(nstark,nesum),rgamec(nstark,nesum), &
   rbrvs(nstark,nvesel),rbzvs(nstark,nvesel),rgamvs(nstark,nvesel),rgamfc(nstark,nfsum), &
   rhsgam(ntime,nstark),rrgam(ntime,nstark),zzgam(ntime,nstark),starkar(ntime,nstark), &
   starkaz(ntime,nstark),a5gam(ntime,nstark),a6gam(ntime,nstark),a7gam(ntime,nstark),a8gam(ntime,nstark), &
   cmgam(nstark,ntime),spatial_fix(nstark,ntime),rgamac(nstark,nacoil),rbrac(nstark,nacoil), &
   rbzac(nstark,nacoil),btgam(nstark),sistark(nstark),qstark(nstark), &
   rmse_gain(nstark),rmse_slope(nstark),rmse_scale(nstark),rmse_offset(nstark),sigam(nstark), &
   bzmse(nstark),bzmsec(nstark),cjmse(nstark),cjmsec(nstark),rhogam(nstark), &
   spatial_avg_gam(nstark,ngam_vars,ngam_u,ngam_w), &
   saisil(nsilop),saimpi(magpri),saipre(mpress),chifcc(nfsum), &
   dfluxc(ntime),cdflux(ntime),edflux(ntime),rspdlc(nffcur), &
   rco2r(nco2r,ntime),rco2v(nco2v,ntime),chordv(nco2v),chordr(nco2r), &
   dco2r(ntime,nco2r),dco2v(ntime,nco2v),erflag(ntime,nflag), &
   cerror(kxiter),csibry(kxiter),csimag(kxiter),cvolp(kxiter),crmaxi(kxiter), &
   czmaxi(kxiter),cemaxi(kxiter),cqmaxi(kxiter),cchisq(kxiter), &
   brfbc(kxiter),tvfbrt(kxiter),cdelz(kxiter), &
   rxray(nangle+ntangle),zxray(nangle+ntangle),xangle(nangle+ntangle), &
   rzeroj(ndata),fwtpre(mpress),vforcep(nvesel),vforcet(nvesel),fwtfcsum(nfsum), &
   limloc(ntime),vsname(nvesel), mpnam2(magpri),lpname(nsilop), &
   qsiw(mqwant),pasmsw(mqwant),fqsiw(mqwant),siwantq(mqwant),fgowsw(npcurn,mqwant), &
   rsi(nsilop),zsi(nsilop),wsi(nsilop),hsi(nsilop),as(nsilop),as2(nsilop), &
   rf(nfcoil),zf(nfcoil),wf(nfcoil),hf(nfcoil),af(nfcoil),af2(nfcoil), &
   turnfc(nfsum),fcid(nfcoil),fcturn(nfcoil), &
   re(necoil),ze(necoil),we(necoil),he(necoil),ecid(necoil), &
   xmp2(magpri),ymp2(magpri),amp2(magpri),smp2(magpri),patmp2(magpri), &
   xlim(nlimit),ylim(nlimit),xlim_180(nlimit),ylim_180(nlimit), &
   xmim(nlimit),ymim(nlimit),xmim_180(nlimit),ymim_180(nlimit), &
   ipsi(ntime),irogw(ntime),imag2(ntime),iplasm(ntime), &
   idlopc(ntime),ifc(ntime),iec(ntime), csumip(kxiter),tratio(kxiter), &
   aveerr(kxiter),iermax(kxiter),jermax(kxiter))
   
ALLOCATE(rbdry(mbdry),zbdry(mbdry),fwtbdry(mbdry),fwtbry(mbdry),sigrbd(mbdry), &
   sigzbd(mbdry),rbdry0(mbdry),zbdry0(mbdry), &
   rsol(mbdry),zsol(mbdry),fwtsol(mbdry),fwtsolw(mbdry), &
   rbdrfc(mbdry,nfsum),rsolfc(mbdry,nfsum),rbdrac(mbdry,nacoil), &
   gbdrpc(mbdry,npcurn),gmp2dz(magpri),gsildz(nsilop),gbrdz(nstark), &
   gbzdz(nstark),rgamdz(nstark),rmlsdz(nmsels),relsdz(nmsels),rpredz(mpress), &
   rprwdz(mpress),gbdrdz(mbdry),erbloc(msbdry), erbsloc(msbdry), &
   fcsum(nfsum),fczero(nfsum),ivbcuse(ndata), &
   pressr(mpress),rpress(mpress),zpress(mpress),tethom(mpress),rteth(mpress), &
   zteth(mpress),tionex(mpress),rion(mpress),zion(mpress),dnethom(mpress), &
   rneth(mpress),zneth(mpress),pbeam(mpress),sibeam(mpress),sgteth(mpress), &
   sigti(mpress),sigpre(mpress),sgneth(mpress),precal(mpress), &
   dnbeam(mpress),dmass(mpress),scalepr(mpress),premea(mpress),saipre2(mpress), &
   pbimth(ndata),bremin(ndata),bremsig(ndata),brmrtan(ndata), &
   brmzelev(ndata),dnbthom(ndata),keebdry(nercur),kee2bdry(nercur), &
   rgamer(nstark,nercur),rmlser(nmsels,nercur),relser(nmsels,nercur), &
   eeknt(nercur),eebdry(nercur),ee2bdry(nercur),cerer(nercur), ermse(nstark), &
   e1rbz(nstark,nercur),e2rbz(nstark,nercur),e3rbr(nstark,nercur), &
   calpa(nppcur,nppcur),xalpa(nppcur),alpax(nppcur),cgama(nffcur,nffcur), &
   xgama(nffcur),gamax(nffcur),comega(nwcurn,nwcurn),xomega(nwcurn), &
   ccoils(nfsum,nfsum),xcoils(nfsum),cloops(nsilop,nsilop),xloops(nsilop), &
   tifit(nppcur),stitho(ndata),xsiion(ndata),tithom(ndata),tefit(nppcur), &
   snitho(ndata),dnitho(ndata),defit(nppcur),zuperts(ntime),rlibim(ntime), &
   xlimbd(nlimbd),ylimbd(nlimbd),fb_plasma(ntime),signn1(magpri), curtn1(ntime), &
   curc139(ntime),curc79(ntime),curc199(ntime),curiu30(ntime), &
   curiu90(ntime),curiu150(ntime),curil30(ntime),curil90(ntime),curil150(ntime), &
   curccoi(ntime,mccoil),curicoi(ntime,micoil))
   
ALLOCATE(xnodef(modef),xnodep(modep),xnodew(modew),kbasef(modef), &
   kbasep(modep),kbasew(modew),omegat(mpress),rpresw(mpress),zpresw(mpress), &
   presw(mpress),sigprw(mpress),rpresws(mpress),fwtprw(mpress),romegat(mpress), &
   zomegat(mpress),sigome(mpress),rpreswv(mpress),scalepw(mpress), &
   saiprw(npress),premew(mpress),saiprw2(mpress),rprwpc(mpress,nwwcur), &
   betapw(ntime),betatw(ntime),wplasw(ntime), &
   fwtec(nesum),cecurr(nesum),chiecc(nesum), &
   ztssym(ntime),ztswid(ntime),ptssym(ntime),ztserr(ntime), &
   psin32(ntime),psin21(ntime),rq32in(ntime),rq21top(ntime), &
   sigfcc(nfsum),sigecc(nesum),sigsil(nsilop),sigmpi(magpri), &
   signbt(magpri))
   
!-------from ecom2-mods.f90--------
ALLOCATE(gridec(nwnh,nesum),gridvs(nwnh,nvesel),gbrpc(nstark,nwnh), &
   gbzpc(nstark,nwnh),gridac(nwnh,nacoil), &
   rsilec(nsilop,nesum),rmp2ec(magpri,nesum),rfcec(nfsum,nesum),ecurrt(nesum), &
   pecur(nesum),rbdrec(mbdry,nesum),rsolec(mbdry,nesum),recec(nesum,nesum), &
   vcurrt(nvesel),rsilvs(nsilop,nvesel),rmp2vs(magpri,nvesel), &
   rfcvs(nfsum,nvesel),rbdrvs(mbdry,nvesel),rsolvs(mbdry,nvesel), &
   rvsfc(nvesel,nfsum),rvsec(nvesel,nesum),rvsvs(nvesel,nvesel), &
   rsilpc(nsilop,npcurn),fgowpc(npcurn),rmp2pc(magpri,npcurn),rprepc(mpress,nppcur), &
   rsilac(nsilop,nacoil),rmp2ac(magpri,nacoil),thetav(nvesel), &
   sinta(nfourier,nvesel),costa(nfourier,nvesel),vecta(2*nfourier+1,nvesel))

ALLOCATE(ppbdry(npcurn),pp2bdry(npcurn), &
   kppbdry(npcurn),kpp2bdry(npcurn), &
   ffbdry(npcurn),ff2bdry(npcurn), &
   kffbdry(npcurn),kff2bdry(npcurn), &
   wwbdry(npcurn),ww2bdry(npcurn), &
   kwwbdry(npcurn),kww2bdry(npcurn), &
   wwknt(npcurn),ffknt(npcurn),ppknt(npcurn), &
   appknt(npcurn),affknt(npcurn),awwknt(npcurn),aeeknt(npcurn), &
   save_gam(ntime,nstark),save_tangam(ntime,nstark), &
   arsp_cw2(ndata,nppcur),wrsp_cw2(nppcur),work_cw2(ndata), &
   bdata_cw2(ndata),ematrix_cw2(nppcur,nppcur),einv_cw2(nppcur,nppcur))

ALLOCATE(rmx(mfila),zmx(mfila),rsilpf(nsilop,mfila), &
   rmp2pf(magpri,mfila), &
   irfila(mfila),jzfila(mfila), &
   wsilpc(nsilop),wmp2pc(magpri),wfcpc(nfsum),wgridpc(nwnh), &
   wecpc(nesum),wvspc(nvesel),npxtra(nxtram),scraps(nxtram), &
   workb_jw4(nsilop))
  
!-------from modules-efit.f90-------   
ALLOCATE(beti(icycred_loopmax,nw-2),abeti(icycred_loopmax,nw-2), &
   wk1(icycred_loopmax,nw-2),rhsdumy1(nwnh),phi(nw),v(nw), &
   wk2(nw),diagl(nw),diagu(nw),psi(nwnh),xpsi(nwnh), &
   vfbrrt(nwnh),psipla(nwnh),xout(npoint),yout(npoint))



ALLOCATE(volecs(nesum),volecc(nesum),rsisec(nesum),volfcs(nfsum),volfcc(nfsum), &
   rvs(nvesel),zvs(nvesel),hvs(nvesel),wvs(nvesel),avs(nvesel),avs2(nvesel), &
   rsisvs(nvesel),alphab(icycred_loopmax),diag1(icycred_loopmax), &
   tempgrid(ncurrt),tempgrid2(ncurrt),rowscale(nrsmat),colscale(mfnpcr), &
   elong(ntime),rout(ntime),zout(ntime),utri(ntime), &
   ltri(ntime),aminor(ntime),volume(ntime),betat(ntime),gaptop(ntime), &
   betap(ntime),li(ntime),gapin(ntime),gapout(ntime),qstar(ntime), &
   rcurrt(ntime),zcurrt(ntime),qout(ntime),sepin(ntime), &
   sepout(ntime),septop(ntime),sibdry(ntime),area(ntime), &
   wmhd(ntime),elongm(ntime),qm(ntime),terror(ntime), &
   rm(ntime),zm(ntime),gapbot(ntime),sepbot(ntime), &
   alpha(ntime),rttt(ntime),dbpli(ntime),delbp(ntime),oring(ntime), &
   sepexp(ntime),shearb(ntime), &
   xtch(ntime),ytch(ntime),q95(ntime),vertn(ntime),aaq1(ntime), &
   aaq2(ntime),aaq3(ntime),btaxp(ntime),btaxv(ntime), &
   psim(ntime),jerror(ntime),dsep(ntime),peak(ntime), &
   wbpol(ntime),taumhd(ntime),betapd(ntime),betatd(ntime), &
   wdia(ntime),taudia(ntime),wbpold(ntime), &
   qmerci(ntime),slantu(ntime),slantl(ntime),zeff(ntime), &
   zeffr(ntime),tave(ntime),rvsin(ntime),zvsin(ntime), &
   rvsout(ntime),zvsout(ntime),wpdot(ntime),wbdot(ntime), &
   vsurfa(ntime),cjor95(ntime),pp95(ntime),drsep(ntime), &
   yyy2(ntime),xnnc(ntime),wtherm(ntime),wfbeam(ntime),taujd3(ntime),tauthn(ntime), &
   li3(ntime),tflux(ntime),twagap(ntime),rseps(2,ntime),zseps(2,ntime))

ALLOCATE(csilop(nsilop,ntime),csilopv(nsilop,ntime),crogow(nrogow,ntime), &
   cmpr2(magpri,ntime),cmpr2v(magpri,ntime),ipmhd(ntime),indent(ntime), &
   ccbrsp(nfsum,ntime),caccurt(ntime,nacoil),qsiwant(ntime), &
   cjorsw(ntime),cjor0(ntime),ssiwant(ntime),ssi95(ntime),cjor99(ntime), &
   cj1ave(ntime),rmidin(ntime),rmidout(ntime),psurfa(ntime), &
   dolubaf(ntime),dolubafm(ntime),diludom(ntime), &
   diludomm(ntime),dminux(ntime),dminlx(ntime), &
   ratsol(ntime),rvsiu(ntime),zvsiu(ntime),rvsou(ntime), &
   zvsou(ntime),rvsid(ntime),zvsid(ntime),rvsod(ntime),zvsod(ntime), &
   vtime(ntime), &
   s1(ntime),s2(ntime),s3(ntime),bpolav(ntime),rfcfc(nfsum,nfsum))
!
!-------from commonblocks module------
ALLOCATE(c(kubicx,lubicx,kubicy,lubicy),wk(nwrk), &
   copy(nw,nh),bkx(lubicx+1),bky(lubicy+1), &
   cw(kubicx,lubicx,kubicy,lubicy),wkw(nwrk), &
   copyw(nw,nh),bwx(lubicx+1),bwy(lubicy+1), &
   sifprw(nw),bwprw(nw),cwprw(nw),dwprw(nw), &
   sfprw(nw),sprwp(nw), &
   rfcpc(nfsum,nwnh), &
   ct(kubicx,lubicx,kubicy,lubicy), &
   wkt(nwrk),bkrt(lubicx+1),bkzt(lubicy+1), &
   psiold(nwnh),psipold(nwnh), &
   work(nwnh),sifpre(nw),bwpre(nw),cwpre(nw), &
   dwpre(nw),sfpre(nw),sprep(nw),worka(nwf),worka2(nw), &
   zeros(nwnh),byringr(nh2),byringz(nh2), &
   cjrf(nwf+5*nxtrap*nxtram+2*nxtram-4*npoint), &
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
   fpxtra(nxtrap,nxtram))

