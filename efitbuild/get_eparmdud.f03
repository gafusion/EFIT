subroutine get_eparmdud_defaults
  use eparmdud129
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

  nwcur2=nwcurn*2
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
  kubics=4
  icycred_loopmax=1290
  nfourier=5

end subroutine

subroutine get_eparmdud_dependents

  use eparmdud129

  nmtark=mse315+mse45+mse15+mse1h+mse315_2+mse210
  nstark=nmtark+libim

  magpol=magpri67+magpri322+magprirdp+magudom
  magpri=magpol+maglds

  nsilop=nsilds+nsilol

  nbwork=nsilop
  msbdry=mbdry+nsilop+nfcoil+1
  msbdr2=2*msbdry

  npcurn=nffcur+nppcur
  necur2=nercur*2
  mfnpcr=nfcoil+npcurn+nvesel+nwwcur+nesum+nfcoil+nercur 
  npcur2=npcurn*2
  nrsmat=nsilop+magpri+nrogow+nffcur+1+npcurn+nwwcur+mpress+nfcoil+nstark+nnece+neceo
  nrsma2=2*nrsmat

  nwcurn=nwwcur+npcurn
  npcur3=npcurn*2

  ncurrt=nvesel+nesum+nfcoil
end subroutine

  subroutine read_eparmdud(filename)
  use eparmdud129
  use var_nio
  implicit none
  integer ::istatus
  character (*) :: filename

  NAMELIST/machinein/nsilds,nsilol,nfcoil,nrogow,nacoil,mfcoil,necoil,nvesel, &
      mpress,nesum,magpri67,magpri322,magprirdp,magudom,maglds,mse315,mse45, &
      mse15,mse1h,mse315_2,mse210,libim,nmsels,nnece,nnecein,neceo,nnnte, &
      ngam_vars,ngam_u,ngam_w,nlimit,nlimbd,nangle,ntangle,nfbcoil,mccoil, &
      micoil,ndata,nwwcur,nffcur,nppcur,nercur,ntime,ndim,kxiter,mqwant, &
      mbdry,mbdry1,nxtram,nxtlim,nco2v,nco2r,modef,modep,modew,kubics, &
      icycred_loopmax,nfourier

  open(unit=nin,status='old',file=filename)
  read (nin,machinein,iostat=istatus)
  close(unit=nin)

  end subroutine

