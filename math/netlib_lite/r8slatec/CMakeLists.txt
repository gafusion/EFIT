######################################################################
#
# CMakeLists.txt for r8slatec
#
# $Id: CMakeLists.txt 53 2011-09-29 18:28:47Z cary $
#
######################################################################

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

set(R8SLATEC_SOURCES
  d1mach.c   j1mach.c
  fmin.f90
  dgbdi.f90 dbdiff.f dbsppp.f dei.f dgbco.f dintrv.f dp1vlu.f dpcoef.f
  dpolft.f fdump.f j4save.f parmsetget.f r8genbun.f r8hwscyl.f
  r8merge.f r8pimach.f r8poisp2.f xerclr.f xercnt.f xerhlt.f
  xermax.f xermsg.f xerprn.f xersve.f xgetf.f xgetua.f
  xsetf.f xsetfp.f xsetun.f xsetunp.f zexp.f zmlt.f zshch.f
  zuchk.f

  d9aimp.F d9atn1.F d9b0mp.F d9b1mp.F d9chu.F d9gmic.F d9gmit.F
  d9knus.F d9lgic.F d9lgit.F d9lgmc.F d9ln2r.F d9pak.F d9upak.F dacosh.F
  dai.F daie.F dasinh.F dasyik.F dasyjy.F datanh.F dbesi.F
  dbesi0.F dbesi1.F dbesj.F dbesj0.F dbesj1.F dbesk.F dbesk0.F
  dbesk1.F dbesks.F dbesy.F dbesy0.F dbesy1.F dbeta.F dbetai.F
  dbi.F dbie.F dbinom.F dbkias.F dbkisr.F dbsgq8.F dbsi0e.F
  dbsi1e.F dbsk0e.F dbsk1e.F dbskes.F dbskin.F dbsknu.F dbspdr.F
  dbspev.F dbspvd.F dbspvn.F dbsqad.F dbsynu.F dbvalu.F dcbrt.F
  dchu.F dcosdg.F dcot.F dcsevl.F ddaws.F de1.F derf.F derfc.F
  dexint.F dexprl.F dfac.F dgami.F dgamic.F dgamit.F dgamlm.F
  dgamln.F dgamma.F dgamr.F dgamrn.F dgbfa.F dgbsl.F dhkseq.F
  dgedi.F dtrco.F dgeco.F90
  djairy.F dlbeta.F dlgams.F dli.F dlngam.F dlnrel.F dpoch.F
  dpoch1.F dpsi.F dpsifn.F dpsixn.F dqmomo.F drc.F drc3jj.F
  drc3jm.F drc6j.F drd.F drf.F drj.F dsindg.F dspenc.F dxadd.F
  dxadj.F dxlegf.F dxnrmp.F dxpmu.F dxpmup.F dxpnrm.F dxpqnu.F
  dxpsi.F dxqmu.F dxqnu.F dxred.F dxset.F dyairy.F initds.F
  r8cosgen.F r8poisd2.F r8poisn2.F r8tri3.F r8trix.F zabs.F
  zacai.F zacon.F zairy.F zasyi.F zbesh.F zbesi.F zbesj.F zbesk.F
  zbesy.F zbinu.F zbiry.F zbknu.F zbuni.F zbunk.F zdiv.F zkscl.F
  zlog.F zmlri.F zrati.F zs1s2.F zseri.F zsqrt.F zunhj.F zuni1.F
  zuni2.F zunik.F zunk1.F zunk2.F zuoik.F zwrsk.F
)

add_library(r8slatec ${R8SLATEC_SOURCES})

install(TARGETS r8slatec
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

