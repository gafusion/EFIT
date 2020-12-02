!
!   $Date: 2008/07/29 23:42:31 $ $Author: radhakri $
!
!   @(#)$RCSfile: basiscomdu.inc,v $ $Revision: 1.1.2.1 $
!
      common/bscom/kppfnc,kppknt,ppknt(npcurn),pptens, &
                   kfffnc,kffknt,ffknt(npcurn),fftens, &
                   kwwfnc,kwwknt,wwknt(npcurn),wwtens, &
                   ppbdry(npcurn),pp2bdry(npcurn), &
                   kppbdry(npcurn),kpp2bdry(npcurn), &
                   ffbdry(npcurn),ff2bdry(npcurn), &
                   kffbdry(npcurn),kff2bdry(npcurn), &
                   wwbdry(npcurn),ww2bdry(npcurn), &
                   kwwbdry(npcurn),kww2bdry(npcurn), &
                   kpphord,kffhord,kwwhord
      common/bscomss/ppbdryss(npcurn),pp2bdryss(npcurn), &
                   kppbdryss(npcurn),kpp2bdryss(npcurn), &
                   ffbdryss(npcurn),ff2bdryss(npcurn), &
                   kffbdryss(npcurn),kff2bdryss(npcurn), &
                   wwbdryss(npcurn),ww2bdryss(npcurn), &
                   kwwbdryss(npcurn),kww2bdryss(npcurn)
      common/autokknot/akchiwt,akerrwt,aktol,akgamwt,akprewt, &
                       kautoknt,kakloop,kakiter
