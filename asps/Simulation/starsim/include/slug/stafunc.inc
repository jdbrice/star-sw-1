* $Id: stafunc.inc,v 1.1.1.1 2004/01/12 23:49:49 potekhin Exp $
* $Name:  $
C Declare types for the things used in the statement function STAFUNC
      INTEGER IIIII,LVPRIN,LVGEOM,LVHIST,LVGRAP,LVDEBU,LVSIMU
     +             ,LVDIGI,LVRECO,LVMFLD,LVANAL,LVBACK,LWPRIN,LWDEBU
C
C Statement Functions for ZEBRA
C datacard value are in LW PRIN/DEBU, current print is in LV PRIN
C
      LVPRIN(IIIII)=IQ(LQ(LKDETM-IIIII)+1)
      LVDEBU(IIIII)=IQ(LQ(LKDETM-IIIII)+2)
      LWPRIN(IIIII)=IQ(LQ(LKDETM-IIIII)+2)     ! For Pavel Nevski to give
      LWDEBU(IIIII)=IQ(LQ(LKDETM-IIIII)+3)     ! a "current" value I guess.
      LVGEOM(IIIII)=IQ(LQ(LKDETM-IIIII)+4)
      LVHIST(IIIII)=IQ(LQ(LKDETM-IIIII)+5)
      LVGRAP(IIIII)=IQ(LQ(LKDETM-IIIII)+6)
      LVSIMU(IIIII)=IQ(LQ(LKDETM-IIIII)+7)
      LVDIGI(IIIII)=IQ(LQ(LKDETM-IIIII)+8)
      LVRECO(IIIII)=IQ(LQ(LKDETM-IIIII)+9)
      LVMFLD(IIIII)=IQ(LQ(LKDETM-IIIII)+10)
      LVANAL(IIIII)=IQ(LQ(LKDETM-IIIII)+11)
      LVBACK(IIIII)=IQ(LQ(LKDETM-IIIII)+12)