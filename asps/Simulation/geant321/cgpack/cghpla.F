*
* $Id: cghpla.F,v 1.1.1.1 2004/01/15 00:11:56 potekhin Exp $
*
* $Log: cghpla.F,v $
* Revision 1.1.1.1  2004/01/15 00:11:56  potekhin
*
* Revision 1.1.1.1  1997/11/03 15:24:27  atlascvs
* Importing CERNLIB version 08.21.
*
* Revision 1.1.1.1  1995/10/24 10:19:44  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.31  by  S.Giani
*-- Author :
      SUBROUTINE CGHPLA(IFACE,XYZ,ABCD)
************************************************************************
*                                                                      *
*     Name: CGHPLA                                                     *
*     Author: E. Chernyaev                       Date:    08.08.88     *
*                                                Revised:              *
*                                                                      *
*     Function: Compute face plane equation coefficients:              *
*                    Ax + By + Cz + D = 0                              *
*                                                                      *
*     References: none                                                 *
*                                                                      *
*     Input: IFACE(*) - face                                           *
*            XYZ(3,*) - node coordinates                               *
*                                                                      *
*     Output: ABCD(4) - plane equation coefficients                    *
*                                                                      *
*     Errors: none                                                     *
*                                                                      *
************************************************************************
      REAL            XYZ(3,*),ABCD(4),GRAV(3)
*SG
      INTEGER         IFACE(*)
*SG
*-
      NEDGE = IFACE(1)
      A      = 0.
      B      = 0.
      C      = 0.
      GRAV(1)= 0.
      GRAV(2)= 0.
      GRAV(3)= 0.
      JF     = 2
      DO 100 NE=1,NEDGE
        N1     = IFACE(JF)
        N2     = IFACE(JF+1)
        JF     = JF + 2
        A      = A + XYZ(2,N1)*XYZ(3,N2) - XYZ(2,N2)*XYZ(3,N1)
        B      = B + XYZ(3,N1)*XYZ(1,N2) - XYZ(3,N2)*XYZ(1,N1)
        C      = C + XYZ(1,N1)*XYZ(2,N2) - XYZ(1,N2)*XYZ(2,N1)
        GRAV(1)= GRAV(1) + XYZ(1,N1) + XYZ(1,N2)
        GRAV(2)= GRAV(2) + XYZ(2,N1) + XYZ(2,N2)
        GRAV(3)= GRAV(3) + XYZ(3,N1) + XYZ(3,N2)
  100   CONTINUE
      HNGRAV = 0.5/NEDGE
      GRAV(1)= GRAV(1) * HNGRAV
      GRAV(2)= GRAV(2) * HNGRAV
      GRAV(3)= GRAV(3) * HNGRAV
      IF (ABS(A) .LT. 0.000001) A=0.
      IF (ABS(B) .LT. 0.000001) B=0.
      IF (ABS(C) .LT. 0.000001) C=0.
      AREAI  = 1./SQRT(A*A + B*B + C*C)
      ABCD(1)= A * AREAI
      ABCD(2)= B * AREAI
      ABCD(3)= C * AREAI
      ABCD(4)=-(A*GRAV(1) + B*GRAV(2) + C*GRAV(3)) * AREAI
      END
