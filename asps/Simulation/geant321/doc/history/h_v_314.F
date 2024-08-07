*
* $Id: h_v_314.F,v 1.1.1.1 2004/01/15 00:11:58 potekhin Exp $
*
* $Log: h_v_314.F,v $
* Revision 1.1.1.1  2004/01/15 00:11:58  potekhin
*
* Revision 1.1.1.1  1997/11/03 15:28:39  atlascvs
* Importing CERNLIB version 08.21.
*
* Revision 1.1.1.1  1995/10/24 10:19:48  cernlib
* Geant
*
*
#include "sys/CERNLIB_machine.h"
#include "pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.37  by  S.Giani
*-- Author :
C   09-FEB-1988 NVE CERN GENEVA
C               NUMBER OF SUBSCRIPTS OF ARRAY "PCM" CHANGED FROM 1
C               TO 2 IN EQUIVALENCE STATEMENT IN ROUTINE PHASP.
C   17-FEB-1988 NVE CERN GENEVA
C               ALL CALLS TO "RN" CHANGED TO "RNDM" AND WITHIN A DO
C               LOOP THE DUMMY ARGUMENT HAS BEEN SET TO THE LOOP
C               INDEX IN ORDER TO ENABLE COMPILATION WITH OPTIMIZATION
C   25-FEB-1988 NVE CERN GENEVA
C               STATEMENT FUNCTION SEQUENCES "ANPOL" AND "ANPOL" ADDED
C   25-FEB-1988 NVE CERN GENEVA
C               SEQUENCE "ANPOL" INCLUDED IN "STEEQ" AND "STPAIR"
C   26-FEB-1988 NVE CERN GENEVA
C               SEQUENCE "LIMITS" ADDED FOR ARGUMENTS OF INT. FUNCTIONS
C   01-MAR-1988 NVE CERN GENEVA
C               ALL CALLS TO "HSTORY" REMOVED
C   04-MAR-1988 NVE CERN GENEVA
C               +SELF MATERIAL ADDED TO ENABLE SUBROUTINE TIMING
C   16-MAR-1988 NVE CERN GENEVA
C               NEW "SAVE" CARD AND INIT. FOR "PMUL" AND "ANORM"
C               IN CASCADE ROUTINES
C   16-MAR-1988 NVE CERN GENEVA
C               DIMENSION OF "ANORM1" CHANGED FROM (4,60) TO (2,60)
C               IN ROUTINE CASASM
C   16-MAR-1988 NVE CERN GENEVA
C               "REAL*8" AND "REAL*4" REPLACED BY "DOUBLE PRECISION"
C               AND "REAL" RESP. IN ROUTINE "PHASP"
C   16-MAR-1988 NVE CERN GENEVA
C               "REAL*8" REMOVED IN ROUTINE "PDK"
C               FOR DOUBLE PRECISION USE "DPDK"
C   16-MAR-1988 NVE CERN GENEVA
C               "SQRT(A)" REPLACED BY "SQRT(ABS(A))" IN "GENXPT"
C   18-MAR-1988 NVE CERN GENEVA
C               GHEISHA UTILITIES ADDED IN PATCH "GHUTILS"
C   18-MAR-1988 NVE CERN GENEVA
C               NEW PATCHES "TIFACE, TATINA, GIFACE AND CASCAD"
C               ADDED TO BE COMPLETED IN FUTURE
C   29-MAR-1988 NVE CERN GENEVA
C               DECK "VECSUB" CLEANED UP AND ALL "ENTRIES" RE-WRITTEN
C               AND CHANGED TO "SUBROUTINES"
C   29-MAR-1988 NVE CERN GENEVA
C               "DCOS" "DSIN" AND "DSQRT" CHANGED TO "COS" "SIN"
C               AND "SQRT" RESP. IN ROUTINE "PHASP"
C   29-MAR-1988 NVE CERN GENEVA
C               "DPDK" "DLPSOR" AND "DOTES2" PUT IN SEPARATE DECKS
C   29-MAR-1988 NVE CERN GENEVA
C               "REAL*8" CHANGED TO "DOUBLE PRECISION" IN SEQ. /GENIO
C   30-MAR-1988 NVE CERN GENEVA
C               PROTECTIONS AGAINST "EXP()" LIMITS ADDED IN REMAINING
C               CASCADE ROUTINES
C   30-MAR-1988 NVE CERN GENEVA
C               PROTECTIONS AGAINST "DIVISION BY ZERO" ADDED IN
C               ROUTINE "TWOCLU"
C   06-APR-1988 NVE CERN GENEVA
C               "GPGHEI" CHANGED IN ORDER TO MATCH GEANT 3.12
C   07-APR-1988 NVE CERN GENEVA
C               SPECIFIC FUNCTIONS "AMAX1" AND "ALOG" REPLACED BY
C               THEIR GENERIC EQUIVALENTS AND
C               GEOMETRICAL CROSS-SECTIONS FOR NUCLEAR FRAGMENTS ADDED
C               IN ROUTINE "GHESIG"
C   11-APR-1988 NVE CERN GENEVA
C               PROTECTIONS AGAINST "EXP()" LIMITS ADDED IN ROUTINE
C               "CINEMA" AND "ALOG" CALLS CHANGED TO CALLS TO
C               GENERIC FUNCTION "LOG"
C   11-APR-1988 NVE CERN GENEVA
C               PROTECTIONS AGAINST NEGATIVE SQRT ADDED IN
C               ROUTINE "TWOCLU"
C   14-APR-1988 NVE CERN GENEVA
C               ARGUMENT OF "RNDM" SET TO LOOP INDEX
C               IN ROUTINE "NORMAL"
C   22-APR-1988 NVE CERN GENEVA
C               ALL CALLS TO "ALOG" "IABS" "DCOS" AND "DEXP"
C               CHANGED TO THE GENERIC EQUIVALENTS
C   02-MAY-1988 NVE CERN GENEVA
C               ROUTINES "GENXPT" AND "TWOCLU" CORRECTED FOR LOW
C               ENERGY BEHAVIOUR
C   02-MAY-1988 R. BRUN/NVE CERN GENEVA
C               STMT. "GO TO 38" CHANGED TO "GO TO 30" IN ROUTINE
C               "GENXPT" TO PREVENT INFINITE LOOPING
C   02-MAY-1988 NVE CERN GENEVA
C               ALL "IPV(..)" STUFF CHANGED TO NORMAL "PV(..)"
C               IN ROUTINE "ADD"
C   04-MAY-1988 NVE CERN GENEVA
C               FUNCTIONS "RMASS" AND "CHARGE" REPLACED BY ARRAYS
C               "RMASS" AND "RCHARG" RESP. AND PUT INTO COMMON /CONSTS
C               MASSES AND CT VALUES OF PARTICLES UPDATED TO LATEST
C               (= GEANT) VALUES AND INITIALIZED IN ROUTINE "GHEINI"
C   04-MAY-1988 NVE CERN GENEVA
C               DISTINCTION INTRODUCED BETWEEN P AND N TARGET IN ALL
C               "CAS..." ROUTINES AND S AND RS COMPUTED ACCORDINGLY
C               IN THE ROUTINES THEMSELVES
C               THE SAME HAS BEEN DONE IN ROUTINE "TWOB"
C   04-MAY-1988 NVE CERN GENEVA
C               CORRECT MASSES INTRODUCED FOR THE VARIOUS PROCESSES
C               IN THE ABSORBTION AND ANNIHILATION ROUTINES
C   06-MAY-1988 NVE CERN GENEVA
C               EXTRA ARGUMENT "NFL" ADDED IN ROUTINE "CORANH"
C               TO INDICATE TARGET TYPE
C   10-MAY-1988 NVE CERN GENEVA
C               OLD VERSION (GEANT 3.11) OF "GHESIG" RESTORED AND
C               ONLY SLIGHTLY MODIFIED TO INCLUDE GEOMETRICAL X-SECT.
C               OF LIGHT MASS FRAGMENTS
C               THE NEW VERSION OF GHESIG IS CONTAINED IN DECK "GHESIX"
C               ON THIS PAM IN ORDER TO BE INVESTIGATED
C   11-MAY-1988 NVE CERN GENEVA
C               K0/K0 BAR CHOICE REMOVED FROM "GHESIG" TO AVOID
C               UNPHYSICAL RESULTS
C               THE K0/K0 BAR CHOICE IS NOW MADE RIGHT AFTER THE
C               CREATION OF THE K0 IN ROUTINE "GHEISH"
C   11-MAY-1988 NVE CERN GENEVA
C               NUCLEAR MASSES OF D, T AND ALPHA RE-COMPUTED
C               TAKING INTO ACCOUNT BINDING ENERGIES AND NEW VALUES
C               USED IN "GHEINI"
C               THESE NEW VALUES WILL ALSO BE USED IN GEANT AND THEY
C               AGREE WITH PARTICLE DATA 1986
C   11-MAY-1988 NVE CERN GENEVA
C               FLAG "NUCFLG" INTRODUCED IN ROUTINES "CASFRG" AND
C               "GNSLWD" TO INDICATE ACTION TAKEN BY "NUCREC" CALL
C   11-MAY-1988 NVE CERN GENEVA
C               CALL TO "COSCAT" ADDED IN ROUTINE "CASFRG" IN CASE
C               THE "NUCREC" CALL DID NOT GIVE ANY ACTION
C   18-MAY-1988 NVE CERN GENEVA
C               DOUBLE COUNTING OF EDEPS FOR STOPPING PARTICLES
C               FIXED (GHEISH AND GHSTOP) AND NEW STEERING FOR
C               LOW ENERGY NUCLEONS IN "GHEISH"
C   18-MAY-1988 NVE CERN GENEVA
C               BUGS FIXED FOR LOW MOMENTUM AND CALL TO "ATOMAS"
C               REPLACED BY "RMASS" ASSIGNMENTS FOR P, N, D, T AND
C               ALPHA IN ROUTINE "NUCREC"
C   19-MAY-1988 NVE CERN GENEVA
C               MORE ACCURATE VALUES OF MASSES OF D AND ALPHA USED
C               IN ROUTINE "ATOMAS"
C   20-MAY-1988 NVE CERN GENEVA
C               PARTICLE MASSES, CHARGES AND CT VALUES FOR GHEISHA
C               DIRECTLY TAKEN FROM THE GEANT ZEBRA BANKS IN
C               ROUTINE "GHEINI"
C   25-MAY-1988 NVE CERN GENEVA
C               COMMON "/CSDATC" REMOVED AND SEQUENCES "/CSDIM"
C               AND "/CSDAT" INTRODUCED FOR CROSS SECTION DATA IN
C               ORDER TO SAVE MEMORY SPACE
C               THESE SEQUENCES ARE INCLUDED IN ROUTINES "GHESIG"
C               AND "PCSDAT"
C   15-JUN-1988 NVE CERN GENEVA
C               SEQUENCES "PCODIM" AND "PCODAT" INTRODUCED FOR
C               CENTRALISED GEANT/GHEISHA PARTICLE CODE CONVERSION
C               THESE SEQUENCES ARE INCLUDED IN ROUTINES "GHESIG"
C               AND "GHEISH"
C   15-JUN-1988 NVE CERN GENEVA
C               FOLLOWING UNUSED SEQUENCES REMOVED
C               /BEAM  /CALSTR /CBRDK /CMUEL  /CPADE /CQINT  /EQUAT
C               /INPEV /MALARR /MUPAC /PARCOM /SCINT /STEPSI
C   15-JUN-1988 NVE CERN GENEVA
C               ARRAY "IBAR(33)" REMOVED FROM ROUTINE "GHEISH"
C   23-JUN-1988 R.BRUN/NVE CERN GENEVA
C               STMT. "KKPART=KPART" ADDED AT THE BEGINNING OF
C               ROUTINE "GHEISH" TO PREVENT INFINITE LOOP FOR
C               PARTICLE OF ZERO MOMENTUM
C   23-JUN-1988 NVE CERN GENEVA
C               ALL VARIABLES "THENVE" RENAMED TO "RTHNVE" TO MATCH
C               NEW IBM COMPILER (= ANSI STANDARD ?)
C   29-JUN-1988 NVE CERN GENEVA
C               NAME OF ROUTINE "LENGTH" AND RELATED CALLS CHANGED
C               TO "LENGTX" TO AVOID CLASH WITH INTERNAL ROUTINE
C               ON THE CRAY COMPUTER
C   12-JUL-1988 NVE CERN GENEVA
C               INELASTICITY FOR HADRONIC INTERACTIONS RE-CALCULATED
C               AND CHANGES MADE IN ROUTINE "CINEMA" ACCORDINGLY
C   13-JUL-1988 NVE CERN GENEVA
C               IPV(...) STUFF REMOVED AND CALL TO "RNDM" CHANGED
C               TO "RANRES" IN ROUTINE "COSCAT"
C               ALSO PROTECTION FOR "SINT" INCLUDED IN "COSCAT"
C   13-JUL-1988 NVE CERN GENEVA
C               FUNCTION "RANRES" ADDED TO PATCH "GHEISHA"
C   01-AUG-1988 NVE CERN GENEVA
C               BUG FIXED IN PARTICLE DISTRIBUTION IN BACKWARD
C               AND FORWARD HEMISPHERE IN ROUTINE "TWOCLU"
C               THIS ALSO AVOIDS INFINITE LOOPING FOR HIGH MULT.
C   15-AUG-1988 NVE CERN GENEVA
C               SETTING OF MECHANISM "STOP" SUPRESSED IN "GHSTOP"
C               TO AVOID DOUBLE SETTING.
C   16-JAN-1989 NVE CERN GENEVA
C               CONDITIONAL CALL TO "GHEINI" ADDED IN ROUTINE "GHEISH"
C               TO ENSURE PROPER INITIALIZATION IN ALL CASES.
C   19-JAN-1989 NVE CERN GENEVA
C               INITIALIZATION OF IPA(...) ARRAY ADDED IN ROUTINES
C               "GHEISH" AND "GNSLWD" BY "CALL VZERO(IPA(1),100)"
C               TO MATCH THE ROUTINE "GENXPT"
C   31-JAN-1989 NVE CERN GENEVA
C               PRINT STATEMENT IN "EXNU" MADE MORE CLEAR
C   01-FEB-1989 NVE CERN GENEVA
C               PARTICLE CODES ARRAY "IKPART" EXTENDED
C               TO SIZE (35) IN SEQUENCE "PCODIM"
C               AND OMEGA AND ANTI-OMEGA ADDED TO GHEISHA CODES
C               AT POS. 33 AND 34 RESP. IN SEQUENCE "PCODAT"
C               ALSO ARRAYS "RMASS" AND "RCHARG" EXTENDED IN SIZE
C               IN "/CONSTS"
C               OMEGA AND ANTI-OMEGA PARAMETERS ALSO ADDED IN "GHEINI"
C   01-FEB-1989 NVE CERN GENEVA
C               ROUTINES "GHESIG", "GHESIX" AND "GHEISH"
C               UPDATED FOR OMEGA AND ANTI-OMEGA PARTICLES
C               ALSO SIZE OF "CSEL" AND "CSIN" ARRAYS EXTENDED
C               IN "/CSDIM" AND "/CSDAT" UPDATED ACCORDINGLY BY
C               CHANGE IN ROUTINE "PCSDAT" OF PATCH "GHUTILS"
C   02-FEB-1989 NVE CERN GENEVA
C               DECKS "CASXM" "CASX0" "CASOM" "CASAXM" "CASAX0"
C               AND "CASAOM" ADDED TO TREAT THE CASCADES OF
C               XI-, XI0, OMEGA- AND THEIR ANTI-PARTICLES.
C               THE CASCADE ROUTINES FOR THE ANTI-PARTICLES ARE
C               FIXUPS FOR THE MOMENT, BUT CORRECT ROUTINES ARE
C               UNDER DEVELOPMENT.
C   16-FEB-1989 RB CERN GENEVA
C               ALL "RNDM" STUFF REPLACED BY NEW GEANT RANDOM NUMBER
C               GENERATOR "GRNDM"
C   13-MAR-1989 NVE CERN GENEVA
C               ALL CHANGES MADE BY RENE BRUN CHECKED AND SOME
C               SERIOUS MISTAKES FIXED
C   17-MAR-1989 NVE CERN GENEVA
C               CALCULATIONS FOR HADRONIC CROSS SECTIONS RE-DONE
C               AND NEW FUNCTION "GHESIG" ADDED IN PATCH "GIFACE"
C               THE OLD VERSION IS STILL AVAILABLE IN DECK "GHESIX"
C               FOR BACKWARD COMPATIBILITY TESTS
C   20-MAR-1989 NVE CERN GENEVA
C               CORRECT INTERPOLATIONS AND EXTRAPOLATIONS IN
C               ROUTINES "STEEQ" AND "STPAIR"
C               ALSO ALL "ANPOL" AND "ANPOL2" STUFF REMOVED
C   31-MAR-1989 NVE CERN GENEVA
C               PARTICLE MASSES FOR PACKAGES OTHER THAN "GEANT" UPDATED
C               ACCORDING TO PART. DATA 1986 IN ROUTINE "GHEINI"
C   23-MAY-1989 NVE CERN GENEVA
C               "DIMENSION PR(50)" CHANGED TO "DIMENSION PR(*)"
C               IN ROUTINE "ROTES2"
C   29-MAY-1989 R.Brun
C               Protection of all SQRT against <0
C               GRNDM changed to 2 arguments only
C   12-JUN-1989 NVE CERN GENEVA
C               "COST" PROTECTION RESTORED IN "GMUNU"
C   13-JUN-1989 H.FESEFELDT
C               Corrections in GHESIG/GHESIX for ACC40
C               R.B, FCA,F.B
C               Protection in GHESIG for IPART>48
C   29-JUN-1989 RB/NVE CERN GENEVA
C               MASSES UPDATED IN GHEINI ACCORDING TO PART. DATA APR-1988
C               AND PROTECTION AGAINST PRECISION PROBLEMS OF EK AND EDEP
C               IN GHEISH
C   29-JUN-1989 NVE CERN GENEVA
C               MASSES TAKEN FROM GEANT PARTICLE DATA BANKS IN
C               GMUNU, GMUNUI AND GMUSIG
C               AND MASSES TAKEN FROM "RMASS" ARRAY IN ATOMAS
C   29-JUN-1989 RB CERN GENEVA
C               "PDK" AND "DPDK" RENAMED TO "GPDK" AND "GDPDK" RESP.
C               BECAUSE OF CONFLICT WITH NEW CERNLIB ROUTINES
C               AND CALLS IN "PHASP" UPDATED ACCORDINGLY
C   29-JUN-1989 NVE CERN GENEVA
C               OBSOLETE ROUTINE "GHESIX" REMOVED FROM PATCH "GIFACE"
C   04-JUL-1989 NVE CERN GENEVA
C               PROTECTION AGAINST ANNIHILATION PROCESSES INTRODUCED
C               IN "STPAIR"
C   18-AUG-1989 NVE CERN GENEVA
C               PROTECTION FOR NGKINE .GT. 100 ADDED IN "GHEISH"
C               AND ALSO PROTECTIONS ADDED IN "SETTRK"
C   14-SEP-1989 HF RWTH AACHEN
C               CORRECTED IPARM SETTING IN "CINEMA"
C   05-MAR-1990 HF RWTH AACHEN
C               DOUBLE PRECISION RESTORED FOR PHASE SPACE PACKAGE
C               AND PHOTO-FISSION AND ABSORBTION PROCESSES MODIFIED
C   05-MAR-1990 RB CERN GENEVA
C               SOME PROTECTIONS INTRODUCED IN MODIFIED CODE
C   11-MAY-1990 NVE CERN GENEVA
C               MODIFICATIONS OF HF AND RB IMPLEMENTED BY RETAINING
C               THE LOGIC OF GEANH 3.13
C   14-MAY-1990 NVE CERN GENEVA
C               PROTECTION FOR ATAN2 INTRODUCED IN "GMUNU"
C   15-MAY-1990 NVE CERN GENEVA
C               COMPLETE "+KEEP" FOR /NUCIO/ UNDER PATCHY FLAG "SINGLE"
C               AND THIS SEQUENCE INTODUCED IN "NUCREC" AND "PHPNUC"
C   17-MAY-1990 HF RWTH AACHEN/NVE CERN GENEVA
C               DOUBLE INVOKATION OF NUCREC REMOVED FROM "GHEISH"
C               FOR SLOW NEUTRONS WHICH RESULTS IN LONGER PATH
C   18-MAY-1990 RB CERN GENEVA
C               CLEANER CONSTRUCTION TO PUT DOUBLE PRECISION OF
C               /NUCIO UNDER CONTROL OF PATCHY FLAG "SINGLE"
C   20-JUN-1990 RB/NVE CERN GENEVA
C               NEW OPTION IHADR=3 INTRODUCED TO SELECT INVOKATION
C               OF THE NUCRIN PACKAGE FROM "GHEISH" FOR ENERGIES
C               BELOW 5 GEV
C   27-SEP-1990 RB/NVE CERN GENEVA
C               NUCRIN PACKAGE REPLACED BY NEW VERSION OF RB
C               AND PROTECTION AGAINST CRAZY P VALUES ADDED IN "GHEISH"
C   27-SEP-1990 MR/NVE CERN GENEVA
C               CORRECTION FOR THETA=PI IN ROUTINE "DEFS1"
C               AS A RESULT OF THE BUG THE ORIGINAL (HIGH ENERGETIC)
C               PARTICLE MIGHT HAVE GOT LOST IN CASE OF AN INTERACTION
C               AT THE BEGINNING OF A DETECTOR
C   27-SEP-1990 HJT/NVE CERN GENEVA
C               AVOIDING NEGATIVE W2 IN "GMUNU" AND CORRECTIONS MADE FOR
C               'SIG' IN "GMUNUI"
C   27-SEP-1990 NVE CERN GENEVA
C               PRINT MESSAGES WITH CONTACT PERSON REMOVED IN "GHEISH"
C   28-SEP-1990 NVE CERN GENEVA
C               GIFACE AND GHEISH CODE CLEANED FROM UNUSED VARIABLES ETC.
C               AS PROPOSED BY RB
C   25-OCT-1990 NVE CERN GENEVA
C               PROTECTION FOR 'ZINTHA' TO BECOME -0. ADDED IN "GHEISH"
C               TO PREVENT INFINITE LOOP DUE TO PRECISION PROBLEM
C   26-OCT-1990 NVE CERN GENEVA
C               LOCAL ARRAYS INITIALISED IN "PHASP" AND "PHPNUC" TO
C               ENABLE RUNNING IN DYNAMIC COMP. MODE UNDER UNIX
C   29-OCT-1990 NVE CERN GENEVA
C               BUG FIXED FOR DOUBLE PRECISION RANDOM NUMBERS IN "PHPNUC"
C               AND FORMAT FOR DEBUG PRINT CHANGED IN "NUCREC"
C
C ********************************************************************
C ***                                                              ***
C ***  IN CASE OF ANY TROUBLES, COMMENTS OR SUGGESTIONS CONTACT :  ***
C ***                                                              ***
C ***                     NICK VAN EIJNDHOVEN                      ***
C ***                     CERN                                     ***
C ***                     B 595 R-004                              ***
C ***                     CH-1211 GENEVA 23                        ***
C ***                     SWITZERLAND                              ***
C ***                     EUROPE                                   ***
C ***                     TEL. +41-22-7674740 OR CERN BEEP 7294    ***
C ***                     BITNET ADDRESS : NICK@CERNVM.CERN.CH     ***
C ***                                                              ***
C ********************************************************************
