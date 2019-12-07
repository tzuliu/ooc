!
!  RCMD SHLIB legislator_cuttingplane_routine.f90 -lRblas -lRlapack
!
! *********************************************************************
!   SUBROUTINE KPLEGIS -- PERFORMS THE LEGISLATIVE PROCEDURE
!                           30 JUNE 1999
!
! *********************************************************************
!
      SUBROUTINE KPLEGIS(NUMMEMBERS,NUMVOTES,&
                         JJJ,NP,NRCALL,NS,NDUAL,XMAT,LLEGERR, &
                        ZVEC,WS,MCUTS,LERROR,LTOTAL,MWRONG, &
                        LDATA,IPRINT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION XMAT(NUMMEMBERS,25),ZVEC(NUMVOTES,25),&
              MCUTS(NUMVOTES,2),WS(NDUAL),LERROR(NUMMEMBERS,NUMVOTES),&
              LLEGERR(NUMMEMBERS,2),LDATA(NUMMEMBERS,NUMVOTES)
!
      DOUBLE PRECISION, ALLOCATABLE :: BB(:,:)
      ALLOCATE(BB(25,NRCALL))
!
!  885 FORMAT(' ERROR CHECK ON ENTRY TO KPLEGIS',I5,I7)
! 1094 FORMAT(' LEG CLASSIFICATION ERROR  ',2I3,2I8,2F10.5)
!
      IF(JJJ.EQ.1) JJJ=1  !hack to get rid of warnings
      CALL ZVECINV(NUMMEMBERS,NUMVOTES,NRCALL,NS,ZVEC,BB,IPRINT)
!
      LTOTAL=0
      LWRONG=0
      MWRONG=0
      DO 111 I=1,NP
      KCHECK4=0
      DO 882 JX=1,NRCALL
      KCHECK4=KCHECK4+LERROR(I,JX)
  882 CONTINUE
      XSAVE1=XMAT(I,1)
      XSAVE2=XMAT(I,2)
      IVOT=I
      CALL KTPXI(NUMMEMBERS,NUMVOTES,&
                 IVOT,NP,NRCALL,NS,NDUAL,MCUTS,BB,XMAT,ZVEC,WS, &
                  LERROR,KTOTAL,KWRONG,LDATA)
!
!      WRITE(21,1492)I,KCHECK4,KWRONG,KTOTAL,XSAVE1,XSAVE2,
!     C                    XMAT(I,1),XMAT(I,2)
! 1492 FORMAT(4I8,4F7.3)
      MWRONG=MWRONG+KWRONG
      LTOTAL=LTOTAL+KTOTAL
      LLEGERR(I,1)=KWRONG
      LLEGERR(I,2)=KTOTAL
!
  111 CONTINUE
      IF(LTOTAL.GT.0)THEN
        XERROR=FLOAT(MWRONG)/FLOAT(LTOTAL)
        YERROR=1.0-XERROR
      ENDIF
!      IF(IPRINT.EQ.1)WRITE(23,1094)JJJ,NS,MWRONG,LTOTAL,XERROR,YERROR
      DEALLOCATE(BB)
      RETURN
      END
!
! **************************************************************************
!  SUBROUTINE ZVECINV -- CALCULATES (Z'Z)-1Z' WHERE Z IS THE NORMAL PLANE
!                        MATRIX
! **************************************************************************
!
      SUBROUTINE ZVECINV(NUMMEMBERS,NUMVOTES,NRCALL,NS,ZVEC,BB,IPRINT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION BB(25,NRCALL),ZVEC(NUMVOTES,25)
      INTEGER, ALLOCATABLE :: IWORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: VVV(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: VVV2(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: FV1(:)
      DOUBLE PRECISION, ALLOCATABLE :: FV2(:)
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: VCOV(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: VCOV2(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: WVEC(:)
      DOUBLE PRECISION, ALLOCATABLE :: UL(:,:)
      ALLOCATE(VVV(25,25))
      ALLOCATE(VVV2(25,25))
      ALLOCATE(FV1(NRCALL))
      ALLOCATE(FV2(NRCALL))
      ALLOCATE(VCOV(25,25))
      ALLOCATE(VCOV2(25,25))
      ALLOCATE(WVEC(25))
      ALLOCATE(UL(25,25))
      ALLOCATE(IWORK(8*25*25+1875))
      ALLOCATE(WORK(8*25*25+1875))
!
      LWORK=8*25*25+1875
      KPDUDE=NUMMEMBERS
!
! 1011 FORMAT(7F10.4)
! 1012 FORMAT(' DECOMPOSITION OF NORMAL VECTOR MATRIX',3I4)
! 1091 FORMAT(' INVERSE MATRIX ERROR',F10.4)
!
!
!    (X'X)
!
      IF(IPRINT.EQ.1) IPRINT=1  !hack to get rid of warnings
      DO 38 K=1,NS
      DO 39 L=1,NS
      SUM=0.0
      DO 40 I=1,NRCALL
      SUM=SUM+ZVEC(I,K)*ZVEC(I,L)
  40  CONTINUE
      VCOV(K,L)=SUM
      VCOV2(K,L)=SUM
  39  CONTINUE
  38  CONTINUE
!
!  EIGENVECTOR-EIGENVALUE DECOMPOSITION OF NORMAL VECTOR MATRIX
!
!
      CALL DGESDD('S',NS,NS,VCOV2,25,WVEC,VVV, &
                 25,VVV2,25,WORK,LWORK,IWORK,IRANK)
!
!      CALL KPRS(25,NS,VCOV,WVEC,1,VVV,FV1,FV2,IER)
!
!  (X'X)-1
!
      DO 83 I=1,NS
      DO 83 K=1,NS
      SUM=0.0
      DO 84 J=1,NS
!      IF(ABS(WVEC(NS+1-J)).GT..0001)THEN
!          SUM=SUM+VVV(K,NS+1-J)*(1.0/WVEC(NS+1-J))*VVV(I,NS+1-J)
      IF(ABS(WVEC(J)).GT..0001)THEN
          SUM=SUM+VVV(K,J)*(1.0/WVEC(J))*VVV(I,J)
      ENDIF
  84  CONTINUE
  83  UL(I,K)=SUM
!
!
!  MATRIX INVERSION CHECK  (X'X)-1(X'X) = I
!
      ASUM=0.0
      DO 933 I=1,NS
      DO 933 J=1,NS
      SUM=0.0
      DO 944 K=1,NS
      SUM=SUM+UL(J,K)*VCOV(K,I)
  944 CONTINUE
      IF(I.EQ.J)ASUM=ASUM+ABS(1.0-SUM)
      IF(I.NE.J)ASUM=ASUM+ABS(SUM)
  933 CONTINUE
!      IF(ASUM.GT..01.AND.IPRINT.EQ.1)WRITE(23,1091)ASUM
!
!  (X'X)-1*X'
!
      DO 85 I=1,NRCALL
      DO 85 J=1,NS
      SUM=0.0
      DO 86 JJ=1,NS
      SUM=SUM+UL(J,JJ)*ZVEC(I,JJ)
  86  CONTINUE
  85  BB(J,I)=SUM
!
      DEALLOCATE(VVV)
      DEALLOCATE(VVV2)
      DEALLOCATE(FV1)
      DEALLOCATE(FV2)
      DEALLOCATE(VCOV)
      DEALLOCATE(VCOV2)
      DEALLOCATE(WVEC)
      DEALLOCATE(UL)
      DEALLOCATE(WORK)
      DEALLOCATE(IWORK)
      RETURN
      END
!
! **************************************************************************
!  SUBROUTINE KTPXI -- PERFORMS LEGISLATOR FITTING
! **************************************************************************
!
      SUBROUTINE KTPXI(NUMMEMBERS,NUMVOTES,&
                       ILEG,NP,NRCALL,NS,NDUAL,MCUTS,BB,XMAT,ZVEC,WS, &
                              LERROR,KTOTAL,KWRONG,LDATA)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION BB(25,NRCALL),XMAT(NUMMEMBERS,25),WS(NDUAL),&
                LERROR(NUMMEMBERS,NUMVOTES), &
                ZVEC(NUMVOTES,25),MCUTS(NUMVOTES,2),&
                LDATA(NUMMEMBERS,NUMVOTES)
!
      INTEGER, ALLOCATABLE :: LL(:)
      INTEGER, ALLOCATABLE :: LSAVE(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: ZWRONG(:)
      DOUBLE PRECISION, ALLOCATABLE :: XDAT(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: XXZ(:)
      DOUBLE PRECISION, ALLOCATABLE :: YWRONG(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: YYWRONG(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: XMAT2(:)
      DOUBLE PRECISION, ALLOCATABLE :: XXY(:)
      DOUBLE PRECISION, ALLOCATABLE :: YYYWRONG(:,:)
      ALLOCATE(LL(NRCALL))
      ALLOCATE(LSAVE(NRCALL,25))
      ALLOCATE(ZWRONG(50))
      ALLOCATE(XDAT(NP,25))
      ALLOCATE(XXZ(NRCALL))
      ALLOCATE(YWRONG(NRCALL,25))
      ALLOCATE(YYWRONG(NRCALL,25))
      ALLOCATE(XMAT2(25))
      ALLOCATE(XXY(NRCALL))
      ALLOCATE(YYYWRONG(NRCALL,25))
!
! 1001 FORMAT(I4,2X,I6,I4)
! 1011 FORMAT(I4,2I3,I4,12F7.3)
! 1012 FORMAT(14X,12F7.3)
! 1013 FORMAT(' LEG',2I3,I4,12F7.3)
! 1014 FORMAT(' WHOA DUDE',I4,I3,I4)
! 1015 FORMAT(' HEY THESE ARE NOT EQUAL!')
! 1016 FORMAT(I4,2I3,I4,12F7.3)
! 1017 FORMAT(14X,I4,12F7.3)
! 1018 FORMAT(10X,2I4,12F7.3)
! 1019 FORMAT(16X,12F7.3)
! 1020 FORMAT(2I4,3F10.4)
! 1021 FORMAT(4I4,12F7.3)
! 1022 FORMAT(' TWO  ',3I4,12F7.3)
! 1023 FORMAT(' THREE',3I4,12F7.3)
! 1024 FORMAT(' FOUR ',3I4,12F7.3)
! 1025 FORMAT(4X,4X,12X,12F7.3)
! 1026 FORMAT(I4,25F7.3)
!
      DO 8999 I=1,50
      ZWRONG(I)=0.0
 8999 CONTINUE
      NTRIAL2=15
      MZZZ=5
!
!      WRITE(21,1019)(XMAT(ILEG,K),K=1,NS)
      DO 203 K=1,NS
      XDAT(20,K)=XMAT(ILEG,K)
  203 CONTINUE
!
      DO 9887 IIII=1,2
      SUM=0.0
      DO 101 K=1,NS
      XDAT(1,K)=0.0
!      IF(IIII.EQ.2)XDAT(1,K)=0.1
      IF(IIII.EQ.2)XDAT(1,K)=XDAT(20,K)
      XDAT(2,K)=XDAT(1,K)
  101 CONTINUE
!
!      WRITE(21,1019)(XDAT(1,K),K=1,NS)
      KWRONG=0
      K3WRONG=0
!
      DO 9888 III=1,MZZZ
!
      DO 988 II=1,NTRIAL2
      NII=II
      DO 987 KL=1,NS
!
      K3WRONG=KWRONG
      XDAT(2,KL)=0.01
!
!
!  CALCULATE FEASIBLE ALPHA VALUE
!
      ASUM=0.0
      AAA=0.0
      BBB=0.0
      DO 902 ML=1,NS
      ASUM=ASUM+XDAT(1,ML)**2
      AAA=AAA+(XDAT(2,ML)-XDAT(1,ML))**2
      BBB=BBB+2.0*(XDAT(2,ML)-XDAT(1,ML))*XDAT(1,ML)
  902 CONTINUE
      CCC=ASUM - 1.0
      RAD=BBB*BBB-4.0*AAA*CCC
      RADSQ=SQRT(ABS(RAD))
      ROOT1=0.0
      ROOT2=0.0
      IF(AAA.GT..00001)THEN
         ROOT1=(-BBB+RADSQ)/(2.0*AAA)
         ROOT2=(-BBB-RADSQ)/(2.0*AAA)
      ENDIF
!
      CALL KTPXIXJ(NUMMEMBERS,NUMVOTES,&
                   NII,ILEG,NP,NRCALL,NS,NDUAL,MCUTS,BB,XDAT,ZVEC,WS, &
                      XXZ,WSSY,XMAT2,ZWRONG,YYWRONG,LERROR,KTOTAL, &
                      KWRONG,ROOT1,ROOT2,LDATA)
!
      LSAVE(II,KL)=KWRONG
      IF(KWRONG.EQ.0)THEN
         LL(1)=II
         DO 668 K=1,NS
         XMAT(ILEG,K)=XMAT2(K)
  668    CONTINUE
         KKNII=NII
         K3WRONG=KWRONG
!         WRITE(21,1021)ILEG,III,II,KWRONG,(XMAT(ILEG,K),K=1,NS)
         GO TO 996
      ENDIF
!
      SUM=0.0
      DO 133 K=1,NS
      SUM=SUM+XMAT2(K)**2
      YWRONG(II,K)=XMAT2(K)
      YYYWRONG(KL,K)=XMAT2(K)
      XDAT(3,K)=XDAT(1,K)
      XDAT(1,K)=XMAT2(K)
      XDAT(2,K)=XMAT2(K)
      XMAT(ILEG,K)=XMAT2(K)
  133 CONTINUE
!
!      WRITE(21,1011)ILEG,II,KL,KWRONG,(XDAT(3,K),K=1,NS)
!      WRITE(21,1012)(XMAT2(K),K=1,NS),(ZWRONG(K),K=1,NS)
!
  987 CONTINUE
!
      SUM=0.0
      KSUM=0
      DO 132 K=1,NS
      IF(II.GT.1)THEN
         SUM=SUM+(YWRONG(II,K)-YWRONG(II-1,K))**2
         KSUM=KSUM+ABS(LSAVE(II,K)-LSAVE(II-1,K))
      ENDIF
  132 CONTINUE
      IF(II.EQ.1)GO TO 988
      IF(SUM.LE..000001)GO TO 986
      IF(II.GT.5.AND.KSUM.EQ.0)GO TO 986
  988 CONTINUE
!
  986 CONTINUE
!      WRITE(21,1021)ILEG,III,II,KWRONG,(XMAT(ILEG,K),K=1,NS)
      KKKNII=NII
      K3WRONG=KWRONG
!
!
!  CALCULATE FEASIBLE ALPHA VALUE
!
      ASUM=0.0
      AAA=0.0
      BBB=0.0
      DO 903 ML=1,NS
      XDAT(1,ML)=XMAT(ILEG,ML)
      XDAT(2,ML)=ZWRONG(K)
      ASUM=ASUM+XDAT(1,ML)**2
      AAA=AAA+(XDAT(2,ML)-XDAT(1,ML))**2
      BBB=BBB+2.0*(XDAT(2,ML)-XDAT(1,ML))*XDAT(1,ML)
  903 CONTINUE
      CCC=ASUM - 1.0
      RAD=BBB*BBB-4.0*AAA*CCC
      RADSQ=SQRT(ABS(RAD))
      ROOT1=0.0
      ROOT2=0.0
      IF(AAA.GT..00001)THEN
         ROOT1=(-BBB+RADSQ)/(2.0*AAA)
         ROOT2=(-BBB-RADSQ)/(2.0*AAA)
      ENDIF
!      ROOT1=(-BBB+RADSQ)/(2.0*AAA)
!      ROOT2=(-BBB-RADSQ)/(2.0*AAA)
!
      CALL KTPXIXJ(NUMMEMBERS,NUMVOTES,&
                   NII,ILEG,NP,NRCALL,NS,NDUAL,MCUTS,BB,XDAT,ZVEC,WS, &
                      XXZ,WSSY,XMAT2,ZWRONG,YYWRONG,LERROR,KTOTAL, &
                      KWRONG,ROOT1,ROOT2,LDATA)
!
!      WRITE(21,1021)ILEG,III,II,KWRONG,(XMAT2(K),K=1,NS),WSSY
!
      DO 44 K=1,NS
      XMAT(ILEG,K)=XMAT2(K)
      XDAT(1,K)=XMAT2(K)
      XDAT(2,K)=XMAT2(K)
      IF(IIII.EQ.1)XDAT(10,K)=XMAT2(K)
  44  CONTINUE
      IF(KWRONG.EQ.0)GO TO 996
 9888 CONTINUE
      IF(IIII.EQ.1)KLWRONG=KWRONG
      IF(IIII.EQ.2)THEN
!
!  CALCULATE FEASIBLE ALPHA VALUE
!
         ASUM=0.0
         AAA=0.0
         BBB=0.0
         DO 905 ML=1,NS
         XDAT(1,ML)=XMAT(ILEG,ML)
         XDAT(2,ML)=XDAT(10,ML)
         ASUM=ASUM+XDAT(1,ML)**2
         AAA=AAA+(XDAT(2,ML)-XDAT(1,ML))**2
         BBB=BBB+2.0*(XDAT(2,ML)-XDAT(1,ML))*XDAT(1,ML)
  905    CONTINUE
         CCC=ASUM - 1.0
         RAD=BBB*BBB-4.0*AAA*CCC
         RADSQ=SQRT(ABS(RAD))
         ROOT1=0.0
         ROOT2=0.0
         IF(AAA.GT..00001)THEN
            ROOT1=(-BBB+RADSQ)/(2.0*AAA)
            ROOT2=(-BBB-RADSQ)/(2.0*AAA)
         ENDIF
!         ROOT1=(-BBB+RADSQ)/(2.0*AAA)
!         ROOT2=(-BBB-RADSQ)/(2.0*AAA)
!
         CALL KTPXIXJ(NUMMEMBERS,NUMVOTES,&
                      NII,ILEG,NP,NRCALL,NS,NDUAL,MCUTS,BB,XDAT,ZVEC, &
                      WS,XXZ,WSSY,XMAT2,ZWRONG,YYWRONG,LERROR,KTOTAL, &
                      KWRONG,ROOT1,ROOT2,LDATA)
!
!         WRITE(21,1021)ILEG,III,II,KWRONG,(XMAT2(K),K=1,NS),WSSY
!
         DO 48 K=1,NS
         XMAT(ILEG,K)=XMAT2(K)
  48     CONTINUE
         IF(KWRONG.EQ.0)GO TO 996
      ENDIF
!
 9887 CONTINUE
!
  996 CONTINUE
!
      SUM=0.0
      DO 12 K=1,NS
      SUM=SUM+XMAT(ILEG,K)**2
  12  CONTINUE
      IF(SUM.GT.1.0)THEN
         DO 13 K=1,NS
         XMAT(ILEG,K)=XMAT(ILEG,K)/SQRT(SUM)
  13     CONTINUE
      ENDIF
!
      KKRITE=0
      KKWRONG=0
      KTOTAL=0
      DO 31 JX=1,NRCALL
      LERROR(ILEG,JX)=0
      SUM=0.0
      DO 32 K=1,NS
      SUM=SUM+XMAT(ILEG,K)*ZVEC(JX,K)
  32  CONTINUE
      XXY(JX)=SUM
      DB2B1=WS(JX)-XXY(JX)
!
!  CALCULATE CLASSIFICATION ERROR
!
      IF(LDATA(ILEG,JX).NE.0)THEN
         KTOTAL=KTOTAL+1
         IF(XXY(JX).LT.WS(JX))THEN
!            KTOTAL=KTOTAL+1
            IF(LDATA(ILEG,JX).EQ.MCUTS(JX,1))THEN
               KKRITE=KKRITE+1
            ENDIF
            IF(LDATA(ILEG,JX).EQ.MCUTS(JX,2))THEN
               LERROR(ILEG,JX)=1
               KKWRONG=KKWRONG+1
            ENDIF
         ENDIF
         IF(XXY(JX).GT.WS(JX))THEN
!            KTOTAL=KTOTAL+1
            IF(LDATA(ILEG,JX).EQ.MCUTS(JX,2))THEN
               KKRITE=KKRITE+1
            ENDIF
            IF(LDATA(ILEG,JX).EQ.MCUTS(JX,1))THEN
               LERROR(ILEG,JX)=1
               KKWRONG=KKWRONG+1
            ENDIF
         ENDIF
      ENDIF
  31  CONTINUE
!      WRITE(21,1001)ILEG,KKRITE,KKWRONG
      IF(KKWRONG.NE.KWRONG)THEN
!         WRITE(21,1015)
      ENDIF
      SUM=0.0
      DO 204 K=1,NS
      SUM=SUM+(XDAT(20,K)-XMAT(ILEG,K))**2
  204 CONTINUE
      SUM=SQRT(SUM)
!      WRITE(23,1026)ILEG,SUM
      DEALLOCATE(LL)
      DEALLOCATE(LSAVE)
      DEALLOCATE(ZWRONG)
      DEALLOCATE(XDAT)
      DEALLOCATE(XXZ)
      DEALLOCATE(YWRONG)
      DEALLOCATE(YYWRONG)
      DEALLOCATE(XMAT2)
      DEALLOCATE(XXY)
      DEALLOCATE(YYYWRONG)
      RETURN
      END
!
! **************************************************************************
!  SUBROUTINE KTPXIXJ -- PERFORMS LEGISLATOR FITTING
! **************************************************************************
!
      SUBROUTINE KTPXIXJ(NUMMEMBERS,NUMVOTES,&
                         NII,ILEG,NP,NRCALL,NS,NDUAL,MCUTS,BB,XDAT,ZVEC, &
                         WS,XXZ,WSSY,XMAT2,ZWRONG,YWRONG,LERROR,KTOTAL, &
                            KKWRONG,ROOT1,ROOT2,LDATA)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION BB(25,NRCALL),XDAT(NP,25),WS(NDUAL),&
                LERROR(NUMMEMBERS,NUMVOTES), &
                ZVEC(NUMVOTES,25),MCUTS(NUMVOTES,2),XXZ(NRCALL), &
                ZWRONG(50),YWRONG(NRCALL,25),XMAT2(25),&
                LDATA(NUMMEMBERS,NUMVOTES)
!
      INTEGER, ALLOCATABLE :: MRITE(:,:)
      INTEGER, ALLOCATABLE :: LALL(:)
      DOUBLE PRECISION, ALLOCATABLE :: YALL(:)
      DOUBLE PRECISION, ALLOCATABLE :: XALL(:)
      DOUBLE PRECISION, ALLOCATABLE :: KALL(:)
      DOUBLE PRECISION, ALLOCATABLE :: XXY(:)
      DOUBLE PRECISION, ALLOCATABLE :: XSAVE(:,:)
      ALLOCATE(MRITE(NRCALL,100))
      ALLOCATE(LALL(NDUAL))
      ALLOCATE(YALL(NDUAL))
      ALLOCATE(XALL(NDUAL))
      ALLOCATE(KALL(NDUAL))
      ALLOCATE(XXY(NRCALL))
      ALLOCATE(XSAVE(NRCALL,100))
!
! 1012 FORMAT(I4,3I6,7F7.3)
! 1013 FORMAT(2I4,3I6,F12.6,6F7.3/8X,8F7.3/8X,8F7.3)
! 1013 FORMAT(2I4,4I6,3F7.3)
! 1014 FORMAT(15X,I6,13X,4F7.3)
! 1015 FORMAT(I5,4F7.3)
! 1016 FORMAT(I5,2F7.3,I3,F7.3)
! 1017 FORMAT(I5,2F10.4,3F7.3)
! 1018 FORMAT(4X,6I5,20F10.4)
! 1019 FORMAT(2I4,4I6)
! 1020 FORMAT(' INITIAL ERROR',I4,3I6,7F7.3)
! 1021 FORMAT(I3,2I5,15F7.3)
! 1023 FORMAT(3I4,F15.10)
!
      KWED1=NII
!
!  INITIALIZE AT ZERO
!
      NTRY=2
!
      DO 99 III=1,NTRY
      KTOTAL=0
      KRITE=0
      KWRONG=0
      DO 11 JX=1,NRCALL
      SUM=0.0
      DO 52 K=1,NS
      SUM=SUM+XDAT(III,K)*ZVEC(JX,K)
  52  CONTINUE
      XXY(JX)=SUM
      XSAVE(JX,III)=SUM
      MRITE(JX,III)=0
      DB2B1=WS(JX)-XXY(JX)
      IF(LDATA(ILEG,JX).NE.0)THEN
         KTOTAL=KTOTAL+1
!
!  CALCULATE CLASSIFICATION ERROR
!
         IF(XXY(JX).LT.WS(JX))THEN
            IF(LDATA(ILEG,JX).EQ.MCUTS(JX,1))THEN
               KRITE=KRITE+1
               MRITE(JX,III)=1
            ENDIF
            IF(LDATA(ILEG,JX).EQ.MCUTS(JX,2))THEN
               KWRONG=KWRONG+1
            ENDIF
         ENDIF
         IF(XXY(JX).GT.WS(JX))THEN
            IF(LDATA(ILEG,JX).EQ.MCUTS(JX,1))THEN
               KWRONG=KWRONG+1
            ENDIF
            IF(LDATA(ILEG,JX).EQ.MCUTS(JX,2))THEN
               KRITE=KRITE+1
               MRITE(JX,III)=1
            ENDIF
         ENDIF
      ENDIF
  11  CONTINUE
!      WRITE(21,1012)III,KRITE,KWRONG,KTOTAL,(XDAT(III,K),K=1,NS)
  99  CONTINUE
!
!
!  CONSTRUCT PROJECTION VECTOR
!
      JJJ=1
      III=2
      ITOT=0
      KRITE=0
      KWRONG=0
      KTOTAL=0
      XERR=0.0
      XKMAX=+99999.0
      XKMIN=-99999.0
      DO 1 JX=1,NRCALL
!
      DB2B1=WS(JX)-XSAVE(JX,III)
      DENOM=XSAVE(JX,III)-XSAVE(JX,JJJ)
      IF(ABS(DENOM).LE.0.00001)THEN
!         WRITE(23,1023)NII,ILEG,JX,DENOM
         GO TO 1
      ENDIF
      XNUM1=+1.0-XSAVE(JX,JJJ)
      XNUM2=-1.0-XSAVE(JX,JJJ)
      XNUM3=WS(JX)-XSAVE(JX,JJJ)
!
!  CALCULATE CLASSIFICATION ERROR
!
      IF(LDATA(ILEG,JX).NE.0)THEN
         IF(XSAVE(JX,III).LT.WS(JX))THEN
            KTOTAL=KTOTAL+1
            IF(LDATA(ILEG,JX).EQ.MCUTS(JX,1))THEN
               KRITE=KRITE+1
!
!  CORRECT TO CORRECT (CASES 3 AND 4)
!
               IF(MRITE(JX,JJJ).EQ.1)THEN
!  CASE 3
                  IF(XSAVE(JX,III).LT.XSAVE(JX,JJJ))THEN
                     ITOT=ITOT+1
                     XALL(ITOT)=XNUM2/DENOM
                     KALL(ITOT)=6
                     ITOT=ITOT+1
                     XALL(ITOT)=XNUM3/DENOM
                     KALL(ITOT)=1
                  ENDIF
!  CASE 4
                  IF(XSAVE(JX,III).GT.XSAVE(JX,JJJ))THEN
                     ITOT=ITOT+1
                     XALL(ITOT)=XNUM3/DENOM
                     KALL(ITOT)=6
                     ITOT=ITOT+1
                     XALL(ITOT)=XNUM2/DENOM
                     KALL(ITOT)=1
                  ENDIF
               ENDIF
!
!  INCORRECT TO CORRECT (CASE 11)
!
               IF(MRITE(JX,JJJ).EQ.0)THEN
                  ITOT=ITOT+1
                  XALL(ITOT)=XNUM2/DENOM
                  KALL(ITOT)=6
                  ITOT=ITOT+1
                  XALL(ITOT)=XNUM3/DENOM
                  KALL(ITOT)=1
               ENDIF
            ENDIF
            IF(LDATA(ILEG,JX).EQ.MCUTS(JX,2))THEN
               KWRONG=KWRONG+1
!
!  CORRECT TO INCORRECT (CASE 9)
!
               IF(MRITE(JX,JJJ).EQ.1)THEN
                  ITOT=ITOT+1
                  XALL(ITOT)=XNUM3/DENOM
                  KALL(ITOT)=6
                  ITOT=ITOT+1
                  XALL(ITOT)=XNUM1/DENOM
                  KALL(ITOT)=1
               ENDIF
!
!  INCORRECT TO INCORRECT (CASES 7 AND 8)
!
               IF(MRITE(JX,JJJ).EQ.0)THEN
!  CASE 7
                  IF(XSAVE(JX,III).LT.XSAVE(JX,JJJ))THEN
                     ITOT=ITOT+1
                     XALL(ITOT)=XNUM3/DENOM
                     KALL(ITOT)=6
                     ITOT=ITOT+1
                     XALL(ITOT)=XNUM1/DENOM
                     KALL(ITOT)=1
                  ENDIF
!  CASE 8
                  IF(XSAVE(JX,III).GT.XSAVE(JX,JJJ))THEN
                     ITOT=ITOT+1
                     XALL(ITOT)=XNUM1/DENOM
                     KALL(ITOT)=6
                     ITOT=ITOT+1
                     XALL(ITOT)=XNUM3/DENOM
                     KALL(ITOT)=1
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
!
!
!
         IF(XSAVE(JX,III).GT.WS(JX))THEN
            KTOTAL=KTOTAL+1
            IF(LDATA(ILEG,JX).EQ.MCUTS(JX,2))THEN
               KRITE=KRITE+1
!
!  CORRECT TO CORRECT (CASES 1 AND 2)
!
               IF(MRITE(JX,JJJ).EQ.1)THEN
!  CASE 2
                  IF(XSAVE(JX,III).LT.XSAVE(JX,JJJ))THEN
                     ITOT=ITOT+1
                     XALL(ITOT)=XNUM3/DENOM
                     KALL(ITOT)=6
                     ITOT=ITOT+1
                     XALL(ITOT)=XNUM1/DENOM
                     KALL(ITOT)=1
                  ENDIF
!  CASE 1
                  IF(XSAVE(JX,III).GT.XSAVE(JX,JJJ))THEN
                     ITOT=ITOT+1
                     XALL(ITOT)=XNUM1/DENOM
                     KALL(ITOT)=6
                     ITOT=ITOT+1
                     XALL(ITOT)=XNUM3/DENOM
                     KALL(ITOT)=1
                  ENDIF
               ENDIF
!
!  INCORRECT TO CORRECT (CASE 12)
!
               IF(MRITE(JX,JJJ).EQ.0)THEN
                  ITOT=ITOT+1
                  XALL(ITOT)=XNUM1/DENOM
                  KALL(ITOT)=6
                  ITOT=ITOT+1
                  XALL(ITOT)=XNUM3/DENOM
                  KALL(ITOT)=1
               ENDIF
            ENDIF
            IF(LDATA(ILEG,JX).EQ.MCUTS(JX,1))THEN
               KWRONG=KWRONG+1
!
!  CORRECT TO INCORRECT (CASE 10)
!
               IF(MRITE(JX,JJJ).EQ.1)THEN
                  ITOT=ITOT+1
                  XALL(ITOT)=XNUM3/DENOM
                  KALL(ITOT)=6
                  ITOT=ITOT+1
                  XALL(ITOT)=XNUM2/DENOM
                  KALL(ITOT)=1
               ENDIF
!
!  INCORRECT TO INCORRECT (CASES 5 AND 6)
!
               IF(MRITE(JX,JJJ).EQ.0)THEN
!  CASE 6
                  IF(XSAVE(JX,III).LT.XSAVE(JX,JJJ))THEN
                     ITOT=ITOT+1
                     XALL(ITOT)=XNUM2/DENOM
                     KALL(ITOT)=6
                     ITOT=ITOT+1
                     XALL(ITOT)=XNUM3/DENOM
                     KALL(ITOT)=1
                  ENDIF
!  CASE 5
                  IF(XSAVE(JX,III).GT.XSAVE(JX,JJJ))THEN
                     ITOT=ITOT+1
                     XALL(ITOT)=XNUM3/DENOM
                     KALL(ITOT)=6
                     ITOT=ITOT+1
                     XALL(ITOT)=XNUM2/DENOM
                     KALL(ITOT)=1
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
  1   CONTINUE
!
!  IMPOSE RANGE CONSTRAINTS ON ALPHA VECTOR -- ROOT2 IS THE LOWER BOUND
!                                              ROOT1 IS THE UPPER BOUND
!
      KK=0
      DO 363 I=1,ITOT
      IF(XALL(I).GT.ROOT2.AND.XALL(I).LT.ROOT1)THEN
         KK=KK+1
         YALL(KK)=XALL(I)
         LALL(KK)=KALL(I)
      ENDIF
  363 CONTINUE
!
      WSSY=0.0
      IF(KK.GT.0)THEN
         CALL KPRSORT(YALL,KK,LALL)
         IROTC=1
         CALL JAN0PT(KK,NP,NDUAL,YALL,LALL,WSSY,JCH,JEH,JCL,JEL,IROTC)
      ENDIF
!
!  (Z'Z)-1*Z'W-HAT
!
!
      DO 3 K=1,NS
      ZWRONG(K)=0.0
      SUMWS=0.0
      DO 4 JJ=1,NRCALL
      SUMWS=SUMWS+BB(K,JJ)*(XSAVE(JJ,JJJ)+ &
                      WSSY*(XSAVE(JJ,III)-XSAVE(JJ,JJJ)))
   4  CONTINUE
!
      XMAT2(K)=SUMWS
   3  CONTINUE
!
      KKRITE=0
      KKWRONG=0
      KTOTAL=0
      XERR2=0.0
      DO 31 JX=1,NRCALL
      LERROR(ILEG,JX)=0
      SUM=0.0
      DO 32 K=1,NS
      SUM=SUM+XMAT2(K)*ZVEC(JX,K)
  32  CONTINUE
      XXY(JX)=SUM
      XXZ(JX)=SUM
      DB2B1=WS(JX)-XXY(JX)
!
!  CALCULATE CLASSIFICATION ERROR
!
      IF(LDATA(ILEG,JX).NE.0)THEN
         IF(XXY(JX).LT.WS(JX))THEN
            KTOTAL=KTOTAL+1
            IF(LDATA(ILEG,JX).EQ.MCUTS(JX,1))THEN
               KKRITE=KKRITE+1
            ENDIF
            IF(LDATA(ILEG,JX).EQ.MCUTS(JX,2))THEN
               LERROR(ILEG,JX)=1
               KKWRONG=KKWRONG+1
               XERR2=XERR2+(XXY(JX)-WS(JX))**2
               XXZ(JX)=+1.0
               DO 62 K=1,NS
               YWRONG(KKWRONG,K)=XMAT2(K)+1.000*DB2B1*ZVEC(JX,K)
               ZWRONG(K)=ZWRONG(K)+XMAT2(K)+ &
                         1.5000*DB2B1*ZVEC(JX,K)
  62           CONTINUE
            ENDIF
         ENDIF
         IF(XXY(JX).GT.WS(JX))THEN
            KTOTAL=KTOTAL+1
            IF(LDATA(ILEG,JX).EQ.MCUTS(JX,2))THEN
               KKRITE=KKRITE+1
            ENDIF
            IF(LDATA(ILEG,JX).EQ.MCUTS(JX,1))THEN
               LERROR(ILEG,JX)=1
               KKWRONG=KKWRONG+1
               XERR2=XERR2+(XXY(JX)-WS(JX))**2
               XXZ(JX)=-1.0
               DO 63 K=1,NS
               YWRONG(KKWRONG,K)=XMAT2(K)+1.000*DB2B1*ZVEC(JX,K)
               ZWRONG(K)=ZWRONG(K)+XMAT2(K)+ &
                         1.5000*DB2B1*ZVEC(JX,K)
  63           CONTINUE
            ENDIF
         ENDIF
      ENDIF
  31  CONTINUE
!
      DO 5 K=1,NS
      IF(KKWRONG.GT.0)ZWRONG(K)=ZWRONG(K)/FLOAT(KKWRONG)
  5   CONTINUE
!      IF(ILEG.EQ.5)THEN
!      WRITE(21,1013)NII,ILEG,KKRITE,KKWRONG,JEH+JEL,KTOTAL,
!     C                  (XMAT2(K),K=1,3)
!      ENDIF
!
      DEALLOCATE(MRITE)
      DEALLOCATE(LALL)
      DEALLOCATE(YALL)
      DEALLOCATE(XALL)
      DEALLOCATE(KALL)
      DEALLOCATE(XXY)
      DEALLOCATE(XSAVE)
      RETURN
      END
!
!  **************************************************************************
!    SUBROUTINE JAN0PT -- FINDS OPTIMAL CUTTING POINT FOR ONE DIMENSION
!  **************************************************************************
!
      SUBROUTINE JAN0PT(KKNP,NP,NDUAL,YSS,KA,WSSY,JCH,JEH,JCL,JEL,IROTC)
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER:: LD=0,LC=0,LA=0,LB=0
      DOUBLE PRECISION:: AB=0.0,AA=0.0
      DIMENSION YSS(NDUAL),KA(NDUAL)
      INTEGER, ALLOCATABLE :: LE(:)
      INTEGER, ALLOCATABLE :: LJEP(:)
      INTEGER, ALLOCATABLE :: LV(:)
      INTEGER, ALLOCATABLE :: LVB(:)
      INTEGER, ALLOCATABLE :: LEB(:)
      INTEGER, ALLOCATABLE :: LAJEP(:)
      INTEGER, ALLOCATABLE :: LBJEP(:)
      INTEGER, ALLOCATABLE :: LCJEP(:)
      INTEGER, ALLOCATABLE :: LDJEP(:)
      INTEGER, ALLOCATABLE :: MJEP(:)
      DOUBLE PRECISION, ALLOCATABLE :: Z(:)
      DOUBLE PRECISION, ALLOCATABLE :: Y(:)
      DOUBLE PRECISION, ALLOCATABLE :: AAJEP(:)
      DOUBLE PRECISION, ALLOCATABLE :: ABJEP(:)
      DOUBLE PRECISION, ALLOCATABLE :: ABABJEP(:)
      ALLOCATE(LE(NDUAL))
      ALLOCATE(LJEP(NDUAL))
      ALLOCATE(LV(NDUAL))
      ALLOCATE(LVB(NDUAL))
      ALLOCATE(LEB(NDUAL))
      ALLOCATE(LAJEP(101))
      ALLOCATE(LBJEP(101))
      ALLOCATE(LCJEP(101))
      ALLOCATE(LDJEP(101))
      ALLOCATE(MJEP(101))
      ALLOCATE(Z(NDUAL))
      ALLOCATE(Y(NDUAL))
      ALLOCATE(AAJEP(101))
      ALLOCATE(ABJEP(101))
      ALLOCATE(ABABJEP(101))
!
      IF(NP.EQ.1) NP=1  !hack to get rid of warnings
      NPN=KKNP+1
      NPP=KKNP-1
      KCUT=1
      LCUT=6
      NOTE=1
      AA1=0.0
      AB1=0.0
      LA1=0
      LB1=0
      LC1=0
      LD1=0
      AA2=999.0
      AB2=0.0
      LA2=0
      LB2=0
      LC2=0
      LD2=0
      DO 999 III=1,1
      IF(III.EQ.2)THEN
         KCUT=6
         LCUT=1
      ENDIF
!
!  CHECK ALL POSSIBLE INTERIOR CUT POINTS  --  THE NP INPUT POINTS
!      ARE HELD FIXED.  THERE ARE NP POSSIBLE CUT-POINTS BEGINNING
!      WITH CUT-POINT 1 WHICH IS .001 UNITS TO THE LEFT OF POINT 1.
!      CUT-POINT 2 IS BETWEEN POINTS 1 AND 2, ETC.
!
!     1   2   3   4   5   6   7   8   9   10   11 ...... NP-1   NP
!    *  *   *   *   *   *   *   *   *   *    *                *
!    1  2   3   4   5   6   7   8   9  10   11  ...........  NP
!
!  IF KCUT=1 AND LCUT=6, THE FOLLOWING NP PATTERNS ARE TESTED
!
! PATTERN
!   1         6666666666666666666666
!   2         1666666666666666666666
!   3         1166666666666666666666
!   4         1116666666666666666666
!   5         1111666666666666666666
!   6         1111166666666666666666
!   7         1111116666666666666666
!   .           .....
!   .           .....
!   .           .....
!  NP-1       1111111111111111111166
!   NP        1111111111111111111116
!
!  BECAUSE THE PROGRAM TRIES BOTH KCUT=1/LCUT=6 AND KCUT=6/LCUT=1, THIS
!  WILL ALSO TEST THE ONE MISSING PATTERN ABOVE, VIZ., ALL "1"S.
!
!
      KSE=0
      KSV=0
      LSV=0
      LSE=0
      KMARK=1
      I=0
  10  I=I+1
      IF((I-KKNP-1).GE.0)GO TO 12
!  61  Z(I)=999.0
      Z(I)=999.0
      IF(I.EQ.1)THEN
         Y(I)=YSS(1)-.001
      ENDIF
      IF(I.GT.1)THEN
         Y(I)=(YSS(I)+YSS(I-1))/2.0
      ENDIF
!      IF(KA(I).EQ.9)GO TO 10
      IF(KMARK.EQ.1)THEN
         DO 3 J=I,KKNP
         IF(KA(J).EQ.9)GO TO 3
         IF((LCUT-KA(J)).EQ.0)GO TO 5
         IF((KCUT-KA(J)).EQ.0)GO TO 6
         IF((KCUT-KA(J)).NE.0)GO TO 3
  5      LSV=LSV+1
         GO TO 3
  6      LSE=LSE+1
  3      CONTINUE
         KMARK=0
         GO TO 31
      ENDIF
      IF(KA(I-1).EQ.KCUT)THEN
         KSV=KSV+1
         LSE=LSE-1
      ENDIF
      IF(KA(I-1).EQ.LCUT)THEN
         KSE=KSE+1
         LSV=LSV-1
      ENDIF
!
  31  CONTINUE
      LJEP(I)=I
      LV(I)=KSV
      LVB(I)=LSV
      LE(I)=KSE
      LEB(I)=LSE
      KT=LV(I)+LE(I)+LVB(I)+LEB(I)
      Z(I)=FLOAT(LE(I)+LEB(I))/FLOAT(KT)
      GO TO 10
  12  CONTINUE
!
!  FIND BEST CUT POINT
!
      CALL KPRSORT(Z,KKNP,LJEP)
      KIN=1
      MJEP(1)=1
      AAJEP(KIN)=Z(1)
      ABJEP(KIN)=Y(LJEP(1))
      ABABJEP(KIN)=ABS(ABJEP(KIN))
      LAJEP(KIN)=LV(LJEP(1))
      LBJEP(KIN)=LE(LJEP(1))
      LCJEP(KIN)=LVB(LJEP(1))
      LDJEP(KIN)=LEB(LJEP(1))
!
!  CHECK IF THERE ARE MULTIPLE CUT-POINTS WITH SAME CLASSIFICATION AND
!    SELECT THAT CUT-POINT CLOSEST TO THE INTERIOR OF THE SPACE
!
      DO 63 I=2,KKNP
      IF(ABS(Z(1)-Z(I)).LE..000001)THEN
         KIN=KIN+1
         MJEP(KIN)=KIN
         AAJEP(KIN)=Z(I)
         ABJEP(KIN)=Y(LJEP(I))
         ABABJEP(KIN)=ABS(ABJEP(KIN))
         LAJEP(KIN)=LV(LJEP(I))
         LBJEP(KIN)=LE(LJEP(I))
         LCJEP(KIN)=LVB(LJEP(I))
         LDJEP(KIN)=LEB(LJEP(I))
         IF(KIN.GT.100)GO TO 633
         GO TO 63
      ENDIF
      IF(Z(1).LT.Z(I))GO TO 633
  63  CONTINUE
  633 CONTINUE
      IF(KIN.EQ.1)THEN
         AA=AAJEP(1)
         AB=ABJEP(1)
         LA=LAJEP(1)
         LB=LBJEP(1)
         LC=LCJEP(1)
         LD=LDJEP(1)
      ENDIF
      IF(KIN.GT.1)THEN
         CALL KPRSORT(ABABJEP,KIN,MJEP)
         AA=AAJEP(MJEP(1))
         AB=ABJEP(MJEP(1))
         LA=LAJEP(MJEP(1))
         LB=LBJEP(MJEP(1))
         LC=LCJEP(MJEP(1))
         LD=LDJEP(MJEP(1))
      ENDIF
!
      IF(III.EQ.1)THEN
         AA1=AA
         AB1=AB
         LA1=LA
         LB1=LB
         LC1=LC
         LD1=LD
      ENDIF
      IF(III.EQ.2)THEN
         AA2=AA
         AB2=AB
         LA2=LA
         LB2=LB
         LC2=LC
         LD2=LD
      ENDIF
!
  999 CONTINUE
!
      IF(AA1.LE.AA2)THEN
         KCCUT=1
         LCCUT=6
         AA=AA1
         AB=AB1
         LA=LA1
         LB=LB1
         LC=LC1
         LD=LD1
      ENDIF
      IF(AA1.GT.AA2)THEN
         KCCUT=6
         LCCUT=1
         AA=AA2
         AB=AB2
         LA=LA2
         LB=LB2
         LC=LC2
         LD=LD2
      ENDIF
      IF(IROTC.EQ.1)THEN
         KCCUT=1
         LCCUT=6
         AA=AA1
         AB=AB1
         LA=LA1
         LB=LB1
         LC=LC1
         LD=LD1
      ENDIF
      WSSY=AB
      JCL=LA
      JEL=LB
      JCH=LC
      JEH=LD
!
      DEALLOCATE(LE)
      DEALLOCATE(LJEP)
      DEALLOCATE(LV)
      DEALLOCATE(LVB)
      DEALLOCATE(LEB)
      DEALLOCATE(LAJEP)
      DEALLOCATE(LBJEP)
      DEALLOCATE(LCJEP)
      DEALLOCATE(LDJEP)
      DEALLOCATE(MJEP)
      DEALLOCATE(Z)
      DEALLOCATE(Y)
      DEALLOCATE(AAJEP)
      DEALLOCATE(ABJEP)
      DEALLOCATE(ABABJEP)
      RETURN
      END
!
!
!  ************************************************************************
!    SUBROUTINE KPRSORT --SORTS A VECTOR 'A' OF REAL ELEMENTS INTO ASCENDING
!    ORDER.  'LA' IS THE NUMBER OF ELEMENTS TO BE SORTED AND 'IR' IS A
!    VECTOR OF INTEGERS THAT RECORDS THE PERMUTATIONS--USUALLY SET TO
!    1,2,3,4,...
!  ************************************************************************
!
!
      SUBROUTINE KPRSORT(A,LA,IR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(LA),IU(21),IL(21),IR(LA)
      IF (LA.LE.0) RETURN
      M = 1
      I = 1
      J = LA
      R = .375
    5 IF (I.EQ.J) GO TO 45
      IF (R.GT..5898437) GO TO 10
      R = R+3.90625E-2
      GO TO 15
   10 R = R-.21875
   15 K = I
!
! SELECT A CENTRAL ELEMENT OF THE
! ARRAY AND SAVE IT IN LOCATION T
!
      IJ = idint(I+(J-I)*R)
      T = A(IJ)
      IT = IR(IJ)
!
! FIRST ELEMENT OF ARRAY IS GREATER
! THAN T, INTERCHANGE WITH T
!
      IF (A(I).LE.T) GO TO 20
      A(IJ) = A(I)
      A(I) = T
      T = A(IJ)
      IR(IJ) = IR(I)
      IR(I) = IT
      IT = IR(IJ)
   20 L = J
!
! IF LAST ELEMENT OF ARRAY IS LESS THAN
! T, INTERCHANGE WITH T
!
      IF (A(J).GE.T) GO TO 30
      A(IJ) = A(J)
      A(J) = T
      T = A(IJ)
      IR(IJ) = IR(J)
      IR(J) = IT
      IT = IR(IJ)
!
! IF FIRST ELEMENT OF ARRAY IS GREATER
! THAN T, INTERCHANGE WITH T
!
      IF (A(I).LE.T) GO TO 30
      A(IJ) = A(I)
      A(I) = T
      T = A(IJ)
      IR(IJ) = IR(I)
      IR(I) = IT
      IT = IR(IJ)
      GO TO 30
   25 IF (A(L).EQ.A(K)) GO TO 30
      TT = A(L)
      A(L) = A(K)
      A(K) = TT
      ITT = IR(L)
      IR(L) = IR(K)
      IR(K) = ITT
!
! FIND AN ELEMENT IN THE SECOND HALF OF
! THE ARRAY WHICH IS SMALLER THAN T
!
   30 L = L-1
      IF (A(L).GT.T) GO TO 30
!
! FIND AN ELEMENT IN THE FIRST HALF OF
! THE ARRAY WHICH IS GREATER THAN T
!
   35 K = K+1
      IF (A(K).LT.T) GO TO 35
!
! INTERCHANGE THESE ELEMENTS
!
      IF (K.LE.L) GO TO 25
!
! SAVE UPPER AND LOWER SUBSCRIPTS OF
! THE ARRAY YET TO BE SORTED
!
      IF (L-I.LE.J-K) GO TO 40
      IL(M) = I
      IU(M) = L
      I = K
      M = M+1
      GO TO 50
   40 IL(M) = K
      IU(M) = J
      J = L
      M = M+1
      GO TO 50
!
! BEGIN AGAIN ON ANOTHER PORTION OF
! THE UNSORTED ARRAY
!
   45 M = M-1
      IF (M.EQ.0) RETURN
      I = IL(M)
      J = IU(M)
   50 IF (J-I.GE.11) GO TO 15
      IF (I.EQ.1) GO TO 5
      I = I-1
   55 I = I+1
      IF (I.EQ.J) GO TO 45
      T = A(I+1)
      IT = IR(I+1)
      IF (A(I).LE.T) GO TO 55
      K = I
   60 A(K+1) = A(K)
      IR(K+1) = IR(K)
      K = K-1
      IF (T.LT.A(K)) GO TO 60
      A(K+1) = T
      IR(K+1) = IT
      GO TO 55
      END
!
! *********************************************************************
!   SUBROUTINE KPCUTPLANE -- FINDS CUTTING LINE USING THE CUTTING
!                            PLANE PROCEDURE
! *********************************************************************
!
!
      SUBROUTINE KPCUTPLANE(NUMMEMBERS,NUMVOTES,&
                            JJJ,NP,NRCALL,NS,NDUAL,XMAT,ZVEC,WS, &
                            MCUTS,LERROR,IFIXX,KTT,KT,LDATA,IPRINT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION XMAT(NUMMEMBERS,25),ZVEC(NUMVOTES,25),WS(NDUAL), &
                LERROR(NUMMEMBERS,NUMVOTES),&
                MCUTS(NUMVOTES,2),LDATA(NUMMEMBERS,NUMVOTES)
!
      DOUBLE PRECISION, ALLOCATABLE :: XJCH(:)
      DOUBLE PRECISION, ALLOCATABLE :: XJEH(:)
      DOUBLE PRECISION, ALLOCATABLE :: XJCL(:)
      DOUBLE PRECISION, ALLOCATABLE :: XJEL(:)
      DOUBLE PRECISION, ALLOCATABLE :: XPROJ(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: XXY(:)
      DOUBLE PRECISION, ALLOCATABLE :: XXX(:)
      DOUBLE PRECISION, ALLOCATABLE :: ZS(:)
      INTEGER, ALLOCATABLE :: LLL(:)
      INTEGER, ALLOCATABLE :: MM(:)
      INTEGER, ALLOCATABLE :: MVOTE(:)
      INTEGER, ALLOCATABLE :: LLV(:)
      INTEGER, ALLOCATABLE :: LLVB(:)
      INTEGER, ALLOCATABLE :: LLE(:)
      INTEGER, ALLOCATABLE :: LLEB(:)
      INTEGER, ALLOCATABLE :: LCERROR(:,:)
      ALLOCATE(XJCH(25))
      ALLOCATE(XJEH(25))
      ALLOCATE(XJCL(25))
      ALLOCATE(XJEL(25))
      ALLOCATE(XPROJ(NUMMEMBERS,NUMVOTES))
      ALLOCATE(XXY(NP))
      ALLOCATE(XXX(NDUAL))
      ALLOCATE(LLL(NDUAL))
      ALLOCATE(MM(NDUAL))
      ALLOCATE(MVOTE(NDUAL))
      ALLOCATE(LLV(NDUAL))
      ALLOCATE(LLVB(NDUAL))
      ALLOCATE(LLE(NDUAL))
      ALLOCATE(LLEB(NDUAL))
      ALLOCATE(LCERROR(NP,NRCALL))
      ALLOCATE(ZS(NDUAL))
!
!  100 FORMAT(5I5)
!  101 FORMAT(7I4,2F10.4)
! 1010 FORMAT(' WHOA DUDE! THESE DO NOT MATCH!')
! 1093 FORMAT(' CLASSIFICATION CHECK  ',I3,4I8)
! 1094 FORMAT(' RC  CLASSIFICATION ERROR  ',2I3,2I8,2F10.5)
! 1118 FORMAT(7I4,F7.3,3I4)
      NCUT=25
      IF(IPRINT.EQ.1) IPRINT=1  !hack to get rid of warnings
      IF(IFIXX.EQ.1) IFIXX=1  !hack to get rid of warnings
!
!   ESTIMATE PROJECTION VECTORS
!
      IF(JJJ.EQ.1)THEN
         DO 286 JX=1,NRCALL
         DO 285 I=1,NP
         LERROR(I,JX)=0
  285    CONTINUE
  286    CONTINUE
      ENDIF
      KCHECK4=0
      DO 284 JX=1,NRCALL
      DO 283 I=1,NP
      LCERROR(I,JX)=LERROR(I,JX)
      IF(LDATA(I,JX).EQ.0)GO TO 283
      KCHECK4=KCHECK4+LERROR(I,JX)
  283 CONTINUE
  284 CONTINUE
!
      KT=0
      KTT=0
      KTTSAVE=0
      KTSAVE=0
      KCHECK=0
      DO 93 JX=1,NRCALL
!
!  GET YES AND NO COUNTS
!
      KYES=0
      KNO=0
      DO 92 I=1,NP
      IF(LDATA(I,JX).EQ.1)KYES=KYES+1
      IF(LDATA(I,JX).EQ.6)KNO=KNO+1
  92  CONTINUE
!
      DO 89 I=1,NP
      SUM=0.0
      DO 90 K=1,NS
      SUM=SUM+XMAT(I,K)*ZVEC(JX,K)
  90  CONTINUE
!
!  SAVE PROJECTION VECTORS -- LEGISLATOR BY ROLL CALL MATRIX
!
      XPROJ(I,JX)=SUM
      XXY(I)=SUM
      LLL(I)=I
      XXX(I)=SUM
      MM(I)=LDATA(I,JX)
      IF(LDATA(I,JX).EQ.0)MM(I)=9
  89  CONTINUE
!
!  SORT PROJECTION VECTOR (Y-HAT)
!
      CALL KPRSORT(XXX,NP,LLL)
      DO 114 I=1,NP
      MVOTE(I)=MM(LLL(I))
  114 CONTINUE
!
!
!  CALCULATE CLASSIFICATION ERRORS OF PROJECTION ONTO NORMAL VECTOR
!
!
      JCH=0
      JEH=0
      JCL=0
      JEL=0
      IROTC=0
      CALL JAN1PT(NUMMEMBERS,NUMVOTES,&
                  NP,NRCALL,NP,NRCALL,NS,NDUAL,JX,XMAT,XXX,MVOTE,WS, &
                  LLV,LLVB,LLE,LLEB,LERROR, &
                  ZS,JCH,JEH,JCL,JEL,IROTC,KCUT,LCUT,LLL, &
                  XJCH,XJEH,XJCL,XJEL)
!
!      IF(IPRINT.EQ.1)WRITE(*,3909)JX,KYES,KNO,JCH,JCL,JEH,JEL, &
!                      KCUT,LCUT,(ZVEC(JX,K),K=1,NS),WS(JX)
! 3909 FORMAT(I3,'***',6I4,2I2,10F7.3)
      IF(JEH+JEL.EQ.0)THEN
         KT=KT+JCH+JEH+JCL+JEL
         KTSAVE=KTSAVE+JCH+JEH+JCL+JEL
         KITTY1=0
         KITTY2=JCH+JEH+JCL+JEL
         IJUST=0
         GO TO 9377
      ENDIF
!
!  SET-UP FOR GRID SEARCH FOR BEST CUTTING LINE
!
      CALL KPSEARCH(NUMMEMBERS,NUMVOTES,&
                    JX,NCUT,NS,NP,NRCALL,NDUAL,KCUT,LCUT,KTT,KT, &
                  XMAT,ZVEC,XPROJ,WS,XXY, &
                  KITTY1,KITTY2,KYES,KNO,LDATA,LERROR,IPRINT)
!
      KTTSAVE=KTTSAVE+KITTY1
      KTSAVE=KTSAVE+KITTY2
!
 9377 CONTINUE
!
!  STORE DIRECTIONALITY OF ROLL CALL
!
      MCUTS(JX,1)=KCUT
      MCUTS(JX,2)=LCUT
!
!
!  LOCATE ERRORS -- WS(.) CONTAINS THE OPTIMAL CUTTING POINT ON THE
!                   PROJECTION VECTOR -- IT CAN BE USED TO CALCULATE THE
!                   CLASSIFICATION ERRORS
!
      KSUM=0
      DO 108 I=1,NP
      LERROR(I,JX)=0
      XXX(I)=XXY(I)
      LLL(I)=I
      IF(LDATA(I,JX).EQ.0)GO TO 108
      IF(XXY(I).LT.WS(JX))THEN
         IF(LDATA(I,JX).NE.KCUT)THEN
            LERROR(I,JX)=1
            KCHECK=KCHECK+1
            KSUM=KSUM+1
         ENDIF
      ENDIF
      IF(XXY(I).GT.WS(JX))THEN
         IF(LDATA(I,JX).NE.LCUT)THEN
            LERROR(I,JX)=1
            KCHECK=KCHECK+1
            KSUM=KSUM+1
         ENDIF
      ENDIF
  108 CONTINUE
!      IF(KSUM.NE.KITTY1)THEN
!         IF(IPRINT.EQ.1)WRITE(11,1010)
!      ENDIF
      KXERROR=KITTY1
      JXERROR=KITTY1
      SAVEWS=WS(JX)
      XINC=0.2
      CALL KPRSEARCH(NUMMEMBERS,NUMVOTES,&
                     NP,NRCALL,NS,NDUAL,XINC,JX,NCUT,KPCUT,LPCUT, &
                   XMAT,ZVEC,WS,KDOWN,KEQUAL,KUP,JXERROR,WSNEW, &
                   LDATA,LERROR)

      IF(JXERROR.EQ.KXERROR)THEN
         WS(JX)=SAVEWS
      ENDIF
!
!  RESET LERROR(,)
!
      IF(JXERROR.LT.KXERROR)THEN
         KTTSAVE=KTTSAVE-KITTY1+JXERROR
         WS(JX)=WSNEW
         SAVEWS=WS(JX)
         MCUTS(JX,1)=KPCUT
         MCUTS(JX,2)=LPCUT
         KCHECK3=0
         DO 191 I=1,NP
         SUMI=0.0
         DO 192 K=1,NS
         SUMI=SUMI+XMAT(I,K)*ZVEC(JX,K)
  192    CONTINUE
         KCUT=MCUTS(JX,1)
         LCUT=MCUTS(JX,2)
         LERROR(I,JX)=0
         IF(LDATA(I,JX).EQ.0)GO TO 191
         IF(SUMI.LT.WS(JX))THEN
            IF(LDATA(I,JX).NE.KCUT)THEN
               LERROR(I,JX)=1
               KCHECK3=KCHECK3+1
            ENDIF
         ENDIF
         IF(SUMI.GT.WS(JX))THEN
            IF(LDATA(I,JX).NE.LCUT)THEN
               LERROR(I,JX)=1
               KCHECK3=KCHECK3+1
            ENDIF
         ENDIF
  191    CONTINUE
         KXERROR=JXERROR
      ENDIF
!
      IF(JXERROR.GE.KXERROR)THEN
         KCHECK33=0
         DO 391 I=1,NP
         SUMI=0.0
         DO 392 K=1,NS
         SUMI=SUMI+XMAT(I,K)*ZVEC(JX,K)
  392    CONTINUE
         KCUT=MCUTS(JX,1)
         LCUT=MCUTS(JX,2)
         LERROR(I,JX)=0
         IF(LDATA(I,JX).EQ.0)GO TO 391
         IF(SUMI.LT.WS(JX))THEN
            IF(LDATA(I,JX).NE.KCUT)THEN
               LERROR(I,JX)=1
               KCHECK33=KCHECK33+1
            ENDIF
         ENDIF
         IF(SUMI.GT.WS(JX))THEN
            IF(LDATA(I,JX).NE.LCUT)THEN
               LERROR(I,JX)=1
               KCHECK33=KCHECK33+1
            ENDIF
         ENDIF
  391    CONTINUE
      ENDIF
!
!      WRITE(38,1118)JX,JJJ,KYES,KNO,KITTY1,JXERROR,KCHECK3,
!     C                XINC,KDOWN,KEQUAL,KUP
  93  CONTINUE
!
      KT=KTSAVE
      KTT=KTTSAVE
      KCHECK2=0
      KCHECK22=0
      DO 282 I=1,NP
      DO 281 JX=1,NRCALL
      IF(LDATA(I,JX).EQ.0)GO TO 281
      KCHECK2=KCHECK2+LERROR(I,JX)
      KCHECK22=KCHECK22+LCERROR(I,JX)
  281 CONTINUE
  282 CONTINUE
!      IF(IPRINT.EQ.1)THEN
!         WRITE(23,1093)NS,KCHECK,KCHECK2,KCHECK22,KCHECK4
!      ENDIF
      IF(KT.GT.0)THEN
        XERROR=FLOAT(KTT)/FLOAT(KT)
        YERROR=1.0-XERROR
      ENDIF
!      IF(IPRINT.EQ.1)WRITE(23,1094)JJJ,NS,KTT,KT,XERROR,YERROR
      DEALLOCATE(XJCH)
      DEALLOCATE(XJEH)
      DEALLOCATE(XJCL)
      DEALLOCATE(XJEL)
      DEALLOCATE(XPROJ)
      DEALLOCATE(XXY)
      DEALLOCATE(XXX)
      DEALLOCATE(LLL)
      DEALLOCATE(MM)
      DEALLOCATE(MVOTE)
      DEALLOCATE(LLV)
      DEALLOCATE(LLVB)
      DEALLOCATE(LLE)
      DEALLOCATE(LLEB)
      DEALLOCATE(LCERROR)
      DEALLOCATE(ZS)
      RETURN
      END
!
!  ************************************************************************
!    SUBROUTINE KPSEARCH
!  ************************************************************************
!
      SUBROUTINE KPSEARCH(NUMMEMBERS,NUMVOTES,&
                        JX,NCUT,NS,NP,NRCALL,NDUAL,KCUT,LCUT, &
                        KTT,KT,XMAT,ZVEC,XPROJ,WS,XXY, &
                        KITTY1,KITTY2,KYES,KNO,LDATA, &
                        LERROR,IPRINT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION XMAT(NUMMEMBERS,25),ZVEC(NUMVOTES,25),&
                LERROR(NUMMEMBERS,NUMVOTES),&
                XPROJ(NUMMEMBERS,NUMVOTES),WS(NDUAL),XXY(NP),&
                LDATA(NUMMEMBERS,NUMVOTES)
      DOUBLE PRECISION SUM
!
      INTEGER, ALLOCATABLE :: KKKCUT(:)
      INTEGER, ALLOCATABLE :: LLLCUT(:)
      INTEGER, ALLOCATABLE :: LLV(:)
      INTEGER, ALLOCATABLE :: LLVB(:)
      INTEGER, ALLOCATABLE :: LLE(:)
      INTEGER, ALLOCATABLE :: LLEB(:)
      INTEGER, ALLOCATABLE :: LWRONG(:)
      INTEGER, ALLOCATABLE :: LLL(:)
      INTEGER, ALLOCATABLE :: MVOTE(:)
      INTEGER, ALLOCATABLE :: LLM(:)
      INTEGER, ALLOCATABLE :: LLN(:)
      INTEGER, ALLOCATABLE :: MM(:)
      DOUBLE PRECISION, ALLOCATABLE :: XJCH(:)
      DOUBLE PRECISION, ALLOCATABLE :: XJEH(:)
      DOUBLE PRECISION, ALLOCATABLE :: XJCL(:)
      DOUBLE PRECISION, ALLOCATABLE :: XJEL(:)
      DOUBLE PRECISION, ALLOCATABLE :: ZS(:)
      DOUBLE PRECISION, ALLOCATABLE :: FV1(:)
      DOUBLE PRECISION, ALLOCATABLE :: FV2(:)
      DOUBLE PRECISION, ALLOCATABLE :: SUMX(:)
      DOUBLE PRECISION, ALLOCATABLE :: UUUU(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: XXX(:)
      DOUBLE PRECISION, ALLOCATABLE :: IWORK(:)
      DOUBLE PRECISION, ALLOCATABLE ::  Y16MIDP(:,:)
      DOUBLE PRECISION, ALLOCATABLE ::  YHAT(:)
      DOUBLE PRECISION, ALLOCATABLE ::  Z16MIDP(:,:)
      DOUBLE PRECISION, ALLOCATABLE ::  VVV(:,:)
      DOUBLE PRECISION, ALLOCATABLE ::  WORK(:)
      DOUBLE PRECISION, ALLOCATABLE ::  X16MIDP(:,:)
      ALLOCATE(KKKCUT(NDUAL))
      ALLOCATE(LLLCUT(NDUAL))
      ALLOCATE(LLV(NDUAL))
      ALLOCATE(LLVB(NDUAL))
      ALLOCATE(LLE(NDUAL))
      ALLOCATE(LLEB(NDUAL))
      ALLOCATE(LWRONG(NDUAL))
      ALLOCATE(LLL(NDUAL))
      ALLOCATE(MVOTE(NDUAL))
      ALLOCATE(LLM(NDUAL))
      ALLOCATE(LLN(NDUAL))
      ALLOCATE(MM(NDUAL))
      ALLOCATE(XJCH(25))
      ALLOCATE(XJEH(25))
      ALLOCATE(XJCL(25))
      ALLOCATE(XJEL(25))
      ALLOCATE(ZS(NDUAL))
      ALLOCATE(FV1(NDUAL))
      ALLOCATE(FV2(NDUAL))
      ALLOCATE(SUMX(NDUAL))
      ALLOCATE(UUUU(NDUAL,25))
      ALLOCATE(XXX(NDUAL))
      ALLOCATE(IWORK(200))
      ALLOCATE(Y16MIDP(NDUAL,25))
      ALLOCATE(YHAT(NDUAL))
      ALLOCATE(Z16MIDP(NDUAL,25))
      ALLOCATE(VVV(25,25))
      ALLOCATE(WORK(2*NUMMEMBERS+1875))
      ALLOCATE(X16MIDP(NDUAL,25))
!
!  104 FORMAT(I5,10F10.4)
!  210 FORMAT(I5,10F12.3)
! 1091 FORMAT(' INVERSE MATRIX ERROR',I4,I5,I8,2F10.4)
! 1099 FORMAT(I3,I5,I3,2I4)
! 1103 FORMAT(' MIDPOINT DECOMPOSITION',5I6)
! 1212 FORMAT(I3,I5,7I4)
! 3909 FORMAT(I5,I3,6I4,2I8,5I5)
!
      IF(KNO.EQ.1) KNO=1  !hack to get rid of warnings
      IF(IPRINT.EQ.1) IPRINT=1  !hack to get rid of warnings
      IF(KYES.EQ.1) KYES=1  !hack to get rid of warnings
      DO 1 I=1,50
      SUMX(I)=0.0
  1   CONTINUE
!
!  PHASE 2
!
!      NCUT2=20
!
      DO 999 IJL=1,NCUT
!
!  SET-UP FOR PHASE 2
!
!
      DO 388 K=1,NS
      UUUU(IJL,K)=ZVEC(JX,K)
  388 CONTINUE
      DO 389 I=1,NP
      SUM=0.0
      DO 390 K=1,NS
      SUM=SUM+XMAT(I,K)*ZVEC(JX,K)
  390 CONTINUE
!
!  SAVE PROJECTION VECTORS -- LEGISLATOR BY ROLL CALL MATRIX
!
      XPROJ(I,JX)=SUM
      XXY(I)=SUM
      LLL(I)=I
      XXX(I)=SUM
      MM(I)=LDATA(I,JX)
      IF(LDATA(I,JX).EQ.0)MM(I)=9
  389 CONTINUE
!
!  SORT PROJECTION VECTOR (Y-HAT)
!
!
      CALL KPRSORT(XXX,NP,LLL)
      DO 314 I=1,NP
      MVOTE(I)=MM(LLL(I))
  314 CONTINUE
!
!
!  CALCULATE CLASSIFICATION ERRORS FOR BEST SOLUTION FROM PHASE 1
!
!
      JCH=0
      JEH=0
      JCL=0
      JEL=0
      IROTC=0
      CALL JAN1PT(NUMMEMBERS,NUMVOTES,&
                  NP,NRCALL,NP,NRCALL,NS,NDUAL,JX,XMAT,XXX,MVOTE,WS, &
                  LLV,LLVB,LLE,LLEB,LERROR, &
                  ZS,JCH,JEH,JCL,JEL,IROTC,KCUT,LCUT,LLL, &
                  XJCH,XJEH,XJCL,XJEL)
!
!      IF(IPRINT.EQ.1)WRITE(11,3909)JX,IJL,KYES,KNO,JCH,JCL,JEH,JEL
!
      LLM(IJL)=IJL
      LLN(IJL)=JEH+JEL
      FV1(IJL)=FLOAT(JEH+JEL)
      FV2(IJL)=WS(JX)
      KKKCUT(IJL)=KCUT
      LLLCUT(IJL)=LCUT
!
      IF(JEH+JEL.EQ.0)THEN
         KT=KT+JCH+JCL+JEH+JEL
         KITTY1=0
         KITTY2=JCH+JCL+JEH+JEL
         IJUST=2
         DEALLOCATE(KKKCUT)
         DEALLOCATE(LLLCUT)
         DEALLOCATE(LLV)
         DEALLOCATE(LLVB)
         DEALLOCATE(LLE)
         DEALLOCATE(LLEB)
         DEALLOCATE(LWRONG)
         DEALLOCATE(LLL)
         DEALLOCATE(MVOTE)
         DEALLOCATE(LLM)
         DEALLOCATE(LLN)
         DEALLOCATE(MM)
         DEALLOCATE(XJCH)
         DEALLOCATE(XJEH)
         DEALLOCATE(XJCL)
         DEALLOCATE(XJEL)
         DEALLOCATE(ZS)
         DEALLOCATE(FV1)
         DEALLOCATE(FV2)
         DEALLOCATE(SUMX)
         DEALLOCATE(UUUU)
         DEALLOCATE(XXX)
         DEALLOCATE(IWORK)
         DEALLOCATE(Y16MIDP)
         DEALLOCATE(YHAT)
         DEALLOCATE(Z16MIDP)
         DEALLOCATE(VVV)
         DEALLOCATE(WORK)
         DEALLOCATE(X16MIDP)
         RETURN
      ENDIF
!
!
      KASTRO=4*(JEH+JEL)
      IF(KASTRO.LT.4*NS)KASTRO=4*NS
      IF(KASTRO.GT.NP)KASTRO=NP
!
      DO 108 I=1,NP
      LWRONG(I)=0
      DB2B1=WS(JX)-XXY(I)
      IF(XXY(I).LT.WS(JX))THEN
!
!  IF CORRECT PLACE LEGISLATOR POINT ON THE CURRENT CUTTING PLANE
!
         IF(LDATA(I,JX).EQ.KCUT)THEN
            DO 109 K=1,NS
            Y16MIDP(I,K)=XMAT(I,K)+DB2B1*ZVEC(JX,K)
  109       CONTINUE
         ENDIF
!
!  IF INCORRECT PUT ACTUAL POINT INTO THE CUTTING CLOUD
!
         IF(LDATA(I,JX).EQ.LCUT)THEN
            LWRONG(I)=1
            DO 110 K=1,NS
            Y16MIDP(I,K)=XMAT(I,K)
  110       CONTINUE
         ENDIF
!
!  IF NOT-VOTING PUT LEGISLATOR POINT ON THE CURRRENT CUTTING PLANE
!
         IF(LDATA(I,JX).EQ.0)THEN
            DO 111 K=1,NS
            Y16MIDP(I,K)=XMAT(I,K)+DB2B1*ZVEC(JX,K)
  111       CONTINUE
         ENDIF
      ENDIF
      IF(XXY(I).GT.WS(JX))THEN
!
!  IF CORRECT PLACE LEGISLATOR POINT ON THE CURRENT CUTTING PLANE
!
         IF(LDATA(I,JX).EQ.LCUT)THEN
            DO 112 K=1,NS
            Y16MIDP(I,K)=XMAT(I,K)+DB2B1*ZVEC(JX,K)
  112       CONTINUE
         ENDIF
!
!  IF INCORRECT PUT ACTUAL POINT INTO THE CUTTING CLOUD
!
         IF(LDATA(I,JX).EQ.KCUT)THEN
            LWRONG(I)=1
            DO 113 K=1,NS
            Y16MIDP(I,K)=XMAT(I,K)
  113       CONTINUE
         ENDIF
!
!  IF NOT-VOTING PUT LEGISLATOR POINT ON THE CURRRENT CUTTING PLANE
!
         IF(LDATA(I,JX).EQ.0)THEN
            DO 214 K=1,NS
            Y16MIDP(I,K)=XMAT(I,K)+DB2B1*ZVEC(JX,K)
  214       CONTINUE
         ENDIF
      ENDIF
!
  108 CONTINUE
!
!  MASS CENTER THE CUTTING PLANE MATRIX (Y16MIDP(,) HAS ALL POINTS)
!
      DO 215 K=1,NS
      SUM=0.0
      DO 216 I=1,NP
      SUM=SUM+Y16MIDP(I,K)
  216 CONTINUE
      DO 217 I=1,NP
      Y16MIDP(I,K)=Y16MIDP(I,K)-SUM/FLOAT(NP)
      SUMX(K)=SUMX(K)+Y16MIDP(I,K)**2
  217 CONTINUE
      SUMX(K)=SUMX(K)/FLOAT(NP)
  215 CONTINUE
!
!  CONSTRUCT PARTIAL CUTTING PLANE MATRIX (X16MIDP(,))
!
      KK=0
      KHIT=0
!
      DO 316 I=1,NP
      IF(LWRONG(I).EQ.1)THEN
         KK=KK+1
         DO 317 K=1,NS
         X16MIDP(KK,K)=Y16MIDP(I,K)
  317    CONTINUE
      ENDIF
  316 CONTINUE
      DO 201 I=1,NP
      IF(LWRONG(I).EQ.0)THEN
         KK=KK+1
         DO 219 K=1,NS
         X16MIDP(KK,K)=Y16MIDP(I,K)
  219    CONTINUE
         IF(KK.EQ.KASTRO)GO TO 203
      ENDIF
  201 CONTINUE
  203 CONTINUE
!
!  MASS CENTER THE PARTIAL CUTTING PLANE MATRIX
!
      DO 815 K=1,NS
      SUM=0.0
      DO 816 I=1,KASTRO
      SUM=SUM+X16MIDP(I,K)
  816 CONTINUE
      DO 817 I=1,KASTRO
      X16MIDP(I,K)=X16MIDP(I,K)-SUM/FLOAT(KASTRO)
      SUMX(K+NS)=SUMX(K+NS)+X16MIDP(I,K)**2
  817 CONTINUE
      SUMX(K+NS)=SUMX(K+NS)/FLOAT(KASTRO)
  815 CONTINUE
!
!  RUN REGRESSION TO ELIMINATE DIMENSION WITH LEAST VARIANCE
!
!
!  CALL SINGULAR VALUE DECOMPOSITION ROUTINE
!
      LWORK=2*NUMMEMBERS+1875
      XTOL=.001
!      CALL LSVRR(NP,NS,Y16MIDP,NP,21,XTOL,IRANK,YHAT,Y16MIDP,
!     C           NP,VVV,25)
      CALL DGESDD('S',NP,NS,Y16MIDP,NDUAL,YHAT,Z16MIDP, &
                 NDUAL,VVV,25,WORK,LWORK,IWORK,IRANK)
!
!      WRITE(23,1094)IRANK
! 1094 FORMAT(' DGESDD ROUTINE',I5)
!      WRITE(23,3908)JX,IJL,(YHAT(K),K=1,NS),(VVV(K,NS),K=1,NS)
      DO 115 K=1,NS
      SUMX(K)=SUMX(K+NS)
!
!  WRONG WAY -- DGESDD RETURNS V_transpose, NOT V
!
!      ZVEC(JX,K)=VVV(K,NS)
!
!  RIGHT WAY -- DGESDD RETURNS V_transpose, NOT V, SO THE
!               S_th ROW IS TRANSFERRED
!
      ZVEC(JX,K)=VVV(NS,K)
  115 CONTINUE
!
! 3908 FORMAT(I5,I3,10F7.3)
!
!
!  RUN REGRESSION TO ELIMINATE DIMENSION WITH LEAST VARIANCE
!
!      CALL LSVRR(KASTRO,NS,X16MIDP,NP,21,XTOL,IRANK,YHAT,X16MIDP,
!     C           NP,VVV,25)
      CALL DGESDD('S',KASTRO,NS,X16MIDP,NDUAL,YHAT,Z16MIDP, &
                 NDUAL,VVV,25,WORK,LWORK,IWORK,IRANK)
!
!      WRITE(23,1094)IRANK
!      WRITE(23,3908)JX,IJL,(YHAT(K),K=1,NS),(VVV(K,NS),K=1,NS)
      IF(IJL.GT.25)THEN
         DO 114 K=1,NS
!
!  WRONG WAY -- DGESDD RETURNS V_transpose, NOT V
!
!         ZVEC(JX,K)=VVV(K,NS)
!
!  RIGHT WAY -- DGESDD RETURNS V_transpose, NOT V, SO THE
!               S_th ROW IS TRANSFERRED
!
         ZVEC(JX,K)=VVV(NS,K)
  114    CONTINUE
      ENDIF
!
!
  999 CONTINUE
!
      CALL KPRSORT(FV1,NCUT,LLM)
!
      DO 281 JJ=1,NCUT
      IF(FV1(1).LT.FV1(JJ))GO TO 282
  281 CONTINUE
  282 KIN=JJ-1
      LLM(1)=LLM(KIN)
!
      DO 387 K=1,NS
      ZVEC(JX,K)=UUUU(LLM(1),K)
  387 CONTINUE
      WS(JX)=FV2(LLM(1))
      KCUT=KKKCUT(LLM(1))
      LCUT=LLLCUT(LLM(1))
      DO 137 I=1,NP
      SUM=0.0
      DO 138 K=1,NS
      SUM=SUM+XMAT(I,K)*ZVEC(JX,K)
  138 CONTINUE
      XPROJ(I,JX)=SUM
      XXY(I)=SUM
  137 CONTINUE
      KTT=KTT+LLN(LLM(1))
      KITTY1=LLN(LLM(1))
      IJUST=3
      KT=KT+JCH+JCL+JEH+JEL
      KITTY2=JCH+JCL+JEH+JEL
      DEALLOCATE(KKKCUT)
      DEALLOCATE(LLLCUT)
      DEALLOCATE(LLV)
      DEALLOCATE(LLVB)
      DEALLOCATE(LLE)
      DEALLOCATE(LLEB)
      DEALLOCATE(LWRONG)
      DEALLOCATE(LLL)
      DEALLOCATE(MVOTE)
      DEALLOCATE(LLM)
      DEALLOCATE(LLN)
      DEALLOCATE(MM)
      DEALLOCATE(XJCH)
      DEALLOCATE(XJEH)
      DEALLOCATE(XJCL)
      DEALLOCATE(XJEL)
      DEALLOCATE(ZS)
      DEALLOCATE(FV1)
      DEALLOCATE(FV2)
      DEALLOCATE(SUMX)
      DEALLOCATE(UUUU)
      DEALLOCATE(XXX)
      DEALLOCATE(IWORK)
      DEALLOCATE(Y16MIDP)
      DEALLOCATE(YHAT)
      DEALLOCATE(Z16MIDP)
      DEALLOCATE(VVV)
      DEALLOCATE(WORK)
      DEALLOCATE(X16MIDP)
      RETURN
      END
!
!  **************************************************************************
!    SUBROUTINE JAN1PT -- FINDS OPTIMAL CUTTING POINT FOR ONE DIMENSION
!  **************************************************************************
!
      SUBROUTINE JAN1PT(NUMMEMBERS,NUMVOTES,&
                        NPZZ,NV,NP,NRCALL,NS,NDUAL,IVOT,XMAT,YSS,KA,WS, &
                        LLV,LLVB,LLE,LLEB, &
                        LERROR,ZS,JCH,JEH,JCL,JEL,IROTC,KCCUT,LCCUT, &
                        LLL,XJCH,XJEH,XJCL,XJEL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER:: LD=0,LC=0,LA=0,LB=0
      DOUBLE PRECISION:: AB=0.0,AA=0.0
      DIMENSION YSS(NDUAL),KA(NDUAL),WS(NDUAL),LV(NDUAL), &
                LEB(NDUAL),Z(NDUAL),Y(NDUAL),LLV(NDUAL), &
                LLVB(NDUAL),LE(NDUAL),LERROR(NUMMEMBERS,NUMVOTES), &
                LLEB(NDUAL),ZS(NDUAL),LLL(NDUAL),XMAT(NUMMEMBERS,25), &
                XJCH(25),XJEH(25),XJCL(25),XJEL(25),AAJEP(101), &
                ABJEP(101),LAJEP(101),LBJEP(101),LCJEP(101), &
                LDJEP(101),ABABJEP(101),MJEP(101), &
                LVB(NDUAL),LLE(NDUAL),LJEP(NDUAL)
!
      JROTC=1
!      open(1,file='dude.txt',POSITION='APPEND')
!      do 1565 i=1,np
!      write(1,1566)i,YSS(i),ka(i)
! 1566 format(i5,f8.3,i3)
! 1565 continue
!      close(1)
      KPDUDE=NP
      KPDUDE2=NRCALL
      IF(IROTC.EQ.2)THEN
         JROTC=0
         IROTC=1
      ENDIF
      NPN=NPZZ+1
      NPP=NPZZ-1
      KCUT=1
      LCUT=6
      NOTE=2
      IF(IROTC.EQ.1)THEN
         NOTE=1
      ENDIF
      AA1=0.0
      AB1=0.0
      LA1=0
      LB1=0
      LC1=0
      LD1=0
      AA2=999.0
      AB2=0.0
      LA2=0
      LB2=0
      LC2=0
      LD2=0
      DO 999 III=1,NOTE
      IF(III.EQ.2)THEN
         KCUT=6
         LCUT=1
      ENDIF
!
!  CHECK ALL POSSIBLE INTERIOR CUT POINTS  --  THE NP INPUT POINTS
!      ARE HELD FIXED.  THERE ARE NP POSSIBLE CUT-POINTS BEGINNING
!      WITH CUT-POINT 1 WHICH IS .001 UNITS TO THE LEFT OF POINT 1.
!      CUT-POINT 2 IS BETWEEN POINTS 1 AND 2, ETC.
!
!     1   2   3   4   5   6   7   8   9   10   11 ...... NP-1   NP
!    *  *   *   *   *   *   *   *   *   *    *                *
!    1  2   3   4   5   6   7   8   9  10   11  ...........  NP
!
!  IF KCUT=1 AND LCUT=6, THE FOLLOWING NP PATTERNS ARE TESTED
!
! PATTERN
!   1         6666666666666666666666
!   2         1666666666666666666666
!   3         1166666666666666666666
!   4         1116666666666666666666
!   5         1111666666666666666666
!   6         1111166666666666666666
!   7         1111116666666666666666
!   .           .....
!   .           .....
!   .           .....
!  NP-1       1111111111111111111166
!   NP        1111111111111111111116
!
!  BECAUSE THE PROGRAM TRIES BOTH KCUT=1/LCUT=6 AND KCUT=6/LCUT=1, THIS
!  WILL ALSO TEST THE ONE MISSING PATTERN ABOVE, VIZ., ALL "1"S.
!
!
      KSE=0
      KSV=0
      LSV=0
      LSE=0
      KMARK=1
      I=0
  10  I=I+1
!      IF(I-NPZZ-1)61,12,12
      IF((I-NPZZ-1).GE.0)GO TO 12
!  61  Z(I)=999.0
      Z(I)=999.0
      IF(I.EQ.1)THEN
         Y(I)=YSS(1)-.001
      ENDIF
      IF(I.GT.1)THEN
         Y(I)=(YSS(I)+YSS(I-1))/2.0
      ENDIF
!      IF(KA(I).EQ.9)GO TO 10
      IF(KMARK.EQ.1)THEN
         DO 3 J=I,NPZZ
         IF(KA(J).EQ.9)GO TO 3
         IF((LCUT-KA(J)).EQ.0)GO TO 5
         IF((KCUT-KA(J)).EQ.0)GO TO 6
         IF((KCUT-KA(J)).NE.0)GO TO 3
  5      LSV=LSV+1
         GO TO 3
  6      LSE=LSE+1
  3      CONTINUE
         KMARK=0
         GO TO 31
      ENDIF
      IF(KA(I-1).EQ.KCUT)THEN
         KSV=KSV+1
         LSE=LSE-1
      ENDIF
      IF(KA(I-1).EQ.LCUT)THEN
         KSE=KSE+1
         LSV=LSV-1
      ENDIF
!
  31  CONTINUE
      LJEP(I)=I
      LV(I)=KSV
      LVB(I)=LSV
      LE(I)=KSE
      LEB(I)=LSE
      KT=LV(I)+LE(I)+LVB(I)+LEB(I)
      Z(I)=FLOAT(LE(I)+LEB(I))/FLOAT(KT)
!
      IF(JROTC.EQ.0)THEN
         ZS(I)=Y(I)
         LLV(I)=LV(I)
         LLE(I)=LE(I)
         LLVB(I)=LVB(I)
         LLEB(I)=LEB(I)
      ENDIF
      GO TO 10
  12  CONTINUE
!
!  FIND BEST CUT POINT
!
      CALL KPRSORT(Z,NPZZ,LJEP)
      KIN=1
      MJEP(1)=1
      AAJEP(KIN)=Z(1)
      ABJEP(KIN)=Y(LJEP(1))
      ABABJEP(KIN)=ABS(ABJEP(KIN))
      LAJEP(KIN)=LV(LJEP(1))
      LBJEP(KIN)=LE(LJEP(1))
      LCJEP(KIN)=LVB(LJEP(1))
      LDJEP(KIN)=LEB(LJEP(1))
!
!  CHECK IF THERE ARE MULTIPLE CUT-POINTS WITH SAME CLASSIFICATION AND
!    SELECT THAT CUT-POINT CLOSEST TO THE INTERIOR OF THE SPACE
!
      DO 63 I=2,NPZZ
      IF(ABS(Z(1)-Z(I)).LE..00001)THEN
         KIN=KIN+1
         MJEP(KIN)=KIN
         AAJEP(KIN)=Z(I)
         ABJEP(KIN)=Y(LJEP(I))
         ABABJEP(KIN)=ABS(ABJEP(KIN))
         LAJEP(KIN)=LV(LJEP(I))
         LBJEP(KIN)=LE(LJEP(I))
         LCJEP(KIN)=LVB(LJEP(I))
         LDJEP(KIN)=LEB(LJEP(I))
         IF(KIN.GT.100)GO TO 633
         GO TO 63
      ENDIF
      IF(Z(1).LT.Z(I))GO TO 633
  63  CONTINUE
  633 CONTINUE
      IF(KIN.EQ.1)THEN
         AA=AAJEP(1)
         AB=ABJEP(1)
         LA=LAJEP(1)
         LB=LBJEP(1)
         LC=LCJEP(1)
         LD=LDJEP(1)
      ENDIF
      IF(KIN.GT.1)THEN
         CALL KPRSORT(ABABJEP,KIN,MJEP)
         AA=AAJEP(MJEP(1))
         AB=ABJEP(MJEP(1))
         LA=LAJEP(MJEP(1))
         LB=LBJEP(MJEP(1))
         LC=LCJEP(MJEP(1))
         LD=LDJEP(MJEP(1))
      ENDIF
!
      IF(III.EQ.1)THEN
         AA1=AA
         AB1=AB
         LA1=LA
         LB1=LB
         LC1=LC
         LD1=LD
      ENDIF
      IF(III.EQ.2)THEN
         AA2=AA
         AB2=AB
         LA2=LA
         LB2=LB
         LC2=LC
         LD2=LD
      ENDIF
!
  999 CONTINUE
!
      IF(AA1.LE.AA2)THEN
         KCCUT=1
         LCCUT=6
         AA=AA1
         AB=AB1
         LA=LA1
         LB=LB1
         LC=LC1
         LD=LD1
      ENDIF
      IF(AA1.GT.AA2)THEN
         KCCUT=6
         LCCUT=1
         AA=AA2
         AB=AB2
         LA=LA2
         LB=LB2
         LC=LC2
         LD=LD2
      ENDIF
      IF(IROTC.EQ.1)THEN
         KCCUT=1
         LCCUT=6
         AA=AA1
         AB=AB1
         LA=LA1
         LB=LB1
         LC=LC1
         LD=LD1
      ENDIF
      WS(IVOT)=AB
      IF(IROTC.EQ.1)WS(IVOT+NV)=AB
      IF(JROTC.EQ.1)THEN
         ZS(IVOT)=AA
         LLV(IVOT)=LA
         LLE(IVOT)=LB
         LLVB(IVOT)=LC
         LLEB(IVOT)=LD
      ENDIF
      JCL=LA
      JEL=LB
      JCH=LC
      JEH=LD
!
      IF(IROTC.EQ.0)THEN
         DO 71 K=1,NS
         XJCH(K)=0.0
         XJEH(K)=0.0
         XJCL(K)=0.0
         XJEL(K)=0.0
  71     CONTINUE
         DO 64 I=1,NPZZ
         IF(LLL(I).LE.NPZZ-1)LERROR(LLL(I),IVOT)=0
         IF(KA(I).EQ.9)GO TO 64
         LERROR(LLL(I),IVOT)=0
         IF(YSS(I).LT.AB)THEN
            IF(KA(I).EQ.KCCUT)THEN
               LERROR(LLL(I),IVOT)=0
               DO 70 K=1,NS
               XJCL(K)=XJCL(K)+XMAT(LLL(I),K)
   70          CONTINUE
            ENDIF
            IF(KA(I).EQ.LCCUT)THEN
               LERROR(LLL(I),IVOT)=1
               DO 72 K=1,NS
               XJEL(K)=XJEL(K)+XMAT(LLL(I),K)
   72          CONTINUE
            ENDIF
         ENDIF
         IF(YSS(I).GT.AB)THEN
            IF(KA(I).EQ.LCCUT)THEN
               LERROR(LLL(I),IVOT)=0
               DO 73 K=1,NS
               XJCH(K)=XJCH(K)+XMAT(LLL(I),K)
   73          CONTINUE
            ENDIF
            IF(KA(I).EQ.KCCUT)THEN
               LERROR(LLL(I),IVOT)=1
               DO 74 K=1,NS
               XJEH(K)=XJEH(K)+XMAT(LLL(I),K)
   74          CONTINUE
            ENDIF
         ENDIF
  64     CONTINUE
         DO 75 K=1,NS
         IF(JCL.GT.0)XJCL(K)=XJCL(K)/FLOAT(JCL)
         IF(JEL.GT.0)XJEL(K)=XJEL(K)/FLOAT(JEL)
         IF(JCH.GT.0)XJCH(K)=XJCH(K)/FLOAT(JCH)
         IF(JEH.GT.0)XJEH(K)=XJEH(K)/FLOAT(JEH)
  75     CONTINUE
      ENDIF
      IF(IROTC.EQ.1)THEN
         DO 65 I=1,NPZZ
         IF(LLL(I).LE.NPZZ-1)LERROR(IVOT,LLL(I))=0
         IF(KA(I).EQ.9)GO TO 65
         LERROR(IVOT,LLL(I))=0
         IF(YSS(I).LT.AB)THEN
            IF(KA(I).EQ.KCCUT)LERROR(IVOT,LLL(I))=0
            IF(KA(I).EQ.LCCUT)LERROR(IVOT,LLL(I))=1
         ENDIF
         IF(YSS(I).GT.AB)THEN
            IF(KA(I).EQ.LCCUT)LERROR(IVOT,LLL(I))=0
            IF(KA(I).EQ.KCCUT)LERROR(IVOT,LLL(I))=1
         ENDIF
  65     CONTINUE
      ENDIF
      RETURN
      END
!
!  ************************************************************************
!    SUBROUTINE KPRSEARCH -- DOES LOCAL SEARCH ON NORMAL VECTORS
!
!  ************************************************************************
!
      SUBROUTINE KPRSEARCH(NUMMEMBERS,NUMVOTES,&
                           NP,NRCALL,NS,NDUAL,XINC,JX,NCUT,KPCUT,LPCUT, &
                      XMAT,ZVEC,WS,KDOWN,KEQUAL,KUP,JXERROR, &
                      WSNEW,LDATA,LERROR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION XMAT(NUMMEMBERS,25),ZVEC(NUMVOTES,25),WS(NDUAL), &
                LERROR(NUMMEMBERS,NUMVOTES),LDATA(NUMMEMBERS,NUMVOTES)
!
      INTEGER, ALLOCATABLE :: LLL(:)
      INTEGER, ALLOCATABLE :: MVOTE(:)
      INTEGER, ALLOCATABLE :: LLV(:)
      INTEGER, ALLOCATABLE :: LLVB(:)
      INTEGER, ALLOCATABLE :: LLE(:)
      INTEGER, ALLOCATABLE :: LLEB(:)
      INTEGER, ALLOCATABLE :: MM(:)
      DOUBLE PRECISION, ALLOCATABLE :: XJCH(:)
      DOUBLE PRECISION, ALLOCATABLE :: XJEH(:)
      DOUBLE PRECISION, ALLOCATABLE :: XJCL(:)
      DOUBLE PRECISION, ALLOCATABLE :: XJEL(:)
      DOUBLE PRECISION, ALLOCATABLE :: ZS(:)
      DOUBLE PRECISION, ALLOCATABLE :: UUU(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: XXX(:)
      DOUBLE PRECISION, ALLOCATABLE :: ZZZ(:)
      DOUBLE PRECISION, ALLOCATABLE :: ZVEC2(:,:)
      ALLOCATE(LLL(NDUAL))
      ALLOCATE(MVOTE(NDUAL))
      ALLOCATE(LLV(NDUAL))
      ALLOCATE(LLVB(NDUAL))
      ALLOCATE(LLE(NDUAL))
      ALLOCATE(LLEB(NDUAL))
      ALLOCATE(MM(NP))
      ALLOCATE(XJCH(25))
      ALLOCATE(XJEH(25))
      ALLOCATE(XJCL(25))
      ALLOCATE(XJEL(25))
      ALLOCATE(ZS(NDUAL))
      ALLOCATE(UUU(NDUAL,25))
      ALLOCATE(XXX(NDUAL))
      ALLOCATE(ZZZ(NP))
      ALLOCATE(ZVEC2(NRCALL,25))
!
!  210 FORMAT(I5,10F12.3)
! 1091 FORMAT(' INVERSE MATRIX ERROR',I4,I5,I8,2F10.4)
! 1099 FORMAT(I3,I5,I3,2I4)
! 1103 FORMAT(' MIDPOINT DECOMPOSITION',5I6)
! 1212 FORMAT(I3,I5,7I4)
! 3909 FORMAT(I5,I3,6I4,2I8,5I5)
!
!      XINC=0.05
!
!
!
      IF(XINC.EQ.1.0) XINC=1.0  !hack to get rid of warnings
      KDOWN=0
      KEQUAL=0
      KUP=0
      DO 998 IJL=1,NCUT
      KQUIT=IJL
!
!  SET-UP FOR PHASE 2
!
      SUM=0.0
      DO 3 K=1,NS
!      ZZZ(K)=(URAND(ISEED)-.50)*0.4 + ZVEC(JX,K)
      CALL RANDOM_NUMBER(RKEITH)
      ZZZ(K)=(RKEITH-.50)*0.4 + ZVEC(JX,K)
!      ZZZ(K)=(RNUNF()-.50)*0.4 + ZVEC(JX,K)
!      ZZZ(K)=(Rand()-.50)*0.4 + ZVEC(JX,K)
!      ZZZ(K)=0.7*0.4 + ZVEC(JX,K)
      SUM=SUM+ZZZ(K)**2
  3   CONTINUE
      SUM2=0.0
      DO 4 K=1,NS
      ZZZ(K)=ZZZ(K)/SQRT(SUM)
      SUM2=SUM2+(ZVEC(JX,K)-ZZZ(K))**2
  4   CONTINUE
      SUM2=SQRT(SUM2)
      SUM3=0.0
      DO 5 K=1,NS
      ZVEC2(JX,K)=ZVEC(JX,K)+(XINC/SUM2)*(ZZZ(K)-ZVEC(JX,K))
!      ZVEC2(JX,K)=(Rand()-.50)
      SUM3=SUM3+ZVEC2(JX,K)**2
  5   CONTINUE
      DO 6 K=1,NS
      ZVEC2(JX,K)=ZVEC2(JX,K)/SQRT(SUM3)
  6   CONTINUE
!
      DO 488 K=1,NS
      UUU(IJL,K)=ZVEC2(JX,K)
  488 CONTINUE
      DO 489 I=1,NP
      SUM=0.0
      DO 490 K=1,NS
      SUM=SUM+XMAT(I,K)*ZVEC2(JX,K)
  490 CONTINUE
!
!  SAVE PROJECTION VECTORS -- LEGISLATOR BY ROLL CALL MATRIX
!
      LLL(I)=I
      XXX(I)=SUM
      MM(I)=LDATA(I,JX)
      IF(LDATA(I,JX).EQ.0)MM(I)=9
  489 CONTINUE
!
!  SORT PROJECTION VECTOR (Y-HAT)
!
!
      CALL KPRSORT(XXX,NP,LLL)
      DO 414 I=1,NP
      MVOTE(I)=MM(LLL(I))
  414 CONTINUE
!
!
!  CALCULATE CLASSIFICATION ERRORS FOR BEST SOLUTION FROM PHASE 1
!
!
      JCH=0
      JEH=0
      JCL=0
      JEL=0
      IROTC=0
      CALL JAN1PT(NUMMEMBERS,NUMVOTES,&
                  NP,NRCALL,NP,NRCALL,NS,NDUAL,JX,XMAT,XXX,MVOTE,WS, &
                  LLV,LLVB,LLE,LLEB,LERROR, &
                  ZS,JCH,JEH,JCL,JEL,IROTC,KCUT,LCUT,LLL, &
                  XJCH,XJEH,XJCL,XJEL)
!
      IF(JEH+JEL.LT.JXERROR)THEN
         KDOWN=KDOWN+1
         JXERROR=JEH+JEL
         DO 997 K=1,NS
         ZVEC(JX,K)=ZVEC2(JX,K)
  997    CONTINUE
         WSNEW=WS(JX)
         KPCUT=KCUT
         LPCUT=LCUT
         GO TO 998
      ENDIF
      IF(JEH+JEL.EQ.JXERROR)KEQUAL=KEQUAL+1
      IF(JEH+JEL.GT.JXERROR)KUP=KUP+1
!
  998 CONTINUE
!
      DEALLOCATE(LLL)
      DEALLOCATE(MVOTE)
      DEALLOCATE(LLV)
      DEALLOCATE(LLVB)
      DEALLOCATE(LLE)
      DEALLOCATE(LLEB)
      DEALLOCATE(MM)
      DEALLOCATE(XJCH)
      DEALLOCATE(XJEH)
      DEALLOCATE(XJCL)
      DEALLOCATE(XJEL)
      DEALLOCATE(ZS)
      DEALLOCATE(UUU)
      DEALLOCATE(XXX)
      DEALLOCATE(ZZZ)
      DEALLOCATE(ZVEC2)
      RETURN
      END
