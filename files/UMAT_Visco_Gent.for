C==========================================================
C   Material subroutine for standard visco-elastomer material
C   (ABAQUS-UMAT-version)
C   ---- Gent hyperelastic model
C   ---- Newtonian fluid model
C  
C   copyright (c), Yanhui Jiang, NJUST, China
C
C   This code may be used solely for education and research. 
C   You may make copies of this code, but the copyright notice must appear in
C   all copies. Any other use of this code requires written permission.
C   Your use of this code means an implicit agreement to the above conditions.
C========================================================== 

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME

      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

      DIMENSION MY_DHASH(9,2)
      INTEGER :: MY_DHASH
      INTEGER :: rs,pq,r,s,p,q,m,n,i,j,k,l,mn,IJ,KL,K1,K2,K3,K4,ITE_STEP

      DOUBLE PRECISION :: P_Mu,P_N,P_LAMDA,P_Mu_V,P_N_V,P_LAMDA_V,
     1 TAO_D,TAO_V,ETA
      DOUBLE PRECISION :: TERM1,TERM2,TERM3,TERM4,TERM5,TM
      DOUBLE PRECISION :: TOLERANCE1

      DIMENSION DFGRDM1(3,3),DFGRDM1_INV(3,3),
     1 BBAR(6),DISTGR(3,3),DCAUSDE_EM(6,6),CAUS_EM(3,3)
      DOUBLE PRECISION :: DFGRDM1,DFGRDM1_INV,BBAR,DISTGR,DCAUSDE_EM,
     1 CAUS_EM
      DOUBLE PRECISION :: DET,TRBBAR,SCALE,EG,EK,PR,EG23,EG4
      DOUBLE PRECISION :: INV1,DWab,DDWab

      DIMENSION  KIRS_V(3,3),CV_INV(3,3),BE_TRI(3,3),
     1 BET_BAR(6),BET_PS(3),BET_AN(3,3),
     2 BET_PS_LOG(3),BE_PS(3),BE_PS_LOG(3),O_BE_PS(3),
     3 KED_PS(3),KEV_PS(3),KK_PS(3,3),KK_PS_INV(3,3)
      DIMENSION DKEDDBE(3,3),DKEVDBE(3,3),RES_EXPR(3),BE_PS_INC(3),
     1 DKIRDBET(3,3),DPKSDE_V(9,9),DCAUSDE_V(9,9),
     2 BE_NEW(3,3),CV_INV_NEW(3,3),CAUS_V(3,3)

      DOUBLE PRECISION :: KIRS_V,CV_INV,BE_TRI,BET_BAR,BET_PS,BET_AN,
     1 BET_PS_LOG,BE_PS,BE_PS_LOG,KED_PS,KEV_PS,KK_PS,KK_PS_INV,O_BE_PS
      DOUBLE PRECISION :: DKEDDBE,DKEVDBE,RES_EXPR,BE_PS_INC,DKIRDBET,
     1 DPKSDE_V,DCAUSDE_V,BE_NEW,CV_INV_NEW,CAUS_V,XI

      DOUBLE PRECISION :: DET_DE,DET_DE1,DET_DE2,TRACE_BE,TERM_ONE,
     1 TRACE_BE_TOR,DET_KK,RES1,RES2,RES3,RES_MAX,
     2 DWab_BE_TOR,DDWab_BE_TOR

      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,FOUR=4.D0,
     1 COF =1.D6)

      MY_DHASH(1,1)=1
      MY_DHASH(1,2)=1
      MY_DHASH(2,1)=2
      MY_DHASH(2,2)=2
      MY_DHASH(3,1)=3
      MY_DHASH(3,2)=3
      MY_DHASH(4,1)=1
      MY_DHASH(4,2)=2
      MY_DHASH(5,1)=1
      MY_DHASH(5,2)=3
      MY_DHASH(6,1)=2
      MY_DHASH(6,2)=3
      MY_DHASH(7,1)=2
      MY_DHASH(7,2)=1
      MY_DHASH(8,1)=3
      MY_DHASH(8,2)=1
      MY_DHASH(9,1)=3
      MY_DHASH(9,2)=2 
C ---------------------------------------------------------------- 
C
C    读入材料参数
C    
C ----------------------------------------------------------------
C    PROPS(1)- Mu_alpha
C    PROPS(2)- J_alpha
C    PROPS(3)- Mu_beta
C    PROPS(4)- J_beta
C    PROPS(5)- ETA  
C ----------------------------------------------------------------
      P_Mu=PROPS(1)
      P_N=PROPS(2) 
      P_LAMDA=P_Mu*COF
      P_Mu_V=PROPS(3)
      P_N_V=PROPS(4) 
      P_LAMDA_V=P_Mu_V*COF
      ETA=PROPS(5)
      TAO_D=ETA/P_MU_V
      TAO_V=ETA/P_MU_V   
C    
      DO I=1,3
      DO J=1,3
      DFGRDM1(I,J)=DFGRD1(I,J)
      END DO
      END DO    
C
      DET=DFGRDM1(1,1)*DFGRDM1(2,2)*DFGRDM1(3,3)
     1   -DFGRDM1(1,2)*DFGRDM1(2,1)*DFGRDM1(3,3)
      IF(NSHR.EQ.3) THEN
        DET=DET+DFGRDM1(1,2)*DFGRDM1(2,3)*DFGRDM1(3,1)
     1         +DFGRDM1(1,3)*DFGRDM1(3,2)*DFGRDM1(2,1)
     2         -DFGRDM1(1,3)*DFGRDM1(3,1)*DFGRDM1(2,2)
     3         -DFGRDM1(2,3)*DFGRDM1(3,2)*DFGRDM1(1,1)
      END IF
      SCALE=DET**(-ONE/THREE)
      DO K1=1, 3
        DO K2=1, 3
          DISTGR(K2,K1)=SCALE*DFGRDM1(K2,K1)
        END DO
      END DO   
C
      BBAR(1)=DISTGR(1,1)**2+DISTGR(1,2)**2+DISTGR(1,3)**2
      BBAR(2)=DISTGR(2,1)**2+DISTGR(2,2)**2+DISTGR(2,3)**2
      BBAR(3)=DISTGR(3,3)**2+DISTGR(3,1)**2+DISTGR(3,2)**2
      BBAR(4)=DISTGR(1,1)*DISTGR(2,1)+DISTGR(1,2)*DISTGR(2,2)
     1       +DISTGR(1,3)*DISTGR(2,3)
      IF(NSHR.EQ.3) THEN
        BBAR(5)=DISTGR(1,1)*DISTGR(3,1)+DISTGR(1,2)*DISTGR(3,2)
     1         +DISTGR(1,3)*DISTGR(3,3)
        BBAR(6)=DISTGR(2,1)*DISTGR(3,1)+DISTGR(2,2)*DISTGR(3,2)
     1         +DISTGR(2,3)*DISTGR(3,3)
      END IF
C
      DFGRDM1_INV(1,1)=DFGRDM1(2,2)*DFGRDM1(3,3)
     1 -DFGRDM1(2,3)*DFGRDM1(3,2) 
      DFGRDM1_INV(1,2)=-DFGRDM1(1,2)*DFGRDM1(3,3)
     1 +DFGRDM1(1,3)*DFGRDM1(3,2)  
      DFGRDM1_INV(1,3)=DFGRDM1(1,2)*DFGRDM1(2,3)
     1 -DFGRDM1(1,3)*DFGRDM1(2,2)
      DFGRDM1_INV(2,1)=-DFGRDM1(2,1)*DFGRDM1(3,3)
     1 +DFGRDM1(2,3)*DFGRDM1(3,1) 
      DFGRDM1_INV(2,2)=DFGRDM1(1,1)*DFGRDM1(3,3)
     1 -DFGRDM1(1,3)*DFGRDM1(3,1) 
      DFGRDM1_INV(2,3)=-DFGRDM1(1,1)*DFGRDM1(2,3)
     1 +DFGRDM1(1,3)*DFGRDM1(2,1)
      DFGRDM1_INV(3,1)=DFGRDM1(2,1)*DFGRDM1(3,2)
     1 -DFGRDM1(2,2)*DFGRDM1(3,1) 
      DFGRDM1_INV(3,2)=-DFGRDM1(1,1)*DFGRDM1(3,2)
     1 +DFGRDM1(1,2)*DFGRDM1(3,1)  
      DFGRDM1_INV(3,3)=DFGRDM1(1,1)*DFGRDM1(2,2)
     1 -DFGRDM1(1,2)*DFGRDM1(2,1)
       
      DO I=1,3
      DO J=1,3
      IF (DET.NE.ZERO) THEN
      DFGRDM1_INV(I,J)=DFGRDM1_INV(I,J)/DET
      END IF
      END DO
      END DO
C
      INV1=BBAR(1)+BBAR(2)+BBAR(3)
      DWab=P_Mu*P_N/TWO/(P_N+THREE-INV1)
      DDWab=P_Mu*P_N/TWO/(P_N+THREE-INV1)**TWO
C
      TRBBAR=INV1/THREE
      EG=TWO*DWab/DET
      EK=TWO*P_LAMDA*DET 
      PR=P_LAMDA*(DET-ONE/DET)
      DO K1=1,NDI
        CAUS_EM(K1,K1)=EG*(BBAR(K1)-TRBBAR)+PR
      END DO
      CAUS_EM(1,2)=EG*BBAR(4)
      CAUS_EM(2,1)=CAUS_EM(1,2)
      IF (NSHR.EQ.3) THEN
      CAUS_EM(1,3)=EG*BBAR(5)
      CAUS_EM(3,1)=CAUS_EM(1,3)
      CAUS_EM(2,3)=EG*BBAR(6)
      CAUS_EM(3,2)=CAUS_EM(2,3)
      END IF 
C
      EG23=EG*TWO/THREE
      EG4=FOUR*DDWab/DET
      DCAUSDE_EM(1,1)= EG23*(BBAR(1)+TRBBAR)+EK
     1+EG4*((BBAR(1)-TRBBAR))**TWO
      DCAUSDE_EM(2,2)= EG23*(BBAR(2)+TRBBAR)+EK
     1+EG4*((BBAR(2)-TRBBAR))**TWO
      DCAUSDE_EM(3,3)= EG23*(BBAR(3)+TRBBAR)+EK
     1+EG4*((BBAR(3)-TRBBAR))**TWO
      DCAUSDE_EM(1,2)=-EG23*(BBAR(1)+BBAR(2)-TRBBAR)+EK
     1+EG4*(BBAR(1)-TRBBAR)*(BBAR(2)-TRBBAR)
      DCAUSDE_EM(1,3)=-EG23*(BBAR(1)+BBAR(3)-TRBBAR)+EK
     1+EG4*(BBAR(1)-TRBBAR)*(BBAR(3)-TRBBAR)
      DCAUSDE_EM(2,3)=-EG23*(BBAR(2)+BBAR(3)-TRBBAR)+EK
     1+EG4*(BBAR(2)-TRBBAR)*(BBAR(3)-TRBBAR)
      DCAUSDE_EM(1,4)= EG23*BBAR(4)/TWO
     1+EG4*(BBAR(1)-TRBBAR)*BBAR(4)
      DCAUSDE_EM(2,4)= EG23*BBAR(4)/TWO
     1+EG4*(BBAR(2)-TRBBAR)*BBAR(4)
      DCAUSDE_EM(3,4)=-EG23*BBAR(4)
     1+EG4*(BBAR(3)-TRBBAR)*BBAR(4)
      DCAUSDE_EM(4,4)= EG*(BBAR(1)+BBAR(2))/TWO
     1+EG4*BBAR(4)*BBAR(4)
      IF(NSHR.EQ.3) THEN
      DCAUSDE_EM(1,5)= EG23*BBAR(5)/TWO
     1+EG4*(BBAR(1)-TRBBAR)*BBAR(5)
      DCAUSDE_EM(2,5)=-EG23*BBAR(5)
     1+EG4*(BBAR(2)-TRBBAR)*BBAR(5)
      DCAUSDE_EM(3,5)= EG23*BBAR(5)/TWO
     1+EG4*(BBAR(3)-TRBBAR)*BBAR(5)
      DCAUSDE_EM(1,6)=-EG23*BBAR(6)
     1+EG4*(BBAR(1)-TRBBAR)*BBAR(6) 
      DCAUSDE_EM(2,6)= EG23*BBAR(6)/TWO
     1+EG4*(BBAR(2)-TRBBAR)*BBAR(6)
      DCAUSDE_EM(3,6)= EG23*BBAR(6)/TWO
     1+EG4*(BBAR(3)-TRBBAR)*BBAR(6) 
      DCAUSDE_EM(5,5)= EG*(BBAR(1)+BBAR(3))/TWO
     1+EG4*BBAR(5)*BBAR(5)
      DCAUSDE_EM(6,6)= EG*(BBAR(2)+BBAR(3))/TWO
     1+EG4*BBAR(6)*BBAR(6)
      DCAUSDE_EM(4,5)= EG*BBAR(6)/TWO
     1+EG4*BBAR(4)*BBAR(5)
      DCAUSDE_EM(4,6)= EG*BBAR(5)/TWO
     1+EG4*BBAR(4)*BBAR(6)
      DCAUSDE_EM(5,6)= EG*BBAR(4)/TWO
     1+EG4*BBAR(5)*BBAR(6)
      END IF
C
      IF ((TIME(2).EQ.ZERO).AND.(TIME(1).EQ.ZERO)) THEN
      DO K1=1,3
         STATEV(K1)=ONE
      END DO
      STATEV(4)=ZERO
      STATEV(5)=ZERO
      STATEV(6)=ZERO
      END IF
C
      DO K1=1,NDI
      CV_INV(K1,K1)=STATEV(K1)
      END DO
      CV_INV(1,2)=STATEV(4)
      CV_INV(2,1)=CV_INV(1,2)
      IF(NSHR.EQ.3) THEN
      CV_INV(1,3)=STATEV(5)
      CV_INV(3,1)=CV_INV(1,3)
      CV_INV(2,3)=STATEV(6)
      CV_INV(3,2)=CV_INV(2,3)
      END IF
C
      DO K1=1,3
      DO K2=K1,3
         BE_TRI(K1,K2)=ZERO
         DO K3=1,3
         DO K4=1,3
      BE_TRI(K1,K2)=BE_TRI(K1,K2)
     1 +DFGRDM1(K1,K3)*CV_INV(K3,K4)*DFGRDM1(K2,K4)
         END DO
         END DO
      END DO
      END DO
C
      BET_BAR(1)=BE_TRI(1,1)
      BET_BAR(2)=BE_TRI(2,2)
      BET_BAR(3)=BE_TRI(3,3)
      BET_BAR(4)=BE_TRI(1,2)
      BET_BAR(5)=BE_TRI(1,3)
      BET_BAR(6)=BE_TRI(2,3)
      CALL SPRIND(BET_BAR,BET_PS,BET_AN,1,NDI,NSHR)
C
      DO K1=1,3
         BE_PS(K1)=BET_PS(K1)
         IF (BET_PS(K1).NE.ZERO) THEN
         BET_PS_LOG(K1)=dlog(dsqrt(BET_PS(K1)))
         ELSE 
         BET_PS_LOG(K1)=ZERO
         END IF
      END DO
C
      TOLERANCE1=1.0D-8
      DO ITE_STEP=1,5000
C
      DET_DE=dsqrt(BE_PS(1)*BE_PS(2)*BE_PS(3))
      DET_DE1=DET_DE**(-ONE/THREE)
      DET_DE2=DET_DE**(-TWO/THREE)
C
      TRACE_BE=BE_PS(1)+BE_PS(2)+BE_PS(3)
      TERM_ONE=TRACE_BE/THREE
C
      TRACE_BE_TOR=DET_DE2*TRACE_BE
C
      DWab_BE_TOR=P_N_V/TWO/(P_N_V+THREE-TRACE_BE_TOR)
C
      DDWab_BE_TOR=P_N_V/TWO/(P_N_V+THREE-TRACE_BE_TOR)**TWO
C    
      DO K1=1,3
         KED_PS(K1)=TWO*DWab_BE_TOR*DET_DE2*(BE_PS(K1)-TERM_ONE)
         KEV_PS(K1)=(DET_DE**TWO-ONE)
      END DO
C
      DO K1=1,3
      DO K2=1,3
         DKEDDBE(K1,K2)=TWO*DDWab_BE_TOR*(DET_DE2**TWO)
     1 *(BE_PS(K1)-TERM_ONE)*(ONE-TERM_ONE/BE_PS(K2))
         IF(K1.EQ.K2) THEN
         DKEDDBE(K1,K2)=DKEDDBE(K1,K2)+FOUR/THREE*DWab_BE_TOR*DET_DE2
         ELSE
         DKEDDBE(K1,K2)=DKEDDBE(K1,K2)-TWO/THREE*DWab_BE_TOR*DET_DE2
         END IF
      DKEDDBE(K1,K2)=DKEDDBE(K1,K2)-TWO/THREE*DWab_BE_TOR*
     1 DET_DE2*(BE_PS(K1)-TERM_ONE)/BE_PS(K2)
      END DO
      END DO 
C
      DO K1=1,3
      DO K2=1,3
      DKEVDBE(K1,K2)=(DET_DE**(TWO))/BE_PS(K2)
      END DO
      END DO
C
      DO K1=1,3
      BE_PS_LOG(K1)=dlog(dsqrt(BE_PS(K1)))
      END DO 
C
      DO K1=1,3
      TERM1=DTIME*(KED_PS(K1)/TWO/TAO_D
     1 +KEV_PS(K1)/THREE/TAO_V)
      RES_EXPR(K1)=BE_PS_LOG(K1)+TERM1-BET_PS_LOG(K1)
      END DO
C
      DO K1=1,3
      DO K2=1,3
         KK_PS(K1,K2)=ZERO
         IF(K1.EQ.K2) THEN
         KK_PS(K1,K2)=ONE
         END IF
      KK_PS(K1,K2)=KK_PS(K1,K2)+DTIME*TWO*dexp(TWO*BE_PS_LOG(K2))
     1 *(DKEDDBE(K1,K2)/TWO/TAO_D
     2 +DKEVDBE(K1,K2)/THREE/TAO_V)
      END DO
      END DO 
C
      DET_KK=KK_PS(1,1)*KK_PS(2,2)*KK_PS(3,3)
     1 -KK_PS(1,2)*KK_PS(2,1)*KK_PS(3,3)
     2 +KK_PS(1,2)*KK_PS(2,3)*KK_PS(3,1)
     3 +KK_PS(1,3)*KK_PS(3,2)*KK_PS(2,1)
     4 -KK_PS(1,3)*KK_PS(3,1)*KK_PS(2,2)
     5 -KK_PS(2,3)*KK_PS(3,2)*KK_PS(1,1)   
      KK_PS_INV(1,1)=KK_PS(2,2)*KK_PS(3,3)
     1 -KK_PS(2,3)*KK_PS(3,2) 
C
      KK_PS_INV(1,2)=-KK_PS(1,2)*KK_PS(3,3)
     1 +KK_PS(1,3)*KK_PS(3,2)  
      KK_PS_INV(1,3)=KK_PS(1,2)*KK_PS(2,3)
     1 -KK_PS(1,3)*KK_PS(2,2)
      KK_PS_INV(2,1)=-KK_PS(2,1)*KK_PS(3,3)
     1 +KK_PS(2,3)*KK_PS(3,1) 
      KK_PS_INV(2,2)=KK_PS(1,1)*KK_PS(3,3)
     1 -KK_PS(1,3)*KK_PS(3,1) 
      KK_PS_INV(2,3)=-KK_PS(1,1)*KK_PS(2,3)
     1 +KK_PS(1,3)*KK_PS(2,1)
      KK_PS_INV(3,1)=KK_PS(2,1)*KK_PS(3,2)
     1 -KK_PS(2,2)*KK_PS(3,1) 
      KK_PS_INV(3,2)=-KK_PS(1,1)*KK_PS(3,2)
     1 +KK_PS(1,2)*KK_PS(3,1)  
      KK_PS_INV(3,3)=KK_PS(1,1)*KK_PS(2,2)
     1 -KK_PS(1,2)*KK_PS(2,1)      
      DO I=1,3
      DO J=1,3
      KK_PS_INV(I,J)=KK_PS_INV(I,J)/DET_KK
      END DO
      END DO
C
      DO K1=1,3
         BE_PS_INC(K1)=ZERO
      DO K2=1,3
         BE_PS_INC(K1)=BE_PS_INC(K1)-KK_PS_INV(K1,K2)*RES_EXPR(K2)
      END DO
      END DO
C
      DO K1=1,3
      BE_PS(K1)=dexp(TWO*(BE_PS_INC(K1)+BE_PS_LOG(K1)))
      END DO 
C
      RES1=ABS(RES_EXPR(1))
      RES2=ABS(RES_EXPR(2))
      RES3=ABS(RES_EXPR(3))
      RES_MAX=MAX(RES1,RES2,RES3)  
      IF(RES_MAX.LT.TOLERANCE1) THEN
      GOTO 10
      END IF
      END DO

      WRITE(*,20) 30000
  20  FORMAT(//,30X,'***WARNING - VISCOELASCITY ALGORITHM DID NOT '
     1 ,'CONVERGE AFTER ',I4,' ITERATIONS')
      STOP
  10  CONTINUE
C
      DO K1=1,3
         KED_PS(K1)=TWO*DWab_BE_TOR*DET_DE2*(BE_PS(K1)-TERM_ONE)
         KEV_PS(K1)=DET_DE**(TWO)-ONE
      END DO 
      DO K1=1,3
      DO K2=1,3
         KIRS_V(K1,K2)=ZERO
      IF(K1.EQ.K2) THEN
      KIRS_V(K1,K1)=P_Mu_V*KED_PS(K1)+P_LAMDA_V*KEV_PS(K1)
      END IF
      END DO
      END DO
C    
      DO K1=1,3
      DO K2=1,3
      CAUS_V(K1,K2)=ZERO
      DO K3=1,3
      DO K4=1,3
      CAUS_V(K1,K2)=CAUS_V(K1,K2)
     1 +BET_AN(K3,K1)*KIRS_V(K3,K4)*BET_AN(K4,K2)  
      END DO
      END DO
      CAUS_V(K1,K2)=CAUS_V(K1,K2)/DET
      END DO
      END DO
C
      DO K1=1,3
      DO K2=1,3
         BE_NEW(K1,K2)=ZERO
        DO K3=1,3
        BE_NEW(K1,K2)=BE_NEW(K1,K2)
     1 +BET_AN(K3,K1)*BE_PS(K3)*BET_AN(K3,K2)
        END DO
      END DO
      END DO
C
      DO K1=1,3
      DO K2=1,3
         CV_INV_NEW(K1,K2)=ZERO
         DO K3=1,3
         DO K4=1,3
      CV_INV_NEW(K1,K2)=CV_INV_NEW(K1,K2)
     1 +DFGRDM1_INV(K1,K3)*BE_NEW(K3,K4)*DFGRDM1_INV(K2,K4)
         END DO
         END DO
      END DO
      END DO

      DO K1=1,3
      DO K2=1,3
         DKIRDBET(K1,K2)=ZERO
         DO K3=1,3
      DKIRDBET(K1,K2)=DKIRDBET(K1,K2)
     1 +TWO*dexp(TWO*BE_PS_LOG(K3))*(P_Mu_V*DKEDDBE(K1,K3)
     2 +P_LAMDA_V*DKEVDBE(K1,K3))*KK_PS_INV(K3,K2)
         END DO
      END DO
      END DO
C
      DO IJ=1,9
      DO KL=1,9 
         DPKSDE_V(IJ,KL)=ZERO
         I=MY_DHASH(IJ,1)
         J=MY_DHASH(IJ,2)
         K=MY_DHASH(KL,1)
         L=MY_DHASH(KL,2)
       IF((I.EQ.J).AND.(K.EQ.L)) THEN
         DPKSDE_V(IJ,KL)=DKIRDBET(I,K)/BET_PS(I)/BET_PS(K)
         IF(I.EQ.K) THEN
       DPKSDE_V(IJ,KL)=DPKSDE_V(IJ,KL)
     1 -TWO*KIRS_V(I,I)/BET_PS(I)/BET_PS(K)
         END IF
       END IF
      IF(I.NE.J) THEN
        IF((I.EQ.K).AND.(J.EQ.L)) THEN
          IF (BET_PS(J).NE.BET_PS(I)) THEN
      DPKSDE_V(IJ,KL)=DPKSDE_V(IJ,KL)+(KIRS_V(J,J)/BET_PS(J)
     1 -KIRS_V(I,I)/BET_PS(I))/(BET_PS(J)-BET_PS(I))
          END IF
        END IF
        IF((I.EQ.L).AND.(J.EQ.K)) THEN
          IF (BET_PS(J).NE.BET_PS(I)) THEN
      DPKSDE_V(IJ,KL)=DPKSDE_V(IJ,KL)+(KIRS_V(J,J)/BET_PS(J)
     1 -KIRS_V(I,I)/BET_PS(I))/(BET_PS(J)-BET_PS(I))
          END IF
        END IF
      END IF
      END DO
      END DO
C  
      DO pq=1,6
      DO mn=pq,6
         p=MY_DHASH(pq,1)
         q=MY_DHASH(pq,2)
         m=MY_DHASH(mn,1)
         n=MY_DHASH(mn,2) 
      DCAUSDE_V(pq,mn)=ZERO
      DO IJ=1,9
      DO KL=1,9
         I=MY_DHASH(IJ,1)
         J=MY_DHASH(IJ,2)
         K=MY_DHASH(KL,1)
         L=MY_DHASH(KL,2) 
         TM=ZERO
         TM=dsqrt(BET_PS(I)*BET_PS(J)*BET_PS(K)*
     1         BET_PS(L))*DPKSDE_V(IJ,KL)
         DCAUSDE_V(pq,mn)=DCAUSDE_V(pq,mn)+
     1         BET_AN(I,p)*BET_AN(J,q)*BET_AN(K,m)*
     2         BET_AN(L,n)*TM/DET
      END DO
      END DO
         IF(p.EQ.m) THEN
          DCAUSDE_V(pq,mn)=DCAUSDE_V(pq,mn)+CAUS_V(q,n)/TWO
         END IF
         IF(p.EQ.n) THEN
          DCAUSDE_V(pq,mn)=DCAUSDE_V(pq,mn)+CAUS_V(q,m)/TWO
         END IF 
         IF(q.EQ.m) THEN
          DCAUSDE_V(pq,mn)=DCAUSDE_V(pq,mn)+CAUS_V(p,n)/TWO
         END IF
         IF(q.EQ.n) THEN
          DCAUSDE_V(pq,mn)=DCAUSDE_V(pq,mn)+CAUS_V(p,m)/TWO
         END IF
      END DO
      END DO
C
      DO K1=1,NDI
         STATEV(K1)=CV_INV_NEW(K1,K1)
      END DO
      STATEV(4)=CV_INV_NEW(1,2)
      IF(NSHR.EQ.3) THEN
      STATEV(5)=CV_INV_NEW(1,3)
      STATEV(6)=CV_INV_NEW(2,3)
      END IF
C
      DDSDDE(1,1)= DCAUSDE_EM(1,1)+DCAUSDE_V(1,1)
      DDSDDE(2,2)= DCAUSDE_EM(2,2)+DCAUSDE_V(2,2)
      DDSDDE(3,3)= DCAUSDE_EM(3,3)+DCAUSDE_V(3,3)
      DDSDDE(1,2)= DCAUSDE_EM(1,2)+DCAUSDE_V(1,2)
      DDSDDE(1,3)= DCAUSDE_EM(1,3)+DCAUSDE_V(1,3)
      DDSDDE(2,3)= DCAUSDE_EM(2,3)+DCAUSDE_V(2,3)
      DDSDDE(1,4)= DCAUSDE_EM(1,4)+DCAUSDE_V(1,4)
      DDSDDE(2,4)= DCAUSDE_EM(2,4)+DCAUSDE_V(2,4)
      DDSDDE(3,4)= DCAUSDE_EM(3,4)+DCAUSDE_V(3,4)
      DDSDDE(4,4)= DCAUSDE_EM(4,4)+DCAUSDE_V(4,4)
      IF(NSHR.EQ.3) THEN
      DDSDDE(1,5)= DCAUSDE_EM(1,5)+DCAUSDE_V(1,5)
      DDSDDE(2,5)= DCAUSDE_EM(2,5)+DCAUSDE_V(2,5)
      DDSDDE(3,5)= DCAUSDE_EM(3,5)+DCAUSDE_V(3,5)
      DDSDDE(1,6)= DCAUSDE_EM(1,6)+DCAUSDE_V(1,6)
      DDSDDE(2,6)= DCAUSDE_EM(2,6)+DCAUSDE_V(2,6)
      DDSDDE(3,6)= DCAUSDE_EM(3,6)+DCAUSDE_V(3,6)
      DDSDDE(5,5)= DCAUSDE_EM(5,5)+DCAUSDE_V(5,5)
      DDSDDE(6,6)= DCAUSDE_EM(6,6)+DCAUSDE_V(6,6)
      DDSDDE(4,5)= DCAUSDE_EM(4,5)+DCAUSDE_V(4,5)
      DDSDDE(4,6)= DCAUSDE_EM(4,6)+DCAUSDE_V(4,6)
      DDSDDE(5,6)= DCAUSDE_EM(5,6)+DCAUSDE_V(5,6)
      END IF

      DO K1=1,NTENS
        DO K2=1,K1-1
          DDSDDE(K1,K2)=DDSDDE(K2,K1)
        END DO
      END DO
C
      STRESS(1)=CAUS_EM(1,1)+CAUS_V(1,1)
      STRESS(2)=CAUS_EM(2,2)+CAUS_V(2,2)
      STRESS(3)=CAUS_EM(3,3)+CAUS_V(3,3)
      STRESS(4)=CAUS_EM(1,2)+CAUS_V(1,2)
      STRESS(5)=CAUS_EM(1,3)+CAUS_V(1,3)
      STRESS(6)=CAUS_EM(2,3)+CAUS_V(2,3)

      RETURN
      END SUBROUTINE UMAT 
