C
C Copyright. Author: Yanhui Jiang 
C 1/23/2018 21:32 
C NEU,Boston,02115
C
C==========================================================
C	Incorporating:Nlgeom, Hypoelasticiy 
C   Parameters: 
C    PROPS(1)- Young's modulus
C    PROPS(2)- Poisson's ratio
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
	 
      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,FOUR=4.D0,
     1 FIVE=5.D0)	 	  
 
      DIMENSION KIRS(3,3),DKIRSDE(3,3),PS1(3),AN1(3,3)
      DOUBLE PRECISION :: KIRS,DKIRSDE,DET1,PS1,AN1	  

      DOUBLE PRECISION :: YOUNG,POISSON,NU,LAMBDA  
	  
      INTEGER KZ,KX,KY
  	     
C----------------------------------------------------------------
C	
C    CALCULATE PRINCIPAL STRETCH 
C
      CALL STRTCH(DFGRD1,NDI,NSHR,DET1,PS1,AN1)	  	  
C ---------------------------------------------------------------- 
C
C    READ MATERIAL PROPERTIES 
C    
      YOUNG=PROPS(1)
      POISSON=PROPS(2)
      NU=YOUNG/TWO/(ONE+POISSON)
      LAMBDA=YOUNG*POISSON/(ONE+POISSON)/(ONE-TWO*POISSON)
C ----------------------------------------------------------------
C
C    LOCAL CALCULATION
C  	  
      KX=1
      KY=2
      KZ=3
      FORALL(I=1:3,J=1:3) KIRS(I,J)=ZERO 
         KIRS(KX,KX)= (TWO*NU*LOG(PS1(KX))+
     1 LAMBDA*LOG(DET1))*DET1
         KIRS(KY,KY)= (TWO*NU*LOG(PS1(KY))+
     1 LAMBDA*LOG(DET1))*DET1
         KIRS(KZ,KZ)= (TWO*NU*LOG(PS1(KZ))+
     1 LAMBDA*LOG(DET1))*DET1

      FORALL(I=1:3,J=1:3) DKIRSDE(I,J)=ZERO 	 
      DKIRSDE(KX,KX)= TWO*NU+LAMBDA
      DKIRSDE(KY,KY)= TWO*NU+LAMBDA
      DKIRSDE(KZ,KZ)= TWO*NU+LAMBDA
      DKIRSDE(KX,KY)= LAMBDA
      DKIRSDE(KY,KX)= LAMBDA
      DKIRSDE(KX,KZ)= LAMBDA
      DKIRSDE(KY,KZ)= LAMBDA
      DKIRSDE(KZ,KX)= LAMBDA
      DKIRSDE(KZ,KY)= LAMBDA 	 
	  
      CALL MLOTOGO(STRESS,DDSDDE,DKIRSDE,KIRS,PS1,AN1,DET1)

      RETURN		  
      END SUBROUTINE UMAT
C----------------------------------------------------------------
C
C     STRETCH AND ITS DIRECTION
C 
      SUBROUTINE STRTCH(DFGRD,NDI,NSHR,DET,PS,AN)
 
      DIMENSION DFGRD(3,3)
      DOUBLE PRECISION :: DFGRD
	  
      INTEGER :: NDI,NSHR
 
      DIMENSION RCG(6),PS(3),AN(3,3),TEMP(3,3)
 
      DOUBLE PRECISION :: RCG,PS,AN,DET,TEMP
C     THE DETERMINANT OF THE DEFORMATION GRADIENT TENSOR
      DET=DFGRD(1,1)*DFGRD(2,2)*DFGRD(3,3)
     1   -DFGRD(1,2)*DFGRD(2,1)*DFGRD(3,3)
      IF(NSHR.EQ.3) THEN
        DET=DET+DFGRD(1,2)*DFGRD(2,3)*DFGRD(3,1)
     1         +DFGRD(1,3)*DFGRD(3,2)*DFGRD(2,1)
     2         -DFGRD(1,3)*DFGRD(3,1)*DFGRD(2,2)
     3         -DFGRD(2,3)*DFGRD(3,2)*DFGRD(1,1)
      END IF  
C     THE RIGHT CAUCHY-GREEN STRAIN TENSOR
      TEMP=MATMUL(TRANSPOSE(DFGRD),DFGRD)
      RCG(1)=TEMP(1,1)
      RCG(2)=TEMP(2,2)
      RCG(3)=TEMP(3,3)
      RCG(4)=TEMP(1,2)
      RCG(5)=TEMP(1,3)
      RCG(6)=TEMP(2,3)
C     THE LOCAL COORDINATE SYS AFTER ROTATION
      CALL SPRIND(RCG,PS,AN,1,NDI,NSHR)	  
C     THE STRETCH RATIO
      FORALL(I=1:NDI) PS(I)=SQRT(PS(I))	  
      END SUBROUTINE STRTCH	  
C----------------------------------------------------------------
C
C	  UPDATED GLOBAL STRESS/DDSDDE 
C     BY USING LOCAL PRINCIPAL KIRCHHOFF STRESS/JACOBIAN
C
      SUBROUTINE MLOTOGO(STRESS,DDSDDE,DKIRSDE,KIRS,PS1,AN1,DET1)
      
      DIMENSION INDHASH(9,2)
      INTEGER INDHASH,pq,p,q,m,n,i,j,k,l,mn,IJ,KL	  
C
      DIMENSION STRESS(6),DPKSDE(9,9),KIRS(3,3),PS(3),DCAUSDE(9,9),
     1 CAUS(3,3),DKIRSDE(3,3),AN(3,3),DDSDDE(6,6),PS1(3),AN1(3,3)
      DOUBLE PRECISION DPKSDE,KIRS,PS,DCAUSDE,CAUS,TM,DKIRSDE,STRESS,AN,
     1 DET,DDSDDE,PS1,AN1,DET1
C	  
      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,FOUR=4.D0)	  
C	  
      INDHASH(1,1)=1
      INDHASH(1,2)=1
      INDHASH(2,1)=2
      INDHASH(2,2)=2
      INDHASH(3,1)=3
      INDHASH(3,2)=3
      INDHASH(4,1)=1
      INDHASH(4,2)=2
      INDHASH(5,1)=1
      INDHASH(5,2)=3
      INDHASH(6,1)=2
      INDHASH(6,2)=3
      INDHASH(7,1)=2
      INDHASH(7,2)=1
      INDHASH(8,1)=3
      INDHASH(8,2)=1
      INDHASH(9,1)=3
      INDHASH(9,2)=2 
C	  BACKUP   
      DET=DET1
      FORALL(I=1:3) PS(I)=PS1(I)**TWO
      FORALL(I=1:3,J=1:3) AN(I,J)=AN1(I,J)	  
C     UPDATE STRESS
      CAUS=MATMUL(MATMUL(TRANSPOSE(AN),KIRS),AN)	  
      FORALL (I=1:3,J=1:3) CAUS(I,J)=CAUS(I,J)/DET  
      STRESS(1)=CAUS(1,1)
      STRESS(2)=CAUS(2,2)	  
      STRESS(3)=CAUS(3,3)
      STRESS(4)=CAUS(1,2)
      STRESS(5)=CAUS(1,3)
      STRESS(6)=CAUS(2,3)
C     CALCULATE DPKSDE 
      DO IJ=1,9
      DO KL=1,9 
         DPKSDE(IJ,KL)=ZERO	 
         I=INDHASH(IJ,1)
         J=INDHASH(IJ,2)
         K=INDHASH(KL,1)
         L=INDHASH(KL,2)
       IF((I.EQ.J).AND.(K.EQ.L)) THEN
         DPKSDE(IJ,KL)=DKIRSDE(I,K)/PS(I)/PS(K)	
         IF(I.EQ.K) THEN
       DPKSDE(IJ,KL)=DPKSDE(IJ,KL)
     1 -TWO*KIRS(I,I)/PS(I)/PS(K)
         END IF
       END IF	 
      IF(I.NE.J) THEN
        IF((I.EQ.K).AND.(J.EQ.L)) THEN
          IF (PS(J).NE.PS(I)) THEN
      DPKSDE(IJ,KL)=DPKSDE(IJ,KL)+(KIRS(J,J)/PS(J)
     1 -KIRS(I,I)/PS(I))/(PS(J)-PS(I))
          ELSE
      DPKSDE(IJ,KL)=DPKSDE(IJ,KL)+TWO*DKIRSDE(J,J)
     1 -TWO*DKIRSDE(I,J)-KIRS(J,J)/PS(J)/PS(J)	  
          END IF
        END IF
        IF((I.EQ.L).AND.(J.EQ.K)) THEN
          IF (PS(J).NE.PS(I)) THEN
      DPKSDE(IJ,KL)=DPKSDE(IJ,KL)+(KIRS(J,J)/PS(J)
     1 -KIRS(I,I)/PS(I))/(PS(J)-PS(I))
	      ELSE
      DPKSDE(IJ,KL)=DPKSDE(IJ,KL)+TWO*DKIRSDE(J,J)
     1 -TWO*DKIRSDE(I,J)-KIRS(J,J)/PS(J)/PS(J)			  
	      END IF
        END IF
      END IF	
      END DO	
      END DO
C     CALCULATE DCAUSDE  
      DO pq=1,6
      DO mn=1,6
         p=INDHASH(pq,1)
         q=INDHASH(pq,2)
         m=INDHASH(mn,1)
         n=INDHASH(mn,2) 
      DCAUSDE(pq,mn)=ZERO
      DO IJ=1,9
      DO KL=1,9
         I=INDHASH(IJ,1)
         J=INDHASH(IJ,2)
         K=INDHASH(KL,1)
         L=INDHASH(KL,2) 
         TM=ZERO
         TM=dsqrt(PS(I)*PS(J)*PS(K)*
     1         PS(L))*DPKSDE(IJ,KL)
         DCAUSDE(pq,mn)=DCAUSDE(pq,mn)+
     1         AN(I,p)*AN(J,q)*AN(K,m)*
     2         AN(L,n)*TM/DET
      END DO
      END DO
          IF(p.EQ.m) THEN
           DCAUSDE(pq,mn)= DCAUSDE(pq,mn)+CAUS(q,n)/TWO
          END IF
          IF(p.EQ.n) THEN
           DCAUSDE(pq,mn)= DCAUSDE(pq,mn)+CAUS(q,m)/TWO
          END IF 
          IF(q.EQ.m) THEN
           DCAUSDE(pq,mn)= DCAUSDE(pq,mn)+CAUS(p,n)/TWO
          END IF			 
          IF(q.EQ.n) THEN
           DCAUSDE(pq,mn)= DCAUSDE(pq,mn)+CAUS(p,m)/TWO
          END IF	 
      END DO
      END DO
C     UPDATE DDSDDE
      FORALL(I=1:6,J=1:6) DDSDDE(I,J)=DCAUSDE(I,J)  		  
      END SUBROUTINE MLOTOGO	  
	  
	  
	 