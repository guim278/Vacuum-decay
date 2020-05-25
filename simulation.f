      MODULE PARAMS
      IMPLICIT NONE
      !number of points in the mesh for the R and t cordinates. DNR and DNT
      !control how often the value of the field is saved in the output file
      INTEGER, PARAMETER :: NR=2000, NT=3000,DNR=5,DNT=10
      !size of the space-time to be simulated
      DOUBLE PRECISION, PARAMETER :: RF=2000.D0, TF=3000.D0
      DOUBLE PRECISION, PARAMETER :: DR0=RF/NR, DT=TF/NT
      DOUBLE PRECISION :: DR=(RF-DR0)/(NR-1)
      DOUBLE PRECISION, DIMENSION(NR) :: GAUGE,M,RHO
      DOUBLE PRECISION, DIMENSION(NR) :: B,R,U,G,K,PHI,DPHI,X
      DOUBLE PRECISION, PARAMETER ::PI=4.D0*ATAN(1.D0)
      DOUBLE PRECISION, PARAMETER :: PHI0=0.029D0
      DOUBLE PRECISION :: DPHI0=-3.2d-4
      DOUBLE PRECISION :: L=520.5D0
      DOUBLE PRECISION :: T=0

      CONTAINS

      FUNCTION V0(XIN) RESULT(V1)
      !Potential to be used for the evolution of the scalar field      
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: XIN
      DOUBLE PRECISION, DIMENSION(SIZE(PHI,1)) :: V1
          V1=1.D-8*((100*XIN)**2-(100*XIN)**3+0.2D0*(100*XIN)**4)+5.D-7
      END FUNCTION V0

      FUNCTION DV0(XIN) RESULT(DV1)
      !Derivative of the potential
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: XIN
      DOUBLE PRECISION, DIMENSION(SIZE(PHI,1)) :: DV1
          DV1=1.D-6*(2*(100*XIN)-3*(100*XIN)**2+0.8D0*(100*XIN)**3)
      END FUNCTION DV0

      FUNCTION DER(F) RESULT(D)
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: F
      DOUBLE PRECISION, DIMENSION(SIZE(F,1)) :: D
      INTEGER I,J
        D(1)=(F(2)-F(1))/DR
        J=SIZE(F,1)
        DO I=2,J-1
            D(I)=(F(I+1)-F(I-1))/(2*DR)
        ENDDO
        D(J)=(F(J)-F(J-1))/DR
      END FUNCTION DER

      FUNCTION DER2(F) RESULT(D)
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: F
      DOUBLE PRECISION, DIMENSION(SIZE(F,1)) :: D
      INTEGER I,J
        D(1)=(F(3)-2*F(2)+F(1))/DR**2
        J=SIZE(F,1)
        DO I=2,J-1
              D(I)=(F(I+1)-2*F(I)+F(I-1))/(DR**2)
        ENDDO
        D(J)=(F(J)-2*F(J-1)+F(J-2))/DR**2
      END FUNCTION DER2

      END MODULE

      PROGRAM MAIN
      USE PARAMS
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(7*NR) :: Y0,Y1
      DOUBLE PRECISION, DIMENSION(NR) :: DB,DP,D2P
      DOUBLE PRECISION T0,T1
      INTEGER I,J
      CALL CPU_TIME(T0)
      CALL INIT(Y0)
      !name of the output file
      OPEN(11,FILE='data.dat')
      
      DO I=1,NT
          CALL RK4(Y0,Y1)

          T=T+DT
          B=Y1(1:NR)
          R=Y1(NR+1:2*NR)
          U=Y1(2*NR+1:3*NR)
          G=Y1(3*NR+1:4*NR)
          K=Y1(4*NR+1:5*NR)
          PHI=Y1(5*NR+1:6*NR)
          DPHI=Y1(6*NR+1:7*NR)
          DP=DER(PHI)

          CALL ADAPT()
          
          Y0=(/B,R,U,G,K,PHI,DPHI/)
          
          Y0(6*NR)=Y0(6*NR-1)
          Y0(5*NR+1)=Y0(5*NR+2)
          
          IF (MOD(I,DNT).EQ.0) THEN
          CALL CHECK()
          DO J=1,NR,DNR
          WRITE(11,*) X(J),T,B(J),R(J),U(J),G(J),K(J),PHI(J),DPHI(J),
     &                GAUGE(J),DP(J)

          ENDDO
          WRITE(11,*)
          CALL CPU_TIME(T1)
          PRINT*, T,T1-T0,X(1),DR
          ENDIF
      ENDDO



      END PROGRAM

     !------------------------------------------------------------------
     !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
     !------------------------------------------------------------------
      SUBROUTINE RK4(Y0,Y1)

      USE PARAMS
      IMPLICIT NONE
      DOUBLE PRECISION Y0(7*NR),Y1(7*NR),KRK(4,7*NR)
    
      CALL F(Y0,KRK(1,:))
      KRK(1,:)=DT*KRK(1,:)

      CALL F(Y0+0.5D0*KRK(1,:),KRK(2,:))
      KRK(2,:)=DT*KRK(2,:)

      CALL F(Y0+0.5D0*KRK(2,:),KRK(3,:))
      KRK(3,:)=DT*KRK(3,:)

      CALL F(Y0+KRK(3,:),KRK(4,:))
      KRK(4,:)=DT*KRK(4,:)

C     Amb les k calculem yyout

      Y1=Y0+(KRK(1,:)+2*KRK(2,:)+2*KRK(3,:)+KRK(4,:))/6

      END SUBROUTINE RK4

      SUBROUTINE F(Y,DY)
          USE PARAMS
          IMPLICIT NONE
          INTEGER I
          DOUBLE PRECISION, DIMENSION(NR) :: BT,RT,UT,GT,KT,PHIT,DPHIT
          DOUBLE PRECISION, DIMENSION(NR) :: DRPT,D2RPT,DRBT,V,DV
          DOUBLE PRECISION, DIMENSION(7*NR) :: Y,DY

          BT=Y(1:NR)
          RT=Y(NR+1:2*NR)
          UT=Y(2*NR+1:3*NR)
          GT=Y(3*NR+1:4*NR)
          KT=Y(4*NR+1:5*NR)
          PHIT=Y(5*NR+1:6*NR)
          DPHIT=Y(6*NR+1:7*NR)
          DRPT=DER(PHIT)
          D2RPT=DER2(PHIT)
          DRBT=DER(B)
          V=V0(PHIT)
          DV=DV0(PHIT)
          
          DY(1:NR)=BT*(KT-2*UT/RT)
          DY(NR+1:2*NR)=UT
          DY(2*NR+1:3*NR)=-(1-GT**2+UT**2)/(2*RT)
     &                       -4*PI*RT*(DPHIT**2/2
     &                        +DRPT**2/(2*BT**2) - V)
          DY(3*NR+1:4*NR)=4*PI*RT*DPHIT*DRPT/BT
          DY(4*NR+1:5*NR)=-(KT-2*UT/RT)**2-2*(UT/RT)**2-4*PI*
     &                        (2*DPHIT**2-2*V)
          DY(5*NR+1:6*NR)=DPHIT
          DY(6*NR+1:7*NR)=-KT*DPHIT+D2RPT/BT**2+2*DRPT*GT/(BT*RT)
     &                      -DRPT*DRBT/BT**3-DV 

      END SUBROUTINE F


      SUBROUTINE INIT(Y)
        USE PARAMS
        IMPLICIT NONE
        DOUBLE PRECISION, DIMENSION(7*NR) :: Y
        INTEGER I,J

      DO I=1,NR
        X(I)=DR*(I-1)+DR0
      ENDDO
    
      B(:)=1.D0

      R(:)=X(:)
      
      G(:)=1.D0      
      
      PHI(:)=PHI0
      
      DPHI(:)=DPHI0*EXP(-X(:)**2/(2*L**2))
      
      RHO(:)=DPHI**2/2*(1+X(:)**2/L**4)+V0(PHI)

      M(:)=4*PI*V0(PHI)*X(:)**3/3+
     &     DPHI0**2/4*SQRT(PI**3)*L*(2*L**2+3)*ERF(X(:)/L) -
     &     DPHI0**2/(2*L**2)*PI*X(:)*(2*L**4+3*L**2+2*X(:)**2)*
     &     EXP(-(X(:)/L)**2)

      U(:)=SQRT(2*M(:)/X(:))

      K(:)=2*PI*RHO(:)*X(:)**(1.5D0)/SQRT(M(:)/2)+
     &             3*SQRT(M(:)/2)*X(:)**(-1.5D0)

      Y=(/B,R,U,G,K,PHI,DPHI/)

      END SUBROUTINE INIT
      

      SUBROUTINE CHECK()
      USE PARAMS
      IMPLICIT NONE 
      DOUBLE PRECISION, DIMENSION(NR) :: V, DV,DRP
      INTEGER I

      DRP=DER(PHI)

      RHO(:)=DPHI**2/2+(DRP/B)**2/2+V0(PHI)

      M=R/2*(1-G**2+U**2)

      GAUGE=ABS(DER(M)-4*PI*R**2*(G*B*RHO+U*DPHI*DRP))
      END SUBROUTINE CHECK

      SUBROUTINE ADAPT()
          USE PARAMS
          IMPLICIT NONE 
          DOUBLE PRECISION,DIMENSION(NR) :: BA,RA,UA,GA,KA,PHIA,DPHIA,XA
          DOUBLE PRECISION DRA
          INTEGER I
C           DRA=RF/R(NR)*DR
C           XA=RF/R(NR)*X
          DRA=(RF-DR0)/(NR-1)*X(NR)/R(NR)
          DO I=1,NR
            XA(I)=DRA*(I-1)+DR0/R(NR)*X(NR)        
          ENDDO
          CALL INTERPOLATE(B,XA,BA)
          CALL INTERPOLATE(R,XA,RA)
          CALL INTERPOLATE(U,XA,UA)
          CALL INTERPOLATE(G,XA,GA)
          CALL INTERPOLATE(K,XA,KA)
          CALL INTERPOLATE(PHI,XA,PHIA)
          CALL INTERPOLATE(DPHI,XA,DPHIA)
          DR=DRA
          X=XA
          B=BA
          R=RA
          U=UA
          G=GA
          K=KA
          PHI=PHIA
          DPHI=DPHIA
      END SUBROUTINE ADAPT

      SUBROUTINE INTERPOLATE(F,XA,FA)
      USE PARAMS
      DOUBLE PRECISION, DIMENSION(NR) :: F,XA,FA
      INTEGER MIN,I
      LOGICAL KYS
        DO I=1,NR
        MIN=INT(XA(I)/DR)
        IF (MIN.GT.0) THEN
        FA(I)=(F(MIN+1)-F(MIN))/(X(MIN+1)-X(MIN))*(XA(I)-X(MIN))+F(MIN)
        ELSE
        FA(I)=(F(2)-F(1))/(X(2)-X(1))*(XA(I)-X(1))+F(1)
        ENDIF

        IF (ISNAN(FA(I))) THEN
          PRINT*,'KYS'
          PRINT*,XA(I),DR,MIN
          PRINT*,B(MIN),R(MIN),G(MIN)
          STOP
        ENDIF

        ENDDO
      END SUBROUTINE INTERPOLATE
