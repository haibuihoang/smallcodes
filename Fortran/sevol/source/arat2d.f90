SUBROUTINE ARAT2D(U,X,Y,NX,NY,PP,A)
!---------------------------------------------------------------------
!
!     Computes for given functional values U(I,J) at points X(I) and
!     Y(J) the coefficients A(I,J,K,L) of a birational spline
!     interpolation. The boundary points are determined through a
!     simple difference scheme.
!
!     THIS IS THE FIRST SUBROUTINE IN THE SEQUENCE; CALLS ARATPE.
!
!  INPUT
!     U         Gridded data matrix (NX,NY).
!     X,Y       Coordinates of the data matrix.
!     NX,NY     Length of the above.
!     PP        Parameter of the birational spline interpolation.
!               NOTE THAT (PP > -1).
!               PP = 0:             Bicubic Spline interpolation.
!               PP --> + Infinity:  Towards bilinear interpolation.
!                                   Highest value used so far: 10.
!
!  OUTPUT
!     A         Coefficient matrix (NX,NY,4,4) of the birational
!               spline interpolation.
!
!---------------------------------------------------------------------
PARAMETER (EPS=1.E-10,NB=4)

REAL U(NX,NY),X(NX),Y(NY),A(NX,NY,4,4)
REAL P(NX,NY),Q(NX,NY),R(NX,NY)
REAL AX(NX),BX(NX),CX(NX),RX(NX)
REAL AY(NY),BY(NY),CY(NY),RY(NY)
REAL DX(NX),DY(NY)
REAL B(NB,NB),C(NB,NB),D(NB,NB),E(NB,NB)
!
!...CHECKS
!
IF ((NX.LT.2).OR.(NY.LT.2)) THEN
   PRINT*,' ARAT2D: LENGTH ERROR, FULL STOP'
   STOP
ENDIF
IF (PP.LT.EPS-1.) THEN
   PRINT*,' ARAT2D: ERROR (PP < -1), FULL STOP'
   STOP
ENDIF
!
!...PARAMETERS AND NONEQUIDISTANT GRID
!
NX1=NX-1
NX2=NX-2
NY1=NY-1
NY2=NY-2
II=MOD(NX1,4)
DO I=1,II
   DX(I)=1./(X(I+1)-X(I))
ENDDO
DO I=1+II,NX1,4
   DX(I)=1./(X(I+1)-X(I))
   DX(I+1)=1./(X(I+2)-X(I+1))
   DX(I+2)=1./(X(I+3)-X(I+2))
   DX(I+3)=1./(X(I+4)-X(I+3))
ENDDO
JJ=MOD(NY1,4)
DO J=1,JJ
   DY(J)=1./(Y(J+1)-Y(J))
ENDDO
DO J=1+JJ,NY1,4
   DY(J)=1./(Y(J+1)-Y(J))
   DY(J+1)=1./(Y(J+2)-Y(J+1))
   DY(J+2)=1./(Y(J+3)-Y(J+2))
   DY(J+3)=1./(Y(J+4)-Y(J+3))
ENDDO
P1=PP+1.
P2=PP+2.
P3=PP+3.
DK=1./P1
B(1,1)=P2*DK
B(2,1)=-DK
B(4,3)=-DK
B(4,1)=DK
B(2,3)=B(1,1)
B(1,3)=-DK
B(3,1)=-DK
B(3,3)=DK
E(1,1)=B(1,1)
E(2,1)=B(2,1)
E(3,1)=B(3,1)
E(4,1)=B(4,1)
E(1,3)=B(1,3)
E(2,3)=B(2,3)
E(3,3)=B(3,3)
E(4,3)=B(4,3)
!
!...COMPUTE PARTIAL DERIVATIVES AT THE BOUNDARIES
!
JJ=MOD(NY,4)
DO J=1,JJ
   P(1,J)=(U(2,J)-U(1,J))*DX(1)
   P(NX,J)=(U(NX,J)-U(NX1,J))*DX(NX1)
ENDDO
TEMP=DX(1)
DO J=1+JJ,NY,4
   P(1,J)=(U(2,J)-U(1,J))*TEMP
   P(1,J+1)=(U(2,J+1)-U(1,J+1))*TEMP
   P(1,J+2)=(U(2,J+2)-U(1,J+2))*TEMP
   P(1,J+3)=(U(2,J+3)-U(1,J+3))*TEMP
ENDDO
TEMP=DX(NX1)
DO J=1+JJ,NY,4
   P(NX,J)=(U(NX,J)-U(NX1,J))*TEMP
   P(NX,J+1)=(U(NX,J+1)-U(NX1,J+1))*TEMP
   P(NX,J+2)=(U(NX,J+2)-U(NX1,J+2))*TEMP
   P(NX,J+3)=(U(NX,J+3)-U(NX1,J+3))*TEMP
ENDDO
II=MOD(NX,4)
DO I=1,II
   Q(I,1)=(U(I,2)-U(I,1))*DY(1)
   Q(I,NY)=(U(I,NY)-U(I,NY1))*DY(NY1)
ENDDO
TEMP=DY(1)
DO I=1+II,NX,4
   Q(I,1)=(U(I,2)-U(I,1))*TEMP
   Q(I+1,1)=(U(I+1,2)-U(I+1,1))*TEMP
   Q(I+2,1)=(U(I+2,2)-U(I+2,1))*TEMP
   Q(I+3,1)=(U(I+3,2)-U(I+3,1))*TEMP
ENDDO
TEMP=DY(NY1)
DO I=1+II,NX,4
   Q(I,NY)=(U(I,NY)-U(I,NY1))*TEMP
   Q(I+1,NY)=(U(I+1,NY)-U(I+1,NY1))*TEMP
   Q(I+2,NY)=(U(I+2,NY)-U(I+2,NY1))*TEMP
   Q(I+3,NY)=(U(I+3,NY)-U(I+3,NY1))*TEMP
ENDDO
R(1,1)=.5*((P(1,2)-P(1,1))*DY(1)+(Q(2,1)-Q(1,1))*DX(1))
R(1,NY)=.5*((P(1,NY)-P(1,NY1))*DY(NY1)+(Q(2,NY)-Q(1,NY))*DX(1))
R(NX,1)=.5*((P(NX,2)-P(NX,1))*DY(1)+(Q(NX,1)-Q(NX1,1))*DX(NX1))
R(NX,NY)=.5*((P(NX,NY)-P(NX,NY1))*DY(NY1)+(Q(NX,NY)-Q(NX1,NY))*DX(NX1))
!
!...PREPARE THE COEFFICIENT MATRIX FOR X
!
70 IF ((NX.EQ.2).AND.(NY.EQ.2)) GOTO 130
IF (NX.EQ.2) GOTO 90
II=MOD(NX2-1,4)
DO I=1,II
   AX(I)=DX(I+1)
   BX(I)=P2*(DX(I+1)+DX(I))
   IF (ABS(BX(I)).LT.EPS) GOTO 82
ENDDO
DO I=1+II,NX2-1,4
   AX(I)=DX(I+1)
   AX(I+1)=DX(I+2)
   AX(I+2)=DX(I+3)
   AX(I+3)=DX(I+4)
ENDDO
DO I=1+II,NX2-1,4
   BX(I)=P2*(DX(I+1)+DX(I))
   BX(I+1)=P2*(DX(I+2)+DX(I+1))
   BX(I+2)=P2*(DX(I+3)+DX(I+2))
   BX(I+3)=P2*(DX(I+4)+DX(I+3))
   IF ((ABS(BX(I)).LT.EPS).OR.(ABS(BX(I+1)).LT.EPS).OR.&
       (ABS(BX(I+2)).LT.EPS).OR.(ABS(BX(I+3)).LT.EPS)) GOTO 82
ENDDO
GOTO 83
82 PRINT*,' ARAT2D: MAIN DIAGONAL ELEMENT OF LU TOO SMALL, FULL STOP'
STOP
83 CONTINUE
BX(NX2)=P2*(DX(NX2+1)+DX(NX2))
!
!...LU- DECOMPOSITION FOR X
!
DO K=2,NX2
   KM1=K-1
   H=AX(KM1)/BX(KM1)
   CX(KM1)=H
   BX(K)=BX(K)-H*AX(KM1)
ENDDO
IVJ=0
IPU=1
CALL ARAT2P(U,NX,NY,P,RX,NX,NY,IVJ,IPU,P3,AX,BX,CX,DX,NX,NY)
90 IF (NY.EQ.2) GOTO 110
!
!...PREPARE THE COEFFICIENT MATRIX FOR Y
!
JJ=MOD(NY2-1,4)
DO J=1,JJ
   AY(J)=DY(J+1)
   BY(J)=P2*(DY(J+1)+DY(J))
   IF (ABS(BY(J)).LT.EPS) GOTO 102
ENDDO
DO J=1+JJ,NY2-1,4
   AY(J)=DY(J+1)
   AY(J+1)=DY(J+2)
   AY(J+2)=DY(J+3)
   AY(J+3)=DY(J+4)
ENDDO
DO J=1+JJ,NY2-1,4
   BY(J)=P2*(DY(J+1)+DY(J))
   BY(J+1)=P2*(DY(J+2)+DY(J+1))
   BY(J+2)=P2*(DY(J+3)+DY(J+2))
   BY(J+3)=P2*(DY(J+4)+DY(J+3))
   IF ((ABS(BY(J)).LT.EPS).OR.(ABS(BY(J+1)).LT.EPS).OR.&
       (ABS(BY(J+2)).LT.EPS).OR.(ABS(BY(J+3)).LT.EPS)) GOTO 102
ENDDO
GOTO 103
102 PRINT*,' ARAT2D: MAIN DIAGONAL ELEMENT OF LU TOO SMALL, FULL STOP'
STOP
103 CONTINUE
BY(NY2)=P2*(DY(NY2+1)+DY(NY2))
!
!...LU- DECOMPOSITION FOR Y
!
DO K=2,NY2
   KM1=K-1
   H=AY(KM1)/BY(KM1)
   CY(KM1)=H
   BY(K)=BY(K)-H*AY(KM1)
ENDDO
IVJ=1
IPU=1
CALL ARAT2P(U,NX,NY,Q,RY,NX,NY,IVJ,IPU,P3,AY,BY,CY,DY,NY,NX)
110 IF (NX.EQ.2) GOTO 120
IVJ=0
IPU=NY1
CALL ARAT2P(Q,NX,NY,R,RX,NX,NY,IVJ,IPU,P3,AX,BX,CX,DX,NX,NY)
120 IF (NY.EQ.2) GOTO 130
IVJ=1
IPU=1
CALL ARAT2P(P,NX,NY,R,RY,NX,NY,IVJ,IPU,P3,AY,BY,CY,DY,NY,NX)
!
!...COMPUTE THE COEFFICIENT MATRIX
!
130 DO I=1,NX1
   I1=I+1
   FA=(X(I1)-X(I))/(PP+3.)
   B(1,2)=B(1,1)*FA
   B(2,2)=B(2,1)*FA
   B(3,2)=-B(1,2)
   B(4,2)=-B(2,2)
   B(1,4)=B(4,2)
   B(2,4)=B(2,2)*(PP+2.)
   B(3,4)=B(2,2)
   B(4,4)=-B(2,4)
   DO J=1,NY1
      J1=J+1
      C(1,1)=U(I,J)
      C(1,2)=Q(I,J)
      C(2,1)=P(I,J)
      C(2,2)=R(I,J)
      C(1,3)=U(I,J1)
      C(1,4)=Q(I,J1)
      C(2,3)=P(I,J1)
      C(2,4)=R(I,J1)
      C(3,1)=U(I1,J)
      C(3,2)=Q(I1,J)
      C(4,1)=P(I1,J)
      C(4,2)=R(I1,J)
      C(3,3)=U(I1,J1)
      C(3,4)=Q(I1,J1)
      C(4,3)=P(I1,J1)
      C(4,4)=R(I1,J1)
      DO K1=1,NB
         DO K2=1,NB
            SUMA=0.
            DO K=1,NB
               SUMA=SUMA+B(K1,K)*C(K,K2)
            ENDDO
            D(K1,K2)=SUMA
         ENDDO
      ENDDO
      FB=(Y(J1)-Y(J))/(PP+3.)
      E(1,2)=E(1,1)*FB
      E(2,2)=E(2,1)*FB
      E(3,2)=-E(1,2)
      E(4,2)=-E(2,2)
      E(1,4)=E(4,2)
      E(2,4)=E(2,2)*(PP+2.)
      E(3,4)=E(2,2)
      E(4,4)=-E(2,4)
      DO K1=1,NB
         DO K2=1,NB
            SUMA=0.
            DO K=1,NB
               SUMA=SUMA+D(K1,K)*E(K2,K)
            ENDDO
            A(I,J,K1,K2)=SUMA
         ENDDO
      ENDDO
   ENDDO
ENDDO

return
END
