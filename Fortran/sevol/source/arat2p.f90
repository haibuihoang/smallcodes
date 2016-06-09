SUBROUTINE ARAT2P(U,MX,MY,P,RX,NNX,NNY,IVJ,IPU,P3,AX,BX,CX,DX,NX,NY)
!---------------------------------------------------------------------
!
!     Solves the equation system for the subroutine ARAT2D.
!
!  INPUT
!     U          First matrix (MX,MY).
!     MX,MY      Physical length of the above.
!     P          Matrix (NNX,NNY) of last sweep.
!     NNX,NNY    Physical length of P.
!     IVJ,IPU    Counting parameters.
!     P3         Parameter of the birational interpolation.
!     AX,BX,CX   Diagonals of L.
!     DX         Grid vector.
!     NX,NY      Length of the calculation.
!
!  OUTPUT
!     P          Second matrix (NNX,NNY).
!     RX         Solution vector (NNX).
!
!---------------------------------------------------------------------
REAL U(MX,MY),P(NNX,NNY),AX(NX),BX(NX),CX(NX),RX(NX),DX(NX)

NX2=NX-2
IF (IVJ.NE.0) THEN
   DO J=1,NY,IPU
      R2=P3*DX(1)*DX(1)*(U(J,2)-U(J,1))
      R1=R2
      R2=P3*DX(2)*DX(2)*(U(J,3)-U(J,2))
      RX(1)=R1+R2
      RX(1)=RX(1)-DX(1)*P(J,1)
      R1=R2
      DO I=3,NX-2
         I1=I-1
         R2=P3*DX(I)*DX(I)*(U(J,I+1)-U(J,I))
         RX(I1)=R1+R2
         R1=R2
      ENDDO
      R2=P3*DX(NX-1)*DX(NX-1)*(U(J,NX)-U(J,NX-1))
      RX(NX2)=R1+R2
      RX(NX2)=RX(NX2)-DX(NX-1)*P(J,NX)
      RX(1)=RX(1)/BX(1)
      DO K=2,NX2
         KM1=K-1
         RX(K)=(RX(K)-AX(KM1)*RX(KM1))/BX(K)
      ENDDO
      DO K=NX2-1,1,-1
         RX(K)=RX(K)-CX(K)*RX(K+1)
      ENDDO
      DO I=2,NX-1
         P(J,I)=RX(I-1)
      ENDDO
   ENDDO
ELSE
   DO J=1,NY,IPU
      R2=P3*DX(1)*DX(1)*(U(2,J)-U(1,J))
      R1=R2
      R2=P3*DX(2)*DX(2)*(U(3,J)-U(2,J))
      RX(1)=R1+R2
      RX(1)=RX(1)-DX(1)*P(1,J)
      R1=R2
      DO I=3,NX-2
         I1=I-1
         R2=P3*DX(I)*DX(I)*(U(I+1,J)-U(I,J))
         RX(I1)=R1+R2
         R1=R2
      ENDDO
      R2=P3*DX(NX-1)*DX(NX-1)*(U(NX,J)-U(NX-1,J))
      RX(NX2)=R1+R2
      RX(NX2)=RX(NX2)-DX(NX-1)*P(NX,J)
      RX(1)=RX(1)/BX(1)
      DO K=2,NX2
         KM1=K-1
         RX(K)=(RX(K)-AX(KM1)*RX(KM1))/BX(K)
      ENDDO
      DO K=NX2-1,1,-1
         RX(K)=RX(K)-CX(K)*RX(K+1)
      ENDDO
      DO I=2,NX-1
         P(I,J)=RX(I-1)
      ENDDO
   ENDDO
ENDIF

return
END
