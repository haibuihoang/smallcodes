SUBROUTINE ARAT2V(A,X,Y,NX,NY,PP,XX,YY,F)
!---------------------------------------------------------------------
!
!     Computation of the interpolated functional value at (XX,YY)
!     using birational splines. The coefficient matrix A is produced
!     by ARAT2D.
!
!  INPUT
!     A        Coefficient matrix (NX,NY,4,4) produced by ARAT2D.
!     X,Y      Coordinate vectors.
!     NX,NY    Length of the above.
!     PP       Parameter for birational spline interpolation as
!              described in ARAT2D.
!     XX,YY    Coordinates of the desired point.
!
!  OUTPUT
!     F        Interpolated functional value at the desired point.
!
!---------------------------------------------------------------------
REAL A(NX,NY,4,4),X(NX),Y(NY),H(4,4),HH(4,4)

! CONSISTENCY CHECK TO MAKE SURE THE x,y POINT IS IN THE DOMAIN

IF (XX.LT.X(1)) THEN
!   PRINT*,' ARAT2V: POINT OUT OF LOWER BOUND X, RESET'
!   PRINT*,'         XX,X(1): ',XX,X(1)
   XX=X(1)
ENDIF
IF (XX.GT.X(NX)) THEN
!   PRINT*,' ARAT2V: POINT OUT OF UPPER BOUND X, RESET'
!   PRINT*,'         XX,X(NX): ',XX,X(NX)
   XX=X(NX)
ENDIF
IF (YY.LT.Y(1)) THEN
!   PRINT*,' ARAT2V: POINT OUT OF LOWER BOUND Y, RESET'
!   PRINT*,'         YY,Y(1): ',YY,Y(1)
   YY=Y(1)
ENDIF
IF (YY.GT.Y(NY)) THEN
!   PRINT*,' ARAT2V: POINT OUT OF UPPER BOUND Y, RESET'
!   PRINT*,'         YY,Y(NY): ',YY,Y(NY)
   YY=Y(NY)
ENDIF
!
!...LOCALIZE THE X INDEX OF THE POINT
!
I=1
J=1
IF (XX.LT.X(I)) GOTO 10
IF (XX.LE.X(I+1)) GOTO 40
L=NX
GOTO 30
10 L=I
I=1
20 K=(I+L)/2
IF (XX.LT.X(K)) THEN
   L=K
ELSE
   I=K
ENDIF
30 IF (L.GT.I+1) GOTO 20
!
!...LOCALIZE THE Y INDEX OF THE POINT
!
40 IF (YY.LT.Y(J)) GOTO 50
IF (YY.LE.Y(J+1)) GOTO 80
L=NY
GOTO 70
50 L=J
J=1
60 K=(J+L)/2
IF (YY.LT.Y(K)) THEN
   L=K
ELSE
   J=K
ENDIF
70 IF (L.GT.J+1) GOTO 60
!
!...COMPUTE THE COEFFICIENTS OF THE INTERPOLATION
!
80 XXX=(XX-X(I))/(X(I+1)-X(I))
YYY=(YY-Y(J))/(Y(J+1)-Y(J))
H1=1.-XXX
H2=1.-YYY
DO L=1,4
   H(1,L)=H1
   H(2,L)=XXX
   H(3,L)=(H1*H1*H1)/(PP*XXX+1.)
   H(4,L)=(XXX*XXX*XXX)/(PP*H1+1.)
   HH(L,1)=H2
   HH(L,2)=YYY
   HH(L,3)=(H2*H2*H2)/(PP*YYY+1.)
   HH(L,4)=(YYY*YYY*YYY)/(PP*H2+1.)
ENDDO
!
!...INTERPOLATE AT XX AND YY
!
F=0.
DO K=1,4
   DO L=1,4
      F=F+A(I,J,K,L)*H(K,L)*HH(K,L)
   ENDDO
ENDDO

return
END
