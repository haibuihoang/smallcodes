!==============================================================================
!																					
!--This module will solve the Sawyer-Eliassen(SE) with Successive Overrelaxation(SOR) method	(Press et.al. 1992)
!																				    
!==============================================================================

MODULE solve

USE data
USE COMM

IMPLICIT NONE

!--Define the terms in the SE equation in elliptic form
DOUBLE PRECISION, dimension(Nr,Nz) :: A,B,C,D,E,F 

!-- Other variables
DOUBLE PRECISION, dimension(Nr,Nz) :: DvDt, Diffusion



!--The gridpoint that violates the eliptic condition the most
integer :: ire,ize
integer :: irec,istep,Time
integer, parameter :: ifile=101

CONTAINS

!--
subroutine solve_init   
   open(ifile,file='evol.dat',access='direct',form='unformatted',recl=Nr*Nz*4)
   irec=1
   istep=1   
   Time=0
end subroutine solve_init
subroutine solve_deinit
   close(ifile)
end subroutine solve_deinit

!--This subroutine predict value for newstep
subroutine new_step
    DOUBLE PRECISION, dimension(Nr,Nz) :: M, LHS, RHS,PR	 ! absolute angular momentum per unit mass
    DOUBLE PRECISION :: Src_mag_c,Rtempx,Rtempz,Src_Ts
    integer :: ir,iz

!0- Update heating source gradually over time
!-- Here we use Potential Radius (PR) defined as
!-- 1/2 PR^2 = fv + 1/2 r^2

    Src_Ts=Src_T*60

    !Src_Speed=2.8  ! The Spead speed of src 0.28m/s ~10km/h
    !test try spread the heating source overtime 
    !to advoid the collapse of it
    Src_sx = Src_sx + Src_speed*DT 
    Src_r = Src_r + Src_speed*DT/2 

    
    do ir=1,Nr
       do iz=1,Nz
          ! Potential radius
          PR(ir,iz) = 2./fc * (R(ir)*V(ir,iz) + 0.5*fc*R(ir)*R(ir) ) 
	  if (PR(ir,iz) .lt. 0 )then
	     open (113,file="out.txt")
	     write(113,*) "Src_speed=",Src_speed
	     write(113,*) "Kr=",Kr
	     write(113,*) "Cd=",Cd
	     write(113,*) "Stop because Pr^2<0"
	     write(113,*) "Pr^2 = ",PR(ir,iz)
	     write(113,*) "At point:", ir, iz
	     write(113,*) "R,Z:", R(ir), Z(iz)
	     write(113,*) "istep, Time:",istep, Time	     
	     close(113)
	     stop "Stop becase Pr^2<0"
	  endif
	  PR(ir,iz) = SQRT(PR(ir,iz))
	  
!	  Update DIABATIC HEATING SOURCE AT (Src_r,Rrc_z)
 	  Rtempx = ABS(PR(ir,iz)-Src_r)
  	  Rtempz = ABS(Z(iz)-Src_z)
	  if (Time < Src_Ts)	then
	     Src_mag_c = Src_mag*Time/Src_Ts
	  else
	     Src_mag_c = Src_mag
	  endif
	  if ((Rtempx<=Src_sx) .and. (Rtempz<=Src_sz)) then
	     Thetadot(ir,iz)=Src_mag_c*cos(pi*(Rtempx/Src_sx)/2)*cos(pi*(Rtempz/Src_sz)/2)
	  endif
       enddo
    enddo

!1- Calculate terms for the SE equation
    call calc_terms
!2- Solve the se equation
    CALL sor(Psi,A,B,C,D,E,F,nr,nz,Dr,Dz)
!3- Compute u,w from stream function
    CALL psi2uw
!4- Calculate new v from the tendency equation

!-- Plot grads
    write(ifile,rec=irec) real (U_se)
	irec=irec+1
    write(ifile,rec=irec) real (V)
	irec=irec+1
    write(ifile,rec=irec) real (W_se)
	irec=irec+1
    write(ifile,rec=irec) real (Theta)
	irec=irec+1
    write(ifile,rec=irec) real (Rho)
	irec=irec+1
    write(ifile,rec=irec) real (P)
	irec=irec+1
    write(ifile,rec=irec) real (Vdot)
	irec=irec+1
    write(ifile,rec=irec) real (Thetadot)
	irec=irec+1
    write(ifile,rec=irec) real (Psi)
	irec=irec+1
    write(ifile,rec=irec) real (Delta)
	irec=irec+1
    write(ifile,rec=irec) real (DvDt)
	irec=irec+1
    write(ifile,rec=irec) real (Diffusion)
	irec=irec+1
    
    !--Absolute Angular Momentum
    do ir=1,Nr
       do iz=1,Nz
          M(ir,iz) = R(ir)*V(ir,iz) + 0.5*fc*R(ir)*R(ir) 
       enddo
    enddo
    write(ifile,rec=irec) real (M)
    irec=irec+1

    !--20100112-Calc LHS and RHS of the thermal wind eq. for checking
    LHS=0.
    RHS=0.
    do ir=1+1,Nr-1
       do iz=1+1,Nz-1
          RHS(ir,iz) = -1/(g*2*Dz)*( V(ir,iz+1)*V(ir,iz+1)/R(ir) + fc*V(ir,iz+1)  & 
	                             -V(ir,iz-1)*V(ir,iz-1)/R(ir) - fc*V(ir,iz-1)    )
          LHS(ir,iz) = 1/(2*Dr)*( log(Rho(ir+1,iz)) - log(Rho(ir-1,iz)) ) &
	              + (V(ir,iz)*V(ir,iz)/R(ir) + fc*V(ir,iz))/g  &
		        *1/(2*Dz)*( log(Rho(ir,iz+1)) - log(Rho(ir,iz-1)) ) 
       enddo
    enddo
    
    write(ifile,rec=irec) real (LHS)
	irec=irec+1
    write(ifile,rec=irec) real (RHS)
	irec=irec+1
    write(ifile,rec=irec) real (PR)
	irec=irec+1
    

    CALL calc_tendency
    V(:,:) = V(:,:) + DvDt(:,:)*DT

    write(ifile,rec=irec) real (V)
    irec=irec+1

!5- Recaculate balanced theta, p and density 
    CALL init_balance_data
!6- Update Friction force
    if (friction) then
 	write(*,*)"Initilizing the friction..."
	do iz=1,Nz
            do ir=1,Nr
		Vdot(ir,iz) = - Cd*(.9*V(ir,1))**2.*exp( -1.*(Z(iz)/z0)**2 ) / H
	    enddo
	enddo
    endif

!	pause
    istep=istep+1
    Time=Time+DT
    write(*,*)"---->",istep,irec,Time
!	pause

end subroutine new_step



!-- This subroutine calculate tendency of tangential wind
subroutine calc_tendency()
	integer :: i,k
	DOUBLE PRECISION, DIMENSION(Nr,Nz) :: DvDr,DvDz
!	DOUBLE PRECISION,Parameter :: Kr=0.  ! Radial Diffusion Coef, Defautl=1.e4 
 	do i=1,Nr
	   do k=1,Nz
	      if (i == Nr) then
	         DvDr(i,k) = (V(i,k)-V(i-2,k))/(2.*Dr)
	      elseif(i==1) then
	         DvDr(i,k) = (V(i+2,k)-V(i,k))/(2.*Dr)
	      else
	         DvDr(i,k) = (V(i+1,k)-V(i-1,k))/(2.*Dr)
	      endif
	    
	      if (k == 1) then
	         DvDz(i,k) = (V(i,k+1)-V(i,k))/Dz
	      elseif (k == Nz) then
	         DvDz(i,k) = (V(i,k)-V(i,k-1))/Dz
	      else
	         DvDz(i,k) = (V(i,k+1)-V(i,k-1))/(2.*Dz)
	      endif
	   enddo
	enddo
	
	do i=2,Nr
	   do k=1,Nz
	      if (i == Nr) then
	         Diffusion(i,k) = 1./R(i) * ( R(i)*Kr*DvDr(i,k) - R(i-1)*Kr*DvDr(i-1,k))/Dr
	      else
	         Diffusion(i,k) = 1./R(i) * ( R(i+1)*Kr*DvDr(i+1,k) - R(i-1)*Kr*DvDr(i-1,k))/(2*Dr)	      
	      endif
	      
	      DvDt(i,k) = - U_se(i,k)*DvDr(i,k) - W_se(i,k)*DvDz(i,k) - U_se(i,k)*V(i,k)/R(i)  &
	                  - fc*U_se(i,k) + Vdot(i,k) + Diffusion(i,k)   
	      
	   enddo
	enddo
	
end subroutine calc_tendency



!--This subroutine calculates the terms of elliptic equations
subroutine calc_terms
	DOUBLE PRECISION :: dKhi_dz,dKhi_dr,  Zeta, R_temp, Min_C
	DOUBLE PRECISION, dimension(Nr,Nz) :: CC, Xi  ! CC=v^2/r + fv, Xi=2*v/r + f 
	DOUBLE PRECISION, dimension(Nr,Nz) :: dCKhi_dZ
	integer :: ir,iz,npoints
	
	DOUBLE PRECISION :: max_delta

	open (3,file='check_ellip.txt')

	A=0.
	B=0.
	C=0.
	npoints=0
    Min_C=9999
	write(*,*) "Calculating the terms of the PDE..."
	do iz=1,Nz
		do ir=2,Nr
		   CC(ir,iz) = V(ir,iz)**2/R(ir) + fc*V(ir,iz)
        enddo
		CC(1,iz)=CC(2,iz)
	enddo

	do iz=1,Nz
		do ir=1,Nr
			R_temp=R(ir)
		    if (iz==1)then
				dKhi_dz = (Khi(ir,iz+1)-Khi(ir,iz))/Dz
				dCKhi_dz(ir,iz) = (Khi(ir,iz+1)*CC(ir,iz+1)-Khi(ir,iz)*CC(ir,iz))/Dz
			else if (iz==Nz) then
				dKhi_dz = (Khi(ir,iz)-Khi(ir,iz-1))/Dz
				dCKhi_dz(ir,iz) = (Khi(ir,iz)*CC(ir,iz)-Khi(ir,iz-1)*CC(ir,iz-1))/Dz
			else
				dKhi_dz = (Khi(ir,iz+1)-Khi(ir,iz-1))/(2.*Dz)
				dCKhi_dz(ir,iz) = (Khi(ir,iz+1)*CC(ir,iz+1)-Khi(ir,iz-1)*CC(ir,iz-1))/(2.*Dz)
			endif

			if (ir==1)then
				Zeta = V(ir+1,iz)/(2*R(ir+1)) + (V(ir+1,iz)-V(ir,iz))/Dr
				DKhi_dr = (Khi(ir+1,iz)-Khi(ir,iz))/Dr
				R_temp=Dr/2
			else if (ir==Nr) then
				Zeta = V(ir,iz)/R(ir) + (V(ir,iz)-V(ir-1,iz))/Dr
				DKhi_dr = (Khi(ir,iz)-Khi(ir-1,iz))/Dr
			else
				Zeta = V(ir,iz)/R(ir) + (V(ir+1,iz)-V(ir-1,iz))/(2.*Dr)
				DKhi_dr = (Khi(ir+1,iz)-Khi(ir-1,iz))/(2.*Dr)
			endif

			Xi(ir,iz) = 2.*V(ir,iz)/R_temp + fc

			A(ir,iz) = - (g*dKhi_dz)/(R_temp*rho(ir,iz))
			B(ir,iz) = -2.*dCKhi_dz(ir,iz)/(R_temp*rho(ir,iz))
			C(ir,iz) = (Xi(ir,iz)*(Zeta+fc)*Khi(ir,iz) + CC(ir,iz)*DKhi_dr)/(R_temp*rho(ir,iz))

			Delta(ir,iz) = B(ir,iz)**2 - 4*A(ir,iz)*C(ir,iz)

			if (Delta(ir,iz)>0) then
				npoints=npoints+1
			endif

            if (Delta(ir,iz)>max_delta) then
			   max_delta=Delta(ir,iz)
			   ire=ir
			   ize=iz
	    endif

            if (C(ir,iz)<Min_C) then
			   Min_C=C(ir,iz)
			endif

			
		enddo
	enddo

	write(*,*) "Number of gridpoints NOT satisfied elliptics condition:",npoints
	write(*,*) "Maximum Delta: ",max_delta
	write(*,*) "...at grid point: ",ire,ize
	write(*,*) "Maximum C: ",Min_C
	!pause

	
	!-Jan 01, 2009 - Adjust the data to overcome ellipticity condition
	if (regularize_field) then

	!-First attempt: increase the I2 GLOBALLY so that it's posivtive everywhere
	!-org: 1.1

	    if (Min_C<0) then	   
	        C(:,:) = C(:,:)	+ abs(1.01*Min_C) 
   	    endif
	
	    Delta(:,:) = B(:,:)**2 - 4*A(:,:)*C(:,:)
	    npoints=0
	    do iz=1,Nz
		    do ir=1,Nr
			    if (Delta(ir,iz)>0) then
				    npoints=npoints+1
			    endif
            enddo 
	    enddo
	    write(*,*) "After 1st regularization:"
	    write(*,*) "Number of gridpoints Still NOT satisfied elliptics condition:",npoints


	!-Second attempt: reduce vertical shear if needed, org: 0.5
	    do iz=1,Nz
		    do ir=1,Nr
			    if (Delta(ir,iz)>=0) then
			        dCKhi_dz(ir,iz) = 0.8*dCKhi_dz(ir,iz)
				    B(ir,iz) = -2.*dCKhi_dz(ir,iz)/(R_temp*rho(ir,iz))
			    endif
            enddo 
	    enddo
	    Delta(:,:) = B(:,:)**2 - 4*A(:,:)*C(:,:)

	    npoints=0
	    do iz=1,Nz
		    do ir=1,Nr
			    if (Delta(ir,iz)>0) then
			 	    npoints=npoints+1
			    endif
            enddo 
	    enddo
	    write(*,*) "After 2nd regularization:"
	    write(*,*) "Number of gridpoints Still NOT satisfied elliptics condition:",npoints
		!pause 'check'
	else
	    write(*,*)"No regularize:!!!",regularize_field
	endif

	do iz=2,Nz-1
		do ir=2,Nr-1
			D(ir,iz) = (A(ir+1,iz)-A(ir-1,iz))/(2.*Dr) + (B(ir,iz+1)-B(ir,iz-1))/(4.*Dz)
			E(ir,iz) = (B(ir+1,iz)-B(ir-1,iz))/(4.*Dr) + (C(ir,iz+1)-C(ir,iz-1))/(2.*Dz)

			F(ir,iz) = g*(Thetadot(ir+1,iz)*Khi(ir+1,iz)**2 - Thetadot(ir-1,iz)*Khi(ir-1,iz)**2)/(2.*Dr)	&
					  +(CC(ir,iz+1)*Thetadot(ir,iz+1)*Khi(ir,iz+1)**2 - CC(ir,iz-1)*Thetadot(ir,iz-1)*Khi(ir,iz-1)**2)/(2.*Dz) &
					  -(Khi(ir,iz+1)*Xi(ir,iz+1)*Vdot(ir,iz+1) - Khi(ir,iz-1)*Xi(ir,iz-1)*Vdot(ir,iz-1))/(2.*Dz)

		enddo
	enddo



end subroutine calc_terms
!--^_^



!--This subroutine use Successive OverRelaxation method the solve an elliptic equation of P
!--  APxx + BPxy + CPyy + DPx + EPy = F
!--
subroutine sor(P,A,B,C,D,E,F,nx,ny,Dx,Dy)
	integer :: Nx, Ny , irec
	DOUBLE PRECISION :: Dx, Dy, P_xx,P_yy,P_xy,P_x,P_y, EIJ
	DOUBLE PRECISION, dimension(Nx,Ny) :: P, A,B,C,D,E,F 

	!--Local variables
	integer ::  max_itt,i,j,itt
	DOUBLE PRECISION :: omega,epsilon,Residual,max_res, Residuale 


	max_itt=6000
	omega=1.8
	epsilon=1.e-20

	write(*,*)"Start SOR subroutine..."

	open(2,file='check.dat',access='direct',form='unformatted',recl=Nr*Nz*4)
	irec=1

	do itt=1,max_itt
	  max_res=-1
 	  do i=2,Nx-1
	    do j=2,Ny-1

		P_xx = ( P(i+1,j) + P(i-1,j) - 2.*P(i,j) ) / Dx**2
		P_yy = ( P(i,j+1) + P(i,j-1) - 2.*P(i,j) ) / Dy**2
		P_xy = ( P(i+1,j+1) - P(i+1,j-1) - P(i-1,j+1) + P(i-1,j-1) ) / (4.*Dx*Dy)
		P_x  = ( P(i+1,j) - P(i-1,j) ) / (2.*Dx)
		P_y  = ( P(i,j+1) - P(i,j-1) ) / (2.*Dy)


		Residual = A(i,j) * P_xx + B(i,j) * P_xy + C(i,j)*P_yy &	
				 + D(i,j)*P_x + E(i,j)*P_y  - F(i,j)

		EIJ = - 2.*A(i,j)/Dx**2 - 2.*C(i,j)/Dy**2



!		!--Update to new value P
		P(i,j)=P(i,j)-omega*Residual/EIJ

		if (max_res < 0) max_res=Residual

		if (Residual>abs(max_res)) max_res=abs(Residual)

		
	    enddo
	  enddo

	  if (max_res<epsilon) then
	    write(*,*) "The itteration converged at itt=",itt 
		return
	  endif
     
   	  if (mod(itt,100)==0) then
	     write(*,*) "Itt, Max_res:",itt,max_res
	  endif


	!-- boundary condition
    !P(Nr,:)=P(Nr-1,:) 	  
  
	enddo

	return
end subroutine sor
!--^_^


!--This subroutine calculate u,w from toroidal function psi
subroutine psi2uw
	integer :: ir,iz
	do ir=2,Nr-1
		do iz=1,Nz-1
			W_se(ir,iz)= (psi(ir+1,iz)-psi(ir-1,iz))/(2.*Dr*R(ir)*Rho(ir,iz))		
			if (iz==1) then
				U_se(ir,iz)=-(psi(ir,iz+1)-psi(ir,iz))/(Dz*R(ir)*Rho(ir,iz))
			else
				U_se(ir,iz)=-(psi(ir,iz+1)-psi(ir,iz-1))/(2.*Dz*R(ir)*Rho(ir,iz))
			endif

		enddo
	enddo
	W_se(1,:)=W_se(2,:)
	U_se(Nr,:)=2*U_se(Nr-1,:) - U_se(Nr-2,:) 

end subroutine psi2uw
!--^_^




END MODULE solve