!==============================================================================
!
! This MODULE contains subroutine to output data for grads and other programs
! The data avaible in MODULE DATA and MODULE SOLVE (for debugging purpose)
!
!==============================================================================

MODULE plot

USE data
USE solve

IMPLICIT NONE

CONTAINS

subroutine plot_ncarg
end subroutine plot_ncarg

!--This subroutine outputs to ASCII files some enssential fields
subroutine write_ascii
	!--Calculate some other fields
	DOUBLE PRECISION, dimension(Nr,Nz) :: S,Xi,Eta,M,I2,N2
	DOUBLE PRECISION :: dtheta_dz
	integer :: ir,iz

	write(*,*)"-->Grads, Nr,Nz:",Nr,Nz

	do iz=1,Nz
		do ir=1,Nr
		    if (iz==1) then
			   dtheta_dz = (theta(ir,iz+1)-theta(ir,iz))/Dz
			   S(ir,iz) = (V(ir,iz+1)-V(ir,iz))/Dz
			else if (iz==Nz) then
			   dtheta_dz = (theta(ir,iz)-theta(ir,iz-1))/Dz
			   S(ir,iz) = (V(ir,iz)-V(ir,iz-1))/Dz
			else
			   dtheta_dz = (theta(ir,iz+1)-theta(ir,iz-1))/(2.*Dz)
			   S(ir,iz) = (V(ir,iz+1)-V(ir,iz-1))/(2.*Dz)
			endif

			
		    if (ir==Nr) then
			   Eta(ir,iz) = V(ir,iz)/R(ir) + (V(ir,iz)-V(ir-1,iz))/Dr  + fc
			   Xi(ir,iz) = 2*V(ir,iz)/R(ir) + fc
			   I2(ir,iz) = Xi(ir,iz)*Eta(ir,iz)/R(ir)
			else if (ir>1) then
			   Eta(ir,iz) = V(ir,iz)/R(ir) + (V(ir+1,iz)-V(ir-1,iz))/(2.*Dr)  + fc
			   Xi(ir,iz) = 2*V(ir,iz)/R(ir) + fc
			   I2(ir,iz) = Xi(ir,iz)*Eta(ir,iz)/R(ir)
			endif
			M(ir,iz) = R(ir)*V(ir,iz) + R(ir)**2*fc/2
			N2(ir,iz) = g/theta(ir,iz)*dtheta_dz
		enddo
	enddo

	Eta(1,:) = 2*Eta(2,:) - Eta(3,:) 
	Xi(1,:) = 2*Xi(2,:) - Xi(3,:) 
	I2(1,:) = 2*I2(2,:) - I2(3,:) 



    open(1,file='mm5_u.txt')
	   write(1,*)'mm5_u',Nr,Nz
	   do iz=1,Nz
	      do ir=1,Nr
		      write(1,*) real (U(ir,iz))
		  enddo
	   enddo
	close(1)
    open(1,file='mm5_v.txt')
	   write(1,*)'mm5_v',Nr,Nz
	   do iz=1,Nz
	      do ir=1,Nr
		      write(1,*) real (V(ir,iz))
		  enddo
	   enddo
	close(1)
    open(1,file='mm5_w.txt')
	   write(1,*)'mm5_w',Nr,Nz
	   do iz=1,Nz
	      do ir=1,Nr
		      write(1,*) real (W(ir,iz))
		  enddo
	   enddo
	close(1)

    open(1,file='mm5_vdot.txt')
	   write(1,*)'mm5_vdot',Nr,Nz
	   do iz=1,Nz
	      do ir=1,Nr
		      write(1,*) real (Vdot(ir,iz))
		  enddo
	   enddo
	close(1)

    open(1,file='mm5_thetadot.txt')
	   write(1,*)'mm5_thetadot',Nr,Nz
	   do iz=1,Nz
	      do ir=1,Nr
		      write(1,*) real (Thetadot(ir,iz))
		  enddo
	   enddo
	close(1)

     open(1,file='se_theta.txt')
	   write(1,*)'se_theta',Nr,Nz
	   do iz=1,Nz
	      do ir=1,Nr
		      write(1,*) real (Theta(ir,iz))
		  enddo
	   enddo
	close(1)

    open(1,file='se_psi.txt')
	   write(1,*)'se_psi',Nr,Nz
	   do iz=1,Nz
	      do ir=1,Nr
		      write(1,*) real (Psi(ir,iz))
		  enddo
	   enddo
	close(1)

    open(1,file='se_u.txt')
	   write(1,*)'U_se',Nr,Nz
	   do iz=1,Nz
	      do ir=1,Nr
		      write(1,*) real (U_se(ir,iz))
		  enddo
	   enddo
	close(1)

    open(1,file='se_w.txt')
	   write(1,*)'W_se',Nr,Nz
	   do iz=1,Nz
	      do ir=1,Nr
		      write(1,*) real (W_se(ir,iz))
		  enddo
	   enddo
	close(1)

    open(1,file='se_delta.txt')
	   write(1,*)'Eliptic condition: delta<0',Nr,Nz
	   do iz=1,Nz
	      do ir=1,Nr
		      write(1,*) real (Delta(ir,iz))
		  enddo
	   enddo
	close(1)

    open(1,file='Isquare.txt')
	   write(1,*)'I square',Nr,Nz
	   do iz=1,Nz
	      do ir=1,Nr
		      write(1,*) real (I2(ir,iz))
		  enddo
	   enddo
	close(1)

    open(1,file='Eta.txt')
	   write(1,*)'Eta=V/r + dV/dr + f',Nr,Nz
	   do iz=1,Nz
	      do ir=1,Nr
		      write(1,*) real (Eta(ir,iz))
		  enddo
	   enddo
	close(1)

    open(1,file='Xi.txt')
	   write(1,*)'Xi=2v/r+f',Nr,Nz
	   do iz=1,Nz
	      do ir=1,Nr
		      write(1,*) real (Xi(ir,iz))
		  enddo
	   enddo
	close(1)

    open(1,file='Rho.txt')
	   write(1,*)'Density',Nr,Nz
	   do iz=1,Nz
	      do ir=1,Nr
		      write(1,*) real (Rho(ir,iz))
		  enddo
	   enddo
	close(1)

    open(1,file='P.txt')
	   write(1,*)'Pressure',Nr,Nz
	   do iz=1,Nz
	      do ir=1,Nr
		      write(1,*) real (P(ir,iz))
		  enddo
	   enddo
	close(1)

end subroutine write_ascii


!--This subroutine calculate some fields and writes to grads datafile
subroutine plot_grads
	!--Calculate some other fields
	DOUBLE PRECISION, dimension(Nr,Nz) :: S,Xi,Eta,M,I2,N2
	DOUBLE PRECISION :: dtheta_dz
	integer :: ir,iz

	write(*,*)"-->Grads, Nr,Nz:",Nr,Nz


	do iz=1,Nz
		do ir=1,Nr
		    if (iz==1) then
			   dtheta_dz = (theta(ir,iz+1)-theta(ir,iz))/Dz
			   S(ir,iz) = (V(ir,iz+1)-V(ir,iz))/Dz
			else if (iz==Nz) then
			   dtheta_dz = (theta(ir,iz)-theta(ir,iz-1))/Dz
			   S(ir,iz) = (V(ir,iz)-V(ir,iz-1))/Dz
			else
			   dtheta_dz = (theta(ir,iz+1)-theta(ir,iz-1))/(2.*Dz)
			   S(ir,iz) = (V(ir,iz+1)-V(ir,iz-1))/(2.*Dz)
			endif

			
		    if (ir==Nr) then
			   Eta(ir,iz) = V(ir,iz)/R(ir) + (V(ir,iz)-V(ir-1,iz))/Dr  + fc
			   Xi(ir,iz) = 2*V(ir,iz)/R(ir) + fc
			   I2(ir,iz) = Xi(ir,iz)*Eta(ir,iz)/R(ir)
			else if (ir>1) then
			   Eta(ir,iz) = V(ir,iz)/R(ir) + (V(ir+1,iz)-V(ir-1,iz))/(2.*Dr)  + fc
			   Xi(ir,iz) = 2*V(ir,iz)/R(ir) + fc
			   I2(ir,iz) = Xi(ir,iz)*Eta(ir,iz)/R(ir)
			endif
			M(ir,iz) = R(ir)*V(ir,iz) + R(ir)**2*fc/2
			N2(ir,iz) = g/theta(ir,iz)*dtheta_dz
		enddo
	enddo

	Eta(1,:) = 2*Eta(2,:) - Eta(3,:) 
	Xi(1,:) = 2*Xi(2,:) - Xi(3,:) 
	I2(1,:) = 2*I2(2,:) - I2(3,:) 


    open(1,file='se.dat',access='direct',form='unformatted',recl=Nr*Nz*4)
 	write(1,rec=1) real (U) 
	write(1,rec=2) real (V) 
	write(1,rec=3) real (W) 
	write(1,rec=4) real (T) 
	write(1,rec=5) real (P) 
	write(1,rec=6) real (Rho) 
	write(1,rec=7) real (Theta)
	write(1,rec=8) real (Thetadot)
	write(1,rec=9) real (Psi)
	write(1,rec=10) real (Delta)
	write(1,rec=11) real (Eta)
	write(1,rec=12) real (Xi)
	write(1,rec=13) real (S)
	write(1,rec=14) real (M)
	write(1,rec=15) real (N2)
	write(1,rec=16) real (I2)
	write(1,rec=17) real (A)
	write(1,rec=18) real (B)
	write(1,rec=19) real (C)
	write(1,rec=20) real (D)
	write(1,rec=21) real (E)
	write(1,rec=22) real (F)
	write(1,rec=23) real (Vdot)
	write(1,rec=24) real (U_se) 
	write(1,rec=25) real (W_se) 
	close(1)

end subroutine plot_grads


END MODULE plot