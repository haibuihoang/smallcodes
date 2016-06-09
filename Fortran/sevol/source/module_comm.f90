!
! The module that define some common variable
!


MODULE COMM
IMPLICIT NONE

!--Parameter for the blance vortex
!DOUBLE PRECISION, parameter :: Rm=40e3, Vm=50., Src_r=80e3, Src_z=5e3, Src_sx=10e3,Src_sz=2000,Src_mag=1e-5


! Heating Source:  Src_r, Src_z, Src_sx,Src_sz,Src_mag
! Momentum Source: 	M_Src_r, M_Src_z, M_Src_sx,M_Src_sz,M_Src_mag
! Diffusion Coeficion: Kr,  default=1.e+4   
! Spread speed of heating source: Src_speed, default=2.8m/s

DOUBLE PRECISION:: Rm, Vm, Src_r, Src_z, Src_sx,Src_sz,Src_mag, Src_T, &
		   M_Src_r, M_Src_z, M_Src_sx,M_Src_sz,M_Src_mag,Kr,Src_speed

logical :: barotropic,friction

INTEGER :: DT, MaxT    ! DT: integration timestep (seconds), MaxT: Max integration Time (minutes)

!--	simplified=0: Full SE equation, =1: Density is a function of z, =2: Density is a constant
integer :: simplified	  


! fc:Coriolis parameter, !hpi:one half pi, g: gravity acc.., pR:1015 mb in Pa, 	  fc=4.97e-5
! Rd:Specific gas constant
DOUBLE PRECISION, parameter :: hpi=1.570796, fc=0.97e-5, g=9.8,pR=101500.0, Rd=287. 

!--For the friction force in the boundary: Cd is the Drag coefiction, z0: vertical scale, H: boundary height
 
DOUBLE PRECISION, parameter :: Cd=2.0e-3, z0=400., H=400.


DOUBLE PRECISION :: alp1,alp2,v1,v2,pi





CONTAINS

subroutine init_comm
	DOUBLE PRECISION :: ch
	namelist /options/ barotropic,friction, simplified, Rm, Vm, Src_r, Src_z, Src_sx,Src_sz,Src_mag, Src_T, &
	                   M_Src_r, M_Src_z, M_Src_sx,M_Src_sz,M_Src_mag, Kr, Src_speed, DT, MaxT
	
	open(1,file='options.nls')
	read(1,options)
	close(1) 
	
	write(*,options)

	ch  = 0.3
	alp2 = 0.15
	alp1 = (1 - alp2*ch*exp(-alp2))/(1 - ch*exp(-alp2))
	v2 = ch*vm
	v1 = vm*exp(alp1)*(1 - ch*exp(-alp2))
    pi=4*atan(1.)

end subroutine init_comm

END MODULE COMM


