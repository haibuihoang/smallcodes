!==============================================================================
!
! The module that define and initialize data for the Sawyer-Eliassen Solver
!
!==============================================================================

MODULE data

USE COMM

IMPLICIT NONE

!--Options
!INTEGER, parameter :: hours=48
INTEGER, parameter :: s_time=48	, e_time=48
LOGICAL, parameter :: vdot_only = .true.
LOGICAL, parameter :: thetadot_only = .false.

LOGICAL, parameter :: balance_field = .true.
LOGICAL, parameter :: regularize_field = .true.
LOGICAL, parameter :: add_eddy = .true.


!--Define the dimensions of domain
INTEGER, parameter :: Nr=101, Nz=31	              ! Nz = 15, 23, 31, 76 
INTEGER, parameter :: Np=23
DOUBLE PRECISION, parameter ::  Dr=10e3, Dz=500.	  ! Dz = 1000, 700, 500, 200  

INTEGER :: time_step


!--Define state variables 
!--we will solve the Sawyer-Eliassen for U and W  (radial and vertical wind component)
!--Some other variables: V (tangential wind component), P: pressure, Rho: air density
!--			 T: Temperature, Theta: Potential T, Psi: Toroidal streamfunction
DOUBLE PRECISION, dimension(Nr,Nz) :: U, V, W, P, Rho, T, Theta, & 
				      Psi, Delta, Khi, Thetadot, Vdot, &
				      Fdif,Fpbl, &
				      U_se,W_se, dvdz_rz 
DOUBLE PRECISION, dimension(Nz) :: Z	! Define vertical coordinates
DOUBLE PRECISION, dimension(Nr) :: R	! Define radial coordinates

DOUBLE PRECISION :: Rdom,Zdom


!-- Data on pressure levels
DOUBLE PRECISION, dimension(Np) :: Plevs

REAL, dimension(Nr,Np) :: Up, Vp, Wp, Thetap, Rhop, GHeight, Heatp, Fdif_p, Fpbl_p, &
                          T2p, T4p, T5p, V2p, V4p, Vdotp        ! Additional eddy term (09 May 30)   

Real, parameter :: PP=0 ! Options for birational interpolation, =0 means bicubic

DATA Plevs /1000, 975, 950, 925, 900, 875, 850, 825, 800, 750, 700,   &
                 650,  600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100/



CONTAINS


!--This subroutine will initialize the domain with an balanced vortex
!--follows Smith 2006 (Tellus).
subroutine init_ideal_data()
    integer, parameter :: m=2
	integer :: i,k, ir, iz

	DOUBLE PRECISION :: x,  Rtempx, Rtempz
	DOUBLE PRECISION :: alogrhoR, alogrho, dvdz,pRz,rhoRz,TRz,xi
	DOUBLE PRECISION, dimension(Nr, Nz) :: dCdz  
	DOUBLE PRECISION, dimension(m) :: y

	write(*,*) "Intializing data..."
	write(*,*) "Rm,Vm",Rm,Vm

	


	!--Init some common variables from COMM module	
	call init_comm

	!--Define some terms in the thermal wind equation
	do k=1,Nz
		Z(k)=(k-1)*Dz
	enddo
	do i=1,Nr
		R(i)=(i-1)*Dr
	enddo
	Rdom=R(Nr)
	Zdom=Z(Nz)

!--Step1: Calc V and Dv/Dz on (r,z) gridpoints
	do k=1,Nz
	   do i=1,Nr
		  call vortex1 (r(i),Z(k),V(i,k),dvdz_rz(i,k))
	   enddo
	enddo
	!V(1,:) = 0

!--Step2:Defining environment sounding at R=r(Nr)   

	write(*,*) "Defining environment sounding..."
	open(1,file='envini.txt')
	write(1,'(5A15)') "Z","P","Rho","T","Theta"

	!--Define the boundary condition at r=0 and r=R(nr)
	do k = 1,Nz
		call pres1 (Z(k),P(Nr,k),Rho(Nr,k),T(Nr,k))
		Theta(Nr,k)  = T(Nr,k)*(1.0e3/P(Nr,k))**0.2865
		Khi(Nr,k)  =  1./Theta(Nr,k)
		write(1,'(5f15.3)') Z(k),P(Nr,k),Rho(Nr,k),T(Nr,k),Theta(Nr,k)
	enddo

	close(1)



!--Step 3: calculate pressure at each gridpoint by integrating outward and density by integrating inward
	write(*,*) "Calculating state variables..."

! Loop through eacg grid points
	do iz = 1,Nz
		write (*,*) 'Z = ',int(Z(iz)),' m'
		do ir = 1,Nr
		
			if (ir > 1) then
				xi = 2.*V(ir,iz)/R(ir) + fc
				dCdz(ir,iz) = xi*dvdz_rz(ir,iz)
			end if      
   
			! START INTEGRATION FROM (r(ir), z(iz)) TO DETERMINE THE HEIGHT y(1) OF THE ISOBARIC 
			! SURFACE THROUGH THIS POINT AT RADIUS r(n) WHERE THE PRESSURE, pRz, TEMPERATURE, 
			! TRz, AND DENSITY, rhoRz, ARE KNOWN. ALSO CALCULATED IS THE LOGARITHM OF DENSITY
			! DIFFERENCE BETWEEN  (r(ir), z) AND (r(n), y(1)

			y(1) = Z(iz)          ! HEIGHT OF THE ISOBARIC SURFACE
			y(2) = 0              ! FIRST GUESS FOR THE LOGARITHM OF DENSITY DIFFERENCE BETWEEN
					      ! r(ir) AND r(Nr)

			do i = ir,Nr 
			    x = r(i) 
			    call intg (x,y,m,Dr)     
			end do ! OVER i

			
			call pres1 (y(1),pRz,rhoRz,TRz)              ! CALCULATES pRz, TRz AND rhoRz
			p(ir,iz)      = pRz                          ! PRESSURE AT (r(ir),z)   
			

			alogrhoR      = log(rhoRz)                   ! CALCULATES ln(rho) AT r = R (= r(n))
			alogrho       = alogrhoR - y(2)              ! CALCULATES ln(rho) AT (r(ir),z)
			Rho(ir,iz)    = exp(alogrho)                 ! CALCULATES rho AT (r(ir),z)
			T(ir,iz)      = pRz/(Rd*rho(ir,iz))          ! TEMPERATURE AT (r(ir),z)
			Theta(ir,iz)  = T(ir,iz)*(1000./pRz)**0.2865 ! POTENTIAL TEMPERATURE AT (r(ir),z)
			Khi(ir,iz)  =  1./Theta(ir,iz)
		
			! INIT A MOMENTUM SOURCE AT (M_Src_r,M_Rrc_z)
			Rtempx = ABS(R(ir)-M_Src_r)
			Rtempz = ABS(Z(iz)-M_Src_z)
			if ((Rtempx<=M_Src_sx) .and. (Rtempz<=M_Src_sz)) then
			    Vdot(ir,iz)=M_Src_mag*cos(pi*(Rtempx/M_Src_sx)/2)*cos(pi*(Rtempz/M_Src_sz)/2)
			endif
			

		end do ! OVER ir
 
	end do   ! OVER iz


!	pause					 


	if (friction) then
   	   write(*,*)"Initilizing the friction..."
 	   do iz=1,Nz
	      do ir=1,Nr
		 Vdot(ir,iz) = Vdot(ir,iz) - Cd*(.9*V(ir,1))**2.*exp( -1.*(Z(iz)/z0)**2 ) / H
	      enddo
	   enddo
	endif


	write(*,*) 'Simply condiftion:',simplified
   	do ir=1,Nr
	    do iz=1,Nz
		!--Now, simplify the air density if neccessary
		if (simplified == 1) then  !--density is a function of z only
			Rho(ir,iz)=Rho(Nr,iz)	
		elseif (simplified==2) then !--density is a constant
			Rho(ir,iz)=Rho(Nr,Nz)	
	 	endif
	    enddo
	enddo


end subroutine init_ideal_data


!==============================================================================
!--This subroutine will initialize calculated data from MM5's simulation
subroutine init_calc_data()
    DOUBLE PRECISION, dimension(Nr) :: Pslv
	real :: temp, temp_srf
    integer :: i,k,ir,iz, dt, itime, startrec1,startrec2

    dt = 15
	itime = (s_time-6)*60/dt    
	startrec1 = itime*(Np+1) + 2  !startrec1 for Np+1 levels
	startrec2 = itime*Np + 1 	  !startrec1 for Np levels
	time_step = 1 + (e_time-s_time)*60/dt
	write(*,*)'itime',itime,startrec1,startrec2,time_step
	!pause

	!--Init some common variables from COMM module	
	call init_comm

	!--Define some terms in the thermal wind equation
	do k=1,Nz
		Z(k)=(k-1)*Dz
	enddo
	do i=1,Nr
		R(i)=(i-1)*Dr
	enddo

!-- First, read data on pressure surface
    write(*,*) "Reading data..."
    
	!--Read Geopotential Height
    CALL read_avaraged_data(GHeight,'Input/gheight.grd   ',startrec1,time_step,1,0)
    !--Read Radial wind
    CALL read_avaraged_data(Up     ,'Input/syuu.grd      ',startrec1,time_step,1,0)
    Up(1,:)=0.0

    !--Read tangential wind
    CALL read_avaraged_data(Vp     ,'Input/syvv.grd      ',startrec1,time_step,1,0)
    !--Read Radial wind
    CALL read_avaraged_data(Wp     ,'Input/syww.grd      ',startrec1,time_step,1,0)


    !--Read Density - Note: zrev
    CALL read_avaraged_data(Rhop   ,'Input/rho.grd       ',startrec2,time_step,0,1)
    !--Read Potential Temperature - Note: zrev
    CALL read_avaraged_data(Thetap ,'Input/theta.grd     ',startrec2,time_step,0,1)
	!--Read  condensational heating (K/hr) - NOTE: zrev
    CALL read_avaraged_data(Heatp  ,'Input/heat.grd      ',startrec2,time_step,0,1)


    !--Read Horizontal Diffusion
    CALL read_avaraged_data(Fdif_p ,'Input/dif.grd       ',startrec2,time_step,0,0)
    !--Read boundary friction
    CALL read_avaraged_data(Fpbl_p ,'Input/pbl.grd       ',startrec2,time_step,0,0)


    !--Read Horizontal Diffusion
    CALL read_avaraged_data(Fdif_p ,'Input/dif.grd       ',startrec2,time_step,0,0)
    !--Read boundary friction
    CALL read_avaraged_data(Fpbl_p ,'Input/pbl.grd       ',startrec2,time_step,0,0)

    !--Read Eddy terms
    CALL read_avaraged_data(T2p ,'Input/secondT2.grd     ',startrec1,time_step,1,0)
    CALL read_avaraged_data(T4p ,'Input/fourthT4.grd     ',startrec1,time_step,1,0)
    CALL read_avaraged_data(T5p ,'Input/fifth.grd     ',startrec1,time_step,1,0)
    CALL read_avaraged_data(V2p ,'Input/secondV2.grd     ',startrec1,time_step,1,0)
    CALL read_avaraged_data(V4p ,'Input/fourthV4.grd     ',startrec1,time_step,1,0)


	open(1,file='checkinput.dat',access='direct',form='unformatted',recl=Nr*Np*4)
	write(1,rec=1) real (Up) 
	write(1,rec=2) real (Vp) 
	write(1,rec=3) real (Wp) 
	write(1,rec=4) real (Rhop) 
	write(1,rec=5) real (Thetap)
	write(1,rec=6) real (Heatp)
	write(1,rec=7) real (Fdif_p)
	write(1,rec=8) real (Fpbl_p)
	write(1,rec=9) real (GHeight)
	write(1,rec=10) real (T2p)
	write(1,rec=11) real (T4p)
	write(1,rec=12) real (V2p)
	write(1,rec=13) real (V4p)
	write(1,rec=14) real (Vdotp)
	write(1,rec=15) real (T5p)
	close(1)



!-- Change unit of fields to SI
!	Heatp(:,:) = Heatp(:,:)/3600.
    Heatp(:,:)=Heatp(:,:)/3600. 

!-- Now, do the interpolation
	Vdotp(:,:)= - ( Fpbl_p(:,:) + Fdif_p(:,:) ) / (96.5*900)   
       	        



! add the eddy terms (09May30)
    if (add_eddy) then
	   Vdotp(:,:) = Vdotp(:,:) - V2p(:,:)/100000. - V4p(:,:)/1000.  
	   Heatp(:,:) = Heatp(:,:) - T2p(:,:) - T4p(:,:)/1000.-T5p(:,:)
	endif




!	Vdot(:,:) =  -1. * (  Fpbl(:,:) + Fdif(:,:)	)/(96.5*900.)
!	Vdot(:,:) =  -1. * Fpbl(:,:)/(96.5*900.)
!	Vdot(:,:) =  -1. * Fdif(:,:)	/(96.5*900.)


    !--First, do a simple reduction for mean sea level pressure
	do i=1,Nr
	   temp  = Thetap(i,1)*(Plevs(1)/1.0e3)**0.2865
	   temp_srf = temp + 0.0065*Gheight(i,1)
	   temp = 0.5*(temp+temp_srf)
	   pslv(i) = Plevs(1)*exp(g*Gheight(i,1)/(Rd*temp))
	enddo


	!--Do some linear interpolation according to Z
	write(*,*)"Field int Vp->V"
	call field_linear_int(DBLE(Vp),V)
	call field_linear_int(DBLE(Up),U)	
	call field_linear_int(DBLE(Wp),W)
	call field_linear_int(DBLE(Rhop),Rho)
	call field_linear_int(DBLE(Thetap),Theta)
	call field_linear_int(DBLE(Heatp),Thetadot)
	call field_linear_int(DBLE(Vdotp),Vdot)

	open(1,file='checkint.dat',access='direct',form='unformatted',recl=Nr*Nz*4)
	write(1,rec=1) real (U) 
	write(1,rec=2) real (V) 
	write(1,rec=3) real (W) 
	write(1,rec=4) real (Rho) 
	write(1,rec=5) real (Theta)
	write(1,rec=6) real (Thetadot)
	write(1,rec=7) real (Fdif)
	write(1,rec=8) real (Fpbl)
	close(1)

    open(1,file='mm5_theta.txt')
	   write(1,*)'mm5_theta',Nr,Nz
	   do iz=1,Nz
	      do ir=1,Nr
		      write(1,*) real (Theta(ir,iz))
		  enddo
	   enddo
	close(1)


   do iz=1,Nz
      do ir=1,Nr
	    call linear_int(DBLE(Gheight(ir,:)),Plevs,Np,Z(iz),P(ir,iz))
        if (iz==1)then
           dvdz_rz(ir,iz)=(v(ir,iz+1)-v(ir,iz))/dz
        else if (iz==Nz) then
	       dvdz_rz(ir,iz)=(v(ir,iz)-v(ir,iz-1))/dz 
        else
	       dvdz_rz(ir,iz)=(v(ir,iz+1)-v(ir,iz-1))/(2.*dz)
        endif
      enddo
    enddo	   



    if (thetadot_only) then
!      Suppress friction    
  	   Vdot = 0.
    endif

    if (vdot_only) then
!      Suppress Heating    
   	   Thetadot = 0.
    endif 

    Khi(:,:) =  1./Theta(:,:)

!Idealized
	if (friction) then
		write(*,*)"Initilizing idealized friction..."
		do iz=1,Nz
			do ir=1,Nr
				Vdot(ir,iz) = Vdot(ir,iz) - Cd*(.9*V(ir,1))**2.*exp( -1.*(Z(iz)/z0)**2 ) / H
			enddo
		enddo
	endif

end  subroutine init_calc_data



!==============================================================================
!--This subroutine calculates balanced vortex following Smith 2006.
subroutine init_balance_data()
	INTEGER, parameter :: m=2
	integer :: i,k, ir, iz
	DOUBLE PRECISION :: x,  Rtempx, Rtempz
	DOUBLE PRECISION :: alogrhoR, alogrho, dvdz,pRz,rhoRz,TRz,xi,pi
	DOUBLE PRECISION, dimension(Nr, Nz) :: dCdz  
	DOUBLE PRECISION, dimension(m) :: y

	pi=4*atan(1.)

	open(1,file='env.txt')
	write(1,'(5A15)') "Z","P","Rho","T","Theta"


    !--Recalculate dv/dz on height gridpoint (dvdz_rz)
	do iz=1,Nz
      do ir=1,Nr
        if (iz==1)then
           dvdz_rz(ir,iz)=(v(ir,iz+1)-v(ir,iz))/dz
        else if (iz==Nz) then
           dvdz_rz(ir,iz)=(v(ir,iz)-v(ir,iz-1))/dz 
        else
           dvdz_rz(ir,iz)=(v(ir,iz+1)-v(ir,iz-1))/(2.*dz)
        endif
      enddo
    enddo	   


    !--No need to define the boundary condition! 


    ! 	Recalculating state variables..."

    ! Loop through each grid points
    do iz = 1,Nz
	!write (*,*) 'Z = ',int(Z(iz)),' m'
        !check!! Rho over time has problems	
	do ir = 1,Nr-1
	    !write (*,*) 'R = ',int(R(ir)/1000.),' km'
	    !we dont need to call this again
	    !call vortex (R(ir),Z(iz),V(ir,iz),dvdz)
   	    dvdz=dVdZ_rz(ir,iz)
	    if (ir > 1) then
		xi = 2.*V(ir,iz)/R(ir) + fc
		dCdz(ir,iz) = xi*dvdz
	    end if      
   
	    ! START INTEGRATION FROM (r(ir), z(iz)) TO DETERMINE THE HEIGHT y(1) OF THE ISOBARIC 
 	    ! SURFACE THROUGH THIS POINT AT RADIUS r(n) WHERE THE PRESSURE, pRz, TEMPERATURE, 
 	    ! TRz, AND DENSITY, rhoRz, ARE KNOWN. ALSO CALCULATED IS THE LOGARITHM OF DENSITY
	    ! DIFFERENCE BETWEEN  (r(ir), z) AND (r(n), y(1)

	    y(1) = Z(iz)          ! HEIGHT OF THE ISOBARIC SURFACE
	    y(2) = 0              ! FIRST GUESS FOR THE LOGARITHM OF DENSITY DIFFERENCE BETWEEN
				  ! r(ir) AND r(Nr)

	    do i = ir,Nr-1 
	       x = r(i) 
	       call intg (x,y,m,Dr)     
	    end do ! OVER i

	
	    call pres1 (y(1),pRz,rhoRz,TRz)               ! CALCULATES pRz, TRz AND rhoRz
	    p(ir,iz)      = pRz                          ! PRESSURE AT (r(ir),z)   
	    alogrhoR      = log(rhoRz)                   ! CALCULATES ln(rho) AT r = R (= r(n))
	    alogrho       = alogrhoR - y(2)              ! CALCULATES ln(rho) AT (r(ir),z)
	    Rho(ir,iz)    = exp(alogrho)                 ! CALCULATES rho AT (r(ir),z)
	    T(ir,iz)      = pRz/(Rd*rho(ir,iz))          ! TEMPERATURE AT (r(ir),z)
	    Theta(ir,iz)  = T(ir,iz)*(1000./pRz)**0.2865 ! POTENTIAL TEMPERATURE AT (r(ir),z)
	    Khi(ir,iz)  =  1./Theta(ir,iz)
		
	    ! INIT A MOMENTUM SOURCE AT (M_Src_r,M_Rrc_z)
	    Rtempx = ABS(R(ir)-M_Src_r)
	    Rtempz = ABS(Z(iz)-M_Src_z)
	    if ((Rtempx<=M_Src_sx) .and. (Rtempz<=M_Src_sz)) then
	        Vdot(ir,iz)=M_Src_mag*cos(pi*(Rtempx/M_Src_sx)/2)*cos(pi*(Rtempz/M_Src_sz)/2)
 	    endif
			
	end do ! OVER ir
     end do   ! OVER iz


     write(*,*) 'Simply condiftion:',simplified
     do ir=1,Nr
        do iz=1,Nz
	!--Now, simplify the air density if neccessary
	    if (simplified == 1) then  !--density is a function of z only
	    Rho(ir,iz)=Rho(Nr,iz)	
	    elseif (simplified==2) then !--density is a constant
	    Rho(ir,iz)=Rho(Nr,Nz)	
	    endif
	enddo
    enddo


end subroutine init_balance_data



!... THIS SUBROUTINE ACCEPTS THE HEIGHT IN METRES AND RETURNS THE
!    PRESSURE p, DENSITY rho, AND VIRTUAL TEMPERATURE T, FOR AN
!    APPROXIMATION TO THE JORDAN SOUNDING
subroutine pres(zz,pp,rho1,T1)
   DOUBLE PRECISION :: zz, pp, rho1, T1, Rd,Kappa,Theta1
   data Rd,Kappa / 287.,.286 /
   call linear_int(Z,P(Nr,:),Nz,zz,pp)
   call linear_int(Z,Theta(Nr,:),Nz,zz,Theta1)
   T1 =  Theta1*(pp/1000.)**Kappa
   rho1  = pp/(Rd*T1)
   return
end subroutine pres




!... This subroutine calculates, v and dvdz in nondimensional units
!    as functions of the nondimensional radius s
!... CALCULATE THE TANGENTIAL VELOCITY AND DERIVATIVES
subroutine vortex (rr,zz,vv,dvdz)
   DOUBLE PRECISION :: rr, zz, vv, dvdz   
   CALL bilinear_int(R,Z,Nr,Nz,V,rr,zz,vv)
   CALL bilinear_int(R,Z,Nr,Nz,dvdz_rz,rr,zz,dvdz)
   return
end subroutine vortex



!... THIS SUBROUTINE ACCEPTS THE HEIGHT IN METRES AND RETURNS THE
!    PRESSURE p, DENSITY rho, AND VIRTUAL TEMPERATURE T, FOR AN
!    APPROXIMATION TO THE JORDAN SOUNDING
subroutine pres1(z,p,rho,T)
   DOUBLE PRECISION :: z, p, rho, T, Rd, Ts, ps, gamma1, c1
   data Rd,Ts,ps,gamma1,c1 / 287., 303., 1015., 0.212151E-04, 5.31593 /
   T    = Ts*(1 - gamma1*z)
   p    = ps*(1 - gamma1*z)**c1
   rho  = p/(Rd*T)
   return
end subroutine pres1


!... This subroutine calculates, v and dvdz in nondimensional units
!    as functions of the nondimensional radius s
!... CALCULATE THE TANGENTIAL VELOCITY AND DERIVATIVES
subroutine vortex1 (rr,zz,vv,dvdz)
DOUBLE PRECISION :: rr, zz, vv, dvdz
DOUBLE PRECISION :: v0, bb,  ss, zdH
   real,parameter :: b	= 0.53
   ss	 = rr/Rm
   !--Original in baroheight.for
   !v0   = v1*ss*exp(-alp1*ss) + v2*ss*exp(-alp2*ss)
	
   !--new one as DeMaria 1987
   v0	 = Vm*ss*exp((1-ss**b)/b)

   bb   = hpi/Zdom
   zdH  = bb*zz
   if (.not. barotropic ) then 
    	vv    = v0*cos(zdH)
	dvdz = -v0*bb*sin(zdH)
   else
        vv = v0
	dvdz = 0.
   endif
return
end subroutine vortex1


! THIS SUBROUTINRE INTEGRATES n FIRST-ORDER DIFFERENTIAL EQUATIONS
! USING KUTTA'S THREE EIGHTH'S RULE.
! *** NOTE THAT x IS UPDATED BY THE ROUTINE ***
subroutine intg (x,y,n,h)
  integer :: i,n
  DOUBLE PRECISION :: x, h, dx
  DOUBLE PRECISION, dimension(n) :: y,yt,y1,y2,y3,y4

     dx = 0.33333333*h
     call deqn (x,y,y1,n)
     do 1 i = 1,n
1    yt(i) = y(i) + dx*y1(i)
     x = x + dx
     call deqn (x,yt,y2,n)
     do 2 i = 1,n
2    yt(i) = y(i) - dx*y1(i) + h*y2(i)
     x = x + dx
     call deqn (x,yt,y3,n)
     do 3 i = 1,n
3    yt(i) = y(i) + h*(y1(i) - y2(i) + y3(i))
     x = x + dx
     call deqn (x,yt,y4,n)
     do 4 i = 1,n
4    y(i) = y(i) + 0.125*h*(y1(i) + 3.0*( y2(i) + y3(i) ) + y4(i))
     x = x - h
  return
end	subroutine intg


subroutine deqn (r,y,f,n)
 integer :: n
 DOUBLE PRECISION :: vv, r, dvdz
 DOUBLE PRECISION,dimension(n) :: y,f
!...  Calculate v at radius r and height z
 call vortex (r,y(1),vv,dvdz)
 if (r > 0) then
  f(1) = vv*(vv/r + fc)/g
  f(2) = -(2*vv/r + fc)*dvdz/g
 else
  f(1) = 0
  f(2) = -fc*dvdz/g 
 end if
 return
end	subroutine deqn


!==============================================================================
!-- This function do a simple linear interpolation from pressure levels to height levels
subroutine field_linear_int(Fp,F)
    DOUBLE PRECISION, dimension(Nr,Nz) :: F
    DOUBLE PRECISION, dimension(Nr,Np) :: Fp
	integer :: i,k
    do i=1,Nr
	   do k=1,Nz
!	      CALL linear_int(DBLE(Gheight(i,:)),Fp(i,:),Np,Z(k),F(i,k))  
	      CALL bilinear_int(R,DBLE(Gheight(i,:)),Nr,Np,Fp,R(i),Z(k),F(i,k))		  
	   enddo
	enddo
end subroutine field_linear_int

!==============================================================================
!-- This subroutine perform a linear interpolation
subroutine linear_int(X,F,Nx,xx,ff)
integer :: Nx,i
DOUBLE PRECISION,dimension(Nx) :: X,F
DOUBLE PRECISION :: xx,ff

if (xx .le. X(1)) then
   ff=F(1)
else if(xx .ge. X(Nx)) then
   ff=F(Nx)
else
   do i=1,Nx
      if ((xx .ge. X(i)) .and. (xx .le. X(i+1))) then
	     ff = F(i) + (F(i+1) - F(i))*(xx - X(i)) / (X(i+1)-X(i))
	     exit
      endif
   enddo
endif
end subroutine linear_int
!==============================================================================

!==============================================================================
!-- This subroutine perform a bilinear interpolation by combined 2 linear one.
!-- This seems to be simple but ... S.L.O.W
subroutine bilinear_int1(X,Y,Nx,Ny,F,xx,yy,ff)
integer :: Nx,Ny,i,j
DOUBLE PRECISION,dimension(Nx) :: X
DOUBLE PRECISION,dimension(Ny) :: Y
DOUBLE PRECISION,dimension(Nx,Ny) :: F 
DOUBLE PRECISION,dimension(Ny) :: Fytemp 
DOUBLE PRECISION :: xx,yy,ff
  
do i=1,Nx
   do j=1,Ny
      call linear_int(X,F(:,j),Nx,xx,Fytemp(j))
   enddo
   do j=1,Ny
      call linear_int(Y,Fytemp,Ny,yy,ff)
   enddo
enddo

end subroutine bilinear_int1
!==============================================================================


!==============================================================================
!-- This subroutine perform a better linear interpolation
!-- Haibuihoang 2008
subroutine bilinear_int(X,Y,Nx,Ny,F,xin,yin,ff)

integer,INTENT(IN) :: Nx,Ny 
DOUBLE PRECISION, INTENT(IN) ::	 xin,yin
DOUBLE PRECISION,dimension(Nx),INTENT(IN) :: X
DOUBLE PRECISION,dimension(Ny),INTENT(IN) :: Y
DOUBLE PRECISION,dimension(Nx,Ny),INTENT(IN) :: F 
DOUBLE PRECISION,INTENT(OUT)::ff

integer :: i1,j1,i2,j2
integer :: L,R,m,RmL
DOUBLE PRECISION :: xx,yy,temp
xx=xin
yy=yin
!Check bound
if (xx < X(1))  xx=X(1)
if (xx > X(Nx)) xx=X(Nx)
if (yy < Y(1))  yy=Y(1)
if (yy > Y(Ny))  yy=Y(Ny)

!Use binary searching to search for the index i
L=1
R=Nx
m=L+(R-L)/2
RmL=R-L
do while (RmL > 1)
   if (xx > X(m)) then
      L = m
   else
      R = m
   endif
   m=L+(R-L)/2
   RmL=R-L
enddo
i1=L
i2=i1+1

!Use binary searching to search for the index j
L=1
R=Ny
m=L+(R-L)/2
RmL=R-L
do while (RmL > 1)
   if (yy > Y(m)) then
      L = m
   else
      R = m
   endif
   m=L+(R-L)/2
   RmL=R-L
enddo
j1=L
j2=j1+1

!http://en.wikipedia.org/wiki/Bilinear_interpolation

temp = 1./( (x(i2)-x(i1))*(y(j2)-y(j1)) )

ff =   F(i1,j1)*temp*( x(i2) - xx )*( y(j2) - yy ) & 
     + F(i2,j1)*temp*( xx - x(i1) )*( y(j2) - yy ) &
     + F(i1,j2)*temp*( x(i2) - xx )*( yy - y(j1) ) &
     + F(i2,j2)*temp*( xx - x(i1) )*( yy - y(j1) ) 

end subroutine bilinear_int
!==============================================================================

subroutine read_avaraged_data(outdata,inputFile,startrec,nstep,hasSurface,zrev)
character(len=20) :: inputFile
integer,intent(in) :: hasSurface,zrev,startrec,nstep
real,dimension(Nr,Np), intent(OUT) :: outdata
real,dimension(Nr,Np) :: tempdata
integer :: irec,istep,k

write(*,*) 'Opening file ',inputFile

!    dt = 15
!	itime = (hours-6)*60/dt    
!	startrec1 = itime*(Np+1) + 2 
!	startrec2 = itime*Np + 1 

open(1,file=inputFile,access='direct',form='unformatted',status='old',recl=Nr*4)
irec = startrec
outdata = 0.
!write(*,*) 'irec=',irec
!pause
do istep=1,nstep
   if (zrev .eq. 1) then
       do k=Np,1,-1
          read(1,rec=irec) tempdata(:,k)
	  irec=irec+1
       enddo
    else
	  do k=1,Np
	     read(1,rec=irec) tempdata(:,k)
	     irec=irec+1
	  enddo
    endif
    if (hasSurface .eq. 1) irec = irec+1  ! Skip the surface level
    outdata(:,:) = outdata(:,:) + tempdata(:,:)/nstep
enddo
close(1)

end subroutine read_avaraged_data


END MODULE data


