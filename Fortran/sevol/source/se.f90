!==============================================================================
! SE -	Sawyer-Eliassen Solver - August 19th, 2008 by Bui Hoang Hai, 
!		Hanoi University of Science
!       MM5 dianogstic version
!
! Synopsis:
!		This program solve the Sawyer-Eliassen  
!		in order to get the toroidal stream function in response to the diabatic heating source 
!
! Modules:
!		module_comm:	contains some common parameters, as well as defines the 
!						ideal vortex,
!						NOTE: solving the equation depends critical on the Elliptic condition
!						which depend a lot on the choice of tangential wind profile
!					 
!		module_data:	defines the data for the program
!		module_solve:	solves the Sawyer-Eliassen equation	using Successive Overrelaxation Method (SOR)
!						The SOR subroutine is to solve a general elliptic 
!		module_plot:	plot data/write data file for plotting
!
!==============================================================================


PROGRAM sawyer_eliassen

USE data
USE solve
USE plot


IMPLICIT NONE

Integer :: Tstep,MaxTime

!CALL init_calc_data
!if (balance_field) then
!   CALL init_balance_data
!endif

CALL init_ideal_data

!--Check initial fields (for dianogstic error)
CALL plot_grads

MaxTime=MaxT*60
Call solve_init
Do While (Time<=MaxTime)
  Call new_step
End Do

Call solve_deinit

!CALL calc_terms
!CALL plot_grads
!CALL sor(Psi,A,B,C,D,E,F,nr,nz,Dr,Dz)
!CALL psi2uw
!--Write the results
!CALL plot_grads
!CALL write_ascii


write(*,*) "Finished!"
read(*,*) 

STOP
END