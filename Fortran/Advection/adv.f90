
!
! Test the pollution propagation using Langarange-Euler hybrid
! Pullutants are single points that is propagate independantly
! Pulltion fields are calculated by counting pollutants in each cells
! Diffusion maybe added simply by gridscale fields
! First trials, simulate a single particle in a circular statis wind field
!
implicit none
integer,parameter :: Nx=100,Ny=100,dt=3600,Np=1  ! ONE hour time step!!
real, parameter :: Dx=5000.         ! 5km grid size,-->500km domain
real, parameter :: Rm=20000.,Vm=10. ! Rankine Vortex
real :: Px, Py    ! Particle positions and velocity
real, dimension(Nx,Ny) :: u,v
real, dimension(Nx)::Rx
real, dimension(Ny)::Ry


integer :: i,j,t,nrec,tmax=3600*24*3
real :: R,Cx,Cy,Drx,Dry,Vt

nrec=tmax/dt+1
!--init a simple rankin vortext
do i=1,Nx
  Rx(i)=(i-1)*Dx
enddo
do i=1,Ny
  Ry(i)=(i-1)*Dx
enddo

Cx=Rx(Nx)/2
Cy=Ry(Ny)/2

u=0
v=0
do i=1,Nx
  Drx=Rx(i)-Cx
  do j=1,Ny
    Dry=Ry(j)-Cy
    R=sqrt( Drx**2+Dry**2 )
    if (R<Rm)then
      Vt=R*Vm/Rm
    else
      Vt=Vm*Rm/R
    endif
    if (R>0)then
      u(i,j)=-Vt*Dry/R
      v(i,j)=Vt*Drx/R
    endif
  enddo
enddo

!write test grads file
open(1,file='uv.dat',access='direct',form='unformatted',recl=Nx*Ny*4)
write(1,rec=1)u
write(1,rec=2)v
close(1)


!Calculate advection
!initial position
Px=Cx+44000
Py=Cy+20000

open(1,file="adv.txt")
write(1,*)nrec
t=0
do while(t<=tmax)
  call new_position()
  write(*,*)t,Px-Cx,Py-Cy
  write(1,*)Px,Py
  t=t+dt
enddo
close(1)


stop

CONTAINS

subroutine new_position()
  real :: u1,v1,u2,v2,um,vm
  real :: x2,y2  
  integer :: itt
  call bilint(Px,Py,u1,v1)
  u2=u1
  v2=v1
  do itt=1,20
    um=(u2+u1)/2
    vm=(v2+v1)/2
    x2 = Px + um*Dt/2
    y2 = Py + vm*Dt/2
    call bilint(x2,y2,u2,v2)
    !write(*,*)"Um,Vm",um,vm
  enddo
  Px=x2
  Py=y2


end subroutine new_position


!BiLinear Interpolate for u and v at position x, y 
subroutine bilint(x,y,uout,vout)
  real:: x,y,uout,vout,X1,X2,Y1,Y2
  integer :: ix,iy
  ix = int(x/Dx)+1
  iy = int(y/Dx)+1
  if (ix<1 .or. ix>Nx-1 .or. iy<1 .or. iy>Ny-1)then
     uout=0
     vout=0
     return
  endif
  X1=Rx(ix)
  X2=Rx(ix+1)
  Y1=Ry(iy)
  Y2=Ry(iy+1)

  uout= ( u(ix,iy)*(X2-x)*(Y2-y) + &
          u(ix+1,iy)*(x-X1)*(Y2-y)+&
          u(ix,iy+1)*(X2-x)*(y-Y1)+&
          u(ix+1,iy+1)*(x-X1)*(y-Y1)) / Dx**2
  vout= ( v(ix,iy)*(X2-x)*(Y2-y) + &
          v(ix+1,iy)*(x-X1)*(Y2-y)+&
          v(ix,iy+1)*(X2-x)*(y-Y1)+&
          v(ix+1,iy+1)*(x-X1)*(y-Y1)) / Dx**2

end subroutine bilint


end




