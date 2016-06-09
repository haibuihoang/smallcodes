! Test Dicret Fourier Transform
implicit none
integer :: i,j,k,u
integer, parameter :: N = 128
real,dimension(N) :: fx, fur,fui
real,parameter :: pi= 3.14159265359

!init function fx with some sin function
do i=1,N/2
  fx(i) = 5 + 5*cos(2*pi*i/N) + 3*sin(2*pi*i*3/N) + &
          2*sin(2*pi*i*6/N) + i/2 + 5*sin(2*pi*20*i/N) 
enddo


!write test file
open(1,file='fx.dat',access='direct',form='unformatted',recl=N*4)
write(1,rec=1)fx

call dft(fx,N,fur,fui)

write(1,rec=2)fur
write(1,rec=3)fui

call restore(fur,fui,N,fx)

write(1,rec=4)fx

!lowpass

fur(10:N)=0
fui(10:N)=0
call restore(fur,fui,N,fx)
write(1,rec=5)fx

close(1)

stop
end

!caculculate fourier function
subroutine dft(fx,n,fur,fui)
  implicit none
  integer :: n,k,l
  real :: x,u
  real, dimension(n) :: fx, fur, fui
  real, parameter :: pi=3.141569265359
  real :: p  ! phases

  do l=1,n ! Determin fourier coeficents
    u=l-1
    fur(l)=0
    fui(l)=0
    do k=1,n  ! loop for each x=ix
      x=k-1  
      p = 2*pi*x*u/n
      fur(l) = fur(l) + fx(k)*cos(p)/N
      fui(l) = fui(l) - fx(k)*sin(p)/N
    enddo
    write(*,*)'u=',u,p,fur(l),fui(l)
  enddo
end subroutine dft

subroutine restore(fur,fui,n,fx)
implicit none
  integer :: n,l,k
  real, dimension(n) :: fx, fur, fui
  real, parameter :: pi=3.141569265359
  real :: p,u,x ,f ! phases
  
  do k=1,n
    x=k-1
    fx(k)=0
    f=0
    do l=1,n
      u=l-1
      p=2*pi*x*u/n
      fx(k)=fx(k)+(fur(l)*cos(p) - fui(l)*sin(p))
    enddo
    write(*,*)k,fx(k),f
  enddo
end subroutine restore


