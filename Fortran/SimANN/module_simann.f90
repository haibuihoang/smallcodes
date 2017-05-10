!
! a simple 3 layers (1 hidden layer) neural network
! (c) 2011 Bui Hoang Hai
! Features  
!    - Autoscale input/target
!     
!
    module simann

	implicit none
    integer ::  Nin, Nhidden, Nout, iin, ihid, iout
	integer :: i,j,k

    Real,allocatable, Dimension(:,:) :: weight_ih, weight_ho, moment_ih, moment_ho
    Real,allocatable, Dimension(:) :: fin, fhid, fout, ftarget,  delta_o, delta_h

! Training
    Real	:: learning_rate,momentum, error

    logical	:: debug


    Contains

!---------------------------------------------
    subroutine ann_init
	real tmp
	write(*,*) 'Initiate the neural network'
	write(*,*) 'Number of input  neurons:',Nin
	write(*,*) 'Number of hidden neurons:',Nhidden
	write(*,*) 'Number of output neurons:',Nout
	
	allocate(weight_ih(0:Nin,Nhidden))
	allocate(weight_ho(0:Nhidden,Nout))
	allocate(moment_ih(0:Nin,Nhidden))
	allocate(moment_ho(0:Nhidden,Nout))
	
	! init weight that takes values randomly from [-1,1]
	do iin=0,Nin
	  do ihid=1,Nhidden
	    call RANDOM_NUMBER(tmp) 
	    weight_ih(iin,ihid)= 2. - tmp
	  enddo
	enddo
	do ihid=0,Nhidden
	  do iout=1,Nout
	    call RANDOM_NUMBER(tmp)
	    weight_ho(ihid,iout)= 2. - tmp
	  enddo
	enddo
	moment_ho = 0.
	moment_ih = 0.

	allocate(fin(Nin))
	allocate(fhid(0:Nhidden))
    allocate(fout(Nout))
    allocate(ftarget(Nout))
    allocate(delta_o(Nout))
    allocate(delta_h(Nhidden))
!---bias
	fin(0)=1.
	fhid(0)=1.
	open(111,file='debug.txt')
    end subroutine ann_init

!---------------------------------------------
    subroutine ann_forward
    real :: sum

! Tinh dau ra cua lop an
	do ihid=1,Nhidden
	   sum=0.
	   do iin=0,Nin
	     sum = sum + fin(iin)*weight_ih(iin,ihid)
	   enddo
       fhid(ihid) = sigma(sum)
	enddo

! Tinh dau ra cua lop dau ra
	do iout=1,Nout
	   sum =0.
	   do ihid=0,Nhidden
	     sum = sum + fhid(ihid)*weight_ho(ihid,iout)
	   enddo
       fout(iout) = sigma(sum)
	enddo

	if (debug) then
	   write(111,*)'FORWARD'
	   write(111,*)'..fin:',fin
	   write(111,*)'..weight_ih:'
	   do j=1,Nhidden
	      write(111,*)(weight_ih(i,j),i=0,Nin)
	   enddo
	   write(111,*)'..fhid',fhid
	   write(111,*)'..weight_ho',weight_ho
	   write(111,*)'..fout',fout
	endif

    end subroutine ann_forward

!---------------------------------------------
!--- NOTE that you have already run forward_prop
!--- Require ftarget
    subroutine ann_compute_error
       error=0.
       do iout=1,Nout
	      delta_o(iout)=(ftarget(iout)-fout(iout))
	      error=error+(delta_o(iout)**2)/Nout
       enddo
    end subroutine ann_compute_error

!---------------------------------------------
!---NOTE that you have already run forward_prop
!---Update weights for output layer
    subroutine ann_backward
    real :: delta_ih,delta_ho

!---Calculate error, delta and update weights for the output layers
    error=0.
    do iout=1,Nout
   	   delta_o(iout)=(ftarget(iout)-fout(iout))
	   error=error+(delta_o(iout)**2)/Nout						!--Mean square error
       delta_o(iout)=delta_o(iout)*fout(iout)*(1-fout(iout))	!--Derivative of Sigma:  y' = y(1-y)  
   	   do ihid=0,Nhidden
	      delta_ho=learning_rate*delta_o(iout)*fhid(ihid)
		  weight_ho(ihid,iout) = weight_ho(ihid,iout) + delta_ho + moment_ho(ihid,iout)
!		  weight_ho(ihid,iout) = weight_ho(ihid,iout) + delta_ho 
	      moment_ho(ihid,iout)=momentum*delta_ho	            !--momentum for next itteration
	   enddo
    enddo



!---Update weights for hidden layer
    do ihid=1,Nhidden
	delta_h(ihid)=0.
	do iout=1,Nout
	    delta_h(ihid) = delta_h(ihid) + weight_ho(ihid,iout)*delta_o(iout)
	enddo
	delta_h(ihid) = delta_h(ihid)*fhid(ihid)*( 1- fhid(ihid))
	do iin=0,Nin
	    delta_ih = learning_rate*delta_h(ihid)*fin(iin)
	    weight_ih(iin,ihid)=weight_ih(iin,ihid) &
				+ delta_ih &
				+ moment_ih(iin,ihid)
!--momentum for next itteration
	    moment_ih(iin,ihid)=momentum*delta_ih
	enddo
    enddo


    if(debug)then
      write(111,*)'weight_ih',weight_ih
    endif

!---
    if(debug)then
      write(111,*)'BACK'
      write(111,*)'--fout',fout
      write(111,*)'--ftarget',ftarget
      write(111,*)'delta_o',delta_o
      write(111,*)'weight_ih',weight_ih
      write(111,*)'weight_ho',weight_ho
    endif
    end subroutine ann_backward
!---------------------------------------------    
!---Sigmoid function
    real function sigma(x)
    real x
	sigma=1./ (1. + exp(-x))
    end function sigma
!---------------------------------------------    

    subroutine ann_deinit
	deallocate(weight_ih)
	deallocate(weight_ho)
	deallocate(fin)
	deallocate(fhid)
	deallocate(fout)
	deallocate(ftarget)
	deallocate(delta_o)
	deallocate(delta_h)
	close(111)
    end subroutine ann_deinit

    End module simann
