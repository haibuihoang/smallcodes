!
! (c) 2011 Bui Hoang Hai
! Features  
!    - Autoscale input/target
!     
!----------------------------------
! Options that will be read from namelist.txt
! The input file is as follows:
!      Ninput  Noutput  Nsample
!	   input1.1 input2.1 ... inputN.1   output1.1 ... outputM.1
!	   input1.2 input2.2 ... inputN.2   output1.2 ... outputM.2
!      ....
! The save_file is created after the trainning process completed, 
!   the save_file contains scaling factors and connection weights
!   if training_mode=False, then save_file  is readed
!   Note: the number of Input, Output in input_file must be equal to those in save_file
!-------------------------------
! When trainning mode is turned on, data is calculated and scale as follows:
!    inputs [min,max] <--> fin [-1,1]
!    outputs[min,max] <--> fout [0.2,0.8]
!    
!	 we will emply a subroutine:  rescale(x,y,x1,x2,y1,y2) that will calc y from x such that [x1,x2]-->[y1,y2]
!
!-------------------------------

module controller

  Use simann

  implicit none

  character (len=30) :: input_file, output_file, save_file
  logical :: training_mode

  ! Use for calculating scaling factors
  Integer :: Nsample, Nepoc

  Real,allocatable, Dimension(:,:) :: Inputs, Outputs
  Real,allocatable, Dimension(:) :: min_in, max_in, min_out, max_out
  Real :: mean_error

 Contains

!---------------------------------------------    
!Read info from namelist
  subroutine controller_init
    namelist /info/ input_file, output_file, save_file, Nhidden, learning_rate,momentum, training_mode,Nepoc

    open(1,file='namelist.txt',status='unknown')
    read(1,info)
    write(*,info)
    close(1)

	!--Init data for input file
	!--Depend on training_mode = T/F
	open(1,file=input_file,status='old')
	read(1,*) Nin, Nout, Nsample

	call ann_init	   ! Init ann first

	allocate(Inputs(Nin,NSample));allocate(Outputs(Nout,NSample))
	allocate(min_in(Nin));allocate(max_in(Nin))
	allocate(min_out(Nout));allocate(max_out(Nout))


	if (training_mode) then
      write(*,*) "Reading data..."     
	  do j=1,NSample
	    read(1,*)(Inputs(i,j),i=1,Nin), (Outputs(i,j),i=1,Nout)
	  enddo       
	  do i=1,Nin
		 min_in(i)=MINVAL(Inputs(i,:))
	     max_in(i)=MAXVAL(Inputs(i,:))
		 write(*,*)i," in:[",min_in(i),'->',max_in(i),']'
	  enddo
	  do i=1,Nout
		 min_out(i)=MINVAL(Outputs(i,:))
	     max_out(i)=MAXVAL(Outputs(i,:))
		 write(*,*)i," out:[",min_out(i),'->',max_out(i),']'
	  enddo
	   	 
	else

	  do j=1,NSample
	    read(1,*)(Inputs(i,j),i=1,Nin)
	  enddo       

	  call controller_load_info
	  call controller_print_info
    endif
    
	close(1)

  end subroutine controller_init

!---------------------------------------------    
  subroutine controller_train_data
     integer epoc

     do epoc=1,Nepoc
       mean_error=0.
	   do j=1, Nsample
	     do i=1,Nin
		    Fin(i) = rescale(Inputs(i,j),min_in(i),max_in(i),-1.,1.)
		 enddo
		 Call ann_forward


	 	 do i=1,Nout
		    Ftarget(i) = rescale(Outputs(i,j),min_out(i),max_out(i),0.2,0.8)
		 enddo	 
		 Call ann_backward
		 mean_error=mean_error+error/Nsample
	   enddo
	   if (mod(epoc,50).eq.0) then
	      write(*,*) "epoc=",epoc,", error=",mean_error
       endif
	 enddo

	 call controller_print_info
	 call controller_save_info 
  end subroutine controller_train_data


!---------------------------------------------    
  subroutine controller_output_data
     
     open(1,file=output_file,status='unknown')

	 write(1,*) Nin, Nout, Nsample
     
	 do j=1, Nsample
	    do i=1,Nin
	       Fin(:) = rescale(Inputs(i,j),min_in(i),max_in(i),-1.,1.)
	    enddo
	    Call ann_forward
	 	do i=1,Nout
		    Ftarget(i) = rescale(Fout(i),0.2,0.8,min_out(i),max_out(i))
		enddo	 
		write(1,*) (Inputs(i,j),i=1,Nin),(Ftarget(i),i=1,Nout)
	 enddo	 

  end subroutine controller_output_data



!---------------------------------------------    
  subroutine controller_print_info
	write(*,*)"Input,Hidden,Output"
	write(*,*) Nin,Nhidden,Nout
	write(*,*)'Weight_input-hidden:'
	do i=0,Nin
	  write(*,*) (weight_ih(i,j),j=1,Nhidden)
	enddo
	write(*,*)'Weight_hidden-output:'
	do i=0,Nhidden
	  write(*,*) (weight_ho(i,j),j=1,Nout)
	enddo
    write(*,*)"min_in:"
    write(*,*) min_in
    write(*,*)"Max_in:"
    write(*,*) max_in
    write(*,*)"min_out:"
    write(*,*) min_out
    write(*,*)"Max_out:"
    write(*,*) max_out
  end subroutine controller_print_info
!---------------------------------------------


!---Save weights and scaling factors
  subroutine controller_save_info
    open(1,file=save_file,status='unknown')
	write(1,*)"Input,Hidden,Output"
	write(1,*) Nin,Nhidden,Nout
	write(1,*)'Weight_input-hidden:'
	do i=0,Nin
	  write(1,*) (weight_ih(i,j),j=1,Nhidden)
	enddo
	write(1,*)'Weight_hidden-output:'
	do i=0,Nhidden
	  write(1,*) (weight_ho(i,j),j=1,Nout)
	enddo
    write(1,*)"min_in:"
    write(1,*) min_in
    write(1,*)"Max_in:"
    write(1,*) max_in
    write(1,*)"min_out:"
    write(1,*) min_out
    write(1,*)"Max_out:"
    write(1,*) max_out
  end subroutine controller_save_info

!---------------------------------------------
  subroutine controller_load_info
    integer :: Nin1,Nout1
    open(1,file=save_file,status='old')
	read(1,*)
	read(1,*) Nin1,Nhidden,Nout1
	if ((Nin.ne.Nin1).or.(Nout.ne.Nout1)) then
	   write(*,*) "Number of inputs,outputs in save file are not equal to those in input file!"
	   stop
	endif

	read(1,*)
	do i=0,Nin
	  read(1,*) (weight_ih(i,j),j=1,Nhidden)
	enddo
	read(1,*)
	do i=0,Nhidden
	  read(1,*) (weight_ho(i,j),j=1,Nout)
	enddo
    read(1,*)
    read(1,*) min_in
    read(1,*)
    read(1,*) max_in
    read(1,*)
    read(1,*) min_out
    read(1,*)
    read(1,*) max_out
  end subroutine controller_load_info

!---------------------------------------------

!---------------------------------------------
 subroutine controller_deinit
	deallocate(Inputs);deallocate(Outputs)
	deallocate(min_in);deallocate(max_in)
	deallocate(min_out);deallocate(max_out)
    call ann_deinit
 end subroutine controller_deinit


!---------------------------------------------
!calc y from x such that [x1,x2]-->[y1,y2]
 real function rescale(x,x1,x2,y1,y2)
    real x,x1,x2,y1,y2
    rescale = y1 + (y2-y1)*(x-x1)/(x2-x1)
 end function	rescale
!---------------------------------------------

 End module controller