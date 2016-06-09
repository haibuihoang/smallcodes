! -------------------------------
! a simple 3 layers (1 hidden layer) neural network
! (c) 2011 Bui Hoang Hai
! -----------------------------------
Program Simple_ANN
use controller

implicit none

debug = .false.


call controller_init

if (training_mode) then
   call controller_train_data
endif

write(*,*) "Writting output..."
call controller_output_data


call controller_deinit	

write(*,*) "All done!"


Stop
End