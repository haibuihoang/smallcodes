use ez_eof
implicit none

integer, parameter :: Nx=2, Nt=10
real, dimension(Nt,Nx) :: Test, Restored

Test = reshape((/ 1,2,4,5,3,2,1,2,3,4, 2,4,5,6,4,5,2,2,4,7/), shape(Test))


call calculate_EOFs(Test,Nt,Nx)

write(*,*)"Data"
write(*,*)Test

write(*,*)"Mean"
write(*,*)e_Xmean

write(*,*)"Anonamlly Field"
write(*,*)e_Xprime

write(*,*)"Covariance matrix:"
write(*,*)e_Cov

write(*,*)"Eigen Values"
write(*,*)e_EiVal

write(*,*)"EOFs"
write(*,*)"1st:",e_EiVec(:,1)
write(*,*)"2nd:",e_EiVec(:,2)

write(*,*)"PCs"
write(*,*)"1st:",e_PC(:,1)
write(*,*)"2nd:",e_PC(:,2)


call restore_Fields(Restored,1)

write(*,*)"Restored:"
write(*,*)Restored

call deinit_EOFs()



stop
end
