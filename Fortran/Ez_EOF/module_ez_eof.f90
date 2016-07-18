!
!  Easy (EZ) EOF/PCA analysis module
!  Requirement: LAPACK
!  Firstversion: 16/7/2016 (c) Bui Hoang Hai, Hanoi University of Science/Kyoto University
!
!  There are just 3 simple subrountines to call with few arguments:
!  Note the we can use the full field, the anomaly will be calculated 
!  
!  1) calculate_EOFs(X,Ns,Nt) will initial all dimension & variables
!  and calculate, and store all neccesary fields which can be use later:
!       e_Mean: climatology/mean field
!       e_Xprime: anomaly = input X - e_Xmean
!       e_Cov: Covariance Matrix of X
!       e_EiVal:  Eigenvalues of e_Cov
!       e_EiVec:  Eigenvectors of e_Cov, i.e. EOFs
!       e_PCs:  Principal Components
!
!   2) Restore_Fields(X,np) will restore to X using np (<=Ns) principal componet
!   3) Deinit_EOFs() will deallocate all fields after use
!
!   


MODULE ez_eof
implicit none

integer :: e_Nt, e_Ns
real, dimension(:,:), allocatable :: e_Cov, e_Xprime, e_EiVec, e_PC
real, dimension(:), allocatable :: e_Xmean, e_EiVal

CONTAINS


! Ns: space dimension
! Nt: time dimension
! X 
subroutine calculate_EOFs(X,Nt,Ns)
integer :: Ns, Nt,i
real, dimension (Nt,Ns) ::  X
real, dimension(:,:),allocatable :: XprimeT, A
real, dimension(:), allocatable :: WORK,W
integer :: LWORK, INFO
e_Nt=Nt
e_Ns=Ns
allocate(e_Xprime(Nt,Ns))  ! Anomaly
allocate(e_PC(Nt,Ns))      ! Principal components
allocate(e_Xmean(Ns))      ! Mean
allocate(e_Cov(Ns,Ns))     ! Covarance Matrix
allocate(e_EiVal(Ns))      ! Eigen value
allocate(e_EiVec(Ns,Ns))   ! Eigen vectors

!caculate mean
do i=1,Ns
   e_Xmean(i)=SUM(X(:,i))/Nt   
enddo

!caculate Xprime
do i=1,Ns
   e_Xprime(:,i) = X(:,i) - e_Xmean(i)
enddo

!calculate the covariance matrix
allocate(XprimeT(Ns,Nt)) 
XprimeT = Transpose(e_Xprime)  !Transpose X: Kxn
e_Cov = Matmul(XPrimeT,e_Xprime)/Real(Nt-1) 
deallocate(XPrimeT)

!Find the eigen vectors and eigen value using Lapack subroutines
allocate(A(Ns,Ns)); 
A = e_Cov  !arg_A will hold the Eigen values!

allocate(W(Ns))

LWORK = 3*Ns - 1
allocate(WORK(1:LWORK))  !just for working array, does need after calling SSYEV

call SSYEV("V","L",Ns, A, Ns, W, WORK, LWORK, INFO)

if (INFO/=0) then
  write(*,*)"SSYEV Failed: ",INFO
endif
deallocate(WORK)

!Reorder the Eigen values and EOFs in Descending order
e_EiVec(:,1:Ns) = A(:,Ns:1:-1)
e_EiVal(1:Ns) = W(Ns:1:-1)
deallocate(A)
deallocate(W)

!Calculate Principle components according to EOFs
e_Pc = Matmul(e_Xprime,e_Eivec)


end subroutine calculate_EOFs


!X --> field to restore
!np --> number of PCs to be process, np<=Ns
subroutine restore_Fields(X,np)
   integer, intent(in) :: np
   integer :: i,npp
   real, dimension(e_Nt,e_Ns) :: X
   real, dimension(:,:), allocatable :: evecT

   npp=np
   if (np .ge. e_Ns) npp=e_Ns

   allocate(evecT(npp,e_Ns))

   evecT = Transpose(e_EiVec(:,1:npp))

   X = matmul(e_PC(:,1:npp),evecT)
   !Add the mean to restored full field
   do i=1,e_Ns
      X(:,i) = X(:,i) + e_Xmean(i)
   enddo
   
   deallocate(evecT)
end subroutine restore_Fields


subroutine deinit_EOFs()
deallocate(e_Xprime,e_PC)
deallocate(e_Xmean)
deallocate(e_Cov)
deallocate(e_EiVal,e_EiVec)
end subroutine 


END MODULE ez_eof
