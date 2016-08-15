module initialise 
use constants 
implicit none 

contains 

! run the ZHEEV routine with the array shapes 
! of the system and left work = -1 to obtain optimal 
! for l work and hence work array 

subroutine inital(N,LDA,LW)

integer,intent(in) :: N,LDA
integer,intent(out) :: LW

complex*16,allocatable :: A(:,:),work(:)
real(kind=dp),allocatable :: E(:),RW(:)

integer :: L=-1,info,status

! allocate arrays based on sizes used in system 
! under consideration 

allocate(A(LDA,N),stat=status)
if (status.ne.0) stop 'cannot allocate eigenvector array in init'

allocate(work(L),stat=status)
if (status.ne.0) stop 'cannot allocate work array in init'

allocate(E(N),stat=status)
if (status.ne.0) stop 'cannot allocate eigenvalue array in init'

allocate(RW(5*N),stat=status) 
if (status.ne.0) stop 'cannot allocate RW array in init'

call zheev('N','U',N,A,LDA,E,work,L,RW,info)

LW = int(work(1))

if (info.ne.0) print*, 'diag didnt exit properly routine'

deallocate (A)
deallocate(work)
deallocate(E)
deallocate(RW)

end subroutine 

end module 
