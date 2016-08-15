module uppertri_solve
use constants

implicit none 

contains

subroutine diag_H(N,LDA,LW,A,E)

complex*16,intent(inout) :: A(:,:)
integer,intent(in) :: N,LDA,LW
real(kind=dp),allocatable,intent(out) :: E(:)

complex*16,allocatable :: work(:)
real(kind=dp),allocatable :: RW(:)

integer :: info,status

allocate(work(LW),stat=status)
if (status.ne.0) stop 'cannot allocate work array in solve'

allocate(E(N),stat=status)
if (status.ne.0) stop 'cannot allocate eigenvalue array in solve'

allocate(RW(5*N),stat=status) 
if (status.ne.0) stop 'cannot allocate RW array in solve'

call zheev('N','U',N,A,LDA,E,work,LW,RW,info)

if (info.ne.0) print*, 'diag didnt exit properly routine'

deallocate(work)
deallocate(RW)

end subroutine 

end module 

