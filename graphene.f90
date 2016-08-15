program tight_bind
use constants
use n_store
use initialise
use uppertri_solve
implicit none 

!--------------------------------
! Physical Variable Decleration
!--------------------------------

! Operator and matrix deccleration
complex*16,dimension(:,:),allocatable :: H
complex*16,dimension(:),allocatable :: nn
integer,dimension(:,:),allocatable :: unit_cell,atom_pos
real(kind=dp),dimension(:),allocatable :: eval

! matrix parameters 
integer :: i,j,imax,jmax,no_cell

! physical variables of type real
real(kind=dp) :: a,k,t


!----------------------------------
!   Compuatational Variables
!----------------------------------

! loop counters + zheev variable
integer :: l,m,lw,row_count,status,unit_form,counter 

!----------------------------------
!      Files for data output
!---------------------------------- 

open(file='eigen.dat',unit=22,status='replace')
open(file='vis_cell.dat',unit=24,status='replace')
open(file='vis_cellempty.dat',unit=25,status='replace')

!----------------------------------
!    system parameters 
!----------------------------------

! Lattice constant (a), Wave Number (k), Hopping energy (t)  
a = 1.0_dp
!k = 0.0_dp
k = -pi
t = 2.7_dp


! unit cell structures  
! the system will be perodic in x direction

! decide lattice type: 1=square, 2=arm, 3=zig-zag

print*, 'enter ribbon type: square(0), graphene armchair (1) zig-zag (2)'
read*, unit_form  

!----------------------------------
! Define unit cell structures
!----------------------------------

! populate unit cell where: 
! 1 = atom in position i,j
! 0 = no atom present 


! Square lattice ribbon

if (unit_form.eq.0) then
  t = 1.0_dp
  imax = 1
  jmax = 100
 
  allocate (unit_cell(0:imax+1,0:jmax+1))
  unit_cell(:,:) = 1

end if


! Armchair Nanoribbon

if (unit_form.eq.1) then 

! m = number of atoms (number of rows) 

print*, 'length of chain (input + 3)'
read*, m

imax = 3 + m
jmax = 4
 
allocate (unit_cell(0:imax+1,0:jmax+1))

unit_cell(:,:) = 0

  ! 1 graphene molecule 
  unit_cell(1,2) = 1
  unit_cell(1,3) = 1
  
  unit_cell(2,1)  = 1
  unit_cell(2,4)  = 1
  
  unit_cell(3,2) = 1
  unit_cell(3,3) = 1
  ! Additional repetitons
  do i = 1,m
  j = 3+i
    if (modulo(i,2).eq.1) then  
      unit_cell(j,1)  = 1
      unit_cell(j,4)  = 1
    else
      unit_cell(j,2) = 1
      unit_cell(j,3) = 1
    end if
  end do

end if


! Zig zag nanoribbon

if (unit_form.eq.2) then

print*, 'input length of chain (integer > 0)'
read*, m

imax = 2*m 
jmax = 2

allocate (unit_cell(0:imax+1,0:jmax+1))

unit_cell(:,:) = 0

! each value of i will produce two carbon atoms whose relative
! positions will depend on i

! for i even : (i,1),(i+1,2)
! for i odd : (i,2),(i+1,1)
j = 1
do i = 1,m
    if (modulo(i,2).eq.1) then
      unit_cell(j,1) = 1 
      unit_cell(j+1,2) = 1       
    else if (modulo(i,2).eq.0) then
      unit_cell(j+1,1) = 1
      unit_cell(j,2) = 1 
    end if
    j = j + 2
end do 
end if

!----------------------------------
!        sea of zero's
!----------------------------------
!top/bottem
unit_cell(0,:) = 0
unit_cell(imax+1,:) = 0

!left/right
unit_cell(:,0) = 0
unit_cell(:,jmax+1) = 0

! no_cell will count the number of atom/unit cell 
! this will determine the dimesnions of the hamiltonian

no_cell = sum(unit_cell)

!----------------------------------
!    Matrix Initilisation
!----------------------------------

! Allocate  Hamiltonian,Nearest Neighbour list and eigenenergy   
allocate (nn(1:no_cell),stat=status)
if (status.ne.0) stop 'Nearest Neigbour array failed to allocate in main'  

allocate (H(1:no_cell,1:no_cell),stat=status)
if (status.ne.0) stop 'Hamiltonian failed to allocate in main'  

allocate (eval(1:no_cell),stat=status)
if (status.ne.0) stop 'Eigenvalue array failed to allocate in main' 

allocate (atom_pos(1:imax,1:jmax),stat=status)
if (status.ne.0) stop 'Atom_pos array failed to allocate in main'

nn(:) = 0
H(:,:) = 0
eval(:) = 0
atom_pos(:,:) = 0
row_count = 0
counter = 0

!--------------------------------
!  run through zheev once to get
!  optimal parameters     
!--------------------------------

call inital(no_cell,no_cell,lw)

!---------------------------------------
! label occupied sites within unit cell
!---------------------------------------

do i = 1,imax
  do j = 1,jmax
    if (unit_cell(i,j).eq.1) then 
      counter = counter + 1
      atom_pos(i,j) = counter
    end if
  end do
end do 

!--------------------------------
!         Main program
!--------------------------------
do l = 1,200+1
  ! run through each occupied site in the unit cell
  do i = 1,imax
    do j = 1,jmax
      if (unit_cell(i,j).eq.1) then
        row_count = row_count + 1 
        call n_list(i,j,imax,jmax,unit_cell,a,k,t,unit_form,atom_pos,nn)
        H(row_count,:) = nn(:)
      end if 
    end do
  end do
  
  ! solve the generated Hamiltonian 

  call diag_H(no_cell,no_cell,lw,H,eval)

  write(22,*) k,eval
  
  k = k + pi/100.0_dp
  
  ! reset system state 
  H(:,:) = 0
  eval(:) = 0
  row_count = 0
  
end do

! output unit cell structure 

do i = 0,imax+1
  do j = 0,jmax+1
    if (unit_cell(i,j).eq.1) then 
      write(24,*) j,i
    else
      write(25,*) j,i
    end if 
  end do 
end do 

!-----------------------------------
! Deallocate arrays and close files
!-----------------------------------
deallocate(H)
deallocate(eval)
deallocate(nn)
deallocate(unit_cell)

close(22)
close(24)
close(25)

print*, 'max num of evals is: ',no_cell

end program 
