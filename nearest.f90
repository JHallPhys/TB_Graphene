module n_store
use constants

implicit none 

contains 

subroutine n_list(i,j,imax,jmax,unit_cell,a,k,t,unit_form,atom_pos,nn)
implicit none 

integer,intent(in) :: imax,jmax,i,j,unit_form
real(kind=dp),intent(in) :: a,k,t

integer,dimension(0:imax+1,0:jmax+1),intent(in) :: unit_cell
integer,dimension(1:imax,1:jmax),intent(in) :: atom_pos
complex*16,intent(inout) :: nn(:) 


complex*16 :: z,t_c,psi,psi_plus
integer :: lbound_check,h_index,rbound_check,comp_check


! complex elements of the hamiltonian for the upper tri region
 
z = cmplx(cos(k*a),-sin(k*a)) ! form of wave moving in neg x direction
t_c = cmplx(-t,0) ! hopping energy
psi = t_c*(1.0_dp+z) ! element for a neg wave leaving unit cell along x
psi_plus = t_c*(1.0_dp+conjg(z)) 
comp_check = 0
! Hamiltonian will contain elements = 0 so set to zero to save time
nn(:) = 0

!---------------------------------------------------
! find nearest neighbours by moving along unit cell
!---------------------------------------------------


  ! left periodic boundary
  !print*, 'left bound'
  lbound_check = j
  if (lbound_check-1.lt.1) then
    lbound_check = jmax
  else 
    lbound_check = j-1
  end if

  ! right periodic boundary

  rbound_check = j
  if (rbound_check+1.gt.jmax) then
    rbound_check = 1
  else 
    rbound_check = j+1
  end if
!------------------------------
!      movement in x
!------------------------------


  ! x - axis  

  ! check atoms at left and right sites exist 
  if ( (unit_cell(i,rbound_check).eq.1).and.(unit_cell(i,lbound_check).eq.1) ) then 
 
    ! Does the electron cross the left boundary 
    if ( lbound_check.eq.jmax) then
      if (j+1.eq.jmax) then
        h_index = atom_pos(i,j+1)
        nn(h_index) =  psi
      else 
        h_index = atom_pos(i,j+1)
        nn(h_index) = t_c
        h_index = atom_pos(i,jmax)
        nn(h_index) = t_c*z
      end if
    end if

    ! only j+1 will lie in the upper triangualr region  
    if ( (j+1.le.jmax).and.(j-1.ge.1) ) then
      h_index = atom_pos(i,j+1)
      nn(h_index) = t_c
    end if
  end if

  ! does one atom exist to the right or left only

  ! right only 
  if ( (unit_cell(i,rbound_check).eq.1).and.(unit_cell(i,lbound_check).eq.0) ) then
    if (rbound_check.gt.j) then 
      h_index = atom_pos(i,rbound_check)
      nn(h_index) = t_c  
    end if

  ! left only
  else if ( (unit_cell(i,rbound_check).eq.0).and.(unit_cell(i,lbound_check).eq.1) ) then
    if (lbound_check.gt.j) then
      h_index = atom_pos(i,lbound_check)
      nn(h_index) = t_c*z
    end if
  end if
  
  ! x' axis

  ! nearest diaganol neigbour within 1 step along  x' axis but lays in upper diag 
  if (unit_form.ge.1) then
    if ( ( j+1.gt.jmax).and.(atom_pos(i+1,rbound_check).gt.atom_pos(i,j)) ) then
      h_index = atom_pos(i+1,rbound_check)
      nn(h_index) = psi_plus
    else if  ( ( j-1.lt.1 ).and.(atom_pos(i+1,lbound_check).gt.atom_pos(i,j)) ) then
      h_index = atom_pos(i+1,lbound_check)
      nn(h_index) = psi 
      comp_check = 1
    else if ( (unit_cell(i+1,j-1).eq.1) ) then
      h_index = atom_pos(i+1,j-1)
      nn(h_index) = t_c   
    end if
  end if

!------------------------------
       !movement in y
!------------------------------

  ! y-axis 

  ! check atom below exists
  if (unit_cell(i+1,j).ne.0) then 
    h_index = atom_pos(i+1,j)
    nn(h_index) = t_c
  end if

  ! y' axis

  ! check atom exists in in neg y' direction exists
 
  if (comp_check.eq.0) then
    if (unit_form.ge.1) then
      if ( (unit_cell(i+1,j+1).eq.1) ) then      
        h_index = atom_pos(i+1,j+1)
        nn(h_index) = t_c
      end if
    end if
  end if 
  
end subroutine 

end module 
