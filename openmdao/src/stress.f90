! Computes the stresses given a displacement field
! Dr. John T. Hwang
! Sicheng He
! June, 2016

subroutine getstresses(num_nodes, num_elems, &
     E, nodes, elems, disp, stress)

  implicit none

  !f2py intent(in) num_nodes, num_elems, E, nodes, elems, disp
  !f2py intent(out) stress
  !f2py depend(num_nodes) nodes, disp
  !f2py depend(num_elems) elems, stress)

  ! Input
  integer, intent(in) :: num_nodes, num_elems
  double precision, intent(in) :: E, nodes(num_nodes, 3)
  integer, intent(in) :: elems(num_elems, 2)
  double precision, intent(in) :: disp(num_nodes, 3)

  ! Output
  double precision, intent(out) :: stress(num_elems)

  ! Working
  double precision :: xyz1(3), xyz2(3), uvw1(3), uvw2(3), u, L
  integer :: ielem

  do ielem = 1, num_elems
     xyz1 = nodes(elems(ielem, 1), :)
     xyz2 = nodes(elems(ielem, 2), :)
     uvw1 = disp(elems(ielem, 1), :)
     uvw2 = disp(elems(ielem, 2), :)
     L = sqrt(dot_product(xyz2 - xyz1, xyz2 - xyz1))

     u = dot_product(xyz2 - xyz1, uvw2 - uvw1) / L
     stress(ielem) = E * u / L
  end do
  
end subroutine getstresses



subroutine getstressder(num_nodes, num_elems, nnz, &
     E, nodes, elems, data, rows, cols)

  implicit none

  !f2py intent(in) num_nodes, num_elems, nnz, E, nodes, elems
  !f2py intent(out) data, rows, cols
  !f2py depend(num_nodes) nodes
  !f2py depend(num_elems) elems
  !f2py depend(nnz) data, rows, cols
  
  ! Input
  integer, intent(in) :: num_nodes, num_elems, nnz
  double precision, intent(in) :: E, nodes(num_nodes, 3)
  integer, intent(in) :: elems(num_elems, 2)

  ! Output
  double precision, intent(out) :: data(nnz)
  integer, intent(out) :: rows(nnz), cols(nnz)

  ! Working
  double precision :: xyz1(3), xyz2(3), L
  integer :: ielem, k, index

  data(:) = 0.
  rows(:) = 0
  cols(:) = 0

  index = 0
  do ielem = 1, num_elems
     xyz1 = nodes(elems(ielem, 1), :)
     xyz2 = nodes(elems(ielem, 2), :)
     L = sqrt(dot_product(xyz2 - xyz1, xyz2 - xyz1))

     do k = 1, 3
        index = index + 1
        data(index) = E / L ** 2 * (xyz1(k) - xyz2(k))
        rows(index) = ielem
        cols(index) = 3 * (elems(ielem, 1) - 1) + k
        
        index = index + 1
        data(index) = E / L ** 2 * (xyz2(k) - xyz1(k))
        rows(index) = ielem
        cols(index) = 3 * (elems(ielem, 2) - 1) + k
     end do
  end do

  if (index .ne. nnz) then
     print *, 'Error in getstressder: did not reach end of nnz vectors'
  end if

  rows(:) = rows(:) - 1
  cols(:) = cols(:) - 1

end subroutine getstressder
