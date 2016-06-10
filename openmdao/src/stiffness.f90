! Computes the matrix for the finite-element linear system
! Dr. John T. Hwang
! Sicheng He
! June, 2016

subroutine getmtx(num_nodes, num_elems, num_cons, nnz, &
     E, nodes, elems, areas, cons, data, rows, cols)

  implicit none

  !f2py intent(in) num_nodes, num_elems, num_cons, nnz, E, nodes, elems, areas, cons
  !f2py intent(out) data, rows, cols
  !f2py depend(num_nodes) nodes
  !f2py depend(num_elems) elems, areas
  !f2py depend(num_cons) cons
  !f2py depend(nnz) data, rows, cols

  ! Input
  integer, intent(in) :: num_nodes, num_elems, num_cons, nnz
  double precision, intent(in) :: E, nodes(num_nodes, 3)
  integer, intent(in) :: elems(num_elems, 2)
  double precision, intent(in) :: areas(num_elems)
  integer, intent(in) :: cons(num_cons)

  ! Output
  double precision, intent(out) :: data(nnz)
  integer, intent(out) :: rows(nnz), cols(nnz)

  ! Working
  double precision :: Telem(2, 6), Kelem(2, 2)
  double precision :: xyz1(3), xyz2(3), A, L, cos_xyz(3)
  double precision :: data6x6(6, 6)
  integer :: rows6x6(6, 6), cols6x6(6, 6)
  integer :: ones11(6, 6), ones12(6, 6), ones21(6, 6), ones22(6, 6)
  integer :: k, k1, k2, ielem, icons, index, ind1, ind2

  Kelem(1, 1) =  1.
  Kelem(1, 2) = -1.
  Kelem(2, 1) = -1.
  Kelem(2, 2) =  1.

  do k2 = 1, 6
     do k1 = 1, 6
        rows6x6(k1, k2) = mod(k1-1, 3)
        cols6x6(k1, k2) = mod(k2-1, 3)
     end do
  end do

  ones11(:, :) = 0
  ones12(:, :) = 0
  ones21(:, :) = 0
  ones22(:, :) = 0
  
  ones11(1:3, 1:3) = 1
  ones12(1:3, 4:6) = 1
  ones21(4:6, 1:3) = 1
  ones22(4:6, 4:6) = 1

  Telem(:, :) = 0.

  data(:) = 0.
  rows(:) = 0
  cols(:) = 0

  index = 0
  do ielem = 1, num_elems
     xyz1 = nodes(elems(ielem, 1), :)
     xyz2 = nodes(elems(ielem, 2), :)
     L = sqrt(dot_product(xyz2 - xyz1, xyz2 - xyz1))
     A = areas(ielem)

     cos_xyz = (xyz2 - xyz1) / L
     Telem(1, 1:3) = cos_xyz
     Telem(2, 4:6) = cos_xyz

     data6x6(:, :) = matmul(matmul(transpose(Telem), Kelem), Telem)
     data6x6 = data6x6 * E * A / L

     ind1 = 3 * (elems(ielem, 1)-1)
     ind2 = 3 * (elems(ielem, 2)-1)

     do k1 = 1, 6
        do k2 = 1, 6
           index = index + 1
           data(index) = data(index) + data6x6(k1, k2)
           rows(index) = rows(index) + rows6x6(k1, k2) + 1
           cols(index) = cols(index) + cols6x6(k1, k2) + 1

           rows(index) = rows(index) + ones11(k1, k2) * ind1
           cols(index) = cols(index) + ones11(k1, k2) * ind1
           
           rows(index) = rows(index) + ones12(k1, k2) * ind1
           cols(index) = cols(index) + ones12(k1, k2) * ind2
           
           rows(index) = rows(index) + ones21(k1, k2) * ind2
           cols(index) = cols(index) + ones21(k1, k2) * ind1
           
           rows(index) = rows(index) + ones22(k1, k2) * ind2
           cols(index) = cols(index) + ones22(k1, k2) * ind2
        end do
     end do
  end do

  do icons = 1, num_cons
     do k = 1, 3
        index = index + 1
        data(index) = 1.
        rows(index) = 3 * (cons(icons)-1) + k
        cols(index) = 3 * num_nodes + 3 * (icons-1) + k
        
        index = index + 1
        data(index) = 1.
        cols(index) = 3 * (cons(icons)-1) + k
        rows(index) = 3 * num_nodes + 3 * (icons-1) + k
     end do
  end do

  if (index .ne. nnz) then
     print *, 'Error in getmtx: did not reach end of nnz vectors'
  end if

  rows(:) = rows(:) - 1
  cols(:) = cols(:) - 1

end subroutine getmtx




subroutine getresder(num_nodes, num_elems, num_aug, nnz, &
     E, nodes, elems, disp_aug, data, rows, cols)

  implicit none

  !f2py intent(in) num_nodes, num_elems, num_aug, nnz, E, nodes, elems, disp_aug
  !f2py intent(out) data, rows, cols
  !f2py depend(num_nodes) nodes
  !f2py depend(num_elems) elems
  !f2py depend(num_aug) disp_aug
  !f2py depend(nnz) data, rows, cols

  ! Input
  integer, intent(in) :: num_nodes, num_elems, num_aug, nnz
  double precision, intent(in) :: E, nodes(num_nodes, 3)
  integer, intent(in) :: elems(num_elems, 2)
  double precision, intent(in) :: disp_aug(num_aug)

  ! Output
  double precision, intent(out) :: data(nnz)
  integer, intent(out) :: rows(nnz), cols(nnz)

  ! Working
  double precision :: Telem(2, 6), Kelem(2, 2)
  double precision :: xyz1(3), xyz2(3), L, cos_xyz(3)
  double precision :: data6x6(6, 6)
  integer :: rows6x6(6, 6), cols6x6(6, 6)
  integer :: ones11(6, 6), ones12(6, 6), ones21(6, 6), ones22(6, 6)
  integer :: k1, k2, ielem, index, ind1, ind2, ind_aug

  Kelem(1, 1) =  1.
  Kelem(1, 2) = -1.
  Kelem(2, 1) = -1.
  Kelem(2, 2) =  1.

  do k2 = 1, 6
     do k1 = 1, 6
        rows6x6(k1, k2) = mod(k1-1, 3)
        cols6x6(k1, k2) = mod(k2-1, 3)
     end do
  end do

  ones11(:, :) = 0
  ones12(:, :) = 0
  ones21(:, :) = 0
  ones22(:, :) = 0
  
  ones11(1:3, 1:3) = 1
  ones12(1:3, 4:6) = 1
  ones21(4:6, 1:3) = 1
  ones22(4:6, 4:6) = 1

  Telem(:, :) = 0.

  data(:) = 0.
  rows(:) = 0

  index = 0
  do ielem = 1, num_elems
     xyz1 = nodes(elems(ielem, 1), :)
     xyz2 = nodes(elems(ielem, 2), :)
     L = sqrt(dot_product(xyz2 - xyz1, xyz2 - xyz1))

     cos_xyz = (xyz2 - xyz1) / L
     Telem(1, 1:3) = cos_xyz
     Telem(2, 4:6) = cos_xyz

     data6x6(:, :) = matmul(matmul(transpose(Telem), Kelem), Telem)
     data6x6 = data6x6 * E / L

     ind1 = 3 * (elems(ielem, 1)-1)
     ind2 = 3 * (elems(ielem, 2)-1)

     do k1 = 1, 6
        do k2 = 1, 6
           index = index + 1
           ind_aug = 0
           data(index) = data(index) + data6x6(k1, k2)
           rows(index) = rows(index) + rows6x6(k1, k2) + 1
           ind_aug = ind_aug + cols6x6(k1, k2) + 1

           rows(index) = rows(index) + ones11(k1, k2) * ind1
           ind_aug = ind_aug + ones11(k1, k2) * ind1
           
           rows(index) = rows(index) + ones12(k1, k2) * ind1
           ind_aug = ind_aug + ones12(k1, k2) * ind2
           
           rows(index) = rows(index) + ones21(k1, k2) * ind2
           ind_aug = ind_aug + ones21(k1, k2) * ind1
           
           rows(index) = rows(index) + ones22(k1, k2) * ind2
           ind_aug = ind_aug + ones22(k1, k2) * ind2

           data(index) = data(index) * disp_aug(ind_aug)

           cols(index) = ielem
        end do
     end do
  end do

  if (index .ne. nnz) then
     print *, 'Error in getmtx: did not reach end of nnz vectors'
  end if

  rows(:) = rows(:) - 1
  cols(:) = cols(:) - 1

end subroutine getresder
