! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   TurboGAP is distributed in the hope that it will be useful for non-commercial
! HND X   academic research, but WITHOUT ANY WARRANTY; without even the implied
! HND X   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! HND X   ASL for more details.
! HND X
! HND X   You should have received a copy of the ASL along with this program
! HND X   (e.g. in a LICENSE.md file); if not, you can write to the original
! HND X   licensor, Miguel Caro (mcaroba@gmail.com). The ASL is also published at
! HND X   http://github.com/gabor1/ASL
! HND X
! HND X   When using this software, please cite the following reference:
! HND X
! HND X   Miguel A. Caro. Phys. Rev. B 100, 024112 (2019)
! HND X
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module decompose

  contains


!**************************************************************************
  subroutine dd_affinity_dense()
    implicit none

    integer, intent(out), allocatable :: ranks
    integer, intent(in) :: dd_grid(3)

    integer :: rank

    allocate( ranks(dd_grid(0):dd_grid(1):dd_grid(2)) )
    rank = 0
    do i = 1, dd_grid(0)
      do j = 1, dd_grid(1)
        do k = 1, dd_grid(2)
          ranks(i,j,k) = rank
          rank = rank + 1
        end do
      end do
    end do
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine dd_assign(grid_root, color, grid, method, ntasks)
    implicit none

    integer, intent(out), allocatable :: grid_root
    integer, intent(out), allocatable :: color
    integer, intent(in) :: grid(3)
    integer, intent(in) :: method
    integer, intent(in) :: ntasks

    integer :: rank
    integer :: step
    integer :: grid_size

    grid_size = grid(0) * grid(1) * grid(2)

    if (method == "cyclic") then
      ! first MPI ranks do MD
      step = 1
      wrap = grid_size
    else (method == "block") then
      ! spread MD ranks uniformly
      step = (ntasks - 1) / grid_size + 1
      wrap = ntasks
    endif

    allocate( color(1:ntasks) )
    do i = 1, ntasks
      color(i) = (((i - 1) / step) * step) % wrap + 1
    end do

    allocate( grid_root(1:grid(1), 1:grid(2), 1:grid(3)) )
    rank = 0
    do i = 1, grid(1)
      do j = 1, grid(2)
        do k = 1, grid(3)
          grid_root(i,j,k) = rank
          rank = rank + step
        end do
      end do
    end do
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine dd_placement(placement, grid, grid_root, grid_cells, ntasks, &
                          n_sites, positions, a_box, b_box, c_box)
    implicit none
    integer, intent(out), allocatable :: placement(:)
    integer, intent(in) :: grid(3)
    integer, intent(in), allocatable :: grid_root(:,:,:)
    real*8, intent(in), allocatable :: grid_cells(:,:,:,:,:)
    integer, intent(in) :: ntasks
    integer, intent(in) :: n_sites
    real*8, intent(in), allocatable :: positions(:)
    real*8, intent(in) :: a_box(3)
    real*8, intent(in) :: b_box(3)
    real*8, intent(in) :: c_box(3)

    integer :: i, j, k
    integer :: cell(3)
    integer, allocatable :: color
    integer :: rank

    allocate( placement(1:n_sites) )
    do i = 1, n_sites
      cell = find_grid_cell(positions(i), grid, surface, &
                            borders_a, borders_b, borders_c)
      placement(i) = grid_root(cell(0), cell(1), cell(2))
    end do
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine find_grid_cell(pos, grid, surface, &
                            borders_a, borders_b, borders_c)
    implicit none

    real*8, intent(in) :: pos(3)
    integer, intent(in) :: grid(3)
    real*8, intent(in) :: surface(3,3)
    real*8, intent(in) :: borders_a(:)
    real*8, intent(in) :: borders_b(:)
    real*8, intent(in) :: borders_c(:)

    integer, intent(out) :: cell(3)

    real*8 :: norm(3)
    integer :: i

    call projection(norm(1), surface(1), pos)
    call projection(norm(2), surface(2), pos)
    call projection(norm(3), surface(3), pos)

    do i = 1, grid(1)
      if (borders_a(i-1) < norm(1) .and. norm(1) <= borders_a(i)) then
        cell(1) = i
      end if
    end do
    do i = 1, grid(2)
      if (borders_b(i-1) < norm(2) .and. norm(2) <= borders_b(i)) then
        cell(2) = i
      end if
    end do
    do i = 1, grid(3)
      if (borders_c(i-1) < norm(3) .and. norm(3) <= borders_c(i)) then
        cell(3) = i
      end if
    end do

    return cell
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine dd_grid_cell(pos, grid, a_box, b_box, c_box)
    implicit none

    real*8, intent(in) :: pos(3)
    integer, intent(in) :: dd_grid(3)
    real*8, intent(in) :: a_box(3)
    real*8, intent(in) :: b_box(3)
    real*8, intent(in) :: c_box(3)

    real*8 :: cell_size(3)
    integer, intent(out) :: cell(3)

    cell_size(0) = a_box(0) / dd_grid(0)
    cell_size(1) = b_box(1) / dd_grid(1)
    cell_size(2) = c_box(2) / dd_grid(2)

    cell = (pos / cell_size) + 1

    return cell
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine surface_vector(s, x, y)
    implicit none

    real*8, intent(out) :: s(3)
    real*8, intent(in) :: x(3)
    real*8, intent(in) :: y(3)

    real*8 :: length

    s(1) = x(2) * y(3) - y(2) * x(3)
    s(2) = x(3) * y(1) - y(3) * x(1)
    s(3) = x(1) * y(2) - y(1) * x(2)

    length = sqrt(s(1)**2 + s(2)**2 + s(3)**2)

    s(:) = s(:) / length
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine dd_surface_vectors(surface, a_box, b_box, c_box)
    implicit none

    real*8, intent(out) :: surface(3,3)
    real*8, intent(in) :: a_box(3)
    real*8, intent(in) :: b_box(3)
    real*8, intent(in) :: c_box(3)

    call surface_vector(surface(1), b_box, c_box)
    call surface_vector(surface(2), c_box, a_box)
    call surface_vector(surface(3), a_box, b_box)
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine projection(norm, s, v)
    implicit none

    real*8, intent(out) :: norm
    real*8, intent(in) :: s(3)
    real*8, intent(in) :: v(3)

    norm = sqrt(s(1) * v(1) + s(2) * v(2) + s(3) * v(3))
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine dd_init_borders(dd_borders_a, dd_borders_b, dd_borders_c, &
                             grid, a_box, b_box, c_box, surface)
    implicit none

    real*8, intent(out), allocatable :: dd_borders_a(:)
    real*8, intent(out), allocatable :: dd_borders_b(:)
    real*8, intent(out), allocatable :: dd_borders_c(:)
    integer, intent(in) :: grid(3)
    real*8, intent(in) :: a_box(3)
    real*8, intent(in) :: b_box(3)
    real*8, intent(in) :: c_box(3)
    real*8, intent(in) :: surface(3,3)

    real*8 :: full(3), step(3)
    integer :: i

    allocate(dd_borders_a(grid(1) + 1))
    allocate(dd_borders_b(grid(2) + 1))
    allocate(dd_borders_c(grid(3) + 1))

    call projection(full(1), surface(1), a_box)
    call projection(full(2), surface(2), b_box)
    call projection(full(3), surface(3), c_box)

    step(:) = full(:) / grid(:)

    dd_borders_a(0) = 0.0
    dd_borders_b(0) = 0.0
    dd_borders_c(0) = 0.0
    dd_borders_a(grid(1) + 1) = full(i)
    dd_borders_b(grid(2) + 1) = full(i)
    dd_borders_c(grid(3) + 1) = full(i)

    do i = 2, grid(1)
      dd_borders_a(i) = step(1) * (i-1)
    end do
    do i = 2, grid(2)
      dd_borders_b(i) = step(2) * (i-1)
    end do
    do i = 2, grid(3)
      dd_borders_c(i) = step(3) * (i-1)
    end do
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine dd_init_cells(grid_cells, grid, a_box, b_box, c_box)
    implicit none

    real*8, intent(out), allocatable :: grid_cells(:,:,:,:,:)
    integer, intent(in) :: grid(3)
    real*8, intent(in) :: a_box(3)
    real*8, intent(in) :: b_box(3)
    real*8, intent(in) :: c_box(3)

    real*8 :: cell(3,3)
    integer :: i, j, k

    allocate(grid_cells(1:grid(1), 1:grid(2), 1:grid(3), 0:3, 1:3))

    cell(1) = a_box / grid(1)
    cell(2) = b_box / grid(2)
    cell(3) = c_box / grid(3)

    do i = 1, grid(1)
      do j = 1, grid(2)
        do k = 1, grid(3)
          grid_cells(i,j,k,0) = (i-1) * cell(1) + (j-1) * cell(2) + (k-1) * cell(3)
          grid_cells(i,j,k,1:3) = cell
        end do
      end do
    end do

    return cell
  end subroutine
!**************************************************************************

end module
