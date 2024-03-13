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
  subroutine check_grid(grid, method)
    use mpi
    implicit none

    integer, intent(in) :: grid(3)
    character*16, intent(in) :: method

    integer :: rank, ntasks, ierr
    integer :: grid_size

    call mpi_comm_size(MPI_COMM_WORLD, ntasks, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)

    grid_size = grid(1) * grid(2) * grid(3)

    if (grid_size > ntasks) then
      if (rank == 0) then
        write(*,'(a,i0,x,i0,x,i0,a,i0,a)') &
          & 'ERROR: domain decomposition grid (', &
          & grid(1), grid(2), grid(3), ') is too large for ', &
          & ntasks, ' MPI tasks.'
      end if
      call mpi_finalize(ierr)
      stop
    end if
    if (mod(ntasks, grid_size) > 0) then
      if (rank == 0) then
        write(*, '(a,i0,a,i0,x,i0,x,i0,a)') &
          & 'WARNING: non-uniform distribution of MPI tasks (', &
          & ntasks, ') in the domain decomposition grid (', &
          & grid(1), grid(2), grid(3), ')'
      end if
      if (method == "block") then
        if (rank == 0) then
          write(*, '(a,a)') &
            & 'ERROR: dd_grid_affinity "block" requires uniform', &
            & 'distribution of MPI tasks'
        end if
        call mpi_finalize(ierr)
        stop
      end if
    end if
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine dd_assign(color, grid, method, ntasks)
    implicit none

    integer, intent(out), allocatable :: color(:)
    integer, intent(in) :: grid(3)
    character*16, intent(in) :: method
    integer, intent(in) :: ntasks

    integer :: i, step, wrap, grid_size

    grid_size = grid(1) * grid(2) * grid(3)

    if (method == "cyclic") then
      ! first MPI ranks do MD
      step = 1
      wrap = grid_size
    else if (method == "block") then
      ! spread MD ranks uniformly
      step = (ntasks - 1) / grid_size + 1
      wrap = ntasks
    endif

    allocate(color(ntasks))
    do i = 1, ntasks
      color(i) = mod((i - 1) / step, wrap)
    end do
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine get_grid_coords(grid_coords, grid_comm, ntasks)
    implicit none
    integer, intent(out), allocatable :: grid_coords(:,:)
    integer, intent(in) :: grid_comm
    integer, intent(in) :: ntasks

    integer :: i, ierr
    integer :: coord(3)

    allocate(grid_coords(ntasks, 3))
    do i = 1, ntasks
      call mpi_cart_coords(grid_comm, i - 1, 3, grid_coords(i,:), ierr)
    end do
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine get_grid_root(grid_root, grid_coords, ntasks, grid)
    implicit none
    integer, intent(out), allocatable :: grid_root(:,:,:)
    integer, intent(in) :: grid_coords(:,:)
    integer, intent(in) :: ntasks
    integer, intent(in) :: grid(3)

    integer :: rank, i, j, k

    allocate(grid_root(grid(1), grid(2), grid(3)))
    do rank = 1, ntasks
      i = grid_coords(rank, 1)
      j = grid_coords(rank, 2)
      k = grid_coords(rank, 3)
      grid_root(i, j, k) = rank
    end do
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine dd_placement(placement, grid, grid_root, n_sites, positions, &
                          surface, borders)
    implicit none
    integer, intent(out), allocatable :: placement(:)
    integer, intent(in) :: grid(3)
    integer, intent(in) :: grid_root(:,:,:)
    integer, intent(in) :: n_sites
    real*8, intent(in) :: positions(:,:)
    real*8, intent(in) :: surface(3,3)
    real*8, intent(in) :: borders(:,:)

    integer :: i
    integer :: cell(3)

    allocate( placement(n_sites) )
    do i = 1, n_sites
      call find_grid_cell(cell, positions(i,:), grid, surface, borders)
      placement(i) = grid_root(cell(1), cell(2), cell(3))
    end do
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine find_grid_cell(cell, pos, grid, surface, borders)
    implicit none

    integer, intent(out) :: cell(3)
    real*8, intent(in) :: pos(3)
    integer, intent(in) :: grid(3)
    real*8, intent(in) :: surface(3,3)
    real*8, intent(in) :: borders(:,:)

    real*8 :: norm(3)
    integer :: i, j

    call projection(norm(1), surface(1,:), pos)
    call projection(norm(2), surface(2,:), pos)
    call projection(norm(3), surface(3,:), pos)

    do i = 1, 3
      do j = 1, grid(i)
        if (borders(i, j-1) < norm(i) .and. norm(i) <= borders(i, j)) then
          cell(i) = j
        end if
      end do
    end do
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

    call surface_vector(surface(1,:), b_box, c_box)
    call surface_vector(surface(2,:), c_box, a_box)
    call surface_vector(surface(3,:), a_box, b_box)
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine vectorised_projection(norm, s, v, n_sites)
    implicit none

    real*8, intent(out) :: norm(n_sites)
    real*8, intent(in) :: s(3)
    real*8, intent(in) :: v(n_sites,3)
    integer, intent(in) :: n_sites

    integer :: i

    do i = 1, n_sites
      norm(i) = sqrt(s(1) * v(i,1) + s(2) * v(i,2) + s(3) * v(i,3))
    end do
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
  subroutine dd_init_borders(borders, grid, a_box, b_box, c_box, surface)
    implicit none

    real*8, intent(out), allocatable :: borders(:,:)
    integer, intent(in) :: grid(3)
    real*8, intent(in) :: a_box(3)
    real*8, intent(in) :: b_box(3)
    real*8, intent(in) :: c_box(3)
    real*8, intent(in) :: surface(3,3)

    real*8 :: full(3), step(3)
    integer :: i, j, m

    m = maxval(grid)
    allocate(borders(3, m + 1))

    call projection(full(1), surface(1,:), a_box)
    call projection(full(2), surface(2,:), b_box)
    call projection(full(3), surface(3,:), c_box)

    step(:) = full(:) / grid(:)

    do i = 1,3
      borders(i,0) = 0.0
      borders(i,grid(i) + 1) = full(i)
      do j = 2, grid(i)
        borders(i,j) = step(i) * (j-1)
      end do
    end do
  end subroutine
!**************************************************************************

end module
