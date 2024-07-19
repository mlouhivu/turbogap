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
  subroutine grid_affinity(color, grid, method, ntasks)
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
    integer, intent(out) :: grid_coords(0:,:)
    integer, intent(in) :: grid_comm
    integer, intent(in) :: ntasks

    integer :: i, ierr

    do i = 0, ntasks - 1
      call mpi_cart_coords(grid_comm, i, 3, grid_coords(i,:), ierr)
    end do
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine get_grid_root(grid_root, grid_coords, ntasks)
    implicit none
    integer, intent(out) :: grid_root(0:,0:,0:)
    integer, intent(in) :: grid_coords(0:,:)
    integer, intent(in) :: ntasks

    integer :: rank, i, j, k

    do rank = 0, ntasks - 1
      i = grid_coords(rank, 1)
      j = grid_coords(rank, 2)
      k = grid_coords(rank, 3)
      grid_root(i, j, k) = rank
    end do
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine grid_coord2id(id, grid, coord)
    implicit none
    integer, intent(out) :: id
    integer, intent(in) :: grid(3)
    integer, intent(in) :: coord(3)

    id = (coord(1) - 1) * grid(2) * grid(3) + &
       & (coord(2) - 1) * grid(3) + coord(3)
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine grid_id2coord(coord, grid, id)
    implicit none
    integer, intent(out) :: coord(3)
    integer, intent(in) :: grid(3)
    integer, intent(in) :: id

    integer tmp

    coord(1) = id / (grid(2) * grid(3)) + 1
    tmp = id - (coord(1) - 1) * grid(2) * grid(3)
    coord(2) = tmp / grid(3) + 1
    coord(3) = tmp - (coord(2) - 1) * grid(3) + 1
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine print_grid_coord(grid, id)
    implicit none
    integer, intent(in) :: grid(3)
    integer, intent(in) :: id
    integer :: coord(3)

    call grid_id2coord(coord, grid, id)
    write(*,"(a,i2,a,i2,x,i2,x,i2,a)") &
      & "Grid ID ", id, " = (", coord(1), coord(2), coord(3), ")"
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine get_sort_order(order, keys, n, bins)
    implicit none
    integer, intent(out) :: order(n)
    integer, intent(in) :: keys(n)
    integer, intent(in) :: n
    integer, intent(in) :: bins

    integer :: slice(bins)
    integer :: offset(bins)
    integer :: tally(bins)
    integer :: i, k

    ! count number of keys in each bin
    slice(:) = 0
    do i = 1, n
      k = keys(i) + 1
      if (k > bins) then
        write(*,*) "ERROR: key is larger than the number of bins"
        stop
      endif
      slice(k) = slice(k) + 1
    end do
    ! calculate position of the first element in each bin
    offset(:) = 1
    do i = 2, bins
      offset(i) = offset(i-1) + slice(i-1)
    end do

    ! sort indeces into bins based on the key values
    order(:) = 0
    tally(:) = 0
    do i = 1, n
      k = keys(i) + 1
      if (.not. order(offset(k) + tally(k)) == 0) then
        write(*,*) "WARN! already set.", k, i
      endif
      order(offset(k) + tally(k)) = i
      tally(k) = tally(k) + 1
    end do
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine grid_placement(placement, grid, grid_root, n_sites, positions, &
                            surface, borders)
    implicit none
    integer, intent(out), allocatable :: placement(:)
    integer, intent(in) :: grid(3)
    integer, intent(in) :: grid_root(0:,0:,0:)
    integer, intent(in) :: n_sites
    real*8, intent(in) :: positions(:,:)
    real*8, intent(in) :: surface(3,3)
    real*8, intent(in) :: borders(:,:)

    integer :: i
    integer :: cell(n_sites,3)

    allocate( placement(n_sites) )
    call find_grid_cell(cell, positions, grid, surface, borders, n_sites)
    do i = 1, n_sites
      placement(i) = grid_root(cell(i,1), cell(i,2), cell(i,3))
    end do
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine find_grid_cell(cell, positions, grid, surface, borders, n_sites)
    implicit none

    integer, intent(out) :: cell(n_sites,3)
    real*8, intent(in) :: positions(3,n_sites)
    integer, intent(in) :: grid(3)
    real*8, intent(in) :: surface(3,3)
    real*8, intent(in) :: borders(:,:)
    integer, intent(in) :: n_sites

    real*8 :: norm(n_sites,3)
    integer :: i, j, n

    call vectorised_projection(norm, surface, positions, n_sites)

    do n = 1, n_sites
      do i = 1, 3
        do j = 1, grid(i)
          if (borders(i, j) < norm(n, i) &
              & .and. norm(n, i) <= borders(i, j + 1)) then
            cell(n, i) = j - 1
            exit
          end if
        end do
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
  subroutine get_surface_vectors(surface, a_box, b_box, c_box)
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

    real*8, intent(out) :: norm(n_sites,3)
    real*8, intent(in) :: s(3,3)
    real*8, intent(in) :: v(3,n_sites)
    integer, intent(in) :: n_sites

    integer :: i, j

    do i = 1, n_sites
      do j = 1, 3
        norm(i,j) = sqrt(s(j,1) * v(1,i) + s(j,2) * v(2,i) + s(j,3) * v(3,i))
      end do
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
  subroutine init_grid_borders(borders, grid, a_box, b_box, c_box, surface)
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
      borders(i,1) = 0.0
      borders(i,grid(i) + 1) = full(i)
      do j = 2, grid(i)
        borders(i,j) = step(i) * (j-1)
      end do
    end do
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine print_grid(rank, ntasks, local_rank, global_rank, color, &
                        grid_coords, grid_root)
    use mpi
    implicit none

    integer, intent(in) :: rank
    integer, intent(in) :: ntasks
    integer, intent(in) :: local_rank
    integer, intent(in) :: global_rank
    integer, intent(in) :: color(:)
    integer, intent(in) :: grid_coords(0:,:)
    integer, intent(in) :: grid_root(0:,0:,0:)

    integer :: coord(3), i, ierr

    integer :: loc(ntasks), glob(ntasks)

    call mpi_gather(local_rank, 1, MPI_INTEGER, loc(:), 1, MPI_INTEGER, 0, &
                    & MPI_COMM_WORLD, ierr)
    call mpi_gather(global_rank, 1, MPI_INTEGER, glob(:), 1, MPI_INTEGER, 0, &
                    & MPI_COMM_WORLD, ierr)
    if (rank == 0) then
      write (*,*) ""
      write (*,*) "Rank  local_rank  global_rank  color  grid_coords  grid_root"
      write (*,*) "------------------------------------------------------------"
      do i = 1, ntasks
        coord = grid_coords(glob(i), :)
        write (*, "(x,i4,x,i11,x,i12,x,i6,a,i2,x,i2,x,i2,a,x,i10)") &
          & i - 1, loc(i), glob(i), color(i), &
          & '    (', coord(1), coord(2), coord(3), ')  ', &
          & grid_root(coord(1), coord(2), coord(3))
      end do
      write (*,*) ""
    endif
  end subroutine
!**************************************************************************

end module
