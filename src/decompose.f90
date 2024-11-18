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
  subroutine grid_dimensions(grid, dims)
    implicit none
    integer, intent(in) :: grid(3)
    integer, intent(out) :: dims

    integer :: i

    dims = 0
    do i = 1, 3
       if (grid(i) > 1) then
          dims = dims + 1
       end if
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
    integer :: shift

    ! shift key values, so that first index is always 1
    shift = 1 - minval(keys)

    ! count number of keys in each bin
    slice(:) = 0
    do i = 1, n
      k = keys(i) + shift
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
      k = keys(i) + shift
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
  subroutine domain_borders(border, grid_borders, grid_coords, &
                          global_rank)
    implicit none
    real*8, intent(out) :: border(6)
    real*8, intent(in) :: grid_borders(:,:)
    integer, intent(in) :: grid_coords(3,:)
    integer, intent(in) :: global_rank

    integer :: i, j, k

    i, j, k = grid_coords(global_rank)
    border(1) = grid_borders(1, i)
    border(2) = grid_borders(1, i+1)
    border(3) = grid_borders(2, j)
    border(4) = grid_borders(2, j+1)
    border(5) = grid_borders(3, k)
    border(6) = grid_borders(3, k+1)
  end subroutine
!**************************************************************************



!**************************************************************************
subroutine migration_mask(mask, border, norm, n_pos)
    implicit none
    logical, intent(out) :: mask(6, n_pos)
    real*8, intent(in) :: border(6)
    real*8, intent(in) :: norm(n_pos, 3)
    integer, intent(in) :: n_pos

    mask(1) = norm(:,1) <= border(1)
    mask(2) = norm(:,1) > border(2)
    mask(3) = norm(:,2) <= border(3)
    mask(4) = norm(:,2) > border(4)
    mask(5) = norm(:,3) <= border(5)
    mask(6) = norm(:,3) > border(6)
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine migrate(grid, surface, borders, neighbors, &
                     grid_comm, local_comm, global_rank, local_rank, &
                     n_sites, n_pos, n_sp, n_sp_sc, &
                     positions, velocities, masses, xyz_species, &
                     species, xyz_species_supercell, species_supercell, &
                     fix_atom, ids)
    implicit none

    integer, intent(in) :: grid(3)
    real*8, intent(in) :: surface(3,3)
    real*8, intent(in) :: borders(:,:)
    integer, intent(in) :: neighbors(6)
    integer, intent(in) :: grid_comm
    integer, intent(in) :: local_comm
    integer, intent(in) :: global_rank
    integer, intent(in) :: local_rank
    integer, intent(in) :: n_sites
    integer, intent(in) :: n_pos
    integer, intent(in) :: n_sp
    integer, intent(in) :: n_sp_sc
    real*8, intent(in) :: positions(3, n_pos)
    real*8, intent(in) :: velocities(3, n_sp)
    real*8, intent(in) :: masses(n_sp)
    integer, intent(in) :: xyz_species(n_sp)
    integer, intent(in) :: species(n_sp)
    integer, intent(in) :: xyz_species_supercell(n_sp_sc)
    integer, intent(in) :: species_supercell(n_sp_sc)
    logical, intent(in) :: fix_atom(3, n_sp)
    integer, intent(in) :: ids(n_sites)

    real*8 :: norm(n_pos, 3)
    integer :: i, j, n
    real*8 :: local_border(6)
    logical :: mask(6, n_pos)
    integer :: send_count(6)
    integer :: recv_count(6)
    integer :: status(MPI_STATUS_SIZE)

    call domain_borders(local_border, grid_borders, grid_coords, global_rank)
    call vectorised_projection(norm, surface, positions, n_sites)
    call migration_mask(mask, local_border, norm, n_pos)

    send_count(:) = count(mask(:))

    ! note: for domain (i, j, k)
    !   neighbors(1) -> (i-1, j, k)
    !   neighbors(2) -> (i+1, j, k)
    !   neighbors(3) -> (i, j-1, k)
    !     ...
    !   neighbors(6) -> (i, j, k+1)

    ! how many sites to migrate?
    mpi_sendrecv(send_count(1), 1, MPI_INTEGER, neighbors(1), 0, &
                 recv_count(2), 1, MPI_INTEGER, neighbors(2), 0, &
                 grid_comm, status, ierr)
    mpi_sendrecv(send_count(2), 1, MPI_INTEGER, neighbors(2), 0, &
                 recv_count(1), 1, MPI_INTEGER, neighbors(1), 0, &
                 grid_comm, status, ierr)
    mpi_sendrecv(send_count(3), 1, MPI_INTEGER, neighbors(3), 0, &
                 recv_count(4), 1, MPI_INTEGER, neighbors(4), 0, &
                 grid_comm, status, ierr)
    mpi_sendrecv(send_count(4), 1, MPI_INTEGER, neighbors(4), 0, &
                 recv_count(3), 1, MPI_INTEGER, neighbors(3), 0, &
                 grid_comm, status, ierr)
    mpi_sendrecv(send_count(5), 1, MPI_INTEGER, neighbors(5), 0, &
                 recv_count(6), 1, MPI_INTEGER, neighbors(6), 0, &
                 grid_comm, status, ierr)
    mpi_sendrecv(send_count(6), 1, MPI_INTEGER, neighbors(6), 0, &
                 recv_count(5), 1, MPI_INTEGER, neighbors(5), 0, &
                 grid_comm, status, ierr)
    ! allocate buffers
    allocate(send_buffer_r(max(send_count)))
    allocate(send_buffer_i(max(send_count)))
    allocate(send_buffer_l(max(send_count)))
    ! migrate sites
    send_buffer_r(:send_count(1)) = positions(mask(1))
    mpi_sendrecv(send_buffer_r, send_count(1), MPI_DOUBLE_PRECISION, neighbors(1), 0, &
                 positions(...), recv_count(1), MPI_DOUBLE_PRECISION, neighbors(2), 0, &
                 grid_comm, status, ierr)
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine exchange(grid, surface, borders, neighbors, &
                     grid_comm, local_comm, global_rank, local_rank, &
                     n_sites, n_pos, n_sp, n_sp_sc, &
                     positions, velocities, masses, xyz_species, &
                     species, xyz_species_supercell, species_supercell, &
                     fix_atom, ids)
    implicit none

    integer, intent(in) :: grid(3)
    real*8, intent(in) :: surface(3,3)
    real*8, intent(in) :: borders(:,:)
    integer, intent(in) :: neighbors(6)
    integer, intent(in) :: grid_comm
    integer, intent(in) :: local_comm
    integer, intent(in) :: global_rank
    integer, intent(in) :: local_rank
    integer, intent(in) :: n_sites
    integer, intent(in) :: n_pos
    integer, intent(in) :: n_sp
    integer, intent(in) :: n_sp_sc
    real*8, intent(in) :: positions(3, n_pos)
    real*8, intent(in) :: velocities(3, n_sp)
    real*8, intent(in) :: masses(n_sp)
    integer, intent(in) :: xyz_species(n_sp)
    integer, intent(in) :: species(n_sp)
    integer, intent(in) :: xyz_species_supercell(n_sp_sc)
    integer, intent(in) :: species_supercell(n_sp_sc)
    logical, intent(in) :: fix_atom(3, n_sp)
    integer, intent(in) :: ids(n_sites)

    real*8 :: norm(n_pos, 3)
    integer :: i, j, n
    real*8 :: local_border(6)
    logical :: mask(6, n_pos)
    integer :: send_count(6)
    integer :: recv_count(6)
    integer :: status(MPI_STATUS_SIZE)

    call domain_borders(local_border, grid_borders, grid_coords, global_rank)
    call vectorised_projection(norm, surface, positions, n_sites)
    call migration_mask(mask, local_border, norm, n_pos) ! FIXME: border shift

    send_count(:) = count(mask(:))

    ! note: for domain (i, j, k)
    !   neighbors(1) -> (i-1, j, k)
    !   neighbors(2) -> (i+1, j, k)
    !   neighbors(3) -> (i, j-1, k)
    !     ...
    !   neighbors(6) -> (i, j, k+1)

    ! how many ghost sites to send / receive?
    mpi_sendrecv(send_count(1), 1, MPI_INTEGER, neighbors(1), 0, &
                 recv_count(2), 1, MPI_INTEGER, neighbors(2), 0, &
                 grid_comm, status, ierr)
    mpi_sendrecv(send_count(2), 1, MPI_INTEGER, neighbors(2), 0, &
                 recv_count(1), 1, MPI_INTEGER, neighbors(1), 0, &
                 grid_comm, status, ierr)
    mpi_sendrecv(send_count(3), 1, MPI_INTEGER, neighbors(3), 0, &
                 recv_count(4), 1, MPI_INTEGER, neighbors(4), 0, &
                 grid_comm, status, ierr)
    mpi_sendrecv(send_count(4), 1, MPI_INTEGER, neighbors(4), 0, &
                 recv_count(3), 1, MPI_INTEGER, neighbors(3), 0, &
                 grid_comm, status, ierr)
    mpi_sendrecv(send_count(5), 1, MPI_INTEGER, neighbors(5), 0, &
                 recv_count(6), 1, MPI_INTEGER, neighbors(6), 0, &
                 grid_comm, status, ierr)
    mpi_sendrecv(send_count(6), 1, MPI_INTEGER, neighbors(6), 0, &
                 recv_count(5), 1, MPI_INTEGER, neighbors(5), 0, &
                 grid_comm, status, ierr)
    ! allocate buffers
    allocate(send_buffer_r(max(send_count)))
    allocate(send_buffer_i(max(send_count)))
    allocate(send_buffer_l(max(send_count)))
    ! halo exchange
    send_buffer_r(:send_count(1)) = positions(mask(1))
    mpi_sendrecv(send_buffer_r, send_count(1), MPI_DOUBLE_PRECISION, neighbors(1), 0, &
                 positions(...), recv_count(1), MPI_DOUBLE_PRECISION, neighbors(2), 0, &
                 grid_comm, status, ierr)
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
        norm(i,j) = s(j,1) * v(1,i) + s(j,2) * v(2,i) + s(j,3) * v(3,i)
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

    norm = s(1) * v(1) + s(2) * v(2) + s(3) * v(3)
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine allocate_grid_borders(borders, borders_size, grid)
    implicit none

    real*8, intent(out), allocatable :: borders(:,:)
    integer, intent(out) :: borders_size
    integer, intent(in) :: grid(3)

    integer :: m

    m = maxval(grid)
    allocate(borders(3, m + 1))
    borders_size = (m + 1) * 3
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine init_grid_borders(borders, grid, a_box, b_box, c_box, surface)
    implicit none

    real*8, intent(out) :: borders(:,:)
    integer, intent(in) :: grid(3)
    real*8, intent(in) :: a_box(3)
    real*8, intent(in) :: b_box(3)
    real*8, intent(in) :: c_box(3)
    real*8, intent(in) :: surface(3,3)

    real*8 :: full(3), step(3)
    integer :: i, j

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
                        grid, grid_dims, grid_coords, grid_root)
    use mpi
    implicit none

    integer, intent(in) :: rank
    integer, intent(in) :: ntasks
    integer, intent(in) :: local_rank
    integer, intent(in) :: global_rank
    integer, intent(in) :: color(:)
    integer, intent(in) :: grid(3)
    integer, intent(in) :: grid_dims
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
      write (*, "(i0,a,i0,a,i0,a,i0)") &
        grid_dims, 'D grid: ', grid(1), ' x ', grid(2), ' x ', grid(3)
      write (*,*) ""
      write (*,*) "Rank  local_rank  global_rank  color  grid_coords  grid_root"
      write (*,*) "------------------------------------------------------------"
      do i = 1, ntasks
        coord = grid_coords(glob(i), :)
        write (*, "(x,i4,x,i11,x,i12,x,i6,a,i0,x,i0,x,i0,a,x,i10)") &
          & i - 1, loc(i), glob(i), color(i), &
          & '    (', coord(1), coord(2), coord(3), ')  ', &
          & grid_root(coord(1), coord(2), coord(3))
      end do
      write (*,*) ""
    endif
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine prepare_grid_distribute(n, placement, elements, grid_size, &
                                     counts, displs)
    implicit none

    integer, intent(in) :: n
    integer, intent(in) :: placement(:)
    integer, intent(in) :: elements
    integer, intent(in) :: grid_size
    integer, intent(out) :: counts(:)
    integer, intent(out) :: displs(:)

    integer :: i, rank
    integer :: j = 1

    counts(:) = 0
    displs(:) = 0

    rank = 0
    do i = 1, size(placement)
      if (placement(i) /= rank) then
        j = j + 1
        rank = placement(i)
      endif
      counts(j) = counts(j) + 1
    end do
    counts(:) = counts(:) * elements
    do i = 2, grid_size
      displs(i) = displs(i - 1) + counts(i - 1)
    end do
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine grid_distribute(sendbuf, recvbuf, counts, displs, datatype, &
                             grid_comm, rank)
    use mpi
    implicit none

    type(*), intent(in) :: sendbuf(..)
    type(*) :: recvbuf(..)
    integer, intent(in) :: counts(:)
    integer, intent(in) :: displs(:)
    integer, intent(in) :: datatype
    integer, intent(in) :: grid_comm
    integer, intent(in) :: rank

    integer :: n
    integer :: ierr

    n = counts(rank + 1)
    call mpi_scatterv(sendbuf, counts, displs, datatype, &
                      recvbuf, n, datatype, &
                      0, grid_comm, ierr)
  end subroutine
!**************************************************************************

end module
