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

    real*8 :: norm(3,n_sites)
    integer :: i, j, n

    call vectorised_projection(norm, surface, positions, n_sites)

    do n = 1, n_sites
      do i = 1, 3
        do j = 1, grid(i)
          if (borders(i, j) < norm(i, n) &
              & .and. norm(i, n) <= borders(i, j + 1)) then
            cell(n, i) = j - 1
            exit
          end if
        end do
      end do
    end do
  end subroutine
!**************************************************************************



!**************************************************************************
  pure integer function get_neighbor_id(i, j, k) result(id)
    integer, intent(in) :: i
    integer, intent(in) :: j
    integer, intent(in) :: k

    integer :: id_grid(-1:1,-1:1,-1:1)
    ! self
    id_grid( 0,  0,  0) = 0
    ! orthogonal neighbors
    id_grid(-1,  0,  0) = 1
    id_grid( 1,  0,  0) = 2
    id_grid( 0, -1,  0) = 3
    id_grid( 0,  1,  0) = 4
    id_grid( 0,  0, -1) = 5
    id_grid( 0,  0,  1) = 6
    ! diagonal neighbors
    id_grid(-1, -1,  0) = 7
    id_grid(-1,  1,  0) = 8
    id_grid(-1,  0, -1) = 9
    id_grid(-1,  0,  1) = 10
    id_grid( 1, -1,  0) = 11
    id_grid( 1,  1,  0) = 12
    id_grid( 1,  0, -1) = 13
    id_grid( 1,  0,  1) = 14
    id_grid( 0, -1, -1) = 15
    id_grid( 0, -1,  1) = 16
    id_grid( 0,  1, -1) = 17
    id_grid( 0,  1,  1) = 18
    ! corners
    id_grid(-1, -1, -1) = 19
    id_grid(-1, -1,  1) = 20
    id_grid( 1, -1, -1) = 21
    id_grid( 1, -1,  1) = 22
    id_grid(-1,  1, -1) = 23
    id_grid(-1,  1,  1) = 24
    id_grid( 1,  1, -1) = 25
    id_grid( 1,  1,  1) = 26

    id = id_grid(i, j, k)
  end function
!**************************************************************************



!**************************************************************************
  subroutine init_grid_neighbors(grid_neighbor, grid_comm, local_comm, &
                                 local_rank)
    use mpi
    implicit none
    integer, intent(out) :: grid_neighbor(26)
    integer, intent(in) :: grid_comm
    integer, intent(in) :: local_comm
    integer, intent(in) :: local_rank

    integer :: i, j, k, ierr
    integer :: send_buffer(6)
    integer :: recv_buffer(6)
    integer :: status(MPI_STATUS_SIZE)

    if (local_rank == 0) then
       call mpi_cart_shift(grid_comm, 0, 1, grid_neighbor(1), grid_neighbor(2), ierr)
       call mpi_cart_shift(grid_comm, 1, 1, grid_neighbor(3), grid_neighbor(4), ierr)
       call mpi_cart_shift(grid_comm, 2, 1, grid_neighbor(5), grid_neighbor(6), ierr)

       send_buffer(1:4) = grid_neighbor(3:6)
       call mpi_sendrecv(send_buffer, 4, MPI_INTEGER, grid_neighbor(2), 0, &
                         recv_buffer, 4, MPI_INTEGER, grid_neighbor(1), 0, &
                         grid_comm, status, ierr)
       grid_neighbor(7:10) = recv_buffer(1:4)
       call mpi_sendrecv(send_buffer, 4, MPI_INTEGER, grid_neighbor(1), 0, &
                         recv_buffer, 4, MPI_INTEGER, grid_neighbor(2), 0, &
                         grid_comm, status, ierr)
       grid_neighbor(11:14) = recv_buffer(1:4)

       send_buffer(1:2) = grid_neighbor(5:6)
       send_buffer(3:4) = grid_neighbor(9:10)
       send_buffer(5:6) = grid_neighbor(13:14)
       call mpi_sendrecv(send_buffer, 6, MPI_INTEGER, grid_neighbor(4), 0, &
                         recv_buffer, 6, MPI_INTEGER, grid_neighbor(3), 0, &
                         grid_comm, status, ierr)
       grid_neighbor(15:16) = recv_buffer(1:2)
       grid_neighbor(19:22) = recv_buffer(3:6)
       call mpi_sendrecv(send_buffer, 4, MPI_INTEGER, grid_neighbor(3), 0, &
                         recv_buffer, 4, MPI_INTEGER, grid_neighbor(4), 0, &
                         grid_comm, status, ierr)
       grid_neighbor(17:18) = recv_buffer(1:2)
       grid_neighbor(23:26) = recv_buffer(3:6)
    end if
    call mpi_bcast(grid_neighbor, 26, MPI_INTEGER, 0, local_comm, ierr)

  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine domain_borders(border, grid_borders, grid_coords, global_rank)
    implicit none
    real*8, intent(out) :: border(6)
    real*8, intent(in) :: grid_borders(:,:)
    integer, intent(in) :: grid_coords(0:,:)
    integer, intent(in) :: global_rank

    integer :: i, j, k

    i = grid_coords(global_rank, 1) + 1
    j = grid_coords(global_rank, 2) + 1
    k = grid_coords(global_rank, 3) + 1
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
    logical, intent(out) :: mask(0:26, n_pos)
    real*8, intent(in) :: border(6)
    real*8, intent(in) :: norm(3, n_pos)
    integer, intent(in) :: n_pos

    logical :: mask_i(n_pos)
    logical :: mask_ij(n_pos)
    logical :: mask_b(3, -1:1, n_pos)
    integer :: i, j, k

    ! true if outside borders
    mask_b(1,-1, :) = norm(1,:) <= border(1)
    mask_b(1, 1, :) = norm(1,:) > border(2)
    mask_b(2,-1, :) = norm(2,:) <= border(3)
    mask_b(2, 1, :) = norm(2,:) > border(4)
    mask_b(3,-1, :) = norm(3,:) <= border(5)
    mask_b(3, 1, :) = norm(3,:) > border(6)
    ! true if inside borders
    mask_b(1, 0, :) = .not. (mask_b(1,-1,:) .or. mask_b(1,1,:))
    mask_b(2, 0, :) = .not. (mask_b(2,-1,:) .or. mask_b(2,1,:))
    mask_b(3, 0, :) = .not. (mask_b(3,-1,:) .or. mask_b(3,1,:))

    ! combine to form a mask for each neighbor (incl. diagonals)
    do i = -1, 1
       mask_i = mask_b(1,i,:)
       do j = -1, 1
          mask_ij = mask_i .and. mask_b(2,j,:)
          do k = -1, 1
             mask(get_neighbor_id(i,j,k),:) = mask_ij .and. mask_b(3,k,:)
          end do
       end do
    end do
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine check_migration_mask(mask, n_pos)
    implicit none
    logical, intent(out) :: mask(0:26, n_pos)
    integer, intent(in) :: n_pos

    integer :: i
    integer :: counts(27)

    do i = 0, 26
       counts(i+1) = count(mask(i,:))
    end do
    if (sum(counts) /= n_pos) then
       write(*,*) "Error: multiple migration targets"
       stop
    end if
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine migration_targets(targets, send_count, neighbors, mask, n_sites, &
                               global_rank)
    implicit none
    integer, intent(out) :: targets(n_sites)
    integer, intent(out) :: send_count(26)
    integer, intent(in) :: neighbors(26)
    logical, intent(in) :: mask(0:26, n_sites)
    integer, intent(in) :: n_sites
    integer, intent(in) :: global_rank

    integer :: i, n
    logical :: self(26)

    self(:) = (neighbors(:) == global_rank)

    targets = 0
    send_count = 0
    do i = 1, n_sites
       do n = 1, 26
          if (mask(n,i)) then
             if (.not. self(n)) then
                targets(i) = n
                send_count(n) = send_count(n) + 1
             end if
             exit
          end if
       end do
    end do
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine migrate(n_alloc, grid, grid_coords, surface, borders, neighbors, &
                     grid_comm, local_comm, global_rank, local_rank, &
                     n_sites, n_pos, n_sp, n_sp_sc, &
                     ids, positions, velocities, masses, xyz_species, &
                     species, xyz_species_supercell, species_supercell, &
                     fix_atom, positions_diff, positions_prev, forces_prev, &
                     debug)
    use mpi
    implicit none

    integer, intent(out) :: n_alloc
    integer, intent(in) :: grid(3)
    integer, intent(in) :: grid_coords(:,:)
    real*8, intent(in) :: surface(3,3)
    real*8, intent(in) :: borders(:,:)
    integer, intent(in) :: neighbors(26)
    integer, intent(in) :: grid_comm
    integer, intent(in) :: local_comm
    integer, intent(in) :: global_rank
    integer, intent(in) :: local_rank
    integer, intent(inout) :: n_sites
    integer, intent(inout) :: n_pos
    integer, intent(inout) :: n_sp
    integer, intent(inout) :: n_sp_sc
    integer, intent(inout), allocatable :: ids(:)
    real*8, intent(inout), allocatable :: positions(:,:)
    real*8, intent(inout), allocatable :: velocities(:,:)
    real*8, intent(inout), allocatable :: masses(:)
    character*8, intent(inout), allocatable :: xyz_species(:)
    integer, intent(inout), allocatable :: species(:)
    character*8, intent(inout), allocatable :: xyz_species_supercell(:)
    integer, intent(inout), allocatable :: species_supercell(:)
    logical, intent(inout), allocatable :: fix_atom(:,:)
    real*8, intent(inout), allocatable :: positions_diff(:,:)
    real*8, intent(inout), allocatable :: positions_prev(:,:)
    real*8, intent(inout), allocatable :: forces_prev(:,:)
    logical, intent(in) :: debug

    real*8 :: norm(3, n_sites)
    integer :: n, s, r, a, o, ierr
    integer :: n_recv, n_send, n_request
    real*8 :: local_border(6)
    logical :: mask(0:26, n_sites)
    integer :: send_count(26)
    integer :: recv_count(26)
    integer :: request(468)
    integer :: status(MPI_STATUS_SIZE, 468)
    !type(mpi_status) :: status(468)
    integer :: targets(n_sites)
    integer, allocatable :: sort_order(:)
    integer, allocatable :: buffer_ids(:)
    real*8, allocatable :: buffer_positions(:,:), buffer_positions_diff(:,:), &
                           buffer_positions_prev(:,:)
    real*8, allocatable :: buffer_velocities(:,:)
    real*8, allocatable :: buffer_masses(:)
    real*8, allocatable :: buffer_forces_prev(:,:)
    integer, allocatable :: buffer_species(:)
    integer, allocatable :: buffer_species_supercell(:)
    logical, allocatable :: buffer_fix_atom(:,:)
    character*8, allocatable :: buffer_xyz_species(:)
    character*8, allocatable :: buffer_xyz_species_supercell(:)


    ! FIXME: assuming n_sites == n_pos == n_sp == n_sp_sc

    call domain_borders(local_border, borders, grid_coords, global_rank)
    call vectorised_projection(norm, surface, positions, n_sites)
    call migration_mask(mask, local_border, norm, n_sites)
    if (debug) then
       call check_migration_mask(mask, n_sites)
    end if
    call migration_targets(targets, send_count, neighbors, mask, n_sites, &
                           global_rank)

    ! how many sites to migrate?
    n_request = 0
    recv_count(:) = 0
    do n = 1, 26
       if (neighbors(n) == global_rank) cycle ! skip self
       call mpi_irecv(recv_count(n), 1, MPI_INTEGER, neighbors(n), 0, &
                      grid_comm, request(n_request+1), ierr)
       call mpi_isend(send_count(n), 1, MPI_INTEGER, neighbors(n), 0, &
                      grid_comm, request(n_request+2), ierr)
       n_request = n_request + 2
    end do
    call mpi_waitall(n_request, request(1:n_request), status, ierr)
    ! sort arrays to move to-be-migrated sites to the end
    allocate(sort_order(n_sites))
    call get_sort_order(sort_order, targets, n_sites, 27)
    targets(1:n_sites) = targets(sort_order)
    ids(1:n_sites) = ids(sort_order)
    positions(1:3, 1:n_sites) = positions(1:3, sort_order)
    velocities(1:3, 1:n_sites) = velocities(1:3, sort_order)
    masses(1:n_sites) = masses(sort_order)
    xyz_species(1:n_sites) = xyz_species(sort_order)
    species(1:n_sites) = species(sort_order)
    xyz_species_supercell(1:n_sites) = xyz_species_supercell(sort_order)
    species_supercell(1:n_sites) = species_supercell(sort_order)
    fix_atom(1:3, 1:n_sites) = fix_atom(1:3, sort_order)
    positions_diff(1:3, 1:n_sites) = positions_diff(1:3, sort_order)
    positions_prev(1:3, 1:n_sites) = positions_prev(1:3, sort_order)
    forces_prev(1:3, 1:n_sites) = forces_prev(1:3, sort_order)
    ! resize receive arrays if needed
    n_send = sum(send_count)
    n_recv = sum(recv_count)
    n_alloc = n_sites - n_send + n_recv
    if (size(ids) < n_alloc) then
       n_alloc = 2 * size(ids)
       allocate(buffer_ids(n_alloc))
       buffer_ids(1:n_sites) = ids(1:n_sites)
       deallocate(ids)
       call move_alloc(buffer_ids, ids)
       allocate(buffer_positions(3, n_alloc))
       buffer_positions(1:3, 1:n_pos) = positions(1:3, 1:n_pos)
       deallocate(positions)
       call move_alloc(buffer_positions, positions)
       allocate(buffer_velocities(3, n_alloc))
       buffer_velocities(1:3, 1:n_pos) = velocities(1:3, 1:n_pos)
       deallocate(velocities)
       call move_alloc(buffer_velocities, velocities)
       allocate(buffer_masses(n_alloc))
       buffer_masses(1:n_sp) = masses(1:n_sp)
       deallocate(masses)
       call move_alloc(buffer_masses, masses)
       allocate(buffer_xyz_species(n_alloc))
       buffer_xyz_species(1:n_sp) = xyz_species(1:n_sp)
       deallocate(xyz_species)
       call move_alloc(buffer_xyz_species, xyz_species)
       allocate(buffer_species(n_alloc))
       buffer_species(1:n_sp) = species(1:n_sp)
       deallocate(species)
       call move_alloc(buffer_species, species)
       allocate(buffer_xyz_species_supercell(n_alloc))
       buffer_xyz_species_supercell(1:n_sp_sc) = xyz_species_supercell(1:n_sp_sc)
       deallocate(xyz_species_supercell)
       call move_alloc(buffer_xyz_species_supercell, xyz_species_supercell)
       allocate(buffer_species_supercell(n_alloc))
       buffer_species_supercell(1:n_sp_sc) = species_supercell(1:n_sp_sc)
       deallocate(species_supercell)
       call move_alloc(buffer_species_supercell, species_supercell)
       allocate(buffer_fix_atom(3, n_alloc))
       buffer_fix_atom(1:3,1:n_sp) = fix_atom(1:3,1:n_sp)
       deallocate(fix_atom)
       call move_alloc(buffer_fix_atom, fix_atom)
       allocate(buffer_positions_diff(3, n_alloc))
       buffer_positions_diff(1:3, 1:n_pos) = positions_diff(1:3, 1:n_pos)
       deallocate(positions_diff)
       call move_alloc(buffer_positions_diff, positions_diff)
       allocate(buffer_positions_prev(3, n_alloc))
       buffer_positions_prev(1:3, 1:n_pos) = positions_prev(1:3, 1:n_pos)
       deallocate(positions_prev)
       call move_alloc(buffer_positions_prev, positions_prev)
       allocate(buffer_forces_prev(3, n_alloc))
       buffer_forces_prev(1:3, 1:n_pos) = forces_prev(1:3, 1:n_pos)
       deallocate(forces_prev)
       call move_alloc(buffer_forces_prev, forces_prev)
    else
       n_alloc = 0
    end if
    ! allocate buffers
    allocate(buffer_ids(n_send))
    allocate(buffer_positions(3, n_send))
    allocate(buffer_velocities(3, n_send))
    allocate(buffer_masses(n_send))
    allocate(buffer_xyz_species(n_send))
    allocate(buffer_species(n_send))
    allocate(buffer_xyz_species_supercell(n_send))
    allocate(buffer_species_supercell(n_send))
    allocate(buffer_fix_atom(3, n_send))
    allocate(buffer_positions_diff(3, n_send))
    allocate(buffer_positions_prev(3, n_send))
    allocate(buffer_forces_prev(3, n_send))
    ! copy to-be-migrated sites to send buffers
    s = 1 + n_sites - n_send
    buffer_ids(1:n_send) = ids(s:n_sites)
    buffer_positions(1:3, 1:n_send) = positions(1:3, s:n_sites)
    buffer_velocities(1:3, 1:n_send) = velocities(1:3, s:n_sites)
    buffer_masses(1:n_send) = masses(s:n_sites)
    buffer_xyz_species(1:n_send) = xyz_species(s:n_sites)
    buffer_species(1:n_send) = species(s:n_sites)
    buffer_xyz_species_supercell(1:n_send) = xyz_species_supercell(s:n_sites)
    buffer_species_supercell(1:n_send) = species_supercell(s:n_sites)
    buffer_fix_atom(1:3, 1:n_send) = fix_atom(1:3, s:n_sites)
    buffer_positions_diff(1:3, 1:n_send) = positions_diff(1:3, s:n_sites)
    buffer_positions_prev(1:3, 1:n_send) = positions_prev(1:3, s:n_sites)
    buffer_forces_prev(1:3, 1:n_send) = forces_prev(1:3, s:n_sites)
    ! migrate sites
    n_request = 0
    s = 1
    r = 1 + n_sites - n_send
    do n = 1, 26
       if (neighbors(n) == global_rank) cycle ! skip self
       a = 12               ! no. of arrays to communicate
       o = n_request        ! request offset
       call mpi_irecv(ids(r:), recv_count(n), &
                      MPI_INTEGER, neighbors(n), 1, &
                      grid_comm, request(o+1), ierr)
       call mpi_irecv(positions(1:3,r:), 3 * recv_count(n), &
                      MPI_DOUBLE_PRECISION, neighbors(n), 2, &
                      grid_comm, request(o+2), ierr)
       call mpi_irecv(velocities(1:3,r:), 3 * recv_count(n), &
                      MPI_DOUBLE_PRECISION, neighbors(n), 3, &
                      grid_comm, request(o+3), ierr)
       call mpi_irecv(masses(r:), recv_count(n), &
                      MPI_DOUBLE_PRECISION, neighbors(n), 4, &
                      grid_comm, request(o+4), ierr)
       call mpi_irecv(xyz_species(r:), 8 * recv_count(n), &
                      MPI_CHARACTER, neighbors(n), 5, &
                      grid_comm, request(o+5), ierr)
       call mpi_irecv(species(r:), recv_count(n), &
                      MPI_INTEGER, neighbors(n), 6, &
                      grid_comm, request(o+6), ierr)
       call mpi_irecv(xyz_species_supercell(r:), 8 * recv_count(n), &
                      MPI_CHARACTER, neighbors(n), 7, &
                      grid_comm, request(o+7), ierr)
       call mpi_irecv(species_supercell(r:), recv_count(n), &
                      MPI_INTEGER, neighbors(n), 8, &
                      grid_comm, request(o+8), ierr)
       call mpi_irecv(fix_atom(1:3,r:), 3 * recv_count(n), &
                      MPI_LOGICAL, neighbors(n), 9, &
                      grid_comm, request(o+9), ierr)
       call mpi_irecv(positions_diff(1:3,r:), 3 * recv_count(n), &
                      MPI_DOUBLE_PRECISION, neighbors(n), 2, &
                      grid_comm, request(o+10), ierr)
       call mpi_irecv(positions_prev(1:3,r:), 3 * recv_count(n), &
                      MPI_DOUBLE_PRECISION, neighbors(n), 2, &
                      grid_comm, request(o+11), ierr)
       call mpi_irecv(forces_prev(1:3,r:), 3 * recv_count(n), &
                      MPI_DOUBLE_PRECISION, neighbors(n), 2, &
                      grid_comm, request(o+12), ierr)
       o = o + a
       call mpi_isend(buffer_ids(s:), send_count(n), &
                      MPI_INTEGER, neighbors(n), 1, &
                      grid_comm, request(o+1), ierr)
       call mpi_isend(buffer_positions(1:3,s:), 3 * send_count(n), &
                      MPI_DOUBLE_PRECISION, neighbors(n), 2, &
                      grid_comm, request(o+2), ierr)
       call mpi_isend(buffer_velocities(1:3,s:), 3 * send_count(n), &
                      MPI_DOUBLE_PRECISION, neighbors(n), 3, &
                      grid_comm, request(o+3), ierr)
       call mpi_isend(buffer_masses(s:), send_count(n), &
                      MPI_DOUBLE_PRECISION, neighbors(n), 4, &
                      grid_comm, request(o+4), ierr)
       call mpi_isend(buffer_xyz_species(s:), 8 * send_count(n), &
                      MPI_CHARACTER, neighbors(n), 5, &
                      grid_comm, request(o+5), ierr)
       call mpi_isend(buffer_species(s:), send_count(n), &
                      MPI_INTEGER, neighbors(n), 6, &
                      grid_comm, request(o+6), ierr)
       call mpi_isend(buffer_xyz_species_supercell(s:), 8 * send_count(n), &
                      MPI_CHARACTER, neighbors(n), 7, &
                      grid_comm, request(o+7), ierr)
       call mpi_isend(buffer_species_supercell(s:), send_count(n), &
                      MPI_INTEGER, neighbors(n), 8, &
                      grid_comm, request(o+8), ierr)
       call mpi_isend(buffer_fix_atom(1:3,s:), 3 * send_count(n), &
                      MPI_LOGICAL, neighbors(n), 9, &
                      grid_comm, request(o+9), ierr)
       call mpi_isend(buffer_positions_diff(1:3,s:), 3 * send_count(n), &
                      MPI_DOUBLE_PRECISION, neighbors(n), 2, &
                      grid_comm, request(o+10), ierr)
       call mpi_isend(buffer_positions_prev(1:3,s:), 3 * send_count(n), &
                      MPI_DOUBLE_PRECISION, neighbors(n), 2, &
                      grid_comm, request(o+11), ierr)
       call mpi_isend(buffer_forces_prev(1:3,s:), 3 * send_count(n), &
                      MPI_DOUBLE_PRECISION, neighbors(n), 2, &
                      grid_comm, request(o+12), ierr)
       n_request = o + a
       s = s + send_count(n)
       r = r + recv_count(n)
    end do
    call mpi_waitall(n_request, request(1:n_request), status, ierr)
    n_sites = n_sites - n_send + n_recv
    n_pos = n_pos - n_send + n_recv
    n_sp = n_sp - n_send + n_recv
    n_sp_sc = n_sp_sc - n_send + n_recv
    ! deallocate buffers
    deallocate(buffer_ids)
    deallocate(buffer_positions)
    deallocate(buffer_velocities)
    deallocate(buffer_masses)
    deallocate(buffer_xyz_species)
    deallocate(buffer_species)
    deallocate(buffer_xyz_species_supercell)
    deallocate(buffer_species_supercell)
    deallocate(buffer_fix_atom)
    deallocate(buffer_positions_diff)
    deallocate(buffer_positions_prev)
    deallocate(buffer_forces_prev)
  end subroutine
!**************************************************************************



!**************************************************************************
  subroutine exchange_mask(mask, border, norm, n, distance)
    implicit none
    logical, intent(out) :: mask(6, n)
    real*8, intent(in) :: border(6)
    real*8, intent(in) :: norm(3, n)
    integer, intent(in) :: n
    real*8, intent(in) :: distance

    mask(1,:) = norm(1,:) <= (border(1) + distance)
    mask(2,:) = norm(1,:) > (border(2) - distance)
    mask(3,:) = norm(2,:) <= (border(3) + distance)
    mask(4,:) = norm(2,:) > (border(4) - distance)
    mask(5,:) = norm(3,:) <= (border(5) + distance)
    mask(6,:) = norm(3,:) > (border(6) - distance)
  end subroutine
!**************************************************************************



!**************************************************************************
  pure logical function in_border_region(n, border, norm, distance) result(res)
    integer, intent(in) :: n
    real*8, intent(in) :: border(6)
    real*8, intent(in) :: norm(3)
    real*8, intent(in) :: distance

    integer :: axis

    axis = (n + 1) / 2
    if (mod(n,2) == 0) then
       res = norm(axis) > (border(n) - distance)
    else
       res = norm(axis) <= (border(n) + distance)
    end if
  end function
!**************************************************************************



!**************************************************************************
  subroutine halo_exchange(n_alloc, grid, grid_coords, surface, borders, neighbors, &
                           grid_comm, local_comm, global_rank, local_rank, &
                           rcut_max, n_sites, n_pos, n_sp, n_sp_sc, &
                           ids, positions, velocities, masses, xyz_species, &
                           species, xyz_species_supercell, species_supercell, &
                           fix_atom, debug)
    use mpi
    implicit none

    integer, intent(out) :: n_alloc
    integer, intent(in) :: grid(3)
    integer, intent(in) :: grid_coords(:,:)
    real*8, intent(in) :: surface(3,3)
    real*8, intent(in) :: borders(:,:)
    integer, intent(in) :: neighbors(26)
    integer, intent(in) :: grid_comm
    integer, intent(in) :: local_comm
    integer, intent(in) :: global_rank
    integer, intent(in) :: local_rank
    real*8, intent(in) :: rcut_max
    integer, intent(inout) :: n_sites
    integer, intent(inout) :: n_pos
    integer, intent(inout) :: n_sp
    integer, intent(inout) :: n_sp_sc
    integer, intent(inout), allocatable :: ids(:)
    real*8, intent(inout), allocatable :: positions(:,:)
    real*8, intent(inout), allocatable :: velocities(:,:)
    real*8, intent(inout), allocatable :: masses(:)
    character*8, intent(inout), allocatable :: xyz_species(:)
    integer, intent(inout), allocatable :: species(:)
    character*8, intent(inout), allocatable :: xyz_species_supercell(:)
    integer, intent(inout), allocatable :: species_supercell(:)
    logical, intent(inout), allocatable :: fix_atom(:,:)
    logical, intent(in) :: debug

    real*8, allocatable :: norm(:,:), tmp_norm(:,:)
    logical, allocatable :: mask(:,:), tmp_mask(:,:)
    integer :: n_alloc_old
    integer :: i, j, n, s, e, ierr
    integer :: src, tgt
    integer :: n_send, n_recv
    integer :: n_send_alloc=100
    real*8 :: local_border(6)
    integer :: status(MPI_STATUS_SIZE)
    integer, allocatable :: buffer_ids(:)
    real*8, allocatable :: buffer_positions(:,:)
    real*8, allocatable :: buffer_velocities(:,:)
    real*8, allocatable :: buffer_masses(:)
    integer, allocatable :: buffer_species(:)
    integer, allocatable :: buffer_species_supercell(:)
    logical, allocatable :: buffer_fix_atom(:,:)
    character*8, allocatable :: buffer_xyz_species(:)
    character*8, allocatable :: buffer_xyz_species_supercell(:)
    integer, allocatable :: tmp_ids(:)
    real*8, allocatable :: tmp_positions(:,:)
    real*8, allocatable :: tmp_velocities(:,:)
    real*8, allocatable :: tmp_masses(:)
    integer, allocatable :: tmp_species(:)
    integer, allocatable :: tmp_species_supercell(:)
    logical, allocatable :: tmp_fix_atom(:,:)
    character*8, allocatable :: tmp_xyz_species(:)
    character*8, allocatable :: tmp_xyz_species_supercell(:)

    call domain_borders(local_border, borders, grid_coords, global_rank)
    n_alloc = size(ids)
    allocate(norm(3, n_alloc))
    call vectorised_projection(norm(1:3,1:n_sites), surface, &
                               positions(1:3,1:n_sites), n_sites)
    allocate(mask(6, n_alloc))
    mask = .false.
    call exchange_mask(mask(1:6,1:n_sites), local_border, &
                       norm(1:3,1:n_sites), n_sites, rcut_max)

    ! allocate buffers
    allocate(buffer_ids(n_send_alloc))
    allocate(buffer_positions(3, n_send_alloc))
    allocate(buffer_velocities(3, n_send_alloc))
    allocate(buffer_masses(n_send_alloc))
    allocate(buffer_xyz_species(n_send_alloc))
    allocate(buffer_species(n_send_alloc))
    allocate(buffer_xyz_species_supercell(n_send_alloc))
    allocate(buffer_species_supercell(n_send_alloc))
    allocate(buffer_fix_atom(3, n_send_alloc))
    if (debug) then
       write(*,*) "send buffers allocated:", n_send_alloc
    end if

    do n = 1, 6
       if (mod(n,2) == 0) then
          ! shift up
          src = neighbors(n-1)
          tgt = neighbors(n)
       else
          ! shift down
          src = neighbors(n+1)
          tgt = neighbors(n)
       end if
       ! how many ghost sites to send / receive?
       n_send = count(mask(n,:))
       call mpi_sendrecv(n_send, 1, MPI_INTEGER, tgt, 0, &
                         n_recv, 1, MPI_INTEGER, src, 0, &
                         grid_comm, status, ierr)
       ! allocate buffers
       if (n_send_alloc < n_send) then
          deallocate(buffer_ids)
          deallocate(buffer_positions)
          deallocate(buffer_velocities)
          deallocate(buffer_masses)
          deallocate(buffer_xyz_species)
          deallocate(buffer_species)
          deallocate(buffer_xyz_species_supercell)
          deallocate(buffer_species_supercell)
          deallocate(buffer_fix_atom)
          n_send_alloc = n_send * 2
          allocate(buffer_ids(n_send_alloc))
          allocate(buffer_positions(3, n_send_alloc))
          allocate(buffer_velocities(3, n_send_alloc))
          allocate(buffer_masses(n_send_alloc))
          allocate(buffer_xyz_species(n_send_alloc))
          allocate(buffer_species(n_send_alloc))
          allocate(buffer_xyz_species_supercell(n_send_alloc))
          allocate(buffer_species_supercell(n_send_alloc))
          allocate(buffer_fix_atom(3, n_send_alloc))
          if (debug) then
             write(*,*) "send buffers allocated:", n_send_alloc
          end if
       end if
       ! copy ghost sites to send buffers
       j = 1
       do i = 1, n_sites
          if (mask(n,i)) then
             buffer_ids(j) = ids(i)
             buffer_positions(1:3, j) = positions(1:3, i)
             buffer_velocities(1:3, j) = velocities(1:3, i)
             buffer_masses(j) = masses(i)
             buffer_xyz_species(j) = xyz_species(i)
             buffer_species(j) = species(i)
             buffer_xyz_species_supercell(j) = xyz_species_supercell(i)
             buffer_species_supercell(j) = species_supercell(i)
             buffer_fix_atom(1:3, j) = fix_atom(1:3, i)
             j = j + 1
          end if
       end do
       ! reallocate arrays if needed
       if (n_alloc < n_sites + n_recv) then
          n_alloc_old = n_alloc
          do while (n_alloc < n_sites + n_recv)
             n_alloc = n_alloc * 2
          end do
          allocate(tmp_ids(n_alloc))
          tmp_ids(1:n_alloc_old) = ids(1:n_alloc_old)
          deallocate(ids)
          call move_alloc(tmp_ids, ids)
          allocate(tmp_positions(3, n_alloc))
          tmp_positions(1:3, 1:n_alloc_old) = positions(1:3, 1:n_alloc_old)
          deallocate(positions)
          call move_alloc(tmp_positions, positions)
          allocate(tmp_velocities(3, n_alloc))
          tmp_velocities(1:3, 1:n_alloc_old) = velocities(1:3, 1:n_alloc_old)
          deallocate(velocities)
          call move_alloc(tmp_velocities, velocities)
          allocate(tmp_masses(n_alloc))
          tmp_masses(1:n_alloc_old) = masses(1:n_alloc_old)
          deallocate(masses)
          call move_alloc(tmp_masses, masses)
          allocate(tmp_xyz_species(n_alloc))
          tmp_xyz_species(1:n_alloc_old) = xyz_species(1:n_alloc_old)
          deallocate(xyz_species)
          call move_alloc(tmp_xyz_species, xyz_species)
          allocate(tmp_species(n_alloc))
          tmp_species(1:n_alloc_old) = species(1:n_alloc_old)
          deallocate(species)
          call move_alloc(tmp_species, species)
          allocate(tmp_xyz_species_supercell(n_alloc))
          tmp_xyz_species_supercell(1:n_alloc_old) = xyz_species_supercell(1:n_alloc_old)
          deallocate(xyz_species_supercell)
          call move_alloc(tmp_xyz_species_supercell, xyz_species_supercell)
          allocate(tmp_species_supercell(n_alloc))
          tmp_species_supercell(1:n_alloc_old) = species_supercell(1:n_alloc_old)
          deallocate(species_supercell)
          call move_alloc(tmp_species_supercell, species_supercell)
          allocate(tmp_fix_atom(3, n_alloc))
          tmp_fix_atom(1:3,1:n_alloc_old) = fix_atom(1:3,1:n_alloc_old)
          deallocate(fix_atom)
          call move_alloc(tmp_fix_atom, fix_atom)
          allocate(tmp_norm(3, n_alloc))
          tmp_norm(1:3, 1:n_alloc_old) = norm(1:3, 1:n_alloc_old)
          deallocate(norm)
          call move_alloc(tmp_norm, norm)
          allocate(tmp_mask(6, n_alloc))
          tmp_mask(1:6, 1:n_alloc_old) = mask(1:6, 1:n_alloc_old)
          deallocate(mask)
          call move_alloc(tmp_mask, mask)
       end if
       ! halo exchange
       s = 1 + n_sites
       e = n_sites + n_recv
       call mpi_sendrecv(buffer_ids, n_send, &
                         MPI_INTEGER, tgt, 0, &
                         ids(s:e), n_recv, &
                         MPI_INTEGER, src, 0, &
                         grid_comm, status, ierr)
       call mpi_sendrecv(buffer_positions, 3 * n_send, &
                         MPI_DOUBLE_PRECISION, tgt, 0, &
                         positions(1:3, s:e), 3 * n_recv, &
                         MPI_DOUBLE_PRECISION, src, 0, &
                         grid_comm, status, ierr)
       call mpi_sendrecv(buffer_velocities, 3 * n_send, &
                         MPI_DOUBLE_PRECISION, tgt, 0, &
                         velocities(1:3, s:e), 3 * n_recv, &
                         MPI_DOUBLE_PRECISION, src, 0, &
                         grid_comm, status, ierr)
       call mpi_sendrecv(buffer_masses, n_send, &
                         MPI_DOUBLE_PRECISION, tgt, 0, &
                         masses(s:e), n_recv, &
                         MPI_DOUBLE_PRECISION, src, 0, &
                         grid_comm, status, ierr)
       call mpi_sendrecv(buffer_xyz_species, 8 * n_send, &
                         MPI_CHARACTER, tgt, 0, &
                         xyz_species(s:e), 8 * n_recv, &
                         MPI_CHARACTER, src, 0, &
                         grid_comm, status, ierr)
       call mpi_sendrecv(buffer_species, n_send, &
                         MPI_INTEGER, tgt, 0, &
                         species(s:e), n_recv, &
                         MPI_INTEGER, src, 0, &
                         grid_comm, status, ierr)
       call mpi_sendrecv(buffer_xyz_species_supercell, 8 * n_send, &
                         MPI_CHARACTER, tgt, 0, &
                         xyz_species_supercell(s:e), 8 * n_recv, &
                         MPI_CHARACTER, src, 0, &
                         grid_comm, status, ierr)
       call mpi_sendrecv(buffer_species_supercell, n_send, &
                         MPI_INTEGER, tgt, 0, &
                         species_supercell(s:e), n_recv, &
                         MPI_INTEGER, src, 0, &
                         grid_comm, status, ierr)
       call mpi_sendrecv(buffer_fix_atom, 3 * n_send, &
                         MPI_LOGICAL, tgt, 0, &
                         fix_atom(1:3, s:e), 3 * n_recv, &
                         MPI_LOGICAL, src, 0, &
                         grid_comm, status, ierr)
       n_sites = n_sites + n_recv
       ! calculate new norms and update masks
       call vectorised_projection(norm(1:3, s:e), surface, &
                                  positions(1:3, s:e), n_recv)
       call exchange_mask(mask(1:6, s:e), local_border, norm(1:3, s:e), &
                          n_recv, rcut_max)
    end do
    ! deallocate buffers
    deallocate(buffer_ids)
    deallocate(buffer_positions)
    deallocate(buffer_velocities)
    deallocate(buffer_masses)
    deallocate(buffer_xyz_species)
    deallocate(buffer_species)
    deallocate(buffer_xyz_species_supercell)
    deallocate(buffer_species_supercell)
    deallocate(buffer_fix_atom)
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

    real*8, intent(out) :: norm(3,n_sites)
    real*8, intent(in) :: s(3,3)
    real*8, intent(in) :: v(3,n_sites)
    integer, intent(in) :: n_sites

    integer :: i, j

    do i = 1, n_sites
      do j = 1, 3
        norm(j,i) = s(j,1) * v(1,i) + s(j,2) * v(2,i) + s(j,3) * v(3,i)
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
