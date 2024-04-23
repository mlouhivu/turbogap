! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   TurboGAP
! HND X
! HND X   TurboGAP is copyright (c) 2019-2023, Miguel A. Caro and others
! HND X
! HND X   TurboGAP is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   This file, exp_interface.f90, is copyright (c) 2019-2023,
! HND X   Miguel A. Caro and Tigany Zarrouk
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

module exp_interface
  use types
  use read_files
#ifdef _MPIF90
  use mpi
#endif
  use exp_utils
  use soap_turbo_functions
contains

  ! This module implements the interfaces for the gradient of experimental functions

  ! subroutine calculate_forces_from_pdf( pair_distribution_der, rjs, xyz, &
  !      & neighbors_list, n_neigh, neighbor_species, species, forces0 )


  ! end subroutine calculate_forces_from_pdf

  subroutine calculate_pair_distribution( params, x_pair_distribution&
       &, y_pair_distribution, y_pair_distribution_temp,&
       & pair_distribution_partial, pair_distribution_partial_temp, &
       & n_species, n_atoms_of_species, n_sites, a_box, b_box, c_box,&
       & indices, md_istep, i_beg, i_end, j_beg, j_end, ierr, rjs, xyz, &
       & neighbors_list, n_neigh, neighbor_species, species, rank,&
       & do_derivatives, pair_distribution_der, pair_distribution_partial_der,&
       & pair_distribution_partial_temp_der, forces_pair_distribution)
    implicit none
    type(input_parameters), intent(inout) :: params
    real*8, allocatable, intent(out) :: x_pair_distribution(:),&
         & y_pair_distribution(:), pair_distribution_partial(:,:),&
         & n_atoms_of_species(:), pair_distribution_partial_temp(:,:),&
         & y_pair_distribution_temp(:), pair_distribution_der(:,:),&
         & pair_distribution_partial_der(:,:,:), &
         & pair_distribution_partial_temp_der(:,:,:), forces_pair_distribution(:,:)
    real*8,  intent(in), allocatable :: rjs(:), xyz(:,:)
    integer, intent(in), allocatable :: neighbors_list(:), n_neigh(:)&
         &, neighbor_species(:), species(:)
    real*8,  intent(in) :: a_box(1:3), b_box(1:3), c_box(1:3)
    real*8 :: v_uc, f, virial(1:3,1:3)
    integer, intent(in) :: n_species, n_sites, i_beg, i_end, j_beg, j_end
    integer, intent(in) :: indices(1:3), md_istep, rank
    integer, intent(inout) :: ierr
    real, allocatable :: factors(:), forces_pair_distribution_temp(:,:)
    integer :: i, j, k, l, i2, n_dim_partial, n_dim_idx
    logical, intent(in) :: do_derivatives
    character*1024 :: filename

    ! Things that are allocated here:
    ! Always:
    !  > x_pair_distribution
    !  > y_pair_distribution
    ! if pair_distribution_partial == .true.
    !  > pair_distribution_partial( n_samples, n_spec * (n_spec + 1)/2 )
    !  if do_derivatives == .true.
    !    > pair_distribution_partial_der( n_samples, n_spec * (n_spec + 1)/2, j_beg : j_end )


    ! first allocate the necessary arrays for the
    ! calculation of the pair correlation function
    if (allocated( x_pair_distribution)) deallocate(x_pair_distribution)
    if (allocated( y_pair_distribution)) deallocate(y_pair_distribution)

    allocate( x_pair_distribution( 1: params%pair_distribution_n_samples) )
    allocate( y_pair_distribution( 1: params%pair_distribution_n_samples) )


    if (params%pair_distribution_partial)then
       n_dim_partial = n_species * ( n_species + 1 ) / 2
       allocate(factors( 1:n_dim_partial ))

       n_dim_idx = 1
       outer: do i = 1, n_species
          do j = 1, n_species
             if (i > j) cycle

             if (i /= j)then
                factors(n_dim_idx) = 2.d0
             else
                factors(n_dim_idx) = 1.d0
             end if

             n_dim_idx = n_dim_idx + 1
             if ( n_dim_idx > n_dim_partial )then
                exit outer
             end if

          end do
       end do outer


       if (.not. allocated(pair_distribution_partial))then   !deallocate(pair_distribution_partial)
          allocate( pair_distribution_partial(1:params%pair_distribution_n_samples,&
               & 1 : n_dim_partial) )
       end if

       pair_distribution_partial = 0.d0


       if (do_derivatives)then
          allocate( pair_distribution_partial_der(1:params%pair_distribution_n_samples,&
            & 1 : n_dim_partial, j_beg : j_end  ))
          pair_distribution_partial_der = 0.d0
       end if
    else
       if (do_derivatives)then
          allocate( pair_distribution_der(1:params%pair_distribution_n_samples,&
            &  j_beg : j_end  ))
          pair_distribution_der = 0.d0
       end if


    end if

    if(allocated(n_atoms_of_species)) deallocate(n_atoms_of_species)
    allocate(n_atoms_of_species(1:n_species))

    do j = 1, n_species
       n_atoms_of_species(j) = 0.d0
       do i2 = 1, n_sites
          if ( species(i2) == j)then
             n_atoms_of_species(j) = n_atoms_of_species(j) + 1.d0
          end if
       end do
    end do


#ifdef _MPIF90
    if (params%pair_distribution_partial)then
       allocate( pair_distribution_partial_temp(1:params%pair_distribution_n_samples, 1 : n_dim_partial) )

       pair_distribution_partial_temp = 0.0d0


    end if

    allocate( y_pair_distribution_temp( 1: params%pair_distribution_n_samples) )
    y_pair_distribution_temp = 0.d0

#endif
    v_uc = dot_product( cross_product(a_box,&
         & b_box), c_box ) / (&
         & dfloat(indices(1)*indices(2)&
         &*indices(3)) )


    !#####################################################################!
    !###---   Calculating the partial pair distribution functions   ---###!
    !#####################################################################!

    if ( params%pair_distribution_partial )then
       n_dim_idx = 1
       outer1: do j = 1, n_species
          do k = 1, n_species

             if (j > k) cycle ! We have already calculated the pair correlation function!

             ! Note that with the calculation of the derivatives here,
             ! this is without the -2 *delta_ik (r_j^alpha - r_i^alpha)
             ! factor, which allows for some freeing of memory

             call get_pair_distribution( n_sites, &
                  & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
                  & neighbor_species(j_beg:j_end), rjs(j_beg:j_end), xyz(1:3,j_beg:j_end), params&
                  &%r_range_min, params%r_range_max, params%pair_distribution_n_samples, x_pair_distribution,&
                  & pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx), params&
                  &%pair_distribution_rcut, .false.,&
                  & params%pair_distribution_partial, j, k,&
                  & params%pair_distribution_kde_sigma,&
                  & dfloat(n_sites)/v_uc,  do_derivatives,&
                  & pair_distribution_partial_der(1:params&
                  &%pair_distribution_n_samples, n_dim_idx, &
                  & j_beg:j_end) )

             n_dim_idx = n_dim_idx + 1

             if ( n_dim_idx > n_dim_partial )then
                exit outer1
             end if

          end do
       end do outer1

    else
       call get_pair_distribution( n_sites, &
            & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
            & neighbor_species(j_beg:j_end), rjs(j_beg:j_end), xyz(1:3,j_beg:j_end),&
            & params%r_range_min, params%r_range_max, params &
            &%pair_distribution_n_samples, x_pair_distribution,&
            & y_pair_distribution, params &
            &%pair_distribution_rcut, .false., .false., 1, 1,&
            & params%pair_distribution_kde_sigma, dfloat(n_sites)&
            &/v_uc, do_derivatives, pair_distribution_der(1:params&
            &%pair_distribution_n_samples,&
            & j_beg:j_end))
    end if


    ! --- MPI communication is here  ---

    if ( params%pair_distribution_partial )then
#ifdef _MPIF90
       call mpi_reduce(pair_distribution_partial,&
            & pair_distribution_partial_temp, params&
            &%pair_distribution_n_samples * n_dim_partial, MPI_DOUBLE_PRECISION, MPI_SUM,&
            & 0, MPI_COMM_WORLD, ierr)

       ! Now store the FULL pair distribution function which comes from these partial pair distribution functions
       ! Note, we have only so far divided by 4 pi r^2 dr
       ! Therefore, we must scale by the density

       pair_distribution_partial =  pair_distribution_partial_temp
       deallocate( pair_distribution_partial_temp )


       call mpi_bcast(pair_distribution_partial, params&
            &%pair_distribution_n_samples * n_dim_partial, MPI_DOUBLE_PRECISION, 0,&
            & MPI_COMM_WORLD, ierr)


       ! Now, we have the derivatives of the partial pair distribution
       ! function with respect to the atom pairs in that rank
       !
       ! We can keep them in the rank and calculate forces


       ! call mpi_reduce(pair_distribution_partial_der,&
       !      & pair_distribution_partial_der_temp, params&
       !      &%pair_distribution_n_samples * n_species *&
       !      & n_species * 3 * n_pairs_tot, MPI_DOUBLE_PRECISION, MPI_SUM,&
       !      & 0, MPI_COMM_WORLD, ierr)

       ! ! Now store the FULL pair distribution function which comes from these partial pair distribution functions
       ! ! Note, we have only so far divided by 4 pi r^2 dr
       ! ! Therefore, we must scale by the density

       ! pair_distribution_partial_der =  pair_distribution_partial_der_temp
       ! deallocate( pair_distribution_partial_der_temp )



#endif


       if ( params%exp_forces )then
          allocate(forces_pair_distribution(1:3,1:n_sites))
          forces_pair_distribution = 0.d0
       end if


       y_pair_distribution = 0.d0
       n_dim_idx = 1
       outer2: do j = 1, n_species
          do k = 1, n_species

             if (j > k) cycle

             pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx) =&
                  & pair_distribution_partial(1:params&
                  &%pair_distribution_n_samples, n_dim_idx) * v_uc &
                  &  /  n_atoms_of_species(j) /  n_atoms_of_species(k) !real(n_sites)


             if (do_derivatives) then
                pair_distribution_partial_der(1:params&
                     &%pair_distribution_n_samples, n_dim_idx, &
                     & j_beg:j_end) =  pair_distribution_partial_der(1:params&
                     & %pair_distribution_n_samples, n_dim_idx, &
                     & j_beg:j_end) * v_uc /  n_atoms_of_species(j) / &
                     & n_atoms_of_species(k) !real(n_sites)
             end if

             

             y_pair_distribution(1:params%pair_distribution_n_samples) = &
                  & y_pair_distribution(1:params%pair_distribution_n_samples)  +  &
                  &  factors(n_dim_idx) * (n_atoms_of_species(j) * n_atoms_of_species(k)) * &
                  & pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx) &
                  &  /  dfloat(n_sites) / dfloat(n_sites)


             if (do_derivatives .and.  params%exp_forces .and. allocated( params%exp_energy_scales ))then
                call get_pair_distribution_forces(  n_sites, params%exp_energy_scales(params%pdf_idx),&
                     & params%exp_data(params%pdf_idx)%y,&
                     & forces_pair_distribution, virial,&
                     & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
                     & neighbor_species(j_beg:j_end), rjs(j_beg:j_end), xyz(1:3,j_beg:j_end), params&
                     &%r_range_min, params%r_range_max, params&
                     &%pair_distribution_n_samples,&
                     & pair_distribution_partial(1:params&
                     &%pair_distribution_n_samples, n_dim_idx), params%pair_distribution_rcut&
                     &, j, k, pair_distribution_partial_der(1:params &
                     &%pair_distribution_n_samples, n_dim_idx,&
                     & j_beg:j_end), params%pair_distribution_partial&
                     & )
             end if


             n_dim_idx = n_dim_idx + 1

             if ( n_dim_idx > n_dim_partial )then
                exit outer2
             end if

          end do
       end do outer2


#ifdef _MPIF90

       ! No need to broadcast tbe pair distribution partials again as each process has done the calculation
       ! call mpi_bcast(pair_distribution_partial, params&
       !      &%pair_distribution_n_samples * n_dim_idx, MPI_DOUBLE_PRECISION, 0,&
       !      & MPI_COMM_WORLD, ierr)

       call mpi_bcast(y_pair_distribution, params&
            &%pair_distribution_n_samples, MPI_DOUBLE_PRECISION, 0,&
            & MPI_COMM_WORLD, ierr)


#endif

       !##################################################################!
       !###---   If not doing partial pair distribution functions   ---###!
       !##################################################################!


    else
#ifdef _MPIF90
       call mpi_reduce(y_pair_distribution,&
            & y_pair_distribution_temp, params&
            &%pair_distribution_n_samples,&
            & MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
            & MPI_COMM_WORLD, ierr)

       y_pair_distribution =  y_pair_distribution_temp
       deallocate( y_pair_distribution_temp )

       call mpi_bcast(y_pair_distribution, params&
            &%pair_distribution_n_samples, MPI_DOUBLE_PRECISION, 0,&
            & MPI_COMM_WORLD, ierr)

       if ( do_derivatives .and. params%exp_forces .and. allocated( params%exp_energy_scales ) )then
          allocate( forces_pair_distribution_temp(1:3, 1:n_sites) )
          forces_pair_distribution_temp = 0.d0

          call mpi_reduce(forces_pair_distribution,&
               & forces_pair_distribution_temp, 3 * n_sites,&
               & MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
               & MPI_COMM_WORLD, ierr)

          forces_pair_distribution = forces_pair_distribution_temp
          deallocate( forces_pair_distribution_temp )

          call mpi_bcast(forces_pair_distribution, 3*n_sites, MPI_DOUBLE_PRECISION, 0,&
               & MPI_COMM_WORLD, ierr)

       end if

#endif
    end if


    if (params%pair_distribution_partial .and. allocated(factors)) deallocate(factors)

    ! Multiplying by the 1/density factor (V/n_sites)
    ! Still need to check that the integral ( 4 * pi * density * ( int dr r^2 g(r) ) ~= N_sites )

    if (.not. params%pair_distribution_partial) then
       y_pair_distribution =y_pair_distribution * &
            & dot_product( cross_product(a_box, b_box),&
            & c_box ) / (dfloat(indices(1)*indices(2)&
            &*indices(3))) / dfloat(n_sites) / dfloat(n_sites)
    end if


    ! Write out the partial pair distribution functions
    if (rank == 0 .and. params%write_pair_distribution) then

       if (params%pair_distribution_partial)then
          n_dim_idx = 1
          outer3: do j = 1, n_species
             do k = 1, n_species

                if (j > k) cycle

                write(filename,'(A)')&
                     & 'pair_distribution_' // trim(params&
                     &%species_types(j)) // '_' // trim(params&
                     &%species_types(k)) //&
                     & "_prediction.dat"
                call write_exp_datan(x_pair_distribution(1:params%pair_distribution_n_samples),&
                     & pair_distribution_partial(1:params%pair_distribution_n_samples, n_dim_idx),&
                     & md_istep <= 0, filename, 'pair_distribution')

                n_dim_idx = n_dim_idx + 1
                if ( n_dim_idx > n_dim_partial )then
                   exit outer3
                end if

             end do
          end do outer3
       end if

       write(filename,'(A)')&
            & "pair_distribution_total.dat"
       call write_exp_datan(x_pair_distribution(1:params%pair_distribution_n_samples),&
            &y_pair_distribution(1:params%pair_distribution_n_samples),&
            & md_istep <= 0, filename, "pair_distribution")

    end if

  end subroutine calculate_pair_distribution


  subroutine finalize_pair_distribution(params, x_pair_distribution&
       &, y_pair_distribution, y_pair_distribution_temp,&
       & pair_distribution_partial, pair_distribution_partial_temp,&
       & do_derivatives, pair_distribution_der, pair_distribution_partial_der,&
       & pair_distribution_partial_temp_der, forces_pair_distribution, n_atoms_of_species, rank)
    implicit none
    type(input_parameters), intent(in) :: params
    integer, intent(in) :: rank
    real*8, allocatable, intent(inout) :: x_pair_distribution(:),&
         & y_pair_distribution(:), pair_distribution_partial(:,:),&
         & n_atoms_of_species(:), pair_distribution_partial_temp(:,:),&
         & y_pair_distribution_temp(:), pair_distribution_der(:,:),&
         & pair_distribution_partial_der(:,:,:), &
         & pair_distribution_partial_temp_der(:,:,:), forces_pair_distribution(:,:)
    logical, intent(in) :: do_derivatives

    ! Naive finalization, include the logic of how things are actually
    ! allocated above rather then allocating and deallocating

    if ( allocated( x_pair_distribution )               ) deallocate(x_pair_distribution)
    if ( allocated( y_pair_distribution )               ) deallocate(y_pair_distribution)
    if ( allocated( y_pair_distribution_temp )          ) deallocate(y_pair_distribution_temp)
    if ( allocated( pair_distribution_partial )         ) deallocate(pair_distribution_partial)
    if ( allocated( pair_distribution_partial_temp )    ) deallocate(pair_distribution_partial_temp)
    if ( allocated( pair_distribution_der )             ) deallocate(pair_distribution_der)
    if ( allocated( pair_distribution_partial_der )     ) deallocate(pair_distribution_partial_der)
    if ( allocated( pair_distribution_partial_temp_der )) deallocate(pair_distribution_partial_temp_der)
    if ( allocated( forces_pair_distribution )          ) deallocate(forces_pair_distribution)
    if ( allocated( n_atoms_of_species )          ) deallocate(n_atoms_of_species)


  end subroutine finalize_pair_distribution



  subroutine calculate_structure_factor( params, x_structure_factor, x_structure_factor_temp,&
       & y_structure_factor, y_structure_factor_temp,&
       & structure_factor_partial, structure_factor_partial_temp,&
       & x_pair_distribution, y_pair_distribution, &
       & pair_distribution_partial, n_species, n_atoms_of_species,&
       & n_sites, a_box, b_box, c_box, indices, md_istep, i_beg,&
       & i_end, j_beg, j_end, ierr, rjs, neighbors_list, n_neigh,&
       & neighbor_species, species, rank , q_beg, q_end, ntasks, sinc_factor_matrix)
    implicit none
    type(input_parameters), intent(in) :: params
    real*8, allocatable, intent(out) :: x_structure_factor(:), x_structure_factor_temp(:), &
         & y_structure_factor(:), structure_factor_partial(:,:),&
         &  structure_factor_partial_temp(:,:),&
         & y_structure_factor_temp(:)
    real*8,  intent(in), allocatable :: rjs(:),&
         & x_pair_distribution(:), y_pair_distribution(:),&
         & pair_distribution_partial(:,:), n_atoms_of_species(:)
    integer, intent(in), allocatable :: neighbors_list(:), n_neigh(:)&
         &, neighbor_species(:), species(:)
    real*8,  intent(in) :: a_box(1:3), b_box(1:3), c_box(1:3)
    real*8, allocatable, intent(inout) :: sinc_factor_matrix(:,:)
    real*8, allocatable :: sinc_factor_matrix_temp(:,:), temp_pdf(:,:)
    real*8 :: v_uc
    integer, intent(in) :: n_species, n_sites, i_beg, i_end, j_beg, j_end, ntasks
    integer, intent(out) :: q_beg, q_end
    integer, intent(in) :: indices(1:3), md_istep, rank
    integer, intent(out) :: ierr
    integer :: i, j, k, l, i2, n_dim_partial, n_dim_idx, n, m
    real*8 :: dq, f, cabh
    real*8, parameter :: pi = acos(-1.0)
    character*1024 :: filename


    if (params%structure_factor_from_rdf) then
       q_beg = 1
       q_end = params%structure_factor_n_samples
#ifdef _MPIF90
       ! We need to do an integral for each q value
       ! This can be split among the processes

       ! Split each q the integrals among each of the
       ! processes, just like we do with the atoms,
       ! and then collect accordingly


       if( rank < mod( params%structure_factor_n_samples, ntasks ) )then
          q_beg = 1 + rank*(params%structure_factor_n_samples / ntasks + 1)
       else
          q_beg = 1 + mod(params%structure_factor_n_samples, ntasks)*(params&
               &%structure_factor_n_samples / ntasks + 1) &
               &+ (rank - mod(params%structure_factor_n_samples, ntasks))*(params&
               &%structure_factor_n_samples / ntasks)
       end if
       if( rank < mod( params%structure_factor_n_samples, ntasks ) )then
          q_end = (rank+1)*(params%structure_factor_n_samples / ntasks + 1)
       else
          q_end = q_beg + params%structure_factor_n_samples/ntasks - 1
       end if


#endif
       if (params%structure_factor_from_rdf .and. params%pair_distribution_partial)then

          n_dim_partial = n_species * ( n_species + 1 ) / 2

          if (allocated(structure_factor_partial))  deallocate(structure_factor_partial)
          allocate( structure_factor_partial(1:params&
               &%structure_factor_n_samples, 1 : n_dim_partial) )
          structure_factor_partial = 0.d0
       end if

    end if

    if (allocated(y_structure_factor))deallocate(y_structure_factor)
    allocate(y_structure_factor(1:params%structure_factor_n_samples))
    y_structure_factor = 0.d0

    if (params%structure_factor_from_rdf .and. params%pair_distribution_partial) then
#ifdef _MPIF90
       allocate( structure_factor_partial_temp(1:params&
            &%structure_factor_n_samples, 1 : n_dim_partial) )

       structure_factor_partial_temp = 0.0d0
#endif

    else

#ifdef _MPIF90
       allocate( y_structure_factor_temp(1:params&
            &%structure_factor_n_samples) )

       y_structure_factor_temp = 0.0d0
#endif

    end if

    if (allocated(x_structure_factor)) deallocate(x_structure_factor)
    if (allocated(x_structure_factor_temp)) deallocate(x_structure_factor_temp)
    call linspace( x_structure_factor, params&
         &%q_range_min, params%q_range_max, params&
         &%structure_factor_n_samples, dq )

    call linspace( x_structure_factor_temp, params&
         &%q_range_min, params%q_range_max, params&
         &%structure_factor_n_samples, dq )

    if( trim(params%q_units) == "xrd" .or. params%q_units == "twotheta")then
       ! assume that theta is given for the Q range
       do i = 1, params%structure_factor_n_samples
          ! This gives s, but Q  = 2 pi * s = 4pi sin(theta) / lambda
          x_structure_factor(i) = 2.d0 * sin( pi *&
               & x_structure_factor(i) /180.d0 / 2.d0 ) /&
               & params%xrd_wavelength
       end do
    elseif (trim(params%q_units) == "saxs" .or. params%q_units == "q")then
       do i = 1, params%structure_factor_n_samples
          x_structure_factor(i) =  x_structure_factor(i)  / 2.d0 / pi
       end do
    end if

    ! Here the units of q range are that called Q (1/A) in
    ! the literature, small q in the literature are =
    ! 2sin(theta)/lambda = Q/2/pi

    if (params%structure_factor_from_rdf .and. params%pair_distribution_partial) then

       v_uc = dot_product( cross_product(a_box,&
            & b_box), c_box ) / (&
            & dfloat(indices(1)*indices(2)&
            &*indices(3)) )


       if ( .not. params%structure_factor_matrix ) then
          call get_partial_structure_factor(q_beg, q_end, &
               & pair_distribution_partial,&
               & x_structure_factor(1:params%structure_factor_n_samples) ,&
               & x_pair_distribution, params%pair_distribution_rcut,&
               & params %pair_distribution_n_samples, params&
               &%structure_factor_n_samples, n_species, n_dim_partial,&
               & n_atoms_of_species, n_sites, dfloat(n_sites)/v_uc, &
               & params%structure_factor_window,&
               & structure_factor_partial)



#ifdef _MPIF90
          call mpi_reduce(structure_factor_partial,&
               & structure_factor_partial_temp, params&
               &%structure_factor_n_samples * n_dim_partial, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
               & MPI_COMM_WORLD, ierr)
          structure_factor_partial = structure_factor_partial_temp
          deallocate(structure_factor_partial_temp)

          call mpi_bcast(structure_factor_partial, params&
               &%structure_factor_n_samples * n_dim_partial, MPI_DOUBLE_PRECISION, 0,&
               & MPI_COMM_WORLD, ierr)

#endif
       else
          if ( .not. allocated(sinc_factor_matrix) )then
             call get_sinc_factor_matrix( q_beg, q_end,&
                  & x_structure_factor, x_pair_distribution,&
                  & params%pair_distribution_rcut, params&
                  &%pair_distribution_n_samples, params&
                  &%structure_factor_n_samples, params%structure_factor_window,&
                  & sinc_factor_matrix )

#ifdef _MPIF90
             allocate( sinc_factor_matrix_temp(  1:params%structure_factor_n_samples, 1:params&
                  &%pair_distribution_n_samples ) )
             sinc_factor_matrix_temp = 0.d0

             call mpi_reduce(sinc_factor_matrix, sinc_factor_matrix_temp&
                  &, params%structure_factor_n_samples * params&
                  &%pair_distribution_n_samples ,&
                  & MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD,&
                  & ierr)

             sinc_factor_matrix = sinc_factor_matrix_temp

             deallocate(sinc_factor_matrix_temp)
             call mpi_bcast(sinc_factor_matrix, params&
                  &%structure_factor_n_samples * params&
                  &%pair_distribution_n_samples, MPI_DOUBLE_PRECISION, 0,&
                  & MPI_COMM_WORLD, ierr)

#endif
          end if

          ! Now make a blas call to perform the matrix multiplication necessary to obtain the structure factor.


          n = params%structure_factor_n_samples
          m = params%pair_distribution_n_samples
          k = n_dim_partial

          ! [sinc_factor_matrix] = n_samples_sf * n_samples_pc  ( N x M )
          ! [g_ab]               = n_samples_pc * n_dim_partial ( M x K )
          ! [S_ab]               = n_samples_sf * n_dim_partial ( N x K )
          ! S_ab = [sinc_factor_matrix] x [g_ab - 1]
          !      = ( N x M ) . ( M x K ) -> ( N x K )
          ! alpha * (A * B) + beta * C
          ! A === sinc_factor_partial
          ! B === [g_ab - 1]
          ! C === S_ab
          !   TRANS_A, TRANS_B  N_ROWS_A  N_COLS_B  N_COLS_A ( == N_ROWS_B ) alpha,
          !                                                   A,  first_dim_A,    B,      first_dim_B,  beta,  C, first_dim_C

          call dgemm("N", "N", n, k, m, 1.d0, sinc_factor_matrix, n,&
               & pair_distribution_partial - 1.d0, m, 0.d0,&
               & structure_factor_partial, n)


          ! All processes do this calculation, so they have the whole of the partial structure factors to work with.
          n_dim_idx = 1
          outersfm: do j = 1, n_species
             do k = 1, n_species

                if (j > k) cycle

                if (j == k) f = 1.d0
                if (j /= k) f = 0.d0

                cabh = ( (n_atoms_of_species(j) / dfloat( n_sites )) &
                     &*  (n_atoms_of_species(k) / dfloat( n_sites )) &
                     & )**(0.5)

                structure_factor_partial(1:params&
                     &%structure_factor_n_samples, n_dim_idx) = f +&
                     & 4.d0 * pi * cabh * (dfloat(n_sites)/v_uc) * &
                     & structure_factor_partial(1:params&
                     &%structure_factor_n_samples, n_dim_idx)

                n_dim_idx = n_dim_idx + 1
                if (n_dim_idx > n_dim_partial) exit outersfm
             end do
          end do outersfm

       end if


       ! using dq as a temp variable
       dq = 0.d0
       n_dim_idx = 1
       outer2: do j = 1, n_species
          dq = dq + (n_atoms_of_species(j) / dfloat(n_sites))
          do k = 1, n_species

             if (j > k) cycle

             if (j == k) f = 1.d0
             if (j /= k) f = 2.d0

             y_structure_factor(1:params%structure_factor_n_samples) = &
                  & y_structure_factor(1:params%structure_factor_n_samples)  +  &
                  & f * ( n_atoms_of_species(j) * n_atoms_of_species(k) )**0.5 * &
                  &  (structure_factor_partial(1:params%structure_factor_n_samples, n_dim_idx) &
                  &  )/ dfloat(n_sites) !/ dfloat(n_sites)

             n_dim_idx = n_dim_idx + 1

             if (n_dim_idx > n_dim_partial) exit outer2

          end do
       end do outer2
       y_structure_factor(1:params%structure_factor_n_samples)&
            & = y_structure_factor(1:params&
            &%structure_factor_n_samples) !/ dq


    elseif (params%structure_factor_from_rdf)then

       v_uc = dot_product( cross_product(a_box,&
            & b_box), c_box ) / (&
            & dfloat(indices(1)*indices(2)&
            &*indices(3)) )

       call get_structure_factor_from_rdf(q_beg, q_end, &
            & y_structure_factor(1:params%structure_factor_n_samples) ,&
            & y_pair_distribution,&
            & x_structure_factor(1:params%structure_factor_n_samples) ,&
            & x_pair_distribution, params%pair_distribution_rcut,&
            & params %pair_distribution_n_samples, params&
            &%structure_factor_n_samples, n_species,&
            & n_atoms_of_species, n_sites, dfloat(n_sites) / v_uc,  params%structure_factor_window)

    else

       call get_structure_factor_explicit( n_sites, &
            & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
            & neighbor_species(j_beg:j_end), rjs(j_beg:j_end),&
            & params%structure_factor_n_samples,&
            & x_structure_factor(1:params&
            &%structure_factor_n_samples), y_structure_factor(1:params&
            &%structure_factor_n_samples),params%r_range_min, params%r_range_max,&
            & params%pair_distribution_rcut, .false., 1&
            &, 1, params%structure_factor_window )
    end if


    if (.not. (params%structure_factor_from_rdf .and. params%pair_distribution_partial))then
#ifdef _MPIF90

       call mpi_reduce(y_structure_factor,&
            & y_structure_factor_temp, params&
            &%structure_factor_n_samples, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
            & MPI_COMM_WORLD, ierr)


       y_structure_factor = y_structure_factor_temp
       deallocate(y_structure_factor_temp)
#endif
       y_structure_factor =  y_structure_factor + 1.d0

#ifdef _MPIF90
       call mpi_bcast(y_structure_factor, params &
            &%structure_factor_n_samples,&
            & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif

    end if

    ! Write out the partial structure functions
    if (rank == 0 .and. params%write_structure_factor) then

       if (params%structure_factor_from_rdf .and. params%pair_distribution_partial)then
          n_dim_idx = 1
          outer3: do j = 1, n_species
             do k = 1, n_species

                if ( j > k ) cycle

                ! write with the temp data
                write(filename,'(A)')&
                     & 'structure_factor_' // trim(params&
                     &%species_types(j)) // '_' // trim(params&
                     &%species_types(k)) //&
                     & "_prediction.dat"
                call write_exp_datan(x_structure_factor_temp(1:params%structure_factor_n_samples)&
                     &,&
                     & structure_factor_partial(1:params&
                     &%structure_factor_n_samples, n_dim_idx),&
                     & md_istep <= 0, filename, "structure_factor: units of "// trim(params%q_units) )

                n_dim_idx = n_dim_idx + 1
                if ( n_dim_idx > n_dim_partial ) exit outer3

             end do
          end do outer3
       end if

       write(filename,'(A)')&
            & 'structure_factor_total.dat'
       call write_exp_datan(x_structure_factor_temp(1:params%structure_factor_n_samples),&
            & y_structure_factor(1:params&
            &%structure_factor_n_samples),&
            & md_istep <= 0, filename, "structure_factor: units of "// trim(params%q_units))


    end if
  end subroutine calculate_structure_factor


  subroutine finalize_structure_factor( params, x_structure_factor, x_structure_factor_temp,&
       & y_structure_factor, y_structure_factor_temp,&
       & structure_factor_partial, structure_factor_partial_temp,&
       & x_pair_distribution, y_pair_distribution, &
       & pair_distribution_partial, sinc_factor_matrix)
    implicit none
    type(input_parameters), intent(in) :: params
    real*8, allocatable, intent(inout) :: x_structure_factor(:), x_structure_factor_temp(:), &
         & y_structure_factor(:), structure_factor_partial(:,:),&
         &  structure_factor_partial_temp(:,:),&
         & y_structure_factor_temp(:)
    real*8,  intent(inout), allocatable :: x_pair_distribution(:), y_pair_distribution(:),&
         & pair_distribution_partial(:,:)
    real*8, allocatable, intent(inout) :: sinc_factor_matrix(:,:)

    if (allocated(x_structure_factor)) deallocate(x_structure_factor)
    if (allocated(x_structure_factor_temp)) deallocate(x_structure_factor_temp)
    if (allocated(y_structure_factor)) deallocate(y_structure_factor)
    if (allocated(structure_factor_partial)) deallocate(structure_factor_partial)
    if (allocated(structure_factor_partial_temp)) deallocate(structure_factor_partial_temp)
    if (allocated(y_structure_factor_temp)) deallocate(y_structure_factor_temp)
    if (allocated(x_pair_distribution)) deallocate(x_pair_distribution)
    if (allocated(y_pair_distribution)) deallocate(y_pair_distribution)
    if (allocated(pair_distribution_partial)) deallocate(pair_distribution_partial)
    if (allocated(sinc_factor_matrix)) deallocate(sinc_factor_matrix)

  end subroutine finalize_structure_factor




  subroutine calculate_xrd( params, x_xrd, x_xrd_temp,&
       & y_xrd, y_xrd_temp, x_structure_factor, x_structure_factor_temp,&
       & structure_factor_partial, structure_factor_partial_temp,&
       & n_species, n_atoms_of_species,&
       & n_sites, a_box, b_box, c_box, indices, md_istep, i_beg,&
       & i_end, j_beg, j_end, ierr, rjs, neighbors_list, n_neigh,&
       & neighbor_species, species, rank , q_beg, q_end, ntasks)
    implicit none
    type(input_parameters), intent(in) :: params
    real*8, allocatable, intent(out) :: x_xrd(:), x_xrd_temp(:), &
         & y_xrd(:), y_xrd_temp(:)
    real*8, allocatable, intent(in) :: structure_factor_partial(:,:),&
         &  structure_factor_partial_temp(:,:), x_structure_factor(:), x_structure_factor_temp(:)
    real*8,  intent(in), allocatable :: rjs(:), n_atoms_of_species(:)
    integer, intent(in), allocatable :: neighbors_list(:), n_neigh(:)&
         &, neighbor_species(:), species(:)
    real*8,  intent(in) :: a_box(1:3), b_box(1:3), c_box(1:3)
    real*8 :: v_uc
    integer, intent(in) :: n_species, n_sites, i_beg, i_end, j_beg, j_end, ntasks
    integer, intent(out) :: q_beg, q_end
    integer, intent(in) :: indices(1:3), md_istep, rank
    integer, intent(out) :: ierr
    integer :: i, j, k, l, i2, n_dim_idx, n_dim_partial
    real*8 :: dq
    real*8, parameter :: pi = acos(-1.0)
    character*1024 :: filename

    n_dim_partial = n_species * (n_species + 1 ) / 2

    ! Get the XRD from the partial structure factors!
    if (allocated(x_xrd)) deallocate(x_xrd)
    allocate(x_xrd(1:params%structure_factor_n_samples))
    if (allocated(x_xrd_temp)) deallocate(x_xrd_temp)
    allocate(x_xrd_temp(1:params%structure_factor_n_samples))

    x_xrd = x_structure_factor
    x_xrd_temp = x_structure_factor_temp


    if (allocated(y_xrd)) deallocate(y_xrd)
    allocate(y_xrd(1:params%structure_factor_n_samples))
    y_xrd = 0.d0

    ! allocate for the structure factor parameters, so we don't have to look through the horrible list



#ifdef _MPIF90
    allocate( y_xrd_temp(1:params&
         &%structure_factor_n_samples ))

    y_xrd_temp = 0.0d0
#endif

    ! Have the same range for the xrd as the structure factor

    if (params%structure_factor_from_rdf .and. params%pair_distribution_partial) then
       call get_xrd_from_partial_structure_factors(q_beg, q_end, &
            & structure_factor_partial(1:params&
            &%structure_factor_n_samples,1:n_dim_partial), n_species, params%species_types&
            &, species, params%xrd_wavelength, params &
            &%xrd_damping, params%xrd_alpha, params &
            &%xrd_method, params%xrd_iwasa, x_xrd(1:params&
            &%structure_factor_n_samples), y_xrd(1:params&
            &%structure_factor_n_samples),&
            & n_atoms_of_species )

    else

       call get_xrd_explicit( n_sites, params%species_types, n_species,  &
            & neighbors_list(j_beg:j_end), n_neigh(i_beg:i_end),&
            & neighbor_species(j_beg:j_end), rjs(j_beg:j_end),&
            & params%structure_factor_n_samples,&
            & x_xrd(1:params&
            &%structure_factor_n_samples), y_xrd(1:params&
            &%structure_factor_n_samples),params%r_range_min, params%r_range_max,&
            & params%pair_distribution_rcut, .false., 1&
            &, 1, params%structure_factor_window )


    end if

    !###################################################################################!
    !###---   Can calculate the Structure factors related to XRD / Neutron here   ---###!
    !###################################################################################!


    ! if (allocated(sf_parameters) ) deallocate( sf_parameters )
    ! allocate( sf_parameters(1:9,1:n_species))
    ! sf_parameters = 0.d0
    ! do i = 1, n_species
    !    call get_scattering_factor_params(params%species_types(i), sf_parameters(1:9,i))
    ! end do



    ! do j = 1, params%structure_factor_n_samples
    !    wfac = 0.d0
    !    do k = 1, n_species
    !       call get_scattering_factor(wfac_temp, sf_parameters(1:9,k), x_xrd(j)/2.d0 )
    !       wfac = wfac + wfac_temp*wfac_temp *n_atoms_of_species(k)  / dfloat(n_sites)
    !    end do
    !    print *, wfac
    !    y_xrd(j) = y_xrd(j) / (wfac)
    ! end do



#ifdef _MPIF90

    call mpi_reduce(y_xrd,&
         & y_xrd_temp, params&
         &%structure_factor_n_samples, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
         & MPI_COMM_WORLD, ierr)
    y_xrd = y_xrd_temp
    deallocate(y_xrd_temp)

    call mpi_bcast(y_xrd, params&
         &%structure_factor_n_samples, MPI_DOUBLE_PRECISION, 0,&
         & MPI_COMM_WORLD, ierr)

#endif

    if (rank == 0 .and. params%write_xrd) then

       write(filename,'(A)')&
            & 'xrd_prediction.dat'
       call write_exp_datan(x_xrd_temp(1:params%structure_factor_n_samples),&
            & y_xrd(1:params&
            &%structure_factor_n_samples),&
            & md_istep <= 0, filename, "xrd: units of "// trim(params%q_units))

    end if

  end subroutine calculate_xrd


  subroutine finalize_xrd( params, x_xrd, x_xrd_temp,&
       & y_xrd, y_xrd_temp, x_structure_factor, x_structure_factor_temp,&
       & structure_factor_partial, structure_factor_partial_temp)
    implicit none
    type(input_parameters), intent(in) :: params
    real*8, allocatable, intent(inout) :: x_xrd(:), x_xrd_temp(:), &
         & y_xrd(:), y_xrd_temp(:)
    real*8, allocatable, intent(inout) :: structure_factor_partial(:,:),&
         &  structure_factor_partial_temp(:,:), x_structure_factor(:), x_structure_factor_temp(:)

    if (allocated(x_xrd)) deallocate(x_xrd)
    if (allocated(x_xrd_temp)) deallocate(x_xrd_temp)
    if (allocated(y_xrd)) deallocate(y_xrd)
    if (allocated(y_xrd_temp)) deallocate(y_xrd_temp)
    if (allocated(x_structure_factor)) deallocate(x_structure_factor)
    if (allocated(x_structure_factor_temp)) deallocate(x_structure_factor_temp)
    if (allocated(structure_factor_partial)) deallocate(structure_factor_partial)
    if (allocated(structure_factor_partial_temp)) deallocate(structure_factor_partial_temp)
  end subroutine finalize_xrd




end module exp_interface
