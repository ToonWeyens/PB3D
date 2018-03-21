!------------------------------------------------------------------------------!
!> Operations considering perturbation quantities.
!------------------------------------------------------------------------------!
module X_ops
#include <PB3D_macros.h>
#include <wrappers.h>
    use str_utilities
    use output_ops
    use messages
    use num_vars, only: dp, iu, max_str_ln, max_name_ln, pi
    use grid_vars, onlY: grid_type, disc_type
    use eq_vars, only: eq_1_type, eq_2_type
    use X_vars, only: X_1_type, X_2_type, modes_type

    implicit none
    private
    public calc_X, calc_magn_ints, print_output_X, resonance_plot, &
        &calc_res_surf, check_X_modes, init_modes, setup_modes, divide_X_jobs, &
        &redistribute_output_X
#if ldebug
    public debug_check_X_modes_2
#endif
    
    ! global variables
#if ldebug
    logical :: debug_check_X_modes_2 = .false.                                  !< plot debug information for check_x_modes_2() \ldebug
#endif
    
    ! interfaces
    
    !> \public Calculates either vectorial or tensorial perturbation variables.
    !!
    !! Optionally, the secondary mode number can  be specified (\c m if poloidal
    !! flux is used as normal coordinate and  c n if toroidal flux). By default,
    !! they are taken from the global \c X_vars variables.
    !!
    !! \return ierr
    interface calc_X
        !> \public
        module procedure calc_X_1
        !> \public
        module procedure calc_X_2
    end interface
    
    !> \public Redistribute the perturbation variables.
    !!
    !! \see redistribute_output_grid()
    !!
    !! \return ierr
    interface redistribute_output_X
        !> \public
        module procedure redistribute_output_X_1
        !> \public
        module procedure redistribute_output_X_2
    end interface
    
    !> \public Print either vectorial or  tensorial perturbation quantities of a
    !! certain order to an output file.
    !!
    !!  - vectorial:
    !!      - U
    !!      - DU
    !!  - tensorial:
    !!      - PV_int
    !!      - KV_int
    !!
    !! (The  non-integrated variables are not  saved because they are  heavy and
    !! not requested.)
    !!
    !! If \c  rich_lvl is provided, \c  "_R_[rich_lvl]" is appended to  the data
    !! name if it is \c >0.
    !!
    !! Optionally,  it  can be  specified that  this is  a divided
    !! parallel grid, corresponding  to the variable \c  eq_jobs_lims with index
    !! \c eq_job_nr. In  this case, the total  grid size is adjusted  to the one
    !! specified by \c eq_jobs_lims and the grid is written as a subset.
    !! \note 
    !!  -# The equilibrium quantities are outputted in Flux coordinates.
    !!  -# The  tensorial perturbation type  can also be used  for field-aligned
    !!  variables, in which case the first  index is assumed to have dimension 1
    !!  only. This can be triggered using \c is_field_averaged.
    !!
    !! \return ierr
    interface print_output_X
        !> \public
        module procedure print_output_X_1
        !> \public
        module procedure print_output_X_2
    end interface
    
contains
    !> \private vectorial version
    integer function calc_X_1(mds,grid_eq,grid_X,eq_1,eq_2,X,lim_sec_X) &
        &result(ierr)
        
        character(*), parameter :: rout_name = 'calc_X_1'
        
        ! input / output
        type(modes_type), intent(in) :: mds                                     !< general modes variables
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid variables
        type(grid_type), intent(in) :: grid_X                                   !< perturbation grid variables
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium
        type(eq_2_type), intent(in) :: eq_2                                     !< metric equilibrium
        type(X_1_type), intent(inout) :: X                                      !< vectorial perturbation variables
        integer, intent(in), optional :: lim_sec_X(2)                           !< limits of \c m_X (pol. flux) or \ n_X (tor. flux)
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Calculate vectorial perturbation variables')
        call lvl_ud(1)
        
        ! create perturbation with modes of current X job
        call X%init(mds,grid_X,lim_sec_X)
        
        ! calculate U and DU
        call writo('Calculating U and DU...')
        call lvl_ud(1)
        ierr = calc_U(grid_eq,grid_X,eq_1,eq_2,X)
        CHCKERR('')
        call lvl_ud(-1)
        
        ! user output
        call lvl_ud(-1)
    end function calc_X_1
    !> \private tensorial version
    integer function calc_X_2(mds,grid_eq,grid_X,eq_1,eq_2,X_a,X_b,X,&
        &lim_sec_X) result(ierr)
        
        character(*), parameter :: rout_name = 'calc_X_2'
        
        ! input / output
        type(modes_type), intent(in) :: mds                                     !< general modes variables
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid variables
        type(grid_type), intent(in) :: grid_X                                   !< perturbation grid variables
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium
        type(eq_2_type), intent(in) :: eq_2                                     !< metric equilibrium
        type(X_1_type), intent(inout) :: X_a                                    !< vectorial perturbation variables (dimension 1)
        type(X_1_type), intent(inout) :: X_b                                    !< vectorial perturbation variables (dimension 2)
        type(X_2_type), intent(inout) :: X                                      !< tensorial perturbation variables
        integer, intent(in), optional :: lim_sec_X(2,2)                         !< limits of \c m_X (pol flux) or \c n_X (tor flux) for both dimensions
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Calculate tensorial perturbation variables')
        call lvl_ud(1)
        
        ! create perturbation with modes of current X job
        call X%init(mds,grid_X,lim_sec_X)
        
        ! Calculate  PV_i  for  all (k,m)  pairs  and n_r  (equilibrium)
        ! values of the normal coordinate
        call writo('Calculating PV...')
        call lvl_ud(1)
        ierr = calc_PV(grid_eq,grid_X,eq_1,eq_2,X_a,X_b,X,lim_sec_X)
        CHCKERR('')
        call lvl_ud(-1)
        
        !  Calculate KV_i  for  all (k,m)  pairs  and n_r  (equilibrium)
        ! values of the normal coordinate
        call writo('Calculating KV...')
        call lvl_ud(1)
        ierr = calc_KV(grid_eq,grid_X,eq_1,eq_2,X_a,X_b,X,lim_sec_X)
        CHCKERR('')
        call lvl_ud(-1)
        
        ! user output
        call lvl_ud(-1)
    end function calc_X_2
    
    !> \private flux version
    integer function redistribute_output_X_1(mds,grid,grid_out,X,X_out) &
        &result(ierr)
        use MPI_utilities, only: redistribute_var
        
        character(*), parameter :: rout_name = 'redistribute_output_X_1'
        
        ! input / output
        type(modes_type), intent(in) :: mds                                     !< general modes variables
        type(grid_type), intent(in) :: grid                                     !< perturbation grid variables
        type(grid_type), intent(in) :: grid_out                                 !< redistributed perturbation grid variables
        type(X_1_type), intent(in) :: X                                         !< vectorial perturbation variables
        type(X_1_type), intent(inout) :: X_out                                  !< vectorial perturbation variables in redistributed grid
        
        ! local variables
        integer :: ld                                                           ! counter
        integer :: lims(2), lims_dis(2)                                         ! limits and distributed limits, taking into account the angular extent
        integer :: siz(3), siz_dis(3)                                           ! size for geometric part of variable
        real(dp), allocatable :: temp_var(:)                                    ! temporary variable
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Redistribute vectorial perturbation variables')
        call lvl_ud(1)
        
        ! test
        if (grid%n(1).ne.grid_out%n(1) .or. grid%n(2).ne.grid_out%n(2)) then
            ierr = 1
            CHCKERR('grid and grid_out are not compatible')
        end if
        
        ! create redistributed vectorial perturbation variables
        call X_out%init(mds,grid_out,lim_sec_X=X%lim_sec_X)
        
        ! set up limits taking into account angular extent and temporary var
        lims(1) = product(grid%n(1:2))*(grid%i_min-1)+1
        lims(2) = product(grid%n(1:2))*grid%i_max
        lims_dis(1) = product(grid%n(1:2))*(grid_out%i_min-1)+1
        lims_dis(2) = product(grid%n(1:2))*grid_out%i_max
        siz = [grid%n(1:2),grid%loc_n_r]
        siz_dis = [grid%n(1:2),grid_out%loc_n_r]
        allocate(temp_var(2*product(siz_dis)))                                  ! factor 2 because complex variable
        
        ! for all mode numbers
        do ld = 1,size(X%U_0,4)
            ! U_0
            ierr = redistribute_var(reshape(&
                &[rp(X%U_0(:,:,:,ld)),ip(X%U_0(:,:,:,ld))],[2*product(siz)]),&
                &temp_var,lims,lims_dis)
            CHCKERR('')
            X_out%U_0(:,:,:,ld) = &
                &reshape(temp_var(1:product(siz_dis)),siz_dis) + iu*&
                &reshape(temp_var(product(siz_dis)+1:2*product(siz_dis)),&
                &siz_dis)
            
            ! U_1
            ierr = redistribute_var(reshape(&
                &[rp(X%U_1(:,:,:,ld)),ip(X%U_1(:,:,:,ld))],[2*product(siz)]),&
                &temp_var,lims,lims_dis)
            CHCKERR('')
            X_out%U_1(:,:,:,ld) = &
                &reshape(temp_var(1:product(siz_dis)),siz_dis) + iu*&
                &reshape(temp_var(product(siz_dis)+1:2*product(siz_dis)),&
                &siz_dis)
            
            ! DU_0
            ierr = redistribute_var(reshape(&
                &[rp(X%DU_0(:,:,:,ld)),ip(X%DU_0(:,:,:,ld))],[2*product(siz)]),&
                &temp_var,lims,lims_dis)
            CHCKERR('')
            X_out%DU_0(:,:,:,ld) = &
                &reshape(temp_var(1:product(siz_dis)),siz_dis) + iu*&
                &reshape(temp_var(product(siz_dis)+1:2*product(siz_dis)),&
                &siz_dis)
            
            ! DU_1
            ierr = redistribute_var(reshape(&
                &[rp(X%DU_1(:,:,:,ld)),ip(X%DU_1(:,:,:,ld))],[2*product(siz)]),&
                &temp_var,lims,lims_dis)
            CHCKERR('')
            X_out%DU_1(:,:,:,ld) = &
                &reshape(temp_var(1:product(siz_dis)),siz_dis) + iu*&
                &reshape(temp_var(product(siz_dis)+1:2*product(siz_dis)),&
                &siz_dis)
        end do
        
        ! user output
        call lvl_ud(-1)
    end function redistribute_output_X_1
    !> \private metric version
    integer function redistribute_output_X_2(mds,grid,grid_out,X,X_out) &
        &result(ierr)
        use MPI_utilities, only: redistribute_var
        
        character(*), parameter :: rout_name = 'redistribute_output_X_2'
        
        ! input / output
        type(modes_type), intent(in) :: mds                                     !< general modes variables
        type(grid_type), intent(in) :: grid                                     !< perturbation grid variables
        type(grid_type), intent(in) :: grid_out                                 !< redistributed perturbation grid variables
        type(X_2_type), intent(in) :: X                                         !< tensorial perturbation variables
        type(X_2_type), intent(inout) :: X_out                                  !< tensorial perturbation variables in redistributed grid
        
        ! local variables
        integer :: ld                                                           ! counters
        integer :: n_loc(2)                                                     ! local n
        integer :: lims(2), lims_dis(2)                                         ! limits and distributed limits, taking into account the angular extent
        integer :: siz(3), siz_dis(3)                                           ! size for geometric part of variable
        logical :: is_field_averaged                                            ! whether field-averaged, so only one dimension for first index
        real(dp), allocatable :: temp_var(:)                                    ! temporary variable
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Redistribute tensorial perturbation variables')
        call lvl_ud(1)
        
        ! test
        is_field_averaged = size(X%PV_0,1).eq.1
        if ((grid%n(1).ne.grid_out%n(1) .or. grid%n(2).ne.grid_out%n(2)) .and. &
            &.not.is_field_averaged) then                                       ! for field-averaged grids, don't do tests
            ierr = 1
            CHCKERR('grid and grid_out are not compatible')
        end if
        
        ! create redistributed tensorial perturbation variables
        call X_out%init(mds,grid_out,lim_sec_X=X%lim_sec_X,&
            &is_field_averaged=is_field_averaged)
        
        ! set up limits taking into account angular extent and temporary var
        n_loc = grid%n(1:2)
        if (is_field_averaged) n_loc(1) = 1
        lims(1) = product(n_loc)*(grid%i_min-1)+1
        lims(2) = product(n_loc)*grid%i_max
        lims_dis(1) = product(n_loc)*(grid_out%i_min-1)+1
        lims_dis(2) = product(n_loc)*grid_out%i_max
        siz = [n_loc,grid%loc_n_r]
        siz_dis = [n_loc,grid_out%loc_n_r]
        allocate(temp_var(2*product(siz_dis)))                                  ! factor 2 because complex variable
        
        ! for all mode number combinations, symmetric
        do ld = 1,size(X%PV_0,4)
            ! PV_0
            ierr = redistribute_var(reshape(&
                &[rp(X%PV_0(:,:,:,ld)),ip(X%PV_0(:,:,:,ld))],[2*product(siz)]),&
                &temp_var,lims,lims_dis)
            CHCKERR('')
            X_out%PV_0(:,:,:,ld) = &
                &reshape(temp_var(1:product(siz_dis)),siz_dis) + iu*&
                &reshape(temp_var(product(siz_dis)+1:2*product(siz_dis)),&
                &siz_dis)
            
            ! PV_2
            ierr = redistribute_var(reshape(&
                &[rp(X%PV_2(:,:,:,ld)),ip(X%PV_2(:,:,:,ld))],[2*product(siz)]),&
                &temp_var,lims,lims_dis)
            CHCKERR('')
            X_out%PV_2(:,:,:,ld) = &
                &reshape(temp_var(1:product(siz_dis)),siz_dis) + iu*&
                &reshape(temp_var(product(siz_dis)+1:2*product(siz_dis)),&
                &siz_dis)
            
            ! KV_0
            ierr = redistribute_var(reshape(&
                &[rp(X%KV_0(:,:,:,ld)),ip(X%KV_0(:,:,:,ld))],[2*product(siz)]),&
                &temp_var,lims,lims_dis)
            CHCKERR('')
            X_out%KV_0(:,:,:,ld) = &
                &reshape(temp_var(1:product(siz_dis)),siz_dis) + iu*&
                &reshape(temp_var(product(siz_dis)+1:2*product(siz_dis)),&
                &siz_dis)
            
            ! KV_2
            ierr = redistribute_var(reshape(&
                &[rp(X%KV_2(:,:,:,ld)),ip(X%KV_2(:,:,:,ld))],[2*product(siz)]),&
                &temp_var,lims,lims_dis)
            CHCKERR('')
            X_out%KV_2(:,:,:,ld) = &
                &reshape(temp_var(1:product(siz_dis)),siz_dis) + iu*&
                &reshape(temp_var(product(siz_dis)+1:2*product(siz_dis)),&
                &siz_dis)
        end do
        
        ! for all mode number combinations, asymmetric
        do ld = 1,size(X%PV_1,4)
            ! PV_1
            ierr = redistribute_var(reshape(&
                &[rp(X%PV_1(:,:,:,ld)),ip(X%PV_1(:,:,:,ld))],[2*product(siz)]),&
                &temp_var,lims,lims_dis)
            CHCKERR('')
            X_out%PV_1(:,:,:,ld) = &
                &reshape(temp_var(1:product(siz_dis)),siz_dis) + iu*&
                &reshape(temp_var(product(siz_dis)+1:2*product(siz_dis)),&
                &siz_dis)
            
            ! KV_1
            ierr = redistribute_var(reshape(&
                &[rp(X%KV_1(:,:,:,ld)),ip(X%KV_1(:,:,:,ld))],[2*product(siz)]),&
                &temp_var,lims,lims_dis)
            CHCKERR('')
            X_out%KV_1(:,:,:,ld) = &
                &reshape(temp_var(1:product(siz_dis)),siz_dis) + iu*&
                &reshape(temp_var(product(siz_dis)+1:2*product(siz_dis)),&
                &siz_dis)
        end do
        
        ! user output
        call lvl_ud(-1)
    end function redistribute_output_X_2
    
    !> \private vectorial version.
    integer function print_output_X_1(grid_X,X,data_name,rich_lvl,par_div,&
        &lim_sec_X) result(ierr)
        
        use num_vars, only: PB3D_name, eq_jobs_lims, eq_job_nr
        use grid_utilities, only: trim_grid
        use HDF5_ops, only: print_HDF5_arrs
        use HDF5_vars, only: dealloc_var_1D, var_1D_type, &
            &max_dim_var_1D
        use X_vars, only: n_mod_X
        
        character(*), parameter :: rout_name = 'print_output_X_1'
        
        ! input / output
        type(grid_type), intent(in) :: grid_X                                   !< perturbation grid variables
        type(X_1_type), intent(in) :: X                                         !< vectorial perturbation variables 
        character(len=*), intent(in) :: data_name                               !< name under which to store
        integer, intent(in), optional :: rich_lvl                               !< Richardson level to print
        logical, intent(in), optional :: par_div                                !< is a parallely divided grid
        integer, intent(in), optional :: lim_sec_X(2)                           !< limits of \c m_X (pol. flux) or \c n_X (tor. flux)
        
        ! local variables
        type(grid_type) :: grid_trim                                            ! trimmed grid
        integer :: n_tot(3)                                                     ! total n
        integer :: par_id(2)                                                    ! parallel interval
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        type(var_1D_type), allocatable, target :: X_1D(:)                       ! 1D equivalent of X variables
        type(var_1D_type), pointer :: X_1D_loc => null()                        ! local element in X_1D
        integer :: id                                                           ! counters
        logical :: par_div_loc                                                  ! local par_div
        integer :: lim_sec_X_loc(2)                                             ! local lim_sec_X
        integer :: loc_size                                                     ! local size
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Write vectorial perturbation variables to output file')
        call lvl_ud(1)
        
        ! trim grid
        ierr = trim_grid(grid_X,grid_trim,norm_id)
        CHCKERR('')
        
        ! set local par_div
        par_div_loc = .false.
        if (present(par_div)) par_div_loc = par_div
        
        ! set total n and parallel interval
        n_tot = grid_trim%n
        par_id = [1,n_tot(1)]
        if (grid_trim%n(1).gt.0 .and. par_div_loc) then                         ! total grid includes all equilibrium jobs
            n_tot(1) = maxval(eq_jobs_lims)-minval(eq_jobs_lims)+1
            par_id = eq_jobs_lims(:,eq_job_nr)
        end if
        
        ! set local size and lim_sec_X
        loc_size = size(X%U_0(:,:,norm_id(1):norm_id(2),:))
        lim_sec_X_loc = [1,n_mod_X]
        if (present(lim_sec_X)) lim_sec_X_loc = lim_sec_X
        
        ! Set up the 1D equivalents of the perturbation variables
        allocate(X_1D(max_dim_var_1D))
        
        ! set up variables X_1D
        id = 1
        
        ! RE_U_0
        X_1D_loc => X_1D(id); id = id+1
        X_1D_loc%var_name = 'RE_U_0'
        allocate(X_1D_loc%tot_i_min(4),X_1D_loc%tot_i_max(4))
        allocate(X_1D_loc%loc_i_min(4),X_1D_loc%loc_i_max(4))
        X_1D_loc%tot_i_min = [1,1,1,1]
        X_1D_loc%tot_i_max = [n_tot,n_mod_X]
        X_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min,lim_sec_X_loc(1)]
        X_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max,&
            &lim_sec_X_loc(2)]
        allocate(X_1D_loc%p(loc_size))
        X_1D_loc%p = reshape(rp(X%U_0(:,:,norm_id(1):norm_id(2),:)),&
            &[loc_size])
        
        ! IM_U_0
        X_1D_loc => X_1D(id); id = id+1
        X_1D_loc%var_name = 'IM_U_0'
        allocate(X_1D_loc%tot_i_min(4),X_1D_loc%tot_i_max(4))
        allocate(X_1D_loc%loc_i_min(4),X_1D_loc%loc_i_max(4))
        X_1D_loc%tot_i_min = [1,1,1,1]
        X_1D_loc%tot_i_max = [n_tot,n_mod_X]
        X_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min,lim_sec_X_loc(1)]
        X_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max,&
            &lim_sec_X_loc(2)]
        allocate(X_1D_loc%p(loc_size))
        X_1D_loc%p = reshape(ip(X%U_0(:,:,norm_id(1):norm_id(2),:)),&
            &[loc_size])
        
        ! RE_U_1
        X_1D_loc => X_1D(id); id = id+1
        X_1D_loc%var_name = 'RE_U_1'
        allocate(X_1D_loc%tot_i_min(4),X_1D_loc%tot_i_max(4))
        allocate(X_1D_loc%loc_i_min(4),X_1D_loc%loc_i_max(4))
        X_1D_loc%tot_i_min = [1,1,1,1]
        X_1D_loc%tot_i_max = [n_tot,n_mod_X]
        X_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min,lim_sec_X_loc(1)]
        X_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max,&
            &lim_sec_X_loc(2)]
        allocate(X_1D_loc%p(loc_size))
        X_1D_loc%p = reshape(rp(X%U_1(:,:,norm_id(1):norm_id(2),:)),&
            &[loc_size])
        
        ! IM_U_1
        X_1D_loc => X_1D(id); id = id+1
        X_1D_loc%var_name = 'IM_U_1'
        allocate(X_1D_loc%tot_i_min(4),X_1D_loc%tot_i_max(4))
        allocate(X_1D_loc%loc_i_min(4),X_1D_loc%loc_i_max(4))
        X_1D_loc%tot_i_min = [1,1,1,1]
        X_1D_loc%tot_i_max = [n_tot,n_mod_X]
        X_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min,lim_sec_X_loc(1)]
        X_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max,&
            &lim_sec_X_loc(2)]
        allocate(X_1D_loc%p(loc_size))
        X_1D_loc%p = reshape(ip(X%U_1(:,:,norm_id(1):norm_id(2),:)),&
            &[loc_size])
        
        ! RE_DU_0
        X_1D_loc => X_1D(id); id = id+1
        X_1D_loc%var_name = 'RE_DU_0'
        allocate(X_1D_loc%tot_i_min(4),X_1D_loc%tot_i_max(4))
        allocate(X_1D_loc%loc_i_min(4),X_1D_loc%loc_i_max(4))
        X_1D_loc%tot_i_min = [1,1,1,1]
        X_1D_loc%tot_i_max = [n_tot,n_mod_X]
        X_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min,lim_sec_X_loc(1)]
        X_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max,&
            &lim_sec_X_loc(2)]
        allocate(X_1D_loc%p(loc_size))
        X_1D_loc%p = reshape(rp(X%DU_0(:,:,norm_id(1):norm_id(2),:)),&
            &[loc_size])
        
        ! IM_DU_0
        X_1D_loc => X_1D(id); id = id+1
        X_1D_loc%var_name = 'IM_DU_0'
        allocate(X_1D_loc%tot_i_min(4),X_1D_loc%tot_i_max(4))
        allocate(X_1D_loc%loc_i_min(4),X_1D_loc%loc_i_max(4))
        X_1D_loc%tot_i_min = [1,1,1,1]
        X_1D_loc%tot_i_max = [n_tot,n_mod_X]
        X_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min,lim_sec_X_loc(1)]
        X_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max,&
            &lim_sec_X_loc(2)]
        allocate(X_1D_loc%p(loc_size))
        X_1D_loc%p = reshape(ip(X%DU_0(:,:,norm_id(1):norm_id(2),:)),&
            &[loc_size])
        
        ! RE_DU_1
        X_1D_loc => X_1D(id); id = id+1
        X_1D_loc%var_name = 'RE_DU_1'
        allocate(X_1D_loc%tot_i_min(4),X_1D_loc%tot_i_max(4))
        allocate(X_1D_loc%loc_i_min(4),X_1D_loc%loc_i_max(4))
        X_1D_loc%tot_i_min = [1,1,1,1]
        X_1D_loc%tot_i_max = [n_tot,n_mod_X]
        X_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min,lim_sec_X_loc(1)]
        X_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max,&
            &lim_sec_X_loc(2)]
        allocate(X_1D_loc%p(loc_size))
        X_1D_loc%p = reshape(rp(X%DU_1(:,:,norm_id(1):norm_id(2),:)),&
            &[loc_size])
        
        ! IM_DU_1
        X_1D_loc => X_1D(id); id = id+1
        X_1D_loc%var_name = 'IM_DU_1'
        allocate(X_1D_loc%tot_i_min(4),X_1D_loc%tot_i_max(4))
        allocate(X_1D_loc%loc_i_min(4),X_1D_loc%loc_i_max(4))
        X_1D_loc%tot_i_min = [1,1,1,1]
        X_1D_loc%tot_i_max = [n_tot,n_mod_X]
        X_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min,lim_sec_X_loc(1)]
        X_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max,&
            &lim_sec_X_loc(2)]
        allocate(X_1D_loc%p(loc_size))
        X_1D_loc%p = reshape(ip(X%DU_1(:,:,norm_id(1):norm_id(2),:)),&
            &[loc_size])
        
        ! write
        ierr = print_HDF5_arrs(X_1D(1:id-1),PB3D_name,trim(data_name),&
            &rich_lvl=rich_lvl,ind_print=.not.grid_trim%divided)
        CHCKERR('')
        
        ! clean up
        call grid_trim%dealloc()
        call dealloc_var_1D(X_1D)
        nullify(X_1D_loc)
        
        ! user output
        call lvl_ud(-1)
    end function print_output_X_1
    !> \private tensorial version
    integer function print_output_X_2(grid_X,X,data_name,rich_lvl,par_div,&
        &lim_sec_X,is_field_averaged) result(ierr)
        use num_vars, only: PB3D_name, eq_jobs_lims, eq_job_nr
        use grid_utilities, only: trim_grid
        use HDF5_ops, only: print_HDF5_arrs
        use HDF5_vars, only: dealloc_var_1D, var_1D_type, &
            &max_dim_var_1D
        use X_utilities, only: is_necessary_X, get_sec_X_range
        use X_vars, only: set_nn_mod
        use num_utilities, only: c
        
        character(*), parameter :: rout_name = 'print_output_X_2'
        
        ! input / output
        type(grid_type), intent(in) :: grid_X                                   !< perturbation grid variables
        type(X_2_type), intent(in) :: X                                         !< tensorial perturbation variables 
        character(len=*), intent(in) :: data_name                               !< name under which to store
        integer, intent(in), optional :: rich_lvl                               !< Richardson level to print
        logical, intent(in), optional :: par_div                                !< is a parallely divided grid
        integer, intent(in), optional :: lim_sec_X(2,2)                         !< limits of \c m_X (pol. flux) or \c n_X (tor. flux)
        logical, intent(in), optional :: is_field_averaged                      !< whether field-averaged, so only one dimension for first index
        
        ! local variables
        type(grid_type) :: grid_trim                                            ! trimmed grid
        integer :: n_tot(3)                                                     ! total n
        integer :: par_id(2)                                                    ! parallel interval
        integer :: par_id_loc(2)                                                ! local parallel interval
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grids
        type(var_1D_type), allocatable, target :: X_1D(:)                       ! 1D equivalent of X variables
        type(var_1D_type), pointer :: X_1D_loc => null()                        ! local element in X_1D
        logical :: print_this(2)                                                ! whether symmetric and asymmetric variables need to be printed
        integer :: nn_mod_tot(2)                                                ! total nr. of modes for symmetric and asymmetric variables
        integer :: nn_mod_loc(2)                                                ! local nr. of modes for symmetric and asymmetric variables
        integer :: id                                                           ! counter
        logical :: par_div_loc                                                  ! local par_div
        integer :: loc_size(2)                                                  ! local size for symmetric and asymmetric variables
        integer :: m, k                                                         ! counters
        integer :: sXr_loc(2,2)                                                 ! local secondary X limits for symmetric and asymmetric variables
        integer :: sXr_tot(2,2)                                                 ! total secondary X limits for symmetric and asymmetric variables
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Write tensorial perturbation variables to output file')
        call lvl_ud(1)
        
        ! trim grid
        ierr = trim_grid(grid_X,grid_trim,norm_id)
        CHCKERR('')
        
        ! set local par_div
        par_div_loc = .false.
        if (present(par_div)) par_div_loc = par_div
        
        ! set total n, parallel interval and total n_mod_X
        n_tot = grid_trim%n
        par_id = [1,n_tot(1)]
        par_id_loc = [1,n_tot(1)]
        if (grid_trim%n(1).gt.0 .and. par_div_loc) then                         ! total grid includes all equilibrium jobs
            n_tot(1) = maxval(eq_jobs_lims)-minval(eq_jobs_lims)+1
            par_id = eq_jobs_lims(:,eq_job_nr)
            par_id_loc = par_id - par_id(1)+1
        else if (present(is_field_averaged)) then
            if (is_field_averaged) then                                         ! only first point
                par_id = [1,1]
                par_id_loc = [1,1]
                n_tot(1) = 1
            end if
        end if
        nn_mod_tot(1) = set_nn_mod(.true.)
        nn_mod_tot(2) = set_nn_mod(.false.)
        
        ! set dimensions, parallel limits and local and total nn_mod_X
        
        ! Set up the 1D equivalents of the perturbation variables
        allocate(X_1D(max_dim_var_1D))
        
        ! set up variables X_1D
        id = 1
        
        ! loop over modes of dimension 2
        do m = 1,X%n_mod(2)
            ! get contiguous range of modes of this m
            call get_sec_X_range(sXr_loc(:,1),sXr_tot(:,1),m,.true.,lim_sec_X)
            call get_sec_X_range(sXr_loc(:,2),sXr_tot(:,2),m,.false.,lim_sec_X)
            nn_mod_loc = sXr_loc(2,:)-sXr_loc(1,:)+1
            print_this = .false.
            do k = 1,2
                if (sXr_loc(1,k).le.sXr_loc(2,k)) print_this(k) = .true.        ! a bound is found
            end do
            
            ! set local size
            loc_size(1) = (par_id(2)-par_id(1)+1)*n_tot(2)*&
                &(norm_id(2)-norm_id(1)+1)*nn_mod_loc(1)
            loc_size(2) = (par_id(2)-par_id(1)+1)*n_tot(2)*&
                &(norm_id(2)-norm_id(1)+1)*nn_mod_loc(2)
            
            if (print_this(1)) then
                ! RE_PV_0
                X_1D_loc => X_1D(id); id = id+1
                X_1D_loc%var_name = 'RE_PV_0'
                allocate(X_1D_loc%tot_i_min(4),X_1D_loc%tot_i_max(4))
                allocate(X_1D_loc%loc_i_min(4),X_1D_loc%loc_i_max(4))
                X_1D_loc%tot_i_min = [1,1,1,1]
                X_1D_loc%tot_i_max = [n_tot,nn_mod_tot(1)]
                X_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min,sXr_tot(1,1)]
                X_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max,&
                    &sXr_tot(2,1)]
                allocate(X_1D_loc%p(loc_size(1)))
                X_1D_loc%p = reshape(rp(X%PV_0(par_id_loc(1):par_id_loc(2),:,&
                    &norm_id(1):norm_id(2),sXr_loc(1,1):sXr_loc(2,1))),&
                    &[loc_size(1)])
                
                ! IM_PV_0
                X_1D_loc => X_1D(id); id = id+1
                X_1D_loc%var_name = 'IM_PV_0'
                allocate(X_1D_loc%tot_i_min(4),X_1D_loc%tot_i_max(4))
                allocate(X_1D_loc%loc_i_min(4),X_1D_loc%loc_i_max(4))
                X_1D_loc%tot_i_min = [1,1,1,1]
                X_1D_loc%tot_i_max = [n_tot,nn_mod_tot(1)]
                X_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min,sXr_tot(1,1)]
                X_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max,&
                    &sXr_tot(2,1)]
                allocate(X_1D_loc%p(loc_size(1)))
                X_1D_loc%p = reshape(ip(X%PV_0(par_id_loc(1):par_id_loc(2),:,&
                    &norm_id(1):norm_id(2),sXr_loc(1,1):sXr_loc(2,1))),&
                    &[loc_size(1)])
                    
                ! RE_PV_2
                X_1D_loc => X_1D(id); id = id+1
                X_1D_loc%var_name = 'RE_PV_2'
                allocate(X_1D_loc%tot_i_min(4),X_1D_loc%tot_i_max(4))
                allocate(X_1D_loc%loc_i_min(4),X_1D_loc%loc_i_max(4))
                X_1D_loc%tot_i_min = [1,1,1,1]
                X_1D_loc%tot_i_max = [n_tot,nn_mod_tot(1)]
                X_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min,sXr_tot(1,1)]
                X_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max,&
                    &sXr_tot(2,1)]
                allocate(X_1D_loc%p(loc_size(1)))
                X_1D_loc%p = reshape(rp(X%PV_2(par_id_loc(1):par_id_loc(2),:,&
                    &norm_id(1):norm_id(2),sXr_loc(1,1):sXr_loc(2,1))),&
                    &[loc_size(1)])
                
                ! IM_PV_2
                X_1D_loc => X_1D(id); id = id+1
                X_1D_loc%var_name = 'IM_PV_2'
                allocate(X_1D_loc%tot_i_min(4),X_1D_loc%tot_i_max(4))
                allocate(X_1D_loc%loc_i_min(4),X_1D_loc%loc_i_max(4))
                X_1D_loc%tot_i_min = [1,1,1,1]
                X_1D_loc%tot_i_max = [n_tot,nn_mod_tot(1)]
                X_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min,sXr_tot(1,1)]
                X_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max,&
                    &sXr_tot(2,1)]
                allocate(X_1D_loc%p(loc_size(1)))
                X_1D_loc%p = reshape(ip(X%PV_2(par_id_loc(1):par_id_loc(2),:,&
                    &norm_id(1):norm_id(2),sXr_loc(1,1):sXr_loc(2,1))),&
                    &[loc_size(1)])
                
                ! RE_KV_0
                X_1D_loc => X_1D(id); id = id+1
                X_1D_loc%var_name = 'RE_KV_0'
                allocate(X_1D_loc%tot_i_min(4),X_1D_loc%tot_i_max(4))
                allocate(X_1D_loc%loc_i_min(4),X_1D_loc%loc_i_max(4))
                X_1D_loc%tot_i_min = [1,1,1,1]
                X_1D_loc%tot_i_max = [n_tot,nn_mod_tot(1)]
                X_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min,sXr_tot(1,1)]
                X_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max,&
                    &sXr_tot(2,1)]
                allocate(X_1D_loc%p(loc_size(1)))
                X_1D_loc%p = reshape(rp(X%KV_0(par_id_loc(1):par_id_loc(2),:,&
                    &norm_id(1):norm_id(2),sXr_loc(1,1):sXr_loc(2,1))),&
                    &[loc_size(1)])
                
                ! IM_KV_0
                X_1D_loc => X_1D(id); id = id+1
                X_1D_loc%var_name = 'IM_KV_0'
                allocate(X_1D_loc%tot_i_min(4),X_1D_loc%tot_i_max(4))
                allocate(X_1D_loc%loc_i_min(4),X_1D_loc%loc_i_max(4))
                X_1D_loc%tot_i_min = [1,1,1,1]
                X_1D_loc%tot_i_max = [n_tot,nn_mod_tot(1)]
                X_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min,sXr_tot(1,1)]
                X_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max,&
                    &sXr_tot(2,1)]
                allocate(X_1D_loc%p(loc_size(1)))
                X_1D_loc%p = reshape(ip(X%KV_0(par_id_loc(1):par_id_loc(2),:,&
                    &norm_id(1):norm_id(2),sXr_loc(1,1):sXr_loc(2,1))),&
                    &[loc_size(1)])
                    
                ! RE_KV_2
                X_1D_loc => X_1D(id); id = id+1
                X_1D_loc%var_name = 'RE_KV_2'
                allocate(X_1D_loc%tot_i_min(4),X_1D_loc%tot_i_max(4))
                allocate(X_1D_loc%loc_i_min(4),X_1D_loc%loc_i_max(4))
                X_1D_loc%tot_i_min = [1,1,1,1]
                X_1D_loc%tot_i_max = [n_tot,nn_mod_tot(1)]
                X_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min,sXr_tot(1,1)]
                X_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max,&
                    &sXr_tot(2,1)]
                allocate(X_1D_loc%p(loc_size(1)))
                X_1D_loc%p = reshape(rp(X%KV_2(par_id_loc(1):par_id_loc(2),:,&
                    &norm_id(1):norm_id(2),sXr_loc(1,1):sXr_loc(2,1))),&
                    &[loc_size(1)])
                
                ! IM_KV_2
                X_1D_loc => X_1D(id); id = id+1
                X_1D_loc%var_name = 'IM_KV_2'
                allocate(X_1D_loc%tot_i_min(4),X_1D_loc%tot_i_max(4))
                allocate(X_1D_loc%loc_i_min(4),X_1D_loc%loc_i_max(4))
                X_1D_loc%tot_i_min = [1,1,1,1]
                X_1D_loc%tot_i_max = [n_tot,nn_mod_tot(1)]
                X_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min,sXr_tot(1,1)]
                X_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max,&
                    &sXr_tot(2,1)]
                allocate(X_1D_loc%p(loc_size(1)))
                X_1D_loc%p = reshape(ip(X%KV_2(par_id_loc(1):par_id_loc(2),:,&
                    &norm_id(1):norm_id(2),sXr_loc(1,1):sXr_loc(2,1))),&
                    &[loc_size(1)])
            end if
            
            if (print_this(2)) then
                ! RE_PV_1
                X_1D_loc => X_1D(id); id = id+1
                X_1D_loc%var_name = 'RE_PV_1'
                allocate(X_1D_loc%tot_i_min(4),X_1D_loc%tot_i_max(4))
                allocate(X_1D_loc%loc_i_min(4),X_1D_loc%loc_i_max(4))
                X_1D_loc%tot_i_min = [1,1,1,1]
                X_1D_loc%tot_i_max = [n_tot,nn_mod_tot(2)]
                X_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min,sXr_tot(1,2)]
                X_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max,&
                    &sXr_tot(2,2)]
                allocate(X_1D_loc%p(loc_size(2)))
                X_1D_loc%p = reshape(rp(X%PV_1(par_id_loc(1):par_id_loc(2),:,&
                    &norm_id(1):norm_id(2),sXr_loc(1,2):sXr_loc(2,2))),&
                    &[loc_size(2)])
                
                ! IM_PV_1
                X_1D_loc => X_1D(id); id = id+1
                X_1D_loc%var_name = 'IM_PV_1'
                allocate(X_1D_loc%tot_i_min(4),X_1D_loc%tot_i_max(4))
                allocate(X_1D_loc%loc_i_min(4),X_1D_loc%loc_i_max(4))
                X_1D_loc%tot_i_min = [1,1,1,1]
                X_1D_loc%tot_i_max = [n_tot,nn_mod_tot(2)]
                X_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min,sXr_tot(1,2)]
                X_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max,&
                    &sXr_tot(2,2)]
                allocate(X_1D_loc%p(loc_size(2)))
                X_1D_loc%p = reshape(ip(X%PV_1(par_id_loc(1):par_id_loc(2),:,&
                    &norm_id(1):norm_id(2),sXr_loc(1,2):sXr_loc(2,2))),&
                    &[loc_size(2)])
                
                ! RE_KV_1
                X_1D_loc => X_1D(id); id = id+1
                X_1D_loc%var_name = 'RE_KV_1'
                allocate(X_1D_loc%tot_i_min(4),X_1D_loc%tot_i_max(4))
                allocate(X_1D_loc%loc_i_min(4),X_1D_loc%loc_i_max(4))
                X_1D_loc%tot_i_min = [1,1,1,1]
                X_1D_loc%tot_i_max = [n_tot,nn_mod_tot(2)]
                X_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min,sXr_tot(1,2)]
                X_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max,&
                    &sXr_tot(2,2)]
                allocate(X_1D_loc%p(loc_size(2)))
                X_1D_loc%p = reshape(rp(X%KV_1(par_id_loc(1):par_id_loc(2),:,&
                    &norm_id(1):norm_id(2),sXr_loc(1,2):sXr_loc(2,2))),&
                    &[loc_size(2)])
                
                ! IM_KV_1
                X_1D_loc => X_1D(id); id = id+1
                X_1D_loc%var_name = 'IM_KV_1'
                allocate(X_1D_loc%tot_i_min(4),X_1D_loc%tot_i_max(4))
                allocate(X_1D_loc%loc_i_min(4),X_1D_loc%loc_i_max(4))
                X_1D_loc%tot_i_min = [1,1,1,1]
                X_1D_loc%tot_i_max = [n_tot,nn_mod_tot(2)]
                X_1D_loc%loc_i_min = [par_id(1),1,grid_trim%i_min,sXr_tot(1,2)]
                X_1D_loc%loc_i_max = [par_id(2),n_tot(2),grid_trim%i_max,&
                    &sXr_tot(2,2)]
                allocate(X_1D_loc%p(loc_size(2)))
                X_1D_loc%p = reshape(ip(X%KV_1(par_id_loc(1):par_id_loc(2),:,&
                    &norm_id(1):norm_id(2),sXr_loc(1,2):sXr_loc(2,2))),&
                    &[loc_size(2)])
            end if
        end do
        
        ! write
        ierr = print_HDF5_arrs(X_1D(1:id-1),PB3D_name,trim(data_name),&
            &rich_lvl=rich_lvl,ind_print=.not.grid_trim%divided)
        CHCKERR('')
        
        ! clean up
        call grid_trim%dealloc()
        call dealloc_var_1D(X_1D)
        nullify(X_1D_loc)
        
        ! user output
        call lvl_ud(-1)
    end function print_output_X_2
    
    !> Initializes some variables concerning the mode numbers.
    !!
    !! Setup  minimum and  maximum  of mode  numbers at  every  flux surface  in
    !! equilibrium coordinates
    !!  - \c min_n_X, \c max_n_X
    !!  - \c min_m_X, \c max_m_X,
    !!
    !! It  functions depending  on  the \c  X_style used:  1  (prescribed) or  2
    !! (fast). For the fast  style, at every flux surface the  range of modes is
    !! sought that resonates most:
    !!  - \f$m = nq \pm \frac{N}{2}\f$      for poloidal flux
    !!  - \f$n = \iota m \pm \frac{N}{2}\f$  for toroidal flux
    !!
    !! where \f$N\f$ = \c n_mod_X.
    !!
    !! However, at  the same  time, both \c  m and  \c n have  to be  larger, in
    !! absolute  value, than  \c min_sec_X.  Therefore,  the range  of width  \c
    !! n_mod_X can be shifted upwards.
    !!
    !! \return ierr
    integer function init_modes(grid_eq,eq) result(ierr)
        use num_vars, only: use_pol_flux_F, X_style
        use X_vars, only: prim_X, min_sec_X, max_sec_X, n_mod_X, min_n_X, &
            &max_n_X, min_m_X, max_m_X, min_nm_X
        use MPI_utilities, only: get_ser_var
        use grid_utilities, only: trim_grid
        
        character(*), parameter :: rout_name = 'init_modes'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(eq_1_type), intent(in) :: eq                                       !< flux equilibrium
        
        ! local variables
        type(grid_type) :: grid_eq_trim                                         ! trimmed version of equilibrium
        real(dp), allocatable :: jq_tot(:)                                      ! saf. fac. or rot. transf. and derivs. in Flux coords.
        integer :: norm_eq_id(2)                                                ! untrimmed normal indices for trimmed grids
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Initializing perturbation modes variables')
        call lvl_ud(1)
        
        ! get trimmed grids
        ierr = trim_grid(grid_eq,grid_eq_trim,norm_eq_id)
        CHCKERR('')
        
        ! (re)allocate variables
        if (allocated(min_n_X)) deallocate(min_n_X)
        if (allocated(max_n_X)) deallocate(max_n_X)
        if (allocated(min_m_X)) deallocate(min_m_X)
        if (allocated(max_m_X)) deallocate(max_m_X)
        allocate(min_n_X(grid_eq_trim%n(3)),max_n_X(grid_eq_trim%n(3)))
        allocate(min_m_X(grid_eq_trim%n(3)),max_m_X(grid_eq_trim%n(3)))
        
        ! get serial version of safety factor or rot. transform
        allocate(jq_tot(grid_eq_trim%n(3)))
        if (use_pol_flux_F) then
            if (grid_eq%divided) then
                ierr = get_ser_var(eq%q_saf_FD(norm_eq_id(1):norm_eq_id(2),0),&
                    &jq_tot,scatter=.true.)
                CHCKERR('')
            else
                jq_tot = eq%q_saf_FD(:,0)
            end if
        else
            if (grid_eq%divided) then
                ierr = get_ser_var(eq%rot_t_FD(norm_eq_id(1):norm_eq_id(2),0),&
                    &jq_tot,scatter=.true.)
                CHCKERR('')
            else
                jq_tot = eq%rot_t_FD(:,0)
            end if
        end if
        
        ! set the limits depending on the X style
        select case (X_style)
            case (1)                                                            ! prescribed
                if (use_pol_flux_F) then
                    min_n_X = prim_X
                    max_n_X = prim_X
                    min_m_X = min_sec_X
                    max_m_X = max_sec_X
                else
                    min_n_X = min_sec_X
                    max_n_X = max_sec_X
                    min_m_X = prim_X
                    max_m_X = prim_X
                end if
            case (2)                                                            ! fast
                if (use_pol_flux_F) then
                    min_n_X = prim_X
                    max_n_X = prim_X
                    if (prim_X*jq_tot(1).gt.0) then                             ! nq > 0: m > 0
                        min_m_X = nint(prim_X*jq_tot-(n_mod_X-1)*0.5)
                        min_m_X = max(min_nm_X,min_m_X)
                        max_m_X = min_m_X + n_mod_X - 1
                    else                                                        ! nq < 0: m < 0
                        max_m_X = nint(prim_X*jq_tot+(n_mod_X-1)*0.5)
                        max_m_X = min(-min_nm_X,max_m_X)
                        min_m_X = max_m_X - n_mod_X + 1
                    end if
                else
                    if (prim_X*jq_tot(1).gt.0) then                             ! m iota > 0: n > 0
                        min_n_X = nint(prim_X*jq_tot-(n_mod_X-1)*0.5)
                        min_n_X = max(min_nm_X,min_n_X)
                        max_n_X = min_n_X + n_mod_X - 1
                    else                                                        ! m iota < 0: n < 0
                        max_n_X = nint(prim_X*jq_tot+(n_mod_X-1)*0.5)
                        max_n_X = min(-min_nm_X,max_n_X)
                        min_n_X = max_n_X - n_mod_X + 1
                    end if
                    min_m_X = prim_X
                    max_m_X = prim_X
                end if
        end select
        
        ! clean up
        call grid_eq_trim%dealloc()
        
        call lvl_ud(-1)
        call writo('Perturbation modes variables initialized')
    end function init_modes
    
    !> Sets up some variables concerning the mode numbers.
    !!
    !! Apart from the mode numbers
    !!  - \c n, \c m,
    !!
    !! also the indices of the secondary modes
    !!  - sec_ind,
    !! 
    !! are set up  in the coordinates of  the grid passed. The  normal values of
    !! this grid are saved as well, in
    !!  - r_F.
    !!
    !! This procedure makes use of the global variables
    !!  - \c min_m_X, max_m_X, min_n_X, max_m_X
    !!
    !! that have to be set up using init_nm_x().
    !!
    !! Optionally, \c n and \c m can be plot.
    !!
    !! \return ierr
    integer function setup_modes(mds,grid_eq,grid,plot_nm) result(ierr)
        use num_vars, only: use_pol_flux_F, rank, norm_disc_prec_X
        use X_vars, only: n_mod_X, min_n_X, max_n_X, min_m_X, max_m_X
        use grid_vars, only: disc_type
        use grid_utilities, only: trim_grid, setup_interp_data, apply_disc
        use eq_vars, only: max_flux_F
        
        character(*), parameter :: rout_name = 'setup_modes'
        
        ! input / output
        type(modes_type), intent(inout) :: mds                                  !< modes variables
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(grid_type), intent(in) :: grid                                     !< grid at which to calculate modes
        logical, intent(in), optional :: plot_nm                                !< plot \c n and \c m
        
        ! local variables
        type(grid_type) :: grid_eq_trim                                         ! trimmed version of equilibrium
        type(grid_type) :: grid_trim                                            ! trimmed version of perturbation grid
        type(disc_type) :: norm_interp_data                                     ! data for normal interpolation
        real(dp), allocatable :: lim_nm_X(:,:)                                  ! bundled mode number limits
        real(dp), allocatable :: lim_nm_X_interp(:,:)                           ! interpolated lim_nm_X
        real(dp), allocatable :: x_plot(:,:)                                    ! x values of plot
        integer :: norm_eq_id(2)                                                ! untrimmed normal indices for trimmed grids
        integer :: id, ld, kd                                                   ! counters
        integer :: n_mod_tot                                                    ! total number of modes, not equal to n_mod_X for X_style 2
        character(len=max_str_ln) :: plot_title                                 ! title for plots
        character(len=max_str_ln) :: plot_name                                  ! file name for plots
        logical :: plot_nm_loc                                                  ! local plot_nm
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Setting up general modes variables')
        call lvl_ud(1)
        
        ! set local plot_nm
        plot_nm_loc = .false.
        if (present(plot_nm)) plot_nm_loc = plot_nm
        
        ! get trimmed grids
        ierr = trim_grid(grid_eq,grid_eq_trim,norm_eq_id)
        CHCKERR('')
        ierr = trim_grid(grid,grid_trim)
        CHCKERR('')
        
        ! bundle the mode limits in one variable and convert to real
        allocate(lim_nm_X(4,grid_eq_trim%n(3)))
        lim_nm_X(1,:) = min_n_X*1._dp
        lim_nm_X(2,:) = max_n_X*1._dp
        lim_nm_X(3,:) = min_m_X*1._dp
        lim_nm_X(4,:) = max_m_X*1._dp
        
        ! setup normal interpolation data
        ierr = setup_interp_data(grid_eq_trim%r_F,grid_trim%r_F,&
            &norm_interp_data,norm_disc_prec_X)
        CHCKERR('')
        
        ! interpolate
        allocate(lim_nm_X_interp(4,grid_trim%n(3)))
        ierr = apply_disc(lim_nm_X,norm_interp_data,lim_nm_X_interp,2)
        CHCKERR('')
        
        ! clean up
        call norm_interp_data%dealloc()
        
        ! set up normal tabulation values
        allocate(mds%r_F(grid_trim%n(3)))
        mds%r_F = grid_trim%r_F
        
        ! set n and m
        allocate(mds%n(grid_trim%n(3),n_mod_X))
        allocate(mds%m(grid_trim%n(3),n_mod_X))
        if (use_pol_flux_F) then
            do kd = 1,grid_trim%n(3)
                mds%n(kd,:) = nint(lim_nm_X_interp(1,kd))
                mds%m(kd,:) = [(ld, ld = nint(lim_nm_X_interp(3,kd)),&
                    &nint(lim_nm_X_interp(4,kd)))]
            end do
        else
            do kd = 1,grid_trim%n(3)
                mds%n(kd,:) = [(ld, ld = nint(lim_nm_X_interp(1,kd)),&
                    &nint(lim_nm_X_interp(2,kd)))]
                mds%m(kd,:) = nint(lim_nm_X_interp(3,kd))
            end do
        end if
        
        ! total number of mds that might occur
        if (use_pol_flux_F) then
            n_mod_tot = maxval(max_m_X)-minval(min_m_X)+1
        else
            n_mod_tot = maxval(max_n_X)-minval(min_n_X)+1
        end if
        
        ! set up indices of n and m
        allocate(mds%sec_ind(grid_trim%n(3),n_mod_tot))
        
        ! loop over all possible mds to find possible indices
        mds%sec_ind = 0
        do kd = 1,grid_trim%n(3)
            do ld = 1,n_mod_tot
                do id = 1,n_mod_X
                    if (use_pol_flux_F) then
                        if (mds%m(kd,id).eq.minval(min_m_X)+ld-1) &
                            &mds%sec_ind(kd,ld) = id
                    else
                        if (mds%n(kd,id).eq.minval(min_n_X)+ld-1) &
                            &mds%sec_ind(kd,ld) = id
                    end if
                end do
            end do
        end do
        
        ! master plots output if requested
        if (rank.eq.0 .and. plot_nm_loc) then
            call writo('Plotting mode numbers')
            call lvl_ud(1)
            
            ! set up x values
            allocate(x_plot(grid_trim%n(3),n_mod_X))
            do ld = 1,n_mod_X
                x_plot(:,ld) = grid_trim%r_F
            end do
            x_plot = x_plot*2*pi/max_flux_F
            ! plot poloidal modes
            plot_title = 'poloidal mode numbers'
            plot_name = 'modes_m_X'
            call print_ex_2D([plot_title],plot_name,mds%m*1._dp,x=x_plot,&
                &draw=.false.)
            call draw_ex([plot_title],plot_name,n_mod_X,1,.false.)
            ! plot toroidal modes
            plot_title = 'toroidal mode numbers'
            plot_name = 'modes_n_X'
            call print_ex_2D([plot_title],plot_name,mds%n*1._dp,x=x_plot,&
                &draw=.false.)
            call draw_ex([plot_title],plot_name,n_mod_X,1,.false.)
            
            call lvl_ud(-1)
        end if
        
        ! clean up
        call grid_eq_trim%dealloc()
        call grid_trim%dealloc()
        
        call lvl_ud(-1)
        call writo('General mode variables set up')
    end function setup_modes
    
    !> Checks whether the high-n approximation is valid:
    !!
    !! This depends on the \c X_style:
    !!
    !! For X_style 1 (prescribed): every mode should resonate at least somewhere
    !! in the whole normal range:
    !!  - \f$\frac{\left|nq-m\right|}{\left|n\right|} < T\f$ and
    !!      \f$\frac{\left|nq-m\right|}{\left|m\right|} < T\f$
    !!      for poloidal flux
    !!  - \f$\frac{\left|q-\iota m\right|}{\left|m\right|} < T\f$ and
    !!      \f$\frac{\left|n-\iota m\right|}{\left|n\right|} < T\f$
    !!      for toroidal flux
    !!
    !! where \f$T\f$ = \c tol << 1.
    !!
    !! This condition is determined by the sign of \f$q\f$ (or \f$\iota\f$) and
    !! given by:
    !!  - \f$\max\left(q_\text{min}-T,\frac{q_\text{min}}{1+T}\right) <
    !!      \frac{m}{n} <
    !!      \min\left(q_\text{max}+T,\frac{q_\text{max}}{1-T}\right)\f$,
    !!      \f$q>0\f$
    !!  - \f$\max\left(q_\text{min}-T,\frac{q_\text{min}}{1-T}\right) <
    !!      \frac{m}{n} <
    !!      \min\left(q_\text{max}+T,\frac{q_\text{max}}{1+T}\right)\f$,
    !!      \f$q<0\f$
    !!
    !! for poloidal flux.
    !!
    !! For  toroidal  flux,  \f$q\f$  should  be  replaced  by  \f$\iota\f$  and
    !! \f$\frac{m}{n}\f$ by \f$\frac{n}{m}\f$.
    !!
    !! For \c  X_style 2 (fast):  the resonance has been  taken care of,  but it
    !! remains to be checked whether the number of modes is efficient.
    !!
    !! \return ierr
    integer function check_X_modes(grid_eq,eq) result(ierr)
        use num_vars, only: X_style, rank, use_pol_flux_F, n_procs
        use MPI_utilities, only: get_ser_var
        
        character(*), parameter :: rout_name = 'check_X_modes'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(eq_1_type), intent(in) :: eq                                       !< flux equilibrium
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        call writo('Checking mode numbers')
        call lvl_ud(1)
        
        ! check the modes depending on the X style
        select case (X_style)
            case (1)                                                            ! prescribed
                ierr = check_X_modes_1(eq)
                CHCKERR('')
            case (2)                                                            ! fast
                ierr = check_X_modes_2(grid_eq,eq)
                CHCKERR('')
        end select
        
        call lvl_ud(-1)
        call writo('Mode numbers checked')
    contains
        ! Version for  X style 1: Checks  whether there exists a  range in which
        ! each of the modes resonates, with a certain tolerance tol_norm.
        !> \private
        integer function check_X_modes_1(eq) result(ierr)
            use X_vars, only: prim_X, min_sec_X, n_mod_X
            use num_vars, only: tol_norm
            
            character(*), parameter :: rout_name = 'check_X_modes_1'
            
            ! input / output
            type(eq_1_type), intent(in) :: eq                                   ! flux equilibrium
            
            ! local variables
            integer :: id                                                       ! counter
            integer :: pmone                                                    ! plus or minus one
            real(dp) :: min_jq, max_jq                                          ! min. and max. values of q or iota
            real(dp) :: lim_lo, lim_hi                                          ! lower and upper limit on n/m (or m/n)
            real(dp), allocatable :: temp_var(:)                                ! temporary variable
            character(len=max_str_ln) :: jq_name                                ! either safety factor or rotational transform
            character(len=1) :: mode_name                                       ! either n or m
            logical :: all_modes_fine                                           ! all modes are fine
            
            ! initialize ierr
            ierr = 0
            
            ! user output
            call writo('The tolerance used is '//trim(r2strt(tol_norm))//'...')
            
            ! set min_jq and max_jq in Flux coordinate system
            if (use_pol_flux_F) then 
                min_jq = minval(eq%q_saf_FD(:,0))
                max_jq = maxval(eq%q_saf_FD(:,0))
            else
                min_jq = minval(eq%rot_t_FD(:,0))
                max_jq = maxval(eq%rot_t_FD(:,0))
            end if
            
            ! combine all processes
            ierr = get_ser_var([min_jq],temp_var)
            CHCKERR('')
            if (rank.eq.0) min_jq = minval(temp_var)
            ierr = get_ser_var([max_jq],temp_var)
            CHCKERR('')
            if (rank.eq.0) max_jq = maxval(temp_var)
            
            ! output if master
            if (rank.eq.0) then
                ! set up jq name and all_modes_fine
                if (use_pol_flux_F) then
                    jq_name = 'safety factor'
                    mode_name = 'm'
                else
                    jq_name = 'rotational transform'
                    mode_name = 'n'
                end if
                all_modes_fine = .true.
                
                ! set up plus minus one, according to the sign of jq
                if (min_jq.lt.0 .and. max_jq.lt.0) then
                    pmone = -1
                else if (min_jq.gt.0 .and. max_jq.gt.0) then
                    pmone = 1
                else
                    err_msg = trim(jq_name)//' cannot change sign'
                    ierr = 1
                    CHCKERR(err_msg)
                end if
                
                ! calculate upper and lower limits
                lim_lo = max(min_jq-tol_norm,min_jq/(1+pmone*tol_norm))
                lim_hi = min(max_jq+tol_norm,max_jq/(1-pmone*tol_norm))
                
                ! for every mode (n,m) check whether  m/n is inside the range of
                ! q values or n/m inside the range of iota values
                do id = 1,n_mod_X
                    ! check if limits are met
                    if (use_pol_flux_F) then
                        if ((min_sec_X+id-1._dp)/prim_X.lt.lim_lo .or. &
                            &(min_sec_X+id-1._dp)/prim_X.gt.lim_hi) then
                            call writo('for (n,m) = ('//trim(i2str(prim_X))//&
                                &','//trim(i2str(min_sec_X+id-1))//&
                                &'), there is no range in the plasma where the &
                                &ratio |n q - m| << |n|,|m| is met',alert=.true.)
                            all_modes_fine = .false.
                        end if
                    else
                        if ((min_sec_X+id-1._dp)/prim_X.lt.lim_lo .or. &
                            &(min_sec_X+id-1._dp)/prim_X.gt.lim_hi) then
                            call writo('for (n,m) = ('//&
                                &trim(i2str(min_sec_X+id-1))//','//&
                                &trim(i2str(prim_X))//'), there is no range in &
                                &the plasma where the ratio |n - iota m| << &
                                &|m|,|n| is met',alert=.true.)
                            all_modes_fine = .false.
                        end if
                    end if
                end do
                
                ! output message
                if (all_modes_fine) then
                    call writo('The modes are all within the allowed range &
                        &of '//trim(i2str(ceiling(prim_X*lim_lo)))//' < '//&
                        &mode_name//' < '//trim(i2str(floor(prim_X*lim_hi)))//&
                        &'...')
                else
                    call writo('Not all modes are within the allowed range &
                        &of '//trim(i2str(ceiling(prim_X*lim_lo)))//' < '//&
                        &mode_name//' < '//trim(i2str(floor(prim_X*lim_hi)))//&
                        &'...')
                end if
            end if
        end function check_X_modes_1
        
        ! version for X style 2: Check  how efficient the chosen number of modes
        ! is.
        !> \private
        integer function check_X_modes_2(grid_eq,eq) result(ierr)
            use X_vars, only: min_n_X, min_m_X, n_mod_X
            use grid_utilities, only: trim_grid
            
            character(*), parameter :: rout_name = 'check_X_modes_2'
            
            ! input / output
            type(grid_type), intent(in) :: grid_eq                              ! equilibrium grid
            type(eq_1_type), intent(in) :: eq                                   ! flux equilibrium
            
            ! local variables
            type(grid_type) :: grid_eq_trim                                     ! trimmed equilibrium grid
            integer :: norm_id(2)                                               ! untrimmed normal indices for trimmed grid
            integer :: jd, kd, ld                                               ! counters
            real(dp), allocatable :: max_frac(:,:)                              ! maximum fraction
            real(dp), allocatable :: max_frac_sum(:)                            ! sum of maximum fraction for all processes
            real(dp), allocatable :: max_frac_max(:)                            ! maximum maximum fraction for all processes
            real(dp), allocatable :: max_frac_temp(:)                           ! auxilliary variable
            real(dp), allocatable :: n(:,:), m(:,:)                             ! n and m
            real(dp), allocatable :: fac_n(:), fac_m(:)                         ! factors to be multiplied with n m m
            character(len=max_str_ln) :: frac_name                              ! name of fraction
#if ldebug
            real(dp), allocatable :: x_plot(:,:)                                ! for plotting
            character(len=max_str_ln) :: plot_title                             ! title for plots
            character(len=max_str_ln) :: plot_name                              ! file name for plots
#endif
            
            ! initialize ierr
            ierr = 0
            
            ! set up helper variables
            allocate(max_frac(grid_eq%loc_n_r,n_mod_X))
            allocate(fac_n(grid_eq%loc_n_r))
            allocate(fac_m(grid_eq%loc_n_r))
            allocate(n(grid_eq%loc_n_r,n_mod_X))
            allocate(m(grid_eq%loc_n_r,n_mod_X))
            if (use_pol_flux_F) then
                do jd = 1,n_mod_X
                    n(:,jd) = min_n_X(grid_eq%i_min:grid_eq%i_max)
                    m(:,jd) = min_m_X(grid_eq%i_min:grid_eq%i_max)+jd-1
                end do
                fac_n = eq%q_saf_FD(:,0)
                fac_m = 1._dp
            else
                do jd = 1,n_mod_X
                    n(:,jd) = min_n_X(grid_eq%i_min:grid_eq%i_max)+jd-1
                    m(:,jd) = min_m_X(grid_eq%i_min:grid_eq%i_max)
                end do
                fac_n = 1._dp
                fac_m = eq%rot_t_FD(:,0)
            end if
            
            ! loop over all modes
            do ld = 1,n_mod_X
                do kd = 1,grid_eq%loc_n_r
                    max_frac(kd,ld) = &
                        &abs(n(kd,ld)*fac_n(kd)-m(kd,ld)*fac_m(kd))/&
                        &min(abs(n(kd,ld)),abs(m(kd,ld)))
                end do
            end do
            
            ! trim grid
            ierr = trim_grid(grid_eq,grid_eq_trim,norm_id=norm_id)
            CHCKERR('')
            
            ! gather all fractions
            if (rank.eq.0) allocate(max_frac_temp(n_procs))
            if (rank.eq.0) allocate(max_frac_sum(n_mod_X))
            if (rank.eq.0) allocate(max_frac_max(n_mod_X))
            do ld = 1,n_mod_X
                ierr = get_ser_var([sum(max_frac(norm_id(1):norm_id(2),ld))],&
                    &max_frac_temp)
                CHCKERR('')
                if (rank.eq.0) max_frac_sum(ld) = sum(max_frac_temp)
                ierr = get_ser_var(&
                    &[maxval(max_frac(norm_id(1):norm_id(2),ld))],max_frac_temp)
                CHCKERR('')
                if (rank.eq.0) max_frac_max(ld) = maxval(max_frac_temp)
            end do
            
            ! output if master
            if (rank.eq.0) then
                if (use_pol_flux_F) then
                    frac_name = '|nq-m|/|n| or |nq-m|/|m|'
                else
                    frac_name = '|n-iota m|/|n| or |n-iota m|/|m|'
                end if
                call writo('The fraction '//trim(frac_name)//' is')
                call lvl_ud(1)
                do ld = 1,n_mod_X
                    call writo('for mode '//trim(i2str(ld))//' maximally '//&
                        &trim(r2strt(max_frac_max(ld)))//&
                        &', average '//&
                        &trim(r2strt(max_frac_sum(ld)/grid_eq%n(3))))
                end do
                if (n_mod_X.gt.1) call writo('so for all modes maximally '//&
                    &trim(r2strt(maxval(max_frac_max)))//', average '//&
                    &trim(r2strt(sum(max_frac_sum)/(grid_eq%n(3)*n_mod_X))))
                call lvl_ud(-1)
            end if
#if ldebug
            if (debug_check_X_modes_2) then
                call writo('Plotting the fraction for all modes')
                allocate(x_plot(grid_eq%n(3),n_mod_X))
                do ld = 1,n_mod_X
                    x_plot(:,ld) = grid_eq%r_F
                end do
                plot_name = 'TEST_max_frac'
                plot_title = 'maximum fraction'
                call print_ex_2D([plot_title],plot_name,max_frac,x=x_plot,&
                    &draw=.false.)
                call draw_ex([plot_title],plot_name,n_mod_X,1,.false.)
            end if
#endif
            
            ! clean up
            call grid_eq_trim%dealloc()
        end function check_X_modes_2
    end function check_X_modes
    
    !> Calculates resonating flux surfaces for the perturbation modes.
    !!
    !! The output  consists of mode  number, resonating normal position  and the
    !! fraction  \f$\frac{m}{n}\f$ or  \f$\frac{n}{m}\f$.  for  those modes  for
    !! which a solution is found that is within the plasma range.
    !!
    !! It contains three pieces of information:
    !!  - <tt>(:,1)</tt>: the mode index
    !!  - <tt>(:,2)</tt>: the radial position in Flux coordinates
    !!  - <tt>(:,3)</tt>: the fraction \f$\frac{m}{n}\f$ or  \f$\frac{n}{m}\f$
    !!
    !! for every single mode in \c sec_ind of \c mds.
    !!
    !! Optionally,  the  total safety  factor  or  rotational transform  can  be
    !! returned to the master.
    !!
    !! Also, information can be displayed.
    !!
    !! \return ierr
    integer function calc_res_surf(mds,grid_eq,eq,res_surf,info,jq) result(ierr)
        use X_vars, only: min_n_X, min_m_X, prim_X
        use num_vars, only: use_pol_flux_F, norm_disc_prec_eq
        use eq_vars, only: max_flux_F, max_flux_F
        use num_ops, only: calc_zero_Zhang
        use grid_utilities, only: trim_grid, setup_interp_data, apply_disc
        use MPI_utilities, only: get_ser_var
        
        character(*), parameter :: rout_name = 'calc_res_surf'
        
        ! input / output
        type(modes_type), intent(in) :: mds                                     !< general modes variables
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid_eq
        type(eq_1_type), intent(in) :: eq                                       !< flux equilibrium
        real(dp), intent(inout), allocatable :: res_surf(:,:)                   !< resonant surface
        logical, intent(in), optional :: info                                   !< if info is displayed
        real(dp), intent(inout), optional, allocatable :: jq(:)                 !< either safety factor or rotational transform
        
        ! local variables (also used in child functions)
        real(dp) :: nmfrac_fun                                                  ! fraction m/n or n/m to determine resonant flux surface
        type(disc_type) :: norm_interp_data                                     ! data for normal interpolation
        
        ! local variables (not to be used in child functions)
        integer :: norm_id(2)                                                   ! untrimmed normal indices for trimmed grid
        integer :: ld                                                           ! counter
        integer :: ld_loc                                                       ! local ld
        logical :: info_loc                                                     ! local version of info
        real(dp), allocatable :: res_surf_loc(:,:)                              ! local copy of res_surf
        real(dp), allocatable :: jq_tot(:)                                      ! saf. fac. or rot. transf. and derivs. in Flux coords.
        real(dp), allocatable :: jq_loc(:)                                      ! local version of jq
        real(dp) :: norm_factor                                                 ! normalization factor to make normal coordinate 0..1
        type(grid_type) :: grid_eq_trim                                         ! trimmed version of equilibrium grid
        integer :: n_loc, m_loc                                                 ! local n and m
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set local info
        info_loc = .false.
        if (present(info)) info_loc = info
        
        ! get trimmed grid
        ierr = trim_grid(grid_eq,grid_eq_trim,norm_id)
        CHCKERR('')
        
        ! get serial version of safety factor or rot. transform
        allocate(jq_tot(grid_eq_trim%n(3)))
        if (use_pol_flux_F) then
            if (grid_eq%divided) then
                ierr = get_ser_var(eq%q_saf_FD(norm_id(1):norm_id(2),0),&
                    &jq_loc,scatter=.true.)
                CHCKERR('')
                jq_tot = jq_loc
            else
                jq_tot = eq%q_saf_FD(:,0)
            end if
        else
            if (grid_eq%divided) then
                ierr = get_ser_var(eq%rot_t_FD(norm_id(1):norm_id(2),0),&
                    &jq_loc,scatter=.true.)
                CHCKERR('')
                jq_tot = jq_loc
            else
                jq_tot = eq%rot_t_FD(:,0)
            end if
        end if
        if (present(jq)) then
            allocate(jq(grid_eq_trim%n(3)))
            jq = jq_tot
        end if
        
        ! initialize local res_surf
        allocate(res_surf_loc(1:size(mds%sec_ind,2),3))
        
        ! calculate normalization factor max_flux / 2pi
        norm_factor = max_flux_F/(2*pi)
        
        ! loop over all modes (and shift the index in m_loc and n_loc by 1)
        ld_loc = 1
        do ld = 1,size(mds%sec_ind,2)
            ! Find place  where q  = m/n or  iota = n/m  in Flux  coordinates by
            ! solving q-m/n  = 0 or  iota-n/m=0, using the functin  jq_fun. Only
            ! consider modes  that appear  somewhere in  the plasma  (e.g. where
            ! sec_ind is nonzero somewhere)
            
            ! set up local m, n and nmfrac for function if mode appears
            if (maxval(mds%sec_ind(:,ld)).gt.0) then
                if (use_pol_flux_F) then
                    n_loc = prim_X
                    m_loc = minval(min_m_X)+ld-1
                    nmfrac_fun = 1.0_dp*m_loc/n_loc
                else
                    n_loc = minval(min_n_X)+ld-1
                    m_loc = prim_X
                    nmfrac_fun = 1.0_dp*n_loc/m_loc
                end if
            else
                cycle
            end if
            
            ! calculate zero using Zhang
            if (use_pol_flux_F) then
                res_surf_loc(ld_loc,1) = m_loc
            else
                res_surf_loc(ld_loc,1) = n_loc
            end if
            res_surf_loc(ld_loc,3) = nmfrac_fun
            err_msg = calc_zero_Zhang(res_surf_loc(ld_loc,2),jq_fun,&
                &[minval(grid_eq_trim%r_F),maxval(grid_eq_trim%r_F)])
            
            ! intercept error
            if (err_msg.ne.'') then
                if (info_loc) then
                    call writo('Mode (n,m) = ('//trim(i2str(n_loc))//','//&
                        &trim(i2str(m_loc))//') is not found')
                    call lvl_ud(1)
                    call writo(trim(err_msg))
                    call lvl_ud(-1)
                end if
            else if (res_surf_loc(ld_loc,2).lt.minval(grid_eq_trim%r_F) .or. &
                &res_surf_loc(ld_loc,2).gt.maxval(grid_eq_trim%r_F)) then
                if (info_loc) call writo('Mode (n,m) = ('//&
                    &trim(i2str(n_loc))//','//trim(i2str(m_loc))//&
                    &') does not resonate in plasma')
            else
                if (info_loc) call writo('Mode (n,m) = ('//&
                    &trim(i2str(n_loc))//','//trim(i2str(m_loc))//&
                    &') resonates in plasma at normalized flux surface '//&
                    &trim(r2str(res_surf_loc(ld_loc,2)/norm_factor)))
                ld_loc = ld_loc + 1                                             ! advance ld_loc
            end if
        end do
        
        ! set res_surf from local copy
        allocate(res_surf(ld_loc-1,3))
        res_surf = res_surf_loc(1:ld_loc-1,:)
        
        ! clean up
        call grid_eq_trim%dealloc()
        call norm_interp_data%dealloc()
    contains
        ! Returns q-m/n or  iota-n/m in Flux coordinates, used to  solve for q =
        ! m/n or iota = n/m.
        !> \private
        real(dp) function jq_fun(pt) result(res)
            ! input / output
            real(dp), intent(in) :: pt                                          ! normal position at which to evaluate
            
            ! local variables
            integer :: i_min, i_max                                             ! index of minimum and maximum value of x
            real(dp) :: x_min, x_max                                            ! minimum and maximum value of x
            real(dp) :: res_loc(1)                                              ! local copy of res
            
            ! initialize res
            res = 0
            
            ! set up min. and max index and value
            x_min = minval(grid_eq_trim%r_F)
            x_max = maxval(grid_eq_trim%r_F)
            i_min = minloc(grid_eq_trim%r_F,1)
            i_max = maxloc(grid_eq_trim%r_F,1)
            
            ! check whether to interpolate or extrapolate
            if (pt.ge.x_min .and. pt.le.x_max) then                             ! point requested within range
                ! setup interpolation data
                ierr = setup_interp_data(grid_eq_trim%r_F,[pt],&
                    &norm_interp_data,norm_disc_prec_eq)
                CHCKERR('')
                ! interpolate
                ierr = apply_disc(jq_tot-nmfrac_fun,norm_interp_data,&
                    &res_loc)
                CHCKERR('')
                ! copy
                res = res_loc(1)
            end if
        end function jq_fun
    end function calc_res_surf
    
    !> plot  \f$q\f$-profile or  \f$\iota\f$-profile  in  flux coordinates  with
    !! \f$nq-m = 0\f$ or \f$n-\iota m = 0\f$ indicated if requested.
    !!
    !! \return ierr
    integer function resonance_plot(mds,grid,eq) result(ierr)
        use num_vars, only: use_pol_flux_F, no_plots, n_theta_plot, &
            &n_zeta_plot, rank, eq_style, min_theta_plot, max_theta_plot, &
            &min_zeta_plot, max_zeta_plot
        use eq_vars, only: max_flux_F
        use grid_utilities, only: calc_XYZ_grid, calc_eqd_grid, coord_F2E, &
            &trim_grid
        use X_vars, only: prim_X
        use VMEC_utilities, only: calc_trigon_factors
        
        character(*), parameter :: rout_name = 'resonance_plot'
        
        ! input / output
        type(modes_type), intent(in) :: mds                                     !< general modes variables
        type(grid_type), intent(in) :: grid                                     !< equilibrium grid
        type(eq_1_type), intent(in) :: eq                                       !< flux equilibrium
        
        ! local variables (not to be used in child functions)
        integer :: ld                                                           ! counter
        integer :: n_mod_loc                                                    ! local n_mod
        real(dp), allocatable :: res_surf(:,:)                                  ! resonant surfaces
        real(dp), allocatable :: x_plot_loc(:,:)                                ! for plotting
        real(dp), allocatable :: y_plot_loc(:,:)                                ! for plotting
        character(len=max_str_ln) :: plot_title, file_name                      ! name of plot, of file
        real(dp), allocatable :: jq(:)                                          ! saf. fac. or rot. transf. in Flux coords.
        integer :: n_r                                                          ! total number of normal points
        integer :: plot_dim(4)                                                  ! plot dimensions (total = local because only masters)
        type(grid_type) :: grid_trim                                            ! trimmed version of grid
        type(grid_type) :: grid_plot                                            ! grid for plotting
        real(dp), allocatable :: r_plot_E(:)                                    ! normal E coordinates of plot (needed to calculate X, Y and Z)
        real(dp), allocatable :: theta_plot(:,:,:), zeta_plot(:,:,:)            ! pol. and tor. angle of plot
        real(dp), allocatable :: X_plot(:,:,:,:), Y_plot(:,:,:,:), &
            &Z_plot(:,:,:,:)                                                    ! X, Y and Z of plot of all surfaces
        real(dp), allocatable :: vars(:,:,:,:)                                  ! variable to plot
        character(len=max_str_ln), allocatable :: plot_titles(:)                ! name of plots
        
        ! initialize ierr
        ierr = 0
        
        ! bypass plots if no_plots
        if (no_plots) return
        
        ! get trimmed grid
        ierr = trim_grid(grid,grid_trim)
        CHCKERR('')
        
        ! initialize variables
        n_r = grid_trim%n(3)
        
        ! print user output
        if (use_pol_flux_F) then
            call writo('Plotting safety factor q and resonant surfaces &
                &q = m/n...')
            plot_title = 'safety factor q'
            file_name = 'q_saf'
        else
            call writo('Plotting rotational transform iota and resonant &
                &surfaces iota = n/m...')
            plot_title = 'rotational transform iota'
            file_name = 'rot_t'
        end if
        
        call lvl_ud(1)
        
        call writo('Calculate resonant surfaces')
        call lvl_ud(1)
        
        ! find resonating surfaces
        ierr = calc_res_surf(mds,grid,eq,res_surf,info=.true.,jq=jq)
        CHCKERR('')
        
        call lvl_ud(-1)
        
        ! only master
        if (rank.eq.0) then
            ! set local n_mod
            n_mod_loc = size(res_surf,1)
            
            ! set up plot titles
            allocate(plot_titles(n_mod_loc))
            if (use_pol_flux_F) then                                            ! n is fixed and m = m/n n
                do ld = 1,n_mod_loc
                    plot_titles(ld) = 'resonance for m,n = '//&
                        &trim(i2str(nint(res_surf(ld,3)*prim_X)))//','//&
                        &trim(i2str(prim_X))
                end do
            else                                                                ! m is fixed and n = n/m m
                do ld = 1,n_mod_loc
                    plot_titles(ld) = 'resonance for m,n = '//&
                        &trim(i2str(prim_X))//','//&
                        &trim(i2str(nint(res_surf(ld,3)*prim_X)))
                end do
            end if
            
            ! initialize x_plot_loc and y_plot_loc
            allocate(x_plot_loc(n_r,n_mod_loc+1)); x_plot_loc = 0
            allocate(y_plot_loc(n_r,n_mod_loc+1)); y_plot_loc = 0
            
            ! set x_plot_loc and y_plot_loc for first column
            x_plot_loc(:,1) = grid_trim%r_F
            y_plot_loc(:,1) = jq(:)
            
            ! set x_plot_loc and y_plot_loc for other columns
            do ld = 1,n_mod_loc
                x_plot_loc(:,ld+1) = res_surf(ld,2)
                y_plot_loc(n_r,ld+1) = res_surf(ld,3)
            end do
            
            ! only if resonant surfaces
            if (size(res_surf,1).gt.0) then
                ! user message
                call writo('Plot resonant flux surfaces using HDF5')
                
                call lvl_ud(1)
                
                ! set up pol. and tor. angle for plot
                allocate(theta_plot(n_theta_plot,n_zeta_plot,1))
                allocate(zeta_plot(n_theta_plot,n_zeta_plot,1))
                ierr = calc_eqd_grid(theta_plot,min_theta_plot*pi,&
                    &max_theta_plot*pi,1)
                CHCKERR('')
                ierr = calc_eqd_grid(zeta_plot,min_zeta_plot*pi,&
                    &max_zeta_plot*pi,2)
                CHCKERR('')
                
                ! set up vars
                allocate(vars(n_theta_plot,n_zeta_plot,1,n_mod_loc))
                do ld = 1,n_mod_loc
                    vars(:,:,:,ld) = y_plot_loc(n_r,ld+1)
                end do
                
                ! set dimensions for single flux surface
                plot_dim = [n_theta_plot,n_zeta_plot,1,n_mod_loc]
                
                ! calculate normal vars in Equilibrium coords.
                allocate(r_plot_E(n_mod_loc))
                ierr = coord_F2E(grid,x_plot_loc(n_r,2:n_mod_loc+1),r_plot_E,&
                    &r_F_array=grid%r_F,r_E_array=grid%r_E)
                CHCKERR('')
                
                ! create plot grid for single flux surface
                ierr = grid_plot%init(plot_dim(1:3))
                CHCKERR('')
                grid_plot%theta_E = theta_plot
                grid_plot%zeta_E = zeta_plot
                
                ! if VMEC, calculate trigonometric factors of plot grid
                if (eq_style.eq.1) then
                    ierr = calc_trigon_factors(grid_plot%theta_E,&
                        &grid_plot%zeta_E,grid_plot%trigon_factors)
                    CHCKERR('')
                end if
                
                ! set up plot X, Y and Z
                allocate(X_plot(n_theta_plot,n_zeta_plot,1,n_mod_loc))
                allocate(Y_plot(n_theta_plot,n_zeta_plot,1,n_mod_loc))
                allocate(Z_plot(n_theta_plot,n_zeta_plot,1,n_mod_loc))
                
                ! loop over all resonant surfaces to calculate X, Y and Z values
                do ld = 1,n_mod_loc
                    ! set loc_r_E of plot grid
                    grid_plot%loc_r_E = r_plot_E(ld)
                    
                    ! calculate X, Y and Z
                    ierr = calc_XYZ_grid(grid,grid_plot,X_plot(:,:,:,ld),&
                        &Y_plot(:,:,:,ld),Z_plot(:,:,:,ld))
                    CHCKERR('')
                end do
                
                ! plot using HDF5
                call plot_HDF5(plot_titles,file_name,vars,X=X_plot,Y=Y_plot,&
                    &Z=Z_plot,col=1,descr='resonant surfaces')
                
                ! deallocate local variables
                deallocate(vars)
                deallocate(theta_plot,zeta_plot,r_plot_E)
                deallocate(X_plot,Y_plot,Z_plot)
                call grid_plot%dealloc()
                
                call lvl_ud(-1)
            end if
            
            call writo('Plot safety factor with resonant flux surfaces using &
                &external program')
            
            call lvl_ud(1)
            
            ! rescale x_plot_loc by max_flux_F/2*pi
            x_plot_loc = x_plot_loc*2*pi/max_flux_F
            
            ! plot using external program
            call print_ex_2D([plot_title,plot_titles],file_name,y_plot_loc,&
                &x=x_plot_loc,draw=.false.)
            call draw_ex([plot_title,plot_titles],file_name,n_mod_loc+1,1,&
                &.false.)
            
            call lvl_ud(-1)
            
            ! only if resonant surfaces
            if (size(res_surf,1).gt.0) then
                call writo('Plot 2D map of resonances using HDF5')
                
                call lvl_ud(1)
                
                ! initialize plot dimensions
                plot_dim = [1,n_mod_loc,1,0]
                
                ! set up X, Y and Z
                allocate(X_plot(plot_dim(1),plot_dim(2),plot_dim(3),1))
                allocate(Y_plot(plot_dim(1),plot_dim(2),plot_dim(3),1))
                allocate(Z_plot(plot_dim(1),plot_dim(2),plot_dim(3),1))
                
                do ld = 1,n_mod_loc
                    X_plot(1,ld,1,1) = x_plot_loc(1,ld+1)
                    Y_plot(1,ld,1,1) = res_surf(ld,1)
                    Z_plot(1,ld,1,1) = 1._dp
                end do
                
                ! plot 2D (to be used with plots of harmonics in sol_ops)
                call plot_HDF5('resonating surfaces','res_surf',&
                    &Z_plot(:,:,:,1),x=X_plot(:,:,:,1),y=Y_plot(:,:,:,1))
                
                ! deallocate local variables
                deallocate(X_plot,Y_plot,Z_plot)
                
                call lvl_ud(-1)
            end if
        end if
        
        ! clean up
        call grid_trim%dealloc()
        
        call lvl_ud(-1)
    end function resonance_plot
    
    !> Calculate \f$U_m^0\f$, \f$U_m^1\f$ or \f$U_n^0\f$, \f$U_n^1\f$.
    !!
    !! This is done at  eq \c loc_n_r values of the  normal coordinate, \c n_par
    !! values of  the parallel coordinate and  \c size_X values of  the poloidal
    !! mode number,  or of  the toroidal  mode number, as  well as  the poloidal
    !! derivatives
    !!
    !! Use is made  of variables \c Ti  that are set up in  the equilibrium grid
    !! and are afterwards  converted to \c Ti_X in the  perturbation grid, which
    !! needs to have the same angular coordinates as the equilibrium grid:
    !!  \f[\begin{aligned}
    !!      T_1 &= \frac{B_\alpha}{B_\theta} \\
    !!      T_2 &= \Theta^\alpha + q' \theta \\
    !!      T_3 &= \frac{B_\alpha}{B_\theta} q' +
    !!          \frac{\mathcal{J}^2}{B_\theta} \mu_0 p' \\
    !!      T_4 &= \frac{B_\alpha}{B_\theta} q' \theta -
    !!          \frac{B_\psi}{B_\theta} \\
    !!      T_5 &= \frac{B_\alpha}{B_\theta}
    !!          \frac{\partial \Theta^\theta}{\partial \theta} -
    !!          \frac{\partial\Theta^\theta}{\partial \alpha} \\
    !!      T_6 &= \frac{B_\alpha}{B_\theta} \Theta^\theta \ ,
    !!  \end{aligned}\f]
    !! with \f$T_i\f$ = \c Ti.
    !!
    !! which is valid for poloidal Flux coordinates.
    !!
    !! For  toroidal   Flux  coordinates,  \f$q'\f$   has  to  be   replaced  by
    !! \f$-\iota'\f$ and \f$\theta\f$ by \f$\zeta\f$.
    !!
    !! The interpolated \c Ti_X are then used to calculate \f$U\f$:
    !!  \f[\begin{aligned}
    !!      U_0 =& -T_2                                      &\text{(order 1)}\\
    !!          &+ \frac{i}{n}\left(T_3 + i(nq-m) T_4\right) &\text{(order 2)}\\
    !!          &+ \left(\frac{i}{n}\right)^2 i(nq-m)
    !!          \left(-T_5 - (nq-m) T_6\right)               &\text{(order 3)}\\
    !!      U_1 =& \frac{i}{n}                               &\text{(order 1)}\\
    !!          &+ \left(\frac{i}{n}\right)^2 i(nq-m)(-T_1)  &\text{(order 2)} .
    !!  \end{aligned}\f]
    !!
    !! This is  valid for poloidal Flux  coordinates and where \f$n\f$  is to be
    !! replaced by  \f$m\f$ and \f$(nq-m)\f$  by \f$(n-\iota m)\f$  for toroidal
    !! Flux coordinates.
    !!
    !! For VMEC, these factors are also derived in the parallel coordinate.
    !!
    !! \see See \cite weyens2014theory.
    !!
    !! \return ierr
    integer function calc_U(grid_eq,grid_X,eq_1,eq_2,X) result(ierr)
        use num_vars, only: use_pol_flux_F, eq_style, U_style, &
            &norm_disc_prec_X, X_grid_style
        use num_utilities, only: c
        use input_utilities, only: get_log, pause_prog
        use eq_vars, only: vac_perm
        use grid_vars, only: disc_type
        use grid_utilities, only: setup_deriv_data, setup_interp_data, &
            &apply_disc
#if ldebug
        use num_vars, only: ltest
#endif
        
        character(*), parameter :: rout_name = 'calc_U'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid variables
        type(grid_type), intent(in) :: grid_X                                   !< perturbation grid variables
        type(eq_1_type), intent(in), target :: eq_1                             !< flux equilibrium
        type(eq_2_type), intent(in), target :: eq_2                             !< metric equilibrium
        type(X_1_type), intent(inout) :: X                                      !< vectorial perturbation variables
        
        ! local variables
        integer :: id, jd, kd, ld                                               ! counters
        type(disc_type) :: par_deriv_data                                       ! data for parallel derivative
        type(disc_type) :: geo_deriv_data                                       ! data for geodesic derivative
        integer :: T_size                                                       ! 2 for VMEC and 1 for HELENA
        type(disc_type) :: norm_interp_data                                     ! data for normal interpolation
        
        ! Jacobian
        real(dp), pointer :: J(:,:,:)                                           ! Jacobian
        real(dp), pointer :: D3J(:,:,:)                                         ! D_theta Jacobian
        ! lower metric factors
        real(dp), pointer :: g13(:,:,:)                                         ! g_alpha,theta
        real(dp), pointer :: D3g13(:,:,:)                                       ! D_theta g_alpha,theta
        real(dp), pointer :: g23(:,:,:)                                         ! g_psi,theta
        real(dp), pointer :: D3g23(:,:,:)                                       ! D_theta g_psi,theta
        real(dp), pointer :: g33(:,:,:)                                         ! g_theta,theta
        real(dp), pointer :: D3g33(:,:,:)                                       ! D_theta g_theta,theta
        ! upper metric factors
        real(dp), pointer :: h12(:,:,:)                                         ! h^alpha,psi
        real(dp), pointer :: D3h12(:,:,:)                                       ! D_theta h^alpha,psi
        real(dp), pointer :: h22(:,:,:)                                         ! h^psi,psi
        real(dp), pointer :: D1h22(:,:,:)                                       ! D_alpha h^psi,psi
        real(dp), pointer :: D3h22(:,:,:)                                       ! D_theta h^psi,psi
        real(dp), pointer :: D13h22(:,:,:)                                      ! D_alpha,theta h^psi,psi
        real(dp), pointer :: D33h22(:,:,:)                                      ! D_theta,theta h^psi,psi
        real(dp), pointer :: h23(:,:,:)                                         ! h^psi,theta
        real(dp), pointer :: D1h23(:,:,:)                                       ! D_alpha h^psi,theta
        real(dp), pointer :: D3h23(:,:,:)                                       ! D_theta h^psi,theta
        real(dp), pointer :: D13h23(:,:,:)                                      ! D_alpha,theta h^psi,theta
        real(dp), pointer :: D33h23(:,:,:)                                      ! D_theta,theta h^psi,theta
        ! helper variables
        real(dp), pointer :: ang_par_F(:,:,:)                                   ! equilibrium parallel angle in flux coordinates
        real(dp), pointer :: ang_geo_F(:,:,:)                                   ! equilibrium geodesical angle in flux coordinates
        real(dp), pointer :: q_saf(:), rot_t(:)                                 ! safety factor and rotational transform in X grid
        real(dp), allocatable :: djq(:)                                         ! either q' (pol. flux) or -iota' (tor. flux)
        real(dp), allocatable :: Theta_3(:,:,:), D1Theta_3(:,:,:), &
            &D3Theta_3(:,:,:)                                                   ! Theta^theta and derivatives
        real(dp), allocatable :: D13Theta_3(:,:,:), D33Theta_3(:,:,:)           ! Theta^theta and derivatives
        real(dp), allocatable :: n_frac(:,:)                                    ! nq-m (pol. flux) or n-iotam (tor. flux) for all modes
        real(dp), allocatable :: nm(:,:)                                        ! either n*A_0 (pol. flux) or m (tor.flux)
        complex(dp), allocatable :: U_corr(:,:,:,:)                             ! correction to U_0 and U_1 for a certain (n,m)
        complex(dp), allocatable :: D3U_corr(:,:,:,:)                           ! D_theta U_corr
        ! U factors
        real(dp), allocatable, target :: T1(:,:,:,:)                            ! B_alpha/B_theta and par. deriv.
        real(dp), allocatable, target :: T2(:,:,:,:)                            ! Theta^alpha + q' theta and par. deriv.
        real(dp), allocatable, target :: T3(:,:,:,:)                            ! B_alpha/B_theta q' + J^2/B_theta mu_0 p' and par. deriv.
        real(dp), allocatable, target :: T4(:,:,:,:)                            ! B_alpha/B_theta q' theta - B_psi/B_theta and par. deriv.
        real(dp), allocatable, target :: T5(:,:,:,:)                            ! B_alpha/B_theta D3Theta^theta - D1Theta^theta and par. deriv.
        real(dp), allocatable, target :: T6(:,:,:,:)                            ! B_alpha/B_theta Theta^theta and par. deriv.
        real(dp), pointer :: T1_X(:,:,:,:)                                      ! T1 in X grid and par. deriv.
        real(dp), pointer :: T2_X(:,:,:,:)                                      ! T2 in X grid and par. deriv.
        real(dp), pointer :: T3_X(:,:,:,:)                                      ! T3 in X grid and par. deriv.
        real(dp), pointer :: T4_X(:,:,:,:)                                      ! T4 in X grid and par. deriv.
        real(dp), pointer :: T5_X(:,:,:,:)                                      ! T5 in X grid and par. deriv.
        real(dp), pointer :: T6_X(:,:,:,:)                                      ! T6 in X grid and par. deriv.
        
        ! initialize ierr
        ierr = 0
        
        ! allocate variables
        ! helper variables
        allocate(nm(grid_X%loc_n_r,X%n_mod))
        allocate(q_saf(grid_X%loc_n_r))
        allocate(rot_t(grid_X%loc_n_r))
        allocate(djq(grid_eq%loc_n_r))
        allocate(Theta_3(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        allocate(D1Theta_3(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        allocate(D3Theta_3(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        allocate(D13Theta_3(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        allocate(D33Theta_3(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        allocate(n_frac(grid_X%loc_n_r,X%n_mod))
        allocate(U_corr(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r,2))
        allocate(D3U_corr(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r,2))
        ! U factors
        select case (eq_style)
            case (1)                                                            ! VMEC
                T_size = 2
            case (2)                                                            ! HELENA
                T_size = 1
        end select
        allocate(T1(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r,T_size))
        allocate(T2(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r,T_size))
        allocate(T3(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r,T_size))
        allocate(T4(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r,T_size))
        allocate(T5(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r,T_size))
        allocate(T6(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r,T_size))
        select case (X_grid_style)
            case (1)                                                            ! equilibrium
                ! do nothing
            case (2)                                                            ! solution
                allocate(q_saf(grid_X%loc_n_r))
                allocate(rot_t(grid_X%loc_n_r))
                allocate(T1_X(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r,T_size))
                allocate(T2_X(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r,T_size))
                allocate(T3_X(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r,T_size))
                allocate(T4_X(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r,T_size))
                allocate(T5_X(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r,T_size))
                allocate(T6_X(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r,T_size))
        end select
        
        ! set pointers
        if (use_pol_flux_F) then
            ang_par_F => grid_eq%theta_F
            ang_geo_F => grid_eq%zeta_F
        else
            ang_par_F => grid_eq%zeta_F
            ang_geo_F => grid_eq%theta_F
        end if
        J => eq_2%jac_FD(:,:,:,0,0,0)
        D3J => eq_2%jac_FD(:,:,:,0,0,1)
        g13 => eq_2%g_FD(:,:,:,c([1,3],.true.),0,0,0)
        D3g13 => eq_2%g_FD(:,:,:,c([1,3],.true.),0,0,1)
        g23 => eq_2%g_FD(:,:,:,c([2,3],.true.),0,0,0)
        D3g23 => eq_2%g_FD(:,:,:,c([2,3],.true.),0,0,1)
        g33 => eq_2%g_FD(:,:,:,c([3,3],.true.),0,0,0)
        D3g33 => eq_2%g_FD(:,:,:,c([3,3],.true.),0,0,1)
        h12 => eq_2%h_FD(:,:,:,c([1,2],.true.),0,0,0)
        D3h12 => eq_2%h_FD(:,:,:,c([1,2],.true.),0,0,1)
        h22 => eq_2%h_FD(:,:,:,c([2,2],.true.),0,0,0)
        D1h22 => eq_2%h_FD(:,:,:,c([2,2],.true.),1,0,0)
        D3h22 => eq_2%h_FD(:,:,:,c([2,2],.true.),0,0,1)
        D13h22 => eq_2%h_FD(:,:,:,c([2,2],.true.),1,0,1)
        D33h22 => eq_2%h_FD(:,:,:,c([2,2],.true.),0,0,2)
        h23 => eq_2%h_FD(:,:,:,c([2,3],.true.),0,0,0)
        D1h23 => eq_2%h_FD(:,:,:,c([2,3],.true.),1,0,0)
        D3h23 => eq_2%h_FD(:,:,:,c([2,3],.true.),0,0,1)
        D13h23 => eq_2%h_FD(:,:,:,c([2,3],.true.),1,0,1)
        D33h23 => eq_2%h_FD(:,:,:,c([2,3],.true.),0,0,2)
        
        ! set up helper variables in eq grid
        if (use_pol_flux_F) then
            nm = X%n
            djq = eq_1%q_saf_FD(:,1)
        else
            djq = -eq_1%rot_t_FD(:,1)
            nm = X%m
        end if
        Theta_3 = h23/h22
        select case (eq_style)
            case (1)                                                            ! VMEC
                D1Theta_3 = D1h23/h22 - h23*D1h22/(h22**2)
                D3Theta_3 = D3h23/h22 - h23*D3h22/(h22**2)
                D13Theta_3 = D13h23/h22 - (D3h23*D1h22+D1h23*D3h22)/(h22**2) &
                    &- h23*D13h22/(h22**2) + 2*h23*D3h22*D1h22/(h22**3)
                D33Theta_3 = D33h23/h22 - 2*D3h23*D3h22/(h22**2) &
                    &- h23*D33h22/(h22**2) + 2*h23*D3h22**2/(h22**3)
            case (2)                                                            ! HELENA
                ! geodesic derivative
                if (grid_eq%n(2).gt.norm_disc_prec_X+3) then                    ! need enough terms
                    do kd = 1,grid_eq%loc_n_r
                        do id = 1,grid_eq%n(1)
                            ierr = setup_deriv_data(ang_geo_F(id,:,kd),&
                                &geo_deriv_data,1,norm_disc_prec_X+1)           ! higher precision because other derivative will be taken later
                            CHCKERR('')
                            ierr = apply_disc(Theta_3(id,:,kd),&
                                &geo_deriv_data,D1Theta_3(id,:,kd))
                            CHCKERR('')
                        end do
                    end do
                    call geo_deriv_data%dealloc()
                else
                    D1Theta_3 = 0._dp
                end if
                ! parallel derivative
                if (grid_eq%n(1).gt.norm_disc_prec_X+3) then                    ! need enough terms
                    do kd = 1,grid_eq%loc_n_r
                        do jd = 1,grid_eq%n(2)
                            ierr = setup_deriv_data(ang_par_F(:,jd,kd),&
                                &par_deriv_data,1,norm_disc_prec_X+1)           ! higher precision because other derivative will be taken later
                            CHCKERR('')
                            ierr = apply_disc(Theta_3(:,jd,kd),&
                                &par_deriv_data,D3Theta_3(:,jd,kd))
                            CHCKERR('')
                        end do
                    end do
                    call par_deriv_data%dealloc()
                else
                    D3Theta_3 = 0._dp
                end if
        end select
        
        ! set up U factors in eq grid
        T1(:,:,:,1) = g13/g33
        do kd = 1,grid_eq%loc_n_r
            T2(:,:,kd,1) = h12(:,:,kd)/h22(:,:,kd) + djq(kd)*ang_par_F(:,:,kd)
            T3(:,:,kd,1) = T1(:,:,kd,1)*djq(kd) + &
                &J(:,:,kd)**2*vac_perm*eq_1%pres_FD(kd,1)/g33(:,:,kd)
            T4(:,:,kd,1) = T1(:,:,kd,1)*djq(kd)*ang_par_F(:,:,kd) - &
                &g23(:,:,kd)/g33(:,:,kd)
        end do
        T5(:,:,:,1) = T1(:,:,:,1)*D3Theta_3 - D1Theta_3
        T6(:,:,:,1) = T1(:,:,:,1)*Theta_3
        
        ! set up parallel derivatives of U in eq grid if VMEC is used
        if (eq_style.eq.1) then
            T1(:,:,:,2) = D3g13/g33 - g13*D3g33/g33**2
            do kd = 1,grid_eq%loc_n_r
                T2(:,:,kd,2) = D3h12(:,:,kd)/h22(:,:,kd) - &
                    &h12(:,:,kd)*D3h22(:,:,kd)/h22(:,:,kd)**2 + djq(kd)
                T3(:,:,kd,2) = T1(:,:,kd,2)*djq(kd) + &
                    &J(:,:,kd)*vac_perm*eq_1%pres_FD(kd,1)/g33(:,:,kd) * &
                    &(2*D3J(:,:,kd)-D3g33(:,:,kd)*J(:,:,kd)/g33(:,:,kd))
                T4(:,:,kd,2) = (T1(:,:,kd,2)*ang_par_F(:,:,kd)+T1(:,:,kd,1))*&
                    &djq(kd) - D3g23(:,:,kd)/g33(:,:,kd) + &
                    &g23(:,:,kd)*D3g33(:,:,kd)/g33(:,:,kd)**2
            end do
            T5(:,:,:,2) = T1(:,:,:,2)*D3Theta_3 + T1(:,:,:,1)*D33Theta_3 - &
                &D13Theta_3
            T6(:,:,:,2) = T1(:,:,:,2)*Theta_3 + T1(:,:,:,1)*D3Theta_3
        end if
        
        select case (X_grid_style)
            case (1)                                                            ! equilibrium
                ! point
                q_saf => eq_1%q_saf_FD(:,0)
                rot_t => eq_1%rot_t_FD(:,0)
                T1_X => T1
                T2_X => T2
                T3_X => T3
                T4_X => T4
                T5_X => T5
                T6_X => T6
            case (2)                                                            ! solution
                ! setup normal interpolation data
                ierr = setup_interp_data(grid_eq%loc_r_F,grid_X%loc_r_F,&
                    &norm_interp_data,norm_disc_prec_X)
                CHCKERR('')
                
                ! interpolate
                ierr = apply_disc(eq_1%q_saf_FD(:,0),norm_interp_data,q_saf)
                CHCKERR('')
                ierr = apply_disc(eq_1%rot_t_FD(:,0),norm_interp_data,rot_t)
                CHCKERR('')
                ierr = apply_disc(T1,norm_interp_data,T1_X,3)
                CHCKERR('')
                ierr = apply_disc(T2,norm_interp_data,T2_X,3)
                CHCKERR('')
                ierr = apply_disc(T3,norm_interp_data,T3_X,3)
                CHCKERR('')
                ierr = apply_disc(T4,norm_interp_data,T4_X,3)
                CHCKERR('')
                ierr = apply_disc(T5,norm_interp_data,T5_X,3)
                CHCKERR('')
                ierr = apply_disc(T6,norm_interp_data,T6_X,3)
                CHCKERR('')
                
                ! clean up
                call norm_interp_data%dealloc()
        end select
        
        ! set up n_frac
        do kd = 1,grid_X%loc_n_r
            if (use_pol_flux_F) then
                n_frac(kd,:) = X%n(kd,:)*q_saf(kd) - X%m(kd,:)
            else
                n_frac(kd,:) = -X%m(kd,:)*rot_t(kd) + X%n(kd,:)
            end if
        end do
        
        ! calculate U_0 and U_1 for all modes
        do ld = 1,X%n_mod
            ! calculate order 1 of U_0 and U_1
            if (U_style.ge.1) then
                X%U_0(:,:,:,ld) = -T2_X(:,:,:,1)
                do kd = 1,grid_X%loc_n_r
                    X%U_1(:,:,kd,ld) = iu/nm(kd,ld)
                end do
            end if
            ! include order 2
            if (U_style.ge.2) then
                do kd = 1,grid_X%loc_n_r
                    U_corr(:,:,kd,1) = T3_X(:,:,kd,1) + &
                        &iu*n_frac(kd,ld)*T4_X(:,:,kd,1)
                    U_corr(:,:,kd,2) = n_frac(kd,ld)/nm(kd,ld)*T1_X(:,:,kd,1)
                    X%U_0(:,:,kd,ld) = X%U_0(:,:,kd,ld) + iu/nm(kd,ld)*&
                        &U_corr(:,:,kd,1)
                    X%U_1(:,:,kd,ld) = X%U_1(:,:,kd,ld) + iu/nm(kd,ld)*&
                        &U_corr(:,:,kd,2)
                end do
            end if
            ! include order 3
            if (U_style.ge.3) then
                do kd = 1,grid_X%loc_n_r
                    U_corr(:,:,kd,1) = iu*n_frac(kd,ld)*&
                        &(-T5_X(:,:,kd,1)-iu*n_frac(kd,ld)*T6_X(:,:,kd,1))
                    X%U_0(:,:,kd,ld) = X%U_0(:,:,kd,ld) + (iu/nm(kd,ld))**2*&
                        &U_corr(:,:,kd,1)
                end do
            end if
            if (U_style.ge.4) then
                call writo('The geodesic perturbation U is implemented up to &
                    &order 3',warning=.true.)
            end if
        end do
        
        ! calculate DU_0 and DU_1, starting with first part
        select case (eq_style)
            case (1)                                                            ! VMEC
                do ld = 1,X%n_mod
                    ! calculate order 1 of DU_0 and DU_1
                    if (U_style.ge.1) then
                        X%DU_0(:,:,:,ld) = -T2_X(:,:,:,2)
                        X%DU_1(:,:,:,ld) = 0._dp
                    end if
                    ! include order 2
                    if (U_style.ge.2) then
                        do kd = 1,grid_X%loc_n_r
                            U_corr(:,:,kd,1) = T3_X(:,:,kd,2) + &
                                &iu*n_frac(kd,ld)*T4_X(:,:,kd,2)
                            U_corr(:,:,kd,2) = n_frac(kd,ld)/nm(kd,ld)*&
                                &T1_X(:,:,kd,2)
                            X%DU_0(:,:,kd,ld) = X%DU_0(:,:,kd,ld) + &
                                &iu/nm(kd,ld)*U_corr(:,:,kd,1)
                            X%DU_1(:,:,kd,ld) = X%DU_1(:,:,kd,ld) + &
                                &iu/nm(kd,ld)*U_corr(:,:,kd,2)
                        end do
                    end if
                    ! include order 3
                    if (U_style.ge.3) then
                        do kd = 1,grid_X%loc_n_r
                            U_corr(:,:,kd,1) = iu*n_frac(kd,ld)*&
                                &(-T5_X(:,:,kd,2)-iu*n_frac(kd,ld)*&
                                &T6_X(:,:,kd,2))
                            X%DU_0(:,:,kd,ld) = X%DU_0(:,:,kd,ld) + &
                                &(iu/nm(kd,ld))**2*U_corr(:,:,kd,1)
                        end do
                    end if
                    if (U_style.ge.4) then
                        call writo('The geodesic perturbation U is implemented &
                            &only up to order 3',warning=.true.)
                    end if
                end do
#if ldebug
                if (ltest) then
                    call writo('Test calculation of DU')
                    if(get_log(.false.,ind=.true.)) then
                        ierr = test_DU()
                        CHCKERR('')
                        call pause_prog(ind=.true.)
                    end if
                end if
#endif
            case (2)                                                            ! HELENA
                ! reset pointers
                if (use_pol_flux_F) then
                    ang_par_F => grid_X%theta_F
                else
                    ang_par_F => grid_X%zeta_F
                end if
                ! set up helper variable for derivative
                ! numerically derive U_0 and U_1
                do kd = 1,grid_X%loc_n_r
                    do jd = 1,grid_X%n(2)
                        ! set up parallel derivative data
                        ierr = setup_deriv_data(ang_par_F(:,jd,kd),&
                            &par_deriv_data,1,norm_disc_prec_X+1)               ! higher precision because previous derivative
                        CHCKERR('')
                        
                        ! calculate DX%U_0
                        ierr = apply_disc(X%U_0(:,jd,kd,:),par_deriv_data,&
                            &X%DU_0(:,jd,kd,:),1)
                        CHCKERR('')
                        
                        ! calculate DX%U_1
                        ierr = apply_disc(X%U_1(:,jd,kd,:),par_deriv_data,&
                            &X%DU_1(:,jd,kd,:),1)
                        CHCKERR('')
                    end do
                end do
                call par_deriv_data%dealloc()
        end select
        
        ! add second part of derivatives ~ n_frac
        do kd = 1,grid_X%loc_n_r
            do ld = 1,X%n_mod
                X%DU_0(:,:,kd,ld) = X%DU_0(:,:,kd,ld) + &
                    &iu*n_frac(kd,ld)*X%U_0(:,:,kd,ld)
                X%DU_1(:,:,kd,ld) = X%DU_1(:,:,kd,ld) + &
                    &iu*n_frac(kd,ld)*X%U_1(:,:,kd,ld)
            end do
        end do
        
        ! clean up
        nullify(ang_par_F,ang_geo_F)
        nullify(J,D3J)
        nullify(g13,D3g13)
        nullify(g23,D3g23)
        nullify(g33,D3g33)
        nullify(h12,D3h12)
        nullify(h22,D1h22,D3h22,D13h22,D33h22)
        nullify(h23,D1h23,D3h23,D13h23,D33h23)
        select case (X_grid_style) 
            case (1)                                                            ! equilibrium
                ! do nothing
            case (2)                                                            ! solution
                deallocate(q_saf,rot_t,T1_X,T2_X,T3_X,T4_X,T5_X,T6_X)
        end select
        nullify(q_saf,rot_t,T1_X,T2_X,T3_X,T4_X,T5_X,T6_X)
#if ldebug
    contains
        ! Test calculation of DU by deriving U numerically.
        integer function test_DU() result(ierr)
            use num_vars, only: norm_disc_prec_X
            
            character(*), parameter :: rout_name = 'test_DU'
            
            ! local variables
            integer :: jd, kd, ld                                               ! counters
            complex(dp), allocatable :: DU_0(:,:,:)                             ! alternative calculation for DU_0
            complex(dp), allocatable :: DU_1(:,:,:)                             ! alternative calculation for DU_1
            character(len=max_str_ln) :: file_name                              ! name of plot file
            character(len=max_str_ln) :: description                            ! description of plot
            
            ! initialize ierr
            ierr = 0
            
            ! warning
            call writo('This test is representable only if there are enough &
                &parallel points in the grid!',warning=.true.)
            
            ! output
            call writo('Going to test whether DU is consistent with U')
            call lvl_ud(1)
            
            ! set up DU_0 and DU_1
            allocate(DU_0(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
            allocate(DU_1(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
            
            ! reset pointers
            if (use_pol_flux_F) then
                ang_par_F => grid_X%theta_F
            else
                ang_par_F => grid_X%zeta_F
            end if
            
            ! loop over all modes
            do ld = 1,X%n_mod
                ! loop over all normal points
                do kd = 1,grid_X%loc_n_r
                    ! derive numerically
                    do jd = 1,grid_X%n(2)
                        ierr = setup_deriv_data(ang_par_F(:,jd,kd),&
                            &par_deriv_data,1,norm_disc_prec_X+1)
                        CHCKERR('')
                        ierr = apply_disc(X%U_0(:,jd,kd,ld),par_deriv_data,&
                            &DU_0(:,jd,kd))
                        CHCKERR('')
                        ierr = apply_disc(X%U_1(:,jd,kd,ld),par_deriv_data,&
                            &DU_1(:,jd,kd))
                        CHCKERR('')
                    end do
                end do
                
                ! set some variables
                file_name = 'TEST_RE_DU_0_'//trim(i2str(ld))
                description = 'Verifying DU_0 by deriving U_0 numerically for &
                    &mode '//trim(i2str(ld)//', real part')
                
                ! plot difference for RE DU_0
                call plot_diff_HDF5(rp(DU_0),rp(X%DU_0(:,:,:,ld)),&
                    &file_name,descr=description,output_message=.true.)
                
                ! set some variables
                file_name = 'TEST_IM_DU_0_'//trim(i2str(ld))
                description = 'Verifying DU_0 by deriving U_0 numerically for &
                    &mode '//trim(i2str(ld)//', imaginary part')
                
                ! plot difference for IM DU_0
                call plot_diff_HDF5(ip(DU_0),ip(X%DU_0(:,:,:,ld)),&
                    &file_name,descr=description,output_message=.true.)
                
                ! set some variables
                file_name = 'TEST_RE_DU_1_'//trim(i2str(ld))
                description = 'Verifying DU_1 by deriving U_1 numerically for &
                    &mode '//trim(i2str(ld)//', real part')
                
                ! plot difference for RE DU_1
                call plot_diff_HDF5(rp(DU_1),rp(X%DU_1(:,:,:,ld)),&
                    &file_name,descr=description,output_message=.true.)
                
                ! set some variables
                file_name = 'TEST_IM_DU_1_'//trim(i2str(ld))
                description = 'Verifying DU_1 by deriving U_1 numerically for &
                    &mode '//trim(i2str(ld)//', imaginary part')
                
                ! plot difference for IM DU_1
                call plot_diff_HDF5(ip(DU_1),ip(X%DU_1(:,:,:,ld)),&
                    &file_name,descr=description,output_message=.true.)
            end do
            call par_deriv_data%dealloc()
            
            ! user output
            call lvl_ud(-1)
            call writo('Test complete')
        end function test_DU
#endif
    end function calc_U
    
    !> calculate      \f$\widetilde{PV}_{k,m}^i\f$      (pol.      flux)      or
    !! \f$\widetilde{PV}_{l,n}^i\f$ (tor.  flux) at  all equilibrium  \c loc_n_r
    !! values.
    !!
    !! Like in calc_U(), use  is made of variables \c Ti that are  set up in the
    !! equilibrium  grid  and  are  afterwards  converted  to  \c  Ti_X  in  the
    !! perturbation grid,  which needs to  have the same angular  coordinates as
    !! the equilibrium grid:
    !!  \f[\begin{aligned}
    !!      T_1 &= \frac{h^{22}}{\mu_0} g_{33} \\
    !!      T_2 &= \mathcal{J}S +
    !!           \mu_0 \sigma \frac{g_{33}}{\mathcal{J}} h^{22} \\
    !!      T_3 &= \frac{\sigma}{\mathcal{J}} T_2 \\
    !!      T_4 &= \frac{1}{\mu_0} \mathcal{J}^2 h^{22} \\
    !!      T_5 &= 2 p' \kappa_n \ ,
    !!  \end{aligned}\f]
    !! with \f$T_i\f$ = \c Ti.
    !!
    !! The     interpolated    Ti_X     are    then     used    to     calculate
    !! \f$\widetilde{PV}_{k,m}^i\f$:
    !!  \f[\begin{aligned}
    !!      \widetilde{PV}_{k,m}^0 &= T_1(DU_k^{0*} - T_2)(DU_m^0 - T_2) -
    !!          T_3 + (nq-m)(nq-k)T_4 - T_5 \\
    !!      \widetilde{PV}_{k,m}^1 &= T_1(DU_k^{0*} - T_2) DU_m^1 \\
    !!      \widetilde{PV}_{k,m}^2 &= T_1 DU_k^{1*} DU_m^1 \ ,
    !!  \end{aligned}\f]
    !! where 
    !!  \f[DU_k^i = i (nq-m) U_k^i + \frac{\partial U_k^i}{\partial\theta} . \f]
    !!
    !! This is  valid for poloidal Flux  coordinates and where \f$n\f$  is to be
    !! replaced by  \f$m\f$ and \f$(nq-m)\f$  by \f$(n-\iota m)\f$  for toroidal
    !! Flux coordinates.
    !!
    !! \see See \cite weyens2014theory .
    !!
    !! \return ierr
    integer function calc_PV(grid_eq,grid_X,eq_1,eq_2,X_a,X_b,X,lim_sec_X) &
        &result(ierr)
        use num_vars, only: use_pol_flux_F, norm_disc_prec_X, X_grid_style
        use eq_vars, only: vac_perm
        use num_utilities, only: c
        use X_utilities, only: is_necessary_X
        use X_vars, only: n_mod_X
        use grid_vars, only: disc_type
        use grid_utilities, only: setup_interp_data, apply_disc
        
        character(*), parameter :: rout_name = 'calc_PV'
        
        ! use input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(grid_type), intent(in) :: grid_X                                   !< perturbation grid
        type(eq_1_type), intent(in), target :: eq_1                             !< flux equilibrium
        type(eq_2_type), intent(in), target :: eq_2                             !< metric equilibrium
        type(X_1_type), intent(in) :: X_a                                       !< vectorial perturbation variables of dimension 2
        type(X_1_type), intent(in) :: X_b                                       !< vectorial perturbation variables of dimension 2
        type(X_2_type), intent(inout) :: X                                      !< tensorial perturbation variables
        integer, intent(in), optional :: lim_sec_X(2,2)                         !< limits of \c m_X (pol flux) or \c n_X (tor flux) for both dimensions
        
        ! local variables
        integer :: m, k, kd                                                     ! counters
        integer :: c_loc(2)                                                     ! local c for symmetric and asymmetric variables
        logical :: calc_this(2)                                                 ! whether this combination needs to be calculated
        type(disc_type) :: norm_interp_data                                     ! data for normal interpolation
        
        ! jacobian
        real(dp), pointer :: J(:,:,:)                                           ! jac
        ! lower metric factors
        real(dp), pointer :: g33(:,:,:)                                         ! h_theta,theta or h_zeta,zeta
        ! upper metric factors
        real(dp), pointer :: h22(:,:,:)                                         ! h^psi,psi
        ! helper variables
        real(dp), pointer :: q_saf(:), rot_t(:)                                 ! safety factor and rotational transform in X grid
        real(dp), allocatable :: fac_n(:), fac_m(:)                             ! multiplicative factors for n and m
        ! PV factors
        real(dp), allocatable, target :: T1(:,:,:)                              ! h22/mu_0 g33
        real(dp), allocatable, target :: T2(:,:,:)                              ! JS + mu_0 sigma g33/J h22
        real(dp), allocatable, target :: T3(:,:,:)                              ! sigma/J T2
        real(dp), allocatable, target :: T4(:,:,:)                              ! 1/mu_0 J^2 h22
        real(dp), allocatable, target :: T5(:,:,:)                              ! 2 p' kappa_n
        real(dp), pointer :: T1_X(:,:,:)                                        ! T1 in X grid and par. deriv.
        real(dp), pointer :: T2_X(:,:,:)                                        ! T2 in X grid and par. deriv.
        real(dp), pointer :: T3_X(:,:,:)                                        ! T3 in X grid and par. deriv.
        real(dp), pointer :: T4_X(:,:,:)                                        ! T4 in X grid and par. deriv.
        real(dp), pointer :: T5_X(:,:,:)                                        ! T5 in X grid and par. deriv.
        
        ! initialize ierr
        ierr = 0
        
        ! allocate variables
        ! helper variables
        allocate(q_saf(grid_X%loc_n_r))
        allocate(rot_t(grid_X%loc_n_r))
        allocate(fac_n(grid_X%loc_n_r),fac_m(grid_X%loc_n_r))
        ! PV factors
        allocate(T1(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        allocate(T2(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        allocate(T3(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        allocate(T4(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        allocate(T5(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        select case (X_grid_style)
            case (1)                                                            ! equilibrium
                ! do nothing
            case (2)                                                            ! solution
                allocate(q_saf(grid_X%loc_n_r))
                allocate(rot_t(grid_X%loc_n_r))
                allocate(T1_X(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
                allocate(T2_X(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
                allocate(T3_X(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
                allocate(T4_X(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
                allocate(T5_X(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
        end select
        
        ! set pointers
        J => eq_2%jac_FD(:,:,:,0,0,0)
        g33 => eq_2%g_FD(:,:,:,c([3,3],.true.),0,0,0)
        h22 => eq_2%h_FD(:,:,:,c([2,2],.true.),0,0,0)
        
        ! set up PV factors in eq grid
        T1 = h22/(vac_perm*g33)
        T2 = J*eq_2%S + vac_perm*eq_2%sigma*g33/(J*h22)
        T3 = eq_2%sigma/J*T2
        T4 = 1._dp/(vac_perm*J**2*h22)
        do kd = 1,grid_eq%loc_n_r
            T5(:,:,kd) = 2*eq_1%pres_FD(kd,1)*eq_2%kappa_n(:,:,kd)
        end do
        
        select case (X_grid_style)
            case (1)                                                            ! equilibrium
                ! point
                q_saf => eq_1%q_saf_FD(:,0)
                rot_t => eq_1%rot_t_FD(:,0)
                T1_X => T1
                T2_X => T2
                T3_X => T3
                T4_X => T4
                T5_X => T5
            case (2)                                                            ! solution
                ! setup normal interpolation data
                ierr = setup_interp_data(grid_eq%loc_r_F,grid_X%loc_r_F,&
                    &norm_interp_data,norm_disc_prec_X)
                CHCKERR('')
                
                ! interpolate
                ierr = apply_disc(eq_1%q_saf_FD(:,0),norm_interp_data,q_saf)
                CHCKERR('')
                ierr = apply_disc(eq_1%rot_t_FD(:,0),norm_interp_data,rot_t)
                CHCKERR('')
                ierr = apply_disc(T1,norm_interp_data,T1_X,3)
                CHCKERR('')
                ierr = apply_disc(T2,norm_interp_data,T2_X,3)
                CHCKERR('')
                ierr = apply_disc(T3,norm_interp_data,T3_X,3)
                CHCKERR('')
                ierr = apply_disc(T4,norm_interp_data,T4_X,3)
                CHCKERR('')
                ierr = apply_disc(T5,norm_interp_data,T5_X,3)
                CHCKERR('')
                
                ! clean up
                call norm_interp_data%dealloc()
        end select
        
        ! set up fac_n and fac_m
        if (use_pol_flux_F) then
            fac_n = q_saf
            fac_m = 1.0_dp
        else
            fac_n = 1.0_dp
            fac_m = rot_t
        end if
        
        ! loop over all modes
        do m = 1,X_b%n_mod
            do k = 1,X_a%n_mod
                ! check whether modes are correct
                do kd = 1,grid_X%loc_n_r
                    if (X%n_1(kd,k).ne.X_a%n(kd,k) .or. &
                        &X%n_2(kd,m).ne.X_b%n(kd,m) .or. &
                        &X%m_1(kd,k).ne.X_a%m(kd,k) .or. &
                        &X%m_2(kd,m).ne.X_b%m(kd,m)) then
                        ierr = 1
                        CHCKERR('Modes do not match')
                    end if
                end do
                
                ! check whether mode combination needs to be calculated
                calc_this(1) = is_necessary_X(.true.,[k,m],lim_sec_X)
                calc_this(2) = is_necessary_X(.false.,[k,m],lim_sec_X)
                
                ! set up c_loc
                c_loc(1) = c([k,m],.true.,n_mod_X,lim_sec_X)
                c_loc(2) = c([k,m],.false.,n_mod_X,lim_sec_X)
                
                ! calculate PV_0
                if (calc_this(1)) then
                    do kd = 1,grid_X%loc_n_r
                        X%PV_0(:,:,kd,c_loc(1)) = T1_X(:,:,kd)*&
                            &(X_b%DU_0(:,:,kd,m) - T2_X(:,:,kd)) * &
                            &(conjg(X_a%DU_0(:,:,kd,k)) - T2_X(:,:,kd)) - &
                            &T3_X(:,:,kd) - T5_X(:,:,kd) + T4_X(:,:,kd)*&
                            &(X_b%n(kd,m)*fac_n(kd)-X_b%m(kd,m)*fac_m(kd))*&
                            &(X_a%n(kd,k)*fac_n(kd)-X_a%m(kd,k)*fac_m(kd))
                    end do
                end if
                
                ! calculate PV_1
                if (calc_this(2)) then
                    X%PV_1(:,:,:,c_loc(2)) = T1_X * X_b%DU_1(:,:,:,m) * &
                        &(conjg(X_a%DU_0(:,:,:,k))-T2_X)
                end if
                
                ! calculate PV_2
                if (calc_this(1)) then
                    X%PV_2(:,:,:,c_loc(1)) = T1_X * X_b%DU_1(:,:,:,m) * &
                        &conjg(X_a%DU_1(:,:,:,k))
                end if
            end do
        end do
        
        ! clean up
        nullify(J)
        nullify(g33)
        nullify(h22)
        select case (X_grid_style) 
            case (1)                                                            ! equilibrium
                ! do nothing
            case (2)                                                            ! solution
                deallocate(q_saf,rot_t,T1_X,T2_X,T3_X,T4_X,T5_X)
        end select
        nullify(q_saf,rot_t,T1_X,T2_X,T3_X,T4_X,T5_X)
    end function calc_PV
    
    !> calculate      \f$\widetilde{KV}_{k,m}^i\f$      (pol.      flux)      or
    !! \f$\widetilde{KV}_{l,n}^i\f$ (tor.  flux) at  all equilibrium  \c loc_n_r
    !! values.
    !!
    !! Like in calc_U(), use  is made of variables \c Ti that are  set up in the
    !! equilibrium  grid  and  are  afterwards  converted  to  \c  Ti_X  in  the
    !! perturbation grid,  which needs to  have the same angular  coordinates as
    !! the equilibrium grid:
    !!  \f[\begin{aligned}
    !!      T_1 &= \rho \mathcal{J}^2 \frac{h^{22}}{\mu_0} g_{33} \\
    !!      T_2 &= \rho \frac{1}{h^{22}} \ ,
    !!  \end{aligned}\f]
    !! with \f$T_i\f$ = \c Ti.
    !!
    !! The     interpolated    Ti_X     are    then     used    to     calculate
    !! \f$\widetilde{KV}_{k,m}^i\f$:
    !!  \f[\begin{aligned}
    !!      \widetilde{KV}_{k,m}^0 &= T_1 U_k^{0*} U_m^0 + T_2 \\
    !!      \widetilde{KV}_{k,m}^1 &= T_1 U_k^{0*} U_m^1 \\
    !!      \widetilde{KV}_{k,m}^2 &= T_1 U_k^{1*} U_m^1 \ ,
    !!  \end{aligned}\f]
    !! where 
    !!  \f[DU_k^i = i (nq-m) U_k^i + \frac{\partial U_k^i}{\partial\theta} . \f]
    !!
    !! This is  valid for poloidal Flux  coordinates and where \f$n\f$  is to be
    !! replaced by  \f$m\f$ and \f$(nq-m)\f$  by \f$(n-\iota m)\f$  for toroidal
    !! Flux coordinates.
    !!
    !! \see See \cite weyens2014theory .
    !!
    !! \return ierr
    integer function calc_KV(grid_eq,grid_X,eq_1,eq_2,X_a,X_b,X,lim_sec_X) &
        &result(ierr)
        use num_vars, only: K_style, norm_disc_prec_X, X_grid_style
        use num_utilities, only: c
        use X_utilities, only: is_necessary_X
        use X_vars, only: n_mod_X
        use grid_vars, only: disc_type
        use grid_utilities, only: setup_interp_data, apply_disc
        
        character(*), parameter :: rout_name = 'calc_KV'
        
        ! use input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(grid_type), intent(in) :: grid_X                                   !< perturbation grid
        type(eq_1_type), intent(in) :: eq_1                                     !< flux equilibrium
        type(eq_2_type), intent(in), target :: eq_2                             !< metric equilibrium
        type(X_1_type), intent(in) :: X_a                                       !< vectorial perturbation variables of dimension 2
        type(X_1_type), intent(in) :: X_b                                       !< vectorial perturbation variables of dimension 2
        type(X_2_type), intent(inout) :: X                                      !< tensorial perturbation variables
        integer, intent(in), optional :: lim_sec_X(2,2)                         !< limits of \c m_X (pol flux) or \c n_X (tor flux) for both dimensions
        
        ! local variables
        integer :: m, k, kd                                                     ! counters
        integer :: c_loc(2)                                                     ! local c for symmetric and asymmetric variables
        logical :: calc_this(2)                                                 ! whether this combination needs to be calculated
        type(disc_type) :: norm_interp_data                                     ! data for normal interpolation
        
        ! jacobian
        real(dp), pointer :: J(:,:,:)                                           ! jac
        ! lower metric factors
        real(dp), pointer :: g33(:,:,:)                                         ! h_theta,theta or h_zeta,zeta
        ! upper metric factors
        real(dp), pointer :: h22(:,:,:)                                         ! h^psi,psi
        ! KV factors
        real(dp), allocatable, target :: T1(:,:,:)                              ! h22/mu_0 J g33
        real(dp), allocatable, target :: T2(:,:,:)                              ! JS + mu_0 sigma g33/h22
        real(dp), pointer :: T1_X(:,:,:)                                        ! T1 in X grid and par. deriv.
        real(dp), pointer :: T2_X(:,:,:)                                        ! T2 in X grid and par. deriv.
        
        ! initialize ierr
        ierr = 0
        
        ! allocate variables
        ! KV factors
        allocate(T1(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        allocate(T2(grid_eq%n(1),grid_eq%n(2),grid_eq%loc_n_r))
        select case (X_grid_style)
            case (1)                                                            ! equilibrium
                ! do nothing
            case (2)                                                            ! solution
                allocate(T1_X(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
                allocate(T2_X(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
        end select
        
        ! set pointers
        J => eq_2%jac_FD(:,:,:,0,0,0)
        g33 => eq_2%g_FD(:,:,:,c([3,3],.true.),0,0,0)
        h22 => eq_2%h_FD(:,:,:,c([2,2],.true.),0,0,0)
        
        ! set up KV factors in eq grid
        do kd = 1,grid_eq%loc_n_r
            T1(:,:,kd) = eq_1%rho(kd)*J(:,:,kd)**2*h22(:,:,kd)/g33(:,:,kd)
            T2(:,:,kd) = eq_1%rho(kd)/h22(:,:,kd)
        end do
        
        select case (X_grid_style)
            case (1)                                                            ! equilibrium
                ! point
                T1_X => T1
                T2_X => T2
            case (2)                                                            ! solution
                ! setup normal interpolation data
                ierr = setup_interp_data(grid_eq%loc_r_F,grid_X%loc_r_F,&
                    &norm_interp_data,norm_disc_prec_X)
                CHCKERR('')
                
                ! interpolate
                ierr = apply_disc(T1,norm_interp_data,T1_X,3)
                CHCKERR('')
                ierr = apply_disc(T2,norm_interp_data,T2_X,3)
                CHCKERR('')
                
                ! clean up
                call norm_interp_data%dealloc()
        end select
        
        ! loop over all modes
        do m = 1,X_b%n_mod
            do k = 1,X_a%n_mod
                ! check whether modes are correct
                do kd = 1,grid_X%loc_n_r
                    if (X%n_1(kd,k).ne.X_a%n(kd,k) .or. &
                        &X%n_2(kd,m).ne.X_b%n(kd,m) .or. &
                        &X%m_1(kd,k).ne.X_a%m(kd,k) .or. &
                        &X%m_2(kd,m).ne.X_b%m(kd,m)) then
                        ierr = 1
                        CHCKERR('Modes do not match')
                    end if
                end do
                
                ! check whether mode combination needs to be calculated
                calc_this(1) = is_necessary_X(.true.,[k,m],lim_sec_X)
                calc_this(2) = is_necessary_X(.false.,[k,m],lim_sec_X)
                
                ! set up c_loc
                c_loc(1) = c([k,m],.true.,n_mod_X,lim_sec_X)
                c_loc(2) = c([k,m],.false.,n_mod_X,lim_sec_X)
                
                select case (K_style)
                    case (1)                                                    ! normalization of full perpendicular component
                        ! calculate KV_0
                        if (calc_this(1)) then
                            X%KV_0(:,:,:,c_loc(1)) = T2_X + T1_X * &
                                &X_b%U_0(:,:,:,m) * conjg(X_a%U_0(:,:,:,k))
                        end if
                        
                        ! calculate KV_1
                        if (calc_thiS(2)) then
                            X%KV_1(:,:,:,c_loc(2)) = T1_X * &
                                &X_b%U_1(:,:,:,m) * conjg(X_a%U_0(:,:,:,k))
                        end if
                        
                        ! calculate KV_2
                        if (calc_this(1)) then
                            X%KV_2(:,:,:,c_loc(1)) = T1_X * &
                                &X_b%U_1(:,:,:,m) * conjg(X_a%U_1(:,:,:,k))
                        end if
                    case (2)                                                    ! normalization of only normal component
                        ! calculate KV_0
                        if (calc_this(1)) then
                            X%KV_0(:,:,:,c_loc(1)) = T2_X
                        end if
                        
                        ! calculate KV_1
                        if (calc_this(2)) then
                            X%KV_1(:,:,:,c_loc(2)) = 0._dp
                        end if
                        
                        ! calculate KV_2
                        if (calc_this(1)) then
                            X%KV_2(:,:,:,c_loc(1)) = 0._dp
                        end if
                end select
            end do
        end do
        
        ! clean up
        nullify(J)
        nullify(g33)
        nullify(h22)
        select case (X_grid_style) 
            case (1)                                                            ! equilibrium
                ! do nothing
            case (2)                                                            ! solution
                deallocate(T1_X,T2_X)
        end select
        nullify(T1_X,T2_X)
    end function calc_KV
    
    !> Calculate  the magnetic  integrals from  \f$\widetilde{PV}_{k,m}^i\f$ and
    !! \f$\widetilde{KV}_{k,m}^i\f$ in  an equidistant grid where  the step size
    !! can vary depending on the normal coordinate.
    !! 
    !! All the variables should thus be field-line oriented. The result is saved
    !! in the first index of the X variables, the other can be ignored.
    !!
    !! Using \c prev_style, the results of this calculation on X can be combined
    !! with the ones already present in \c X_int.
    !!  - \c prev_style=0: Overwrite [def].
    !!  - \c prev_style=1: Add to current integral.
    !!  - \c prev_style=2: Divide by  2 and add to current integral.  Also modify the
    !!  indices of the current integral:
    !!      - for \c magn_int_style = 1: <tt>1 2 2 2 2 2 1 -> 2 2 2 2 2 2</tt>,
    !!      - for \c magn_int_style = 2: <tt>1 3 3 2 3 3 1 -> 3 2 3 3 2 3</tt>.
    !!  - \c prev_style=3:  Add to current integral. Also modify  the indices of
    !!  the current integral as in \c prev_style = 2.
    !!  - else: ignore previous magnetic integrals.
    !!
    !! Therefore, it is to be used in order. For example:
    !!  -# \f$\phantom{\rightarrow}\f$  (R=1,E=1): style 0,
    !!  -# \f$\rightarrow\f$            (R=1,E=2): style 1,
    !!  -# \f$\rightarrow\f$            (R=2,E=1): style 2,
    !!  -# \f$\rightarrow\f$            (R=2,E=2): style 3.
    !!
    !! Afterwards, it would cycle through 2 and 3,  as style 0 and 1 are only for
    !! the first Richardson level.
    !!
    !! \note
    !!  -# The variable type for the integrated tensorial perturbation variables
    !!  is also \c X_2, but it is  assumed that the first index has dimension 1,
    !!  the rest is ignored.
    !!  -# Simpson's 3/8  rule converges faster than the  trapezoidal rule, but
    !!  normally needs a better starting point (i.e. higher \c min_n_par_X)
    !!  -#  For alpha  style 1,  Weyl's  lemma is  used to  convert the  surface
    !!  integral  into a  line integral  along  the magnetic  field. This  means
    !!  that  the  resulting line  integral  still  needs  to be  multiplied  by
    !!  \f$\frac{2\pi}{M}\f$, where \f$M\f$ is the  number of times the poloidal
    !!  or  toroidal angle  completes a  full \f$2\pi\f$  tour, depending  on \c
    !!  use_pol_flux_F \cite Helander2014.
    !!  -# For alpha style 2, this factor is just the step size in alpha.
    !!
    !! \see See x_vars.x_2_type.
    !!
    !! \return ierr
    integer function calc_magn_ints(grid_eq,grid_X,eq,X,X_int,prev_style,&
        &lim_sec_X) result(ierr)
        use num_vars, only: use_pol_flux_F, norm_disc_prec_X, magn_int_style, &
            &alpha_style, X_grid_style
        use X_utilities, only: is_necessary_X
        use X_vars, only: n_mod_X
        use num_utilities, only: c
        use grid_vars, only: min_par_X, max_par_X, n_alpha, &
            &disc_type
        use grid_utilities, only: setup_interp_data, apply_disc
        use rich_vars, only: rich_lvl
        
        character(*), parameter :: rout_name = 'calc_magn_ints'
        
        ! input / output
        type(grid_type), intent(in) :: grid_eq                                  !< equilibrium grid
        type(grid_type), intent(in) :: grid_X                                   !< perturbation grid
        type(eq_2_type), intent(in), target :: eq                               !< metric equilibrium
        type(X_2_type), intent(in) :: X                                         !< tensorial perturbation variables
        type(X_2_type), intent(inout) :: X_int                                  !< interpolated tensorial perturbation variables, containing all mode combinations
        integer, intent(in), optional :: prev_style                             !< style to treat X_prev
        integer, intent(in), optional :: lim_sec_X(2,2)                         !< limits of \c m_X (pol flux) or \c n_X (tor flux) for both dimensions
        
        ! local variables, not used in subroutines
        logical :: calc_this(2)                                                 ! whether this combination needs to be calculated
        integer :: k, m                                                         ! counters
        integer :: kd                                                           ! counter
        integer :: c_loc(2)                                                     ! local c for symmetric and asymmetric variables
        integer :: c_tot(2)                                                     ! total c for symmetric and asymmetric variables
        integer :: prev_style_loc                                               ! local prev_style
        real(dp) :: alpha_fac                                                   ! factor for alpha (style 1: Weyl factor, style 2: step size)
        real(dp) :: prev_mult_fac                                               ! multiplicative factor for previous style
        real(dp), allocatable :: step_size(:,:)                                 ! step size for every field line and normal point
        real(dp), pointer :: ang_par_F(:,:,:) => null()                         ! parallel angle in flux coordinates
        complex(dp), allocatable :: V_int_work(:,:)                             ! work V_int
        type(disc_type) :: norm_interp_data                                     ! data for normal interpolation
        
        ! local variables, also used in subroutines
        integer :: nr_int_regions                                               ! nr. of integration regions
        integer, allocatable :: int_dims(:,:)                                   ! dimension information of integration regions
        real(dp), allocatable :: int_facs(:)                                    ! integration factors for each of the regions
        complex(dp), allocatable :: J_exp_ang(:,:,:)                            ! J * exponential of Flux parallel angle
        
        ! Makes use of nr_int_regions, int_dims, int_facs and J_exp_ang
        
        ! jacobian
        real(dp), pointer :: J(:,:,:)                                       ! jac
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Calculate field-line averages')
        call lvl_ud(1)
        
        ! set up local prev_style
        prev_style_loc = 0
        if (present(prev_style)) prev_style_loc = prev_style
        
        ! allocate variables
        allocate(J_exp_ang(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
        select case (X_grid_style)
            case (1)                                                            ! equilibrium
                ! point
                J => eq%jac_FD(:,:,:,0,0,0)
            case (2)                                                            ! solution
                allocate(J(grid_X%n(1),grid_X%n(2),grid_X%loc_n_r))
                
                ! setup normal interpolation data
                ierr = setup_interp_data(grid_eq%loc_r_F,grid_X%loc_r_F,&
                    &norm_interp_data,norm_disc_prec_X)
                CHCKERR('')
                
                ! interpolate
                ierr = apply_disc(eq%jac_FD(:,:,:,0,0,0),norm_interp_data,J,3)
                CHCKERR('')
                
                ! clean up
                call norm_interp_data%dealloc()
        end select
        
        ! set up parallel angle in flux coordinates
        if (use_pol_flux_F) then
            ang_par_F => grid_X%theta_F
        else
            ang_par_F => grid_X%zeta_F
        end if
        
        ! set up alpha factor
        select case (alpha_style)
            case (1)                                                            ! single field line, multiple turns
                alpha_fac = (2*pi)**2/((max_par_X-min_par_X)*pi)                ! Weyl factor
            case (2)                                                            ! multiple field lines, single turns
                alpha_fac = (2*pi)/n_alpha                                      ! grid size
        end select
        
        ! set up integration variables
        
        ! set nr_int_regions
        select case (magn_int_style)
            case (1)                                                            ! Trapezoidal rule
                select case (prev_style_loc)
                    case (2,3)                                                  ! change indices
                        nr_int_regions = 1
                    case default                                                ! don't change indices
                        nr_int_regions = 2
                end select
            case (2)                                                            ! Simpson's 3/8 rule
                select case (prev_style_loc)
                    case (2,3)                                                  ! change indices
                        nr_int_regions = 3
                    case default                                                ! don't change indices
                        nr_int_regions = 4
                end select
        end select
        
        ! set up integration factors and dimensions
        allocate(int_facs(nr_int_regions))
        allocate(int_dims(nr_int_regions,3))
        select case (magn_int_style)
            case (1)                                                            ! Trapezoidal rule
                select case (prev_style_loc)
                    case (2,3)                                                  ! change indices
                        int_facs = 0.5_dp*[2]
                        int_dims(1,:) = [1,grid_X%n(1),1]                       ! region 1: all points
                    case default                                                ! don't change indices
                        int_facs = 0.5_dp*[1,2]
                        int_dims(1,:) = [1,grid_X%n(1),grid_X%n(1)-1]           ! region 1: first and last points
                        int_dims(2,:) = [2,grid_X%n(1)-1,1]                     ! region 2: intermediate points
                end select
            case (2)                                                            ! Simpson's 3/8 rule
                select case (prev_style_loc)
                    case (2,3)                                                  ! change indices
                        int_facs = 0.375_dp*[3,2,3]
                        int_dims(1,:) = [1,grid_X%n(1)-2,3]                     ! region 1: intermediate points ~ 3
                        int_dims(2,:) = [2,grid_X%n(1)-1,3]                     ! region 2: intermediate points ~ 2
                        int_dims(3,:) = [3,grid_X%n(1),3]                       ! region 3: intermediate points ~ 3
                    case default                                                ! don't change indices
                        int_facs = 0.375_dp*[1,3,3,2]
                        int_dims(1,:) = [1,grid_X%n(1),grid_X%n(1)-1]           ! region 1: first and last points
                        int_dims(2,:) = [2,grid_X%n(1)-2,3]                     ! region 2: intermediate points ~ 3
                        int_dims(3,:) = [3,grid_X%n(1)-1,3]                     ! region 3: intermediate points ~ 3
                        int_dims(4,:) = [4,grid_X%n(1)-3,3]                     ! region 4: intermediate points ~ 2
                end select
        end select
        
        ! set up multiplicative factor for previous level
        select case (prev_style_loc)
            case (1,3)                                                          ! add to integral of previous equilibrium job
                prev_mult_fac = 1.0_dp
            case (2)                                                            ! add to integral of previous Richardson level / 2
                prev_mult_fac = 0.5_dp
            case default                                                        ! don't add, overwrite
                prev_mult_fac = 0.0_dp
        end select
        
        ! set up step size
        allocate(step_size(grid_X%n(2),grid_X%loc_n_r))
        step_size = ang_par_F(2,:,:)-ang_par_F(1,:,:)
        if (rich_lvl.gt.1) step_size = step_size*0.5_dp                         ! only half of points present: step size is actually half
        
        ! initialize local V_int
        allocate(V_int_work(grid_X%n(2),grid_X%loc_n_r))
        
        ! loop over all modes
        do m = 1,X%n_mod(2)
            do k = 1,X%n_mod(1)
                ! check whether mode combination needs to be calculated
                calc_this(1) = is_necessary_X(.true.,[k,m],lim_sec_X)
                calc_this(2) = is_necessary_X(.false.,[k,m],lim_sec_X)
                
                ! set up local and total indices
                c_loc(1) = c([k,m],.true.,n_mod_X,lim_sec_X)
                c_loc(2) = c([k,m],.false.,n_mod_X,lim_sec_X)
                c_tot(1) = c([lim_sec_X(1,1)-1+k,lim_sec_X(1,2)-1+m],.true.,&
                    &n_mod_X)
                c_tot(2) = c([lim_sec_X(1,1)-1+k,lim_sec_X(1,2)-1+m],.false.,&
                    &n_mod_X)
                
                ! calculate J_exp_ang
                do kd = 1,grid_X%loc_n_r
                    if (use_pol_flux_F) then
                        J_exp_ang(:,:,kd) = J(:,:,kd)*&
                            &exp(iu*(X%m_1(kd,k)-X%m_2(kd,m))*&
                            &ang_par_F(:,:,kd))
                    else
                        J_exp_ang(:,:,kd) = J(:,:,kd)*&
                            &exp(-iu*(X%n_1(kd,k)-X%n_2(kd,m))*&
                            &ang_par_F(:,:,kd))
                    end if
                end do
                
                ! integrate PV_0 and KV_0
                if (calc_this(1)) then
                    X_int%PV_0(1,:,:,c_tot(1)) = &
                        &prev_mult_fac*X_int%PV_0(1,:,:,c_tot(1))
                    X_int%KV_0(1,:,:,c_tot(1)) = &
                        &prev_mult_fac*X_int%KV_0(1,:,:,c_tot(1))
                    call calc_magn_int_loc(X%PV_0(:,:,:,c_loc(1))*alpha_fac,&
                        &X_int%PV_0(1,:,:,c_tot(1)),V_int_work,step_size)
                    call calc_magn_int_loc(X%KV_0(:,:,:,c_loc(1))*alpha_fac,&
                        &X_int%KV_0(1,:,:,c_tot(1)),V_int_work,step_size)
                end if
                
                ! integrate PV_1 and KV_1
                if (calc_this(2)) then
                    X_int%PV_1(1,:,:,c_tot(2)) = &
                        &prev_mult_fac*X_int%PV_1(1,:,:,c_tot(2))
                    X_int%KV_1(1,:,:,c_tot(2)) = &
                        &prev_mult_fac*X_int%KV_1(1,:,:,c_tot(2))
                    call calc_magn_int_loc(X%PV_1(:,:,:,c_loc(2))*alpha_fac,&
                        &X_int%PV_1(1,:,:,c_tot(2)),V_int_work,step_size)
                    call calc_magn_int_loc(X%KV_1(:,:,:,c_loc(2))*alpha_fac,&
                        &X_int%KV_1(1,:,:,c_tot(2)),V_int_work,step_size)
                end if
                
                ! integrate PV_2 and KV_2
                if (calc_this(1)) then
                    X_int%PV_2(1,:,:,c_tot(1)) = &
                        &prev_mult_fac*X_int%PV_2(1,:,:,c_tot(1))
                    X_int%KV_2(1,:,:,c_tot(1)) = &
                        &prev_mult_fac*X_int%KV_2(1,:,:,c_tot(1))
                    call calc_magn_int_loc(X%PV_2(:,:,:,c_loc(1))*alpha_fac,&
                        &X_int%PV_2(1,:,:,c_tot(1)),V_int_work,step_size)
                    call calc_magn_int_loc(X%KV_2(:,:,:,c_loc(1))*alpha_fac,&
                        &X_int%KV_2(1,:,:,c_tot(1)),V_int_work,step_size)
                end if
            end do
        end do
        
        ! clean up
        nullify(ang_par_F)
        select case (X_grid_style) 
            case (1)                                                            ! equilibrium
                ! do nothing
            case (2)                                                            ! solution
                deallocate(J)
        end select
        nullify(J)
        
        ! user output
        call lvl_ud(-1)
    contains
        ! Integrate local magnetic integral, adding to the previous result.
        !
        ! Makes use of nr_int_regions, int_dims, int_facs and J_exp_ang.
        !> \private
        subroutine calc_magn_int_loc(V,V_int,V_int_work,step_size)
            ! input / output
            complex(dp), intent(in) :: V(:,:,:)                                 ! V
            complex(dp), intent(inout) :: V_int(:,:)                            ! integrated V
            complex(dp), intent(inout) :: V_int_work(:,:)                       ! workint variable
            real(dp) :: step_size(:,:)                                          ! step size for every normal point
            
            ! local variables
            integer :: id, kd, ld                                               ! counters
            
            ! loop over regions
            do ld = 1,nr_int_regions
                ! initialize
                V_int_work = 0
                
                ! integrate
                do kd = 1,size(step_size,2)
                    do id = int_dims(ld,1),int_dims(ld,2),int_dims(ld,3)
                        V_int_work(:,kd) = V_int_work(:,kd) + &
                            &J_exp_ang(id,:,kd)*V(id,:,kd)*step_size(:,kd)
                    end do
                end do
                
                ! scale by integration factor and update V_int
                V_int = V_int + V_int_work*int_facs(ld)
            end do
        end subroutine calc_magn_int_loc
    end function calc_magn_ints
    
    !> Divides the perturbation jobs.
    !!
    !! This  concerns the  calculation of  the magnetic  integrals of  blocks of
    !! tensorial perturbation variables. These are  set up using the equilibrium
    !! and vectorial perturbation variables that  are stored in memory, but only
    !! the  integrated tensorial  result is  stored in  memory, with  negligible
    !! variables size.
    !!
    !! The size of the (k,m) pairs to  be calculated is determined by looking at
    !! what  fits in  memory  when the  equilibrium  and vectorial  perturbation
    !! variables are stored in memory first.
    !!
    !! \return ierr
    integer function divide_X_jobs(arr_size) result(ierr)
        use num_vars, only: max_X_mem, n_procs, X_jobs_lims
        use X_vars, only: n_mod_X
        use X_utilities, only: calc_memory_X
        
        character(*), parameter :: rout_name = 'divide_X_jobs'
        
        ! input / output
        integer, intent(in) :: arr_size                                         !< array size (using loc_n_r)
        
        ! local variables
        integer :: n_div                                                        ! factor by which to divide the total size
        real(dp) :: mem_size                                                    ! approximation of memory required for X variables
        integer :: n_mod_block                                                  ! nr. of modes in block
        integer, allocatable :: n_mod_loc(:)                                    ! number of modes per block
        character(len=max_str_ln) :: block_message                              ! message about how many blocks
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Dividing the perturbation jobs')
        call lvl_ud(1)
        
        ! calculate largest possible block of (k,m) values
        n_div = 0
        mem_size = huge(1._dp)
        do while (mem_size.gt.max_X_mem)
            n_div = n_div + 1
            n_mod_block = ceiling(n_mod_X*1._dp/n_div)
            ierr = calc_memory_X(2,arr_size,n_mod_block,mem_size)
            CHCKERR('')
            if (n_div.gt.n_mod_X) then
                ierr = 1
                err_msg = 'The memory limit is too low'
                CHCKERR(err_msg)
            end if
        end do
        if (n_div.gt.1) then
            block_message = 'The '//trim(i2str(n_mod_X))//&
                &' Fourier modes are split into '//trim(i2str(n_div))//&
                &' and '//trim(i2str(n_div**2))//' jobs are done separately'
            if (n_procs.lt.n_div**2) then
                block_message = trim(block_message)//', '//&
                    &trim(i2str(n_procs))//' at a time'
            else
                block_message = trim(block_message)//', but simultaneously'
            end if
        else
            block_message = 'The '//trim(i2str(n_mod_X))//' Fourier modes &
                &can be done without splitting them'
        end if
        call writo(block_message)
        call writo('The total memory for all processes together is estimated &
            &to be about '//trim(i2str(ceiling(mem_size)))//'MB')
        call lvl_ud(1)
        call writo('(maximum: '//trim(i2str(ceiling(max_X_mem)))//&
            &'MB left over from equilibrium variables)')
        call lvl_ud(-1)
        
        ! set up jobs data as illustrated below for 3 divisions, order 1:
        !   [1,2,3]
        ! or order 2:
        !   [1,4,7]
        !   [2,5,8]
        !   [3,6,9]
        ! etc.
        ! Also initialize the jobs taken to false.
        allocate(n_mod_loc(n_div))
        n_mod_loc = n_mod_X/n_div                                               ! number of radial points on this processor
        n_mod_loc(1:mod(n_mod_X,n_div)) = n_mod_loc(1:mod(n_mod_X,n_div)) + 1   ! add a mode to if there is a remainder
        call calc_X_jobs_lims(X_jobs_lims,n_mod_loc,2)
        
        ! user output
        call lvl_ud(-1)
        call writo('Perturbation jobs divided')
    contains
        ! Calculate X_jobs_lims.
        !
        ! These are  stored as follows: The  job index is the  second one, while
        ! the first index  ranges from 1..2. For  every job in the  range of job
        ! indices, thee first index shows the limits of
        ! the different orders sequentially. For example:
        !   - for order 1: [3,4],
        !       which corresponds to the 1-D range 3..4,
        !   - for order 2: [3,4,1,5],
        !       which corresponds to the 2-D range 3..4 x 1..5,
        !
        ! etc.
        !> \private
        recursive subroutine calc_X_jobs_lims(res,n_mod,ord)
            ! input / output
            integer, intent(inout), allocatable :: res(:,:)                     ! result
            integer, intent(in) :: n_mod(:)                                     ! X jobs data
            integer, intent(in) :: ord                                          ! order of data
            
            ! local variables
            integer, allocatable :: res_loc(:,:)                                ! local result
            integer :: n_div                                                    ! nr. of divisions of modes
            integer :: id                                                       ! counter
            
            ! set up nr. of divisions
            n_div = size(n_mod)
            
            ! (re)allocate result
            if (allocated(res)) deallocate(res)
            allocate(res(2*ord,n_div**ord))
            
            ! loop over divisions
            do id = 1,n_div
                if (ord.gt.1) then                                              ! call lower order
                    call calc_X_jobs_lims(res_loc,n_mod,ord-1)
                    res(1:2*ord-2,(id-1)*n_div**(ord-1)+1:id*n_div**(ord-1)) = &
                        &res_loc
                end if                                                          ! set for own order
                res(2*ord-1,(id-1)*n_div**(ord-1)+1:id*n_div**(ord-1)) = &
                    &sum(n_mod(1:id-1))+1
                res(2*ord,(id-1)*n_div**(ord-1)+1:id*n_div**(ord-1)) = &
                    &sum(n_mod(1:id))
            end do
        end subroutine calc_X_jobs_lims
    end function divide_X_jobs
end module
