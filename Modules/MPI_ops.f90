!------------------------------------------------------------------------------!
!   Operations related to MPI                                                  !
!------------------------------------------------------------------------------!
module MPI_ops
#include <PB3D_macros.h>
    use str_ops
    use messages
    use MPI
    use num_vars, only: dp, max_str_ln, pi
    
    implicit none
    private
    public start_MPI, stop_MPI, split_MPI, split_MPI_POST, abort_MPI, &
        &broadcast_input_vars, merge_MPI, get_next_job, divide_X_grid
    
contains
    ! start MPI and gather information
    ! [MPI] Collective call
    integer function start_MPI() result(ierr)
        use num_vars, only: glb_rank, glb_n_procs, grp_rank, n_groups, plt_rank
        
        character(*), parameter :: rout_name = 'start_MPI'
        
        ! initialize ierr
        ierr = 0
        
        ! start MPI
        call MPI_init(ierr)                                                     ! initialize MPI
        CHCKERR('MPI init failed')
        call MPI_Comm_rank(MPI_Comm_world,glb_rank,ierr)                        ! global rank
        CHCKERR('MPI rank failed')
        call MPI_Comm_size(MPI_Comm_world,glb_n_procs,ierr)                     ! nr. processes
        CHCKERR('MPI size failed')
        
        ! no alpha groups yet, so set group and plot rank to global rank
        grp_rank = glb_rank
        plt_rank = glb_rank
        
        ! initialize n_groups to 0
        n_groups = 0
    end function start_MPI
    
    ! determine   how   to   split    the   communicator   MPI_Comm_world   into
    ! subcommunicators,  according to  n_procs_per_alpha,  which determines  how
    ! many processors to use per field line
    ! [MPI] Collective call
    integer function split_MPI(n_r_eq,eq_limits) result(ierr)
        use X_vars, only: min_n_r_X
        use num_vars, only: n_procs_per_alpha, n_procs, n_alpha, &
            &MPI_Comm_groups, glb_rank, glb_n_procs, grp_rank, next_job, &
            &grp_n_procs, grp_nr, n_groups, next_job_win, plt_rank
        use files_ops, only: open_output
        use MPI_utilities, only: calc_n_groups
        !use num_vars, only: MPI_Comm_masters
        
        character(*), parameter :: rout_name = 'split_MPI'
        
        ! input / output
        integer, intent(in) :: n_r_eq                                           ! tot. nr. normal points in equilibrium grid
        integer, intent(inout) :: eq_limits(2)                                  ! min. and max. index of eq. grid for this process
        
        ! local variables
        integer :: remainder                                                    ! remainder of division
        integer :: id                                                           ! counter
        integer :: sum_groups                                                   ! sum of previous colros
        logical :: grp_found                                                    ! when searching for group of local process
        character(len=max_str_ln) :: str_1, str_2, str_3                        ! strings used in user messages
        !integer, allocatable :: glb_grp_master_rank(:)                          ! global ranks of the group masters
        !integer :: world_group                                                  ! MPI group associated to MPI_Comm_world
        !integer :: master_group                                                 ! MPI group associated to alpha group masters
        integer :: intsize                                                      ! size of MPI real
        integer(kind=MPI_ADDRESS_KIND) :: size_one = 1                          ! size equal to 1
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! determine how many  subcommunicators can be made with  n_procs. If the
        ! remainder of  the division is not  zero, assign an extra  processor to
        ! the first groups
        
        ! get n_groups
        call calc_n_groups(n_groups)
        
        ! limit n_procs_per_alpha to glb_n_procs
        if (n_procs_per_alpha.gt.glb_n_procs) n_procs_per_alpha = glb_n_procs   ! too few processors for even one group
        
        ! allocate n_procs for sorting into groups
        allocate(n_procs(n_groups))                                             ! allocate n_procs
        n_procs = n_procs_per_alpha
        
        ! add an extra process if the remainder is not zero
        remainder = glb_n_procs - n_groups*n_procs_per_alpha
        id = 1
        do while (remainder.gt.0)
            n_procs(id) = n_procs(id)+1
            id = mod(id,n_groups)+1
            remainder = remainder-1
        end do
        
        ! user messages
        if (glb_n_procs.eq.1) then
            str_1 = '1 Process'
        else
            str_1 = trim(i2str(glb_n_procs))//' Processes'
        end if
        if (n_groups.eq.1) then
            str_2 = '1 group'
        else
            str_2 = trim(i2str(n_groups))//' groups'
        end if
        if (n_procs_per_alpha.eq.1) then
            str_3 = ''
        else
            str_3 = ' of at least '//trim(i2str(n_procs_per_alpha))//&
                &' processes'
        end if
        if (grp_nr.eq.0) call writo(trim(str_1)//' in total, for '//&
            &trim(str_2)//trim(str_3))                                          ! how many groups
        if (n_groups.eq.1) then
            str_1 = 'This group dynamically solves'
        else
            str_1 = 'These '//trim(i2str(n_groups))//' groups &
                &dynamically solve'
        end if
        if (n_alpha.eq.1) then
            str_2 = '1 field line'
        else
            str_2 = trim(i2str(n_alpha))//' field lines'
        end if
        if (glb_n_procs.eq.1) then
            str_3 = ' in serial'
        else
            str_3 = ' in parallel'
        end if
        call writo(trim(str_1)//' a total of '//trim(str_2)//trim(str_3))       ! how many field lines to solve
        
        ! determine  group nr. of  local process, group rank and  grp_n_procs
        grp_nr = 0                                                              ! start group nr at 0, in MPI spirit
        grp_found = .false.                                                     ! group not yet found
        sum_groups = 0                                                          ! start sum of ranks at 0
        do while (.not.grp_found)
            sum_groups = sum_groups + n_procs(grp_nr+1)                         ! upper bound of nr. processes of groups from 0 to grp_nr+1
            if (glb_rank+1.le.sum_groups) then                                  ! group found
                grp_found = .true.
                grp_n_procs = n_procs(grp_nr+1)
            else                                                                ! group not found: try next group
                grp_rank = grp_rank - n_procs(grp_nr+1)                         ! subtract n_procs in this group from grp_rank
                grp_nr = grp_nr + 1                                             ! increment group nr.
            end if
        end do
        
        ! split MPI_Comm_world according to n_procs
        call MPI_Comm_split(MPI_Comm_world,grp_nr,grp_rank,MPI_Comm_groups,&
            &ierr)
        CHCKERR('Failed to split in groups')
        
        ! set plot rank equal to grp_rank
        plt_rank = grp_rank
        
        ! increment n_r_X if lower than grp_n_procs
        if (grp_n_procs.gt.min_n_r_X) then
            call writo('WARNING: The starting nr. of normal points of the &
                &perturbation is increased from '//trim(i2str(min_n_r_X))//&
                &' to '//trim(i2str(grp_n_procs))//' because there are too &
                &many processes')
            min_n_r_X = grp_n_procs
        end if
        
        !! take subset of MPI_Comm_world  containing all group masters with their
        !! ranks according to the group rank
        !allocate(glb_grp_master_rank(0:n_groups-1))                             ! allocate glb_grp_master_rank
        !glb_grp_master_rank(0) = 0                                              ! master of first group is also global master
        !do id = 1,n_groups-1
            !glb_grp_master_rank(id) = glb_grp_master_rank(id-1) + n_procs(id)
        !end do
        !call MPI_Comm_group(MPI_Comm_world,world_group,ierr)                    ! get group of MPI_Comm_world
        !CHCKERR('Failed to get global group')
        !call MPI_group_incl(world_group,n_groups,glb_grp_master_rank,&
            !&master_group,ierr)                                                 ! take master subset
        !CHCKERR('Failed to create group of group masters')
        !call MPI_Comm_create(MPI_Comm_world,master_group,MPI_Comm_masters,ierr) ! create communicator for master subset
        !err_msg = 'Failed to create communicator to group of group masters'
        !CHCKERR(err_msg)
        
        ! set starting next_job to 1 on global master
        if (glb_rank.eq.0) then
            next_job = 1
        else
            next_job = 0
        end if
        
        ! create  a window  to the  variable next_job in  the global  master for
        ! asynchronous access
        call MPI_Type_extent(MPI_INTEGER,intsize,ierr) 
        err_msg = 'Couldn''t determine the extent of MPI_REAL'
        CHCKERR(err_msg)
        if (grp_rank.eq.0) then
            call MPI_Win_create(next_job,1*size_one,intsize,MPI_INFO_NULL,&
                &MPI_Comm_world,next_job_win,ierr)
        else
            call MPI_Win_create(0,0*size_one,intsize,MPI_INFO_NULL,&
                &MPI_Comm_world,next_job_win,ierr)
        end if
        CHCKERR('Couldn''t create window to next_job')
        
        ! calculate the r range for the equilibrium calculations
        ierr = calc_eq_r_range(eq_limits)
        CHCKERR('')
        
        ! set fence so that the global master holds next_job = 1 for sure
        call MPI_Win_fence(0,next_job_win,ierr)
        CHCKERR('Couldn''t set fence')
        
        ! open output files  for the groups that are not  the master group (i.e.
        ! the group containing the global master)
        if (grp_nr.ne.0) then
            ierr = open_output()
            CHCKERR('')
        end if
    contains
        ! calculate group limits for equilibrium grid for this rank in the alpha
        ! group so that  the equilibrium quantities are calculated  for only the
        ! normal range that is relevant for the perturbation calculations, which
        ! are to be anticipated.
        ! This is done as follows, assuming an equidistant perturbation grid (in
        ! the perturbation  normal variable):  if the  nr. of  perturbation grid
        ! points n_r_X goes  to infinity, the division under  the processes will
        ! approach the  ideal value of  1/grp_n_procs each. For lower  values of
        ! n_r_X, the first  processes can carry an  additional perturbation grid
        ! point if the  remainder of the division is not  zero. The situation is
        ! at its most asymetric for the lowest value of n_r_X.
        ! Therefore,  to be  able  to cover  all the  levels  in the  Richardson
        ! extrapolation  (when n_r_X  is  subsequently  increased), the  minimum
        ! of  the   equilibrium  range  of   each  processor  is  to   be  given
        ! by  the  grid   point  that  includes  the  lowest   of  the  possible
        ! perturbation  points,  at  n_r_X  going  to  infinity  (i.e.  r_min  +
        ! grp_rank*(r_max-r_min)/grp_n_procs) and the maximum of the equilibrium
        ! range is to  be given by the  grid point that includes  the highest of
        ! the possible perturbation  points, at n_r_X at its  lowest value (i.e.
        ! given by divide_X_grid) plus 1,  because the plotting routines need an
        ! assymetric ghost region.
        ! Furthermore, for the conversion between points on the continuous range
        ! (0..1) and  discrete grid points  in the equilbrium  grid (1..n_r_eq),
        ! the subroutine  con2dis is used.
        ! [MPI] Collective call
        integer function calc_eq_r_range(eq_limits) result(ierr)
            use num_vars, only: grp_n_procs, grp_rank, use_pol_flux_E, &
                &use_pol_flux_F, eq_style, tol_norm_r
            use utilities, only: con2dis, dis2con, calc_int, interp_fun, &
                &calc_deriv, round_with_tol
            use VMEC, only: flux_t_V, Dflux_t_V, rot_t_V
            use HELENA, only: flux_p_H, qs
            use X_vars, only: min_n_r_X, min_r_X, max_r_X
            
            character(*), parameter :: rout_name = 'calc_eq_r_range'
            
            ! input / output
            integer, intent(inout) :: eq_limits(2)                              ! min. and max. index of eq. grid for this process
            
            ! local variables
            real(dp), allocatable :: flux_F(:), flux_E(:)                       ! either pol. or tor. flux in F and E
            real(dp), allocatable :: Dflux_p_H(:)                               ! normal derivative of flux_p_H
            integer :: X_limits(2)                                              ! min. and max. index of X grid for this process
            real(dp) :: grp_min_r_eq_F_con                                      ! grp_min_r_eq in continuous F coords.
            real(dp) :: grp_min_r_eq_E_con                                      ! grp_min_r_eq in continuous E coords.
            real(dp) :: grp_min_r_eq_E_dis                                      ! grp_min_r_eq in discrete E coords., unrounded
            integer :: grp_max_r_eq_F_dis                                       ! grp_max_r_eq in discrete F coords.
            real(dp) :: grp_max_r_eq_F_con                                      ! grp_max_r_eq in continuous F coords.
            real(dp) :: grp_max_r_eq_E_con                                      ! grp_max_r_eq in continuous E coords.
            real(dp) :: grp_max_r_eq_E_dis                                      ! grp_max_r_eq in discrete E coords., unrounded
            
            ! initialize ierr
            ierr = 0
            
            ! set up flux_F and flux_E,  depending on which equilibrium style is
            ! being used:
            !   1:  VMEC
            !   2:  HELENA
            allocate(flux_F(n_r_eq),flux_E(n_r_eq))
            select case (eq_style)
                case (1)                                                        ! VMEC
                    ! set up F flux
                    if (use_pol_flux_F) then
                        ierr = calc_int(-rot_t_V*Dflux_t_V,&
                            &1.0_dp/(n_r_eq-1.0_dp),flux_F)                     ! conversion VMEC LH -> RH coord. system
                        CHCKERR('')
                    else
                        flux_F = flux_t_V
                    end if
                    ! set up E flux
                    if (use_pol_flux_E) then
                    write(*,*) '!!! CHANGED rot_t_V to -rot_t_V BUT NOT SURE !!'
                        ierr = calc_int(-rot_t_V*Dflux_t_V,&
                            &1.0_dp/(n_r_eq-1.0_dp),flux_E)                     ! conversion VMEC LH -> RH coord. system
                        CHCKERR('')
                    else
                        flux_E = flux_t_V
                    end if
                case (2)                                                        ! HELENA
                    ! calculate normal derivative of flux_p_H
                    allocate(Dflux_p_H(n_r_eq))
                    ierr = calc_deriv(flux_p_H,Dflux_p_H,flux_p_H,1,1)
                    CHCKERR('')
                    ! set up F flux
                    if (use_pol_flux_F) then
                        flux_F = flux_p_H
                    else
                        ierr = calc_int(qs*Dflux_p_H,flux_p_H,flux_F)
                        CHCKERR('')
                    end if
                    ! set up E flux
                    if (use_pol_flux_E) then
                        flux_E = flux_p_H
                    else
                        ierr = calc_int(qs*Dflux_p_H,flux_p_H,flux_E)
                        CHCKERR('')
                    end if
                case default
                    err_msg = 'No equilibrium style associated with '//&
                        &trim(i2str(eq_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            ! normalize flux_F and flux_E to (0..1)
            ! Note: max_flux  is not yet  determined, so use  flux_F(n_r_eq) and
            ! flux_E(n_r_eq), knowing that this is the full global range
            flux_F = flux_F/flux_F(n_r_eq)
            flux_E = flux_E/flux_E(n_r_eq)
            
            ! use min_r_X and max_r_X, with grp_n_procs to get the minimum bound
            ! grp_min_r_eq for this rank
            ! 1. continuous F coords (min_r_X..max_r_X)
            grp_min_r_eq_F_con = min_r_X + &
                &grp_rank*(max_r_X-min_r_X)/grp_n_procs
            ! 2. include tolerance for lowermost bound
            if (grp_rank.eq.0) grp_min_r_eq_F_con = &
                &max(grp_min_r_eq_F_con-tol_norm_r,0._dp)                       ! include tolerance for lowermost bound
            ! 3 round with tolerance
            ierr = round_with_tol(grp_min_r_eq_F_con,0._dp,1._dp)
            CHCKERR('')
            ! 4. continuous E coords
            ierr = interp_fun(grp_min_r_eq_E_con,flux_E,grp_min_r_eq_F_con,&
                &flux_F)
            CHCKERR('')
            ! 5. discrete E index, unrounded
            ierr = con2dis(grp_min_r_eq_E_con,grp_min_r_eq_E_dis,flux_E)
            CHCKERR('')
            ! 6. discrete E index, rounded down
            eq_limits(1) = floor(grp_min_r_eq_E_dis)
            
            ! use min_r_X and max_r_X to calculate grp_max_r_eq
            ! 1. divide_X_grid for min_n_r_X normal points
            ierr = divide_X_grid(min_n_r_X,X_limits)                            ! divide the grid for the min_n_r_X tot. normal points
            CHCKERR('')
            ! 2. discrete F coords (1..min_n_r_X)
            grp_max_r_eq_F_dis = X_limits(2)
            ! 3. continous F coords
            ierr = dis2con(grp_max_r_eq_F_dis,grp_max_r_eq_F_con,[1,min_n_r_X],&
                &[min_r_X,max_r_X])                                             ! the perturbation grid is equidistant
            CHCKERR('')
            ! 4. include normal tolerance for uppermost bound
            if (grp_rank.eq.grp_n_procs-1) grp_max_r_eq_F_con = &
                min(grp_max_r_eq_F_con+tol_norm_r,1._dp)
            ! 5. round with tolerance
            ierr = round_with_tol(grp_max_r_eq_F_con,0._dp,1._dp)
            CHCKERR('')
            ! 6. continuous E coords
            ierr = interp_fun(grp_max_r_eq_E_con,flux_E,&
                &grp_max_r_eq_F_con,flux_F)
            CHCKERR('')
            ! 7. discrete E index, unrounded
            ierr = con2dis(grp_max_r_eq_E_con,grp_max_r_eq_E_dis,flux_E)
            CHCKERR('')
            ! 8. discrete E index, rounded up
            eq_limits(2) = ceiling(grp_max_r_eq_E_dis)
            
            deallocate(flux_F,flux_E)
        end function calc_eq_r_range
    end function split_MPI
    
    ! Divides a perturbation grid over all  the processes and use this to divide
    ! the equilibrium grid  as well. From now  on, there is supposed  to be only
    ! one group so  group nr. is equal  to global nr. This makes  using the PB3D
    ! routines easier.
    integer function split_MPI_POST(r_F_eq,r_F_X,i_lim_eq,i_lim_X) result(ierr)
        use num_vars, only: grp_nr, n_groups, grp_n_procs, glb_n_procs, &
            &grp_rank, ghost_width_POST
        use MPI_utilities, only: get_ser_var
        
        character(*), parameter :: rout_name = 'split_MPI_POST'
        
        ! input / output
        real(dp), intent(in) :: r_F_eq(:), r_F_X(:)                             ! equilibrium and perturbation r_F
        integer, intent(inout) :: i_lim_eq(2), i_lim_X(2)                       ! i_lim of equilibrium and perturbation
        
        ! local variables
        integer :: n_r_eq, n_r_X                                                ! total nr. of normal points in eq and X grid
        integer, allocatable :: grp_n_r_eq(:), grp_n_r_X(:)                     ! group nr. of normal points in eq and X grid
        integer :: remainder                                                    ! remainder of division
        integer :: id                                                           ! counter
        real(dp) :: min_X, max_X                                                ! min. and max. of r_F_X in range of this process
        real(dp), parameter :: tol = 1.E-6                                      ! tolerance for grids
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! set some variables
        n_groups = 1
        grp_nr = 1
        grp_n_procs = glb_n_procs
        n_r_eq = size(r_F_eq)
        n_r_X = size(r_F_X)
        allocate(grp_n_r_eq(glb_n_procs),grp_n_r_X(glb_n_procs))
        
        ! divide the perturbation grid equally over all the processes
        grp_n_r_X = n_r_X/glb_n_procs
        
        ! add an extra process if the remainder is not zero
        remainder = mod(n_r_X,glb_n_procs)
        id = 1
        do while (remainder.gt.0)
            grp_n_r_X(id) = grp_n_r_X(id)+1
            id = mod(id,glb_n_procs)+1
            remainder = remainder - 1
        end do
        
        ! set i_lim_X
        i_lim_X = [sum(grp_n_r_X(1:grp_rank))+1,sum(grp_n_r_X(1:grp_rank+1))]
        if (grp_rank.gt.0) i_lim_X(1) = i_lim_X(1)-ghost_width_POST             ! ghost region for num. deriv.
        if (grp_rank.lt.glb_n_procs-1) i_lim_X(2) = &
            &i_lim_X(2)+ghost_width_POST+1                                      ! ghost region for num. deriv. and overlap for int_vol
        min_X = minval(r_F_X(i_lim_X(1):i_lim_X(2)))
        max_X = maxval(r_F_X(i_lim_X(1):i_lim_X(2)))
        
        ! determine i_lim_eq: smallest eq range comprising X range
        i_lim_eq = [0,size(r_F_eq)+1]
        if (r_F_eq(1).lt.r_F_eq(size(r_F_eq))) then                             ! ascending r_F_eq
            do id = 1,size(r_F_eq)
                if (r_F_eq(id).le.min_X+tol) i_lim_eq(1) = id                   ! move lower limit up
                if (r_F_eq(size(r_F_eq)-id+1).ge.max_X-tol) &
                    &i_lim_eq(2) = size(r_F_eq)-id+1                            ! move upper limit down
            end do
        else                                                                    ! descending r_F_eq
            do id = 1,size(r_F_eq)
                if (r_F_eq(size(r_F_eq)-id+1).le.min_X+tol) &
                    &i_lim_eq(1) = size(r_F_eq)-id+1                            ! move lower limit up
                if (r_F_eq(id).ge.max_X-tol) i_lim_eq(2) = id                   ! move upper limit down
            end do
        end if
        
        ! check if valid limits found
        if (i_lim_eq(1).lt.1 .or. i_lim_eq(2).gt.size(r_F_eq)) then
            ierr = 1
            err_msg = 'Perturbation range not contained in equilibrium range'
            CHCKERR(err_msg)
        end if
    end function split_MPI_POST
    
    ! merge the MPI groups back to MPI_Comm_world
    ! [MPI] Collective call
    integer function merge_MPI() result(ierr)
        use num_vars, only: MPI_Comm_groups, grp_rank, next_job_win, glb_rank, &
            &plt_rank
        !use num_vars, only: n_groups, MPI_Comm_masters
        
        character(*), parameter :: rout_name = 'merge_MPI'
        
        ! local variables
        !character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! free window for next_job
        call MPI_Win_free(next_job_win,ierr)
        CHCKERR('Unable to free window for next job')
        
        ! free groups and masters communicators
        call MPI_Comm_free(MPI_Comm_groups,ierr)
        CHCKERR('Unable to free communicator for groups')
        !if (n_groups.gt.1) then
            !if (grp_rank.eq.0) then
                !call MPI_Comm_free(MPI_Comm_masters,ierr)
                !err_msg = 'Unable to free communicator for masters'
                !CHCKERR(err_msg)
            !end if
        !end if
        
        ! reset meaningless grp_rank and plt_rank to glb_rank
        grp_rank = glb_rank
        plt_rank = glb_rank
        
        ! barrier
        call MPI_Barrier(MPI_Comm_world,ierr)
        CHCKERR('Coulnd''t set barrier')
    end function merge_MPI
    
    ! stop MPI
    ! [MPI] Collective call
    integer function stop_MPI() result(ierr)
        character(*), parameter :: rout_name = 'stop_MPI'
        
        ! initialize ierr
        ierr = 0
        
        call writo('Stopping MPI')
        call MPI_finalize(ierr)
        CHCKERR('MPI stop failed')
    end function stop_MPI
    
    ! abort MPI
    ! [MPI] Collective call
    integer function abort_MPI() result(ierr)
        character(*), parameter :: rout_name = 'abort_MPI'
        
        ! initialize ierr
        ierr = 0
        
        call MPI_Abort(MPI_Comm_world,0,ierr)
        CHCKERR('MPI abort failed')
    end function abort_MPI
    
    ! compares  own next_job with next_job  of other group masters  and gets the
    ! maximum value. Returns this value
    ! [MPI] Collective call
    integer function get_next_job(next_job) result(ierr)
        use num_vars, only: grp_rank, next_job_win, MPI_Comm_groups, &
            &n_alpha, grp_nr
        
        character(*), parameter :: rout_name = 'get_next_job'
        
        ! input / output
        integer, intent(inout) :: next_job
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer(kind=MPI_ADDRESS_KIND) :: null_disp = 0                         ! zero target displacement
        
        ! initialize ierr
        ierr = 0
        
        ! test whether it is a master that calls this routine
        if (grp_rank.eq.0) then
            err_msg = 'Group '//trim(i2str(grp_nr+1))//&
                &' coulnd''t lock window on global master'
            call MPI_Win_lock(MPI_LOCK_EXCLUSIVE,0,0,next_job_win,ierr)
            CHCKERR(err_msg)
            err_msg = 'Group '//trim(i2str(grp_nr+1))//&
                &' coulnd''t increment next_job on global master'
            !call MPI_Get_accumulate(1,1,MPI_INTEGER,next_job,1,MPI_INTEGER,0,&
                !&null_disp,1,MPI_INTEGER,MPI_SUM,next_job_win,ierr)             ! ONLY WORKS FOR OPENMPI 1.8 or higher
            !CHCKERR(err_msg)
            call MPI_Get(next_job,1,MPI_INTEGER,0,null_disp,1,MPI_INTEGER,&
                &next_job_win,ierr)
            CHCKERR(err_msg)
            call MPI_accumulate(1,1,MPI_INTEGER,0,null_disp,1,MPI_INTEGER,&
                &MPI_SUM,next_job_win,ierr)
            CHCKERR(err_msg)
            err_msg = 'Group '//trim(i2str(grp_nr+1))//&
                &' coulnd''t unlock window on global master'
            call MPI_Win_unlock(0,next_job_win,ierr)
            CHCKERR(err_msg)
            
            ! if all jobs reached, output -1
            if (next_job.ge.n_alpha+1) then
                next_job = -1
            end if
        end if
        
        ! broadcast next_job to whole group
        call MPI_Bcast(next_job,1,MPI_INTEGER,0,MPI_Comm_groups,ierr)
        CHCKERR('MPI broadcast failed')
    end function get_next_job
    
    ! divides a  grid of  n_r_X points  under the  ranks of  MPI_Comm_groups and
    ! assigns grp_n_r_X and grp_min_r_X and  grp_max_r_X to each rank. Also sets
    ! up r_X, which contains the normal  variable in the perturbation grid (with
    ! global range (min_r_X..max_r_X) translated to max_flux/2pi)
    ! Note: for  the first ranks,  the upper index is  one higher than  might be
    ! expected because  the plotting  routines need  information about  the next
    ! perturbation point (so this is an asymetric ghost region)
    integer function divide_X_grid(n_r_X,X_limits,r_X) result(ierr)
        use num_vars, only: MPI_Comm_groups, use_pol_flux_F, grp_rank, &
            &grp_n_procs
        use X_vars, only: min_r_X, max_r_X
        use utilities, only: round_with_tol
        use eq_vars, only: max_flux_p_F, max_flux_t_F
        
        character(*), parameter :: rout_name = 'divide_X_grid'
        
        ! input / output
        integer, intent(in) :: n_r_X                                            ! tot. nr. of normal points in pert. grid
        integer, intent(inout) :: X_limits(2)                                   ! min. and max. index of X grid for this process
        real(dp), intent(inout), allocatable, optional :: r_X(:)                ! normal points in Flux coords., globally normalized to (min_r_X..max_r_X)
        
        ! local variables
        integer :: n_procs                                                      ! nr. of procs.
        integer :: rank                                                         ! rank
        integer :: id                                                           ! counter
        integer :: grp_n_r_X                                                    ! nr. of points in group normal X grid
        real(dp) :: max_flux_F                                                  ! either max_flux_p_F or max_flux_t_F
        
        ! initialize ierr
        ierr = 0
        
        ! set n_procs and rank
        call MPI_Comm_size(MPI_Comm_groups,n_procs,ierr)
        CHCKERR('Failed to get MPI size')
        call MPI_Comm_rank(MPI_Comm_groups,rank,ierr)
        CHCKERR('Failed to get MPI rank')
        
        ! calculate n_loc for this rank
        grp_n_r_X = divide_X_grid_ind(rank,n_r_X,n_procs)
        
        ! Note: ghost region is only needed for output
        ! add ghost region if not last process
        if (grp_rank.lt.grp_n_procs-1) grp_n_r_X = grp_n_r_X + 1
        
        ! calculate the starting index of this rank
        X_limits(1) = 1
        do id = 0,rank-1
            X_limits(1) = X_limits(1) + divide_X_grid_ind(id,n_r_X,n_procs)
        end do
        ! calculate the end index of this rank
        X_limits(2) = X_limits(1) - 1 + grp_n_r_X
        
        ! limit grp_n_r_X and X_limits(2) to total n_r_X
        X_limits(2) = min(X_limits(2),n_r_X)
        grp_n_r_X = min(grp_n_r_X,X_limits(2)-X_limits(1)+1)
        
        ! set up r_X if present (equidistant grid in Flux coordinates)
        if (present(r_X)) then
            ! set up max_flux
            if (use_pol_flux_F) then
                max_flux_F = max_flux_p_F
            else
                max_flux_F = max_flux_t_F
            end if
            
            ! allocate r_X
            if (allocated(r_X)) deallocate(r_X)
            allocate(r_X(n_r_X))
            
            ! calculate r_X in range from 0 to 1
            r_X = [(min_r_X + (id-1.0_dp)/(n_r_X-1.0_dp)*(max_r_X-min_r_X),&
                &id=1,n_r_X)]
            
            ! round with standard tolerance
            ierr = round_with_tol(r_X,0.0_dp,1.0_dp)
            CHCKERR('')
            
            ! translate to the real normal variable in range from 0..flux/2pi
            r_X = r_X*max_flux_F/(2*pi)
        end if
    contains 
        integer function divide_X_grid_ind(rank,n,n_procs) result(n_loc)
            ! input / output
            integer, intent(in) :: rank                                         ! rank for which to be calculate divide_X_grid_ind
            integer, intent(in) :: n                                            ! tot. nr. of points to be divided under n_procs
            integer, intent(in) :: n_procs                                      ! nr. of processes
            
            n_loc = n/n_procs                                                   ! number of radial points on this processor
            if (mod(n,n_procs).gt.0) then                                       ! check if there is a remainder
                if (mod(n,n_procs).gt.rank) n_loc = n_loc + 1                   ! add a point to the first ranks if there is a remainder
            end if
        end function divide_X_grid_ind
    end function divide_X_grid
    
    ! Broadcasts all  the relevant variable that have been  determined so far in
    ! the global master process using the inputs to the other processes
    ! [MPI] Collective call
    integer function broadcast_input_vars() result(ierr)
        use num_vars, only: max_str_ln, ltest, EV_style, max_it_NR, max_it_r, &
            &n_alpha, n_procs_per_alpha, minim_style, max_alpha, min_alpha, &
            &tol_NR, glb_rank, glb_n_procs, no_guess, n_sol_requested, &
            &nyq_fac, tol_r, use_pol_flux_F, use_pol_flux_E, retain_all_sol, &
            &plot_flux_q, plot_grid, no_plots, eq_style, use_normalization, &
            &n_sol_plotted, n_theta_plot, n_zeta_plot, plot_resonance, EV_BC, &
            &tol_slepc, rho_style, prog_style, max_it_inv, norm_disc_ord, &
            &BC_style, tol_norm_r, max_n_it_slepc
        use VMEC, only: mpol, ntor, lasym, lfreeb, nfp, rot_t_V, gam, R_V_c, &
            &R_V_s, Z_V_c, Z_V_s, L_V_c, L_V_s, flux_t_V, Dflux_t_V, pres_V
        use HELENA, only: pres_H, qs, flux_p_H, nchi, chi_H, ias, h_H_11, &
            &h_H_12, h_H_33, RBphi, R_H, Z_H
        use eq_vars, only: R_0, pres_0, B_0, psi_0, rho_0, T_0, vac_perm
        use X_vars, only: min_m_X, max_m_X, min_n_X, max_n_X, min_n_r_X, &
            &min_r_X, max_r_X
        use grid_vars, only: n_r_eq, n_par_X, min_par_X, max_par_X
        use HDF5_ops, only: var_1D
        use PB3D_vars, only: vars_1D_eq, vars_1D_eq_B, vars_1D_X
#if ldebug
        use VMEC, only: B_V_sub_c, B_V_sub_s, B_V_c, B_V_s, jac_V_c, jac_V_s
#endif
        
        character(*), parameter :: rout_name = 'broadcast_input_vars'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: id                                                           ! counter
        
        ! initialize ierr
        ierr = 0
        
        ! set err_msg
        err_msg = 'MPI broadcast failed'
        
        if (glb_n_procs.gt.1) then                                              ! need to broadcast from global rank 0 to other processes
            ! print message
            call writo('Broadcasting variables determined by input')
            call lvl_ud(1)
            
            ! broadcast eq_style to determine what else to broadcast
            call MPI_Bcast(eq_style,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            
            ! variables that are sent for every program style:
            call MPI_Bcast(use_normalization,1,MPI_LOGICAL,0,MPI_Comm_world,&
                &ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(min_m_X,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(max_m_X,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(min_n_X,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(max_n_X,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(min_r_X,1,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(max_r_X,1,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(vac_perm,1,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,&
                &ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(rho_style,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(plot_resonance,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(plot_flux_q,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            call MPI_Bcast(plot_grid,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
            CHCKERR(err_msg)
            
            ! select according to program style
            select case (prog_style)
                case(1)                                                         ! PB3D
                    call MPI_Bcast(ltest,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(use_pol_flux_F,1,MPI_LOGICAL,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(use_pol_flux_E,1,MPI_LOGICAL,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(no_guess,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(no_plots,1,MPI_LOGICAL,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(max_it_NR,1,MPI_INTEGER,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(max_it_r,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(max_it_inv,1,MPI_INTEGER,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(minim_style,1,MPI_INTEGER,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(norm_disc_ord,1,MPI_INTEGER,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(BC_style,2,MPI_INTEGER,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(n_par_X,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(n_r_eq,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(EV_style,1,MPI_INTEGER,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(n_procs_per_alpha,1,MPI_INTEGER,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(n_alpha,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(n_sol_requested,1,MPI_INTEGER,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(min_n_r_X,1,MPI_INTEGER,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(nyq_fac,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(retain_all_sol,1,MPI_INTEGER,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(n_theta_plot,1,MPI_INTEGER,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(n_zeta_plot,1,MPI_INTEGER,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(max_n_it_slepc,1,MPI_INTEGER,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(min_alpha,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(tol_norm_r,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(max_alpha,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(tol_NR,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(tol_r,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(min_par_X,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(max_par_X,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(gam,1,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(R_0,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(pres_0,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(B_0,1,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(psi_0,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(rho_0,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(T_0,1,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(EV_BC,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(tol_slepc,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                case(2)                                                         ! PB3D_POST
                    call MPI_Bcast(n_sol_plotted,4,MPI_INTEGER,0,&
                        &MPI_Comm_world,ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(n_theta_plot,1,MPI_INTEGER,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call MPI_Bcast(n_zeta_plot,1,MPI_INTEGER,0,MPI_Comm_world,&
                        &ierr)
                    CHCKERR(err_msg)
                    call bcast_size_1_var_1D(vars_1D_X)
                    CHCKERR(err_msg)
                    do id = 1,size(vars_1D_X)
                        call bcast_size_1_R(vars_1D_X(id)%p)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_X(id)%p,size(vars_1D_X(id)%p),&
                            &MPI_DOUBLE_PRECISION,&
                            &0,MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                        call bcast_size_1_I(vars_1D_X(id)%tot_i_min)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_X(id)%tot_i_min,&
                            &size(vars_1D_X(id)%tot_i_min),MPI_INTEGER,0,&
                            &MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                        call bcast_size_1_I(vars_1D_X(id)%tot_i_max)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_X(id)%tot_i_max,&
                            &size(vars_1D_X(id)%tot_i_max),MPI_INTEGER,0,&
                            &MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_X(id)%var_name,max_str_ln,&
                            &MPI_CHARACTER,0,MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                    end do
                    call bcast_size_1_var_1D(vars_1D_eq)
                    CHCKERR(err_msg)
                    do id = 1,size(vars_1D_eq)
                        call bcast_size_1_R(vars_1D_eq(id)%p)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_eq(id)%p,size(vars_1D_eq(id)%p),&
                            &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                        call bcast_size_1_I(vars_1D_eq(id)%tot_i_min)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_eq(id)%tot_i_min,&
                            &size(vars_1D_eq(id)%tot_i_min),MPI_INTEGER,0,&
                            &MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                        call bcast_size_1_I(vars_1D_eq(id)%tot_i_max)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_eq(id)%tot_i_max,&
                            &size(vars_1D_eq(id)%tot_i_max),MPI_INTEGER,0,&
                            &MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_eq(id)%var_name,max_str_ln,&
                            &MPI_CHARACTER,0,MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                    end do
                    call bcast_size_1_var_1D(vars_1D_eq_B)
                    CHCKERR(err_msg)
                    do id = 1,size(vars_1D_eq_B)
                        call bcast_size_1_R(vars_1D_eq_B(id)%p)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_eq_B(id)%p,&
                            &size(vars_1D_eq_B(id)%p),MPI_DOUBLE_PRECISION,0,&
                            &MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                        call bcast_size_1_I(vars_1D_eq_B(id)%tot_i_min)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_eq_B(id)%tot_i_min,&
                            &size(vars_1D_eq_B(id)%tot_i_min),MPI_INTEGER,0,&
                            &MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                        call bcast_size_1_I(vars_1D_eq_B(id)%tot_i_max)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_eq_B(id)%tot_i_max,&
                            &size(vars_1D_eq_B(id)%tot_i_max),MPI_INTEGER,0,&
                            &MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                        call MPI_Bcast(vars_1D_eq_B(id)%var_name,max_str_ln,&
                            &MPI_CHARACTER,0,MPI_Comm_world,ierr)
                        CHCKERR(err_msg)
                    end do
                case default
                    err_msg = 'No program style associated with '//&
                        &trim(i2str(prog_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            ! For  specific variables, choose  which equilibrium style  is being
            ! used:
            !   1:  VMEC
            !   2:  HELENA
            select case (eq_style)
                case (1)                                                        ! VMEC
                    ! select according to program style
                    select case (prog_style)
                        case(1)                                                 ! PB3D
                            call bcast_size_1_R(flux_t_V)
                            CHCKERR(err_msg)
                            call MPI_Bcast(flux_t_V,size(flux_t_V),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_1_R(Dflux_t_V)
                            CHCKERR(err_msg)
                            call MPI_Bcast(Dflux_t_V,size(Dflux_t_V),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_1_R(rot_t_V)
                            CHCKERR(err_msg)
                            call MPI_Bcast(rot_t_V,size(rot_t_V),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_1_R(pres_V)
                            CHCKERR(err_msg)
                            call MPI_Bcast(pres_V,size(pres_V),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call MPI_Bcast(lasym,1,MPI_LOGICAL,0,&
                                &MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call MPI_Bcast(lfreeb,1,MPI_LOGICAL,0,&
                                &MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call MPI_Bcast(mpol,1,MPI_INTEGER,0,&
                                &MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call MPI_Bcast(ntor,1,MPI_INTEGER,0,MPI_Comm_world,&
                                &ierr)
                            CHCKERR(err_msg)
                            call MPI_Bcast(nfp,1,MPI_INTEGER,0,MPI_Comm_world,&
                                    &ierr)
                            CHCKERR(err_msg)
                            call bcast_size_4_R(R_V_c,0)
                            CHCKERR(err_msg)
                            call MPI_Bcast(R_V_c,size(R_V_c),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_4_R(R_V_s,0)
                            CHCKERR(err_msg)
                            call MPI_Bcast(R_V_s,size(R_V_s),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_4_R(Z_V_c,0)
                            CHCKERR(err_msg)
                            call MPI_Bcast(Z_V_c,size(Z_V_c),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_4_R(Z_V_s,0)
                            CHCKERR(err_msg)
                            call MPI_Bcast(Z_V_s,size(Z_V_s),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_4_R(L_V_c,0)
                            CHCKERR(err_msg)
                            call MPI_Bcast(L_V_c,size(L_V_c),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_4_R(L_V_s,0)
                            CHCKERR(err_msg)
                            call MPI_Bcast(L_V_s,size(L_V_s),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
#if ldebug
                            if (ltest) then                                     ! ltest has already been broadcast
                                call bcast_size_4_R(B_V_sub_c,1)
                                CHCKERR(err_msg)
                                call MPI_Bcast(B_V_sub_c,size(B_V_sub_c),&
                                    &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                                CHCKERR(err_msg)
                                call bcast_size_4_R(B_V_sub_s,1)
                                CHCKERR(err_msg)
                                call MPI_Bcast(B_V_sub_s,size(B_V_sub_s),&
                                    &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                                CHCKERR(err_msg)
                                call bcast_size_3_R(B_V_c)
                                CHCKERR(err_msg)
                                call MPI_Bcast(B_V_c,size(B_V_c),&
                                    &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                                CHCKERR(err_msg)
                                call bcast_size_3_R(B_V_s)
                                CHCKERR(err_msg)
                                call MPI_Bcast(B_V_s,size(B_V_s),&
                                    &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                                CHCKERR(err_msg)
                                call bcast_size_3_R(jac_V_c)
                                CHCKERR(err_msg)
                                call MPI_Bcast(jac_V_c,size(jac_V_c),&
                                    &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                                CHCKERR(err_msg)
                                call bcast_size_3_R(jac_V_s)
                                CHCKERR(err_msg)
                                call MPI_Bcast(jac_V_s,size(jac_V_s),&
                                    &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                                CHCKERR(err_msg)
                            end if
#endif
                        case(2)                                                 ! PB3D_POST
                            ! do nothing extra
                        case default
                            err_msg = 'No program style associated with '//&
                                &trim(i2str(prog_style))
                            ierr = 1
                            CHCKERR(err_msg)
                    end select
                case (2)                                                        ! HELENA
                    ! select according to program style
                    select case (prog_style)
                        case(1)                                                 ! PB3D
                            call bcast_size_2_R(R_H)
                            CHCKERR(err_msg)
                            call MPI_Bcast(R_H,size(R_H),MPI_DOUBLE_PRECISION,&
                                &0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_2_R(Z_H)
                            CHCKERR(err_msg)
                            call MPI_Bcast(Z_H,size(Z_H),MPI_DOUBLE_PRECISION,&
                                &0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_1_R(chi_H)
                            CHCKERR(err_msg)
                            call MPI_Bcast(chi_H,size(chi_H),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_1_R(flux_p_H)
                            CHCKERR(err_msg)
                            call MPI_Bcast(flux_p_H,size(flux_p_H),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call MPI_Bcast(nchi,1,MPI_INTEGER,0,MPI_Comm_world,&
                                &ierr)
                            CHCKERR(err_msg)
                            call MPI_Bcast(ias,1,MPI_INTEGER,0,MPI_Comm_world,&
                                &ierr)
                            CHCKERR(err_msg)
                            call bcast_size_1_R(qs)
                            CHCKERR(err_msg)
                            call MPI_Bcast(qs,size(qs),MPI_DOUBLE_PRECISION,0,&
                                &MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_1_R(pres_H)
                            CHCKERR(err_msg)
                            call MPI_Bcast(pres_H,size(pres_H),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_1_R(RBphi)
                            CHCKERR(err_msg)
                            call MPI_Bcast(RBphi,size(RBphi),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_2_R(h_H_11)
                            CHCKERR(err_msg)
                            call MPI_Bcast(h_H_11,size(h_H_11),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_2_R(h_H_12)
                            CHCKERR(err_msg)
                            call MPI_Bcast(h_H_12,size(h_H_12),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                            call bcast_size_2_R(h_H_33)
                            CHCKERR(err_msg)
                            call MPI_Bcast(h_H_33,size(h_H_33),&
                                &MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
                            CHCKERR(err_msg)
                        case(2)                                                 ! PB3D_POST
                            ! do nothing extra
                        case default
                            err_msg = 'No program style associated with '//&
                                &trim(i2str(prog_style))
                            ierr = 1
                            CHCKERR(err_msg)
                    end select
                case default
                    err_msg = 'No equilibrium style associated with '//&
                        &trim(i2str(eq_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            call lvl_ud(-1)
            call writo('Variables broadcasted')
        end if
    contains
        ! broadcasts the size of an array, so this array can be allocated in the
        ! slave processes
        ! The index of this array is (1:)
        subroutine bcast_size_1_I(arr)                                          ! version with 1 integer argument
            ! input / output
            integer, intent(inout), allocatable :: arr(:)
            
            ! local variables
            integer :: arr_size                                                 ! sent ahead so arrays can be allocated
            
            if (glb_rank.eq.0) arr_size = size(arr)
            
            call MPI_Bcast(arr_size,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            if (glb_rank.ne.0) allocate(arr(1:arr_size))
        end subroutine bcast_size_1_I
        ! The index of this array is (1:)
        subroutine bcast_size_1_R(arr)                                          ! version with 1 real argument
            ! input / output
            real(dp), intent(inout), allocatable :: arr(:)
            
            ! local variables
            integer :: arr_size                                                 ! sent ahead so arrays can be allocated
            
            if (glb_rank.eq.0) arr_size = size(arr)
            
            call MPI_Bcast(arr_size,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            if (glb_rank.ne.0) allocate(arr(1:arr_size))
        end subroutine bcast_size_1_R
        ! The index of this array is (1:)
        subroutine bcast_size_1_var_1D(arr)                                     ! version with 1D var argument (see HDF5_ops)
            ! input / output
            type(var_1D), intent(inout), allocatable :: arr(:)
            
            ! local variables
            integer :: arr_size                                                 ! sent ahead so arrays can be allocated
            
            if (glb_rank.eq.0) arr_size = size(arr)
            
            call MPI_Bcast(arr_size,1,MPI_INTEGER,0,MPI_Comm_world,ierr)
            if (glb_rank.ne.0) allocate(arr(1:arr_size))
        end subroutine bcast_size_1_var_1D
        ! The array index is (1:,1:)
        subroutine bcast_size_2_R(arr)                                          ! version with 2 real arguments
            ! input / output
            real(dp), intent(inout), allocatable :: arr(:,:)
            
            ! local variables
            integer :: arr_size(2)                                              ! sent ahead so arrays can be allocated
            
            if (glb_rank.eq.0) arr_size = [size(arr,1),size(arr,2)]
            
            call MPI_Bcast(arr_size,2,MPI_INTEGER,0,MPI_Comm_world,ierr)
            if (glb_rank.ne.0) allocate(arr(1:arr_size(1),1:arr_size(2)))
        end subroutine bcast_size_2_R
        ! The array index is (0:mpol-1,-ntor:ntor,1:,start_id)
        subroutine bcast_size_4_R(arr,start_id)                                 ! version with 4 real arguments
            ! input / output
            real(dp), intent(inout), allocatable :: arr(:,:,:,:)
            integer, intent(in) :: start_id
            
            ! local variables
            integer :: arr_size(4)                                              ! sent ahead so arrays can be allocated
            
            if (glb_rank.eq.0) arr_size = [size(arr,1),size(arr,2),&
                &size(arr,3),size(arr,4)]
            
            call MPI_Bcast(arr_size,4,MPI_INTEGER,0,MPI_Comm_world,ierr)
            if (glb_rank.ne.0) allocate(arr(0:arr_size(1)-1,&
                &-(arr_size(2)-1)/2:(arr_size(2)-1)/2,1:arr_size(3),&
                &start_id:arr_size(4)+start_id-1))
        end subroutine bcast_size_4_R
        ! The array index is (0:mpol-1,-ntor:ntor,1:)
        subroutine bcast_size_3_R(arr)                                          ! version with 3 real arguments
            ! input / output
            real(dp), intent(inout), allocatable :: arr(:,:,:)
            
            ! local variables
            integer :: arr_size(3)                                              ! sent ahead so arrays can be allocated
            
            if (glb_rank.eq.0) arr_size = [size(arr,1),size(arr,2),size(arr,3)]
            
            call MPI_Bcast(arr_size,3,MPI_INTEGER,0,MPI_Comm_world,ierr)
            if (glb_rank.ne.0) allocate(arr(0:arr_size(1)-1,&
                &-(arr_size(2)-1)/2:(arr_size(2)-1)/2,1:arr_size(3)))
        end subroutine bcast_size_3_R
    end function broadcast_input_vars
end module MPI_ops
