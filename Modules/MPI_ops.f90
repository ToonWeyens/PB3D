!------------------------------------------------------------------------------!
!   Routines related to MPI                                                    !
!------------------------------------------------------------------------------!
module MPI_ops
#include <PB3D_macros.h>
    use MPI
    use str_ops, only: i2str
    use message_ops, only: writo, lvl_ud, print_ar_1
    use num_vars, only: dp, max_str_ln
    use output_ops, only: print_GP_2D, draw_GP
    
    implicit none
    private
    public start_MPI, stop_MPI, split_MPI, abort_MPI, broadcast_vars, &
        &merge_MPI, get_next_job, divide_grid, get_ser_var, get_ghost_X_vec, &
        &wait_MPI
    
    ! interfaces
    interface get_ser_var
        module procedure get_ser_var_complex, get_ser_var_real
    end interface
    
contains
    ! start MPI and gather information
    ! [MPI] Collective call
    integer function start_MPI() result(ierr)
        use num_vars, only: glb_rank, glb_n_procs, grp_rank
        character(*), parameter :: rout_name = 'start_MPI'
        
        ! initialize ierr
        ierr = 0
        
        ! start MPI
        call MPI_init(ierr)                                                     ! initialize MPI
        CHCKERR('MPI init failed')
        call MPI_Comm_rank(MPI_COMM_WORLD,glb_rank,ierr)                        ! global rank
        CHCKERR('MPI rank failed')
        call MPI_Comm_size(MPI_COMM_WORLD,glb_n_procs,ierr)                     ! nr. processes
        CHCKERR('MPI size failed')
        
        ! no alpha groups yet, so set group rank to global rank
        grp_rank = glb_rank
    end function start_MPI
    
    ! determine   how   to   split    the   communicator   MPI_COMM_WORLD   into
    ! subcommunicators,  according to  n_procs_per_alpha,  which determines  how
    ! many processors to use per field line
    ! [MPI] Collective call
    integer function split_MPI() result(ierr)
        use num_vars, only: n_procs_per_alpha, n_procs, n_alpha, min_n_r_X, &
            &MPI_Comm_groups, glb_rank, glb_n_procs, grp_rank, next_job, &
            &grp_n_procs, grp_nr, n_groups, next_job_win
        !use num_vars, only: MPI_Comm_masters
        use file_ops, only: open_output
        
        character(*), parameter :: rout_name = 'split_MPI'
        
        ! local variables
        integer :: remainder                                                    ! remainder of division
        integer :: id                                                           ! counter
        integer :: sum_groups                                                   ! sum of previous colros
        logical :: grp_found                                                    ! when searching for group of local process
        character(len=max_str_ln) :: str_1, str_2, str_3                        ! strings used in user messages
        !integer, allocatable :: glb_grp_master_rank(:)                          ! global ranks of the group masters
        !integer :: world_group                                                  ! MPI group associated to MPI_COMM_WORLD
        !integer :: master_group                                                 ! MPI group associated to alpha group masters
        integer :: intsize                                                      ! size of MPI real
        integer(kind=MPI_ADDRESS_KIND) :: size_one = 1                          ! size equal to 1
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! determine how many  subcommunicators can be made with  n_procs. If the
        ! remainder of  the division is not  zero, assign an extra  processor to
        ! the first groups
        
        ! set n_procs to n_procs per alpha
        if (glb_n_procs.lt.n_procs_per_alpha) n_procs_per_alpha = glb_n_procs   ! too few processors for even one group
        
        n_groups = glb_n_procs/n_procs_per_alpha                                ! how many groups of alpha
        if (n_groups.gt.n_alpha) n_groups = n_alpha                             ! too many groups for alpha's -> limit
        
        allocate(n_procs(n_groups))                                             ! allocate n_procs
        n_procs = n_procs_per_alpha
        
        ! add an extra process if the remainder is not zero
        remainder = glb_n_procs - n_groups*n_procs_per_alpha
        id = 1
        do while (remainder.gt.0)
            n_procs(id) = n_procs(id)+1
            if (id.lt.n_groups) then                                            ! go to next group
                id = id+1
            else                                                                ! back to first group
                id = 1
            end if
            remainder = remainder - 1
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
        call writo(trim(str_1)//' in total, for '//trim(str_2)//trim(str_3))    ! how many groups
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
        grp_rank = glb_rank                                                     ! start grp_rank at global rank
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
        
        ! split MPI_COMM_WORLD according to n_procs
        call MPI_Comm_split(MPI_COMM_WORLD,grp_nr,grp_rank,MPI_Comm_groups,&
            &ierr)
        CHCKERR('Failed to split in groups')
        
        ! increment n_r_X if lower than grp_n_procs
        if (grp_n_procs.gt.min_n_r_X) then
            call writo('WARNING: The starting nr. of normal points of the &
                &perturbation is increased from '//trim(i2str(min_n_r_X))//&
                &' to '//trim(i2str(grp_n_procs))//' because there are too &
                &many processes')
            min_n_r_X = grp_n_procs
        end if
        
        !! take subset of MPI_COMM_WORLD  containing all group masters with their
        !! ranks according to the group rank
        !allocate(glb_grp_master_rank(0:n_groups-1))                             ! allocate glb_grp_master_rank
        !glb_grp_master_rank(0) = 0                                              ! master of first group is also global master
        !do id = 1,n_groups-1
            !glb_grp_master_rank(id) = glb_grp_master_rank(id-1) + n_procs(id)
        !end do
        !call MPI_Comm_group(MPI_COMM_WORLD,world_group,ierr)                    ! get group of MPI_COMM_WORLD
        !CHCKERR('Failed to get global group')
        !call MPI_group_incl(world_group,n_groups,glb_grp_master_rank,&
            !&master_group,ierr)                                                 ! take master subset
        !CHCKERR('Failed to create group of group masters')
        !call MPI_Comm_create(MPI_COMM_WORLD,master_group,MPI_Comm_masters,ierr) ! create communicator for master subset
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
                &MPI_COMM_WORLD,next_job_win,ierr)
        else
            call MPI_Win_create(0,0*size_one,intsize,MPI_INFO_NULL,&
                &MPI_COMM_WORLD,next_job_win,ierr)
        end if
        CHCKERR('Couldn''t create window to next_job')
        
        ! calculate the r range for the equilibrium calculations
        ierr = calc_eq_r_range()
        CHCKERR('')
        
        ! set fence so that the global master holds next_job = 1 for sure
        call MPI_Win_fence(0,next_job_win,ierr)
        CHCKERR('Couldn''t set fence')
        
        ! open an output file for the groups that are not the master group (i.e.
        ! the group containing the global master)
        if (grp_nr.ne.0) then
            ierr = open_output()
            CHCKERR('')
        end if
    contains
        ! calculate grp_min_r_eq  and grp_max_r_eq  for this  rank in  the alpha
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
        ! the  possible  perturbation  points,  at n_r_X  at  its  lowest  value
        ! (i.e.  given by  divide_grid with  cumul .true.)  plus 1,  because the
        ! perturbation  quantity of  perturbation  normal point  (i) depends  on
        ! perturbation normal point (i+1)
        ! Furthermore,  the conversion  between points  r_X on  the perturbation
        ! range  (0..1) and  discrete grid  points r_eq  on the  equilbrium grid
        ! (1..n_r_eq), the  subroutine con2dis  is used with  equilibrium values
        ! tabulated in flux_eq (normalized)
        ! [MPI] Collective call
        integer function calc_eq_r_range() result(ierr)
            use num_vars, only: min_n_r_X, grp_n_procs, grp_rank, min_r_X, &
                &max_r_X, use_pol_flux, eq_style
            use utilities, only: con2dis, dis2con, calc_int, interp_fun_1D, &
                &calc_deriv, round_with_tol
            use eq_vars, only: grp_min_r_eq, grp_max_r_eq, n_r_eq, &
                &eq_use_pol_flux
            use VMEC_vars, only: phi, phi_r, iotaf
            use HEL_vars, only: flux_H, qs
            use X_vars, only: grp_max_r_X
            
            character(*), parameter :: rout_name = 'calc_eq_r_range'
            
            ! local variables
            real(dp), allocatable :: flux(:), flux_eq(:)                        ! either pol. or tor. flux, in VMEC coord.
            real(dp), allocatable :: flux_H_r(:)                                ! normal derivative of flux_H
            real(dp) :: grp_min_r_eq_X_con                                      ! grp_min_r_eq in continuous perturbation grid
            real(dp) :: grp_min_r_eq_eq_con                                     ! grp_min_r_eq in continuous equilibrium grid
            real(dp) :: grp_min_r_eq_eq_dis                                     ! grp_min_r_eq in discrete equilibrium grid, unrounded
            integer :: grp_max_r_eq_X_dis                                       ! grp_max_r_eq in discrete perturbation grid
            real(dp) :: grp_max_r_eq_X_con                                      ! grp_max_r_eq in continuous perturbation grid
            real(dp) :: grp_max_r_eq_eq_con                                     ! grp_max_r_eq in continuous equilibrium grid
            real(dp) :: grp_max_r_eq_dis                                        ! grp_max_r_eq in discrete equilibrium grid, unrounded
            
            ! initialize ierr
            ierr = 0
            
            ! set up  flux and flux_eq, depending on which  equilibrium style is
            ! being used:
            !   1:  VMEC
            !   2:  HELENA
            allocate(flux(n_r_eq),flux_eq(n_r_eq))
            select case (eq_style)
                case (1)                                                        ! VMEC
                    ! set up perturbation flux
                    if (use_pol_flux) then
                        ierr = calc_int(-iotaf*phi_r,1.0_dp/(n_r_eq-1.0_dp),&
                            &flux)
                        CHCKERR('')
                    else
                        flux = phi
                    end if
                    ! set up equilibrium flux
                    if (eq_use_pol_flux) then
                        ierr = calc_int(-iotaf*phi_r,1.0_dp/(n_r_eq-1.0_dp),&
                            &flux_eq)
                        CHCKERR('')
                    else
                        flux_eq = phi
                    end if
                case (2)                                                        ! HELENA
                    ! calculate normal derivative of flux_H
                    allocate(flux_H_r(n_r_eq))
                    ierr = calc_deriv(flux_H,flux_H_r,flux_H,1,1)
                    CHCKERR('')
                    ! set up perturbation flux
                    if (use_pol_flux) then
                        flux = flux_H
                    else
                        ierr = calc_int(qs*flux_H_r,flux_H,flux)
                        CHCKERR('')
                    end if
                    ! set up equilibrium flux
                    if (eq_use_pol_flux) then
                        flux_eq = flux_H
                    else
                        ierr = calc_int(qs*flux_H_r,flux_H,flux_eq)
                        CHCKERR('')
                    end if
                case default
                    err_msg = 'No equilibrium style associated with '//&
                        &trim(i2str(eq_style))
                    ierr = 1
                    CHCKERR(err_msg)
            end select
            
            ! normalize flux and flux_eq to (0..1)
            ! Note:  max_flux is  not yet  determined, so  use flux(n_r_eq)  and
            ! flux_eq(n_r_eq), knowing that this is the full global range
            flux = flux/flux(n_r_eq)
            flux_eq = flux_eq/flux_eq(n_r_eq)
            
            ! use min_r_X and max_r_X, with grp_n_procs to get the minimum bound
            ! grp_min_r_eq for this rank
            ! 1. continuous perturbation grid (0..1)
            grp_min_r_eq_X_con = min_r_X + &
                &grp_rank*(max_r_X-min_r_X)/grp_n_procs
            ! 2 round with tolerance
            ierr = round_with_tol(grp_min_r_eq_X_con,0._dp,1._dp)
            CHCKERR('')
            ! 3. continuous equilibrium grid (0..1)
            ierr = interp_fun_1D(grp_min_r_eq_eq_con,flux_eq,&
                &grp_min_r_eq_X_con,flux)
            CHCKERR('')
            ! 4. discrete equilibrium grid, unrounded
            call con2dis(grp_min_r_eq_eq_con,grp_min_r_eq_eq_dis,flux_eq)
            ! 5. discrete equilibrium grid, rounded down
            grp_min_r_eq = floor(grp_min_r_eq_eq_dis)
            
            ! use min_r_X and max_r_X to calculate grp_max_r_eq
            ! 1. divide_grid for min_n_r_X normal points
            ierr = divide_grid(min_n_r_X)                                       ! divide the grid for the min_n_r_X tot. normal points
            CHCKERR('')
            ! 2. discrete perturbation grid (1..min_n_r_X)
            grp_max_r_eq_X_dis = grp_max_r_X
            ! 3. add one to max if not last global point
            if (grp_rank.ne.grp_n_procs-1) &
                &grp_max_r_eq_X_dis = grp_max_r_eq_X_dis + 1
            ! 4. continous perturbation grid (0..1)
            call dis2con(grp_max_r_eq_X_dis,grp_max_r_eq_X_con,[1,min_n_r_X],&
                &[min_r_X,max_r_X])                                             ! the perturbation grid is equidistant
            ! 5 round with tolerance
            ierr = round_with_tol(grp_max_r_eq_X_con,0._dp,1._dp)
            CHCKERR('')
            ! 6. continuous equilibrium grid (0..1)
            ierr = interp_fun_1D(grp_max_r_eq_eq_con,flux_eq,&
                &grp_max_r_eq_X_con,flux)
            CHCKERR('')
            ! 7. discrete equilibrium grid, unrounded
            call con2dis(grp_max_r_eq_eq_con,grp_max_r_eq_dis,flux_eq)
            ! 8. discrete equlibrium grid, rounded up
            grp_max_r_eq = ceiling(grp_max_r_eq_dis)
            
            deallocate(flux,flux_eq)
        end function calc_eq_r_range
    end function split_MPI
    
    ! merge the MPI groups back to MPI_COMM_WORLD
    ! [MPI] Collective call
    integer function merge_MPI() result(ierr)
        use num_vars, only: MPI_Comm_groups, MPI_Comm_masters, n_groups, &
            &grp_rank, next_job_win, glb_rank
        
        character(*), parameter :: rout_name = 'merge_MPI'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! free window for next_job
        call MPI_Win_free(next_job_win,ierr)
        CHCKERR('Unable to free window for next job')
        
        ! free groups and masters communicators
        call MPI_Comm_free(MPI_Comm_groups,ierr)
        CHCKERR('Unable to free communicator for groups')
        if (n_groups.gt.1) then
            if (grp_rank.eq.0) then
                call MPI_Comm_free(MPI_Comm_masters,ierr)
                err_msg = 'Unable to free communicator for masters'
                CHCKERR(err_msg)
            end if
        end if
        
        ! reset meaningless grp_rank to glb_rank (for output, etc...)
        grp_rank = glb_rank
        
        ! barrier
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        CHCKERR('Coulnd''t set barrier')
    end function merge_MPI
    
    ! MPI Barrier
    integer function wait_MPI() result(ierr)
        use num_vars, only: MPI_Comm_groups
        
        character(*), parameter :: rout_name = 'wait_MPI'
        
        ! initialize ierr
        ierr = 0
        
        call MPI_Barrier(MPI_Comm_groups,ierr)
        CHCKERR('MPI Barrier failed')
    end function wait_MPI
    
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
        
        call MPI_Abort(MPI_COMM_WORLD,0,ierr)
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
            call MPI_Win_lock(MPI_LOCK_EXCLUSIVE,0,0,next_job_win,ierr)
            err_msg = 'Group '//trim(i2str(grp_nr))//&
                &' coulnd''t lock window on global master'
            CHCKERR(err_msg)
            call MPI_Get_accumulate(1,1,MPI_INTEGER,next_job,1,&
                &MPI_INTEGER,0,null_disp,1,MPI_INTEGER,MPI_SUM,&
                &next_job_win,ierr)
            err_msg = 'Group '//trim(i2str(grp_nr))//&
                &' coulnd''t increment next_job on global master'
            CHCKERR(err_msg)
            call MPI_Win_unlock(0,next_job_win,ierr)
            err_msg = 'Group '//trim(i2str(grp_nr))//&
                &' coulnd''t unlock window on global master'
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
    
    ! Gather parallel variable in serial version on group master
    ! Note: the serial variable has to be allocatable and unallocated
    ! [MPI] Collective call
    integer function get_ser_var_complex(var,ser_var) result(ierr)              ! complex version
        use num_vars, only: MPI_Comm_groups, grp_rank, grp_n_procs
        
        character(*), parameter :: rout_name = 'get_ser_var'
        
        ! input / output
        complex(dp), intent(in) :: var(:)                                       ! parallel vector
        complex(dp), allocatable, intent(inout) :: ser_var(:)                   ! serial vector
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer, allocatable :: recvcounts(:)                                   ! counts of nr. of elements received from each processor
        integer, allocatable :: displs(:)                                       ! displacements elements received from each processor
        integer :: id                                                           ! counter
        
        ! initialize ierr
        ierr = 0
        
        ! gather local size  of var of all groups onto  main processor, to serve
        ! as receive counts on group master
        if (grp_rank.eq.0) then
            allocate(recvcounts(grp_n_procs))
            allocate(displs(grp_n_procs))
        else
            allocate(recvcounts(0))
            allocate(displs(0))
        end if
        call MPI_Gather(size(var),1,MPI_INTEGER,recvcounts,1,MPI_INTEGER,0,&
            &MPI_Comm_groups,ierr)
        err_msg = 'Failed to gather size of parallel variable'
        CHCKERR(err_msg)
        
        ! allocate serial variable
        allocate(ser_var(sum(recvcounts)),stat=ierr)
        err_msg = 'Serial variable was already allocated'
        CHCKERR(err_msg)
        
        ! deduce displacements by summing recvcounts
        if (grp_rank.eq.0) then
            displs(1) = 0
            do id = 2,grp_n_procs
                displs(id) = displs(id-1) + recvcounts(id-1)
            end do
        end if
        
        call MPI_Gatherv(var,size(var),MPI_DOUBLE_COMPLEX,ser_var,&
            &recvcounts,displs,MPI_DOUBLE_COMPLEX,0,MPI_Comm_groups,ierr)
        err_msg = 'Failed to gather parallel variable'
        CHCKERR(err_msg)
    end function get_ser_var_complex
    integer function get_ser_var_real(var,ser_var) result(ierr)                 ! real version
        use num_vars, only: MPI_Comm_groups, grp_rank, grp_n_procs
        
        character(*), parameter :: rout_name = 'get_ser_var_real'
        
        ! input / output
        real(dp), intent(in) :: var(:)                                          ! parallel vector
        real(dp), allocatable, intent(inout) :: ser_var(:)                      ! serial vector
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer, allocatable :: recvcounts(:)                                   ! counts of nr. of elements received from each processor
        integer, allocatable :: displs(:)                                       ! displacements elements received from each processor
        integer :: id                                                           ! counter
        
        ! initialize ierr
        ierr = 0
        
        ! gather local size  of var of all groups onto  main processor, to serve
        ! as receive counts on group master
        if (grp_rank.eq.0) then
            allocate(recvcounts(grp_n_procs))
            allocate(displs(grp_n_procs))
        else
            allocate(recvcounts(0))
            allocate(displs(0))
        end if
        call MPI_Gather(size(var),1,MPI_INTEGER,recvcounts,1,MPI_INTEGER,0,&
            &MPI_Comm_groups,ierr)
        err_msg = 'Failed to gather size of parallel variable'
        CHCKERR(err_msg)
        
        ! allocate serial variable
        allocate(ser_var(sum(recvcounts)),stat=ierr)
        err_msg = 'Serial variable was already allocated'
        CHCKERR(err_msg)
        
        ! deduce displacements by summing recvcounts
        if (grp_rank.eq.0) then
            displs(1) = 0
            do id = 2,grp_n_procs
                displs(id) = displs(id-1) + recvcounts(id-1)
            end do
        end if
        
        call MPI_Gatherv(var,size(var),MPI_DOUBLE_PRECISION,ser_var,&
            &recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_Comm_groups,ierr)
        err_msg = 'Failed to gather parallel variable'
        CHCKERR(err_msg)
    end function get_ser_var_real
    
    ! fill the  ghost regions  in  X_vec. Every  message  is  identified by  its
    ! sending process
    ! [MPI] Collective call
    integer function get_ghost_X_vec(X_vec,size_ghost) result(ierr)
        use num_vars, only: MPI_Comm_groups, grp_rank, grp_n_procs
        
        character(*), parameter :: rout_name = 'get_ghost_X_vec'
        
        ! input / output
        complex(dp), intent(in) :: X_vec(:,:)                                   ! parallel vector
        integer, intent(in) :: size_ghost                                       ! width of ghost region
        
        ! local variables
        integer :: n_modes                                                      ! number of modes
        integer :: tot_size                                                     ! total size (including ghost region)
        integer :: istat(MPI_STATUS_SIZE)                                       ! status of send-receive
        
        ! initialize ierr
        ierr = 0
        
        ! initialize number of modes and total size
        n_modes = size(X_vec,1)
        tot_size = size(X_vec,2)
        
        ! ghost regions only make sense if there is more than 1 process
        if (grp_n_procs.gt.1) then
            if (grp_rank.eq.0) then                                             ! first rank only receives
                call MPI_Recv(X_vec(:,tot_size-size_ghost+1:tot_size),&
                    &size_ghost*n_modes,MPI_DOUBLE_COMPLEX,grp_rank+1,&
                    &grp_rank+1,MPI_Comm_groups,istat,ierr)
                CHCKERR('Failed to receive')
            else if (grp_rank+1.eq.grp_n_procs) then                            ! last rank only sends
                call MPI_Send(X_vec(:,1:size_ghost),size_ghost*n_modes,&
                    &MPI_DOUBLE_COMPLEX,grp_rank-1,grp_rank,MPI_Comm_groups,&
                    &ierr)
                CHCKERR('Failed to send')
            else                                                                ! middle ranks send and receive
                call MPI_Sendrecv(X_vec(:,1:size_ghost),size_ghost*n_modes,&
                    &MPI_DOUBLE_COMPLEX,grp_rank-1,grp_rank,&
                    &X_vec(:,tot_size-size_ghost+1:tot_size),size_ghost*n_modes,&
                    &MPI_DOUBLE_COMPLEX,grp_rank+1,grp_rank+1,MPI_Comm_groups,&
                    &istat,ierr)
                CHCKERR('Failed to send and receive')
            end if
        end if
    end function get_ghost_X_vec
    
    ! divides a  grid of  n_r_X points  under the  ranks of  MPI_Comm_groups and
    ! assigns grp_n_r_X and grp_min_r_X and  grp_max_r_X to each rank. Also sets
    ! up  grp_r_X, which contains the  normal variable in the  perturbation grid
    ! for this rank (global range (min_r_X..max_r_X))
    integer function divide_grid(n_r_X_in) result(ierr)
        use num_vars, only: MPI_Comm_groups, min_r_X, max_r_X, grp_rank, &
            &grp_n_procs
        use X_vars, only: n_r_X, grp_n_r_X, grp_min_r_X, grp_max_r_X, grp_r_X
        use utilities, only: round_with_tol
        
        character(*), parameter :: rout_name = 'divide_grid'
        
        ! input / output
        integer, intent(in), optional :: n_r_X_in                               ! custom user-provided number n_r_X
        
        ! local variables
        integer :: n_procs                                                      ! nr. of procs.
        integer :: rank                                                         ! rank
        integer :: id                                                           ! counter
        integer :: n_r_X_loc                                                    ! local value of n_r_X, either user-provided or from X_vars
        integer :: grp_n_r_X_one                                                ! either grp_n_r_X + 1 (first ranks) or grp_n_r_X (last rank)
        
        ! initialize ierr
        ierr = 0
        
        ! set n_r_X_loc
        if (present(n_r_X_in)) then
            n_r_X_loc = n_r_X_in
        else
            n_r_X_loc = n_r_X
        end if
        
        ! set n_procs and rank
        call MPI_Comm_size(MPI_Comm_groups,n_procs,ierr)
        CHCKERR('Failed to get MPI size')
        call MPI_Comm_rank(MPI_Comm_groups,rank,ierr)
        CHCKERR('Failed to get MPI rank')
        
        ! calculate n_loc for this rank
        grp_n_r_X = divide_grid_ind(rank,n_r_X_loc,n_procs)
        
        ! calculate the starting index of this rank
        grp_min_r_X = 1
        do id = 0,rank-1
            grp_min_r_X = grp_min_r_X + divide_grid_ind(id,n_r_X_loc,n_procs)
        end do
        ! calculate the end index of this rank
        grp_max_r_X = grp_min_r_X - 1 + grp_n_r_X
        
        ! set up grp_r_X
        ! Note: for the first ranks, the upper index is one higher than might be
        ! expected because  the routine fill_matrix needs  information about the
        ! next perturbation point (so this is an asymetric ghost region)
        if (allocated(grp_r_X)) deallocate(grp_r_X)
        if (grp_rank+1.lt.grp_n_procs) then
            grp_n_r_X_one = grp_n_r_X + 1
        else
            grp_n_r_X_one = grp_n_r_X
        end if
        allocate(grp_r_X(grp_n_r_X_one))
        grp_r_X = [(min_r_X + (grp_min_r_X+id-2.0_dp)/(n_r_X_loc-1.0_dp)*&
            &(max_r_X-min_r_X),id=1,grp_n_r_X_one)]
        
        ! round with standard tolerance
        ierr = round_with_tol(grp_r_X,0.0_dp,1.0_dp)
        CHCKERR('')
    contains 
        integer function divide_grid_ind(rank,n,n_procs) result(n_loc)
            ! input / output
            integer, intent(in) :: rank                                         ! rank for which to be calculate divide_grid_ind
            integer, intent(in) :: n                                            ! tot. nr. of points to be divided under n_procs
            integer, intent(in) :: n_procs                                      ! nr. of processes
            
            n_loc = n/n_procs                                                   ! number of radial points on this processor
            if (mod(n,n_procs).gt.0) then                                       ! check if there is a remainder
                if (mod(n,n_procs).gt.rank) &
                    &n_loc = n_loc + 1                                          ! add a point to the first ranks if there is a remainder
            end if
        end function divide_grid_ind
    end function divide_grid
    
    ! Broadcasts all  the relevant variable that have been  determined so far in
    ! the global master process using the inputs to the other processes
    ! [MPI] Collective call
    integer function broadcast_vars() result(ierr)
        use VMEC_vars, only: mpol, ntor, lasym, lfreeb, nfp, iotaf, gam, R_c, &
            &R_s, Z_c, Z_s, L_c, L_s, phi, phi_r, presf
        use num_vars, only: max_str_ln, output_name, ltest, EV_style, &
            &max_it_NR, max_it_r, n_alpha, n_procs_per_alpha, minim_style, &
            &max_alpha, min_alpha, tol_NR, glb_rank, glb_n_procs, no_guess, &
            &n_sol_requested, min_n_r_X, min_r_X, max_r_X, nyq_fac, tol_r, &
            &use_pol_flux, max_n_plots, plot_grid, no_plots, output_style, &
            &eq_style, use_normalization, n_sol_plotted, n_theta_plot, &
            &n_zeta_plot, max_deriv
        use X_vars, only: min_m_X, max_m_X, min_n_X, max_n_X
        use eq_vars, only: n_par, max_par, min_par, grp_min_r_eq, n_r_eq, &
            &grp_max_r_eq, R_0, pres_0, B_0, psi_0, rho_0, eq_use_pol_flux
        use HEL_vars, only: R_0_H, B_0_H, p0, qs, flux_H, nchi, chi_H, ias, &
            &h_H_11, h_H_12, h_H_33, RBphi, R_H, Z_H
        
        character(*), parameter :: rout_name = 'broadcast_vars'
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        if (glb_n_procs.gt.1) then                                              ! need to broadcast from global rank 0 to other processes
            ! print message
            call writo('Broadcasting variables determined by input')
            call lvl_ud(1)
            
            ! broadcast eq_style to determine what else to broadcast
            call MPI_Bcast(eq_style,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            
            ! variables that are sent for every equilibrium style:
            call MPI_Bcast(output_name,max_str_ln,MPI_CHARACTER,0,&
                &MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(ltest,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(use_pol_flux,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(eq_use_pol_flux,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(no_guess,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(no_plots,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(plot_grid,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(use_normalization,1,MPI_LOGICAL,0,MPI_COMM_WORLD,&
                &ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(max_it_NR,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(max_it_r,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(minim_style,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(n_par,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(n_r_eq,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(EV_style,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(n_procs_per_alpha,1,MPI_INTEGER,0,MPI_COMM_WORLD,&
                &ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(n_alpha,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(min_m_X,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(max_m_X,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(min_n_X,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(max_n_X,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(n_sol_requested,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(min_n_r_X,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(grp_min_r_eq,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(grp_max_r_eq,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(nyq_fac,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(max_n_plots,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(output_style,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(n_theta_plot,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(n_zeta_plot,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(max_deriv,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(n_sol_plotted,4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(min_alpha,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,&
                &ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(max_alpha,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,&
                &ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(min_r_X,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(max_r_X,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(tol_NR,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(tol_r,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(min_par,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(max_par,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(gam,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(R_0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(pres_0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(B_0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(psi_0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(rho_0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            
            ! For  specific variables, choose  which equilibrium style  is being
            ! used:
            !   1:  VMEC
            !   2:  HELENA
            select case (eq_style)
                case (1)                                                        ! VMEC
                    call MPI_Bcast(lasym,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
                    call MPI_Bcast(lfreeb,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
                    call MPI_Bcast(mpol,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
                    call MPI_Bcast(ntor,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
                    call MPI_Bcast(nfp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
                    call bcast_size_1_R(phi)
                    CHCKERR('MPI broadcast failed')
                    call MPI_Bcast(phi,size(phi),MPI_DOUBLE_PRECISION,0,&
                        &MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
                    call bcast_size_1_R(phi_r)
                    CHCKERR('MPI broadcast failed')
                    call MPI_Bcast(phi_r,size(phi_r),MPI_DOUBLE_PRECISION,0,&
                        &MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
                    call bcast_size_1_R(iotaf)
                    CHCKERR('MPI broadcast failed')
                    call MPI_Bcast(iotaf,size(iotaf),MPI_DOUBLE_PRECISION,0,&
                        &MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
                    call bcast_size_4_R(R_c)
                    CHCKERR('MPI broadcast failed')
                    call MPI_Bcast(R_c,size(R_c),MPI_DOUBLE_PRECISION,0,&
                        &MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
                    call bcast_size_4_R(R_s)
                    CHCKERR('MPI broadcast failed')
                    call MPI_Bcast(R_s,size(R_s),MPI_DOUBLE_PRECISION,0,&
                        &MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
                    call bcast_size_4_R(Z_c)
                    CHCKERR('MPI broadcast failed')
                    call MPI_Bcast(Z_c,size(Z_c),MPI_DOUBLE_PRECISION,0,&
                        &MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
                    call bcast_size_4_R(Z_s)
                    CHCKERR('MPI broadcast failed')
                    call MPI_Bcast(Z_s,size(Z_s),MPI_DOUBLE_PRECISION,0,&
                        &MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
                    call bcast_size_4_R(L_c)
                    CHCKERR('MPI broadcast failed')
                    call MPI_Bcast(L_c,size(L_c),MPI_DOUBLE_PRECISION,0,&
                        &MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
                    call bcast_size_4_R(L_s)
                    CHCKERR('MPI broadcast failed')
                    call MPI_Bcast(L_s,size(L_s),MPI_DOUBLE_PRECISION,0,&
                        &MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
                    call bcast_size_1_R(presf)
                    CHCKERR('MPI broadcast failed')
                    call MPI_Bcast(presf,size(presf),MPI_DOUBLE_PRECISION,0,&
                        &MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
                case (2)                                                        ! HELENA
                    call MPI_Bcast(nchi,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
                    call MPI_Bcast(ias,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
                    call bcast_size_1_R(flux_H)
                    CHCKERR('MPI broadcast failed')
                    call MPI_Bcast(R_0_H,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
                    call MPI_Bcast(B_0_H,1,MPI_DOUBLE_PRECISION,0,&
                        &MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
                    call MPI_Bcast(flux_H,size(flux_H),MPI_DOUBLE_PRECISION,0,&
                        &MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
                    call bcast_size_1_R(qs)
                    CHCKERR('MPI broadcast failed')
                    call MPI_Bcast(qs,size(qs),MPI_DOUBLE_PRECISION,0,&
                        &MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
                    call bcast_size_1_R(p0)
                    CHCKERR('MPI broadcast failed')
                    call MPI_Bcast(p0,size(p0),MPI_DOUBLE_PRECISION,0,&
                        &MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
                    call bcast_size_1_R(chi_H)
                    CHCKERR('MPI broadcast failed')
                    call MPI_Bcast(chi_H,size(chi_H),MPI_DOUBLE_PRECISION,0,&
                        &MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
                    call bcast_size_1_R(RBphi)
                    CHCKERR('MPI broadcast failed')
                    call MPI_Bcast(RBphi,size(RBphi),MPI_DOUBLE_PRECISION,0,&
                        &MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
                    call bcast_size_2_R(R_H)
                    CHCKERR('MPI broadcast failed')
                    call MPI_Bcast(R_H,size(R_H),MPI_DOUBLE_PRECISION,0,&
                        &MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
                    call bcast_size_2_R(Z_H)
                    CHCKERR('MPI broadcast failed')
                    call MPI_Bcast(Z_H,size(Z_H),MPI_DOUBLE_PRECISION,0,&
                        &MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
                    call bcast_size_2_R(h_H_11)
                    CHCKERR('MPI broadcast failed')
                    call MPI_Bcast(h_H_11,size(h_H_11),MPI_DOUBLE_PRECISION,0,&
                        &MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
                    call bcast_size_2_R(h_H_12)
                    CHCKERR('MPI broadcast failed')
                    call MPI_Bcast(h_H_12,size(h_H_12),MPI_DOUBLE_PRECISION,0,&
                        &MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
                    call bcast_size_2_R(h_H_33)
                    CHCKERR('MPI broadcast failed')
                    call MPI_Bcast(h_H_33,size(h_H_33),MPI_DOUBLE_PRECISION,0,&
                        &MPI_COMM_WORLD,ierr)
                    CHCKERR('MPI broadcast failed')
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
        ! The index of this array is (1:n_r_eq)
        subroutine bcast_size_1_I(arr)                                          ! version with 1 integer argument
            ! input / output
            integer, intent(inout), allocatable :: arr(:)
            
            ! local variables
            integer :: arr_size                                                 ! sent ahead so arrays can be allocated
            
            if (glb_rank.eq.0) arr_size = size(arr)
            
            call MPI_Bcast(arr_size,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            if (glb_rank.ne.0) allocate(arr(1:arr_size))
        end subroutine bcast_size_1_I
        ! The index of this array is (1:n_r_eq)
        subroutine bcast_size_1_R(arr)                                          ! version with 1 real argument
            ! input / output
            real(dp), intent(inout), allocatable :: arr(:)
            
            ! local variables
            integer :: arr_size                                                 ! sent ahead so arrays can be allocated
            
            if (glb_rank.eq.0) arr_size = size(arr)
            
            call MPI_Bcast(arr_size,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            if (glb_rank.ne.0) allocate(arr(1:arr_size))
        end subroutine bcast_size_1_R
        ! The array index is (1:n_par,1:n_r_eq)
        subroutine bcast_size_2_R(arr)                                          ! version with 2 real arguments
            ! input / output
            real(dp), intent(inout), allocatable :: arr(:,:)
            
            ! local variables
            integer :: arr_size(2)                                              ! sent ahead so arrays can be allocated
            
            if (glb_rank.eq.0) arr_size = [size(arr,1),size(arr,2)]
            
            call MPI_Bcast(arr_size,2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            if (glb_rank.ne.0) allocate(arr(1:arr_size(1),1:arr_size(2)))
        end subroutine bcast_size_2_R
        ! The array index is (0:mpol-1,-ntor:ntor,1:n_r_eq,0:max_deriv)
        subroutine bcast_size_4_R(arr)                                          ! version with 4 real arguments
            ! input / output
            real(dp), intent(inout), allocatable :: arr(:,:,:,:)
            
            ! local variables
            integer :: arr_size(4)                                              ! sent ahead so arrays can be allocated
            
            if (glb_rank.eq.0) arr_size = [size(arr,1),size(arr,2),&
                &size(arr,3),size(arr,4)]
            
            call MPI_Bcast(arr_size,4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            if (glb_rank.ne.0) allocate(arr(0:arr_size(1)-1,&
                &-(arr_size(2)-1)/2:(arr_size(2)-1)/2,1:arr_size(3),&
                &0:arr_size(4)-1))
        end subroutine bcast_size_4_R
    end function broadcast_vars
end module
