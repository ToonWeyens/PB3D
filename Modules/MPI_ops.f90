!------------------------------------------------------------------------------!
!   Routines related to MPI                                                    !
!------------------------------------------------------------------------------!
module MPI_ops
#include <PB3D_macros.h>
    use MPI
    use str_ops, only: i2str
    use output_ops, only: writo, lvl_ud
    use num_vars, only: dp
    
    implicit none
    private
    public start_MPI, stop_MPI, split_MPI, abort_MPI, broadcast_vars, &
        &merge_MPI, get_next_job, divide_grid, broadcast_logical
    
contains
    ! start MPI and gather information
    ! [MPI] Collective call
    integer function start_MPI() result(ierr)
        use num_vars, only: glob_rank, glob_n_procs, group_rank
        character(*), parameter :: rout_name = 'start_MPI'
        
        ! initialize ierr
        ierr = 0
        
        ! start MPI
        call MPI_init(ierr)                                                     ! initialize MPI
        CHCKERR('MPI init failed')
        call MPI_Comm_rank(MPI_COMM_WORLD,glob_rank,ierr)                       ! global rank
        CHCKERR('MPI rank failed')
        call MPI_Comm_size(MPI_COMM_WORLD,glob_n_procs,ierr)                    ! nr. processes
        CHCKERR('MPI size failed')
        
        ! no alpha groups yet, so set group rank to global rank
        group_rank = glob_rank
    end function start_MPI
    
    ! determine   how   to   split    the   communicator   MPI_COMM_WORLD   into
    ! subcommunicators,  according to  n_procs_per_alpha,  which determines  how
    ! many processors to use per field line
    ! [MPI] Collective call
    integer function split_MPI() result(ierr)
        use num_vars, only: n_procs_per_alpha, n_procs, n_alpha, max_str_ln, &
            &MPI_Comm_groups, glob_rank, glob_n_procs, group_rank, next_job, &
            &group_n_procs, group_nr, n_groups,  MPI_Comm_masters, &
            &next_job_win
        use file_ops, only: open_output
        
        character(*), parameter :: rout_name = 'split_MPI'
        
        ! local variables
        integer :: remainder                                                    ! remainder of division
        integer :: id                                                           ! counter
        integer :: sum_groups                                                   ! sum of previous colros
        logical :: group_found                                                  ! when searching for group of local process
        character(len=max_str_ln) :: str_1, str_2, str_3                        ! strings used in user messages
        integer, allocatable :: glob_group_master_rank(:)                       ! global ranks of the group masters
        integer :: world_group                                                  ! MPI group associated to MPI_COMM_WORLD
        integer :: master_group                                                 ! MPI group associated to alpha group masters
        integer :: intsize                                                      ! size of MPI real
        integer(kind=MPI_ADDRESS_KIND) :: size_one = 1                          ! size equal to 1
        character(len=max_str_ln) :: err_msg                                    ! error message
        
        ! initialize ierr
        ierr = 0
        
        ! determine how many  subcommunicators can be made with  n_procs. If the
        ! remainder of  the division is not  zero, assign an extra  processor to
        ! the first groups
        
        ! set n_procs to n_procs per alpha
        if (glob_n_procs.lt.n_procs_per_alpha) n_procs_per_alpha = glob_n_procs ! too few processors for even one group
        
        n_groups = glob_n_procs/n_procs_per_alpha                               ! how many groups of alpha
        if (n_groups.gt.n_alpha) n_groups = n_alpha                             ! too many groups for alpha's -> limit
        
        allocate(n_procs(n_groups))                                             ! allocate n_procs
        n_procs = n_procs_per_alpha
        
        ! add an extra process if the remainder is not zero
        remainder = glob_n_procs - n_groups*n_procs_per_alpha
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
        if (glob_n_procs.eq.1) then
            str_1 = '1 process'
        else
            str_1 = trim(i2str(glob_n_procs))//' processes'
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
        if (glob_n_procs.eq.1) then
            str_3 = ' in serial'
        else
            str_3 = ' in parallel'
        end if
        call writo(trim(str_1)//' a total of '//trim(str_2)//trim(str_3))       ! how many field lines to solve
        
        ! determine  group nr. of  local process, group rank and  group_n_procs
        group_nr = 0                                                            ! start group nr at 0, in MPI spirit
        group_found = .false.                                                   ! group not yet found
        sum_groups = 0                                                          ! start sum of ranks at 0
        group_rank = glob_rank                                                  ! start group_rank at global rank
        do while (.not.group_found)
            sum_groups = sum_groups + n_procs(group_nr+1)                       ! upper bound of nr. processes of groups from 0 to group_nr+1
            if (glob_rank+1.le.sum_groups) then                                 ! group found
                group_found = .true.
                group_n_procs = n_procs(group_nr+1)
            else                                                                ! group not found: try next group
                group_rank = group_rank - n_procs(group_nr+1)                   ! subtract n_procs in this group from group_rank
                group_nr = group_nr + 1                                         ! increment group nr.
            end if
        end do
        
        ! split MPI_COMM_WORLD according to n_procs
        call MPI_Comm_split(MPI_COMM_WORLD,group_nr,group_rank,MPI_Comm_groups,&
            &ierr)
        CHCKERR('Failed to split in groups')
        
        ! take subset of MPI_COMM_WORLD  containing all group masters with their
        ! ranks according to the group rank
        allocate(glob_group_master_rank(0:n_groups-1))                          ! allocate glob_group_master_rank
        glob_group_master_rank(0) = 0                                           ! master of first group is also global master
        do id = 1,n_groups-1
            glob_group_master_rank(id) = glob_group_master_rank(id-1) + &
                &n_procs(id)
        end do
        call MPI_Comm_group(MPI_COMM_WORLD,world_group,ierr)                    ! get group of MPI_COMM_WORLD
        CHCKERR('Failed to get global group')
        call MPI_Group_incl(world_group,n_groups,glob_group_master_rank,&
            &master_group,ierr)                                                 ! take master subset
        CHCKERR('Failed to create group of group masters')
        call MPI_Comm_create(MPI_COMM_WORLD,master_group,MPI_Comm_masters,ierr) ! create communicator for master subset
        err_msg = 'Failed to create communicator to group of group masters'
        CHCKERR(err_msg)
        !!!! MPI_Comm_masters NOT NECESSARY ANYMORE !!!!
        
        ! set starting next_job to 1 on global master
        if (glob_rank.eq.0) then
            next_job = 1
        else
            next_job = 0
        end if
        
        ! create a window to the variable next_job in the global master
        call MPI_Type_extent(MPI_INTEGER,intsize,ierr) 
        err_msg = 'Couldn''t determine the extent of MPI_REAL'
        CHCKERR(err_msg)
        if (group_rank.eq.0) then
            call MPI_Win_create(next_job,1*size_one,intsize,MPI_INFO_NULL,&
                &MPI_COMM_WORLD,next_job_win,ierr)
        else
            call MPI_Win_create(0,0*size_one,intsize,MPI_INFO_NULL,&
                &MPI_COMM_WORLD,next_job_win,ierr)
        end if
        CHCKERR('Couldn''t create window to next_job')
        
        ! calculate the r range for the equilibrium calculations
        ierr = calc_r_range()
        CHCKERR('')
        
        ! set fence so that the global master holds next_job = 1 for sure
        call MPI_Win_fence(0,next_job_win,ierr)
        CHCKERR('Couldn''t set fence')
        
        ! open an output file for the groups that are not the master group (i.e.
        ! the group containing the global master)
        if (group_nr.ne.0) then
            ierr = open_output()
            CHCKERR('')
        end if
    contains
        ! calculate min_r_eq  and max_r_eq for this  rank in the alpha  group so
        ! that the  equilibrium quantities  are calculated  for only  the normal
        ! range that is relevant for the perturbation calculations, which are to
        ! be anticipated.
        ! This is done as follows: if  the nr. of perturbation grid points n_r_X
        ! goes to infinity,  the division under the processes  will approach the
        ! ideal value  of 1/group_n_procs each.  For lower values of  n_r_X, the
        ! first processes can carry an additional perturbation grid point if the
        ! remainder of  the division is not  zero. The situation is  at its most
        ! assymetric for the lowest value of n_r_X.
        ! Therefore,  to be  able  to cover  all the  levels  in the  Richardson
        ! extrapolation  (when n_r_X  is  subsequently  increased), the  minimum
        ! of  the   equilibrium  range  of   each  processor  is  to   be  given
        ! by  the  grid   point  that  includes  the  lowest   of  the  possible
        ! perturbation  points,  at  n_r_X  going  to  infinity  (i.e.  r_min  +
        ! group_rank*(r_max-r_min)/group_n_procs)   and  the   maximum  of   the
        ! equilibrium range is  to be given by the grid  point that includes the
        ! highest of  the possible perturbation  points, at n_r_X at  its lowest
        ! value (i.e. given by divide_grid with  cumul true) plus 1, because the
        ! perturbation  quantity of  perturbation  normal point  (i) depends  on
        ! perturbation normal point (i+1)
        ! Furthermore,  the conversion  between points  r_X on  the perturbation
        ! range  (0..1) and  discrete grid  points r_eq  on the  equilbrium grid
        ! (1..n_r), the subroutine  con2dis is used with  equilibrium limits set
        ! to [1..n_r]
        ! [MPI] Collective call
        integer function calc_r_range() result(ierr)
            use num_vars, only: min_n_r_X, group_n_procs, group_rank, min_r_X, &
                &max_r_X
            use utilities, only: con2dis, dis2con
            use eq_vars, only: min_r_eq, max_r_eq
            use VMEC_vars, only: n_r
            
            character(*), parameter :: rout_name = 'calc_r_range'
            
            ! local variables
            real(dp) :: min_r_eq_X_con                                          ! min_r_eq in continuous perturbation grid
            real(dp) :: min_r_eq_dis                                            ! min_r_eq in discrete equilibrium grid, unrounded
            integer :: max_r_eq_X_dis                                           ! max_r_eq in discrete perturbation grid
            real(dp) :: max_r_eq_X_con                                          ! max_r_eq in continuous perturbation grid
            real(dp) :: max_r_eq_dis                                            ! max_r_eq in discrete equilibrium grid, unrounded
            
            ! initialize ierr
            ierr = 0
            
            ! use min_r_X  and max_r_X,  with group_n_procs  to get  the minimum
            ! bound for this rank
            min_r_eq_X_con = min_r_X + &
                &group_rank*(max_r_X-min_r_X)/group_n_procs                     ! min_r_eq in continuous perturbation grid (0..1)
            call con2dis(min_r_eq_X_con,[0._dp,1._dp],min_r_eq_dis,[1,n_r])     ! min_r_eq in discrete equilibrium grid, unrounded
            min_r_eq = floor(min_r_eq_dis)                                      ! rounded up
            
            ! use grid_max with min_n_r_X to calculate max_r_eq
            ierr = divide_grid(MPI_Comm_groups,min_n_r_X,max_r_eq_X_dis,&
                &cumul=.true.)                                                  ! max_r_eq in discrete perturbation grid (1..min_n_r_X)
            CHCKERR('')
            if (group_rank.ne.group_n_procs-1) &
                &max_r_eq_X_dis = max_r_eq_X_dis + 1                            ! only add 1 if not last global point
            call dis2con(max_r_eq_X_dis,[1,min_n_r_X],max_r_eq_X_con,&
                &[0._dp,1._dp])                                                 ! max_r_eq in continuous perturbation grid (min_r_X..max_r_X)
            max_r_eq_X_con = min_r_X + (max_r_X-min_r_X)*max_r_eq_X_con         ! max_r_eq in continuous perturbation grid (0..1)
            call con2dis(max_r_eq_X_con,[0._dp,1._dp],max_r_eq_dis,[1,n_r])     ! max_r_eq in discrete equilibrium grid, unrounded
            max_r_eq = ceiling(max_r_eq_dis)                                    ! round down
        end function calc_r_range
    end function split_MPI
    
    ! merge the MPI groups back to MPI_COMM_WORLD
    ! [MPI] Collective call
    integer function merge_MPI() result(ierr)
        use num_vars, only: MPI_Comm_groups, MPI_Comm_masters, n_groups, &
            &group_rank, next_job_win, glob_rank, max_str_ln
        
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
            if (group_rank.eq.0) then
                call MPI_Comm_free(MPI_Comm_masters,ierr)
                err_msg = 'Unable to free communicator for masters'
                CHCKERR(err_msg)
            end if
        end if
        
        ! reset meaningless group_rank to glob_rank (for output, etc...)
        group_rank = glob_rank
        
        ! barrier
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        CHCKERR('Coulnd''t set barrier')
        
        !!! ALSO CLOSE THEIR OUTPUT FILES, MAYBE DELETE?!?!?!
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
        
        call MPI_Abort(MPI_COMM_WORLD,0,ierr)
        CHCKERR('MPI abort failed')
    end function abort_MPI
    
    ! compares  own next_job with next_job  of other group masters  and gets the
    ! maximum value. Returns this value
    ! [MPI] Collective call
    integer function get_next_job(next_job) result(ierr)
        use num_vars, only: group_rank, max_str_ln, next_job_win, &
            &MPI_Comm_groups, n_alpha, group_nr
        
        character(*), parameter :: rout_name = 'get_next_job'
        
        ! input / output
        integer, intent(inout) :: next_job
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer(kind=MPI_ADDRESS_KIND) :: null_disp = 0                         ! zero target displacement
        
        ! initialize ierr
        ierr = 0
        
        ! test whether it is a master that calls this routine
        if (group_rank.eq.0) then
            call MPI_Win_lock(MPI_LOCK_EXCLUSIVE,0,0,next_job_win,ierr)
            err_msg = 'Group '//trim(i2str(group_nr))//&
                &' coulnd''t lock window on global master'
            CHCKERR(err_msg)
            call MPI_Get_accumulate(1,1,MPI_INTEGER,next_job,1,&
                &MPI_INTEGER,0,null_disp,1,MPI_INTEGER,MPI_SUM,&
                &next_job_win,ierr)
            err_msg = 'Group '//trim(i2str(group_nr))//&
                &' coulnd''t increment next_job on global master'
            CHCKERR(err_msg)
            call MPI_Win_unlock(0,next_job_win,ierr)
            err_msg = 'Group '//trim(i2str(group_nr))//&
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
    
    ! divides a grid  of n points under  the ranks of MPI  Communicator comm and
    ! assigns n_loc to each rank
    integer function divide_grid(comm,n,n_loc,cumul) result(ierr)
        character(*), parameter :: rout_name = 'divide_grid'
        
        ! input / output
        integer, intent(in) :: comm                                             ! MPI communicator for the group
        integer, intent(in) :: n                                                ! number of grid points to be divided
        integer, intent(inout) :: n_loc                                         ! number of grid points for this rank
        logical, intent(in), optional :: cumul                                  ! if .true., n_loc is given cumulative
        
        ! local variables
        integer :: n_procs                                                      ! nr. of procs.
        integer :: rank                                                         ! rank
        logical :: cumul_loc                                                    ! local version of cumul
        integer :: id                                                           ! counter
        
        ! initialize ierr
        ierr = 0
        
        ! initialize cumul_loc
        if (present(cumul)) then
            cumul_loc = cumul
        else
            cumul_loc = .false.
        end if
        
        ! set n_procs and rank
        call MPI_Comm_size(comm,n_procs,ierr)
        CHCKERR('Failed to get MPI size')
        call MPI_Comm_rank(comm,rank,ierr)
        CHCKERR('Failed to get MPI rank')
        
        if (cumul_loc) then
            n_loc = 0
            do id = 0,rank
                n_loc = n_loc + divide_grid_ind(id)
            end do
        else
            n_loc = divide_grid_ind(rank)
        end if
    contains 
        integer function divide_grid_ind(rank) result(n_loc)
            ! input / output
            integer, intent(in) :: rank                                         ! rank for which to be calculate divide_grid_ind
            
            n_loc = n/n_procs                                                   ! number of radial points on this processor
            if (mod(n,n_procs).gt.0) then                                       ! check if there is a remainder
                if (mod(n,n_procs).gt.rank) &
                    &n_loc = n_loc + 1                                          ! add a point to the first ranks if there is a remainder
            end if
        end function divide_grid_ind
    end function divide_grid
    
    ! broadcasts a logical variable inside an alpha group
    ! [MPI] Collective call
    integer function broadcast_logical(comm,flag) result(ierr)
        character(*), parameter :: rout_name = 'divide_grid'
        
        ! input / output
        integer, intent(in) :: comm                                             ! MPI communicator for the group
        logical, intent(in) :: flag                                             ! logical to be broadcast
        
        ! initialize ierr
        ierr = 0
        
        call MPI_Bcast(flag,1,MPI_LOGICAL,0,comm,ierr)
        CHCKERR('MPI broadcast failed')
    end function broadcast_logical
        
    ! THIS SHOULD BE DONE WITH A STRUCT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! BUT DON'T FORGET TO ALLOCATE BEFORE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Broadcasts all  the relevant variable that have been  determined so far in
    ! the global master process using the inputs to the other processes
    ! The messages consist of (from open_input, default_input and read_VMEC):
    !   1   character(len=max_str_ln)   output_name
    !   2   logical                     ltest
    !   3   logical                     theta_var_along_B
    !   4   logical                     lasym
    !   5   logical                     lrfp
    !   6   logical                     lfreeb
    !   7   integer                     max_it_NR
    !   8   integer                     max_it_r
    !   9   integer                     format_out
    !   10  integer                     style
    !   11  integer                     n_par
    !   12  integer                     n_r
    !   13  integer                     mpol
    !   14  integer                     ntor
    !   15  integer                     nfp
    !   16  integer                     EV_style
    !   17  integer                     n_procs_per_alpha
    !   18  integer                     n_alpha
    !   19  integer                     n_X
    !   20  integer                     min_m_X
    !   21  integer                     max_m_X
    !   22  integer                     n_sol_requested
    !   23  integer                     min_n_r_X
    !   24  integer                     min_r_eq
    !   25  integer                     max_r_eq
    !   26  real_dp                     min_alpha
    !   27  real_dp                     max_alpha
    !   28  real_dp                     min_r_X
    !   29  real_dp                     max_r_X
    !   30  real_dp                     tol_NR
    !   31  real_dp                     min_par
    !   32  real_dp                     max_par
    !   33  real_dp                     version
    !   34  real_dp(n_r)                phi(n_r) 
    !   35  real_dp(n_r)                phi_r(n_r) 
    !   36  real_dp(n_r)                iotaf(n_r) 
    !   37  real_dp(n_r)                presf(n_r) 
    !   38  real_dp(*)                  R_c(*)
    !   39  real_dp(*)                  R_s(*)
    !   40  real_dp(*)                  Z_c(*)
    !   41  real_dp(*)                  Z_s(*)
    !   42  real_dp(*)                  L_c(*)
    !   43  real_dp(*)                  L_s(*)
    !   with (*) = (0:mpol-1,-ntor:ntor,1:n_r,0:max_deriv(3))
    ! [MPI] Collective call
    integer function broadcast_vars() result(ierr)
        use VMEC_vars, only: mpol, ntor, lasym, lrfp, lfreeb, nfp, iotaf, &
            &R_c, R_s, Z_c, Z_s, L_c, L_s, phi, phi_r, presf, version, n_r
        use num_vars, only: max_str_ln, output_name, ltest, &
            &theta_var_along_B, EV_style, max_it_NR, max_it_r, n_alpha, &
            &n_procs_per_alpha, style, max_alpha, min_alpha, tol_NR, glob_rank, &
            &glob_n_procs, n_sol_requested, min_n_r_X, min_r_X, max_r_X
        use output_ops, only: format_out
        use X_vars, only: n_X, min_m_X, max_m_X
        use eq_vars, only: n_par, max_par, min_par, min_r_eq, max_r_eq
        
        character(*), parameter :: rout_name = 'broadcast_vars'
        
        ! initialize ierr
        ierr = 0
        
        if (glob_n_procs.gt.1) then                                             ! need to broadcast from global rank 0 to other processes
            ! print message
            call writo('Broadcasting variables determined by input')
            call lvl_ud(1)
            
            call MPI_Bcast(output_name,max_str_ln,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(ltest,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(theta_var_along_B,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(lasym,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(lrfp,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(lfreeb,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
            CHCKERR('MPI broadcast failed')
            call MPI_Bcast(max_it_NR,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(max_it_r,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(format_out,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(style,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(n_par,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(n_r,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(mpol,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(ntor,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(nfp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(EV_style,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(n_procs_per_alpha,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(n_alpha,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(n_X,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(min_m_X,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(max_m_X,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(n_sol_requested,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(min_n_r_X,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(min_r_eq,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(max_r_eq,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(min_alpha,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(max_alpha,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(min_r_X,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(max_r_X,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(tol_NR,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(min_par,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(max_par,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_Bcast(version,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call bcast_size_1_R(phi)
            call MPI_Bcast(phi,size(phi),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call bcast_size_1_R(phi_r)
            call MPI_Bcast(phi_r,size(phi_r),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call bcast_size_1_R(iotaf)
            call MPI_Bcast(iotaf,size(iotaf),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call bcast_size_1_R(presf)
            call MPI_Bcast(presf,size(presf),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call bcast_size_4_R(R_c)
            call MPI_Bcast(R_c,size(R_c),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call bcast_size_4_R(R_s)
            call MPI_Bcast(R_s,size(R_s),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call bcast_size_4_R(Z_c)
            call MPI_Bcast(Z_c,size(Z_c),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call bcast_size_4_R(Z_s)
            call MPI_Bcast(Z_s,size(Z_s),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call bcast_size_4_R(L_c)
            call MPI_Bcast(L_c,size(L_c),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call bcast_size_4_R(L_s)
            call MPI_Bcast(L_s,size(L_s),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            
            call lvl_ud(-1)
            call writo('Variables broadcasted')
        end if
    contains
        ! broadcasts the size of an array, so this array can be allocated in the
        ! slave processes
        ! The index of this array is (1:n_r)
        subroutine bcast_size_1_I(arr)                                          ! version with 1 integer argument
            ! input / output
            integer, intent(inout), allocatable :: arr(:)
            
            ! local variables
            integer :: arr_size                                                 ! sent ahead so arrays can be allocated
            
            if (glob_rank.eq.0) arr_size = size(arr)
            
            call MPI_Bcast(arr_size,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            if (glob_rank.ne.0) allocate(arr(1:arr_size))
        end subroutine bcast_size_1_I
        ! The index of this array is (1:n_r)
        subroutine bcast_size_1_R(arr)                                            ! version with 1 real argument
            ! input / output
            real(dp), intent(inout), allocatable :: arr(:)
            
            ! local variables
            integer :: arr_size                                                 ! sent ahead so arrays can be allocated
            
            if (glob_rank.eq.0) arr_size = size(arr)
            
            call MPI_Bcast(arr_size,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            if (glob_rank.ne.0) allocate(arr(1:arr_size))
        end subroutine bcast_size_1_R
        ! The index of this array is (0:mpol-1,-ntor:ntor,1:n_r,0:max_deriv(3))
        subroutine bcast_size_4_R(arr)                                          ! version with 4 real arguments
            ! input / output
            real(dp), intent(inout), allocatable :: arr(:,:,:,:)
            
            ! local variables
            integer :: arr_size(4)                                              ! sent ahead so arrays can be allocated
            
            if (glob_rank.eq.0) arr_size = [size(arr,1),size(arr,2),&
                &size(arr,3),size(arr,4)]
            
            call MPI_Bcast(arr_size,4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            if (glob_rank.ne.0) allocate(arr(0:arr_size(1)-1,&
                &-(arr_size(2)-1)/2:(arr_size(2)-1)/2,1:arr_size(3),&
                &0:arr_size(4)-1))
        end subroutine bcast_size_4_R
    end function broadcast_vars
end module
