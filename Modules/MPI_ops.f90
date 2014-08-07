!------------------------------------------------------------------------------!
!   Routines related to MPI                                                    !
!------------------------------------------------------------------------------!
module MPI_ops
#include <PB3D_macros.h>
    use MPI
    use num_vars, only: glob_rank, glob_n_procs, MPI_comm_groups
    use str_ops, only: i2str
    use output_ops, only: writo, lvl_ud
    
    implicit none
    private
    public start_MPI, stop_MPI, split_MPI, abort_MPI
    
contains
    ! start MPI and gather information
    subroutine start_MPI(ierr)
        character(*), parameter :: rout_name = 'start_MPI'
        
        ! input / output
        integer, intent(inout) :: ierr                                          ! error variable
        
        ! initialize ierr
        ierr = 0
        
        call MPI_init(ierr)                                                     ! initialize MPI
        CHCKERR('')
        call MPI_Comm_rank(MPI_COMM_WORLD,glob_rank,ierr)                       ! rank
        CHCKERR('')
        call MPI_Comm_size(MPI_COMM_WORLD,glob_n_procs,ierr)                    ! nr. processes
        CHCKERR('')
    end subroutine
    
    ! determine   how   to   split    the   communicator   MPI_COMM_WORLD   into
    ! subcommunicators,  according to  n_procs_per_alpha,  which determines  how
    ! many processors to use per field line
    subroutine split_MPI(ierr)
        use num_vars, only: n_procs_per_alpha, n_procs, n_alpha, max_str_ln, &
            &MPI_comm_groups, glob_rank, glob_n_procs
        
        character(*), parameter :: rout_name = 'split_MPI'
        
        ! input / output
        integer, intent(inout) :: ierr                                          ! error
        
        ! local variables
        integer :: n_groups                                                     ! how many subgroups are formed for alpha's
        integer :: remainder                                                    ! remainder of division
        integer :: id                                                           ! counter
        integer :: group                                                        ! group of local process, used to split MPI_COMM_WORLD
        integer :: sum_groups                                                   ! sum of previous colros
        logical :: group_found                                                  ! when searching for group of local process
        character(len=max_str_ln) :: str_procs, str_groups, str_min_procs       ! strings used in user messages
        
        !integer :: loc_group_size, loc_rank
        
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
            str_procs = '1 process'
        else
            str_procs = trim(i2str(glob_n_procs))//' processes'
        end if
        if (n_groups.eq.1) then
            str_groups = '1 group'
        else
            str_groups = trim(i2str(n_groups))//' groups'
        end if
        if (n_procs_per_alpha.eq.1) then
            str_min_procs = ''
        else
            str_min_procs = ' of at least '//trim(i2str(n_procs_per_alpha))//&
                &' processes'
        end if
        call writo(trim(str_procs)//' in total, for '//trim(str_groups)//&
            &trim(str_min_procs))                                               ! how many groups
        if (n_groups.eq.1) then
            str_groups = 'This group dynamically solves'
        else
            str_groups = 'These '//trim(i2str(n_groups))//' groups &
                &dynamically solve'
        end if
        if (n_alpha.eq.1) then
            str_procs = '1 field line'
        else
            str_procs = trim(i2str(n_alpha))//' field lines'
        end if
        call writo(trim(str_groups)//' a total of '//trim(str_procs)//&
            &' in parallel')                                                    ! how many field lines to solve
        
        ! determine group of local process
        group = 1                                                               ! start group on 1
        group_found = .false.
        sum_groups = 0
        do while (.not.group_found)
            sum_groups = sum_groups + n_procs(group)
            if (glob_rank+1.le.sum_groups) then                                 ! group found
                group_found = .true.
            else                                                                ! group not found: try next group
                group = group + 1
            end if
        end do
        
        ! split MPI_COMM_WORLD according to n_procs
        call MPI_Comm_split(MPI_COMM_WORLD,group,0,MPI_comm_groups,ierr)
        CHCKERR('')
    end subroutine
    
    ! stop MPI
    subroutine stop_MPI(ierr)
        character(*), parameter :: rout_name = 'stop_MPI'
        
        ! input / output
        integer, intent(inout) :: ierr                                          ! error variable
        
        call MPI_finalize(ierr)
        CHCKERR('')
    end subroutine
    
    ! abort MPI
    subroutine abort_MPI(ierr)
        character(*), parameter :: rout_name = 'abort_MPI'
        
        ! input / output
        integer, intent(inout) :: ierr                                          ! error variable
        
        call MPI_Abort(MPI_COMM_WORLD,0,ierr)
        CHCKERR('')
    end subroutine
end module
