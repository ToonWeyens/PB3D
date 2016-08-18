!------------------------------------------------------------------------------!
!   Numerical utilities related to MPI                                         !
!------------------------------------------------------------------------------!
module MPI_utilities
#include <PB3D_macros.h>
    use str_ops
    use messages
    use MPI
    use num_vars, only: dp, max_str_ln, pi
    
    implicit none
    private
    public get_ser_var, get_ghost_arr, broadcast_var, wait_MPI, mutex_req_acc, &
        &mutex_return_acc
#if ldebug
    public debug_mutex
#endif
    
#if ldebug
    ! global variables
    logical :: debug_mutex = .false.                                            ! print debug information about mutex operations
#endif
    
    ! interfaces
    interface get_ser_var
        module procedure get_ser_var_complex, get_ser_var_real, get_ser_var_int
    end interface
    interface get_ghost_arr
        module procedure get_ghost_arr_3D_complex, get_ghost_arr_3D_real, &
            &get_ghost_arr_2D_complex, get_ghost_arr_1D_real
    end interface
    interface broadcast_var
        module procedure broadcast_var_real, broadcast_var_int, &
            &broadcast_var_log
    end interface
    
contains
    ! Gather parallel variable in serial version on group master or, optionally,
    ! to all the processes, using the variable "scatter"
    ! Note:  the serial variable  has to be  allocatable and if  unallocated, it
    ! will be allocated.
    ! [MPI] Collective call
    integer function get_ser_var_complex(var,ser_var,scatter) result(ierr)      ! complex version
        use num_vars, only: rank, n_procs
        
        character(*), parameter :: rout_name = 'get_ser_var'
        
        ! input / output
        complex(dp), intent(in) :: var(:)                                       ! parallel vector
        complex(dp), allocatable, intent(inout) :: ser_var(:)                   ! serial vector
        logical, intent(in), optional :: scatter                                ! optionally scatter the result to all the processes
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer, allocatable :: recvcounts(:)                                   ! counts of nr. of elements received from each processor
        integer, allocatable :: displs(:)                                       ! displacements elements received from each processor
        integer :: id                                                           ! counter
        logical :: scatter_loc                                                  ! local copy of scatter
        
        ! initialize ierr
        ierr = 0
        
        ! set local copy of scatter
        scatter_loc = .false.
        if (present(scatter)) scatter_loc = scatter
        
        ! gather local size  of var of all groups onto  main processor, to serve
        ! as receive counts on group master
        if (rank.eq.0 .or. scatter_loc) then
            allocate(recvcounts(n_procs))
            allocate(displs(n_procs))
        else
            allocate(recvcounts(0))
            allocate(displs(0))
        end if
        if (scatter_loc) then
            call MPI_Allgather(size(var),1,MPI_INTEGER,recvcounts,1,&
                &MPI_INTEGER,MPI_Comm_world,ierr)
        else
            call MPI_Gather(size(var),1,MPI_INTEGER,recvcounts,1,&
                &MPI_INTEGER,0,MPI_Comm_world,ierr)
        end if
        err_msg = 'Failed to gather size of parallel variable'
        CHCKERR(err_msg)
        
        ! allocate serial variable
        if (allocated(ser_var)) then
            if (size(ser_var).ne.sum(recvcounts)) then
                ierr = 1
                err_msg = 'ser_var has wrong dimensions'
                CHCKERR(err_msg)
            end if
        else
            allocate(ser_var(sum(recvcounts)))
        end if
        
        ! deduce displacements by summing recvcounts
        if (rank.eq.0 .or. scatter_loc) then
            displs(1) = 0
            do id = 2,n_procs
                displs(id) = displs(id-1) + recvcounts(id-1)
            end do
        end if
        
        if (scatter_loc) then
            call MPI_Allgatherv(var,size(var),MPI_DOUBLE_COMPLEX,ser_var,&
                &recvcounts,displs,MPI_DOUBLE_COMPLEX,MPI_Comm_world,ierr)
        else
            call MPI_Gatherv(var,size(var),MPI_DOUBLE_COMPLEX,ser_var,&
                &recvcounts,displs,MPI_DOUBLE_COMPLEX,0,MPI_Comm_world,ierr)
        end if
        err_msg = 'Failed to gather parallel variable'
        CHCKERR(err_msg)
    end function get_ser_var_complex
    integer function get_ser_var_real(var,ser_var,scatter) result(ierr)         ! real version
        use num_vars, only: rank, n_procs
        
        character(*), parameter :: rout_name = 'get_ser_var_real'
        
        ! input / output
        real(dp), intent(in) :: var(:)                                          ! parallel vector
        real(dp), allocatable, intent(inout) :: ser_var(:)                      ! serial vector
        logical, intent(in), optional :: scatter                                ! optionally scatter the result to all the processes
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer, allocatable :: recvcounts(:)                                   ! counts of nr. of elements received from each processor
        integer, allocatable :: displs(:)                                       ! displacements elements received from each processor
        integer :: id                                                           ! counter
        logical :: scatter_loc                                                  ! local copy of scatter
        
        ! initialize ierr
        ierr = 0
        
        ! set local copy of scatter
        scatter_loc = .false.
        if (present(scatter)) scatter_loc = scatter
        
        ! gather local size  of var of all groups onto  main processor, to serve
        ! as receive counts on group master
        if (rank.eq.0 .or. scatter_loc) then
            allocate(recvcounts(n_procs))
            allocate(displs(n_procs))
        else
            allocate(recvcounts(0))
            allocate(displs(0))
        end if
        if (scatter_loc) then
            call MPI_Allgather(size(var),1,MPI_INTEGER,recvcounts,1,&
                &MPI_INTEGER,MPI_Comm_world,ierr)
        else
            call MPI_Gather(size(var),1,MPI_INTEGER,recvcounts,1,&
                &MPI_INTEGER,0,MPI_Comm_world,ierr)
        end if
        err_msg = 'Failed to gather size of parallel variable'
        CHCKERR(err_msg)
        
        ! allocate serial variable
        if (allocated(ser_var)) then
            if (size(ser_var).ne.sum(recvcounts)) then
                ierr = 1
                err_msg = 'ser_var has wrong dimensions'
                CHCKERR(err_msg)
            end if
        else
            allocate(ser_var(sum(recvcounts)))
        end if
        
        ! deduce displacements by summing recvcounts
        if (rank.eq.0 .or. scatter_loc) then
            displs(1) = 0
            do id = 2,n_procs
                displs(id) = displs(id-1) + recvcounts(id-1)
            end do
        end if
        
        if (scatter_loc) then
            call MPI_Allgatherv(var,size(var),MPI_DOUBLE_PRECISION,ser_var,&
                &recvcounts,displs,MPI_DOUBLE_PRECISION,MPI_Comm_world,ierr)
        else
            call MPI_Gatherv(var,size(var),MPI_DOUBLE_PRECISION,ser_var,&
                &recvcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_Comm_world,ierr)
        end if
        err_msg = 'Failed to gather parallel variable'
        CHCKERR(err_msg)
    end function get_ser_var_real
    integer function get_ser_var_int(var,ser_var,scatter) result(ierr)          ! integer version
        use num_vars, only: rank, n_procs
        
        character(*), parameter :: rout_name = 'get_ser_var_int'
        
        ! input / output
        integer, intent(in) :: var(:)                                           ! parallel vector
        integer, allocatable, intent(inout) :: ser_var(:)                       ! serial vector
        logical, intent(in), optional :: scatter                                ! optionally scatter the result to all the processes
        
        ! local variables
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer, allocatable :: recvcounts(:)                                   ! counts of nr. of elements received from each processor
        integer, allocatable :: displs(:)                                       ! displacements elements received from each processor
        integer :: id                                                           ! counter
        logical :: scatter_loc                                                  ! local copy of scatter
        
        ! initialize ierr
        ierr = 0
        
        ! set local copy of scatter
        scatter_loc = .false.
        if (present(scatter)) scatter_loc = scatter
        
        ! gather local size  of var of all groups onto  main processor, to serve
        ! as receive counts on group master
        if (rank.eq.0 .or. scatter_loc) then
            allocate(recvcounts(n_procs))
            allocate(displs(n_procs))
        else
            allocate(recvcounts(0))
            allocate(displs(0))
        end if
        if (scatter_loc) then
            call MPI_Allgather(size(var),1,MPI_INTEGER,recvcounts,1,&
                &MPI_INTEGER,MPI_Comm_world,ierr)
        else
            call MPI_Gather(size(var),1,MPI_INTEGER,recvcounts,1,&
                &MPI_INTEGER,0,MPI_Comm_world,ierr)
        end if
        err_msg = 'Failed to gather size of parallel variable'
        CHCKERR(err_msg)
        
        ! allocate serial variable
        if (allocated(ser_var)) then
            if (size(ser_var).ne.sum(recvcounts)) then
                ierr = 1
                err_msg = 'ser_var has wrong dimensions'
                CHCKERR(err_msg)
            end if
        else
            allocate(ser_var(sum(recvcounts)))
        end if
        
        ! deduce displacements by summing recvcounts
        if (rank.eq.0 .or. scatter_loc) then
            displs(1) = 0
            do id = 2,n_procs
                displs(id) = displs(id-1) + recvcounts(id-1)
            end do
        end if
        
        if (scatter_loc) then
            call MPI_Allgatherv(var,size(var),MPI_INTEGER,ser_var,&
                &recvcounts,displs,MPI_INTEGER,MPI_Comm_world,ierr)
        else
            call MPI_Gatherv(var,size(var),MPI_INTEGER,ser_var,&
                &recvcounts,displs,MPI_INTEGER,0,MPI_Comm_world,ierr)
        end if
        err_msg = 'Failed to gather parallel variable'
        CHCKERR(err_msg)
    end function get_ser_var_int
    
    ! Fill the ghost regions in an array  by sending the first normal point of a
    ! process to  the left process. Every  message is identified by  its sending
    ! process. The array should have the extended size, including ghost regions.
    ! [MPI] Collective call
    integer function get_ghost_arr_3D_complex(arr,size_ghost) result(ierr)      ! 3D complex version
        use num_vars, only: rank, n_procs
        
        character(*), parameter :: rout_name = 'get_ghost_arr_3D_complex'
        
        ! input / output
        complex(dp), intent(inout) :: arr(:,:,:)                                ! divided array
        integer, intent(in) :: size_ghost                                       ! width of ghost region
        
        ! local variables
        integer :: n_modes(2)                                                   ! number of modes
        integer :: tot_size                                                     ! total size (including ghost region)
        integer :: istat(MPI_STATUS_SIZE)                                       ! status of send-receive
        
        ! initialize ierr
        ierr = 0
        
        ! initialize number of modes and total size
        n_modes = [size(arr,1),size(arr,2)]
        tot_size = size(arr,3)
        
        ! ghost regions only make sense if there is more than 1 process
        if (n_procs.gt.1) then
            if (rank.eq.0) then                                                 ! first rank only receives
                call MPI_Recv(arr(:,:,tot_size-size_ghost+1:tot_size),&
                    &size_ghost*product(n_modes),MPI_DOUBLE_COMPLEX,rank+1,&
                    &rank+1,MPI_Comm_world,istat,ierr)
                CHCKERR('Failed to receive')
            else if (rank+1.eq.n_procs) then                                    ! last rank only sends
                call MPI_Send(arr(:,:,1:size_ghost),&
                    &size_ghost*product(n_modes),MPI_DOUBLE_COMPLEX,rank-1,&
                    &rank,MPI_Comm_world,ierr)
                CHCKERR('Failed to send')
            else                                                                ! middle ranks send and receive
                call MPI_Sendrecv(arr(:,:,1:size_ghost),&
                    &size_ghost*product(n_modes),MPI_DOUBLE_COMPLEX,rank-1,&
                    &rank,arr(:,:,tot_size-size_ghost+1:tot_size),&
                    &size_ghost*product(n_modes),MPI_DOUBLE_COMPLEX,&
                    &rank+1,rank+1,MPI_Comm_world,istat,ierr)
                CHCKERR('Failed to send and receive')
            end if
        end if
    end function get_ghost_arr_3D_complex
    integer function get_ghost_arr_3D_real(arr,size_ghost) result(ierr)         ! 3D real version
        use num_vars, only: rank, n_procs
        
        character(*), parameter :: rout_name = 'get_ghost_arr_3D_real'
        
        ! input / output
        real(dp), intent(inout) :: arr(:,:,:)                                   ! divided array
        integer, intent(in) :: size_ghost                                       ! width of ghost region
        
        ! local variables
        integer :: n_modes(2)                                                   ! number of modes
        integer :: tot_size                                                     ! total size (including ghost region)
        integer :: istat(MPI_STATUS_SIZE)                                       ! status of send-receive
        
        ! initialize ierr
        ierr = 0
        
        ! initialize number of modes and total size
        n_modes = [size(arr,1),size(arr,2)]
        tot_size = size(arr,3)
        
        ! ghost regions only make sense if there is more than 1 process
        if (n_procs.gt.1) then
            if (rank.eq.0) then                                                 ! first rank only receives
                call MPI_Recv(arr(:,:,tot_size-size_ghost+1:tot_size),&
                    &size_ghost*product(n_modes),MPI_DOUBLE_PRECISION,&
                    &rank+1,rank+1,MPI_Comm_world,istat,ierr)
                CHCKERR('Failed to receive')
            else if (rank+1.eq.n_procs) then                                    ! last rank only sends
                call MPI_Send(arr(:,:,1:size_ghost),&
                    &size_ghost*product(n_modes),MPI_DOUBLE_PRECISION,&
                    &rank-1,rank,MPI_Comm_world,ierr)
                CHCKERR('Failed to send')
            else                                                                ! middle ranks send and receive
                call MPI_Sendrecv(arr(:,:,1:size_ghost),&
                    &size_ghost*product(n_modes),MPI_DOUBLE_PRECISION,&
                    &rank-1,&
                    &rank,arr(:,:,tot_size-size_ghost+1:tot_size),&
                    &size_ghost*product(n_modes),MPI_DOUBLE_PRECISION,&
                    &rank+1,rank+1,MPI_Comm_world,istat,ierr)
                CHCKERR('Failed to send and receive')
            end if
        end if
    end function get_ghost_arr_3D_real
    integer function get_ghost_arr_2D_complex(arr,size_ghost) result(ierr)      ! 2D complex version
        use num_vars, only: rank, n_procs
        
        character(*), parameter :: rout_name = 'get_ghost_arr_2D_complex'
        
        ! input / output
        complex(dp), intent(inout) :: arr(:,:)                                  ! divided array
        integer, intent(in) :: size_ghost                                       ! width of ghost region
        
        ! local variables
        integer :: n_modes                                                      ! number of modes
        integer :: tot_size                                                     ! total size (including ghost region)
        integer :: istat(MPI_STATUS_SIZE)                                       ! status of send-receive
        
        ! initialize ierr
        ierr = 0
        
        ! initialize number of modes and total size
        n_modes = size(arr,1)
        tot_size = size(arr,2)
        
        ! ghost regions only make sense if there is more than 1 process
        if (n_procs.gt.1) then
            if (rank.eq.0) then                                                 ! first rank only receives
                call MPI_Recv(arr(:,tot_size-size_ghost+1:tot_size),&
                    &size_ghost*n_modes,MPI_DOUBLE_COMPLEX,rank+1,&
                    &rank+1,MPI_Comm_world,istat,ierr)
                CHCKERR('Failed to receive')
            else if (rank+1.eq.n_procs) then                                    ! last rank only sends
                call MPI_Send(arr(:,1:size_ghost),size_ghost*n_modes,&
                    &MPI_DOUBLE_COMPLEX,rank-1,rank,MPI_Comm_world,&
                    &ierr)
                CHCKERR('Failed to send')
            else                                                                ! middle ranks send and receive
                call MPI_Sendrecv(arr(:,1:size_ghost),size_ghost*n_modes,&
                    &MPI_DOUBLE_COMPLEX,rank-1,rank,&
                    &arr(:,tot_size-size_ghost+1:tot_size),size_ghost*n_modes,&
                    &MPI_DOUBLE_COMPLEX,rank+1,rank+1,MPI_Comm_world,&
                    &istat,ierr)
                CHCKERR('Failed to send and receive')
            end if
        end if
    end function get_ghost_arr_2D_complex
    integer function get_ghost_arr_1D_real(arr,size_ghost) result(ierr)         ! 1D real version
        use num_vars, only: rank, n_procs
        
        character(*), parameter :: rout_name = 'get_ghost_arr_1D_real'
        
        ! input / output
        real(dp), intent(in) :: arr(:)                                          ! divided array
        integer, intent(in) :: size_ghost                                       ! width of ghost region
        
        ! local variables
        integer :: tot_size                                                     ! total size (including ghost region)
        integer :: istat(MPI_STATUS_SIZE)                                       ! status of send-receive
        
        ! initialize ierr
        ierr = 0
        
        ! initialize number of modes and total size
        tot_size = size(arr)
        
        ! ghost regions only make sense if there is more than 1 process
        if (n_procs.gt.1) then
            if (rank.eq.0) then                                                 ! first rank only receives
                call MPI_Recv(arr(tot_size-size_ghost+1:tot_size),&
                    &size_ghost,MPI_DOUBLE_PRECISION,rank+1,&
                    &rank+1,MPI_Comm_world,istat,ierr)
                CHCKERR('Failed to receive')
            else if (rank+1.eq.n_procs) then                                    ! last rank only sends
                call MPI_Send(arr(1:size_ghost),size_ghost,&
                    &MPI_DOUBLE_PRECISION,rank-1,rank,MPI_Comm_world,&
                    &ierr)
                CHCKERR('Failed to send')
            else                                                                ! middle ranks send and receive
                call MPI_Sendrecv(arr(1:size_ghost),size_ghost,&
                    &MPI_DOUBLE_PRECISION,rank-1,rank,&
                    &arr(tot_size-size_ghost+1:tot_size),size_ghost,&
                    &MPI_DOUBLE_PRECISION,rank+1,rank+1,&
                    &MPI_Comm_world,istat,ierr)
                CHCKERR('Failed to send and receive')
            end if
        end if
    end function get_ghost_arr_1D_real
    
    ! wrapper function to broadcast a single variable
    integer function broadcast_var_real(var,source) result(ierr)                ! version for reals
        character(*), parameter :: rout_name = 'broadcast_var_real'
        
        ! input / output
        real(dp), intent(in) :: var                                             ! variable to be broadcast
        integer, intent(in), optional :: source                                 ! process that sends
        
        ! local variables
        integer :: source_loc = 0                                               ! local value for source
        
        ! initialize ierr
        ierr = 0
        
        ! set local source if given
        if (present(source)) source_loc = source
        
        call MPI_Bcast(var,1,MPI_DOUBLE_PRECISION,source_loc,MPI_COMM_WORLD,&
            &ierr)
        CHCKERR('MPI broadcast failed')
    end function broadcast_var_real
    integer function broadcast_var_int(var,source) result(ierr)                 ! version for integers
        character(*), parameter :: rout_name = 'broadcast_var_int'
        
        ! input / output
        integer, intent(in) :: var                                              ! variable to be broadcast
        integer, intent(in), optional :: source                                 ! process that sends
        
        ! local variables
        integer :: source_loc = 0                                               ! local value for source
        
        ! initialize ierr
        ierr = 0
        
        ! set local source if given
        if (present(source)) source_loc = source
        
        call MPI_Bcast(var,1,MPI_INTEGER,source_loc,MPI_COMM_WORLD,ierr)
        CHCKERR('MPI broadcast failed')
    end function broadcast_var_int
    integer function broadcast_var_log(var,source) result(ierr)                 ! version for logicals
        character(*), parameter :: rout_name = 'broadcast_var_log'
        
        ! input / output
        logical, intent(in) :: var                                              ! variable to be broadcast
        integer, intent(in), optional :: source                                 ! process that sends
        
        ! local variables
        integer :: source_loc = 0                                               ! local value for source
        
        ! initialize ierr
        ierr = 0
        
        ! set local source if given
        if (present(source)) source_loc = source
        
        call MPI_Bcast(var,1,MPI_LOGICAL,source_loc,MPI_COMM_WORLD,ierr)
        CHCKERR('MPI broadcast failed')
    end function broadcast_var_log
    
    ! MPI wait
    integer function wait_MPI() result(ierr)
        character(*), parameter :: rout_name = 'wait_MPI'
        
        ! initialize ierr
        ierr = 0
        
        ! set barrier
        call MPI_Barrier(MPI_Comm_world,ierr)
        CHCKERR('MPI Barrier failed')
    end function wait_MPI
    
    ! Request access to mutex.
    ! (based on http://www.mcs.anl.gov/~thakur/papers/atomic-mode.pdf)
    integer function mutex_req_acc(mutex_win,wakeup_tag) result(ierr)
#if ldebug
        use num_vars, only: rank
#endif
        
        character(*), parameter :: rout_name = 'mutex_req_acc'
        
        ! input / output
        integer, intent(in) :: mutex_win                                        ! window to mutex
        integer, intent(in) :: wakeup_tag                                       ! wakeup tag
        
        ! local variables
        integer, allocatable :: mutex_loc(:)                                    ! local copy of waiting to list
        integer :: dum_buf                                                      ! dummy buffer
        
        ! initialize ierr
        ierr = 0
        
        ! add current process to mutex waiting list
        ierr = mutex_wlist_add_remove(.true.,mutex_win,mutex_loc)
        
        ! wait for notification if mutex already held
        if (sum(mutex_loc).gt.0) then
#if ldebug
            if (debug_mutex) &
                &write(*,*) '>>>', rank, 'waiting for mutex on', mutex_loc
#endif
            call MPI_Recv(dum_buf,1,MPI_INTEGER,MPI_ANY_SOURCE,&
                &wakeup_tag,MPI_Comm_world,MPI_STATUS_IGNORE,ierr)
            CHCKERR('Failed to receive notification')
        end if
#if ldebug
            if (debug_mutex) write(*,*) '>>>', rank, 'got mutex'
#endif
    end function mutex_req_acc
    
    ! Returns access to a mutex.
    ! (based on http://www.mcs.anl.gov/~thakur/papers/atomic-mode.pdf)
    integer function mutex_return_acc(mutex_win,wakeup_tag) result(ierr)
        use num_vars, only: rank, n_procs
        
        character(*), parameter :: rout_name = 'mutex_return_acc'
        
        ! input / output
        integer, intent(in) :: mutex_win                                        ! window to mutex
        integer, intent(in) :: wakeup_tag                                       ! wakeup tag
        
        ! local variables
        integer, allocatable :: mutex_loc(:)                                    ! local copy of waiting list to mutex
        integer :: id                                                           ! counter
        integer :: next_rank                                                    ! next rank to be running
        integer :: dum_buf                                                      ! dummy buffer
        
        ! initialize ierr
        ierr = 0
        
#if ldebug
        if(debug_mutex) write(*,*) '>>>', rank, 'returning access'
#endif
        ! remove current process to mutex waiting list
        ierr = mutex_wlist_add_remove(.false.,mutex_win,mutex_loc)
        CHCKERR('')
        
        ! give mutex to next proces that requests it
        notify_next: if (sum(mutex_loc).gt.0) then
            ! find the next rank, wrapping around zero
            do id = 1,n_procs-1
                next_rank = mod(rank+id-1,n_procs-1)                            ! wrap around n_procs-1 to zero
                if (mutex_loc(next_rank+1).eq.1) then                           ! arrays start at 1, ranks at 0
                    dum_buf = 1
                    ! next rank possibly needs to be incremented for full arrays
                    if (next_rank.ge.rank) next_rank = next_rank+1
                    ! send to next rank
                    call MPI_Send(dum_buf,1,MPI_INTEGER,next_rank,&
                        &wakeup_tag,MPI_Comm_world,ierr)
                    CHCKERR('Failed to send notification')
#if ldebug
                    if (debug_mutex) &
                        &write(*,*) '>>>', rank, 'notified', next_rank
#endif
                    exit notify_next
                end if
            end do
            ! no next process found: impossible situation
            ierr = 1
            CHCKERR('No next process found')
        end if notify_next
    end function mutex_return_acc
    
    ! Adds or removes a rank from the waiting list for a mutex and returns the
    ! previous mutex waiting list.
    ! (based on http://www.mcs.anl.gov/~thakur/papers/atomic-mode.pdf)
    integer function mutex_wlist_add_remove(add,mutex_win,mutex_excl) &
        &result(ierr)
        use num_vars, only: n_procs, rank
        
        character(*), parameter :: rout_name = 'mutex_wlist_add_remove'
        
        ! input / output
        logical, intent(in) :: add                                              ! true if add, false if remove from mutex waiting list
        integer, intent(in) :: mutex_win                                        ! window to mutex
        integer, intent(inout), allocatable :: mutex_excl(:)                    ! previous mutex waiting list, excluding local proc.
        
        ! local variables
        integer :: lens(2), disps(2)                                            ! lengths and displacements of excluded type
        integer :: mutex_type                                                   ! indexed type to access mutex, excluding the actual process
        integer :: onezero                                                      ! one or zero
        integer(kind=MPI_ADDRESS_KIND) :: one = 1                               ! one
        
        ! initialize ierr
        ierr = 0
        
        ! set plus or minus one and allocate mutex_excl
        if (add) then
            onezero = 1
        else
            onezero = 0
        end if
        allocate(mutex_excl(n_procs-1))
        
        ! set lengths and displacements of extended type to exclude current rank
        lens = [rank, n_procs-rank-1]
        disps = [0, rank+1]
        
        ! create type
        call MPI_Type_indexed(2,lens,disps,MPI_INTEGER,mutex_type,ierr)
        CHCKERR('Failed to create indexed type')
        call MPI_Type_commit(mutex_type,ierr)
        CHCKERR('Failed to commit type')
        
        ! get window to mutex
        call MPI_Win_lock(MPI_LOCK_EXCLUSIVE,0,0,mutex_win,ierr)
        CHCKERR('Failed to lock window')
        
        ! get previous list and put value for current proc
        call MPI_Get(mutex_excl,n_procs-1,MPI_INTEGER,0,one*0,1,&
            &mutex_type,mutex_win,ierr)
        CHCKERR('Failed to get waiting list')
        call MPI_Put(onezero,1,MPI_INTEGER,0,rank*one,1,MPI_INTEGER,mutex_win,&
            &ierr)
        CHCKERR('Failed to add to waiting list')
        
        ! unlock window
        call MPI_Win_unlock(0,mutex_win,ierr)
        CHCKERR('Failed to lock window')
        
        ! free mutex type
        call MPI_Type_free(mutex_type,ierr)
        CHCKERR('Failed to free type')
    end function mutex_wlist_add_remove
end module MPI_utilities
