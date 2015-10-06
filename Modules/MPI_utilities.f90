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
    public get_ser_var, get_ghost_arr, broadcast_var, wait_MPI
    
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
end module MPI_utilities
