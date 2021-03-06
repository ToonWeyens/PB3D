!------------------------------------------------------------------------------!
!> Numerical utilities related to MPI.
!!
!! This  includes  a  lock  system,  which  can be  both  BL  (blocking)  or  NB
!! (non-blocking). It  is based on the  implementation of an MPI-IO  atomic mode
!! without file support, described in \cite RossAtomicIO.
!!
!! \see See mpi_vars.
!!
!! The reason  for this was the fact  that using a simple lock file  can lead to
!! crashes.
!!
!! \note A  downside of this method  is that in  some rare cases a  deadlock may
!! occur as the master process, which contains the shared variable with a window
!! that other  processes may use,  is idle and  waiting, whereas the  others are
!! still performing  lock options. As  the master is  idle and waiting,  its MPI
!! asynchronous  communication  is not  performed.  To  remedy this,  just  call
!! wait_MPI() after procedures where lock operations are performed.
!------------------------------------------------------------------------------!
module MPI_utilities
#include <PB3D_macros.h>
    use str_utilities
    use messages
    use MPI
    use num_vars, only: dp, max_str_ln, pi
    use MPI_vars, only: lock_type
    
    implicit none
    private
    public get_ser_var, redistribute_var, get_ghost_arr, broadcast_var, &
        &wait_MPI, lock_req_acc, lock_return_acc
#if ldebug
    public lock_header, lock_wl_change, &
        &debug_lock, n_waits
#endif
    
#if ldebug
    ! global variables
    ! Note: use "sort -nk1 [output file]" to get output sorted by time.
    logical :: debug_lock = .false.                                             !< print debug information about lock operations \ldebug
    integer :: n_waits = 0                                                      !< number of waits \ldebug
#endif
    
    ! interfaces
    
    !> \public Gather parallel variable in serial version on group master.
    !!
    !! Optionally,  all the  processes receive  the parallel  variable using  \c
    !! scatter.
    !!
    !! \note The  serial variable has to  be allocatable and if  unallocated, it
    !! will be allocated.
    !!
    !! \return ierr
    interface get_ser_var
        !> \public
        module procedure get_ser_var_complex
        !> \public
        module procedure get_ser_var_real
        !> \public
        module procedure get_ser_var_int
    end interface
    
    !> \public Fill the ghost regions in an array.
    !!
    !! This is done by  sending the first normal point of a  process to the left
    !! process.
    !!
    !! Every MPI message is identified by  its sending process. The array should
    !! have the extended size, including ghost regions.
    !!
    !! \return ierr
    interface get_ghost_arr
        !> \public
        module procedure get_ghost_arr_3D_complex
        !> \public
        module procedure get_ghost_arr_3D_real
        !> \public
        module procedure get_ghost_arr_2D_complex
        !> \public
        module procedure get_ghost_arr_1D_real
    end interface
    
    !> \public Wrapper function to broadcast a single variable using MPI.
    !!
    !! \return ierr
    interface broadcast_var
        !> \public
        module procedure broadcast_var_real
        !> \public
        module procedure broadcast_var_int
        !> \public
        module procedure broadcast_var_log
        !> \public
        module procedure broadcast_var_complex_arr
        !> \public
        module procedure broadcast_var_real_arr
        !> \public
        module procedure broadcast_var_int_arr
        !> \public
        module procedure broadcast_var_log_arr
    end interface
    
contains
    !> \private complex version
    integer function get_ser_var_complex(var,ser_var,scatter) result(ierr)
        use num_vars, only: rank, n_procs
        
        character(*), parameter :: rout_name = 'get_ser_var'
        
        ! input / output
        complex(dp), intent(in) :: var(:)                                       !< parallel vector
        complex(dp), allocatable, intent(inout) :: ser_var(:)                   !< serial vector
        logical, intent(in), optional :: scatter                                !< optionally scatter the result to all the processes
        
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
                write(*,*) size(ser_var), sum(recvcounts)
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
    !> \private real version
    integer function get_ser_var_real(var,ser_var,scatter) result(ierr)
        use num_vars, only: rank, n_procs
        
        character(*), parameter :: rout_name = 'get_ser_var_real'
        
        ! input / output
        real(dp), intent(in) :: var(:)                                          !< parallel vector
        real(dp), allocatable, intent(inout) :: ser_var(:)                      !< serial vector
        logical, intent(in), optional :: scatter                                !< optionally scatter the result to all the processes
        
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
                write(*,*) size(ser_var), sum(recvcounts)
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
    !> \private integer version
    integer function get_ser_var_int(var,ser_var,scatter) result(ierr)
        use num_vars, only: rank, n_procs
        
        character(*), parameter :: rout_name = 'get_ser_var_int'
        
        ! input / output
        integer, intent(in) :: var(:)                                           !< parallel vector
        integer, allocatable, intent(inout) :: ser_var(:)                       !< serial vector
        logical, intent(in), optional :: scatter                                !< optionally scatter the result to all the processes
        
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
                write(*,*) size(ser_var), sum(recvcounts)
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
    
    !> Redistribute variables according to new limits.
    integer function redistribute_var(var,dis_var,lims,lims_dis) result(ierr)
        use num_vars, only: n_procs, rank
        
        character(*), parameter :: rout_name = 'redistribute_var'
        
        ! input / output
        real(dp), intent(in) :: var(:)                                          !< parallel vector
        real(dp), intent(inout) :: dis_var(:)                                   !< redistributed vector
        integer, intent(in) :: lims(2)                                          !< indices of parallel vector
        integer, intent(in) :: lims_dis(2)                                      !< indices of redistributed parallel vector
        
        ! local variables
        integer :: id, jd                                                       ! counters
        integer, allocatable :: lims_tot(:,:)                                   ! limits of all processes
        integer, allocatable :: lims_dis_tot(:,:)                               ! redistributed limits of all processes
        integer, allocatable :: temp_lim(:)                                     ! temporary variable
        integer, allocatable :: n_vars(:,:)                                     ! number of values to be received and sent from each process (row: to, column: from)
        integer, allocatable :: nr_sen(:)                                       ! number of values to be sent to each process from this process
        integer, allocatable :: nr_rec(:)                                       ! number of values to be received from each process to this process
        integer, allocatable :: id_rec(:)                                       ! starting position of variable to be received in local memory (starting at 0)
        integer, allocatable :: id_sen(:)                                       ! starting position of variable to be sent in local memory (starting at 0)
        integer :: max_loc                                                      ! local maximum
        
        ! initialize ierr
        ierr = 0
        
        ! get limits and redistributed limits of all procs
        allocate(lims_tot(2,n_procs))
        allocate(lims_dis_tot(2,n_procs))
        do id = 1,2
            ierr = get_ser_var(lims(id:id),temp_lim,scatter=.true.)
            CHCKERR('')
            lims_tot(id,:) = temp_lim
            ierr = get_ser_var(lims_dis(id:id),temp_lim,scatter=.true.)
            CHCKERR('')
            lims_dis_tot(id,:) = temp_lim
        end do
        
        ! find  out how many  variables to send and  receive. A row  indicates a
        ! receive and a column a send. E.g.:
        !   (a b)
        !   (c d)
        ! means that  process 0  gets a  from process  0 and  b from  process 1.
        ! Process 1 gets c from process 0 and d from process 1.
        allocate(n_vars(n_procs,n_procs))
        n_vars = 0
        do jd = 1,n_procs
            max_loc = lims_dis_tot(1,jd)
            do id = 1,n_procs
                if (lims_tot(2,id).ge.max_loc) then
                    n_vars(jd,id) = min(lims_tot(2,id),lims_dis_tot(2,jd)) - &
                        &max_loc + 1
                    max_loc = max_loc + n_vars(jd,id)
                end if
            end do
        end do
        
        ! set up how many are to be sent by each process and which memory index
        !   - nr_rec is the row of n_vars corresponding to this process,
        !   - nr_sen is the column of n_vars corresponding to this process,
        !   -  id_rec  counts  the  number of  variables  gotten  from  previous
        !   processes,
        !   - id_sen counts the number  of variables sent to previous processes,
        !   taking into account  the difference betweeen the lower  limit of the
        !   process sending to and the current process.
        allocate(nr_rec(n_procs))
        allocate(nr_sen(n_procs))
        allocate(id_rec(n_procs))
        allocate(id_sen(n_procs))
        do id = 1,n_procs
            nr_rec(id) = n_vars(rank+1,id)
            nr_sen(id) = n_vars(id,rank+1)
            id_rec(id) = sum(n_vars(rank+1,1:id-1))
            id_sen(id) = sum(n_vars(id,1:rank)) + lims_dis_tot(1,id) - lims(1)
        end do
        
        call MPI_Alltoallv(var,nr_sen,id_sen,MPI_DOUBLE_PRECISION,&
            &dis_var,nr_rec,id_rec,MPI_DOUBLE_PRECISION,MPI_Comm_world,ierr)
        CHCKERR('')
        
        !call sleep(rank)
        !write(*,*) 'rank', rank
        !write(*,*) '    nr_rec', nr_rec
        !write(*,*) '    id_rec', id_rec
        !write(*,*) '    nr_sen', nr_sen
        !write(*,*) '    id_sen', id_sen
        !do id = 1,n_procs
            !write(*,*) '    ', n_vars(id,:)
        !end do
    end function redistribute_var
    
    !> \private 3-D complex version
    integer function get_ghost_arr_3D_complex(arr,size_ghost) result(ierr)
        use num_vars, only: rank, n_procs
        
        character(*), parameter :: rout_name = 'get_ghost_arr_3D_complex'
        
        ! input / output
        complex(dp), intent(inout) :: arr(:,:,:)                                !< divided array
        integer, intent(in) :: size_ghost                                       !< width of ghost region
        
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
    !> \private 3-D real version
    integer function get_ghost_arr_3D_real(arr,size_ghost) result(ierr)
        use num_vars, only: rank, n_procs
        
        character(*), parameter :: rout_name = 'get_ghost_arr_3D_real'
        
        ! input / output
        real(dp), intent(inout) :: arr(:,:,:)                                   !< divided array
        integer, intent(in) :: size_ghost                                       !< width of ghost region
        
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
    !> \private 2-D complex version
    integer function get_ghost_arr_2D_complex(arr,size_ghost) result(ierr)
        use num_vars, only: rank, n_procs
        
        character(*), parameter :: rout_name = 'get_ghost_arr_2D_complex'
        
        ! input / output
        complex(dp), intent(inout) :: arr(:,:)                                  !< divided array
        integer, intent(in) :: size_ghost                                       !< width of ghost region
        
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
    !> \private 1-D real version
    integer function get_ghost_arr_1D_real(arr,size_ghost) result(ierr)
        use num_vars, only: rank, n_procs
        
        character(*), parameter :: rout_name = 'get_ghost_arr_1D_real'
        
        ! input / output
        real(dp), intent(in) :: arr(:)                                          !< divided array
        integer, intent(in) :: size_ghost                                       !< width of ghost region
        
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
    
    !> \private real version
    integer function broadcast_var_real(var,source) result(ierr)
        character(*), parameter :: rout_name = 'broadcast_var_real'
        
        ! input / output
        real(dp), intent(in) :: var                                             !< variable to be broadcast
        integer, intent(in), optional :: source                                 !< process that sends
        
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
    !> \private integer version
    integer function broadcast_var_int(var,source) result(ierr)
        character(*), parameter :: rout_name = 'broadcast_var_int'
        
        ! input / output
        integer, intent(in) :: var                                              !< variable to be broadcast
        integer, intent(in), optional :: source                                 !< process that sends
        
        ! local variables
        integer :: source_loc = 0                                               ! local value for source
        
        ! initialize ierr
        ierr = 0
        
        ! set local source if given
        if (present(source)) source_loc = source
        
        call MPI_Bcast(var,1,MPI_INTEGER,source_loc,MPI_COMM_WORLD,ierr)
        CHCKERR('MPI broadcast failed')
    end function broadcast_var_int
    !> \private logical version
    integer function broadcast_var_log(var,source) result(ierr)
        character(*), parameter :: rout_name = 'broadcast_var_log'
        
        ! input / output
        logical, intent(in) :: var                                              !< variable to be broadcast
        integer, intent(in), optional :: source                                 !< process that sends
        
        ! local variables
        integer :: source_loc = 0                                               ! local value for source
        
        ! initialize ierr
        ierr = 0
        
        ! set local source if given
        if (present(source)) source_loc = source
        
        call MPI_Bcast(var,1,MPI_LOGICAL,source_loc,MPI_COMM_WORLD,ierr)
        CHCKERR('MPI broadcast failed')
    end function broadcast_var_log
    !> \private complex array version
    integer function broadcast_var_complex_arr(var,source) result(ierr)
        character(*), parameter :: rout_name = 'broadcast_var_complex_arr'
        
        ! input / output
        complex(dp), intent(in) :: var(:)                                       !< variable to be broadcast
        integer, intent(in), optional :: source                                 !< process that sends
        
        ! local variables
        integer :: source_loc = 0                                               ! local value for source
        
        ! initialize ierr
        ierr = 0
        
        ! set local source if given
        if (present(source)) source_loc = source
        
        call MPI_Bcast(var,size(var),MPI_DOUBLE_COMPLEX,source_loc,&
            &MPI_COMM_WORLD,ierr)
        CHCKERR('MPI broadcast failed')
    end function broadcast_var_complex_arr
    !> \private real array version
    integer function broadcast_var_real_arr(var,source) result(ierr)
        character(*), parameter :: rout_name = 'broadcast_var_real_arr'
        
        ! input / output
        real(dp), intent(in) :: var(:)                                          !< variable to be broadcast
        integer, intent(in), optional :: source                                 !< process that sends
        
        ! local variables
        integer :: source_loc = 0                                               ! local value for source
        
        ! initialize ierr
        ierr = 0
        
        ! set local source if given
        if (present(source)) source_loc = source
        
        call MPI_Bcast(var,size(var),MPI_DOUBLE_PRECISION,source_loc,&
            &MPI_COMM_WORLD,ierr)
        CHCKERR('MPI broadcast failed')
    end function broadcast_var_real_arr
    !> \private integer array version
    integer function broadcast_var_int_arr(var,source) result(ierr)
        character(*), parameter :: rout_name = 'broadcast_var_int_arr'
        
        ! input / output
        integer, intent(in) :: var(:)                                           !< variable to be broadcast
        integer, intent(in), optional :: source                                 !< process that sends
        
        ! local variables
        integer :: source_loc = 0                                               ! local value for source
        
        ! initialize ierr
        ierr = 0
        
        ! set local source if given
        if (present(source)) source_loc = source
        
        call MPI_Bcast(var,size(var),MPI_INTEGER,source_loc,MPI_COMM_WORLD,ierr)
        CHCKERR('MPI broadcast failed')
    end function broadcast_var_int_arr
    !> \private logical array version
    integer function broadcast_var_log_arr(var,source) result(ierr)
        character(*), parameter :: rout_name = 'broadcast_var_log_arr'
        
        ! input / output
        logical, intent(in) :: var(:)                                           !< variable to be broadcast
        integer, intent(in), optional :: source                                 !< process that sends
        
        ! local variables
        integer :: source_loc = 0                                               ! local value for source
        
        ! initialize ierr
        ierr = 0
        
        ! set local source if given
        if (present(source)) source_loc = source
        
        call MPI_Bcast(var,size(var),MPI_LOGICAL,source_loc,MPI_COMM_WORLD,ierr)
        CHCKERR('MPI broadcast failed')
    end function broadcast_var_log_arr
    
    !> Wait for all processes, wrapper to MPI barrier.
    !!
    !! \return ierr
    integer function wait_MPI() result(ierr)
        character(*), parameter :: rout_name = 'wait_MPI'
        
        ! initialize ierr
        ierr = 0
        
        ! set barrier
        call MPI_Barrier(MPI_Comm_world,ierr)
        CHCKERR('MPI Barrier failed')
        
#if ldebug
        n_waits = n_waits + 1
#endif
    end function wait_MPI
    
    !>  Request  access  to  lock  of   a  BL  (blocking)  or  optionally  a  NB
    !! (non-blocking) type.
    !!
    !! \note Based on \cite RossAtomicIO.
    !!
    !! \return ierr
    integer function lock_req_acc(lock,blocking) result(ierr)
        use num_vars, only: rank
        use num_utilities, only: bubble_sort
        
        character(*), parameter :: rout_name = 'lock_req_acc'
        
        ! input / output
        type(lock_type), intent(inout) :: lock                                  !< lock
        
        ! local variables
        integer :: id                                                           ! counter
        integer, allocatable :: wl_loc(:)                                       ! local copy of waiting list
        integer, allocatable :: next_nb_procs(:)                                ! next NB process(es)
        integer, allocatable :: ranks_to_activate(:)                            ! ranks to be set active
        logical, intent(in), optional :: blocking                               ! is a blocking process
        logical :: next_nb_proc_exists                                          ! a next NB process exists
        logical :: direct_receipt                                               ! receipt was direct
#if ldebug
        integer :: istat                                                        ! status
#endif
        
        ! initialize ierr
        ierr = 0
        
#if ldebug
        ! check for initialized lock
        if (.not.allocated(lock%wl)) then
            ierr = 1
            CHCKERR('lock not intialized')
        end if
#endif
        
        ! set blocking property
        lock%blocking = .true.
        if (present(blocking)) lock%blocking = blocking
        
#if ldebug
        if(debug_lock) write(*,*,IOSTAT=istat) trim(lock_header(lock)), &
            &'requesting access'
#endif
        
        ! add current process to waiting list
        ierr = lock_wl_change(1,lock%blocking,lock,wl_loc)
        CHCKERR('')
        
        ! set direct_receipt when wait list empty
        direct_receipt = wl_empty(wl_loc,[-2,-1,1,2])
        
        ! wait for notification if other processes
        if (direct_receipt) then
#if ldebug
            if(debug_lock) write(*,*,IOSTAT=istat) trim(lock_header(lock)), &
                &'and got it right away'
#endif
            
            ! find next NB if this one is NB also
            next_nb_proc_exists = .false.
            if (.not.lock%blocking) then
                ! find all NB procs
                next_nb_proc_exists =  &
                    &.not.wl_empty(wl_loc,[-1],next_procs=next_nb_procs)
#if ldebug
                if (debug_lock) write(*,*,IOSTAT=istat) &
                    &trim(lock_header(lock)),&
                    &'    NB proc, next NB procs found:', next_nb_procs
#endif
            end if
            
            ! set ranks to activate
            if (.not.lock%blocking .and. next_nb_proc_exists) then              ! more processes
                allocate(ranks_to_activate(size(next_nb_procs)+1))
                ranks_to_activate = [rank,next_nb_procs]
            else                                                                ! only this process
                allocate(ranks_to_activate(1))
                ranks_to_activate = rank
            end if
            call bubble_sort(ranks_to_activate)                                 ! sort
            
            ! activate
#if ldebug
            if(debug_lock) write(*,*,IOSTAT=istat) trim(lock_header(lock)), &
                &'    setting status to activate:', ranks_to_activate
#endif
            ierr = lock_wl_change(2,.false.,lock,wl_loc,ranks=ranks_to_activate)
            CHCKERR('')
            
            ! notify other ranks to activate (NB)
            if (next_nb_proc_exists) then
                do id = 1,size(next_nb_procs)
                    ierr = lock_notify(lock,next_nb_procs(id))
                    CHCKERR('')
                end do
            end if
        else
            ! wait to get notified
            ierr = lock_get_notified(lock)
            CHCKERR('')
        end if
    end function lock_req_acc
    
    !> Returns  access  to a  lock. 
    !!
    !! The blocking property has been set when requesting the lock.
    !!
    !! \note Based on \cite RossAtomicIO.
    !!
    !! \return ierr
    integer function lock_return_acc(lock) result(ierr)
        use num_utilities, only: bubble_sort
        
        character(*), parameter :: rout_name = 'lock_return_acc'
        
        ! input / output
        type(lock_type), intent(inout) :: lock                                  !< lock
        
        ! local variables
        integer :: id                                                           ! counter
        integer, allocatable :: wl_loc(:)                                       ! local copy of waiting list
        integer, allocatable :: next_procs(:)                                   ! next process(es)
        integer, allocatable :: ranks_to_activate(:)                            ! ranks to be set active
        logical :: next_proc_exists                                             ! a next process exists
        logical :: next_proc_Bl                                                 ! next process is BL
#if ldebug
        integer :: istat                                                        ! status
#endif
        
        ! initialize ierr
        ierr = 0
        
#if ldebug
        if (.not.allocated(lock%wl)) then
            ierr = 1
            CHCKERR('lock not intialized')
        end if
        
        if(debug_lock) write(*,*,IOSTAT=istat) trim(lock_header(lock)), &
            &'returning access'
#endif
        
        ! remove current process to lock waiting list
        ierr = lock_wl_change(0,lock%blocking,lock,wl_loc)
        CHCKERR('')
        
        ! find next processes if:
        !   - BL: always
        !   - NB: was last NB running
        next_proc_Bl = .false.
        next_proc_exists = .false.
        if (lock%blocking .or. wl_empty(wl_loc,[-2])) then
            ! find all BL procs
            next_proc_exists = .not.wl_empty(wl_loc,[1],next_procs=next_procs)
            if (next_proc_exists) next_proc_Bl = .true.
#if ldebug
            if (debug_lock) write(*,*,IOSTAT=istat) trim(lock_header(lock)),&
                &'    next BL procs found:', next_procs
#endif
            
            ! if none found, try NB procs
            if (.not. next_proc_exists) then
                next_proc_exists = .not.wl_empty(wl_loc,[-1],&
                    &next_procs=next_procs)
                if (next_proc_exists) next_proc_Bl = .false.
#if ldebug
                if (debug_lock) write(*,*,IOSTAT=istat) &
                    &trim(lock_header(lock)),&
                    &'    next NB procs found:', next_procs
#endif
            end if
#if ldebug
        else
            if(debug_lock) write(*,*,IOSTAT=istat) trim(lock_header(lock)), &
                &'but has no notification rights'
#endif
        end if
        
        if (next_proc_exists) then
            ! set ranks to activate
            if (next_proc_Bl) then                                              ! only the BL process
                allocate(ranks_to_activate(1))
                ranks_to_activate = next_procs(1)
            else                                                                ! multiple NB processes
                allocate(ranks_to_activate(size(next_procs)))
                ranks_to_activate = next_procs
            end if
            call bubble_sort(ranks_to_activate)                                 ! sort
            
            ! activate
#if ldebug
            if(debug_lock) write(*,*,IOSTAT=istat) trim(lock_header(lock)), &
                &'    setting status to activate:', ranks_to_activate
#endif
            ierr = lock_wl_change(2,next_proc_BL,lock,wl_loc,&
                &ranks=ranks_to_activate)
            CHCKERR('')
            
            ! notify ranks to activate (NB or BL)
            if (next_proc_exists) then
                do id = 1,size(ranks_to_activate)
                    ierr = lock_notify(lock,ranks_to_activate(id))
                    CHCKERR('')
                end do
            end if
        end if
    end function lock_return_acc
    
    !> Decides whether a waiting list is empty.
    !!
    !! The type of process to find is  indicated by an array of possible values.
    !!
    !! \see See lock_wl_change() for an explanation of the process type.
    !!
    !! Additionally, for NB processes, the  negative inverse of these values are
    !! used.
    !!
    !! If the waiting list is not  empty, the next process(es) can optionally be
    !! returned.
    !!
    !! \note Based on \cite RossAtomicIO.
    !!
    !! \return ierr
    logical function wl_empty(wl,proc_type,next_procs)
        use num_vars, only: n_procs, rank
        
        ! input / output
        integer, intent(in) :: wl(:)                                            !< waiting list
        integer, intent(in) :: proc_type(:)                                     !< types of processes accepted
        integer, intent(inout), optional, allocatable :: next_procs(:)          !< next process(es) if not empty
        
        ! local variables
        integer :: id, jd                                                       ! counters
        integer :: nr_np                                                        ! nr. of next processes found
        integer, allocatable :: next_procs_loc(:)                               ! local next process
        
        ! initialize
        nr_np = 0
        wl_empty = .true.
        allocate(next_procs_loc(n_procs))
        
        ! find the next rank
        do id = 1,n_procs-1
            ! wrap around n_procs-1 to zero
            next_procs_loc(nr_np+1) = mod(rank+id,n_procs)
            
            proc_types: do jd = 1,size(proc_type)
                if (wl(next_procs_loc(nr_np+1)+1).eq.proc_type(jd)) then        ! arrays start at 1, ranks at 0
                    wl_empty = .false.
                    nr_np = nr_np + 1
                    exit proc_types
                end if
            end do proc_types
        end do
        
        if (.not.wl_empty .and. present(next_procs)) then
            if (allocated(next_procs)) deallocate(next_procs)
            allocate(next_procs(nr_np))
            next_procs = next_procs_loc(1:nr_np)
        end if
    end function wl_empty
    
    !> Notifies a rank that they can get the lock.
    !!
    !! The signal sent is the rank + 1.
    !!
    !! \note Based on \cite RossAtomicIO.
    !!
    !! \return ierr
    integer function lock_notify(lock_loc,rec_rank) result(ierr)
        use num_vars, only: rank
        
        character(*), parameter :: rout_name = 'lock_notify'
        
        ! input / output
        type(lock_type), intent(in) :: lock_loc                                 !< lock
        integer, intent(in) :: rec_rank                                         !< receiving rank
        
        ! local variables
#if ldebug
        integer :: istat                                                        ! status
#endif
        
        ! initialize ierr
        ierr = 0
        
        call MPI_Send(rank+1,1,MPI_INTEGER,rec_rank,lock_loc%wu_tag,&
            &MPI_Comm_world,ierr)
        CHCKERR('Failed to send notification')
        
#if ldebug
        if (debug_lock) write(*,*,IOSTAT=istat) trim(lock_header(lock_loc)), &
            &'    notified', rec_rank
#endif
    end function lock_notify
    
    !> Get notified that the rank can  get the lock.
    !!
    !! \note Based on \cite RossAtomicIO.
    !!
    !! \return ierr
    integer function lock_get_notified(lock_loc) result(ierr)
        character(*), parameter :: rout_name = 'lock_get_notified'
        
        ! input / output
        type(lock_type), intent(in) :: lock_loc                                 !< lock
        
        ! local variables
        integer :: dum_buf                                                      ! dummy buffer
#if ldebug
        integer :: istat                                                        ! status
#endif
        
        ! initialize ierr
        ierr = 0
        
#if ldebug
        if (debug_lock) write(*,*,IOSTAT=istat) trim(lock_header(lock_loc)), &
            &'    but needs to wait for lock'
#endif
        call MPI_Recv(dum_buf,1,MPI_INTEGER,MPI_ANY_SOURCE,&
            &lock_loc%wu_tag,MPI_Comm_world,MPI_STATUS_IGNORE,ierr)
        CHCKERR('Failed to receive notification')
#if ldebug
        if (debug_lock) write(*,*,IOSTAT=istat) trim(lock_header(lock_loc)), &
            &'    got notified by ',dum_buf-1
#endif
    end function lock_get_notified
    
    !> Adds, removes or sets  to active a rank from the waiting  list for a lock
    !! and returns the lock waiting list:
    !!
    !! Actions:
    !!  - \c wl_action = 0: remove
    !!  - \c wl_action = 1: add
    !!  - \c wl_action = 2: active
    !!
    !! Or negative equivalents for non-blocking (NB) procs.
    !!
    !! Optionally, the rank(s)  of the process for which to  perform this action
    !! can  be passed.  This is  useful for  doing the  same action  on multiple
    !! processes.
    !!
    !! \note Based on \cite RossAtomicIO.
    !!
    !! \ldebug
    !!
    !! \return ierr
    integer function lock_wl_change(wl_action,blocking,lock_loc,wl,ranks) &
        &result(ierr)
        use num_vars, only: n_procs, rank
        
        character(*), parameter :: rout_name = 'lock_wl_change'
        
        ! input / output
        integer, intent(in) :: wl_action                                        !< action to perform
        logical, intent(in) :: blocking                                         !< the ranks to be changed are blocking
        type(lock_type), intent(inout) :: lock_loc                              !< lock
        integer, intent(inout), allocatable :: wl(:)                            !< waiting list
        integer, intent(in), optional :: ranks(:)                               !< rank(s) for which to perform option
        
        ! local variables
        integer :: id                                                           ! counter
        integer :: n_ranks                                                      ! number of ranks to set
        integer :: put_val                                                      ! value to put
        integer, allocatable :: ranks_loc(:)                                    ! local ranks
        integer(kind=MPI_ADDRESS_KIND) :: one = 1                               ! one
        integer(kind=MPI_ADDRESS_KIND) :: disp                                  ! local displacement
        integer :: ln                                                           ! local length
#if ldebug
        integer :: istat                                                        ! status
        integer(kind=8) :: window_time(2)                                       ! time that a window is locked
#endif
        
        ! initialize ierr
        ierr = 0
        
        ! set value to put
        put_val = wl_action
        if (.not.blocking) put_val = -put_val                                   ! NB indicated in waiting list using a negative value
        !allocate(wl_excl(n_procs-1))
        if (allocated(wl)) deallocate(wl)
        allocate(wl(n_procs))
        
        ! set local ranks if not present
        if (.not.present(ranks)) then
            allocate(ranks_loc(1))
            ranks_loc(1) = rank
        else
            allocate(ranks_loc(size(ranks)))
            ranks_loc = ranks
        end if
        n_ranks = size(ranks_loc)
        
        ! get window to lock
#if ldebug
        if (debug_lock) call system_clock(window_time(1))
#endif
        call MPI_Win_lock(MPI_LOCK_EXCLUSIVE,0,0,lock_loc%wl_win,ierr)
        CHCKERR('Failed to lock window')
        
        ! get previous list and put value for current ranks
        ! Do this by looking at every rank that is to be put, whether there are
        ! values to be gotten between the previous put and the current put.
        ! Finally, do one more get for the ranks between  the last put  and the
        ! end.
        do id = 1,n_ranks+1
            ! set local displacement and length
            if (id.eq.1) then
                disp = 0
                ln = ranks_loc(1)
            else
                disp = disp + ln + 1                                            ! previous step used ln gets and 1 put
                if (id.lt.n_ranks+1) then
                    ln = ranks_loc(id)-ranks_loc(id-1)-1
                else
                    ln = n_procs-ranks_loc(id-1)-1
                end if
            end if
            
            ! perform the gets between this rank and previous rank
            if (ln.gt.0) then                                                   ! if there is something to get between last put and this one
                call MPI_Get(wl(int(disp)+1:int(disp)+ln),ln,MPI_INTEGER,0,&
                    &disp,ln,MPI_INTEGER,lock_loc%wl_win,ierr)
                CHCKERR('Failed to get waiting list')
            end if
            
            ! perform the single put up to the current rank if not last piece
            if (id.lt.n_ranks+1) then
                call MPI_Put(put_val,1,MPI_INTEGER,0,ranks_loc(id)*one,1,&
                    &MPI_INTEGER,lock_loc%wl_win,ierr)
                CHCKERR('Failed to add to waiting list')
                wl(ranks_loc(id)+1) = put_val
            end if
        end do
        
        ! unlock window
        call MPI_Win_unlock(0,lock_loc%wl_win,ierr)
        CHCKERR('Failed to lock window')
#if ldebug
        if (debug_lock) call system_clock(window_time(2))
#endif

#if ldebug
        if(debug_lock) write(*,*,IOSTAT=istat) trim(lock_header(lock_loc)),&
            &'    time: '//&
            &trim(r2strt(1.E-9_dp*(window_time(2)-window_time(1)))),&
            &' waiting list:', wl
#endif
    end function lock_wl_change

#if ldebug
    !> Returns the header for lock debug messages.
    !!
    !! \ldebug
    character(len=max_str_ln) function lock_header(lock_loc) result(header)
        use num_vars, only: rank
        
        ! input / output
        type(lock_type), intent(in) :: lock_loc                                 !< lock
        
        ! local variables
        integer(kind=8) :: clock                                                ! current clock
        character(len=2) :: block_char                                          ! BL for blocking and NB for non-blocking processes
        
        ! get clock
        call system_clock(clock)
        
        ! set block_char
        if (lock_loc%blocking) then
            block_char = 'BL'
        else
            block_char = 'NB'
        end if
        
        header = trim(ii2str(clock))//' '//trim(i2str(lock_loc%wu_tag))//' '//&
            &block_char//' '//trim(i2str(rank))//'-'
    end function lock_header
#endif
end module MPI_utilities
