!------------------------------------------------------------------------------!
!> Variables pertaining to MPI.
!------------------------------------------------------------------------------!
module MPI_vars
#include <PB3D_macros.h>
    use num_vars, only: max_str_ln, dp, plot_dir, data_dir, script_dir
    use str_utilities, only: i2str, r2str, r2strt
    use MPI
    
    implicit none
    private
    public init_lock, dealloc_lock, &
        &HDF5_lock

    !> lock type
    !!
    !! There is a blocking (BL) and  a nonblocking (NB) version where the former
    !! requires an exclusive lock and the latter  a shared one. This is saved in
    !! the variable \c blocking.
    !!
    !! NB processes  that get the lock  directly on request (meaning  that there
    !! were no  other processes in  the queue) notify  directly all the  next NB
    !! processes after gaining access. It also sets their status to active. When
    !! a NB process gains the lock when notified after waiting, it does not have
    !! to check for other  NB processes, as this has been  done by the notifying
    !! process.
    !!
    !! A BL process retains exclusive access upon receipt of the lock. Similarly
    !! to NB processes, if the receipt was  direct on request, the status is set
    !! to active, but only of this NB process.
    !!
    !! When returning the lock, all BL  processes and NB that find themselves to
    !! be the last  active NB process, scan  the waiting list and  pass the lock
    !! preferably to another BL process to notify. If not available, it searches
    !! for all the NB processes to notify together.
    !!
    !! The advantage of prefering BL processes after finishing a process is that
    !! this  way  NB  processes  are  accumulated,  and  then  quickly  finished
    !! afterwards.
    !!
    !! \note  Every process  in  the  waiting queue  will  eventually receive  a
    !! notification.
    !!
    !! Scheme:
    !! <table>
    !!  <tr> <th> <th> request access <th> gain acces <th> return access
    !!  <tr> <th> BL
    !!      <td> <ul><li> add to queue <li> wait </ul>
    !!      <td> if direct: <ul><li> activate </ul>
    !!      <td> <ul><li> remove from queue <li> find next BL/NB
    !!          <li> notify <li> activate all </ul>
    !!  <tr> <th> NB
    !!      <td> <ul><li> add to queue <li> wait </ul>
    !!      <td> if direct: <ul><li> find next NB </ul>
    !!          for all next NB: <ul> <li> notify <li> activate </ul>
    !!      <td> always: <ul><li> remove from queue </ul>
    !!          if last NB: <ul><li> find next BL/NB
    !!          <li> notify <li> activate all </ul>
    !! </table>
    !! with preference BL > NB(s).
    !!
    !! \see This is based on an extension of the ideas in \c RossAtomicIO.
    type, public :: lock_type
        integer, allocatable :: wl(:)                                           !< waiting list
        integer :: wl_win                                                       !< window to waiting list
        integer :: wu_tag                                                       !< wakeup tag
        logical :: blocking                                                     !< is a normal blocking process
    contains
        procedure :: init => init_lock                                          !< initialize
        procedure :: dealloc => dealloc_lock                                    !< deallocate
    end type lock_type
    
    ! global variables
    type(lock_type) :: HDF5_lock                                                ! HDF5 lock
    
contains
    !> Initializes a lock.
    !!
    !! \note
    !!  -# Should be called collectively.
    !!  -# Every lock should have a unique wakeup tag.
    !!
    !! \return ierr
    integer function init_lock(lock_loc,wu_tag) result(ierr)
        use num_vars, only: n_procs, rank
        
        character(*), parameter :: rout_name = 'init_lock'
        
        ! input / output
        class(lock_type), intent(inout) :: lock_loc                             !< lock
        integer, intent(in) :: wu_tag                                           !< wakeup tag
        
        ! local variables
        integer(kind=MPI_ADDRESS_KIND) :: intlb                                 ! lower bound of int type
        integer(kind=MPI_ADDRESS_KIND) :: intex                                 ! extent of int type
        
        ! initialize ierr
        ierr = 0
        
        call MPI_Type_get_extent(MPI_INTEGER,intlb,intex,ierr)
        CHCKERR('Failed to get extent of int')
        
        if (rank.eq.0) then                                                     ! master
            allocate(lock_loc%wl(n_procs))
            
            lock_loc%wl = 0
            call MPI_Win_create(lock_loc%wl,n_procs*intex,int(intex),&
                &MPI_INFO_NULL,MPI_Comm_world,lock_loc%wl_win,ierr)
            CHCKERR('Failed to create window to lock_loc')
        else
            allocate(lock_loc%wl(0))
            
            call MPI_Win_create(lock_loc%wl,0*intex,int(intex),MPI_INFO_NULL,&
                &MPI_Comm_world,lock_loc%wl_win,ierr)
            CHCKERR('Failed to create window to lock_loc')
        end if
        
        lock_loc%wu_tag = wu_tag
        
        ! set fences
        call MPI_Win_fence(0,lock_loc%wl_win,ierr) 
        CHCKERR('Couldn''t set fence') 
    end function init_lock
    
    !> Deallocates a lock.
    !!
    !! \note Should be called collectively.
    !!
    !! \return ierr
    integer function dealloc_lock(lock_loc) result(ierr)
        character(*), parameter :: rout_name = 'dealloc_lock'
        
        ! input / output
        class(lock_type), intent(inout) :: lock_loc                             !< lock
        
        ! initialize ierr
        ierr = 0
        
        ! free lock window
        call MPI_Win_free(lock_loc%wl_win,ierr)
        CHCKERR('Failed to free window to HDF5_lock')
        
        ! deallocate lock
        deallocate(lock_loc%wl)
    end function dealloc_lock
end module MPI_vars
