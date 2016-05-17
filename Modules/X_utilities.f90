!------------------------------------------------------------------------------!
!   Numerical utilities related to perturbation operations                     !
!------------------------------------------------------------------------------!
module X_utilities
#include <PB3D_macros.h>
    use str_ops
    use messages
    use num_vars, only: dp, max_name_ln, iu, max_str_ln
    use grid_vars, only: grid_type
    use X_vars, only: n_mod_X
    
    implicit none
    
    private
    public get_suffix, is_necessary_X, divide_X_jobs
    
    ! interfaces
    interface get_suffix
        module procedure get_suffix_1, get_suffix_2
    end interface
    
contains
    ! Sets the suffix used to refer to a perturbation quantity.
    function get_suffix_1(id,lim_sec_X) result(res)                             ! vectorial version
        ! input / output
        integer, intent(in) :: id                                               ! mode index
        integer, intent(in), optional :: lim_sec_X(2)                           ! limits of m_X (pol. flux) or n_X (tor. flux)
        integer :: res                                                          ! output
        
        ! local variables
        integer :: lim_sec_X_loc(2)                                             ! local lim_sec_X
        
        ! set local lim_sec_X
        lim_sec_X_loc = [1,n_mod_X]
        if (present(lim_sec_X)) lim_sec_X_loc = lim_sec_X
        
        ! set suffix
        res = lim_sec_X_loc(1)-1+id
    end function get_suffix_1
    function get_suffix_2(id,jd,lim_sec_X) result(res)                          ! tensorial version
        ! input / output
        integer, intent(in) :: id, jd                                           ! mode indices
        integer, intent(in), optional :: lim_sec_X(2,2)                         ! limits of m_X (pol flux) or n_X (tor flux) for both dimensions
        integer :: res(2)                                                       ! output
        
        ! local variables
        integer :: lim_sec_X_loc(2,2)                                           ! local lim_sec_X
        
        ! set local lim_sec_X
        lim_sec_X_loc(:,1) = [1,n_mod_X]
        lim_sec_X_loc(:,2) = [1,n_mod_X]
        if (present(lim_sec_X)) lim_sec_X_loc = lim_sec_X
        
        ! set suffix
        res = [lim_sec_X_loc(1,1)-1+id,lim_sec_X_loc(1,2)-1+jd]
    end function get_suffix_2
    
    ! Determines whether a variable needs to be  considered: Only if it is on or
    ! below the diagonal for symmetric quantities.
    logical function is_necessary_X(sym,sec_X_id,lim_sec_X) result(res)
        ! input / output
        logical, intent(in) :: sym                                              ! whether the variable is symmetric
        integer, intent(in) :: sec_X_id(2)                                      ! mode indices
        integer, intent(in), optional :: lim_sec_X(2,2)                         ! limits of m_X (pol flux) or n_X (tor flux) for both dimensions
        
        ! local variables
        integer :: lim_sec_X_loc(2,2)                                           ! local version of lim_sec_X
        
        ! initialize res
        res = .true.
        
        ! set local lim_sec_X
        lim_sec_X_loc(:,1) = [1,n_mod_X]
        lim_sec_X_loc(:,2) = [1,n_mod_X]
        if (present(lim_sec_X)) lim_sec_X_loc = lim_sec_X
        
        ! modify res depending on symmetry
        if (sym) then
            if (lim_sec_X_loc(1,1)+sec_X_id(1).lt.&
                &lim_sec_X_loc(1,2)+sec_X_id(2)) res = .false.
        end if
    end function is_necessary_X
    
    ! divides the perturbation jobs
    ! The (k,m) pairs have to be calculated, but might have to be broken up into
    ! pieces  due  to memory  constraints.  In  its  most  extreme case,  for  a
    ! certain mode  number k, therefore, all  the pairs (k,m) can  be calculated
    ! sequentially, saving  the values  for this  k in  the process  and cycling
    ! through m.
    ! This  can be  directly  scaled  up to  a  block of  k  and  m values,  and
    ! ultimately, all the values simultaneously.
    ! Therefore, the whole  load is divided into jobs depending  on the sizes of
    ! the blocks of k  and m values in memory. These jobs  start with the vector
    ! phase  by calculating  U and  DU, followed  in the  tensor phase  by their
    ! combinations in  PV and  KV. Then,  these are  integrated in  the parallel
    ! coordinate,  with  a  possible  interpolation  in  between.  Finally,  the
    ! integrated values are saved and the next jobs starts.
    ! This  function does  the  job of  dividing the  grids  setting the  global
    ! variables 'X_jobs_lims'  and 'X_jobs_taken' for  data of a  certain order,
    ! given by div_ord.  E.g. for the vector  phase, the order is 1  and for the
    ! tensorial phase it is 2.
    integer function divide_X_jobs(arr_size,div_ord) result(ierr)
        use num_vars, only: max_X_mem_per_proc, n_procs, X_jobs_lims, rank, &
            &X_jobs_file_name, X_jobs_taken, X_jobs_lock_file_name, &
            &mem_scale_fac
        use X_vars, only: n_mod_X
        use files_utilities, only: nextunit
        use MPI_utilities, only: wait_MPI
        
        character(*), parameter :: rout_name = 'divide_X_jobs'
        
        ! input / output
        integer, intent(in) :: arr_size                                         ! array size (using loc_n_r)
        integer, intent(in) :: div_ord                                          ! division order
        
        ! local variables
        integer :: n_div                                                        ! factor by which to divide the total size
        real(dp) :: mem_size                                                    ! approximation of memory required for X variables
        integer :: n_mod_block                                                  ! nr. of modes in block
        integer, allocatable :: n_mod_loc(:)                                    ! number of modes per block
        character(len=max_str_ln) :: block_message                              ! message about how many blocks
        character(len=max_str_ln) :: err_msg                                    ! error message
        integer :: file_i                                                       ! file number
        
        ! initialize ierr
        ierr = 0
        
        ! user output
        call writo('Dividing the perturbation jobs of order '//&
            &trim(i2str(div_ord)))
        call lvl_ud(1)
        
        ! calculate largest possible block of (k,m) values
        n_div = 0
        mem_size = huge(1._dp)
        do while (mem_size.gt.max_X_mem_per_proc)
            n_div = n_div + 1
            n_mod_block = ceiling(n_mod_X*1._dp/n_div)
            ierr = calc_memory(div_ord,arr_size,n_mod_block,mem_size)
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
                &' and '//trim(i2str(n_div**div_ord))//&
                &' jobs are done separately'
            if (n_procs.lt.n_div**div_ord) then
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
        call writo('The memory per process is estimated to be about '//&
            &trim(i2str(ceiling(mem_size)))//'MB (maximum: '//&
            &trim(i2str(ceiling(max_X_mem_per_proc)))//'MB)')
        
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
        X_jobs_lims = calc_X_jobs_lims(n_mod_loc,div_ord)
        if (allocated(X_jobs_taken)) deallocate(X_jobs_taken)
        allocate(X_jobs_taken(n_div**div_ord))
        X_jobs_taken = .false.
        
        ! initialize file with global variable
        if (rank.eq.0) then
            open(STATUS='REPLACE',unit=nextunit(file_i),file=X_jobs_file_name,&
                &iostat=ierr)
            CHCKERR('Cannot open perturbation jobs file')
            write(file_i,*) '# Process, X job'
            close(file_i,iostat=ierr)
            CHCKERR('Failed to close file')
        end if
        
        ! delete possible lock file
        if (rank.eq.0) then
            open(unit=nextunit(file_i),file=X_jobs_lock_file_name,iostat=ierr)
            CHCKERR('Failed to open lock file')
            close(file_i,status='DELETE',iostat=ierr)
            CHCKERR('Failed to delete lock file')
        end if
        
        ! synchronize MPI
        ierr = wait_MPI()
        CHCKERR('')
        
        ! user output
        call lvl_ud(-1)
        call writo('Perturbation jobs divided')
    contains
        ! Calculate memory in MB necessary for X variables of a certain order
        !   - order 1: 4x n_par_X x n_geo x loc_n_r x n_mod
        !   - order 2: 2x n_par_X x n_geo x loc_n_r x n_mod^2
        !              4x n_par_X x n_geo x loc_n_r x n_mod(n_mod+1)/2
        !   - higher order: not used
        ! where n_par_X  x n_geo x  loc_n_r should  be passed as  'arr_size' and
        ! n_mod as well
        integer function calc_memory(ord,arr_size,n_mod,mem_size) result(ierr)
            use ISO_C_BINDING
            
            character(*), parameter :: rout_name = 'calc_memory'
            
            ! input / output
            integer, intent(in) :: ord                                          ! order of data
            integer, intent(in) :: arr_size                                     ! size of part of X array
            integer, intent(in) :: n_mod                                        ! number of modes
            real(dp), intent(inout) :: mem_size                                 ! total size
            
            ! local variables
            integer(C_SIZE_T) :: dp_size                                        ! size of dp
            character(len=max_str_ln) :: err_msg                                ! error message
            
            ! initialize ierr
            ierr = 0
            
            call lvl_ud(1)
            
            ! get size of complex variable
            dp_size = 2*sizeof(1._dp)                                           ! complex variable
            
            ! calculate memory depending on order
            select case(ord)
                case (1)                                                        ! vectorial data: U, DU
                    ! set memory size
                    mem_size = 4*arr_size
                    mem_size = mem_size*n_mod**ord
                    mem_size = mem_size*dp_size
                case (2)                                                        ! tensorial data: PV, KV
                    ! set memory size
                    mem_size = arr_size
                    mem_size = mem_size*(2*n_mod**ord+4*n_mod*(n_mod+1)/2)
                    mem_size = mem_size*dp_size
                case default
                    ierr = 1
                    err_msg = 'Orders > 2 are not implemented'
                    CHCKERR(err_msg)
            end select
            
            ! convert B to MB
            mem_size = mem_size*1.E-6_dp
            
            ! scale memory to account for rough estimation
            mem_size = mem_size*mem_scale_fac
            
            ! test overflow
            if (mem_size.lt.0) then
                ierr = 1
                CHCKERR('Overflow occured')
            end if
            
            call lvl_ud(-1)
        end function calc_memory
        
        ! Calculate X_jobs_lims.
        recursive function calc_X_jobs_lims(n_mod,ord) result(res)
            ! input / output
            integer, intent(in) :: n_mod(:)                                     ! X jobs data
            integer, intent(in) :: ord                                          ! order of data
            integer, allocatable :: res(:,:)                                    ! result
            
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
                    res_loc = calc_X_jobs_lims(n_mod,ord-1)
                    res(1:2*ord-2,(id-1)*n_div**(ord-1)+1:id*n_div**(ord-1)) = &
                        &res_loc
                end if                                                          ! set for own order
                res(2*ord-1,(id-1)*n_div**(ord-1)+1:id*n_div**(ord-1)) = &
                    &sum(n_mod(1:id-1))+1
                res(2*ord,(id-1)*n_div**(ord-1)+1:id*n_div**(ord-1)) = &
                    &sum(n_mod(1:id))
            end do
        end function calc_X_jobs_lims
    end function divide_X_jobs
end module X_utilities
