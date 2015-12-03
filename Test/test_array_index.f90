program test_array_index
    use ISO_FORTRAN_ENV
    
    implicit none
    
    ! local variables
    integer, parameter :: dp = REAL64                                           ! double precision
    real(dp), allocatable :: var_out(:,:)                                       ! variable to be passed
    integer :: var_size(2)                                                      ! size of variable
    integer :: ind_shift                                                        ! shift of first index
    integer :: id                                                               ! counter
    integer :: ierr                                                             ! error variable
    real(dp), allocatable :: var_in(:,:)                                        ! variable to be allocated in other function
    
    ! set up variable
    ind_shift = 3
    var_size = [2,3]
    allocate(var_out(1-ind_shift:var_size(1)-ind_shift,var_size(2)))
    var_out = reshape([(id,id=1,product(var_size))],var_size)
    
    ! pass variable
    write(*,*) 'passing variable'
    ierr = export_variable(var_out)
    
    ! allocate variable in other function
    write(*,*) 'allocate variable in other function'
    ierr = allocate_variable(var_in,var_size,ind_shift)
    do id = 1-ind_shift-1,size(var_in,1)-ind_shift+1
        write(*,*) 'var_in(',id,') = ',var_in(id,:)
    end do
contains
    integer function export_variable(var_in) result(ierr)
        ! input / output
        real(dp), intent(inout) :: var_in(:,:)                                  ! variable to be passed
        
        ! local variables
        integer :: id, jd                                                       ! counters
        
        ! initialize ierr
        ierr = 0
        
        write(*,*) 'size of var_in = ', shape(var_in)
        
        do id = 1-1,size(var_in,1)+1
            write(*,*) 'var_in(',id,') = ',var_in(id,:)
        end do
    end function export_variable
    
    integer function allocate_variable(var_in,var_size,ind_shift) result(ierr)
        ! input / output
        real(dp), intent(inout), allocatable :: var_in(:,:)                     ! variable to be allocated
        integer, intent(in) :: var_size(2)                                      ! size of variable
        integer, intent(in) :: ind_shift                                        ! shift of first index
        
        ! initialize ierr
        ierr = 0
        
        ! set up variable
        allocate(var_in(1-ind_shift:var_size(1)-ind_shift,var_size(2)))
        var_in = reshape([(id,id=1,product(var_size))],var_size)
    end function allocate_variable
end program test_array_index
