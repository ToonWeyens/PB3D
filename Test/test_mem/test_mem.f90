module constants
    use ISO_FORTRAN_ENV
    
    ! constants
    integer, parameter :: dp = REAL64                                           ! double precision
end module constants

module array_class
    use constants
    
    implicit none
    
    type array_type
        real(dp), allocatable :: p(:)
    contains
        procedure :: dealloc => dealloc_array_type
    end type
contains
    subroutine dealloc_array_type(array)
        class(array_type), intent(out) :: array
    end subroutine dealloc_array_type
end module array_class

program test_mem
    use array_class
    use constants
    
    implicit none
    
    ! global variables
    integer :: mem_usage_i
    integer :: mem_id = 0
    integer, parameter :: max_str_ln = 300
    real(dp), allocatable :: large_array_1(:)
    integer :: id, jd
    integer :: arr_size
    real(dp) :: arr_mem_size
    integer :: dp_size
    type(array_type), allocatable :: large_array_2(:)
    integer :: nr_arrs
    
    ! open mem_usage file
    mem_usage_i = 10
    open(UNIT=mem_usage_i,FILE='mem_usage.dat',STATUS='replace')
    call writo('memory usage file opened')
    
    ! allocate a large array
    call writo('1. allocating and deallocating array')
    write(*,*)
    call writo('Array of size 10^n, give n?')
    read(*,*) arr_size
    arr_size = 10**arr_size
    dp_size = sizeof(1._dp)
    arr_mem_size = arr_size*dp_size*1.E-6_dp
    allocate(large_array_1(arr_size))
    call writo('large array allocated with size '//&
        &trim(r2strt(arr_mem_size))//' MB')
    
    ! write fifths
    do id = 1,5
        large_array_1((id-1)*arr_size/5+1:id*arr_size/5) = id*1._dp
        call writo('writing '//trim(i2str(id))//'st fifth')
    end do
    
    ! deallocate
    deallocate(large_array_1)
    call writo('large array deallocated')
    write(*,*)
    
    ! use user type
    call writo('2. allocating and deallocating custom type')
    write(*,*)
    
    call writo('How many arrays?')
    read(*,*) nr_arrs
    allocate(large_array_2(nr_arrs))
    call writo('Allocated large array with '//trim(i2str(nr_arrs))//' elements')
    
    do jd = 1,nr_arrs
        allocate(large_array_2(jd)%p(arr_size))
        call writo('large array '//trim(i2str(jd))//' allocated with size '//&
            &trim(r2strt(arr_mem_size))//' MB')
        
        ! write fifths
        do id = 1,5
            large_array_2(jd)%p((id-1)*arr_size/5+1:id*arr_size/5) = id*1._dp
            call writo('writing '//trim(i2str(id))//'st fifth')
        end do
    end do
    
    ! deallocate
    !do jd = 1,nr_arrs
        !call large_array_2(jd)%dealloc
        !call writo('large array '//trim(i2str(jd))//' deallocated')
        !write(*,*)
    !end do
    deallocate(large_array_2)
    call writo('large array deallocated all at once')
    
    write(*,*)
contains
    subroutine writo(input_str)
        ! input / output
        character(len=*), intent(in) :: input_str                               ! the name that is searched for
        
        ! local variables
        character(len=max_str_ln) :: output_str
        integer :: mem_usage
        
        ! increment mem_id and get memory usage
        mem_id = mem_id + 1
        mem_usage = get_mem_usage()
        
        ! set output string
        output_str = trim(input_str)//' - ['//trim(i2str(mem_id))//': '//&
            &trim(i2str(mem_usage))//'MB]'
        
        write(*,*) trim(output_str)
        write(mem_usage_i,*) mem_id, mem_usage
    end subroutine writo
    
    elemental character(len=max_str_ln) function i2str(k)
        ! input / output
        integer, intent(in) :: k
        
        write (i2str, *) k
        i2str = adjustl(i2str)
    end function i2str
    elemental character(len=max_str_ln) function r2str(k)
        ! input / output
        real(dp), intent(in) :: k
        
        write (r2str, '(ES23.16)') k
        r2str = adjustl(r2str)
    end function r2str
    elemental character(len=max_str_ln) function r2strt(k)
        ! input / output
        real(dp), intent(in) :: k
        
        write (r2strt, '(ES9.2)') k
        r2strt = adjustl(r2strt)
    end function r2strt
    
    integer function get_mem_usage() result(mem)
        ! local variables
        character(len=200):: filename=' '                                       ! name of file where memory stored
        character(len=80) :: line                                               ! line of memory file
        character(len=8)  :: pid_char=' '                                       ! process ID in string
        logical :: exists                                                       ! whether memory file exists
        integer :: pid                                                          ! process ID
        integer :: istat                                                        ! status
        
        ! initiazlie mem to negative number
        mem = -1
        
        ! get process ID and write it in string
        pid = getpid()
        write(pid_char,'(I8)') pid
        
        ! set up file name
        filename='/proc/'//trim(adjustl(pid_char))//'/status'
        
        ! inquire about memory file
        inquire (file=filename,exist=exists)
        if (.not.exists) return
        
        ! open and read memory file
        open(UNIT=mem_usage_i-1,FILE=filename,ACTION='read',IOSTAT=istat)
        if (istat.ne.0) return
        do
            read (mem_usage_i-1,'(A)',IOSTAT=istat) line
            if (istat.ne.0) return
            if (line(1:6).eq.'VmRSS:') then
                read (line(7:),*) mem
                exit
            endif
        end do
        
        mem = mem/1000                                                          ! convert to MB
        
        ! close memory file
        close(mem_usage_i-1)
    end function get_mem_usage
end program test_mem
