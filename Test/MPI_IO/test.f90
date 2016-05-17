program test
    use MPI
    implicit none
    
    ! local variables
    integer :: ierr, rank, n_procs
    integer :: max_str_ln = 20
    
    write(*,*) '!! THIS IS NOT WORKING !!'
    
    call MPI_init(ierr)
    if (ierr.ne.0) write(*,*) rank, 'ERROR: '//trim(i2str(ierr))
    call MPI_Comm_rank(MPI_Comm_world,rank,ierr)
    if (ierr.ne.0) write(*,*) rank, 'ERROR: '//trim(i2str(ierr))
    call MPI_Comm_size(MPI_Comm_world,n_procs,ierr)
    if (ierr.ne.0) write(*,*) 'ERROR: '//trim(i2str(ierr))
    write(*,*) 'rank', rank, rank, 'n_procs', n_procs
    
    call run_test()
    
    call MPI_finalize(ierr)
    if (ierr.ne.0) write(*,*) 'ERROR: '//trim(i2str(ierr))
contains
    subroutine run_test()
        character(len=max_str_ln) :: file_name
        integer :: fh
        integer :: id
        integer(kind=MPI_OFFSET_KIND) disp
        character(len=max_str_ln) :: buf
        
        file_name = 'test_file.txt'
        buf = 'This is a test string'
        
        call MPI_File_open(MPI_COMM_WORLD,trim(file_name),MPI_MODE_CREATE&
            &+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
        if (ierr.ne.0) write(*,*) rank, 'ERROR: '//trim(i2str(ierr))
        
        call MPI_File_set_view(fh,disp,MPI_CHARACTER,MPI_CHARACTER,'native',&
            &MPI_INFO_NULL,ierr) 
        if (ierr.ne.0) write(*,*) rank, 'ERROR: '//trim(i2str(ierr))
        
        call MPI_File_write(fh,buf,max_str_ln,MPI_CHARACTER,MPI_STATUS_IGNORE,&
            &ierr)
        if (ierr.ne.0) write(*,*) rank, 'ERROR: '//trim(i2str(ierr))
    end subroutine run_test
    
    elemental character(len=max_str_ln) function i2str(k)
        ! input / output
        integer, intent(in) :: k
        
        write (i2str, *) k
        i2str = adjustl(i2str)
    end function i2str
end program
