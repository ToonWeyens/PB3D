! This module contains operations on variables
module var_ops
    use num_vars, only: dp, max_str_ln
    implicit none
    private
    public i2str, r2str, r2strt, strh2l, mat_mult, mat_sub, det

contains
    ! Convert an integer to string 
    ! (from http://stackoverflow.com/questions/1262695/converting-integers-to-strings-in-fortran)
    character(len=max_str_ln) function i2str(k)
        integer, intent(in) :: k
        write (i2str, *) k
        i2str = adjustl(i2str)
    end function i2str

    ! Convert a real (double) to string
    character(len=max_str_ln) function r2str(k)
        real(dp), intent(in) :: k
        write (r2str, *) k
        r2str = adjustl(r2str)
    end function r2str

    ! Convert a real (double) to string and truncate it
    character(len=max_str_ln) function r2strt(k)
        real(dp), intent(in) :: k
        write (r2strt, '(E9.3)') k
        r2strt = adjustl(r2strt)
    end function r2strt

    ! Convert a string to lowercase 
    ! (from Figure 3.5B, pg 80, "Upgrading to Fortran 90", by Cooper Redwine,
    ! 1995 Springer-Verlag, New York. ) 
    function strh2l(input_string) result(output_string)
        character(*), intent(in)     :: input_string
        character(len(input_string)) :: output_string
        integer :: i, n
        character(*), parameter :: lower_case = 'abcdefghijklmnopqrstuvwxyz'
        character(*), parameter :: upper_case = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
            
        ! copy input string
        output_string = input_string
        
        ! convert case character by character
        do i = 1, len(output_string)
          n = index(upper_case, output_string(i:i))
          if ( n /= 0 ) output_string(i:i) = lower_case(n:n)
        end do
    end function strh2l

    ! multipy two matrices
    function mat_mult(A,B)
        use output_ops, only: writo
        real(dp) :: A(:,:), B(:,:)
        real(dp), allocatable :: mat_mult(:,:)
        integer :: id, jd, kd
        
        if (size(A,2).ne.size(B,1)) then
            call writo('ERROR: matrices A and B not compatible')
            stop
        end if
        
        allocate(mat_mult(size(A,1),size(B,2)))
        mat_mult = 0.0_dp
        do jd = 1,size(B,2)
            do id = 1,size(A,1)
                do kd = 1, size(A,2)
                    mat_mult(id,jd) = mat_mult(id,jd) + A(id,kd)*B(kd,jd)
                end do
            end do
        end do
    end function mat_mult

    ! subtract two matrices
    function mat_sub(A,B)
        use output_ops, only: writo
        real(dp) :: A(:,:), B(:,:)
        real(dp), allocatable :: mat_sub(:,:)
        integer :: id, jd
        
        if (size(A,1).ne.size(B,1) .or. size(A,2).ne.size(B,2)) then
            call writo('ERROR: matrices A and B not compatible')
            stop
        end if
        
        allocate(mat_sub(size(A,1),size(A,2)))
        mat_sub = 0.0_dp
        do jd = 1,size(B,2)
            do id = 1,size(A,1)
                mat_sub(id,jd) = A(id,jd)-B(id,jd)
            end do
        end do
    end function mat_sub

    ! calculate determinant of a matrix
    ! (adapted from http://dualm.wordpress.com/2012/01/06/computing-determinant-in-fortran/)
    real(dp) function det(N, mat)
        integer, intent(in) :: N 
        real(dp), intent(inout) :: mat(:,:)
        integer :: i, info
        integer, allocatable :: ipiv(:)
        real(dp) :: sgn
        
        allocate(ipiv(N))
        
        ipiv = 0
        
        call dgetrf(N, N, mat, N, ipiv, info)
        
        det = 1.0_dp
        do i = 1, N
            det = det*mat(i, i)
        end do
        
        sgn = 1.0_dp
        do i = 1, N
            if(ipiv(i) /= i) then
                sgn = -sgn
            end if
        end do
        det = sgn*det   
    end function det
end module var_ops
