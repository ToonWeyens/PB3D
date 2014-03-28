!------------------------------------------------------------------------------!
!   This module contains operations on variables other than strings            !
!------------------------------------------------------------------------------!
module var_ops
    use num_vars, only: dp, max_str_ln
    use output_ops, only: writo
    implicit none
    private
    public mat_mult, matvec_mult, mat_sub, det

contains
    ! multipy two matrices
    function mat_mult(A,B)
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
    
    ! multiply a matrix with a vector
    function matvec_mult(A,b)
        real(dp) :: A(:,:), b(:)
        real(dp), allocatable :: matvec_mult(:)
        integer :: id, kd
        
        if (size(A,2).ne.size(b)) then
            call writo('ERROR: Matrix A and vector b not compatible')
            stop
        end if
        
        allocate(matvec_mult(size(A,1)))
        matvec_mult = 0.0_dp
        do id = 1, size(A,1)
            do kd = 1, size(A,2)
                matvec_mult(id) = matvec_mult(id) + A(id,kd)*b(kd)
            end do
        end do
    end function matvec_mult

    ! subtract two matrices
    function mat_sub(A,B)
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
